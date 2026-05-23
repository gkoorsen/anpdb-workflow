"""Generate ligand drug-likeness and ADMET-triage summary tables.

This is a local, reproducible triage from standardized SMILES. It uses RDKit
descriptors, rule-of-thumb filters, and structural alert catalogs. It does not
replace experimental ADMET or a validated endpoint-specific ML predictor.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, Lipinski, QED, rdMolDescriptors


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "md_inputs" / "anpdb_truly_novel_std.csv"
OUT_DIR = ROOT / "results" / "ligand_triage"


@dataclass(frozen=True)
class Candidate:
    target: str
    compound_id: str
    md_status: str


CANDIDATES = [
    Candidate("CYP1B1", "Mol_11315", "3 x 100 ns complete"),
    Candidate("SGLT2", "Mol_13144", "rep1/rep2 complete; rep3 pending"),
    Candidate("MAO-B", "Mol_14056", "3 x 100 ns complete"),
]


def yesno(value: bool) -> str:
    return "yes" if value else "no"


def filter_matches(mol: Chem.Mol, catalog_name: int) -> list[str]:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(catalog_name)
    catalog = FilterCatalog.FilterCatalog(params)
    entries = catalog.GetMatches(mol)
    return sorted({entry.GetDescription() for entry in entries})


def esol_log_s(mol: Chem.Mol, mw: float, logp: float, rot_bonds: int) -> float:
    """Delaney ESOL approximation from simple RDKit descriptors."""

    heavy = mol.GetNumHeavyAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    aromatic_proportion = aromatic_atoms / heavy if heavy else 0.0
    return 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rot_bonds) - (0.74 * aromatic_proportion)


def solubility_class(log_s: float) -> str:
    if log_s >= -2:
        return "soluble"
    if log_s >= -4:
        return "moderately soluble"
    if log_s >= -6:
        return "poorly soluble"
    return "very poorly soluble"


def triage_label(row: dict[str, object]) -> str:
    if row["rdkit_parse_ok"] != "yes":
        return "fail: invalid structure"
    if row["pains_alert_count"] > 0:
        return "deprioritize: PAINS alert"
    if row["lipinski_pass"] == "no" or row["veber_pass"] == "no":
        return "caution: oral rule violation"
    if row["brenk_alert_count"] >= 3:
        return "caution: multiple structural alerts"
    if row["qed"] >= 0.50 and row["esol_class"] in {"soluble", "moderately soluble"}:
        return "favorable"
    return "acceptable: review alerts/properties"


def analyse_candidate(candidate: Candidate, source: pd.DataFrame) -> dict[str, object]:
    hit = source.loc[source["molecule_id"] == candidate.compound_id]
    if hit.empty:
        return {
            "target": candidate.target,
            "compound_id": candidate.compound_id,
            "rdkit_parse_ok": "no",
            "triage_label": "fail: missing source row",
        }

    record = hit.iloc[0].to_dict()
    smiles = str(record.get("std_smiles", "")).strip()
    mol = Chem.MolFromSmiles(smiles)

    base = {
        "target": candidate.target,
        "compound_id": candidate.compound_id,
        "compound_name": record.get("mol_name", ""),
        "md_status": candidate.md_status,
        "std_smiles": smiles,
        "truly_novel_std": record.get("truly_novel_std", ""),
        "nearest_coconut_tanimoto": record.get("max_tanimoto", ""),
    }

    if mol is None:
        base.update({"rdkit_parse_ok": "no", "triage_label": "fail: invalid structure"})
        return base

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot_bonds = Lipinski.NumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    frac_csp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    mr = Crippen.MolMR(mol)
    qed = QED.qed(mol)
    log_s = esol_log_s(mol, mw, logp, rot_bonds)
    pains = filter_matches(mol, FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    brenk = filter_matches(mol, FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)

    lipinski_violations = sum(
        [
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10,
        ]
    )
    row: dict[str, object] = {
        **base,
        "rdkit_parse_ok": "yes",
        "mw_Da": round(mw, 2),
        "clogp": round(logp, 2),
        "tpsa_A2": round(tpsa, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot_bonds,
        "heavy_atoms": heavy_atoms,
        "rings": rings,
        "aromatic_rings": aromatic_rings,
        "fraction_csp3": round(frac_csp3, 2),
        "formal_charge": formal_charge,
        "molar_refractivity": round(mr, 2),
        "qed": round(qed, 3),
        "esol_logS": round(log_s, 2),
        "esol_class": solubility_class(log_s),
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": yesno(lipinski_violations == 0),
        "veber_pass": yesno(tpsa <= 140 and rot_bonds <= 10),
        "egan_pass": yesno(tpsa <= 131.6 and logp <= 5.88),
        "ghose_pass": yesno(160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 40 <= mr <= 130 and 20 <= heavy_atoms <= 70),
        "lead_like_pass": yesno(mw <= 350 and logp <= 3.5 and rot_bonds <= 7),
        "rough_high_gi_absorption": yesno(mw <= 500 and logp <= 5 and tpsa <= 140 and rot_bonds <= 10),
        "rough_bbb_permeable": yesno(mw <= 400 and 0 <= logp <= 6 and tpsa <= 90 and hbd <= 3),
        "pains_alert_count": len(pains),
        "pains_alerts": "; ".join(pains),
        "brenk_alert_count": len(brenk),
        "brenk_alerts": "; ".join(brenk),
        "endpoint_prediction_note": "No validated CYP/hERG/DILI/clearance predictor run; table uses descriptors and structural alerts only.",
    }
    row["triage_label"] = triage_label(row)
    return row


def write_markdown(summary: pd.DataFrame, path: Path) -> None:
    cols = [
        "target",
        "compound_id",
        "compound_name",
        "mw_Da",
        "clogp",
        "tpsa_A2",
        "hbd",
        "hba",
        "rotatable_bonds",
        "qed",
        "esol_class",
        "lipinski_pass",
        "veber_pass",
        "pains_alert_count",
        "brenk_alert_count",
        "rough_high_gi_absorption",
        "rough_bbb_permeable",
        "triage_label",
    ]
    visible = summary[cols].copy()
    table_rows = [
        "| " + " | ".join(visible.columns) + " |",
        "| " + " | ".join(["---"] * len(visible.columns)) + " |",
    ]
    for _, row in visible.iterrows():
        values = [str(row[col]).replace("|", "/") for col in visible.columns]
        table_rows.append("| " + " | ".join(values) + " |")

    text = [
        "# Ligand Triage Summary",
        "",
        "Descriptor and structural-alert triage for the MD candidate ligands.",
        "This table is not a validated endpoint-specific ADMET prediction; CYP inhibition, hERG, DILI, clearance, and transporter liability still need dedicated prediction or experimental follow-up.",
        "",
        "\n".join(table_rows),
        "",
        "Files:",
        "",
        "- `ligand_triage_summary.csv`: compact table for manuscript/supporting summaries.",
        "- `ligand_triage_full.csv`: full descriptor/rule/alert table including SMILES and alert names.",
        "",
    ]
    path.write_text("\n".join(text))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(INPUT, keep_default_na=False)
    rows = [analyse_candidate(candidate, source) for candidate in CANDIDATES]
    full = pd.DataFrame(rows)

    full_path = OUT_DIR / "ligand_triage_full.csv"
    summary_path = OUT_DIR / "ligand_triage_summary.csv"
    md_path = OUT_DIR / "README.md"

    full.to_csv(full_path, index=False, quoting=csv.QUOTE_MINIMAL)

    summary_cols = [
        "target",
        "compound_id",
        "compound_name",
        "md_status",
        "mw_Da",
        "clogp",
        "tpsa_A2",
        "hbd",
        "hba",
        "rotatable_bonds",
        "qed",
        "esol_logS",
        "esol_class",
        "lipinski_pass",
        "lipinski_violations",
        "veber_pass",
        "egan_pass",
        "ghose_pass",
        "lead_like_pass",
        "rough_high_gi_absorption",
        "rough_bbb_permeable",
        "pains_alert_count",
        "brenk_alert_count",
        "triage_label",
    ]
    full[summary_cols].to_csv(summary_path, index=False)
    write_markdown(full, md_path)

    print(f"Wrote {summary_path.relative_to(ROOT)}")
    print(f"Wrote {full_path.relative_to(ROOT)}")
    print(f"Wrote {md_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
