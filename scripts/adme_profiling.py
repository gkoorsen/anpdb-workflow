"""SwissADME-equivalent ADME profiling for the top-50 consensus compounds.

Computes Lipinski, Veber, Egan, BOILED-Egg (GI/BBB) rules locally with RDKit
descriptors. No web API needed.

BOILED-Egg model (Daina & Zoete 2016):
  GI absorption: TPSA <= 142 and -0.4 <= WLOGP <= 5.6
  BBB permeant:  TPSA <= 79 and -1.0 <= WLOGP <= 3.5

Outputs (under output/adme/)
----------------------------
- adme_top50.csv          full descriptor table
- adme_summary.txt        pass/fail counts for each rule set
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
SHORTLIST = ROOT / "output" / "chembl_nn" / "consensus_top50_shortlist.tsv"
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"
OUT = ROOT / "output" / "adme"
OUT.mkdir(parents=True, exist_ok=True)


def compute_adme(smi: str) -> dict | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotb = Descriptors.NumRotatableBonds(mol)
    mrefr = Descriptors.MolMR(mol)
    n_atoms = mol.GetNumHeavyAtoms()
    n_rings = Descriptors.RingCount(mol)
    n_arom_rings = Lipinski.NumAromaticRings(mol)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)

    lipinski_violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10,
    ])
    veber_ok = rotb <= 10 and tpsa <= 140
    egan_ok = tpsa <= 131.6 and -1.0 <= logp <= 5.88
    gi_absorb = tpsa <= 142 and -0.4 <= logp <= 5.6
    bbb_permeant = tpsa <= 79 and -1.0 <= logp <= 3.5

    return {
        "MW": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "TPSA": round(tpsa, 2),
        "RotBonds": rotb,
        "MolRefractivity": round(mrefr, 2),
        "HeavyAtoms": n_atoms,
        "Rings": n_rings,
        "AromaticRings": n_arom_rings,
        "Fsp3": round(fsp3, 3),
        "Lipinski_violations": lipinski_violations,
        "Lipinski_pass": lipinski_violations <= 1,
        "Veber_pass": veber_ok,
        "Egan_pass": egan_ok,
        "GI_absorption": gi_absorb,
        "BBB_permeant": bbb_permeant,
    }


def main() -> int:
    shortlist = pd.read_csv(SHORTLIST, sep="\t", dtype=str, keep_default_na=False)
    compound_ids = shortlist["compound_id"].unique().tolist()
    print(f"Shortlist compounds: {len(compound_ids)}")

    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    smiles_map = dict(zip(novel["molecule_id"], novel["std_smiles"].fillna("")))

    rows: list[dict] = []
    for cid in compound_ids:
        smi = smiles_map.get(cid, "")
        if not smi:
            continue
        props = compute_adme(smi)
        if props is None:
            continue
        rows.append({"compound_id": cid, "smiles": smi, **props})

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "adme_top50.csv", index=False)

    n = len(df)
    summary = (
        f"ADME Profiling — Top 50 Consensus Compounds\n"
        f"============================================\n\n"
        f"Compounds profiled: {n}\n\n"
        f"Lipinski Ro5 (<=1 violation):  {df['Lipinski_pass'].sum()}/{n} pass "
        f"({100*df['Lipinski_pass'].mean():.0f}%)\n"
        f"Veber (RotB<=10, TPSA<=140):   {df['Veber_pass'].sum()}/{n} pass "
        f"({100*df['Veber_pass'].mean():.0f}%)\n"
        f"Egan (TPSA<=131.6, LogP ok):   {df['Egan_pass'].sum()}/{n} pass "
        f"({100*df['Egan_pass'].mean():.0f}%)\n"
        f"GI absorption (BOILED-Egg):    {df['GI_absorption'].sum()}/{n} "
        f"({100*df['GI_absorption'].mean():.0f}%)\n"
        f"BBB permeant (BOILED-Egg):     {df['BBB_permeant'].sum()}/{n} "
        f"({100*df['BBB_permeant'].mean():.0f}%)\n\n"
        f"Descriptor statistics:\n"
        f"  MW:    median={df['MW'].median():.0f}  range=[{df['MW'].min():.0f}, {df['MW'].max():.0f}]\n"
        f"  LogP:  median={df['LogP'].median():.1f}  range=[{df['LogP'].min():.1f}, {df['LogP'].max():.1f}]\n"
        f"  TPSA:  median={df['TPSA'].median():.0f}  range=[{df['TPSA'].min():.0f}, {df['TPSA'].max():.0f}]\n"
        f"  HBD:   median={df['HBD'].median():.0f}  range=[{df['HBD'].min()}, {df['HBD'].max()}]\n"
        f"  HBA:   median={df['HBA'].median():.0f}  range=[{df['HBA'].min()}, {df['HBA'].max()}]\n"
        f"  Fsp3:  median={df['Fsp3'].median():.2f}  range=[{df['Fsp3'].min():.2f}, {df['Fsp3'].max():.2f}]\n"
    )
    with open(OUT / "adme_summary.txt", "w") as fh:
        fh.write(summary)
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
