"""Bemis-Murcko scaffold analysis of the 1,012 truly-novel ANPDB compounds.

Outputs (under output/scaffolds/)
---------------------------------
- scaffold_counts.csv           unique scaffolds ranked by frequency
- scaffold_per_compound.csv     per-compound scaffold + generic scaffold
- scaffold_summary.txt          top-level stats
- fig_scaffold_freq.png         top-20 most frequent scaffolds
- fig_scaffold_diversity.png    cumulative scaffold coverage curve
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

PYTHON = sys.executable
ROOT = Path("/Users/gerritkoorsen/anpdb-novelty")
NOVEL_STD = ROOT / "output" / "anpdb_truly_novel_std.csv"
OUT = ROOT / "output" / "scaffolds"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    novel = pd.read_csv(NOVEL_STD, dtype=str, keep_default_na=False)
    print(f"Loaded {len(novel):,} novel compounds")

    rows: list[dict] = []
    n_fail = 0
    for _, r in novel.iterrows():
        smi = r["std_smiles"] or r["smiles"]
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            n_fail += 1
            rows.append({
                "molecule_id": r["molecule_id"],
                "smiles": smi,
                "scaffold_smiles": "",
                "generic_scaffold_smiles": "",
            })
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(mol)
            scaffold_smi = Chem.MolToSmiles(core)
        except Exception:
            scaffold_smi = ""
        try:
            generic = MurckoScaffold.MakeScaffoldGeneric(core) if core else None
            generic_smi = Chem.MolToSmiles(generic) if generic else ""
        except Exception:
            generic_smi = ""
        rows.append({
            "molecule_id": r["molecule_id"],
            "smiles": smi,
            "scaffold_smiles": scaffold_smi,
            "generic_scaffold_smiles": generic_smi,
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "scaffold_per_compound.csv", index=False)
    print(f"  parsed: {len(df) - n_fail:,}  failed: {n_fail}")

    valid = df[df["scaffold_smiles"] != ""]
    scaffold_counts = Counter(valid["scaffold_smiles"])
    generic_counts = Counter(valid["generic_scaffold_smiles"])

    n_unique_scaffolds = len(scaffold_counts)
    n_unique_generic = len(generic_counts)
    n_singleton = sum(1 for c in scaffold_counts.values() if c == 1)
    n_singleton_gen = sum(1 for c in generic_counts.values() if c == 1)

    sc_df = pd.DataFrame(
        scaffold_counts.most_common(),
        columns=["scaffold_smiles", "n_compounds"],
    )
    sc_df["fraction"] = sc_df["n_compounds"] / len(valid)
    sc_df["cumulative_fraction"] = sc_df["n_compounds"].cumsum() / len(valid)
    sc_df.to_csv(OUT / "scaffold_counts.csv", index=False)

    gc_df = pd.DataFrame(
        generic_counts.most_common(),
        columns=["generic_scaffold_smiles", "n_compounds"],
    )
    gc_df["fraction"] = gc_df["n_compounds"] / len(valid)
    gc_df["cumulative_fraction"] = gc_df["n_compounds"].cumsum() / len(valid)
    gc_df.to_csv(OUT / "generic_scaffold_counts.csv", index=False)

    summary = (
        f"Bemis-Murcko Scaffold Analysis\n"
        f"==============================\n\n"
        f"Novel compounds analysed:       {len(valid):,} / {len(df):,}\n"
        f"Unique Murcko scaffolds:        {n_unique_scaffolds:,}\n"
        f"Unique generic frameworks:      {n_unique_generic:,}\n"
        f"Scaffold diversity (scaffolds/compounds): {n_unique_scaffolds/len(valid):.3f}\n"
        f"Generic diversity:              {n_unique_generic/len(valid):.3f}\n"
        f"Singleton scaffolds (appear once): {n_singleton:,} ({100*n_singleton/n_unique_scaffolds:.1f}%)\n"
        f"Singleton generic frameworks:   {n_singleton_gen:,} ({100*n_singleton_gen/n_unique_generic:.1f}%)\n\n"
        f"Top 10 scaffolds:\n"
    )
    for _, r in sc_df.head(10).iterrows():
        summary += f"  {r['scaffold_smiles'][:70]:<70}  n={r['n_compounds']:>3}  ({100*r['fraction']:.1f}%)\n"
    summary += f"\nTop 10 generic frameworks:\n"
    for _, r in gc_df.head(10).iterrows():
        summary += f"  {r['generic_scaffold_smiles'][:70]:<70}  n={r['n_compounds']:>3}  ({100*r['fraction']:.1f}%)\n"

    with open(OUT / "scaffold_summary.txt", "w") as fh:
        fh.write(summary)
    print(summary)

    # Save data for plotting separately (pidgin4_env lacks matplotlib)
    print(f"\nAll scaffold outputs in {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
