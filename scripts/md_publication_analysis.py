"""Generate publication-oriented MD analysis figures for completed MD systems.

It reuses the per-run outputs from scripts/md_analyze_production.py for RMSD,
RMSF, contact occupancy, and pose-retention descriptors, and computes additional
summary descriptors from strided trajectories to keep runtime manageable.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd

from md_analyze_production import LIPID_RESNAMES, heavy_indices, ligand_pose_retention


@dataclass(frozen=True)
class RunSpec:
    system: str
    replicate: str
    run_dir: Path
    ligand_resname: str


RUNS = [
    RunSpec("CYP1B1", "rep1", Path("md_runs/production/cyp1b1_mol11315_amber/rep1"), "LIG"),
    RunSpec("CYP1B1", "rep2", Path("md_runs/production/cyp1b1_mol11315_amber/rep2"), "LIG"),
    RunSpec("CYP1B1", "rep3", Path("md_runs/production/cyp1b1_mol11315_amber/rep3"), "LIG"),
    RunSpec("MAO-B", "rep1_rerun2", Path("md_runs/production/maob_mol14056_amber/rep1_rerun2"), "LIG"),
    RunSpec("MAO-B", "rep2", Path("md_runs/production/maob_mol14056_amber/rep2"), "LIG"),
    RunSpec("MAO-B", "rep3", Path("md_runs/production/maob_mol14056_amber/rep3"), "LIG"),
    RunSpec("SGLT2 mol13144", "rep1", Path("md_runs/production/sglt2_mol13144/rep1"), "UNK"),
    RunSpec("SGLT2 mol13144", "rep2", Path("md_runs/production/sglt2_mol13144/rep2"), "UNK"),
    RunSpec("SGLT2 mol13144", "rep3", Path("md_runs/production/sglt2_mol13144/rep3"), "UNK"),
    RunSpec("SGLT2 mol13733", "rep1", Path("md_runs/production/sglt2_mol13733/rep1"), "UNK"),
    RunSpec("SGLT2 mol13733", "rep2", Path("md_runs/production/sglt2_mol13733/rep2"), "UNK"),
    RunSpec("SGLT2 mol13733", "rep3", Path("md_runs/production/sglt2_mol13733/rep3"), "UNK"),
    RunSpec("SGLT2 mol15088", "rep1", Path("md_runs/production/sglt2_mol15088/rep1"), "UNK"),
    RunSpec("SGLT2 mol15088", "rep2", Path("md_runs/production/sglt2_mol15088/rep2"), "UNK"),
    RunSpec("SGLT2 mol15088", "rep3", Path("md_runs/production/sglt2_mol15088/rep3"), "UNK"),
    RunSpec("OPRK1 mol16614", "rep1", Path("md_runs/production/oprk1_mol16614/rep1"), "UNL"),
    RunSpec("OPRK1 mol16614", "rep2", Path("md_runs/production/oprk1_mol16614/rep2"), "UNL"),
    RunSpec("OPRK1 mol16614", "rep3", Path("md_runs/production/oprk1_mol16614/rep3"), "UNL"),
]

COLORS = {
    "CYP1B1": "#1f77b4",
    "MAO-B": "#d95f02",
    "SGLT2 mol13144": "#1b9e77",
    "SGLT2 mol13733": "#7570b3",
    "SGLT2 mol15088": "#e7298a",
    "OPRK1 mol16614": "#66a61e",
}


SYSTEMS = list(dict.fromkeys(run.system for run in RUNS))


def ensure_run_ready(run: RunSpec) -> None:
    required = [
        run.run_dir / "equilibrated.pdb",
        run.run_dir / "production.dcd",
        run.run_dir / "run_manifest.json",
        run.run_dir / "analysis" / "analysis_summary.json",
        run.run_dir / "analysis" / "rmsd_timeseries.csv",
        run.run_dir / "analysis" / "rmsf_ca.csv",
        run.run_dir / "analysis" / "contact_occupancy.csv",
        run.run_dir / "analysis" / "pose_retention_timeseries.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing required analysis inputs:\n" + "\n".join(missing))


def residue_number(label: str) -> int:
    match = re.search(r"(\d+)$", str(label))
    return int(match.group(1)) if match else -1


def load_rmsd(runs: list[RunSpec]) -> pd.DataFrame:
    frames = []
    for run in runs:
        df = pd.read_csv(run.run_dir / "analysis" / "rmsd_timeseries.csv")
        df["system"] = run.system
        df["replicate"] = run.replicate
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_rmsf(runs: list[RunSpec]) -> pd.DataFrame:
    frames = []
    for run in runs:
        df = pd.read_csv(run.run_dir / "analysis" / "rmsf_ca.csv")
        df["system"] = run.system
        df["replicate"] = run.replicate
        df["residue_index"] = df["residue"].map(residue_number)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_contacts(runs: list[RunSpec]) -> pd.DataFrame:
    frames = []
    for run in runs:
        df = pd.read_csv(run.run_dir / "analysis" / "contact_occupancy.csv")
        df["system"] = run.system
        df["replicate"] = run.replicate
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_summaries(runs: list[RunSpec]) -> pd.DataFrame:
    rows = []
    for run in runs:
        data = json.loads((run.run_dir / "analysis" / "analysis_summary.json").read_text())
        data["system"] = run.system
        data["replicate"] = run.replicate
        data["run_dir"] = str(run.run_dir)
        rows.append(data)
    return pd.DataFrame(rows)


def compute_extra_metrics(run: RunSpec, stride: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    traj = md.load(
        str(run.run_dir / "production.dcd"),
        top=str(run.run_dir / "equilibrated.pdb"),
        stride=stride,
    )
    top = traj.topology
    ligand = np.array([atom.index for atom in top.atoms if atom.residue.name == run.ligand_resname], dtype=int)
    ligand_set = set(int(index) for index in ligand)
    protein = np.array([int(index) for index in top.select("protein") if int(index) not in ligand_set], dtype=int)
    backbone = np.array([int(index) for index in top.select("protein and backbone") if int(index) not in ligand_set], dtype=int)
    if len(protein) == 0 or len(backbone) == 0 or len(ligand) == 0:
        raise SystemExit(f"Could not select protein/backbone/ligand atoms for {run.run_dir}")

    traj.superpose(traj, frame=0, atom_indices=backbone)
    times_ns = np.arange(traj.n_frames, dtype=float) * stride * 0.05

    protein_traj = traj.atom_slice(protein)
    rg_A = md.compute_rg(protein_traj) * 10.0

    sasa_traj = traj.atom_slice(np.concatenate([protein, ligand]))
    sasa_atom = md.shrake_rupley(sasa_traj, mode="atom")
    n_protein = len(protein)
    protein_sasa_nm2 = sasa_atom[:, :n_protein].sum(axis=1)
    ligand_sasa_nm2 = sasa_atom[:, n_protein:].sum(axis=1)

    ligand_local = np.arange(n_protein, n_protein + len(ligand), dtype=int)
    protein_heavy_local = np.array(
        [
            i
            for i, atom in enumerate(sasa_traj.topology.atoms)
            if i < n_protein and atom.element is not None and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    ligand_heavy_local = np.array(
        [
            i
            for i in ligand_local
            if sasa_traj.topology.atom(int(i)).element is not None
            and sasa_traj.topology.atom(int(i)).element.symbol != "H"
        ],
        dtype=int,
    )
    neighbors = md.compute_neighbors(
        sasa_traj,
        0.4,
        ligand_heavy_local,
        haystack_indices=protein_heavy_local,
    )
    contact_counts = np.array([len(frame_neighbors) for frame_neighbors in neighbors], dtype=int)
    protein_heavy = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.index not in ligand_set
            and atom.residue.is_protein
            and atom.element is not None
            and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    ligand_heavy = heavy_indices(top, ligand)
    lipid_indices = np.array(
        [
            atom.index
            for atom in top.atoms
            if atom.residue.name in LIPID_RESNAMES and atom.element is not None and atom.element.symbol != "H"
        ],
        dtype=int,
    )
    pose = ligand_pose_retention(traj, ligand_heavy, protein_heavy, lipid_indices)

    metrics = pd.DataFrame(
        {
            "system": run.system,
            "replicate": run.replicate,
            "time_ns": times_ns,
            "protein_rg_A": rg_A,
            "protein_sasa_nm2": protein_sasa_nm2,
            "ligand_sasa_nm2": ligand_sasa_nm2,
            "protein_ligand_contact_atoms": contact_counts,
        }
    )
    metrics = pd.concat([metrics, pose], axis=1)

    hb_rows = []
    try:
        hbonds = md.baker_hubbard(sasa_traj, freq=0.1, periodic=False)
    except Exception:
        hbonds = np.empty((0, 3), dtype=int)
    for donor, hydrogen, acceptor in hbonds:
        donor_atom = sasa_traj.topology.atom(int(donor))
        acceptor_atom = sasa_traj.topology.atom(int(acceptor))
        donor_lig = donor_atom.residue.name == run.ligand_resname
        acceptor_lig = acceptor_atom.residue.name == run.ligand_resname
        if donor_lig == acceptor_lig:
            continue
        hb_rows.append(
            {
                "system": run.system,
                "replicate": run.replicate,
                "donor": str(donor_atom),
                "acceptor": str(acceptor_atom),
                "ligand_role": "donor" if donor_lig else "acceptor",
            }
        )
    hbonds_df = pd.DataFrame(hb_rows)

    extra_summary = pd.DataFrame(
        [
            {
                "system": run.system,
                "replicate": run.replicate,
                "extra_stride": stride,
                "extra_frames": traj.n_frames,
                "protein_rg_A_mean": float(rg_A.mean()),
                "protein_rg_A_sd": float(rg_A.std(ddof=1)),
                "protein_sasa_nm2_mean": float(protein_sasa_nm2.mean()),
                "protein_sasa_nm2_sd": float(protein_sasa_nm2.std(ddof=1)),
                "ligand_sasa_nm2_mean": float(ligand_sasa_nm2.mean()),
                "ligand_sasa_nm2_sd": float(ligand_sasa_nm2.std(ddof=1)),
                "contact_atoms_mean": float(contact_counts.mean()),
                "contact_atoms_sd": float(contact_counts.std(ddof=1)),
                "ligand_com_displacement_A_mean": float(pose["ligand_com_displacement_A"].mean()),
                "ligand_com_displacement_A_max": float(pose["ligand_com_displacement_A"].max()),
                "ligand_pocket_com_distance_A_mean": float(pose["ligand_pocket_com_distance_A"].mean()),
                "protein_ligand_min_heavy_distance_A_mean": float(pose["protein_ligand_min_heavy_distance_A"].mean()),
                "contact_atoms_4A_mean": float(pose["contact_atoms_4A"].mean()),
                "contact_atoms_6A_mean": float(pose["contact_atoms_6A"].mean()),
                "contact_atoms_8A_mean": float(pose["contact_atoms_8A"].mean()),
                "ligand_z_from_membrane_center_A_mean": float(pose["ligand_z_from_membrane_center_A"].mean()),
                "ligand_protein_hbonds_unique_freq_ge_0p1": int(len(hbonds_df)),
            }
        ]
    )
    return metrics, hbonds_df, extra_summary


def savefig(fig: plt.Figure, out_dir: Path, name: str) -> None:
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(out_dir / f"{name}.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rmsd(rmsd: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, len(SYSTEMS), figsize=(3.4 * len(SYSTEMS), 7), sharex=True, squeeze=False)
    for col, title, row in [
        ("backbone_rmsd_A", "Backbone RMSD", 0),
        ("ligand_rmsd_A", "Ligand RMSD", 1),
    ]:
        for ax, system in zip(axes[row], SYSTEMS, strict=True):
            sub = rmsd[rmsd["system"] == system]
            for rep, rep_df in sub.groupby("replicate"):
                ax.plot(rep_df["time_ns"], rep_df[col], lw=0.9, alpha=0.85, label=rep)
            ax.set_title(system)
            ax.set_ylabel(f"{title} (A)")
            ax.grid(alpha=0.25)
            if row == 1:
                ax.set_xlabel("Time (ns)")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Protein backbone and ligand RMSD across 100 ns production replicates")
    savefig(fig, out_dir, "fig01_backbone_ligand_rmsd")


def plot_rmsd_distribution(rmsd: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(max(10, 1.6 * len(SYSTEMS)), 4))
    positions = np.arange(1, len(SYSTEMS) + 1)
    for ax, col, title in [
        (axes[0], "backbone_rmsd_A", "Backbone RMSD"),
        (axes[1], "ligand_rmsd_A", "Ligand RMSD"),
    ]:
        values = [rmsd.loc[rmsd["system"] == system, col].to_numpy() for system in SYSTEMS]
        parts = ax.violinplot(values, positions=positions, showmeans=True, showextrema=False)
        for body, system in zip(parts["bodies"], SYSTEMS, strict=True):
            body.set_facecolor(COLORS[system])
            body.set_alpha(0.35)
        parts["cmeans"].set_color("#222222")
        ax.set_xticks(positions, SYSTEMS, rotation=35, ha="right")
        ax.set_ylabel(f"{title} (A)")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)
    savefig(fig, out_dir, "fig02_rmsd_distributions")


def plot_rmsf(rmsf: pd.DataFrame, out_dir: Path) -> None:
    ncols = 2
    nrows = int(np.ceil(len(SYSTEMS) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.6 * nrows), sharey=True, squeeze=False)
    flat_axes = axes.ravel()
    for ax, system in zip(flat_axes, SYSTEMS, strict=False):
        sub = rmsf[rmsf["system"] == system]
        pivot = sub.pivot_table(index="residue_index", columns="replicate", values="rmsf_A")
        mean = pivot.mean(axis=1)
        sd = pivot.std(axis=1)
        x = mean.index.to_numpy(dtype=float)
        ax.plot(x, mean.to_numpy(), color=COLORS[system], lw=1.5)
        ax.fill_between(x, (mean - sd).to_numpy(), (mean + sd).to_numpy(), color=COLORS[system], alpha=0.2, lw=0)
        ax.set_title(system)
        ax.set_xlabel("Residue number")
        ax.grid(alpha=0.25)
    for ax in flat_axes[len(SYSTEMS):]:
        ax.axis("off")
    axes[0, 0].set_ylabel("C-alpha RMSF (A)")
    fig.suptitle("Per-residue C-alpha RMSF, mean +/- SD across replicates")
    savefig(fig, out_dir, "fig03_ca_rmsf")


def plot_contacts(contacts: pd.DataFrame, out_dir: Path) -> None:
    for system in SYSTEMS:
        sub = contacts[contacts["system"] == system]
        top = (
            sub.groupby("residue", as_index=False)["occupancy"]
            .mean()
            .sort_values("occupancy", ascending=False)
            .head(15)
            .sort_values("occupancy")
        )
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.barh(top["residue"], top["occupancy"], color=COLORS[system], alpha=0.8)
        ax.set_xlabel("Mean contact occupancy")
        ax.set_title(f"{system}: top ligand-contact residues")
        ax.set_xlim(0, 1)
        ax.grid(axis="x", alpha=0.25)
        savefig(fig, out_dir, f"fig04_{system.lower().replace('-', '')}_top_contacts")


def plot_extra_timeseries(extra: pd.DataFrame, out_dir: Path) -> None:
    plot_defs = [
        ("protein_rg_A", "Protein radius of gyration (A)", "fig05_protein_radius_of_gyration"),
        ("protein_sasa_nm2", "Protein SASA (nm2)", "fig06_protein_sasa"),
        ("ligand_sasa_nm2", "Ligand SASA (nm2)", "fig07_ligand_sasa"),
        ("protein_ligand_contact_atoms", "Protein atoms within 4 A of ligand", "fig08_ligand_contact_count"),
        ("ligand_com_displacement_A", "Ligand COM displacement from start (A)", "fig09_ligand_com_displacement"),
        ("ligand_pocket_com_distance_A", "Ligand COM to starting pocket COM (A)", "fig10_ligand_pocket_com_distance"),
        ("protein_ligand_min_heavy_distance_A", "Minimum protein-ligand heavy-atom distance (A)", "fig11_min_heavy_atom_distance"),
        ("ligand_z_from_membrane_center_A", "Ligand z from membrane center (A)", "fig12_ligand_membrane_z"),
    ]
    for column, ylabel, name in plot_defs:
        ncols = 2
        nrows = int(np.ceil(len(SYSTEMS) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, 3.4 * nrows), sharey=False, squeeze=False)
        flat_axes = axes.ravel()
        for ax, system in zip(flat_axes, SYSTEMS, strict=False):
            sub = extra[extra["system"] == system]
            for rep, rep_df in sub.groupby("replicate"):
                ax.plot(rep_df["time_ns"], rep_df[column], lw=0.9, alpha=0.85, label=rep)
            ax.set_title(system)
            ax.set_xlabel("Time (ns)")
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.25)
            ax.legend(frameon=False, fontsize=8)
        for ax in flat_axes[len(SYSTEMS):]:
            ax.axis("off")
        savefig(fig, out_dir, name)


def plot_hbond_candidates(hbonds_df: pd.DataFrame, out_dir: Path) -> None:
    if hbonds_df.empty:
        return
    counts = (
        hbonds_df.groupby(["system", "replicate"], as_index=False)
        .size()
        .rename(columns={"size": "hbond_candidates"})
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = [f"{row.system}\n{row.replicate}" for row in counts.itertuples()]
    colors = [COLORS.get(row.system, "#666666") for row in counts.itertuples()]
    ax.bar(labels, counts["hbond_candidates"], color=colors, alpha=0.8)
    ax.set_ylabel("Unique ligand-protein H-bond candidates")
    ax.set_title("H-bond candidates with >=10% trajectory frequency")
    ax.grid(axis="y", alpha=0.25)
    savefig(fig, out_dir, "fig13_hbond_candidate_counts")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("output/md_publication_all_systems"))
    parser.add_argument("--extra-stride", type=int, default=10)
    args = parser.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = out_dir / "tables"
    tables_dir.mkdir(exist_ok=True)

    for run in RUNS:
        ensure_run_ready(run)

    rmsd = load_rmsd(RUNS)
    rmsf = load_rmsf(RUNS)
    contacts = load_contacts(RUNS)
    summaries = load_summaries(RUNS)

    extra_frames = []
    hb_frames = []
    extra_summary_frames = []
    for run in RUNS:
        print(f"Computing extra metrics for {run.system} {run.replicate}...", flush=True)
        metrics, hbonds, extra_summary = compute_extra_metrics(run, args.extra_stride)
        extra_frames.append(metrics)
        hb_frames.append(hbonds)
        extra_summary_frames.append(extra_summary)

    extra = pd.concat(extra_frames, ignore_index=True)
    hbonds_df = pd.concat(hb_frames, ignore_index=True) if hb_frames else pd.DataFrame()
    extra_summary = pd.concat(extra_summary_frames, ignore_index=True)

    summaries.to_csv(tables_dir / "per_replicate_rmsd_summary.csv", index=False)
    rmsd.to_csv(tables_dir / "rmsd_timeseries_all.csv", index=False)
    rmsf.to_csv(tables_dir / "rmsf_ca_all.csv", index=False)
    contacts.to_csv(tables_dir / "contact_occupancy_all.csv", index=False)
    extra.to_csv(tables_dir / "extra_metrics_timeseries_stride.csv", index=False)
    extra_summary.to_csv(tables_dir / "extra_metrics_summary.csv", index=False)
    hbonds_df.to_csv(tables_dir / "ligand_protein_hbond_candidates.csv", index=False)

    system_summary = (
        summaries.groupby("system")
        .agg(
            n=("replicate", "count"),
            backbone_rmsd_A_mean=("backbone_rmsd_A_mean", "mean"),
            backbone_rmsd_A_sd=("backbone_rmsd_A_mean", "std"),
            ligand_rmsd_A_mean=("ligand_rmsd_A_mean", "mean"),
            ligand_rmsd_A_sd=("ligand_rmsd_A_mean", "std"),
            backbone_rmsd_A_final_mean=("backbone_rmsd_A_final", "mean"),
            ligand_rmsd_A_final_mean=("ligand_rmsd_A_final", "mean"),
            contact_residues_mean=("n_contact_residues", "mean"),
        )
        .reset_index()
    )
    system_summary.to_csv(tables_dir / "system_level_summary.csv", index=False)

    contact_summary = (
        contacts.groupby(["system", "residue"], as_index=False)
        .agg(mean_occupancy=("occupancy", "mean"), sd_occupancy=("occupancy", "std"), n_reps=("replicate", "nunique"))
        .sort_values(["system", "mean_occupancy"], ascending=[True, False])
    )
    contact_summary.to_csv(tables_dir / "contact_residue_summary.csv", index=False)

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    plot_rmsd(rmsd, out_dir)
    plot_rmsd_distribution(rmsd, out_dir)
    plot_rmsf(rmsf, out_dir)
    plot_contacts(contacts, out_dir)
    plot_extra_timeseries(extra, out_dir)
    plot_hbond_candidates(hbonds_df, out_dir)

    print(f"Wrote publication analysis to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
