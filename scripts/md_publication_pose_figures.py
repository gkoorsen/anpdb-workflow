"""Generate publication pose figures and LigPlot-style interaction schematics."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import AllChem


@dataclass(frozen=True)
class PoseSpec:
    system: str
    compound_id: str
    replicate: str
    run_dir: Path
    ligand_resname: str


POSES = [
    PoseSpec("MAO-B", "Mol_14056", "rep1_rerun2", Path("md_runs/production/maob_mol14056_amber/rep1_rerun2"), "LIG"),
    PoseSpec("SGLT2", "Mol_13144", "rep1", Path("md_runs/production/sglt2_mol13144/rep1"), "UNK"),
    PoseSpec("SGLT2", "Mol_13733", "rep1", Path("md_runs/production/sglt2_mol13733/rep1"), "UNK"),
    PoseSpec("SGLT2", "Mol_15088", "rep2", Path("md_runs/production/sglt2_mol15088/rep2"), "UNK"),
    PoseSpec("OPRK1", "Mol_16614", "rep3", Path("md_runs/production/oprk1_mol16614/rep3"), "UNK"),
]


ELEMENT_COLORS = {
    "C": "#4a4a4a",
    "N": "#2b6cb0",
    "O": "#c53030",
    "S": "#b7791f",
    "P": "#805ad5",
    "F": "#2f855a",
    "CL": "#2f855a",
    "BR": "#975a16",
}

COVALENT_RADII_NM = {
    "C": 0.076,
    "N": 0.071,
    "O": 0.066,
    "S": 0.105,
    "P": 0.107,
    "F": 0.057,
    "CL": 0.102,
    "BR": 0.120,
}


def element_symbol(atom) -> str:
    if atom.element is None:
        return atom.name[0].upper()
    return atom.element.symbol.upper()


def nonhydrogen(atom) -> bool:
    return atom.element is None or atom.element.symbol.upper() != "H"


def read_ligand_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df.set_index("molecule_id")


def ligand_name(ligands: pd.DataFrame, compound_id: str) -> str:
    if compound_id in ligands.index:
        return str(ligands.loc[compound_id, "mol_name"])
    return compound_id


def ligand_smiles(ligands: pd.DataFrame, compound_id: str) -> str | None:
    if compound_id in ligands.index:
        return str(ligands.loc[compound_id, "std_smiles"])
    return None


def heavy_indices(topology, indices: np.ndarray) -> np.ndarray:
    selected = set(int(index) for index in indices)
    return np.array([atom.index for atom in topology.atoms if atom.index in selected and nonhydrogen(atom)], dtype=int)


def ligand_bonds(topology, ligand_heavy: np.ndarray, xyz: np.ndarray) -> list[tuple[int, int]]:
    heavy_set = set(int(index) for index in ligand_heavy)
    bonds = []
    for bond in topology.bonds:
        a = int(bond.atom1.index)
        b = int(bond.atom2.index)
        if a in heavy_set and b in heavy_set:
            bonds.append((a, b))
    if bonds:
        return bonds

    ligand = list(map(int, ligand_heavy))
    for pos, atom_i in enumerate(ligand):
        sym_i = element_symbol(topology.atom(atom_i))
        for atom_j in ligand[pos + 1 :]:
            sym_j = element_symbol(topology.atom(atom_j))
            cutoff = COVALENT_RADII_NM.get(sym_i, 0.077) + COVALENT_RADII_NM.get(sym_j, 0.077) + 0.045
            if np.linalg.norm(xyz[atom_i] - xyz[atom_j]) <= cutoff:
                bonds.append((atom_i, atom_j))
    return bonds


def residue_contact_rows(traj, ligand_heavy: np.ndarray, protein_heavy: np.ndarray, occupancy: dict[str, float]) -> list[dict]:
    xyz = traj.xyz[0]
    rows_by_residue: dict[str, dict] = {}
    ligand_xyz = xyz[ligand_heavy]
    for atom_index in protein_heavy:
        atom = traj.topology.atom(int(atom_index))
        residue = str(atom.residue)
        deltas = ligand_xyz - xyz[int(atom_index)]
        dists = np.sqrt((deltas * deltas).sum(axis=1)) * 10.0
        nearest = int(np.argmin(dists))
        distance = float(dists[nearest])
        if distance > 4.5:
            continue
        ligand_atom = traj.topology.atom(int(ligand_heavy[nearest]))
        polar = element_symbol(atom) in {"N", "O", "S"} and element_symbol(ligand_atom) in {"N", "O", "S"}
        hbond_like = polar and distance <= 3.5
        current = rows_by_residue.get(residue)
        if current is None or distance < current["min_distance_A"]:
            rows_by_residue[residue] = {
                "residue": residue,
                "min_distance_A": distance,
                "nearest_ligand_atom_index": int(ligand_heavy[nearest]),
                "nearest_ligand_atom_name": ligand_atom.name,
                "nearest_protein_atom_name": atom.name,
                "hbond_like": hbond_like,
                "occupancy": occupancy.get(residue, np.nan),
            }
        elif hbond_like:
            current["hbond_like"] = True
    rows = list(rows_by_residue.values())
    rows.sort(key=lambda row: (-(row["occupancy"] if not math.isnan(row["occupancy"]) else 0.0), row["min_distance_A"]))
    return rows


def load_occupancy(run_dir: Path) -> dict[str, float]:
    path = run_dir / "analysis" / "contact_occupancy.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    return {str(row.residue): float(row.occupancy) for row in df.itertuples()}


def set_equal_3d(ax, coords: np.ndarray) -> None:
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    radius = max(float((maxs - mins).max()) / 2.0, 0.5)
    for setter, value in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center, strict=False):
        setter(value - radius, value + radius)


def pca_project(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = coords.mean(axis=0)
    centered = coords - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axes = vh[:2].T
    projected = centered @ axes
    return projected, center, axes


def project_coords(coords: np.ndarray, center: np.ndarray, axes: np.ndarray) -> np.ndarray:
    return (coords - center) @ axes


def plot_binding_pose(spec: PoseSpec, traj, ligand_heavy: np.ndarray, contact_rows: list[dict], out_dir: Path, title: str) -> Path:
    xyz_A = traj.xyz[0] * 10.0
    top = traj.topology
    ligand_set = set(int(index) for index in ligand_heavy)
    contact_residues = {row["residue"] for row in contact_rows[:12]}
    ligand_xyz = xyz_A[ligand_heavy]
    pocket_atoms = []
    for atom in top.atoms:
        if (
            str(atom.residue) in contact_residues
            and atom.index not in ligand_set
            and nonhydrogen(atom)
            and atom.residue.name != spec.ligand_resname
        ):
            dists = np.sqrt(((ligand_xyz - xyz_A[atom.index]) ** 2).sum(axis=1))
            if float(dists.min()) <= 5.0:
                pocket_atoms.append(atom.index)

    all_indices = list(ligand_heavy) + pocket_atoms if pocket_atoms else list(ligand_heavy)
    projected, center, axes = pca_project(xyz_A[all_indices])
    projected_by_atom = {atom_index: projected[pos] for pos, atom_index in enumerate(all_indices)}

    fig, ax = plt.subplots(figsize=(6.7, 5.3))
    ax.set_facecolor("white")
    ax.axis("off")
    if pocket_atoms:
        coords = np.array([projected_by_atom[int(atom_index)] for atom_index in pocket_atoms])
        ax.scatter(coords[:, 0], coords[:, 1], s=18, c="#b8b8b8", alpha=0.38, edgecolor="none", zorder=1)
    for atom_i, atom_j in ligand_bonds(top, ligand_heavy, traj.xyz[0]):
        coords = np.array([projected_by_atom[int(atom_i)], projected_by_atom[int(atom_j)]])
        ax.plot(coords[:, 0], coords[:, 1], color="#111111", lw=2.8, solid_capstyle="round", zorder=4)
    for atom_index in ligand_heavy:
        atom = top.atom(int(atom_index))
        sym = element_symbol(atom)
        color = ELEMENT_COLORS.get(sym, "#4a4a4a")
        xy = projected_by_atom[int(atom_index)]
        ax.scatter(xy[0], xy[1], s=92, c=color, edgecolor="white", linewidth=0.8, zorder=5)
        if sym != "C":
            ax.text(xy[0], xy[1], sym, fontsize=7, color="white", ha="center", va="center", zorder=6)

    ligand_center = project_coords(xyz_A[ligand_heavy].mean(axis=0, keepdims=True), center, axes)[0]
    for row in contact_rows[:8]:
        residue_atoms = []
        for atom in top.atoms:
            if str(atom.residue) == row["residue"] and nonhydrogen(atom):
                dists = np.sqrt(((ligand_xyz - xyz_A[atom.index]) ** 2).sum(axis=1))
                if float(dists.min()) <= 5.0:
                    residue_atoms.append(atom.index)
        if not residue_atoms:
            continue
        residue_center = project_coords(xyz_A[residue_atoms].mean(axis=0, keepdims=True), center, axes)[0]
        start = projected_by_atom[row["nearest_ligand_atom_index"]]
        ax.plot(
            [start[0], residue_center[0]],
            [start[1], residue_center[1]],
            color="#2b6cb0" if row["hbond_like"] else "#dd6b20",
            lw=1.2,
            ls="--",
            alpha=0.72,
            zorder=2,
        )
        direction = residue_center - ligand_center
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
        label_pos = ligand_center + 0.98 * (residue_center - ligand_center) + 0.12 * direction
        ax.text(label_pos[0], label_pos[1], row["residue"], fontsize=8.5, color="#222222", ha="center", va="center")

    all_plot = np.array(list(projected_by_atom.values()))
    mins = all_plot.min(axis=0)
    maxs = all_plot.max(axis=0)
    span = np.maximum(maxs - mins, 1.0)
    pad = 0.18 * span
    ax.set_xlim(mins[0] - pad[0], maxs[0] + pad[0])
    ax.set_ylim(mins[1] - pad[1], maxs[1] + pad[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=12, pad=10)
    fig.tight_layout()
    path = out_dir / f"binding_pose_{spec.compound_id.lower()}_{spec.replicate}.png"
    for ext in ("png", "svg", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=450, bbox_inches="tight")
    plt.close(fig)
    return path


def rdkit_2d_coords(smiles: str | None) -> tuple[Chem.Mol | None, np.ndarray | None]:
    if not smiles:
        return None, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
    coords -= coords.mean(axis=0)
    scale = max(float(np.abs(coords).max()), 1.0)
    coords /= scale
    return mol, coords


def plot_interaction_diagram(
    spec: PoseSpec,
    mol,
    coords_2d: np.ndarray | None,
    ligand_heavy: np.ndarray,
    contact_rows: list[dict],
    out_dir: Path,
    title: str,
) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.set_facecolor("white")
    ax.axis("off")
    if mol is None or coords_2d is None:
        ax.text(0.0, 0.0, spec.compound_id, fontsize=18, ha="center", va="center")
        ligand_atom_positions = {int(atom_index): np.array([0.0, 0.0]) for atom_index in ligand_heavy}
    else:
        for bond in mol.GetBonds():
            begin = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            xy = coords_2d[[begin, end]]
            ax.plot(xy[:, 0], xy[:, 1], color="#202020", lw=2.0, solid_capstyle="round")
            if bond.GetBondTypeAsDouble() == 2:
                offset = np.array([-(xy[1, 1] - xy[0, 1]), xy[1, 0] - xy[0, 0]])
                norm = np.linalg.norm(offset)
                if norm > 0:
                    offset = 0.025 * offset / norm
                    ax.plot(xy[:, 0] + offset[0], xy[:, 1] + offset[1], color="#202020", lw=1.0)
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            sym = atom.GetSymbol().upper()
            color = ELEMENT_COLORS.get(sym, "#4a4a4a")
            ax.scatter(coords_2d[idx, 0], coords_2d[idx, 1], s=120 if sym != "C" else 42, c=color, edgecolor="white", zorder=4)
            if sym != "C":
                ax.text(coords_2d[idx, 0], coords_2d[idx, 1], sym, fontsize=8, color="white", ha="center", va="center", zorder=5)
        heavy_count = min(len(ligand_heavy), len(coords_2d))
        ligand_atom_positions = {int(ligand_heavy[i]): coords_2d[i] for i in range(heavy_count)}

    shown = contact_rows[:8]
    if shown:
        angles = np.linspace(0, 2 * np.pi, len(shown), endpoint=False)
        angles += np.pi / 8.0
        radius_x, radius_y = 1.78, 1.32
        for angle, row in zip(angles, shown, strict=False):
            anchor = ligand_atom_positions.get(row["nearest_ligand_atom_index"], np.array([0.0, 0.0]))
            label_xy = np.array([radius_x * np.cos(angle), radius_y * np.sin(angle)])
            color = "#2b6cb0" if row["hbond_like"] else "#dd6b20"
            ax.plot([anchor[0], label_xy[0]], [anchor[1], label_xy[1]], color=color, lw=1.2, ls="--", alpha=0.75)
            occupancy = row["occupancy"]
            occ_text = f"{occupancy:.2f}" if not math.isnan(occupancy) else "n/a"
            text = f"{row['residue']}\n{row['min_distance_A']:.1f} A | occ {occ_text}"
            ax.text(
                label_xy[0],
                label_xy[1],
                text,
                ha="center",
                va="center",
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": color, "linewidth": 1.0},
            )
    ax.text(-1.95, -1.55, "blue: polar contact <=3.5 A; orange: close contact <=4.5 A", fontsize=8, color="#333333")
    ax.set_xlim(-2.1, 2.1)
    ax.set_ylim(-1.7, 1.7)
    ax.set_title(title, fontsize=12, pad=8)
    fig.tight_layout()
    path = out_dir / f"interaction_diagram_{spec.compound_id.lower()}_{spec.replicate}.png"
    for ext in ("png", "svg", "pdf"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=450, bbox_inches="tight")
    plt.close(fig)
    return path


def trim_whitespace(image: Image.Image, pad: int = 42) -> Image.Image:
    arr = np.asarray(image.convert("RGB"))
    nonwhite = np.any(arr < 248, axis=2)
    if not nonwhite.any():
        return image
    ys, xs = np.where(nonwhite)
    left = max(int(xs.min()) - pad, 0)
    upper = max(int(ys.min()) - pad, 0)
    right = min(int(xs.max()) + pad + 1, image.width)
    lower = min(int(ys.max()) + pad + 1, image.height)
    return image.crop((left, upper, right, lower))


def make_panel(image_paths: list[Path], out_path: Path, columns: int = 2) -> None:
    images = [trim_whitespace(Image.open(path).convert("RGB")) for path in image_paths]
    width = max(image.width for image in images)
    height = max(image.height for image in images)
    rows = math.ceil(len(images) / columns)
    gutter = 60
    canvas = Image.new("RGB", (columns * width + (columns + 1) * gutter, rows * height + (rows + 1) * gutter), "white")
    labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 48)
    except OSError:
        font = ImageFont.load_default()
    draw = ImageDraw.Draw(canvas)
    for idx, image in enumerate(images):
        x0 = gutter + (idx % columns) * (width + gutter)
        y0 = gutter + (idx // columns) * (height + gutter)
        x = x0 + (width - image.width) // 2
        y = y0 + (height - image.height) // 2
        canvas.paste(image, (x, y))
        draw.text((x0 + 10, y0 + 4), labels[idx], fill="black", font=font)
    canvas.save(out_path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("results/md_analysis/publication_pose_figures"))
    parser.add_argument("--ligands", type=Path, default=Path("data/md_inputs/anpdb_truly_novel_std.csv"))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ligands = read_ligand_table(args.ligands)
    rows = []
    pose_paths = []
    diagram_paths = []
    for spec in POSES:
        pdb = spec.run_dir / "final.pdb"
        if not pdb.exists():
            raise SystemExit(f"Missing final structure: {pdb}")
        traj = md.load(str(pdb))
        top = traj.topology
        ligand = np.array([atom.index for atom in top.atoms if atom.residue.name == spec.ligand_resname], dtype=int)
        ligand_heavy = heavy_indices(top, ligand)
        protein_heavy = np.array(
            [
                atom.index
                for atom in top.atoms
                if atom.residue.is_protein and atom.residue.name != spec.ligand_resname and nonhydrogen(atom)
            ],
            dtype=int,
        )
        occupancy = load_occupancy(spec.run_dir)
        contacts = residue_contact_rows(traj, ligand_heavy, protein_heavy, occupancy)
        name = ligand_name(ligands, spec.compound_id)
        title = f"{spec.system} / {spec.compound_id} ({spec.replicate})"
        pose_paths.append(plot_binding_pose(spec, traj, ligand_heavy, contacts, args.out_dir, title))
        mol, coords_2d = rdkit_2d_coords(ligand_smiles(ligands, spec.compound_id))
        diagram_paths.append(plot_interaction_diagram(spec, mol, coords_2d, ligand_heavy, contacts, args.out_dir, title))
        for rank, contact in enumerate(contacts[:12], start=1):
            rows.append(
                {
                    "system": spec.system,
                    "compound_id": spec.compound_id,
                    "compound_name": name,
                    "replicate": spec.replicate,
                    "rank": rank,
                    **contact,
                    "ligand_heavy_atoms": len(ligand_heavy),
                    "rdkit_heavy_atoms": mol.GetNumAtoms() if mol is not None else np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(args.out_dir / "interaction_diagram_contacts.csv", index=False)
    make_panel(pose_paths, args.out_dir / "binding_pose_panel.png", columns=2)
    make_panel(diagram_paths, args.out_dir / "interaction_diagram_panel.png", columns=2)
    print(f"Wrote pose figures to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
