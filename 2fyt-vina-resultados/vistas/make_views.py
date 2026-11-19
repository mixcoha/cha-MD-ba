#!/usr/bin/env python3
"""Vistas del receptor 2FYT con la pose cristalina y las poses de Vina."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem

ROOT = Path("/agent/repos/cha-md-ba/2fyt-vina-resultados")
SRC = Path("/cursor/stores/self/internal/2fyt-vina")
OUT = ROOT / "vistas"
OUT.mkdir(parents=True, exist_ok=True)

TYPE_ELEM = {
    "SA": "S",
    "OA": "O",
    "NA": "N",
    "HD": "H",
    "HS": "H",
    "A": "C",
    "C": "C",
    "N": "N",
    "O": "O",
    "S": "S",
    "P": "P",
    "F": "F",
    "Cl": "Cl",
    "CL": "Cl",
    "Br": "Br",
    "BR": "Br",
    "I": "I",
}


def is_hydrogen(line: str) -> bool:
    name = line[12:16].strip()
    atype = line.split()[-1]
    return atype in {"H", "HD", "HS"} or name.startswith("H")


def models_of(path: Path):
    lines = path.read_text().splitlines()
    if not any(line.startswith("MODEL") for line in lines):
        atoms = [line for line in lines if line.startswith(("ATOM", "HETATM"))]
        return [(None, atoms)]
    models = []
    affinity = None
    atoms = []
    inside = False
    for line in lines:
        if line.startswith("MODEL"):
            inside = True
            affinity = None
            atoms = []
        elif line.startswith("ENDMDL"):
            models.append((affinity, atoms))
            inside = False
        elif inside:
            if line.startswith("REMARK VINA RESULT:"):
                affinity = float(line.split()[3])
            elif line.startswith(("ATOM", "HETATM")):
                atoms.append(line)
    return models


def atom_fields(line: str):
    name = line[12:16]
    resn = line[17:20]
    chain = line[21] if line[21] not in (" ", "") else "A"
    resi = line[22:26]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atype = line.split()[-1]
    elem = TYPE_ELEM.get(atype, name.strip()[:1])
    return name, resn, chain, resi, np.array([x, y, z]), elem


def to_pdb(atom_lines, hetatm=False, resname=None) -> str:
    rows = []
    serial = 0
    for line in atom_lines:
        if is_hydrogen(line):
            continue
        serial += 1
        name, resn, chain, resi, xyz, elem = atom_fields(line)
        if resname:
            resn = resname
        record = "HETATM" if hetatm else "ATOM  "
        el = elem.rjust(2)
        rows.append(
            f"{record}{serial:5d} {name} {resn} {chain}{resi}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {el}"
        )
    rows.append("END")
    return "\n".join(rows) + "\n"


def heavy_coords(atom_lines):
    return np.array([atom_fields(line)[4] for line in atom_lines if not is_hydrogen(line)])


def ca_coords(atom_lines):
    pts = []
    for line in atom_lines:
        if is_hydrogen(line):
            continue
        name = line[12:16].strip()
        if name == "CA":
            pts.append(atom_fields(line)[4])
    return np.array(pts)


def ligand_bonds():
    mol = Chem.MolFromMolFile(str(SRC / "prep" / "SAH_pH74_named.sdf"), removeHs=True)
    conf = mol.GetConformer()
    rd = np.array(
        [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    )
    pdbqt = heavy_coords(models_of(ROOT / "prep" / "SAH.pdbqt")[0][1])
    order = []
    used = set()
    for point in pdbqt:
        dist = np.linalg.norm(rd - point, axis=1)
        dist[list(used)] = np.inf
        k = int(np.argmin(dist))
        used.add(k)
        order.append(k)
    inv = np.empty(len(order), dtype=int)
    inv[np.array(order)] = np.arange(len(order))
    bonds = []
    for bond in mol.GetBonds():
        a = int(inv[bond.GetBeginAtomIdx()])
        b = int(inv[bond.GetEndAtomIdx()])
        bonds.append((a, b))
    return bonds


def draw_ligand(ax, coords, bonds, color, label, lw=1.6, s=18):
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], s=s, c=color, depthshade=False, label=label)
    for i, j in bonds:
        ax.plot(
            [coords[i, 0], coords[j, 0]],
            [coords[i, 1], coords[j, 1]],
            [coords[i, 2], coords[j, 2]],
            color=color,
            lw=lw,
        )


def save_overview(ca, ligands, bonds, path, title):
    fig = plt.figure(figsize=(11.5, 8.2), dpi=140)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color="#b0b0b0", lw=0.7, label="Receptor (traza CA)")
    for coords, color, label in ligands:
        draw_ligand(ax, coords, bonds, color, label)
    mins = ca.min(axis=0)
    maxs = ca.max(axis=0)
    center = (mins + maxs) / 2
    span = (maxs - mins).max() / 2
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    bonds = ligand_bonds()
    receptor_lines = models_of(ROOT / "prep" / "receptor.pdbqt")[0][1]
    ca = ca_coords(receptor_lines)
    receptor_pdb = to_pdb(receptor_lines, hetatm=False)
    (OUT / "receptor_heavy.pdb").write_text(receptor_pdb)

    crystal = heavy_coords(models_of(ROOT / "prep" / "SAH.pdbqt")[0][1])
    (OUT / "cristal.pdb").write_text(to_pdb(models_of(ROOT / "prep" / "SAH.pdbqt")[0][1], hetatm=True, resname="SAH"))

    selected = []
    # mode 1 of each exhaustiveness, plus the closest mode of exhaustiveness 64
    picks = [
        (16, 1, "#d62728", "Exh. 16, modo 1"),
        (256, 1, "#ff7f0e", "Exh. 256, modo 1"),
        (1024, 1, "#17becf", "Exh. 1024, modo 1"),
    ]
    pdb_blocks = {"receptor": receptor_pdb, "cristal": (OUT / "cristal.pdb").read_text()}
    ligands_for_plot = [(crystal, "#2ca02c", "Cristal")]
    for exh, mode, color, label in picks:
        aff, atoms = models_of(ROOT / "poses" / f"exh{exh}.pdbqt")[mode - 1]
        coords = heavy_coords(atoms)
        pdb = to_pdb(atoms, hetatm=True, resname="SAH")
        key = f"exh{exh}_modo{mode}"
        (OUT / f"{key}.pdb").write_text(pdb)
        pdb_blocks[key] = pdb
        rmsd = None
        import json

        data = json.loads((ROOT / "rmsd" / f"exh{exh}.json").read_text())
        for row in data["modes"]:
            if row["mode"] == mode:
                rmsd = row["rmsd_symmetry_A"]
        text = f"{label}: {aff:.2f} kcal/mol, RMSD {rmsd:.1f} Å"
        ligands_for_plot.append((coords, color, text))
        selected.append((key, label, color, aff, rmsd))

    save_overview(
        ca,
        ligands_for_plot,
        bonds,
        OUT / "receptor_y_poses.png",
        "2FYT: receptor, SAH cristalino y poses de Vina",
    )

    # All 20 modes of the recommended run, exhaustiveness 16.
    modes16 = models_of(ROOT / "poses" / "exh16.pdbqt")
    all16 = [(crystal, "#2ca02c", "Cristal")]
    combined = ["MODEL        0\n" + pdb_blocks["cristal"].replace("\nEND", "") + "ENDMDL\n"]
    for i, (aff, atoms) in enumerate(modes16, start=1):
        coords = heavy_coords(atoms)
        all16.append((coords, "#d62728", "Modos de exh. 16" if i == 1 else None))
        block = to_pdb(atoms, hetatm=True, resname="SAH")
        combined.append(f"MODEL    {i:4d}\n" + block.replace("\nEND", "") + "ENDMDL\n")
    (OUT / "exh16_todos.pdb").write_text("".join(combined) + "END\n")
    # drop None labels so legend has one entry
    plot16 = []
    labeled = False
    for coords, color, label in all16:
        if label is None and labeled:
            plot16.append((coords, color, None))
        else:
            plot16.append((coords, color, label))
            if label and label.startswith("Modos"):
                labeled = True
    # custom draw to avoid repeating legend
    fig = plt.figure(figsize=(11.5, 8.2), dpi=140)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(ca[:, 0], ca[:, 1], ca[:, 2], color="#b0b0b0", lw=0.7, label="Receptor (traza CA)")
    draw_ligand(ax, crystal, bonds, "#2ca02c", "Cristal", lw=2.2, s=22)
    first = True
    for aff, atoms in modes16:
        coords = heavy_coords(atoms)
        draw_ligand(ax, coords, bonds, "#d62728", "20 modos, exhaustiveness 16" if first else None, lw=0.8, s=8)
        first = False
    mins = ca.min(axis=0)
    maxs = ca.max(axis=0)
    center = (mins + maxs) / 2
    span = (maxs - mins).max() / 2
    ax.set_xlim(center[0] - span, center[0] + span)
    ax.set_ylim(center[1] - span, center[1] + span)
    ax.set_zlim(center[2] - span, center[2] + span)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=18, azim=-55)
    ax.set_xlabel("x (Å)")
    ax.set_ylabel("y (Å)")
    ax.set_zlabel("z (Å)")
    ax.set_title("Exhaustiveness 16: las 20 poses sobre el receptor")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "exh16_todas.png", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    write_html(pdb_blocks, selected)
    print("wrote", OUT)


def write_html(pdb_blocks, selected):
    js = Path("/tmp/3dmol/3Dmol-min.js").read_text()
    # Escape for a JS template literal.
    def js_string(text):
        return text.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")

    models_js = ",\n".join(
        f'"{key}": `{js_string(pdb_blocks[key])}`'
        for key in ["receptor", "cristal"] + [row[0] for row in selected]
    )
    legend = [
        ("#2ca02c", "Cristal, SAH A 549"),
    ]
    for key, label, color, aff, rmsd in selected:
        legend.append((color, f"{label}: {aff:.3f} kcal/mol, RMSD {rmsd:.2f} Å"))
    legend_html = "\n".join(
        f'<div class="item"><span class="swatch" style="background:{color}"></span>{text}</div>'
        for color, text in legend
    )
    adds = "\n".join(
        f"""
  viewer.addModel(models["{key}"], "pdb");
  viewer.setStyle({{model: {i+1}}}, {{stick: {{color: "{color}", radius: 0.18}}}});
"""
        for i, (key, label, color, aff, rmsd) in enumerate(selected)
    )
    html = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<title>2FYT, poses de Vina en el receptor</title>
<style>
  body {{ margin: 0; font-family: sans-serif; color: #222; }}
  #bar {{ padding: 10px 16px; background: #f6f6f6; border-bottom: 1px solid #ddd; }}
  h1 {{ font-size: 18px; margin: 0 0 6px; }}
  #legend {{ display: flex; flex-wrap: wrap; gap: 10px 18px; font-size: 13px; }}
  .item {{ display: flex; align-items: center; gap: 6px; }}
  .swatch {{ width: 12px; height: 12px; border-radius: 2px; display: inline-block; }}
  #viewer {{ width: 100vw; height: calc(100vh - 92px); }}
</style>
</head>
<body>
<div id="bar">
  <h1>2FYT: receptor, pose cristalina y poses de AutoDock Vina</h1>
  <div id="legend">
    {legend_html}
  </div>
</div>
<div id="viewer"></div>
<script>{js}</script>
<script>
const models = {{
{models_js}
}};
const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
viewer.addModel(models["receptor"], "pdb");
viewer.setStyle({{model: 0}}, {{cartoon: {{color: "#d0d0d0"}}}});
viewer.addModel(models["cristal"], "pdb");
viewer.setStyle({{model: 1}}, {{stick: {{color: "#2ca02c", radius: 0.22}}}});
{adds}
viewer.zoomTo();
viewer.render();
window.viewerReady = true;
</script>
</body>
</html>
"""
    (OUT / "poses.html").write_text(html)


if __name__ == "__main__":
    main()
