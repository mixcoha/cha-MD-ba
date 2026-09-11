"""Empaquetado de mutantes 6M03 para el nodo LARCAD."""

import importlib.util
from pathlib import Path

from test_benchmark_6m03 import HIS41_CYS145_PDB


def _load_pack():
    root = Path(__file__).resolve().parents[2]
    spec = importlib.util.spec_from_file_location(
        "pack_larcad_6m03", root / "scripts" / "pack_larcad_6m03.py"
    )
    pack = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pack)
    return pack


def test_pack_larcad_writes_three_mutants(tmp_path):
    pack = _load_pack()
    source = tmp_path / "mini.pdb"
    source.write_text(HIS41_CYS145_PDB)
    cluster = tmp_path / "larcad"
    (cluster / "mdp").mkdir(parents=True)
    (cluster / "mdp" / "em.mdp").write_text("title = em\n")
    (cluster / "run_model.sh").write_text("#!/bin/bash\n")

    written = pack.prepare_models(cluster, source_pdb=source)
    assert set(written["models"]) == {"6M03_H41A", "6M03_C145A", "6M03_H41A_C145A"}

    h41 = (cluster / "models" / "6M03_H41A" / "protein.pdb").read_text()
    c145 = (cluster / "models" / "6M03_C145A" / "protein.pdb").read_text()
    both = (cluster / "models" / "6M03_H41A_C145A" / "protein.pdb").read_text()
    assert "ALA" in h41 and "CYS" in h41
    assert " ND1 " not in h41
    assert " SG " not in c145
    assert "ALA" in c145
    assert "HIS" not in both and "HIE" not in both and "CYS" not in both

    tarball = tmp_path / "bundle.tar.gz"
    pack.make_tarball(cluster, tarball)
    assert tarball.exists() and tarball.stat().st_size > 0


def test_larcad_scripts_parse_and_nacl_is_0_15():
    import subprocess

    root = Path(__file__).resolve().parents[2] / "cluster" / "larcad"
    for name in (
        "run_model.sh",
        "submit_all.sh",
        "submit_one.sh",
        "submit_model.slurm",
        "diagnose_jobs.sh",
    ):
        subprocess.check_call(["bash", "-n", str(root / name)])

    run_model = (root / "run_model.sh").read_text()
    assert "-conc 0.15" in run_model
    assert "GPU_FLAGS=(-nb gpu -pme gpu)" in run_model
    assert "run_mdrun em" in run_model

    npt = (root / "mdp" / "npt.mdp").read_text()
    npt_active = [
        line for line in npt.splitlines() if line.strip() and not line.lstrip().startswith(";")
    ]
    assert "-DPOSRES" in npt
    assert any("C-rescale" in line for line in npt_active)
    assert not any("Parrinello-Rahman" in line for line in npt_active)

    md = (root / "mdp" / "md.mdp").read_text()
    assert "Parrinello-Rahman" in md
    assert "0.15 M" in md

    slurm = (root / "submit_model.slurm").read_text()
    assert "--output=%x_%j.out" in slurm
    assert "--output=logs/" not in slurm
    for name in ("submit_all.sh", "submit_one.sh"):
        active = [
            line
            for line in (root / name).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        assert not any("--export=ALL" in line for line in active)
