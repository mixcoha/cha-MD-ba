"""Pruebas del benchmark 6M03 y de las condiciones NaCl 0.5 M / 310 K."""

from unittest.mock import Mock

from cha_md_ba.benchmark import (
    Benchmark6M03Config,
    parse_system_composition,
    write_protocol_mdps,
)
from cha_md_ba.gmx_utils import find_gmx
from cha_md_ba.npt import NPTEquilibrator
from cha_md_ba.nvt import NVTEquilibrator
from cha_md_ba.pdb_utils import clean_protein_pdb
from cha_md_ba.prepare import MDSystemPreparator


SAMPLE_PDB = """HEADER    TEST
TITLE     MINI PROTEIN
ATOM      1  N   SER A   1      11.000  12.000  13.000  1.00 10.00           N
ATOM      2  CA  SER A   1      12.000  12.000  13.000  1.00 10.00           C
ATOM      3  N   HIS A   2      13.000  12.000  13.000  1.00 10.00           N
ATOM      4  CA  HIS A   2      14.000  12.000  13.000  1.00 10.00           C
HETATM    5  O   HOH A 101      20.000  20.000  20.000  1.00 30.00           O
TER       6      HIS A   2
END
"""


def test_benchmark_default_conditions():
    config = Benchmark6M03Config()
    assert config.pdb_id == "6M03"
    assert config.temperature == 310.0
    assert config.ion_concentration == 0.5
    assert config.pressure == 1.0
    assert config.forcefield == "amber99sb-ildn"
    assert config.water_model == "tip3p"


def test_clean_protein_pdb_strips_water_and_renames_his(tmp_path):
    raw = tmp_path / "raw.pdb"
    clean = tmp_path / "clean.pdb"
    raw.write_text(SAMPLE_PDB)

    stats = clean_protein_pdb(str(raw), str(clean))
    text = clean.read_text()

    assert stats["atom_kept"] == 4
    assert stats["hetatm_skipped"] == 1
    assert stats["n_residues"] == 2
    assert "HOH" not in text
    assert "HIE" in text
    assert "HIS" not in text


def test_genion_includes_salt_concentration(tmp_path, monkeypatch):
    pdb = tmp_path / "prot.pdb"
    pdb.write_text(
        "HEADER TEST\n"
        "ATOM      1  N   SER A   1      1.000   1.000   1.000  1.00  0.00           N\n"
        "END\n"
    )
    captured = []

    def fake_run(cmd, **kwargs):
        captured.append(list(cmd))
        return Mock(returncode=0)

    monkeypatch.setattr("cha_md_ba.prepare.subprocess.run", fake_run)

    preparator = MDSystemPreparator(
        pdb_path=str(pdb),
        ion_concentration=0.5,
        gmx="gmx_mpi",
    )
    preparator.prepare_system(output_dir=str(tmp_path), minimize=False, ions=True)

    genion_cmds = [cmd for cmd in captured if len(cmd) > 1 and cmd[1] == "genion"]
    assert genion_cmds, "No se invocó genion"
    genion = genion_cmds[0]
    assert "-conc" in genion
    assert "0.5" in genion
    assert "-neutral" in genion
    assert genion[0] == "gmx_mpi"

    pdb2gmx_cmds = [cmd for cmd in captured if len(cmd) > 1 and cmd[1] == "pdb2gmx"]
    assert "-ignh" in pdb2gmx_cmds[0]
    assert "-i" in pdb2gmx_cmds[0]


def test_nvt_and_npt_mdp_use_310k(tmp_path):
    gro = tmp_path / "x.gro"
    top = tmp_path / "t.top"
    gro.write_text("")
    top.write_text("")

    nvt = NVTEquilibrator(str(gro), str(top), temperature=310, gmx="gmx")
    nvt_mdp = nvt.create_mdp_file(tmp_path / "nvt.mdp")
    nvt_text = nvt_mdp.read_text()
    assert "ref_t               = 310       310" in nvt_text
    assert "gen_temp            = 310" in nvt_text
    assert "-DPOSRES" in nvt_text

    npt = NPTEquilibrator(str(gro), str(top), temperature=310, gmx="gmx")
    npt_mdp = npt.create_mdp_file(tmp_path / "npt.mdp", gen_vel=False, continuation=True)
    npt_text = npt_mdp.read_text()
    assert "ref_t               = 310       310" in npt_text
    assert "gen_temp            = 310" in npt_text
    assert "ref_p               = 1" in npt_text


def test_write_protocol_mdps(tmp_path):
    paths = write_protocol_mdps(tmp_path, Benchmark6M03Config())
    md_text = paths["md"].read_text()
    assert "310" in md_text
    assert "0.5 M NaCl" in md_text
    assert "nsteps              = 5000000" in md_text


def test_parse_system_composition(tmp_path):
    top = tmp_path / "topol.top"
    top.write_text(
        """[ molecules ]
; Compound        #mols
Protein_chain_A     1
SOL             18420
NA                175
CL                171
"""
    )
    composition = parse_system_composition(top)
    assert composition["SOL"] == 18420
    assert composition["NA"] == 175
    assert composition["CL"] == 171


def test_find_gmx_respects_env(tmp_path, monkeypatch):
    fake = tmp_path / "gmx"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("CHA_MD_BA_GMXBIN", str(fake))
    assert find_gmx() == str(fake)
