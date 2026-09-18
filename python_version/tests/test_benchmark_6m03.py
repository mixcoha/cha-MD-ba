"""Pruebas del benchmark 6M03 y de las condiciones NaCl 0.15 M / 310 K."""

from unittest.mock import Mock

from cha_md_ba.benchmark import (
    DEFAULT_DATA_DIR,
    DEFAULT_OUTPUT_DIR,
    Benchmark6M03Config,
    config_for_model,
    parse_system_composition,
    resolve_model_keys,
    resume_stages,
    write_protocol_mdps,
)
from cha_md_ba.gmx_utils import detect_gpu_ids, find_gmx
from cha_md_ba.npt import NPTEquilibrator
from cha_md_ba.nvt import NVTEquilibrator
from cha_md_ba.pdb_utils import clean_protein_pdb, mutate_to_alanine, parse_mutation
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
    assert config.ion_concentration == 0.15
    assert config.pressure == 1.0
    assert config.forcefield == "amber99sb-ildn"
    assert config.water_model == "tip3p"


def test_benchmark_defaults_to_local_work_dir():
    assert DEFAULT_OUTPUT_DIR == "work"
    assert DEFAULT_DATA_DIR == "work/data"


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
        ion_concentration=0.15,
        gmx="gmx_mpi",
    )
    preparator.prepare_system(output_dir=str(tmp_path), minimize=False, ions=True)

    genion_cmds = [cmd for cmd in captured if len(cmd) > 1 and cmd[1] == "genion"]
    assert genion_cmds, "No se invocó genion"
    genion = genion_cmds[0]
    assert "-conc" in genion
    assert "0.15" in genion
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
    assert "0.15 M NaCl" in md_text
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


def test_resume_stages_from_empty_dir(tmp_path):
    stages = resume_stages(tmp_path, Benchmark6M03Config())
    assert stages[0] == "download"
    assert "nvt" in stages
    assert "npt" in stages


def test_resume_stages_continues_nvt_after_minimize(tmp_path):
    base = tmp_path / "6M03"
    prep = base / "1_preparation"
    prep.mkdir(parents=True)
    (prep / "topol.top").write_text("[ molecules ]\nSOL 1\n")
    (prep / "6M03_ions.gro").write_text("dummy\n")
    min_dir = base / "2_minimization"
    min_dir.mkdir()
    (min_dir / "minimized.gro").write_text("dummy\n")
    nvt_1000 = base / "3_nvt" / "posre_constante" / "1000"
    nvt_1000.mkdir(parents=True)
    (nvt_1000 / "nvt.gro").write_text("dummy\n")
    (nvt_1000 / "nvt.edr").write_bytes(b"")

    stages = resume_stages(tmp_path, Benchmark6M03Config())
    assert stages == ["nvt", "npt"]


def test_resume_stages_npt_when_nvt_done(tmp_path):
    base = tmp_path / "6M03"
    prep = base / "1_preparation"
    prep.mkdir(parents=True)
    (prep / "topol.top").write_text("[ molecules ]\nSOL 1\n")
    (prep / "6M03_ions.gro").write_text("dummy\n")
    min_dir = base / "2_minimization"
    min_dir.mkdir()
    (min_dir / "minimized.gro").write_text("dummy\n")
    for fc in (1000, 800, 600, 400, 200):
        fc_dir = base / "3_nvt" / "posre_constante" / str(fc)
        fc_dir.mkdir(parents=True)
        (fc_dir / "nvt.gro").write_text("dummy\n")
        (fc_dir / "nvt.edr").write_bytes(b"")

    stages = resume_stages(tmp_path, Benchmark6M03Config())
    assert stages == ["npt"]


def test_nvt_skips_completed_force_constant(tmp_path, monkeypatch):
    gro = tmp_path / "min.gro"
    top = tmp_path / "topol.top"
    gro.write_text("dummy\n")
    top.write_text("[ molecules ]\n")
    nvt_dir = tmp_path / "3_nvt"
    done = nvt_dir / "posre_constante" / "1000"
    done.mkdir(parents=True)
    (done / "nvt.gro").write_text("done\n")
    (done / "nvt.edr").write_bytes(b"")
    (done / "tmp.gro").write_text("centered\n")

    captured = []

    def fake_run(cmd, **kwargs):
        captured.append(list(cmd))
        return Mock(returncode=0)

    monkeypatch.setattr("cha_md_ba.nvt.subprocess.run", fake_run)
    nvt = NVTEquilibrator(str(gro), str(top), temperature=310, gmx="gmx")
    results = nvt.equilibrate(str(nvt_dir), force_constants=[1000])
    assert results[1000]["gro"].endswith("tmp.gro")
    assert captured == []


def test_detect_gpu_ids_none_and_explicit():
    assert detect_gpu_ids("none") is None
    assert detect_gpu_ids("cpu") is None
    assert detect_gpu_ids("0") == "0"


HIS41_CYS145_PDB = """\
HEADER    TEST MPRO DYAD
ATOM      1  N   SER A   1      11.000  12.000  13.000  1.00 10.00           N
ATOM      2  CA  SER A   1      12.000  12.000  13.000  1.00 10.00           C
ATOM     10  N   HIS A  41      20.000  21.000  22.000  1.00 10.00           N
ATOM     11  CA  HIS A  41      21.000  21.000  22.000  1.00 10.00           C
ATOM     12  C   HIS A  41      22.000  21.000  22.000  1.00 10.00           C
ATOM     13  O   HIS A  41      23.000  21.000  22.000  1.00 10.00           O
ATOM     14  CB  HIS A  41      21.000  22.000  22.000  1.00 10.00           C
ATOM     15  CG  HIS A  41      21.000  23.000  22.000  1.00 10.00           C
ATOM     16  ND1 HIS A  41      20.000  24.000  22.000  1.00 10.00           N
ATOM     17  CD2 HIS A  41      22.000  24.000  22.000  1.00 10.00           C
ATOM     18  CE1 HIS A  41      21.000  25.000  22.000  1.00 10.00           C
ATOM     19  NE2 HIS A  41      22.000  25.000  22.000  1.00 10.00           N
ATOM     30  N   CYS A 145      30.000  31.000  32.000  1.00 10.00           N
ATOM     31  CA  CYS A 145      31.000  31.000  32.000  1.00 10.00           C
ATOM     32  C   CYS A 145      32.000  31.000  32.000  1.00 10.00           C
ATOM     33  O   CYS A 145      33.000  31.000  32.000  1.00 10.00           O
ATOM     34  CB  CYS A 145      31.000  32.000  32.000  1.00 10.00           C
ATOM     35  SG  CYS A 145      31.000  33.000  32.000  1.00 10.00           S
END
"""


def test_parse_mutation_h41a_and_c145a():
    h41a = parse_mutation("H41A")
    assert h41a.from_aa == "H"
    assert h41a.resseq == 41
    assert h41a.to_aa == "A"
    assert h41a.chain == "A"
    assert h41a.label == "H41A"
    c145a = parse_mutation("C145A")
    assert c145a.resseq == 145
    assert c145a.from_aa == "C"


def test_mutate_model1_h41a_drops_imidazole(tmp_path):
    raw = tmp_path / "clean.pdb"
    out = tmp_path / "H41A.pdb"
    raw.write_text(HIS41_CYS145_PDB)
    stats = mutate_to_alanine(str(raw), str(out), ["H41A"])
    text = out.read_text()
    assert stats["atoms_dropped"] == 5  # CG, ND1, CD2, CE1, NE2
    assert "H41A" in text
    assert " ND1 " not in text
    assert " SG  CYS A 145" in text
    assert " CB  ALA A  41" in text or " CB  ALA A 41" in text
    assert "HIS A  41" not in text
    assert "CYS A 145" in text


def test_mutate_model2_c145a_drops_sg(tmp_path):
    raw = tmp_path / "clean.pdb"
    out = tmp_path / "C145A.pdb"
    raw.write_text(HIS41_CYS145_PDB)
    stats = mutate_to_alanine(str(raw), str(out), ["C145A"])
    text = out.read_text()
    assert stats["atoms_dropped"] == 1
    assert " SG " not in text
    assert " ND1 HIS A  41" in text
    assert " CB  ALA A 145" in text or " CB  ALA A145" in text


def test_mutate_model3_double(tmp_path):
    raw = tmp_path / "clean.pdb"
    out = tmp_path / "double.pdb"
    raw.write_text(HIS41_CYS145_PDB)
    stats = mutate_to_alanine(str(raw), str(out), ["H41A", "C145A"])
    text = out.read_text()
    assert stats["atoms_dropped"] == 6
    assert "HIS" not in text
    assert "CYS" not in text
    assert text.count("ALA") >= 2


def test_mutate_hie_after_clean(tmp_path):
    raw = tmp_path / "raw.pdb"
    clean = tmp_path / "clean.pdb"
    mutated = tmp_path / "H41A.pdb"
    raw.write_text(HIS41_CYS145_PDB)
    clean_protein_pdb(str(raw), str(clean))
    assert "HIE" in clean.read_text()
    mutate_to_alanine(str(clean), str(mutated), ["H41A"])
    text = mutated.read_text()
    assert "HIE A  41" not in text
    assert "ALA" in text


def test_mutate_wrong_residue_raises(tmp_path):
    raw = tmp_path / "clean.pdb"
    raw.write_text(HIS41_CYS145_PDB)
    try:
        mutate_to_alanine(str(raw), str(tmp_path / "x.pdb"), ["A41A"])
        assert False, "debía fallar"
    except ValueError as exc:
        assert "H41A" not in str(exc) or "ALA" in str(exc) or "HIS" in str(exc)


def test_catalytic_models_and_resume_includes_mutate(tmp_path):
    assert resolve_model_keys("all") == ["1", "2", "3"]
    m1 = config_for_model("1")
    m2 = config_for_model("2")
    m3 = config_for_model("3")
    assert m1.mutations == ("H41A",)
    assert m1.run_id == "6M03_H41A"
    assert m2.mutations == ("C145A",)
    assert m2.run_id == "6M03_C145A"
    assert m3.mutations == ("H41A", "C145A")
    assert m3.run_id == "6M03_H41A_C145A"
    stages = resume_stages(tmp_path, m1)
    assert "mutate" in stages
    assert stages[0] == "download"
    wt_stages = resume_stages(tmp_path, Benchmark6M03Config())
    assert "mutate" not in wt_stages

