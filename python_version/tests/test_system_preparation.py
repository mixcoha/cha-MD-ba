import os
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, Mock
import pytest
from cha_md_ba.prepare import MDSystemPreparator
from cha_md_ba.minimize import EnergyMinimizer

@pytest.fixture
def temp_dir():
    """Fixture que crea un directorio temporal y lo limpia después"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)

@pytest.fixture
def mock_gmx():
    """Fixture que mockea los comandos de GROMACS"""
    with patch("subprocess.run") as mock_run:
        # Configurar el mock para simular éxito en los comandos
        mock_run.return_value = Mock(returncode=0)
        yield mock_run

def test_system_preparation_initialization():
    """Test para la inicialización de MDSystemPreparator"""
    preparator = MDSystemPreparator(
        pdb_path="test.pdb",
        forcefield="amber99sb-ildn",
        water_model="tip3p"
    )
    
    assert preparator.pdb_path == Path("test.pdb")
    assert preparator.forcefield == "amber99sb-ildn"
    assert preparator.water_model == "tip3p"
    assert preparator.gmx == "gmx_mpi"
    assert preparator.system_name == "test"

def test_directory_structure_creation(temp_dir, mock_gmx):
    """Test para la creación de la estructura de directorios"""
    # Crear un archivo PDB de prueba
    test_pdb = temp_dir / "test.pdb"
    test_pdb.write_text("HEADER TEST PDB")
    
    preparator = MDSystemPreparator(pdb_path=str(test_pdb))
    result = preparator.prepare_system(output_dir=str(temp_dir), minimize=False)
    
    # Verificar que se crearon los directorios correctos
    assert (temp_dir / "test").exists()
    assert (temp_dir / "test" / "1_preparation").exists()
    assert (temp_dir / "test" / "2_minimization").exists()
    assert (temp_dir / "test" / "3_nvt").exists()
    assert (temp_dir / "test" / "4_npt").exists()
    
    # Verificar que el PDB se copió correctamente
    assert (temp_dir / "test" / "1_preparation" / "test.pdb").exists()

def test_minimizer_initialization():
    """Test para la inicialización de EnergyMinimizer"""
    minimizer = EnergyMinimizer(
        input_gro="test.gro",
        topol_top="topol.top",
        mdp_file="min.mdp"
    )
    
    assert minimizer.input_gro == Path("test.gro")
    assert minimizer.topol_top == Path("topol.top")
    assert minimizer.mdp_file == Path("min.mdp")
    assert minimizer.gmx == "gmx_mpi"

def test_minimizer_mdp_creation(temp_dir):
    """Test para la creación del archivo .mdp"""
    minimizer = EnergyMinimizer(
        input_gro="test.gro",
        topol_top="topol.top"
    )
    
    mdp_path = minimizer.create_mdp_file(temp_dir / "min.mdp")
    
    assert mdp_path.exists()
    content = mdp_path.read_text()
    assert "title               = Energy Minimization" in content
    assert "integrator          = steep" in content
    assert "emtol               = 1000.0" in content

def test_gmx_mpi_usage(temp_dir, mock_gmx):
    """Test para verificar que se usa gmx_mpi en los comandos"""
    # Crear archivos de prueba
    test_pdb = temp_dir / "test.pdb"
    test_pdb.write_text("HEADER TEST PDB")
    
    preparator = MDSystemPreparator(pdb_path=str(test_pdb))
    preparator.prepare_system(output_dir=str(temp_dir), minimize=False)
    
    # Verificar que se usó gmx_mpi en los comandos
    for call in mock_gmx.call_args_list:
        cmd = call[0][0]  # El primer argumento de subprocess.run es la lista de comandos
        assert cmd[0] == "gmx_mpi" 