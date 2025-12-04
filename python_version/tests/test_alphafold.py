"""
Pruebas para el módulo de integración de AlphaFold.
"""

import os
import pytest
from pathlib import Path
from cha_md_ba.alphafold import AlphaFoldPredictor

@pytest.fixture
def test_sequence():
    return "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"

@pytest.fixture
def test_output_dir(tmp_path):
    return tmp_path

@pytest.fixture
def predictor():
    return AlphaFoldPredictor(
        device="cpu",
        model_preset="monomer",
        num_ensemble=1,
        max_seq_len=2048
    )

def test_predictor_initialization(predictor):
    """Prueba la inicialización del predictor."""
    assert predictor.device == "cpu"
    assert predictor.model_path is None
    assert predictor.model_preset == "monomer"
    assert predictor.num_ensemble == 1
    assert predictor.max_seq_len == 2048
    assert predictor.model is not None
    assert predictor.model_config is not None

def test_predictor_initialization_with_custom_model(tmp_path):
    """Prueba la inicialización del predictor con un modelo personalizado."""
    model_path = tmp_path / "model"
    model_path.mkdir()
    
    predictor = AlphaFoldPredictor(
        model_path=str(model_path),
        device="cpu",
        model_preset="monomer_ptm"
    )
    
    assert predictor.model_path == str(model_path)
    assert predictor.model_preset == "monomer_ptm"
    assert predictor.model is not None
    assert predictor.model_config is not None

def test_predictor_initialization_invalid_preset():
    """Prueba la inicialización del predictor con un preset inválido."""
    with pytest.raises(ValueError):
        AlphaFoldPredictor(model_preset="invalid_preset")

def test_predict_structure(predictor, test_sequence, test_output_dir):
    """Prueba la predicción de estructura."""
    result = predictor.predict_structure(
        sequence=test_sequence,
        output_dir=str(test_output_dir)
    )
    
    assert isinstance(result, dict)
    assert "pdb_path" in result
    assert "confidence" in result
    assert "metrics" in result
    assert os.path.exists(result["pdb_path"])

def test_prepare_for_md(predictor, test_output_dir):
    """Prueba la preparación para simulación MD."""
    # Crear un archivo PDB de prueba
    test_pdb = test_output_dir / "test.pdb"
    test_pdb.write_text("ATOM      1  N   ALA A   1      27.268  24.339   4.298  1.00 10.00\n")
    
    result = predictor.prepare_for_md(
        pdb_path=str(test_pdb),
        output_dir=str(test_output_dir)
    )
    
    assert isinstance(result, dict)
    assert "prepared_pdb" in result
    assert "topology" in result
    assert "parameters" in result
    assert os.path.exists(result["prepared_pdb"]) 