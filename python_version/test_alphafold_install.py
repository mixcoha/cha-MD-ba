"""
Script de prueba para verificar la instalación de AlphaFold3.
"""

import os
import torch
import logging
from pathlib import Path
from cha_md_ba.alphafold import AlphaFoldPredictor

# Configurar logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    # Verificar disponibilidad de GPU
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Usando dispositivo: {device}")
    
    # Crear directorio de prueba
    test_dir = Path("test_output")
    test_dir.mkdir(exist_ok=True)
    
    # Secuencia de prueba (fragmento de una proteína pequeña)
    test_sequence = "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"
    
    try:
        # Inicializar predictor
        predictor = AlphaFoldPredictor(
            device=device,
            model_preset="monomer",
            num_ensemble=1
        )
        logger.info("Predictor inicializado correctamente")
        
        # Intentar una predicción
        result = predictor.predict_structure(
            sequence=test_sequence,
            output_dir=str(test_dir)
        )
        
        logger.info(f"Predicción completada. Resultados guardados en: {result['pdb_path']}")
        logger.info(f"Confianza de la predicción: {result['confidence']}")
        
    except Exception as e:
        logger.error(f"Error durante la prueba: {e}")
        raise

if __name__ == "__main__":
    main() 