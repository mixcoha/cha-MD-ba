"""
Script de prueba para verificar la instalación de ESMFold.
"""

import os
import torch
import logging
from pathlib import Path
from cha_md_ba.protein_predictor import ProteinStructurePredictor

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
    
    # Secuencia de prueba (insulina humana)
    test_sequence = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    
    try:
        # Inicializar el predictor
        predictor = ProteinStructurePredictor()
        
        # Predecir estructura
        output_path = test_dir / "test_structure.pdb"
        result = predictor.predict_structure(test_sequence, str(output_path))
        
        logger.info(f"Predicción completada exitosamente")
        logger.info(f"Archivo de estructura guardado en: {output_path}")
        logger.info(f"Métricas de confianza: {result}")
        
    except Exception as e:
        logger.error(f"Error durante la predicción: {str(e)}")
        raise

if __name__ == "__main__":
    main() 