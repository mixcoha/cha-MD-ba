"""
Script para probar el predictor de estructura de proteínas.
"""

import os
import logging
from cha_md_ba import ProteinStructurePredictor

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    # Crear directorio de salida
    os.makedirs("test_output", exist_ok=True)
    
    # Secuencia de prueba (insulina humana)
    sequence = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    
    try:
        # Inicializar predictor
        logger.info("Inicializando predictor de estructura...")
        predictor = ProteinStructurePredictor(device="cpu")
        
        # Predecir estructura
        logger.info("Prediciendo estructura...")
        confidence = predictor.predict_structure(
            sequence=sequence,
            output_path="test_output/insulin.pdb",
            visualize=True  # Activar visualización
        )
        
        logger.info(f"Predicción completada con métricas de confianza: {confidence}")
        
        # Preparar para MD
        logger.info("Preparando estructura para MD...")
        predictor.prepare_for_md(
            input_pdb="test_output/insulin.pdb",
            output_pdb="test_output/insulin_md.pdb",
            visualize=True  # Activar visualización
        )
        
        logger.info("Prueba completada exitosamente")
        
    except Exception as e:
        logger.error(f"Error durante la prueba: {str(e)}")
        raise

if __name__ == "__main__":
    main() 