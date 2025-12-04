"""
Script para probar todas las funcionalidades del predictor de proteínas.
"""

import os
import torch
import logging
from cha_md_ba.protein_predictor import ProteinStructurePredictor

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    """Función principal para probar el predictor de proteínas."""
    try:
        # Verificar GPU
        device = "cuda" if torch.cuda.is_available() else "cpu"
        logger.info(f"Usando dispositivo: {device}")
        
        # Crear directorio de salida
        output_dir = "test_output"
        os.makedirs(output_dir, exist_ok=True)
        
        # Secuencia de prueba (proteína con dominios de unión a ARN)
        test_sequence = "MVKLAKAGKNQGDPKKMAPPPKEVEEDSEDEEMSEDEEDDSSGEEVVIPQKKGKKAAATSAKKVVVSPTKKVAVATPAKKAAVTPGKKAAATPAKKTVTPAKAVTTPGKKGATPGKALVATPGKKGAAIPAKGAKNGKNAKKEDSDEEEDDDSEEDEEDDEDEDEDEDEIEPAAMKAAAAAPASEDEDDEDDEDDEDDDDDEEDDSEEEAMETTPAKGKKAAKVVPVKAKNVAEDEDEEEDDEDEDDDDDEDDEDDDDEDDEEEEEEEEEEPVKEAPGKRKKEMAKQKAAPEAKKQKVEGTEPTTAFNLFVGNLNFNKSAPELKTGISDVFAKNDLAVVDVRIGMTRKFGYVDFESAEDLEKALELTGLKVFGNEIKLEKPKGKDSKKERDARTLLAKNLPYKVTQDELKEVFEDAAEIRLVSKDGKSKGIAYIEFKTEADAEKTFEEKQGTEIDGRSISLYYTGEKGQNQDYRGGKNSTWSGESKTLVLSNLSYSATEETLQEVFEKATFIKVPQNQNGKSKGYAFIEFASFEDAKEALNSCNKREIEGRAIRLELQGPRGSPNARSQPSKTLFVKGLSEDTTEETLKESFDGSVRARIVTDRETGSSKGFGFVDFNSEEDAKAAKEAMEDGEIDGNKVTLDWAKPKGEGGFGGRGGGRGGFGGRGGGRGGRGGFGGRGRGGFGGRGGFRGGRGGGGDHKPQGKKTKFE"
        
        # Inicializar predictor
        logger.info("Inicializando predictor de proteínas...")
        predictor = ProteinStructurePredictor()
        
        # Predecir estructura
        logger.info("Prediciendo estructura...")
        output_path = os.path.join(output_dir, "rna_binding_protein.pdb")
        predictor.predict_structure(test_sequence, output_path)
        
        # Verificar archivo de salida
        if os.path.exists(output_path):
            logger.info(f"Estructura guardada en: {output_path}")
            
            # Calcular métricas de confianza
            logger.info("Calculando métricas de confianza...")
            confidence_metrics = predictor._calculate_confidence(output_path)
            logger.info(f"Métricas de confianza: {confidence_metrics}")
            
            # Preparar para MD
            logger.info("Preparando estructura para dinámica molecular...")
            md_path = os.path.join(output_dir, "rna_binding_protein_md.pdb")
            predictor.prepare_for_md(output_path, md_path)
            
            if os.path.exists(md_path):
                logger.info(f"Estructura preparada para MD guardada en: {md_path}")
            else:
                logger.error("Error al preparar estructura para MD")
        else:
            logger.error("Error al generar estructura")
            
    except Exception as e:
        logger.error(f"Error durante la prueba: {str(e)}")
        raise

if __name__ == "__main__":
    main() 