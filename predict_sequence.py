"""
Script para predecir la estructura de una secuencia de proteína ingresada por el usuario.
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
    os.makedirs("output", exist_ok=True)
    
    try:
        # Obtener secuencia del usuario
        print("\n=== Predicción de Estructura de Proteínas ===")
        print("Ingresa la secuencia de aminoácidos (formato de una letra)")
        print("Ejemplo: FVNQHLCGSHLVEALYLVCGERGFFYTPKT")
        sequence = input("\nSecuencia: ").strip().upper()
        
        # Validar secuencia
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        if not all(aa in valid_aa for aa in sequence):
            print("Error: La secuencia contiene caracteres no válidos.")
            print("Solo se permiten los siguientes aminoácidos:")
            print("A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y")
            return
        
        # Inicializar predictor
        logger.info("Inicializando predictor de estructura...")
        predictor = ProteinStructurePredictor(device="cpu")
        
        # Predecir estructura
        logger.info(f"Prediciendo estructura para secuencia de longitud {len(sequence)}...")
        output_path = f"output/{sequence[:10]}_predicted.pdb"
        confidence = predictor.predict_structure(
            sequence=sequence,
            output_path=output_path,
            visualize=True
        )
        
        # Mostrar resultados
        print("\n=== Resultados ===")
        print(f"Estructura guardada en: {output_path}")
        print(f"Imagen guardada en: {output_path.replace('.pdb', '_view.png')}")
        print("\nMétricas de confianza:")
        print(f"- pLDDT: {confidence['plddt']:.2f}")
        print(f"- PAE: {confidence['pae']:.2f}")
        print(f"- Confianza general: {confidence['confidence']:.2%}")
        
        # Preparar para MD
        logger.info("Preparando estructura para MD...")
        md_output = output_path.replace('.pdb', '_md.pdb')
        predictor.prepare_for_md(
            input_pdb=output_path,
            output_pdb=md_output,
            visualize=True
        )
        
        print(f"\nEstructura preparada para MD guardada en: {md_output}")
        print(f"Imagen de MD guardada en: {md_output.replace('.pdb', '_view.png')}")
        
    except Exception as e:
        logger.error(f"Error durante la predicción: {str(e)}")
        raise

if __name__ == "__main__":
    main() 