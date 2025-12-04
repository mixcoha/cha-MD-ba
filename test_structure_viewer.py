"""
Script para probar el visualizador de estructuras PDB.
"""

import os
import logging
from cha_md_ba.visualization.structure_viewer import StructureViewer

# Configurar logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    # Crear directorio de salida
    os.makedirs("test_output", exist_ok=True)
    
    # Ruta al archivo PDB
    pdb_path = "test_output/insulin.pdb"
    
    try:
        # Inicializar visualizador
        logger.info("Inicializando visualizador de estructuras...")
        viewer = StructureViewer()
        
        # Visualizar estructura
        logger.info("Visualizando estructura...")
        view = viewer.view_structure(pdb_path)
        
        # Guardar imagen
        logger.info("Guardando imagen...")
        viewer.save_image("test_output/insulin_structure.png")
        
        logger.info("Visualización completada exitosamente")
        
    except Exception as e:
        logger.error(f"Error durante la prueba: {str(e)}")
        raise

if __name__ == "__main__":
    main() 