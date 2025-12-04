"""
Módulo para visualización de estructuras de proteínas.
"""

import os
import logging
import nglview as nv
from typing import Optional

logger = logging.getLogger(__name__)

class StructureViewer:
    """Clase para visualizar estructuras de proteínas."""
    
    def __init__(self):
        """Inicializa el visualizador."""
        self.current_view = None
        
    def view_structure(self, pdb_path: str) -> nv.NGLWidget:
        """
        Visualiza una estructura PDB de forma interactiva.
        
        Args:
            pdb_path: Ruta al archivo PDB
            
        Returns:
            Widget de NGLView para visualización interactiva
        """
        try:
            if not os.path.exists(pdb_path):
                raise FileNotFoundError(f"No se encontró el archivo PDB: {pdb_path}")
                
            # Crear vista
            view = nv.show_structure_file(pdb_path)
            
            # Configurar representación
            view.add_representation('cartoon', selection='protein', color='residueindex')
            view.add_representation('ball+stick', selection='protein')
            
            self.current_view = view
            return view
            
        except Exception as e:
            logger.error(f"Error al visualizar estructura: {str(e)}")
            raise
            
    def save_image(self, output_path: str) -> None:
        """
        Guarda una imagen de la estructura actual.
        
        Args:
            output_path: Ruta donde guardar la imagen
        """
        try:
            if self.current_view is None:
                raise ValueError("No hay una estructura visualizada para guardar")
                
            # Guardar imagen
            self.current_view.render_image(filename=output_path)
            logger.info(f"Imagen guardada en: {output_path}")
            
        except Exception as e:
            logger.error(f"Error al guardar imagen: {str(e)}")
            raise 