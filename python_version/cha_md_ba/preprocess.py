"""
Módulo para preprocesamiento de trayectorias de dinámica molecular
"""

import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple
from rich.console import Console
from rich.progress import Progress
from MDAnalysis import Universe
from MDAnalysis.analysis import align
import numpy as np

console = Console()

class TrajectoryPreprocessor:
    """Clase para preprocesar trayectorias de dinámica molecular"""
    
    def __init__(self, trajectory_path: str, topology_path: str):
        """
        Inicializa el preprocesador de trayectorias
        
        Args:
            trajectory_path: Ruta al archivo de trayectoria (.xtc)
            topology_path: Ruta al archivo de topología (.tpr)
        """
        self.trajectory_path = Path(trajectory_path)
        self.topology_path = Path(topology_path)
        self.universe = Universe(str(topology_path), str(trajectory_path))
        
    def center_and_align(self, output_path: str, selection: str = "protein", 
                        ref_frame: int = 0, center_selection: str = "protein") -> Path:
        """
        Centra y alinea la trayectoria
        
        Args:
            output_path: Ruta para guardar la trayectoria procesada
            selection: Selección de átomos para alineación
            ref_frame: Frame de referencia para alineación
            center_selection: Selección de átomos para centrado
            
        Returns:
            Ruta al archivo de trayectoria procesada
        """
        output_path = Path(output_path)
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Centrando y alineando trayectoria...", total=1)
            
            # Seleccionar átomos para alineación y centrado
            align_atoms = self.universe.select_atoms(selection)
            center_atoms = self.universe.select_atoms(center_selection)
            
            # Crear frame de referencia
            ref = self.universe.trajectory[ref_frame]
            ref_pos = align_atoms.positions
            
            # Preparar archivo de salida
            with self.universe.trajectory.Writer(str(output_path), align_atoms.n_atoms) as W:
                for ts in self.universe.trajectory:
                    # Alinear
                    mobile_pos = align_atoms.positions
                    R, COM = align.rotation_matrix(mobile_pos, ref_pos)
                    align_atoms.positions = np.dot(align_atoms.positions - COM, R.T) + COM
                    
                    # Centrar
                    center_com = center_atoms.center_of_mass()
                    self.universe.atoms.positions -= center_com
                    
                    # Escribir frame
                    W.write(align_atoms)
                    
            progress.advance(task)
            
        return output_path
        
    def remove_rotations_translations(self, output_path: str, selection: str = "protein") -> Path:
        """
        Elimina rotaciones y traslaciones globales de la trayectoria
        
        Args:
            output_path: Ruta para guardar la trayectoria procesada
            selection: Selección de átomos para el procesamiento
            
        Returns:
            Ruta al archivo de trayectoria procesada
        """
        output_path = Path(output_path)
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Eliminando rotaciones y traslaciones...", total=1)
            
            # Seleccionar átomos
            atoms = self.universe.select_atoms(selection)
            
            # Crear archivo de salida
            with self.universe.trajectory.Writer(str(output_path), atoms.n_atoms) as W:
                for ts in self.universe.trajectory:
                    # Calcular centro de masa
                    com = atoms.center_of_mass()
                    
                    # Eliminar traslación
                    atoms.positions -= com
                    
                    # Calcular matriz de inercia
                    positions = atoms.positions
                    masses = atoms.masses
                    inertia = np.zeros((3, 3))
                    for i in range(3):
                        for j in range(3):
                            inertia[i,j] = np.sum(masses * positions[:,i] * positions[:,j])
                    
                    # Diagonalizar matriz de inercia
                    eigenvals, eigenvecs = np.linalg.eigh(inertia)
                    
                    # Rotar al sistema de ejes principales
                    atoms.positions = np.dot(atoms.positions, eigenvecs)
                    
                    # Escribir frame
                    W.write(atoms)
                    
            progress.advance(task)
            
        return output_path
        
    def process_trajectory(self, output_dir: str, selection: str = "protein") -> dict:
        """
        Procesa completamente una trayectoria
        
        Args:
            output_dir: Directorio para guardar los archivos procesados
            selection: Selección de átomos para el procesamiento
            
        Returns:
            Diccionario con rutas a los archivos procesados
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 1. Centrar y alinear
        aligned_path = output_dir / "aligned.xtc"
        self.center_and_align(str(aligned_path), selection)
        
        # 2. Eliminar rotaciones y traslaciones
        processed_path = output_dir / "processed.xtc"
        self.remove_rotations_translations(str(processed_path), selection)
        
        return {
            "aligned": aligned_path,
            "processed": processed_path
        } 