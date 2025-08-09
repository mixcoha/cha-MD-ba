"""
Módulo para análisis de trayectorias de dinámica molecular
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Tuple, Dict
from MDAnalysis import Universe
from rich.console import Console
from rich.progress import Progress
from .preprocess import TrajectoryPreprocessor

console = Console()

class MDTrajectoryAnalyzer:
    """Clase para analizar trayectorias de dinámica molecular"""
    
    def __init__(self, trajectory_path: str, topology_path: str):
        """
        Inicializa el analizador de trayectorias
        
        Args:
            trajectory_path: Ruta al archivo de trayectoria (.xtc)
            topology_path: Ruta al archivo de topología (.tpr)
        """
        self.trajectory_path = Path(trajectory_path)
        self.topology_path = Path(topology_path)
        self.universe = Universe(str(topology_path), str(trajectory_path))
        
    def preprocess_trajectory(self, output_dir: str, selection: str = "protein") -> Dict[str, Path]:
        """
        Preprocesa la trayectoria antes del análisis
        
        Args:
            output_dir: Directorio para guardar los archivos procesados
            selection: Selección de átomos para el procesamiento
            
        Returns:
            Diccionario con rutas a los archivos procesados
        """
        preprocessor = TrajectoryPreprocessor(str(self.trajectory_path), str(self.topology_path))
        return preprocessor.process_trajectory(output_dir, selection)
        
    def calculate_rmsd(self, selection: str = "protein", ref_frame: int = 0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calcula el RMSD de una selección respecto a un frame de referencia
        
        Args:
            selection: Selección de átomos para el cálculo
            ref_frame: Frame de referencia para el RMSD
            
        Returns:
            Tuple con tiempos y valores de RMSD
        """
        with Progress() as progress:
            task = progress.add_task("[cyan]Calculando RMSD...", total=len(self.universe.trajectory))
            
            rmsd_values = []
            times = []
            
            ref = self.universe.select_atoms(selection)
            ref_pos = ref.positions
            
            for ts in self.universe.trajectory:
                mobile = self.universe.select_atoms(selection)
                rmsd = np.sqrt(np.mean(np.sum((mobile.positions - ref_pos) ** 2, axis=1)))
                rmsd_values.append(rmsd)
                times.append(ts.time)
                progress.advance(task)
                
        return np.array(times), np.array(rmsd_values)
    
    def calculate_rog(self, selection: str = "protein") -> Tuple[np.ndarray, np.ndarray]:
        """
        Calcula el radio de giro de una selección
        
        Args:
            selection: Selección de átomos para el cálculo
            
        Returns:
            Tuple con tiempos y valores de radio de giro
        """
        with Progress() as progress:
            task = progress.add_task("[cyan]Calculando radio de giro...", total=len(self.universe.trajectory))
            
            rog_values = []
            times = []
            
            for ts in self.universe.trajectory:
                group = self.universe.select_atoms(selection)
                rog = group.radius_of_gyration()
                rog_values.append(rog)
                times.append(ts.time)
                progress.advance(task)
                
        return np.array(times), np.array(rog_values)
    
    def plot_rmsd(self, output_path: Optional[str] = None):
        """
        Genera un gráfico de RMSD vs tiempo
        
        Args:
            output_path: Ruta para guardar el gráfico
        """
        times, rmsd = self.calculate_rmsd()
        
        plt.figure(figsize=(10, 6))
        plt.plot(times/1000, rmsd, 'k-', label='RMSD')
        plt.xlabel('Tiempo (ns)')
        plt.ylabel('RMSD (nm)')
        plt.title('RMSD vs Tiempo')
        plt.grid(True)
        plt.legend()
        
        if output_path:
            plt.savefig(output_path)
        else:
            plt.show()
            
    def plot_rog(self, output_path: Optional[str] = None):
        """
        Genera un gráfico de radio de giro vs tiempo
        
        Args:
            output_path: Ruta para guardar el gráfico
        """
        times, rog = self.calculate_rog()
        
        plt.figure(figsize=(10, 6))
        plt.plot(times/1000, rog, 'k-', label='Radio de giro')
        plt.xlabel('Tiempo (ns)')
        plt.ylabel('Radio de giro (nm)')
        plt.title('Radio de giro vs Tiempo')
        plt.grid(True)
        plt.legend()
        
        if output_path:
            plt.savefig(output_path)
        else:
            plt.show()
            
    def analyze_hbonds(self, selection1: str, selection2: str) -> dict:
        """
        Analiza puentes de hidrógeno entre dos selecciones
        
        Args:
            selection1: Primera selección de átomos
            selection2: Segunda selección de átomos
            
        Returns:
            Diccionario con estadísticas de puentes de hidrógeno
        """
        # Implementar análisis de puentes de hidrógeno
        pass
        
    def cluster_analysis(self, selection: str = "protein", cutoff: float = 0.2) -> dict:
        """
        Realiza análisis de clustering de la trayectoria
        
        Args:
            selection: Selección de átomos para el clustering
            cutoff: Distancia de corte para el clustering
            
        Returns:
            Diccionario con resultados del clustering
        """
        # Implementar análisis de clustering
        pass 