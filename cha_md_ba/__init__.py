"""
CHA-MD-BA: Herramientas avanzadas para automatizar simulaciones de dinámica molecular
"""

# Importar módulos principales (versión avanzada)
from .prepare import MDSystemPreparator
from .minimize import EnergyMinimizer
from .nvt import NVTEquilibrator
from .npt import NPTEquilibrator
from .analysis import MDTrajectoryAnalyzer
from .preprocess import TrajectoryPreprocessor
from .cli import main as cli_main

# Importar módulos básicos (compatibilidad)
from .core.system import SystemPreparation
from .analysis.trajectory import TrajectoryAnalysis
from .visualization.plotter import ResultsVisualizer

__version__ = "0.2.0"
__all__ = [
    # Módulos avanzados
    'MDSystemPreparator', 'EnergyMinimizer', 'NVTEquilibrator', 
    'NPTEquilibrator', 'MDTrajectoryAnalyzer', 'TrajectoryPreprocessor',
    'cli_main',
    # Módulos básicos (compatibilidad)
    'SystemPreparation', 'TrajectoryAnalysis', 'ResultsVisualizer'
] 