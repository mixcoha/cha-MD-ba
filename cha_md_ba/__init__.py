"""CHA-MD-BA: Herramienta para análisis de dinámica molecular."""

from .core.system import SystemPreparation
from .analysis.trajectory import TrajectoryAnalysis
from .visualization.plotter import ResultsVisualizer

__version__ = '0.1.0'
__all__ = ['SystemPreparation', 'TrajectoryAnalysis', 'ResultsVisualizer']
