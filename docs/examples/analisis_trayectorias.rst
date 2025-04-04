Análisis de Trayectorias
=====================

Ejemplo de análisis de trayectorias:

.. code-block:: python

   from cha_md_ba.analysis import TrajectoryAnalysis
   
   # Crear analizador
   analyzer = TrajectoryAnalysis()
   
   # Cargar trayectoria
   analyzer.load_trajectory("trajectory.xtc")
   
   # Calcular RMSD
   rmsd = analyzer.calculate_rmsd()
   
   # Calcular RMSF
   rmsf = analyzer.calculate_rmsf()
   
   # Generar gráficos
   analyzer.plot_rmsd()
   analyzer.plot_rmsf() 