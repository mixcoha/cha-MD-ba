Ejemplos de Uso
=============

.. toctree::
   :maxdepth: 2
   :caption: Contenido:

   preparacion_sistema
   simulacion_basica
   analisis_trayectorias
   visualizacion_resultados

Preparación del Sistema
---------------------

Este ejemplo muestra cómo preparar un sistema para simulación de dinámica molecular:

.. code-block:: python

   from cha_md_ba.core import SystemPreparation
   
   # Crear instancia del preparador
   prep = SystemPreparation()
   
   # Cargar estructura PDB
   prep.load_pdb("proteina.pdb")
   
   # Añadir caja de agua
   prep.add_water_box()
   
   # Generar topología
   prep.generate_topology()
   
   # Guardar archivos
   prep.save_files()

Simulación Básica
---------------

Ejemplo de cómo ejecutar una simulación básica:

.. code-block:: python

   from cha_md_ba.core import Simulation
   
   # Crear instancia de simulación
   sim = Simulation()
   
   # Configurar parámetros
   sim.set_parameters(
       temperature=300,
       pressure=1,
       time_step=0.002,
       n_steps=1000000
   )
   
   # Ejecutar simulación
   sim.run()

Análisis de Trayectorias
----------------------

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

Visualización de Resultados
-------------------------

Ejemplo de visualización de resultados:

.. code-block:: python

   from cha_md_ba.visualization import ResultsVisualizer
   
   # Crear visualizador
   viz = ResultsVisualizer()
   
   # Cargar datos
   viz.load_data("results.pkl")
   
   # Generar gráficos
   viz.plot_energy()
   viz.plot_temperature()
   viz.plot_pressure()
   
   # Exportar reporte
   viz.export_report("report.pdf") 