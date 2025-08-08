Visualización de Resultados
========================

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