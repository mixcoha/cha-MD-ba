Análisis Avanzado
==============

Análisis de Estructura Secundaria
------------------------------

.. code-block:: python

   from cha_md_ba.analysis import SecondaryStructure
   
   analyzer = SecondaryStructure()
   analyzer.analyze("trajectory.xtc")
   analyzer.plot()

Análisis de Energía
----------------

.. code-block:: python

   from cha_md_ba.analysis import EnergyAnalysis
   
   analyzer = EnergyAnalysis()
   analyzer.analyze("ener.edr")
   analyzer.plot() 