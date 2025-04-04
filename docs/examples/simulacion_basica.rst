Simulación Básica
==============

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