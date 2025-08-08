Preparación del Sistema
====================

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