Despliegue
=========

Proceso de Despliegue
------------------

1. Actualizar versión en setup.py
2. Actualizar CHANGELOG.md
3. Crear tag de versión
4. Construir y subir a PyPI

.. code-block:: bash

   # Construir paquete
   python setup.py sdist bdist_wheel
   
   # Subir a PyPI
   twine upload dist/* 