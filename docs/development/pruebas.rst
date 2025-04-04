Pruebas
======

Ejecutar Pruebas
--------------

.. code-block:: bash

   # Instalar dependencias de desarrollo
   pip install -e ".[dev]"
   
   # Ejecutar pruebas
   pytest
   
   # Ejecutar con cobertura
   pytest --cov=cha_md_ba

Estructura de Pruebas
------------------

::

   tests/
   ├── unit/              # Pruebas unitarias
   ├── integration/       # Pruebas de integración
   └── fixtures/         # Datos de prueba 