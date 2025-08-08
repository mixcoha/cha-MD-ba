Guía de Desarrollo
===============

.. toctree::
   :maxdepth: 2
   :caption: Contenido:

   estructura_proyecto
   contribucion
   pruebas
   despliegue

Estructura del Proyecto
--------------------

Organización de Directorios
~~~~~~~~~~~~~~~~~~~~~~~

::

   cha-MD-ba/
   ├── cha_md_ba/           # Código fuente principal
   │   ├── cli/            # Interfaz de línea de comandos
   │   ├── core/           # Funcionalidad principal
   │   ├── analysis/       # Análisis de datos
   │   └── visualization/  # Visualización de resultados
   ├── tests/              # Pruebas unitarias y de integración
   ├── docs/               # Documentación
   ├── examples/           # Ejemplos de uso
   └── scripts/            # Scripts de utilidad

Estándares de Código
~~~~~~~~~~~~~~~~~

* Seguir PEP 8 para estilo de código
* Documentar todas las funciones y clases
* Incluir pruebas unitarias
* Mantener la cobertura de código > 80%

Ejemplo de Documentación
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def prepare_system(pdb_file: str, output_dir: str) -> None:
       """Prepara un sistema para simulación de dinámica molecular.
       
       Args:
           pdb_file (str): Ruta al archivo PDB de entrada
           output_dir (str): Directorio de salida para los archivos generados
           
       Raises:
           FileNotFoundError: Si el archivo PDB no existe
           ValueError: Si el archivo PDB no es válido
       """
       pass

Contribución
----------

Proceso de Contribución
~~~~~~~~~~~~~~~~~~~~

1. Fork el repositorio
2. Crear una rama para tu feature
3. Hacer commit de tus cambios
4. Crear un Pull Request

Requisitos para Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~~

* Código documentado
* Pruebas unitarias
* Sin errores de linter
* Actualización de documentación

Pruebas
------

Ejecutar Pruebas
~~~~~~~~~~~~~

.. code-block:: bash

   # Instalar dependencias de desarrollo
   pip install -e ".[dev]"
   
   # Ejecutar pruebas
   pytest
   
   # Ejecutar con cobertura
   pytest --cov=cha_md_ba

Estructura de Pruebas
~~~~~~~~~~~~~~~~~~

::

   tests/
   ├── unit/              # Pruebas unitarias
   ├── integration/       # Pruebas de integración
   └── fixtures/         # Datos de prueba

Despliegue
--------

Proceso de Despliegue
~~~~~~~~~~~~~~~~~~

1. Actualizar versión en setup.py
2. Actualizar CHANGELOG.md
3. Crear tag de versión
4. Construir y subir a PyPI

.. code-block:: bash

   # Construir paquete
   python setup.py sdist bdist_wheel
   
   # Subir a PyPI
   twine upload dist/* 