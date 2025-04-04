Estructura del Proyecto
====================

Organización de Directorios
------------------------

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
------------------

* Seguir PEP 8 para estilo de código
* Documentar todas las funciones y clases
* Incluir pruebas unitarias
* Mantener la cobertura de código > 80%

Ejemplo de Documentación
---------------------

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