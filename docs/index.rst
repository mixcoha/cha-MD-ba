CHA-MD-BA: Pipeline de Simulación de Dinámica Molecular
====================================================

Bienvenido a la documentación de CHA-MD-BA, un pipeline para simulación de dinámica molecular.

.. toctree::
   :maxdepth: 2
   :caption: Contenido:

   api/api
   examples/examples
   guides/guides
   development/development

Características Principales
------------------------

* Preparación automática de sistemas para simulación
* Integración con GROMACS
* Análisis avanzado de trayectorias
* Visualización de resultados
* Interfaz de línea de comandos intuitiva

Instalación Rápida
----------------

.. code-block:: bash

   pip install cha-md-ba

Para más detalles, consulta la :doc:`guía de instalación <guides/guides>`.

Primer Paso
---------

1. Prepara tu estructura PDB:

   .. code-block:: bash

      cha-md-ba prepare -i proteina.pdb -o sistema_preparado

2. Ejecuta la simulación:

   .. code-block:: bash

      cha-md-ba simulate -i sistema_preparado -o simulacion

3. Analiza los resultados:

   .. code-block:: bash

      cha-md-ba analyze -i simulacion -o resultados

Para más ejemplos, consulta la sección de :doc:`ejemplos <examples/examples>`.

Documentación
-----------

* :doc:`API Reference <api/api>` - Documentación detallada de la API
* :doc:`Guías de Usuario <guides/guides>` - Guías paso a paso
* :doc:`Ejemplos <examples/examples>` - Ejemplos de código
* :doc:`Guía de Desarrollo <development/development>` - Información para desarrolladores

Soporte
------

* `GitHub Issues <https://github.com/tu-usuario/cha-MD-ba/issues>`_
* `Documentación en línea <https://cha-md-ba.readthedocs.io>`_
* Email: tu-email@ejemplo.com

Índices y Tablas
--------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 