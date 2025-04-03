Guía de Usuario
==============

Esta guía proporciona una introducción al uso de CHA-MD-BA para el análisis de dinámica molecular de proteínas.

Conceptos Básicos
----------------

CHA-MD-BA está diseñado para trabajar con trayectorias de dinámica molecular y estructuras de proteínas. Los conceptos principales incluyen:

* **Estructura**: Representación de la proteína en un momento dado
* **Trayectoria**: Serie de estructuras a lo largo del tiempo
* **Análisis**: Cálculo de propiedades y métricas de la proteína

Uso Básico
----------

Importar el paquete:

.. code-block:: python

   import cha_md_ba as cmb

Cargar una estructura:

.. code-block:: python

   structure = cmb.load_structure("protein.pdb")

Cargar una trayectoria:

.. code-block:: python

   trajectory = cmb.load_trajectory("trajectory.xtc", structure)

Análisis de Estructura
---------------------

Calcular el radio de giro:

.. code-block:: python

   rg = cmb.calculate_radius_of_gyration(structure)

Calcular la matriz de contactos:

.. code-block:: python

   contacts = cmb.calculate_contact_matrix(structure)

Análisis de Trayectoria
----------------------

Calcular RMSD:

.. code-block:: python

   rmsd = cmb.calculate_rmsd(trajectory)

Calcular RMSF:

.. code-block:: python

   rmsf = cmb.calculate_rmsf(trajectory)

Visualización
------------

Visualizar la estructura:

.. code-block:: python

   cmb.visualize_structure(structure)

Visualizar la trayectoria:

.. code-block:: python

   cmb.visualize_trajectory(trajectory)

Ejemplos Avanzados
-----------------

Análisis de correlación:

.. code-block:: python

   correlation = cmb.calculate_correlation(trajectory)

Análisis de clusters:

.. code-block:: python

   clusters = cmb.cluster_structures(trajectory)

Solución de Problemas
--------------------

Si encuentra algún problema, consulte la sección de :doc:`faq` o abra un issue en el repositorio de GitHub. 