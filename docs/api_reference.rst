Referencia de API
================

Este documento proporciona una referencia detallada de todas las funciones y clases disponibles en CHA-MD-BA.

Módulos Principales
------------------

.. toctree::
   :maxdepth: 2

   api/io
   api/analysis
   api/visualization
   api/utils

Módulo io
---------

.. automodule:: cha_md_ba.io
   :members:
   :undoc-members:
   :show-inheritance:

Módulo analysis
--------------

.. automodule:: cha_md_ba.analysis
   :members:
   :undoc-members:
   :show-inheritance:

Módulo visualization
-------------------

.. automodule:: cha_md_ba.visualization
   :members:
   :undoc-members:
   :show-inheritance:

Módulo utils
-----------

.. automodule:: cha_md_ba.utils
   :members:
   :undoc-members:
   :show-inheritance:

Tipos de Datos
-------------

Estructura
~~~~~~~~~

.. autoclass:: cha_md_ba.Structure
   :members:
   :undoc-members:
   :show-inheritance:

Trayectoria
~~~~~~~~~~

.. autoclass:: cha_md_ba.Trajectory
   :members:
   :undoc-members:
   :show-inheritance:

Funciones de Utilidad
--------------------

Carga de Datos
~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.load_structure
.. autofunction:: cha_md_ba.load_trajectory

Análisis Estructural
~~~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.calculate_radius_of_gyration
.. autofunction:: cha_md_ba.calculate_contact_matrix

Análisis de Trayectoria
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.calculate_rmsd
.. autofunction:: cha_md_ba.calculate_rmsf
.. autofunction:: cha_md_ba.calculate_correlation
.. autofunction:: cha_md_ba.cluster_structures

Visualización
~~~~~~~~~~~~

.. autofunction:: cha_md_ba.visualize_structure
.. autofunction:: cha_md_ba.visualize_trajectory 