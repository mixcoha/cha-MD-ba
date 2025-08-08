Módulo analysis
=============

Este módulo proporciona funciones para el análisis de estructuras y trayectorias en CHA-MD-BA.

Clases
------

Analyzer
~~~~~~~~

.. autoclass:: cha_md_ba.analysis.Analyzer
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__

   .. automethod:: calculate_rmsd

   .. automethod:: calculate_rmsf

   .. automethod:: analyze_secondary_structure

   .. automethod:: analyze_hbonds

   .. automethod:: cluster_structures

   .. automethod:: calculate_properties

Funciones
--------

calculate_radius_of_gyration
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_radius_of_gyration

calculate_contact_matrix
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_contact_matrix

calculate_correlation
~~~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_correlation

calculate_sasa
~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_sasa

calculate_dihedrals
~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_dihedrals

calculate_angles
~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_angles

calculate_distances
~~~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.analysis.calculate_distances 