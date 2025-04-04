Módulo io
========

Este módulo proporciona funciones para la entrada y salida de datos en CHA-MD-BA.

Clases
------

SystemPreparator
~~~~~~~~~~~~~~

.. autoclass:: cha_md_ba.io.SystemPreparator
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__

   .. automethod:: prepare

   .. automethod:: generate_topology

   .. automethod:: solvate

   .. automethod:: add_ions

   .. automethod:: add_custom_parameters

Funciones
--------

load_structure
~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.load_structure

load_trajectory
~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.load_trajectory

save_structure
~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.save_structure

save_trajectory
~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.save_trajectory

create_directory
~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.create_directory

remove_directory
~~~~~~~~~~~~~~

.. autofunction:: cha_md_ba.io.remove_directory

copy_file
~~~~~~~~

.. autofunction:: cha_md_ba.io.copy_file 