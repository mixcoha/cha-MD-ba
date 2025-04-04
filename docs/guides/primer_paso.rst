Primer Paso
==========

Configuración Inicial
------------------

1. Asegúrate de tener GROMACS instalado y configurado en tu sistema
2. Configura las variables de entorno necesarias
3. Verifica la instalación ejecutando:

.. code-block:: bash

   cha-md-ba --version

Ejecución Básica
--------------

1. Prepara tu estructura PDB
2. Ejecuta el pipeline básico:

.. code-block:: bash

   cha-md-ba prepare -i proteina.pdb -o sistema_preparado
   cha-md-ba simulate -i sistema_preparado -o simulacion
   cha-md-ba analyze -i simulacion -o resultados 