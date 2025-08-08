Guías de Usuario
==============

.. toctree::
   :maxdepth: 2
   :caption: Contenido:

   instalacion
   primer_paso
   configuracion
   analisis_avanzado

Instalación
----------

Requisitos del Sistema
~~~~~~~~~~~~~~~~~~~~

* Python 3.8 o superior
* GROMACS 2020 o superior
* NumPy
* Pandas
* Matplotlib
* MDAnalysis

Instalación desde PyPI
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   pip install cha-md-ba

Instalación desde el Código Fuente
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/tu-usuario/cha-MD-ba.git
   cd cha-MD-ba
   pip install -e .

Primer Paso
----------

Configuración Inicial
~~~~~~~~~~~~~~~~~~

1. Asegúrate de tener GROMACS instalado y configurado en tu sistema
2. Configura las variables de entorno necesarias
3. Verifica la instalación ejecutando:

.. code-block:: bash

   cha-md-ba --version

Ejecución Básica
~~~~~~~~~~~~~~

1. Prepara tu estructura PDB
2. Ejecuta el pipeline básico:

.. code-block:: bash

   cha-md-ba prepare -i proteina.pdb -o sistema_preparado
   cha-md-ba simulate -i sistema_preparado -o simulacion
   cha-md-ba analyze -i simulacion -o resultados

Configuración
-----------

Archivo de Configuración
~~~~~~~~~~~~~~~~~~~~~

El archivo de configuración se encuentra en ``~/.cha-md-ba/config.yaml``:

.. code-block:: yaml

   gromacs:
     path: /usr/local/gromacs/bin/gmx
     force_field: amber99sb-ildn
   
   simulation:
     temperature: 300
     pressure: 1
     time_step: 0.002
   
   analysis:
     rmsd: true
     rmsf: true
     sasa: true

Análisis Avanzado
--------------

Análisis de Estructura Secundaria
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba.analysis import SecondaryStructure
   
   analyzer = SecondaryStructure()
   analyzer.analyze("trajectory.xtc")
   analyzer.plot()

Análisis de Energía
~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba.analysis import EnergyAnalysis
   
   analyzer = EnergyAnalysis()
   analyzer.analyze("ener.edr")
   analyzer.plot() 