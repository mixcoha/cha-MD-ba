Configuración
=============

Esta sección describe la configuración del sistema CHA-MD-BA.

Archivo de Configuración
----------------------

El archivo de configuración se encuentra en ``~/.cha-md-ba/config.yaml`` y contiene los siguientes parámetros:

.. code-block:: yaml

    gromacs:
        executable: /usr/local/gromacs/bin/gmx
        force_field: amber99sb-ildn

    simulation:
        temperature: 300  # K
        pressure: 1      # atm
        time_step: 0.002 # ps

    analysis:
        rmsd: true
        rmsf: true
        sasa: true

Parámetros de GROMACS
--------------------

- ``executable``: Ruta al ejecutable de GROMACS
- ``force_field``: Campo de fuerzas a utilizar

Parámetros de Simulación
-----------------------

- ``temperature``: Temperatura del sistema en Kelvin
- ``pressure``: Presión del sistema en atmósferas
- ``time_step``: Paso de tiempo en picosegundos

Parámetros de Análisis
---------------------

- ``rmsd``: Habilitar cálculo de RMSD
- ``rmsf``: Habilitar cálculo de RMSF
- ``sasa``: Habilitar cálculo de SASA 