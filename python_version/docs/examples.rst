Ejemplos
========

Esta sección proporciona ejemplos prácticos de uso de CHA-MD-BA.

Uso Básico
----------

Preparación del Sistema
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import prepare

   # Inicializar el preparador
   preparator = prepare.SystemPreparator(
       input_file="protein.pdb",
       force_field="amber99sb-ildn",
       water_model="tip3p"
   )

   # Preparar el sistema
   preparator.prepare()

Minimización de Energía
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import minimize

   # Inicializar el minimizador
   minimizer = minimize.Minimizer(
       max_steps=50000,
       emtol=1000.0
   )

   # Realizar minimización
   minimizer.minimize()

Equilibración NVT
~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import nvt

   # Inicializar el equilibrador NVT
   nvt_equilibrator = nvt.NVTEquilibrator(
       force_constant=1000.0,
       temperature=300.0
   )

   # Realizar equilibración NVT
   nvt_equilibrator.equilibrate()

Equilibración NPT
~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import npt

   # Inicializar el equilibrador NPT
   npt_equilibrator = npt.NPTEquilibrator(
       pressure=1.0,
       temperature=300.0
   )

   # Realizar equilibración NPT
   npt_equilibrator.equilibrate()

Producción
~~~~~~~~~

.. code-block:: python

   from cha_md_ba import production

   # Inicializar el ejecutor de producción
   prod = production.ProductionRunner(
       time=100.0,  # 100 ns
       dt=0.002
   )

   # Ejecutar la simulación de producción
   prod.run()

Análisis
~~~~~~~~

.. code-block:: python

   from cha_md_ba import analysis

   # Inicializar el analizador
   analyzer = analysis.Analyzer(
       reference="reference.pdb",
       trajectory="traj.xtc"
   )

   # Calcular RMSD
   rmsd = analyzer.calculate_rmsd()

   # Calcular RMSF
   rmsf = analyzer.calculate_rmsf()

   # Analizar estructura secundaria
   ss = analyzer.analyze_secondary_structure()

   # Analizar enlaces de hidrógeno
   hbonds = analyzer.analyze_hbonds()

   # Realizar análisis de clusters
   clusters = analyzer.cluster_analysis()

Uso Avanzado
-----------

Campo de Fuerza Personalizado
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import prepare

   # Inicializar con campo de fuerza personalizado
   preparator = prepare.SystemPreparator(
       input_file="protein.pdb",
       force_field="charmm36",
       water_model="tip3p"
   )

   # Agregar parámetros personalizados
   preparator.add_custom_parameters("custom.itp")
   preparator.prepare()

Simulación con Restricciones
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import nvt

   # Inicializar con restricciones de posición
   nvt_equilibrator = nvt.NVTEquilibrator(
       force_constant=1000.0,
       temperature=300.0,
       restraints=["backbone", "CA"]
   )

   # Realizar equilibración con restricciones
   nvt_equilibrator.equilibrate()

Rampa de Temperatura
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import nvt

   # Inicializar con rampa de temperatura
   nvt_equilibrator = nvt.NVTEquilibrator(
       force_constant=1000.0,
       temperature=300.0,
       temperature_ramp={
           "start": 100.0,
           "end": 300.0,
           "steps": 50000
       }
   )

   # Realizar rampa de temperatura
   nvt_equilibrator.equilibrate()

Análisis Personalizado
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from cha_md_ba import analysis

   # Inicializar con selecciones personalizadas
   analyzer = analysis.Analyzer(
       reference="reference.pdb",
       trajectory="traj.xtc",
       selections={
           "protein": "protein",
           "backbone": "backbone",
           "sidechains": "not backbone"
       }
   )

   # Calcular propiedades personalizadas
   properties = analyzer.calculate_properties(
       properties=["radius_of_gyration", "sasa"],
       selections=["protein", "backbone"]
   )

Ejemplos de Línea de Comandos
---------------------------

Preparación del Sistema
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cha-md-ba prepare \
       --input protein.pdb \
       --force-field amber99sb-ildn \
       --water-model tip3p \
       --box-type dodecahedron \
       --box-distance 1.0

Minimización de Energía
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cha-md-ba minimize \
       --max-steps 50000 \
       --emtol 1000.0 \
       --emstep 0.01

Equilibración NVT
~~~~~~~~~~~~~~~~

.. code-block:: bash

   cha-md-ba nvt \
       --force-constant 1000.0 \
       --temperature 300.0 \
       --dt 0.002 \
       --nsteps 50000

Equilibración NPT
~~~~~~~~~~~~~~~~

.. code-block:: bash

   cha-md-ba npt \
       --pressure 1.0 \
       --temperature 300.0 \
       --dt 0.002 \
       --nsteps 50000

Producción
~~~~~~~~~

.. code-block:: bash

   cha-md-ba production \
       --time 100.0 \
       --dt 0.002 \
       --temperature 300.0 \
       --pressure 1.0

Análisis
~~~~~~~~

.. code-block:: bash

   cha-md-ba analyze \
       --type rmsd \
       --reference reference.pdb \
       --trajectory traj.xtc \
       --selection protein 