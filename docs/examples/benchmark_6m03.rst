Benchmark 6M03 en agua con NaCl 0.5 M a 310 K
=============================================

Este ejemplo prepara un sistema de referencia de la proteasa principal de
SARS-CoV-2 en forma apo (PDB ``6M03``), disuelta en agua TIP3P con **NaCl 0.5 M**
a **310 K**.

Condiciones
-----------

* Campo de fuerzas: AMBER99SB-ILDN
* Modelo de agua: TIP3P
* Caja dodecaédrica con 1.2 nm de margen
* Neutralización + NaCl 0.5 M (``gmx genion -neutral -conc 0.5``)
* Temperatura: 310 K
* Presión: 1 bar

Ejecución
---------

.. code-block:: bash

   python scripts/run_benchmark_6m03.py \
       --output-dir simulations \
       --data-dir data

Desde Python:

.. code-block:: python

   from cha_md_ba.benchmark import Benchmark6M03Config, run_benchmark

   config = Benchmark6M03Config(
       temperature=310.0,
       ion_concentration=0.5,
   )
   result = run_benchmark(config=config, stages=["download", "clean", "prepare", "mdps"])

El reporte de composición y las rutas generadas quedan en
``simulations/6M03/benchmark_report.json`` (directorio local, gitignored).
Una copia portable de la corrida de referencia (PDB, gro, topología, logs,
``tpr`` y reportes JSON; **sin** trayectorias ``.trr``) está en
``benchmarks/6M03_nacl_0.5M_310K/run/``. La descripción completa del protocolo
está en ``benchmarks/6M03_nacl_0.5M_310K/README.md``.
