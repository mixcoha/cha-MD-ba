Benchmark 6M03 en agua con NaCl 0.5 M a 310 K
=============================================

Este ejemplo prepara un sistema de referencia de la proteasa principal de
SARS-CoV-2 en forma apo (PDB ``6M03``), disuelta en agua TIP3P con **NaCl 0.5 M**
a **310 K**.

Las corridas se escriben en ``work/`` (gitignored) y no se publican en GitHub.
La descripción del protocolo está en ``benchmarks/6M03_nacl_0.5M_310K/README.md``.

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

   python scripts/run_benchmark_6m03.py --resume
   python scripts/run_benchmark_6m03.py --model 1 --resume
   python scripts/run_benchmark_6m03.py --model 2 --resume
   python scripts/run_benchmark_6m03.py --model 3 --resume
   python scripts/run_benchmark_6m03.py --model all --resume

Desde Python:

.. code-block:: python

   from cha_md_ba.benchmark import Benchmark6M03Config, run_benchmark

   config = Benchmark6M03Config(
       temperature=310.0,
       ion_concentration=0.5,
   )
   result = run_benchmark(config=config, stages=["download", "clean", "prepare", "mdps"])

El reporte de composición y las rutas generadas quedan en
``work/6M03/benchmark_report.json`` (directorio local, gitignored).
