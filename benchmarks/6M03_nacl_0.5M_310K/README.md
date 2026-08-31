# Benchmark 6M03 — agua + NaCl 0.5 M a 310 K

Protocolo de referencia de CHA-MD-BA para la proteasa principal de SARS-CoV-2
en forma apo ([PDB 6M03](https://www.rcsb.org/structure/6M03)), disuelta en
agua TIP3P con **NaCl 0.5 M** a **310 K** y 1 bar.

## Condiciones

| Parámetro | Valor |
| --- | --- |
| Estructura | 6M03, Mpro apo, cadena A (306 residuos) |
| Campo de fuerzas | AMBER99SB-ILDN |
| Agua | TIP3P |
| Caja | dodecaedro, 1.2 nm al borde |
| Iones | NaCl 0.5 M + neutralización |
| Temperatura | 310 K (V-rescale, grupos Protein / Non-Protein) |
| Presión | 1 bar (Parrinello–Rahman, NPT y producción) |
| Paso de tiempo | 2 fs |

## Cómo ejecutarlo

Desde la raíz del repositorio:

```bash
python scripts/run_benchmark_6m03.py \
    --output-dir simulations \
    --data-dir data
```

Etapas por defecto: `download`, `clean`, `prepare`, `mdps`, `minimize`.

Para generar solo los `.mdp` a 310 K (sin GROMACS):

```bash
python scripts/run_benchmark_6m03.py --stages mdps
```

Para continuar con equilibración (requiere GPU o un cluster):

```bash
python scripts/run_benchmark_6m03.py --stages download,clean,prepare,mdps,minimize,nvt,npt --gpu-ids 0
```

## Pipeline

1. Descarga de 6M03 desde RCSB.
2. Limpieza: se conservan ATOM de proteína; se descartan HOH cristalográficas.
3. `pdb2gmx` (AMBER99SB-ILDN, TIP3P, `-ignh`).
4. Caja dodecaédrica centrada (`-d 1.2`).
5. Solvatación y `genion -neutral -conc 0.5`.
6. Minimización (steepest descent, `emtol = 1000 kJ mol⁻¹ nm⁻¹`).
7. NVT a 310 K con POSRES decrecientes (1000 → 200 kJ mol⁻¹ nm⁻²).
8. NPT a 310 K y 1 bar.
9. Producción NPT (10 ns en el `.mdp` de referencia; ampliar según el recurso).

Los archivos `.mdp` generados quedan en `simulations/6M03/protocol/` y una copia de
referencia en `mdp/` de este directorio.
La composición del sistema (n.º de SOL, NA, CL) se guarda en
`simulations/6M03/benchmark_report.json`.

## Resultados de la corrida (tras `git pull`)

`data/` y `simulations/` están en `.gitignore`: una corrida en el Cloud Agent
**no** aparece en tu laptop aunque hagas pull. Los artefactos portables de esta
corrida (preparación + minimización + NVT 100 ps, POSRES 1000) están versionados
en este directorio:

```
benchmarks/6M03_nacl_0.5M_310K/run/
```

Tras `git pull` de la rama `cursor/6m03-benchmark-nacl-310k-c08b` verás, entre otros:

* `run/pdb/6M03_rcsb.pdb` y `run/pdb/6M03.pdb`
* `run/1_preparation/6M03_ions.gro`, `topol.top`, `posre.itp`
* `run/2_minimization/minimized.gro`, `em.log`, `em.edr`
* `run/3_nvt/posre_constante/1000/nvt.gro`, `nvt.log`, `nvt.edr`, `topol.tpr`
* `run/benchmark_report.json` y `run/timing_local.json`

**No** está en git la trayectoria `nvt.trr` (~158 MB; GitHub rechaza archivos
≥ 100 MB). El detalle de lo omitido está en `run/README.md`.

## Sistema de referencia verificado

Con GROMACS 2023.3, AMBER99SB-ILDN y TIP3P, el protocolo produce:

* 306 residuos (cadena A), 4682 átomos de proteína con hidrógenos
* 21 061 moléculas de agua
* 215 Na⁺ y 211 Cl⁻ (carga neta 0; ~0.5 M en una caja de ~700 nm³)
* `grompp` NVT/NPT con `ref_t = 310 K` y POSRES (`posre.itp` junto a `topol.top`)

