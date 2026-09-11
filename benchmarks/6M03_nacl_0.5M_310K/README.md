# Benchmark 6M03 — agua + NaCl 0.5 M a 310 K

Protocolo de referencia de CHA-MD-BA para la proteasa principal de SARS-CoV-2
en forma apo ([PDB 6M03](https://www.rcsb.org/structure/6M03)), disuelta en
agua TIP3P con **NaCl 0.5 M** a **310 K** y 1 bar.

Las **corridas** (coordenadas, trayectorias, logs) viven en `work/` en tu
clon local y **no se suben a GitHub**. Aquí solo está el protocolo.

## Modelos de la díada catalítica

His41 y Cys145 son la díada de Mpro. Se truncan a alanina (se conservan
N, CA, C, O y CB; `pdb2gmx -ignh` completa los hidrógenos):

| Modelo | Mutación | Carpeta local |
| --- | --- | --- |
| silvestre | ninguna | `work/6M03/` |
| 1 | H41A | `work/6M03_H41A/` |
| 2 | C145A | `work/6M03_C145A/` |
| 3 | H41A + C145A | `work/6M03_H41A_C145A/` |

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

## Cómo ejecutarlo (local)

Desde la raíz del repositorio (GROMACS en PATH; GPU opcional):

```bash
python scripts/run_benchmark_6m03.py --resume
python scripts/run_benchmark_6m03.py --model 1 --resume    # H41A
python scripts/run_benchmark_6m03.py --model 2 --resume    # C145A
python scripts/run_benchmark_6m03.py --model 3 --resume    # H41A + C145A
python scripts/run_benchmark_6m03.py --model all --resume  # los tres mutantes
```

El silvestre escribe en `work/6M03/`; cada mutante tiene su carpeta.
En el clúster (LARCAD) el paquete a copiar está en `cluster/larcad/`
(`scripts/pack_larcad_6m03.py` y `scripts/upload_larcad.sh`).
`--gpu-ids auto` (por defecto) usa la GPU 0 si hay NVIDIA; `--gpu-ids none` fuerza CPU.

Para generar solo los `.mdp` a 310 K (sin GROMACS):

```bash
python scripts/run_benchmark_6m03.py --stages mdps
```

## Pipeline

1. Descarga de 6M03 desde RCSB.
2. Limpieza: se conservan ATOM de proteína; se descartan HOH cristalográficas.
3. Mutación a ALA si el modelo no es silvestre (H41A y/o C145A).
4. `pdb2gmx` (AMBER99SB-ILDN, TIP3P, `-ignh`).
5. Caja dodecaédrica centrada (`-d 1.2`).
6. Solvatación y `genion -neutral -conc 0.5`.
7. Minimización (steepest descent, `emtol = 1000 kJ mol⁻¹ nm⁻¹`).
8. NVT a 310 K con POSRES decrecientes (1000 → 200 kJ mol⁻¹ nm⁻²).
9. NPT a 310 K y 1 bar.
10. Producción NPT (10 ns en el `.mdp` de referencia; ampliar según el recurso).

Los archivos `.mdp` generados quedan en `work/6M03/protocol/` y una copia de
referencia en `mdp/` de este directorio.
La composición del sistema (n.º de SOL, NA, CL) se guarda en
`work/6M03/benchmark_report.json` (local, gitignored).

## Sistema de referencia esperado

Con GROMACS 2023.3, AMBER99SB-ILDN y TIP3P, el protocolo produce:

* 306 residuos (cadena A), 4682 átomos de proteína con hidrógenos
* ~21 000 moléculas de agua
* Na⁺ y Cl⁻ con carga neta 0 a ~0.5 M en una caja de ~700 nm³
* `grompp` NVT/NPT con `ref_t = 310 K` y POSRES (`posre.itp` junto a `topol.top`)
