# Espacio de trabajo local

Esta carpeta es para **corridas de dinámica molecular en tu máquina**.
Git la ignora por completo (salvo este README): nada de lo que escribas
aquí se sube a GitHub.

## Qué va aquí

| Ruta | Contenido |
| --- | --- |
| `work/data/` | PDB descargados (p. ej. `6M03.pdb`) |
| `work/6M03/` | Preparación, minimización, NVT, NPT, producción |

No copies `.gro`, `.trr`, `.xtc`, `.log` ni reportes de corrida a
`benchmarks/` ni a ningún otro directorio versionado.

## Cómo clonar y trabajar en local

```bash
git clone https://github.com/mixcoha/cha-MD-ba.git
cd cha-MD-ba
git fetch origin cursor/6m03-benchmark-nacl-310k-c08b
git checkout cursor/6m03-benchmark-nacl-310k-c08b
```

Instala el paquete (o usa `PYTHONPATH=python_version`) y GROMACS, y lanza:

```bash
python scripts/run_benchmark_6m03.py --gpu-ids 0
```

Por defecto escribe en `work/` y `work/data/`. Para continuar una
equilibración ya empezada:

```bash
python scripts/run_benchmark_6m03.py \
    --stages nvt,npt \
    --gpu-ids 0
```

`data/` y `simulations/` siguen ignorados por compatibilidad; el
pipeline nuevo usa `work/`.
