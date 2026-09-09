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

Necesitas GROMACS (p. ej. 2023.3) y, si hay NVIDIA, el driver. Luego:

```bash
python scripts/run_benchmark_6m03.py --resume
```

`--resume` mira `work/6M03/` y continúa desde lo que falte (preparación,
minimización, NVT POSRES 1000→200, NPT). Si no hay nada, empieza de cero.
La GPU 0 se usa sola si `nvidia-smi` ve una tarjeta; si no, corre en CPU.
Para forzar CPU: `--gpu-ids none`. Para una GPU concreta: `--gpu-ids 0`.

Por defecto escribe en `work/` y `work/data/`.

Si interrumpes un `mdrun`, vuelve a lanzar el mismo comando `--resume`:
las constantes NVT que ya tengan `nvt.gro` + `nvt.edr` no se repiten.

`data/` y `simulations/` siguen ignorados por compatibilidad; el
pipeline nuevo usa `work/`.
