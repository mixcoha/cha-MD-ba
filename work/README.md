# Espacio de trabajo local

Esta carpeta es para **corridas de dinámica molecular en tu máquina**.
Git la ignora por completo (salvo este README): nada de lo que escribas
aquí se sube a GitHub.

## Qué va aquí

| Ruta | Contenido |
| --- | --- |
| `work/data/` | PDB descargados y mutados (p. ej. `6M03.pdb`, `6M03_H41A.pdb`) |
| `work/6M03/` | Silvestre: preparación, minimización, NVT, NPT, producción |
| `work/6M03_H41A/` | Modelo 1: H41A |
| `work/6M03_C145A/` | Modelo 2: C145A |
| `work/6M03_H41A_C145A/` | Modelo 3: H41A + C145A |

No copies `.gro`, `.trr`, `.xtc`, `.log` ni reportes de corrida a
`benchmarks/` ni a ningún otro directorio versionado.

## Cómo clonar y trabajar en local

```bash
git clone https://github.com/mixcoha/cha-MD-ba.git
cd cha-MD-ba
git fetch origin cursor/6m03-mutantes-h41a-c145a-3810
git checkout cursor/6m03-mutantes-h41a-c145a-3810
```

Necesitas GROMACS (p. ej. 2023.3) y, si hay NVIDIA, el driver. Luego:

```bash
python scripts/run_benchmark_6m03.py --resume
python scripts/run_benchmark_6m03.py --model 1 --resume
python scripts/run_benchmark_6m03.py --model 2 --resume
python scripts/run_benchmark_6m03.py --model 3 --resume
```

`--resume` mira la carpeta del modelo (`work/6M03/` o `work/6M03_H41A/`, etc.)
y continúa desde lo que falte.
La GPU 0 se usa sola si `nvidia-smi` ve una tarjeta; si no, corre en CPU.
Para forzar CPU: `--gpu-ids none`. Para una GPU concreta: `--gpu-ids 0`.

Por defecto escribe en `work/` y `work/data/`.

Si interrumpes un `mdrun`, vuelve a lanzar el mismo comando `--resume`:
las constantes NVT que ya tengan `nvt.gro` + `nvt.edr` no se repiten.

`data/` y `simulations/` siguen ignorados por compatibilidad; el
pipeline nuevo usa `work/`.
