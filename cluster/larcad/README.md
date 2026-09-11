# Paquete LARCAD — mutantes 6M03 (H41A, C145A, doble)

Carpeta para **copiar al nodo** del Laboratorio Regional de Cómputo de Alto
Desempeño (UNACH) y encolar los tres modelos con GROMACS.

No sube trayectorias a GitHub. En el nodo las corridas quedan en `runs/`.

## Qué hay que subir

| Ruta | Contenido |
| --- | --- |
| `models/6M03_H41A/protein.pdb` | Modelo 1: His41→Ala |
| `models/6M03_C145A/protein.pdb` | Modelo 2: Cys145→Ala |
| `models/6M03_H41A_C145A/protein.pdb` | Modelo 3: ambas |
| `mdp/` | `em.mdp`, `nvt.mdp`, `nvt_cont.mdp`, `npt.mdp`, `md.mdp` |
| `run_model.sh` | Preparación + EM + NVT (1000→200) + NPT + 10 ns |
| `submit_model.slurm` / `submit_all.sh` | Envío SLURM |

## Desde tu laptop

```bash
git checkout cursor/6m03-mutantes-h41a-c145a-3810
python3 scripts/pack_larcad_6m03.py
bash scripts/upload_larcad.sh
```

El login del nodo es `ssh -p 212 mixcoha@148.222.27.130`. El script ya usa ese host y puerto.

Si no usas el script:

```bash
rsync -avz -e "ssh -p 212" cluster/larcad/ mixcoha@148.222.27.130:~/cha-md-ba-6m03/
```

El tar.gz (sin `runs/`) queda en `work/larcad_bundle/cha-md-ba-6m03-larcad.tar.gz`.

## En el nodo

```bash
ssh -p 212 mixcoha@148.222.27.130
cd ~/cha-md-ba-6m03
cp env.sh.example env.sh
# Colas: larcad (n1–n4 CPU) y gpu_rtxA5000 (gpu1, RTX A5000).
# Módulos: gromacs-mpi-2026.2  /  gromacs-mpi-cuda-2026.2
cp env.sh.example env.sh
bash submit_all.sh                 # los tres mutantes (GPU por defecto)
squeue -u "$USER"
```

Un solo modelo: `MODEL=6M03_H41A sbatch --export=ALL,MODEL=6M03_H41A submit_model.slurm`

Si el scheduler no es SLURM, ejecuta a mano en el nodo asignado:

```bash
export OMP_NUM_THREADS=8
./run_model.sh 6M03_H41A
```

`run_model.sh` reanuda etapas ya hechas (útil si el job se corta).

## Protocolo

Igual que el benchmark local: AMBER99SB-ILDN, TIP3P, dodecaedro 1.2 nm,
NaCl 0.15 M, 310 K, POSRES 1000→200, NPT 100 ps, producción 10 ns.
Los resultados **no** se versionan; cópialos de vuelta con `rsync` cuando
terminen (`runs/` y `logs/`).
