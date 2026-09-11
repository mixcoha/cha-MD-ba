#!/bin/bash
#SBATCH --job-name=6M03
#SBATCH --partition=larcad
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
# Los .out/.err van al directorio de envío. Si apuntan a logs/ y esa
# carpeta no existe, Slurm acepta el job y lo mata al arrancar (PD → desaparece).
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Encola con: bash submit_one.sh 6M03_H41A
# No uses --export=ALL ni sed sobre este archivo.

set -u
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs
echo "=== job ${SLURM_JOB_ID:-n/a} host=$(hostname) date=$(date) pwd=$(pwd) ==="

if [[ -f env.sh ]]; then
  # shellcheck disable=SC1091
  source env.sh
fi

echo "MODEL=${MODEL:-}  partition=${SLURM_JOB_PARTITION:-}  gres=${SLURM_JOB_GRES:-${SLURM_JOB_GPUS:-}}"
echo "ntasks=${SLURM_NTASKS:-}  cpus=${SLURM_CPUS_PER_TASK:-}"

# En jobs batch `module` a veces no existe hasta cargar el init.
if ! type module >/dev/null 2>&1; then
  for _modinit in \
    /etc/profile \
    /usr/share/Modules/init/bash \
    /usr/share/lmod/lmod/init/bash \
    /opt/ohpc/admin/lmod/lmod/init/bash
  do
    if [[ -f "$_modinit" ]]; then
      # shellcheck disable=SC1090
      source "$_modinit" && type module >/dev/null 2>&1 && break
    fi
  done
fi

load_modules() {
  local mods="$1"
  local mod
  for mod in $mods; do
    echo "module load $mod"
    if ! module load "$mod"; then
      echo "ERROR: falló module load $mod" >&2
      return 1
    fi
  done
}

GPU_STACK="${LARCAD_MODULES:-gcc-15.1.1 mpich-cuda-4.2.3 gromacs-mpi-cuda-2026.2}"
CPU_STACK="${LARCAD_CPU_MODULES:-gcc-15.1.1 mpich-4.2.3 gromacs-mpi-2026.2}"

if type module >/dev/null 2>&1; then
  module purge 2>/dev/null || true
  if ! load_modules "$GPU_STACK"; then
    echo "WARN: pila GPU falló; se intenta CPU: $CPU_STACK"
    module purge 2>/dev/null || true
    load_modules "$CPU_STACK" || true
    export LARCAD_USE_GPU=0
  fi
  module list || true
else
  echo "WARN: no hay comando module; se usa PATH tal cual"
fi

if [[ -z "${GMX:-}" ]]; then
  if command -v gmx >/dev/null 2>&1; then
    GMX=gmx
  elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX=gmx_mpi
  fi
fi
export GMX="${GMX:-gmx_mpi}"

echo "GMX=$GMX -> $(command -v "$GMX" || true)"
if ! command -v "$GMX" >/dev/null 2>&1; then
  echo "ERROR: no hay $GMX después de module load" >&2
  echo "PATH=$PATH" >&2
  type module >&2 || true
  module avail gromacs 2>&1 || true
  exit 1
fi

"$GMX" --version 2>&1 | head -25 || true
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi -L || true
else
  echo "nvidia-smi no está en PATH"
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
if [[ -z "${MODEL:-}" ]]; then
  echo "ERROR: MODEL vacío. Encola con: bash submit_one.sh 6M03_H41A" >&2
  exit 1
fi
unset DISPLAY || true

exec ./run_model.sh "$MODEL"
