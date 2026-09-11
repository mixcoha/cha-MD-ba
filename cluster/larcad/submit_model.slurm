#!/bin/bash
#SBATCH --job-name=6M03
#SBATCH --partition=larcad
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Encola con submit_all.sh (ahí se fijan partición, cuenta y MODEL).
# Ajusta módulos en env.sh después de: module avail gromacs && sinfo

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs

if [[ -f env.sh ]]; then
  # shellcheck disable=SC1091
  source env.sh
fi

module purge 2>/dev/null || true
if [[ -n "${LARCAD_MODULES:-}" ]]; then
  # gcc + MPI + GROMACS, en ese orden (module avail de LARCAD).
  for mod in ${LARCAD_MODULES}; do
    module load "$mod"
  done
elif [[ -n "${LARCAD_GMX_MODULE:-}" ]]; then
  module load ${LARCAD_GMX_MODULE}
fi
module list || true

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
MODEL="${MODEL:?export MODEL=6M03_H41A (o usa submit_all.sh)}"

echo "host=$(hostname)  job=${SLURM_JOB_ID:-local}  model=$MODEL  ntomp=$OMP_NUM_THREADS"
command -v gmx || command -v gmx_mpi
./run_model.sh "$MODEL"
