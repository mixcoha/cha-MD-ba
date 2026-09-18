#!/bin/bash
# Encola uno o más modelos. Crea el directorio de trabajo y NO usa --export=ALL.
# Uso: bash submit_one.sh 6M03_H41A
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p logs

if [[ -f env.sh ]]; then
  # shellcheck disable=SC1091
  source env.sh
fi

if [[ $# -eq 0 ]]; then
  echo "Uso: bash submit_one.sh 6M03_H41A [6M03_C145A ...]" >&2
  exit 1
fi

EXTRA=()
[[ -n "${LARCAD_PARTITION:-}" ]] && EXTRA+=(--partition="$LARCAD_PARTITION")
[[ -n "${LARCAD_ACCOUNT:-}" ]] && EXTRA+=(--account="$LARCAD_ACCOUNT")
[[ -n "${LARCAD_QOS:-}" ]] && EXTRA+=(--qos="$LARCAD_QOS")
[[ -n "${LARCAD_TIME:-}" ]] && EXTRA+=(--time="$LARCAD_TIME")
[[ -n "${LARCAD_CPUS:-}" ]] && EXTRA+=(--cpus-per-task="$LARCAD_CPUS")
[[ -n "${LARCAD_NODES:-}" ]] && EXTRA+=(--nodes="$LARCAD_NODES")
[[ -n "${LARCAD_GRES:-}" ]] && EXTRA+=(--gres="$LARCAD_GRES")

for model in "$@"; do
  if [[ ! -f "models/${model}/protein.pdb" ]]; then
    echo "Falta models/${model}/protein.pdb" >&2
    exit 1
  fi
  echo "sbatch ${EXTRA[*]:-} --job-name=$model --export=MODEL=$model"
  sbatch "${EXTRA[@]}" \
    --job-name="$model" \
    --export=MODEL="$model" \
    --output="${model}_%j.out" \
    --error="${model}_%j.err" \
    submit_model.slurm
done
