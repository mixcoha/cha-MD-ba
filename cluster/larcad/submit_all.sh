#!/bin/bash
# Encola los tres mutantes de 6M03 en LARCAD.
# Antes: copia env.sh.example a env.sh y rellena partición/módulo.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
mkdir -p logs

if [[ -f env.sh ]]; then
  # shellcheck disable=SC1091
  source env.sh
fi

EXTRA=()
[[ -n "${LARCAD_PARTITION:-}" ]] && EXTRA+=(--partition="$LARCAD_PARTITION")
[[ -n "${LARCAD_ACCOUNT:-}" ]] && EXTRA+=(--account="$LARCAD_ACCOUNT")
[[ -n "${LARCAD_QOS:-}" ]] && EXTRA+=(--qos="$LARCAD_QOS")
[[ -n "${LARCAD_TIME:-}" ]] && EXTRA+=(--time="$LARCAD_TIME")
[[ -n "${LARCAD_CPUS:-}" ]] && EXTRA+=(--cpus-per-task="$LARCAD_CPUS")
[[ -n "${LARCAD_NODES:-}" ]] && EXTRA+=(--nodes="$LARCAD_NODES")
[[ -n "${LARCAD_GRES:-}" ]] && EXTRA+=(--gres="$LARCAD_GRES")

MODELS=(6M03_H41A 6M03_C145A 6M03_H41A_C145A)
if [[ "${1:-}" == "wt" || "${1:-}" == "all-plus-wt" ]]; then
  MODELS=(6M03 "${MODELS[@]}")
fi
if [[ $# -gt 0 && "$1" != "wt" && "$1" != "all-plus-wt" && "$1" != "all" ]]; then
  MODELS=("$@")
fi

for model in "${MODELS[@]}"; do
  if [[ ! -f "models/${model}/protein.pdb" ]]; then
    echo "Falta models/${model}/protein.pdb — corre scripts/pack_larcad_6m03.py" >&2
    exit 1
  fi
  echo "sbatch ${EXTRA[*]} --job-name=$model  $model"
  sbatch "${EXTRA[@]}" --job-name="$model" --export=ALL,MODEL="$model" submit_model.slurm
done
