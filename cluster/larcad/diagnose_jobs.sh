#!/bin/bash
# En el nodo:  bash diagnose_jobs.sh
# Imprime cola, sacct y las colas de logs/.err que suelen explicar por qué
# se cayeron los jobs (módulos, gmx_mpi, GPU, LINCS).
set -u
cd "$(dirname "$0")"

echo "=== host $(hostname)  $(date -Is)  pwd=$(pwd) ==="
echo
echo "=== squeue ==="
squeue -u "${USER}" -o '%.10i %.12P %.18j %.8u %.8T %.10M %.6D %R' 2>/dev/null || echo "(squeue no disponible)"

echo
echo "=== sacct (48 h) ==="
sacct -u "${USER}" --starttime=now-2days \
  --format=JobID,JobName,State,ExitCode,Elapsed,NodeList,Partition,ReqTRES -P 2>/dev/null \
  || echo "(sacct no disponible)"

echo
echo "=== logs/ ==="
if [[ -d logs ]]; then
  ls -lt logs | head -30
else
  echo "no existe logs/ — Slurm no pudo escribir stdout/stderr"
fi

echo
for f in logs/*.err logs/*.out; do
  [[ -e "$f" ]] || continue
  echo "---------- $f (últimas 80 líneas) ----------"
  tail -80 "$f"
  echo
done

echo "=== runs/ (logs de GROMACS) ==="
if [[ -d runs ]]; then
  find runs -type f \( -name '*.log' -o -name 'mdout.mdp' \) | sort
  echo
  while IFS= read -r f; do
    echo "---------- $f (errores) ----------"
    grep -E -i 'fatal|error|segmentation|LINCS|cancelled|GPU|not found' "$f" | tail -25 || true
    echo
  done < <(find runs -type f -name '*.log' | sort)
else
  echo "no existe runs/ — el pipeline no llegó a escribir etapas"
fi
