#!/bin/bash
# Pipeline GROMACS para un modelo 6M03 en un nodo (preparación → EM → NVT → NPT → 10 ns).
# Uso: ./run_model.sh 6M03_H41A
set -euo pipefail

MODEL="${1:?indica el modelo: 6M03_H41A, 6M03_C145A o 6M03_H41A_C145A}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
RUNS="${RUNS_DIR:-$ROOT/runs}"
MDP="$ROOT/mdp"
PDB="$ROOT/models/${MODEL}/protein.pdb"
GMX="${GMX:-}"

if [[ -z "$GMX" ]]; then
  if command -v gmx >/dev/null 2>&1; then
    GMX=gmx
  elif command -v gmx_mpi >/dev/null 2>&1; then
    GMX=gmx_mpi
  else
    echo "No se encontró gmx ni gmx_mpi. Carga el módulo de GROMACS o exporta GMX=..." >&2
    exit 1
  fi
fi

if [[ ! -f "$PDB" ]]; then
  echo "No existe $PDB. Genera los PDB con scripts/pack_larcad_6m03.py" >&2
  exit 1
fi

NT="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-8}}"
export OMP_NUM_THREADS="$NT"
export GMX_MAXBACKUP="${GMX_MAXBACKUP:-0}"

mdrun() {
  local deffnm="$1"
  shift
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v srun >/dev/null 2>&1; then
    srun "$GMX" mdrun -deffnm "$deffnm" -ntomp "$NT" "$@"
  else
    "$GMX" mdrun -deffnm "$deffnm" -ntomp "$NT" "$@"
  fi
}

BASE="$RUNS/$MODEL"
PREP="$BASE/1_preparation"
EM="$BASE/2_minimization"
NVT="$BASE/3_nvt"
NPT="$BASE/4_npt"
PROD="$BASE/5_production"
mkdir -p "$PREP" "$EM" "$NVT/posre_constante" "$NPT" "$PROD"

echo "==> $MODEL  GMX=$GMX  hilos=$NT  $BASE"

if [[ ! -f "$PREP/${MODEL}_ions.gro" || ! -f "$PREP/topol.top" ]]; then
  echo "--> preparación (pdb2gmx, caja, solvente, NaCl 0.5 M)"
  cp "$PDB" "$PREP/${MODEL}.pdb"
  (
    cd "$PREP"
    "$GMX" pdb2gmx -f "${MODEL}.pdb" -o "${MODEL}.gro" -p topol.top -i posre.itp \
      -ff amber99sb-ildn -water tip3p -ignh
    "$GMX" editconf -f "${MODEL}.gro" -o "${MODEL}_box.gro" -c -bt dodecahedron -d 1.2
    "$GMX" solvate -cp "${MODEL}_box.gro" -cs spc216.gro -o "${MODEL}_solv.gro" -p topol.top
    "$GMX" grompp -f "$MDP/em.mdp" -c "${MODEL}_solv.gro" -p topol.top -o ions.tpr -maxwarn 1
    printf 'SOL\n' | "$GMX" genion -s ions.tpr -o "${MODEL}_ions.gro" -p topol.top \
      -pname NA -nname CL -neutral -conc 0.5
  )
else
  echo "--> preparación ya existe, se reutiliza"
fi

if [[ ! -f "$EM/minimized.gro" ]]; then
  echo "--> minimización"
  (
    cd "$EM"
    "$GMX" grompp -f "$MDP/em.mdp" -c "$PREP/${MODEL}_ions.gro" -p "$PREP/topol.top" \
      -o em.tpr -maxwarn 1
    mdrun em
    cp em.gro minimized.gro
  )
else
  echo "--> minimización ya existe, se reutiliza"
fi

PREV_GRO="$EM/minimized.gro"
PREV_CPT=""
FIRST=1
for FC in 1000 800 600 400 200; do
  FC_DIR="$NVT/posre_constante/$FC"
  mkdir -p "$FC_DIR"
  if [[ -f "$FC_DIR/nvt.gro" ]]; then
    echo "--> NVT POSRES $FC ya existe"
    PREV_GRO="$FC_DIR/nvt.gro"
    if [[ -f "$FC_DIR/nvt.cpt" ]]; then
      PREV_CPT="$FC_DIR/nvt.cpt"
    fi
    FIRST=0
    continue
  fi
  echo "--> NVT POSRES $FC (310 K, 100 ps)"
  (
    cd "$FC_DIR"
    cp "$PREP/topol.top" .
    shopt -s nullglob
    for itp in "$PREP"/*.itp; do
      name="$(basename "$itp")"
      if [[ "$name" == posre.itp ]]; then
        continue
      fi
      cp "$itp" "$name"
    done
    sed "s/1000/${FC}/g" "$PREP/posre.itp" > posre.itp
    if [[ "$FIRST" -eq 1 ]]; then
      "$GMX" grompp -f "$MDP/nvt.mdp" -c "$PREV_GRO" -r "$PREV_GRO" -p topol.top \
        -o nvt.tpr -maxwarn 1
    else
      "$GMX" grompp -f "$MDP/nvt_cont.mdp" -c "$PREV_GRO" -r "$PREV_GRO" -t "$PREV_CPT" \
        -p topol.top -o nvt.tpr -maxwarn 1
    fi
    mdrun nvt
  )
  PREV_GRO="$FC_DIR/nvt.gro"
  PREV_CPT="$FC_DIR/nvt.cpt"
  FIRST=0
done

if [[ ! -f "$NPT/npt.gro" ]]; then
  echo "--> NPT 310 K, 1 bar, 100 ps"
  (
    cd "$NPT"
    "$GMX" grompp -f "$MDP/npt.mdp" -c "$PREV_GRO" -t "$PREV_CPT" -p "$PREP/topol.top" \
      -o npt.tpr -maxwarn 1
    mdrun npt
  )
else
  echo "--> NPT ya existe, se reutiliza"
fi

if [[ ! -f "$PROD/md.gro" ]]; then
  echo "--> producción 10 ns"
  (
    cd "$PROD"
    "$GMX" grompp -f "$MDP/md.mdp" -c "$NPT/npt.gro" -t "$NPT/npt.cpt" -p "$PREP/topol.top" \
      -o md.tpr -maxwarn 1
    mdrun md
  )
else
  echo "--> producción ya existe"
fi

echo "==> listo $MODEL"
