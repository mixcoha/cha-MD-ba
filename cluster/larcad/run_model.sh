#!/bin/bash
# Pipeline GROMACS para un modelo 6M03 en un nodo (preparación → EM → NVT → NPT → 10 ns).
# Uso: ./run_model.sh 6M03_H41A
set -euo pipefail

MODEL="${1:?indica el modelo: 6M03_H41A, 6M03_C145A o 6M03_H41A_C145A}"
ROOT="$(cd "$(dirname "$0")" && pwd)"
RUNS="${RUNS_DIR:-$ROOT/runs}"
MDP="$ROOT/mdp"
PDB="$ROOT/models/${MODEL}/protein.pdb"
MAXWARN="${LARCAD_MAXWARN:-3}"

if [[ -z "${GMX:-}" ]]; then
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
unset DISPLAY || true

# GPU solo en dinámica (NVT/NPT/MD). EM con steep + -bonded gpu suele ser fatal.
# Tampoco usamos -bonded gpu: choca con POSRES en varios builds de GROMACS 2026.
GPU_FLAGS=()
if [[ "${LARCAD_USE_GPU:-0}" == "1" ]]; then
  if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    GPU_FLAGS=(-nb gpu -pme gpu)
  else
    echo "WARN: LARCAD_USE_GPU=1 pero no hay GPU visible; se usa CPU"
  fi
fi

gmx_cmd() {
  echo "+ $GMX $*"
  "$GMX" "$@"
}

run_mdrun() {
  local deffnm="$1"
  shift
  local flags=()
  if [[ "$deffnm" != "em" && ${#GPU_FLAGS[@]} -gt 0 ]]; then
    flags=("${GPU_FLAGS[@]}")
  fi

  local launcher=()
  case "${LARCAD_MPI_LAUNCH:-auto}" in
    srun) launcher=(srun --ntasks=1 --cpu-bind=none) ;;
    mpirun) launcher=(mpirun -np 1) ;;
    none|direct|auto) launcher=() ;;
    *)
      echo "LARCAD_MPI_LAUNCH inválido: ${LARCAD_MPI_LAUNCH}" >&2
      exit 1
      ;;
  esac

  echo "+ ${launcher[*]:-} $GMX mdrun -deffnm $deffnm -ntomp $NT ${flags[*]:-} $*"
  # gmx_mpi en un solo nodo suele arrancar como proceso único; srun/mpirun
  # se activan con LARCAD_MPI_LAUNCH si el binario exige launcher MPI.
  "${launcher[@]}" "$GMX" mdrun -deffnm "$deffnm" -ntomp "$NT" "${flags[@]}" "$@"
}

copy_topol() {
  local dest="$1"
  local fc="${2:-}"
  cp "$PREP/topol.top" "$dest/topol.top"
  local itp name
  shopt -s nullglob
  for itp in "$PREP"/*.itp; do
    name="$(basename "$itp")"
    if [[ "$name" == posre.itp ]]; then
      continue
    fi
    cp "$itp" "$dest/$name"
  done
  if [[ -n "$fc" ]]; then
    sed "s/1000/${fc}/g" "$PREP/posre.itp" > "$dest/posre.itp"
  else
    cp "$PREP/posre.itp" "$dest/posre.itp"
  fi
}

BASE="$RUNS/$MODEL"
PREP="$BASE/1_preparation"
EM="$BASE/2_minimization"
NVT="$BASE/3_nvt"
NPT="$BASE/4_npt"
PROD="$BASE/5_production"
mkdir -p "$PREP" "$EM" "$NVT/posre_constante" "$NPT" "$PROD"

echo "==> $MODEL  GMX=$GMX ($(command -v "$GMX"))  hilos=$NT  gpu=${GPU_FLAGS[*]:-cpu}  $BASE"

if [[ ! -f "$PREP/${MODEL}_ions.gro" || ! -f "$PREP/topol.top" ]]; then
  echo "--> preparación (pdb2gmx, caja, solvente, NaCl 0.15 M)"
  cp "$PDB" "$PREP/${MODEL}.pdb"
  (
    cd "$PREP"
    gmx_cmd pdb2gmx -f "${MODEL}.pdb" -o "${MODEL}.gro" -p topol.top -i posre.itp \
      -ff amber99sb-ildn -water tip3p -ignh
    gmx_cmd editconf -f "${MODEL}.gro" -o "${MODEL}_box.gro" -c -bt dodecahedron -d 1.2
    gmx_cmd solvate -cp "${MODEL}_box.gro" -cs spc216.gro -o "${MODEL}_solv.gro" -p topol.top
    gmx_cmd grompp -f "$MDP/em.mdp" -c "${MODEL}_solv.gro" -p topol.top -o ions.tpr -maxwarn "$MAXWARN"
    printf 'SOL\n' | gmx_cmd genion -s ions.tpr -o "${MODEL}_ions.gro" -p topol.top \
      -pname NA -nname CL -neutral -conc 0.15
  )
else
  echo "--> preparación ya existe, se reutiliza"
fi

if [[ ! -f "$EM/minimized.gro" ]]; then
  echo "--> minimización (CPU; steep no lleva flags GPU)"
  (
    cd "$EM"
    gmx_cmd grompp -f "$MDP/em.mdp" -c "$PREP/${MODEL}_ions.gro" -p "$PREP/topol.top" \
      -o em.tpr -maxwarn "$MAXWARN"
    run_mdrun em
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
    copy_topol "$FC_DIR" "$FC"
    if [[ "$FIRST" -eq 1 ]]; then
      gmx_cmd grompp -f "$MDP/nvt.mdp" -c "$PREV_GRO" -r "$PREV_GRO" -p topol.top \
        -o nvt.tpr -maxwarn "$MAXWARN"
    else
      local_t=()
      if [[ -n "$PREV_CPT" && -f "$PREV_CPT" ]]; then
        local_t=(-t "$PREV_CPT")
      fi
      gmx_cmd grompp -f "$MDP/nvt_cont.mdp" -c "$PREV_GRO" -r "$PREV_GRO" "${local_t[@]}" \
        -p topol.top -o nvt.tpr -maxwarn "$MAXWARN"
    fi
    run_mdrun nvt
  )
  PREV_GRO="$FC_DIR/nvt.gro"
  PREV_CPT="$FC_DIR/nvt.cpt"
  FIRST=0
done

if [[ ! -f "$NPT/npt.gro" ]]; then
  echo "--> NPT 310 K, 1 bar, 100 ps (C-rescale + POSRES 200)"
  (
    cd "$NPT"
    copy_topol "$NPT" 200
    npt_t=()
    if [[ -n "$PREV_CPT" && -f "$PREV_CPT" ]]; then
      npt_t=(-t "$PREV_CPT")
    fi
    gmx_cmd grompp -f "$MDP/npt.mdp" -c "$PREV_GRO" -r "$PREV_GRO" "${npt_t[@]}" \
      -p topol.top -o npt.tpr -maxwarn "$MAXWARN"
    run_mdrun npt
  )
else
  echo "--> NPT ya existe, se reutiliza"
fi

if [[ ! -f "$PROD/md.gro" ]]; then
  echo "--> producción 10 ns"
  (
    cd "$PROD"
    gmx_cmd grompp -f "$MDP/md.mdp" -c "$NPT/npt.gro" -t "$NPT/npt.cpt" -p "$PREP/topol.top" \
      -o md.tpr -maxwarn "$MAXWARN"
    run_mdrun md
  )
else
  echo "--> producción ya existe"
fi

echo "==> listo $MODEL"
