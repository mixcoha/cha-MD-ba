#!/bin/bash
# Sube el paquete de 6M03 al nodo de LARCAD.
#   bash scripts/upload_larcad.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HOST="${LARCAD_HOST:-mixcoha@148.222.27.130}"
PORT="${LARCAD_PORT:-212}"
DEST="${LARCAD_REMOTE_DIR:-~/cha-md-ba-6m03}"
CLUSTER="$ROOT/cluster/larcad"
SSH=(ssh -p "$PORT" -o StrictHostKeyChecking=accept-new)

if [[ ! -f "$CLUSTER/models/6M03_H41A/protein.pdb" ]]; then
  echo "Generando PDB mutados..."
  python3 "$ROOT/scripts/pack_larcad_6m03.py" --no-tarball
fi

echo "rsync $CLUSTER/ -> $HOST:$DEST/  (puerto $PORT)"
"${SSH[@]}" "$HOST" "mkdir -p $DEST"
rsync -avz --progress \
  -e "ssh -p $PORT -o StrictHostKeyChecking=accept-new" \
  --exclude '.cache/' \
  --exclude 'runs/' \
  --exclude 'logs/' \
  --exclude 'env.sh' \
  "$CLUSTER/" "$HOST:$DEST/"

echo
echo "Listo. En el nodo:"
echo "  ssh -p $PORT $HOST"
echo "  cd $DEST"
echo "  cp env.sh.example env.sh && nano env.sh"
echo "  sinfo"
echo "  module avail gromacs"
echo "  bash submit_all.sh"
