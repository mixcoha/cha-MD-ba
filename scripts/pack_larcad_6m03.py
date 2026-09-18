#!/usr/bin/env python3
"""Prepara PDB mutados y el paquete que se copia al nodo de LARCAD."""

from __future__ import annotations

import argparse
import shutil
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_version"))

from cha_md_ba.benchmark import MUTANT_MODEL_KEYS, config_for_model  # noqa: E402
from cha_md_ba.pdb_utils import clean_protein_pdb, mutate_to_alanine  # noqa: E402
from cha_md_ba.prepare import download_pdb  # noqa: E402

CLUSTER = ROOT / "cluster" / "larcad"
DEFAULT_TARBALL = ROOT / "work" / "larcad_bundle" / "cha-md-ba-6m03-larcad.tar.gz"


def prepare_models(
    cluster_dir: Path = CLUSTER,
    source_pdb: Path | None = None,
    include_wt: bool = False,
) -> dict:
    """Descarga 6M03 (o usa un PDB local), limpia y escribe los tres mutantes."""
    models_dir = cluster_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    cache = cluster_dir / ".cache"
    cache.mkdir(exist_ok=True)

    if source_pdb is None:
        raw = cache / "6M03_rcsb.pdb"
        if not raw.exists():
            ok = download_pdb("6M03", str(raw))
            if not ok:
                raise RuntimeError("No se pudo descargar 6M03 desde RCSB.")
        source_pdb = raw

    clean = cache / "6M03_clean.pdb"
    stats = clean_protein_pdb(str(source_pdb), str(clean))
    written = {"clean": str(clean), "pdb_stats": stats, "models": {}}

    keys = list(MUTANT_MODEL_KEYS)
    if include_wt:
        keys = ["wt"] + keys

    for key in keys:
        config = config_for_model(key)
        dest_dir = models_dir / config.run_id
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / "protein.pdb"
        if config.mutations:
            mut_stats = mutate_to_alanine(str(clean), str(dest), config.mutations)
            written["models"][config.run_id] = mut_stats
        else:
            shutil.copy2(clean, dest)
            written["models"][config.run_id] = {"mutations": [], "output": str(dest)}
    return written


def make_tarball(cluster_dir: Path, tarball: Path) -> Path:
    """Empaqueta cluster/larcad (sin corridas ni caché) para scp/rsync."""
    tarball.parent.mkdir(parents=True, exist_ok=True)
    exclude = {".cache", "runs", "logs", "env.sh"}
    with tarfile.open(tarball, "w:gz") as tar:
        for path in sorted(cluster_dir.rglob("*")):
            if not path.is_file():
                continue
            rel = path.relative_to(cluster_dir)
            if rel.parts and rel.parts[0] in exclude:
                continue
            tar.add(path, arcname=f"cha-md-ba-6m03-larcad/{rel.as_posix()}")
    return tarball


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Genera PDB mutados y un tar.gz para subir a LARCAD."
    )
    parser.add_argument("--source-pdb", type=Path, default=None, help="PDB 6M03 local (opcional)")
    parser.add_argument("--cluster-dir", type=Path, default=CLUSTER)
    parser.add_argument("--tarball", type=Path, default=DEFAULT_TARBALL)
    parser.add_argument("--include-wt", action="store_true", help="Incluye también el silvestre")
    parser.add_argument("--no-tarball", action="store_true")
    args = parser.parse_args(argv)

    written = prepare_models(args.cluster_dir, args.source_pdb, include_wt=args.include_wt)
    print("Modelos:")
    for name, info in written["models"].items():
        print(f"  {name}: {info.get('mutations') or 'silvestre'}")
    if not args.no_tarball:
        tar = make_tarball(args.cluster_dir, args.tarball)
        print(f"Paquete: {tar} ({tar.stat().st_size} bytes)")
        print("Súbelo con: bash scripts/upload_larcad.sh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
