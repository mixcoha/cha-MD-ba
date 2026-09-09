"""Utilidades para descargar y limpiar estructuras PDB."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, Optional, Set

STANDARD_AMINOACIDS: Set[str] = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASX",
    "GLX",
    "SEC",
    "PYL",
}


def _resname(line: str) -> str:
    return line[17:20].strip().upper()


def _altloc(line: str) -> str:
    if len(line) < 17:
        return " "
    return line[16]


def _chain_id(line: str) -> str:
    if len(line) < 22:
        return " "
    return line[21]


def clean_protein_pdb(
    input_path: str,
    output_path: str,
    chains: Optional[Iterable[str]] = None,
    keep_hetatm: bool = False,
    histidine: str = "HIE",
) -> Dict[str, object]:
    """Conserva átomos de proteína y descarta cristalizantes, iones y ligandos.

    Se mantienen registros ATOM (y TER asociados) de aminoácidos estándar.
    Las localizaciones alternativas distintas de blanco o A se eliminan.
    Las aguas cristalográficas (HOH) y el resto de HETATM se descartan por
    defecto: el solvente se reconstruye al solvatar.
    Los residuos HIS se renombran a HIE por defecto para evitar la selección
    interactiva de pdb2gmx.

    Returns:
        Resumen con recuentos de líneas leídas y escritas.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    allowed_chains = {c.strip().upper() for c in chains} if chains else None
    stats = {
        "input": str(in_path),
        "output": str(out_path),
        "atom_kept": 0,
        "atom_skipped_altloc": 0,
        "atom_skipped_nonprotein": 0,
        "hetatm_skipped": 0,
        "hetatm_kept": 0,
        "ter_kept": 0,
        "chains": set(),
        "residues": set(),
    }

    kept_lines = []
    with in_path.open() as handle:
        for line in handle:
            record = line[:6].strip()
            if record in {"HEADER", "TITLE", "CRYST1", "END"}:
                kept_lines.append(line if line.endswith("\n") else line + "\n")
                continue

            if record == "ATOM":
                if allowed_chains and _chain_id(line).upper() not in allowed_chains:
                    continue
                alt = _altloc(line)
                if alt not in {" ", "A"}:
                    stats["atom_skipped_altloc"] += 1
                    continue
                res = _resname(line)
                if res not in STANDARD_AMINOACIDS:
                    stats["atom_skipped_nonprotein"] += 1
                    continue
                if res == "HIS" and histidine:
                    line = f"{line[:17]}{histidine:<3}{line[20:]}"
                    res = histidine
                stats["atom_kept"] += 1
                stats["chains"].add(_chain_id(line))
                stats["residues"].add((_chain_id(line), line[22:26].strip(), res))
                kept_lines.append(line if line.endswith("\n") else line + "\n")
                continue

            if record == "TER":
                if allowed_chains and _chain_id(line).upper() not in allowed_chains:
                    continue
                if _resname(line) == "HIS" and histidine:
                    line = f"{line[:17]}{histidine:<3}{line[20:]}"
                stats["ter_kept"] += 1
                kept_lines.append(line if line.endswith("\n") else line + "\n")
                continue

            if record == "HETATM":
                if keep_hetatm:
                    stats["hetatm_kept"] += 1
                    kept_lines.append(line if line.endswith("\n") else line + "\n")
                else:
                    stats["hetatm_skipped"] += 1
                continue

    if not any(line.startswith("END") for line in kept_lines):
        kept_lines.append("END\n")

    out_path.write_text("".join(kept_lines))
    stats["chains"] = sorted(stats["chains"])
    stats["n_residues"] = len(stats["residues"])
    del stats["residues"]
    return stats
