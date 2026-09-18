"""Utilidades para descargar, limpiar y mutar estructuras PDB."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

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
    "HIE",
    "HID",
    "HIP",
    "HSD",
    "HSE",
    "HSP",
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


def _resseq(line: str) -> int:
    return int(line[22:26])


def _icode(line: str) -> str:
    if len(line) < 27:
        return " "
    return line[26]


def _atom_name(line: str) -> str:
    return line[12:16].strip()


# Alias de pdb2gmx / AMBER para el aminoácido de una letra.
AMINO_ALIASES: Dict[str, Set[str]] = {
    "A": {"ALA"},
    "C": {"CYS", "CYX", "CYM"},
    "D": {"ASP", "ASH"},
    "E": {"GLU", "GLH"},
    "F": {"PHE"},
    "G": {"GLY"},
    "H": {"HIS", "HIE", "HID", "HIP", "HSD", "HSE", "HSP"},
    "I": {"ILE"},
    "K": {"LYS", "LYN"},
    "L": {"LEU"},
    "M": {"MET"},
    "N": {"ASN"},
    "P": {"PRO"},
    "Q": {"GLN"},
    "R": {"ARG"},
    "S": {"SER"},
    "T": {"THR"},
    "V": {"VAL"},
    "W": {"TRP"},
    "Y": {"TYR"},
}

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

# Átomos de alanina (esqueleto + CB). pdb2gmx -ignh reconstruye hidrógenos.
_ALA_ATOM_CORE = {"N", "C", "O", "CA", "CB", "H", "HA", "HN", "HB"}
_MUTATION_RE = re.compile(
    r"^(?:(?P<chain>[A-Za-z]):)?(?P<from_aa>[A-Z])(?P<resseq>\d+)(?P<to_aa>[A-Z])$"
)


def _atom_core(name: str) -> str:
    """Quita dígitos de un nombre PDB (HB1 → HB, 2HB → HB)."""
    return re.sub(r"\d+", "", name.upper().replace(" ", ""))


def is_alanine_keep_atom(atom_name: str) -> bool:
    """True si el átomo puede quedarse al truncar a ALA."""
    return _atom_core(atom_name) in _ALA_ATOM_CORE


@dataclass(frozen=True)
class PointMutation:
    """Mutación puntual en numeración PDB (p. ej. H41A de la cadena A)."""

    from_aa: str
    resseq: int
    to_aa: str
    chain: str = "A"
    icode: str = " "

    @property
    def label(self) -> str:
        chain = self.chain.strip()
        prefix = f"{chain}:" if chain and chain != "A" else ""
        return f"{prefix}{self.from_aa}{self.resseq}{self.to_aa}"

    def matches_residue(self, line: str) -> bool:
        chain = _chain_id(line).strip() or " "
        want = (self.chain.strip() or "A").upper()
        if chain.upper() != want:
            return False
        if _resseq(line) != self.resseq:
            return False
        return _icode(line) == self.icode

    def allowed_from_resnames(self) -> Set[str]:
        return set(AMINO_ALIASES[self.from_aa])


def parse_mutation(spec: str, default_chain: str = "A") -> PointMutation:
    """Interpreta 'H41A', 'C145A' o 'A:H41A'."""
    token = spec.strip().replace(" ", "").upper()
    match = _MUTATION_RE.match(token)
    if not match:
        raise ValueError(
            f"Mutación no reconocida: {spec!r}. Usa el formato H41A o A:C145A."
        )
    from_aa = match.group("from_aa")
    to_aa = match.group("to_aa")
    if from_aa not in ONE_TO_THREE or to_aa not in ONE_TO_THREE:
        raise ValueError(f"Aminoácido no estándar en {spec!r}.")
    if to_aa != "A":
        raise ValueError(
            f"Solo se admite truncamiento a alanina (got {spec!r}). "
            "pdb2gmx reconstruye la cadena lateral de ALA a partir de N/CA/C/O/CB."
        )
    chain = match.group("chain") or default_chain
    return PointMutation(
        from_aa=from_aa,
        resseq=int(match.group("resseq")),
        to_aa=to_aa,
        chain=chain.upper(),
    )


def parse_mutations(
    specs: Sequence[str] | str,
    default_chain: str = "A",
) -> Tuple[PointMutation, ...]:
    """Acepta 'H41A,C145A' o una secuencia de especificaciones."""
    if isinstance(specs, str):
        parts = [part for part in specs.replace(";", ",").split(",") if part.strip()]
    else:
        parts = [str(item) for item in specs]
    return tuple(parse_mutation(part, default_chain=default_chain) for part in parts)


def run_id_from_mutations(pdb_id: str, mutations: Sequence[PointMutation | str]) -> str:
    """Nombre de carpeta local: 6M03, 6M03_H41A, 6M03_H41A_C145A."""
    parsed: List[PointMutation] = []
    for item in mutations:
        if isinstance(item, PointMutation):
            parsed.append(item)
        else:
            parsed.extend(parse_mutations(item))
    if not parsed:
        return pdb_id
    return pdb_id + "_" + "_".join(mut.label.replace(":", "") for mut in parsed)


def mutate_to_alanine(
    input_path: str,
    output_path: str,
    mutations: Sequence[PointMutation | str],
    default_chain: str = "A",
) -> Dict[str, object]:
    """Trunca residuos a ALA conservando N, CA, C, O y CB.

    Pensado para knockout de cadena lateral (H41A, C145A) antes de pdb2gmx
    con ``-ignh``. Comprueba que el aminoácido original coincide.
    """
    muts = tuple(
        item if isinstance(item, PointMutation) else parse_mutation(item, default_chain)
        for item in mutations
    )
    if not muts:
        raise ValueError("Se necesita al menos una mutación.")

    in_path = Path(input_path)
    out_path = Path(output_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    found: Dict[str, str] = {}
    atoms_kept = 0
    atoms_dropped = 0
    lines_out: List[str] = []

    remarks = [
        "REMARK   4 CHA-MD-BA alanine truncation\n",
        "REMARK   4 Mutations: " + ", ".join(m.label for m in muts) + "\n",
    ]

    with in_path.open() as handle:
        source_lines = handle.readlines()

    header_done = False
    for line in source_lines:
        record = line[:6].strip()
        if not header_done and record in {"ATOM", "HETATM", "TER"}:
            lines_out.extend(remarks)
            header_done = True

        if record not in {"ATOM", "HETATM", "TER"}:
            lines_out.append(line if line.endswith("\n") else line + "\n")
            continue

        matched = next((mut for mut in muts if mut.matches_residue(line)), None)
        if matched is None:
            lines_out.append(line if line.endswith("\n") else line + "\n")
            continue

        res = _resname(line)
        allowed = matched.allowed_from_resnames() | {ONE_TO_THREE[matched.to_aa]}
        if res not in allowed:
            raise ValueError(
                f"En {matched.label} se esperaba "
                f"{'/'.join(sorted(matched.allowed_from_resnames()))}, "
                f"pero el PDB tiene {res!r} (cadena {_chain_id(line)!r}, "
                f"residuo {matched.resseq})."
            )
        found[matched.label] = res

        if record == "TER":
            line = f"{line[:17]}{ONE_TO_THREE[matched.to_aa]:<3}{line[20:]}"
            lines_out.append(line if line.endswith("\n") else line + "\n")
            continue

        if matched.to_aa == "A" and not is_alanine_keep_atom(_atom_name(line)):
            atoms_dropped += 1
            continue

        line = f"{line[:17]}{ONE_TO_THREE[matched.to_aa]:<3}{line[20:]}"
        atoms_kept += 1
        lines_out.append(line if line.endswith("\n") else line + "\n")

    missing = [mut.label for mut in muts if mut.label not in found]
    if missing:
        raise ValueError(
            "No se encontraron en el PDB los residuos: " + ", ".join(missing)
        )

    if not any(line.startswith("END") for line in lines_out):
        lines_out.append("END\n")

    out_path.write_text("".join(lines_out))
    return {
        "input": str(in_path),
        "output": str(out_path),
        "mutations": [mut.label for mut in muts],
        "found_resnames": found,
        "atoms_kept": atoms_kept,
        "atoms_dropped": atoms_dropped,
    }


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
