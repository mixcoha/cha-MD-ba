"""Benchmark de dinámica molecular para 6M03 en condiciones estándar."""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional

from rich.console import Console

from .gmx_utils import find_gmx, gmx_is_available
from .minimize import EnergyMinimizer
from .npt import NPTEquilibrator
from .nvt import NVTEquilibrator
from .pdb_utils import clean_protein_pdb
from .prepare import MDSystemPreparator, download_pdb

console = Console()

DEFAULT_OUTPUT_DIR = "work"
DEFAULT_DATA_DIR = "work/data"


@dataclass
class Benchmark6M03Config:
    """Condiciones del benchmark 6M03 en agua con NaCl a temperatura fisiológica."""

    pdb_id: str = "6M03"
    description: str = (
        "SARS-CoV-2 Mpro apo (6M03) disuelta en agua TIP3P con NaCl 0.5 M a 310 K"
    )
    temperature: float = 310.0
    ion_concentration: float = 0.5
    pressure: float = 1.0
    forcefield: str = "amber99sb-ildn"
    water_model: str = "tip3p"
    box_type: str = "dodecahedron"
    box_distance: float = 1.2
    dt: float = 0.002
    minimization_nsteps: int = 50000
    nvt_nsteps: int = 50000
    npt_nsteps: int = 50000
    production_nsteps: int = 5_000_000  # 10 ns
    nvt_force_constants: tuple = (1000, 800, 600, 400, 200)
    ion_positive: str = "NA"
    ion_negative: str = "CL"


def parse_system_composition(topol_path: Path) -> Dict[str, int]:
    """Lee el recuento de moléculas del bloque [ molecules ] de topol.top."""
    text = Path(topol_path).read_text()
    match = re.search(r"\[ molecules \]\s*(.*)$", text, flags=re.DOTALL | re.IGNORECASE)
    if not match:
        return {}
    composition: Dict[str, int] = {}
    for line in match.group(1).splitlines():
        stripped = line.split(";")[0].strip()
        if not stripped:
            continue
        if stripped.startswith("["):
            break
        parts = stripped.split()
        if len(parts) >= 2 and parts[-1].lstrip("-").isdigit():
            composition[parts[0]] = int(parts[-1])
    return composition


def write_protocol_mdps(output_dir: Path, config: Benchmark6M03Config) -> Dict[str, Path]:
    """Escribe los .mdp del protocolo a 310 K para NVT, NPT y producción."""
    output_dir.mkdir(parents=True, exist_ok=True)
    dummy_gro = output_dir / ".placeholder.gro"
    dummy_top = output_dir / ".placeholder.top"
    dummy_gro.write_text("")
    dummy_top.write_text("")

    nvt = NVTEquilibrator(
        str(dummy_gro),
        str(dummy_top),
        temperature=config.temperature,
        nsteps=config.nvt_nsteps,
    )
    npt = NPTEquilibrator(
        str(dummy_gro),
        str(dummy_top),
        temperature=config.temperature,
        nsteps=config.npt_nsteps,
        pressure=config.pressure,
    )

    nvt_mdp = nvt.create_mdp_file(
        output_dir / "nvt.mdp",
        title=f"NVT {config.pdb_id} {config.temperature} K",
        use_posres=True,
        gen_vel=True,
        continuation=False,
    )
    npt_mdp = npt.create_mdp_file(
        output_dir / "npt.mdp",
        title=f"NPT {config.pdb_id} {config.temperature} K",
        is_production=False,
        gen_vel=False,
        continuation=True,
    )
    md_mdp = npt.create_mdp_file(
        output_dir / "md.mdp",
        title=f"Production {config.pdb_id} {config.temperature} K {config.ion_concentration} M NaCl",
        is_production=True,
        nsteps=config.production_nsteps,
        gen_vel=False,
        continuation=True,
    )

    dummy_gro.unlink(missing_ok=True)
    dummy_top.unlink(missing_ok=True)
    return {"nvt": nvt_mdp, "npt": npt_mdp, "md": md_mdp}


def run_benchmark(
    config: Optional[Benchmark6M03Config] = None,
    output_dir: str = DEFAULT_OUTPUT_DIR,
    data_dir: str = DEFAULT_DATA_DIR,
    stages: Optional[List[str]] = None,
    gpu_ids: Optional[str] = None,
    gmx: Optional[str] = None,
) -> Dict[str, object]:
    """Ejecuta el benchmark 6M03: descarga, limpia, solvata con NaCl 0.5 M y 310 K."""
    config = config or Benchmark6M03Config()
    stages = stages or ["download", "clean", "prepare", "mdps", "minimize"]
    gmx_cmd = gmx or find_gmx()
    os.environ.setdefault("GMX_MAXBACKUP", "0")
    base_output = Path(output_dir)
    data_path = Path(data_dir)
    data_path.mkdir(parents=True, exist_ok=True)

    raw_pdb = data_path / f"{config.pdb_id}_raw.pdb"
    clean_pdb = data_path / f"{config.pdb_id}.pdb"
    result: Dict[str, object] = {
        "config": asdict(config),
        "gmx": gmx_cmd,
        "stages": stages,
    }

    if "download" in stages:
        console.print(f"[cyan]Descargando {config.pdb_id} desde RCSB...[/cyan]")
        raw_download = data_path / f"{config.pdb_id}_rcsb.pdb"
        if not raw_download.exists():
            ok = download_pdb(config.pdb_id, str(raw_download))
            if not ok:
                raise RuntimeError(f"No se pudo descargar {config.pdb_id}")
        raw_pdb = raw_download
        result["raw_pdb"] = str(raw_pdb)

    if "clean" in stages:
        console.print("[cyan]Limpiando PDB (solo proteína, sin aguas cristalográficas)...[/cyan]")
        source = raw_pdb if raw_pdb.exists() else Path(str(result.get("raw_pdb", raw_pdb)))
        stats = clean_protein_pdb(str(source), str(clean_pdb))
        result["clean_pdb"] = str(clean_pdb)
        result["pdb_stats"] = stats
        console.print(
            f"[green]Átomos de proteína: {stats['atom_kept']} "
            f"({stats['n_residues']} residuos, cadenas {stats['chains']})[/green]"
        )

    protocol_dir = base_output / config.pdb_id / "protocol"
    if "mdps" in stages:
        console.print(f"[cyan]Escribiendo MDP del protocolo a {config.temperature} K...[/cyan]")
        mdps = write_protocol_mdps(protocol_dir, config)
        result["mdps"] = {key: str(path) for key, path in mdps.items()}

    prepared = None
    if "prepare" in stages:
        if not gmx_is_available(gmx_cmd):
            raise RuntimeError(
                "GROMACS no está disponible. Instálalo o define CHA_MD_BA_GMXBIN. "
                "Los MDP del protocolo ya pueden generarse con --stages mdps."
            )
        pdb_for_prep = clean_pdb if clean_pdb.exists() else raw_pdb
        console.print(
            f"[cyan]Preparando sistema: {config.forcefield}, {config.water_model}, "
            f"NaCl {config.ion_concentration} M...[/cyan]"
        )
        preparator = MDSystemPreparator(
            pdb_path=str(pdb_for_prep),
            forcefield=config.forcefield,
            water_model=config.water_model,
            ion_concentration=config.ion_concentration,
            gmx=gmx_cmd,
            ion_positive=config.ion_positive,
            ion_negative=config.ion_negative,
        )
        prepared = preparator.prepare_system(
            output_dir=str(base_output),
            box_type=config.box_type,
            box_size=config.box_distance,
            ions=True,
            minimize=False,
        )
        result["preparation"] = {key: str(value) if value is not None else None for key, value in prepared.items()}
        composition = parse_system_composition(prepared["topology"])
        result["composition"] = composition
        console.print(f"[green]Composición: {composition}[/green]")

    if "minimize" in stages:
        if prepared is None:
            raise RuntimeError("La minimización requiere la etapa 'prepare'.")
        console.print("[cyan]Minimización de energía...[/cyan]")
        minimizer = EnergyMinimizer(
            input_gro=str(prepared["ions"]),
            topol_top=str(prepared["topology"]),
            gmx=gmx_cmd,
            nsteps=config.minimization_nsteps,
        )
        min_files = minimizer.minimize(prepared["minimization_dir"], gpu_ids=gpu_ids)
        result["minimization"] = {key: str(path) for key, path in min_files.items()}

    if "nvt" in stages:
        if prepared is None:
            raise RuntimeError("NVT requiere la etapa 'prepare'.")
        min_gro = Path(str(result.get("minimization", {}).get("final") or prepared["minimization_dir"] / "minimized.gro"))
        if not min_gro.exists():
            raise RuntimeError("No hay estructura minimizada; ejecuta la etapa 'minimize'.")
        console.print(f"[cyan]Equilibración NVT a {config.temperature} K...[/cyan]")
        nvt = NVTEquilibrator(
            input_gro=str(min_gro),
            topol_top=str(prepared["topology"]),
            temperature=config.temperature,
            nsteps=config.nvt_nsteps,
            gmx=gmx_cmd,
        )
        nvt_files = nvt.equilibrate(
            prepared["nvt_dir"],
            force_constants=list(config.nvt_force_constants),
            gpu_ids=gpu_ids,
        )
        result["nvt"] = nvt_files

    if "npt" in stages:
        nvt_result = result.get("nvt")
        if not nvt_result:
            raise RuntimeError("NPT requiere la etapa 'nvt'.")
        last_fc = list(config.nvt_force_constants)[-1]
        last_gro = nvt_result[last_fc]["gro"]
        console.print(f"[cyan]Equilibración NPT a {config.temperature} K y {config.pressure} bar...[/cyan]")
        npt = NPTEquilibrator(
            input_gro=str(last_gro),
            topol_top=str(prepared["topology"]),
            temperature=config.temperature,
            nsteps=config.npt_nsteps,
            pressure=config.pressure,
            gmx=gmx_cmd,
        )
        npt_files = npt.equilibrate(prepared["npt_dir"], gpu_ids=gpu_ids)
        result["npt"] = {key: str(path) for key, path in npt_files.items()}

    report_dir = base_output / config.pdb_id
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "benchmark_report.json"
    serializable = _jsonable(result)
    report_path.write_text(json.dumps(serializable, indent=2, ensure_ascii=False))
    result["report"] = str(report_path)
    console.print(f"[bold green]Benchmark 6M03 listo. Reporte: {report_path}[/bold green]")
    return result


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark de 6M03 (Mpro apo de SARS-CoV-2) en agua TIP3P "
            "con NaCl 0.5 M a 310 K."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directorio local de simulaciones (gitignored; por defecto work/)",
    )
    parser.add_argument(
        "--data-dir",
        default=DEFAULT_DATA_DIR,
        help="Directorio local para PDB descargados (gitignored; por defecto work/data/)",
    )
    parser.add_argument(
        "--stages",
        default="download,clean,prepare,mdps,minimize",
        help="Etapas separadas por coma: download,clean,prepare,mdps,minimize,nvt,npt",
    )
    parser.add_argument("--gpu-ids", default=None, help='IDs de GPU, p. ej. "0". Vacío = solo CPU')
    parser.add_argument("--gmx", default=None, help="Ejecutable de GROMACS (gmx / gmx_mpi)")
    parser.add_argument("--em-nsteps", type=int, default=None, help="Pasos máximos de minimización")
    parser.add_argument("--box-distance", type=float, default=None, help="Distancia proteína-caja en nm")
    return parser


def run_6m03_benchmark(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = Benchmark6M03Config()
    if args.em_nsteps:
        config.minimization_nsteps = args.em_nsteps
    if args.box_distance:
        config.box_distance = args.box_distance
    stages = [item.strip() for item in args.stages.split(",") if item.strip()]
    run_benchmark(
        config=config,
        output_dir=args.output_dir,
        data_dir=args.data_dir,
        stages=stages,
        gpu_ids=args.gpu_ids,
        gmx=args.gmx,
    )
    return 0
