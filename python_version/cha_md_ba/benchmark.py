"""Benchmark de dinámica molecular para 6M03 en condiciones estándar."""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rich.console import Console

from .gmx_utils import detect_gpu_ids, find_gmx, gmx_is_available
from .minimize import EnergyMinimizer
from .npt import NPTEquilibrator
from .nvt import NVTEquilibrator
from .pdb_utils import (
    PointMutation,
    clean_protein_pdb,
    mutate_to_alanine,
    parse_mutations,
    run_id_from_mutations,
)
from .prepare import MDSystemPreparator, download_pdb

console = Console()

DEFAULT_OUTPUT_DIR = "work"
DEFAULT_DATA_DIR = "work/data"


# Díada catalítica de Mpro (6M03, cadena A, numeración PDB).
# Modelo 1: H41A; modelo 2: C145A; modelo 3: ambas.
MPRO_CATALYTIC_MODELS = {
    "wt": {
        "run_id": "6M03",
        "mutations": (),
        "title": "silvestre (6M03)",
    },
    "1": {
        "run_id": "6M03_H41A",
        "mutations": ("H41A",),
        "title": "Modelo 1: H41A",
    },
    "2": {
        "run_id": "6M03_C145A",
        "mutations": ("C145A",),
        "title": "Modelo 2: C145A",
    },
    "3": {
        "run_id": "6M03_H41A_C145A",
        "mutations": ("H41A", "C145A"),
        "title": "Modelo 3: H41A + C145A",
    },
}
MPRO_CATALYTIC_MODELS["h41a"] = MPRO_CATALYTIC_MODELS["1"]
MPRO_CATALYTIC_MODELS["c145a"] = MPRO_CATALYTIC_MODELS["2"]
MPRO_CATALYTIC_MODELS["double"] = MPRO_CATALYTIC_MODELS["3"]
MPRO_CATALYTIC_MODELS["h41a_c145a"] = MPRO_CATALYTIC_MODELS["3"]
MUTANT_MODEL_KEYS = ("1", "2", "3")


@dataclass
class Benchmark6M03Config:
    """Condiciones del benchmark 6M03 en agua con NaCl a temperatura fisiológica."""

    pdb_id: str = "6M03"
    run_id: str = "6M03"
    mutations: tuple = ()
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

    def parsed_mutations(self) -> Tuple[PointMutation, ...]:
        if not self.mutations:
            return ()
        return parse_mutations(self.mutations)

    @property
    def mutation_labels(self) -> List[str]:
        return [mut.label for mut in self.parsed_mutations()]


def config_for_model(model: str, base: Optional[Benchmark6M03Config] = None) -> Benchmark6M03Config:
    """Devuelve la config de un modelo (wt, 1/H41A, 2/C145A, 3/doble)."""
    key = str(model).strip().lower()
    if key not in MPRO_CATALYTIC_MODELS:
        raise ValueError(
            f"Modelo desconocido: {model!r}. Usa wt, 1 (H41A), 2 (C145A), 3 (ambas) o all."
        )
    spec = MPRO_CATALYTIC_MODELS[key]
    config = Benchmark6M03Config() if base is None else Benchmark6M03Config(**{
        field: getattr(base, field)
        for field in (
            "pdb_id",
            "temperature",
            "ion_concentration",
            "pressure",
            "forcefield",
            "water_model",
            "box_type",
            "box_distance",
            "dt",
            "minimization_nsteps",
            "nvt_nsteps",
            "npt_nsteps",
            "production_nsteps",
            "nvt_force_constants",
            "ion_positive",
            "ion_negative",
        )
    })
    config.run_id = spec["run_id"]
    config.mutations = spec["mutations"]
    title = spec["title"]
    if config.mutations:
        config.description = (
            f"{title}. SARS-CoV-2 Mpro (6M03) en agua TIP3P con NaCl 0.5 M a 310 K"
        )
    return config


def resolve_model_keys(model: str) -> List[str]:
    """Expande 'all' a los tres mutantes; el resto queda como un solo modelo."""
    key = str(model).strip().lower()
    if key in {"all", "mutants", "mutantes"}:
        return list(MUTANT_MODEL_KEYS)
    if key not in MPRO_CATALYTIC_MODELS:
        raise ValueError(
            f"Modelo desconocido: {model!r}. Usa wt, 1, 2, 3 o all."
        )
    canonical = {
        "h41a": "1",
        "c145a": "2",
        "double": "3",
        "h41a_c145a": "3",
    }
    return [canonical.get(key, key)]


def system_paths(output_dir: Path, run_id: str) -> Dict[str, Path]:
    """Rutas estándar de una corrida local bajo ``output_dir/run_id``."""
    base = Path(output_dir) / run_id
    prep = base / "1_preparation"
    return {
        "base_dir": base,
        "preparation_dir": prep,
        "minimization_dir": base / "2_minimization",
        "nvt_dir": base / "3_nvt",
        "npt_dir": base / "4_npt",
        "production_dir": base / "5_production",
        "topology": prep / "topol.top",
        "ions": prep / f"{run_id}_ions.gro",
        "minimized": base / "2_minimization" / "minimized.gro",
    }


def nvt_stage_complete(nvt_dir: Path, force_constant: int) -> bool:
    """True si esa constante de POSRES ya tiene estructura y energías NVT."""
    fc_dir = Path(nvt_dir) / "posre_constante" / str(force_constant)
    has_gro = (fc_dir / "nvt.gro").exists()
    has_edr = (fc_dir / "nvt.edr").exists() or (fc_dir / "ener.edr").exists()
    return has_gro and has_edr


def nvt_output_gro(nvt_dir: Path, force_constant: int) -> Path:
    """GRO de salida de un NVT (centrado si existe ``tmp.gro``)."""
    fc_dir = Path(nvt_dir) / "posre_constante" / str(force_constant)
    centered = fc_dir / "tmp.gro"
    return centered if centered.exists() else fc_dir / "nvt.gro"


def last_completed_nvt_fc(nvt_dir: Path, force_constants: tuple) -> Optional[int]:
    """Última constante NVT completada en orden; no usa max() numérico."""
    last: Optional[int] = None
    for force_constant in force_constants:
        if nvt_stage_complete(nvt_dir, force_constant):
            last = force_constant
        else:
            break
    return last


def load_prepared_paths(output_dir: Path, run_id: str) -> Dict[str, Path]:
    """Reconstruye las rutas de un sistema ya preparado en ``work/``."""
    paths = system_paths(output_dir, run_id)
    missing = [name for name in ("topology", "ions") if not paths[name].exists()]
    if missing:
        raise RuntimeError(
            f"No hay sistema preparado en {paths['base_dir']} "
            f"(faltan: {', '.join(missing)}). Ejecuta --resume o la etapa prepare."
        )
    return paths


def initial_stages(config: Benchmark6M03Config) -> List[str]:
    """Etapas de una corrida desde cero, con mutación si el modelo la pide."""
    stages = ["download", "clean"]
    if config.mutations:
        stages.append("mutate")
    stages.extend(["prepare", "mdps", "minimize", "nvt", "npt"])
    return stages


def resume_stages(output_dir: Path, config: Benchmark6M03Config) -> List[str]:
    """Elige las etapas que faltan para continuar una corrida local."""
    paths = system_paths(output_dir, config.run_id)
    if not paths["topology"].exists() or not paths["ions"].exists():
        return initial_stages(config)
    if not paths["minimized"].exists():
        return ["mdps", "minimize", "nvt", "npt"]
    if last_completed_nvt_fc(paths["nvt_dir"], config.nvt_force_constants) != config.nvt_force_constants[-1]:
        return ["nvt", "npt"]
    if not (paths["npt_dir"] / "npt.gro").exists():
        return ["npt"]
    return []


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
        title=f"NVT {config.run_id} {config.temperature} K",
        use_posres=True,
        gen_vel=True,
        continuation=False,
    )
    npt_mdp = npt.create_mdp_file(
        output_dir / "npt.mdp",
        title=f"NPT {config.run_id} {config.temperature} K",
        is_production=False,
        gen_vel=False,
        continuation=True,
    )
    md_mdp = npt.create_mdp_file(
        output_dir / "md.mdp",
        title=f"Production {config.run_id} {config.temperature} K {config.ion_concentration} M NaCl",
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
    mutated_pdb = data_path / f"{config.run_id}.pdb"
    result: Dict[str, object] = {
        "config": asdict(config),
        "gmx": gmx_cmd,
        "stages": stages,
        "run_id": config.run_id,
        "mutations": list(config.mutations),
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

    if "mutate" in stages or (config.mutations and "prepare" in stages):
        if not config.mutations:
            raise RuntimeError("La etapa mutate requiere --model o --mutations.")
        source = clean_pdb if clean_pdb.exists() else raw_pdb
        if not source.exists():
            raise RuntimeError(
                f"No hay PDB limpio para mutar ({source}). Ejecuta download+clean o --resume."
            )
        console.print(
            f"[cyan]Mutando {', '.join(config.mutation_labels)} → {mutated_pdb.name}[/cyan]"
        )
        mut_stats = mutate_to_alanine(str(source), str(mutated_pdb), config.mutations)
        result["mutated_pdb"] = str(mutated_pdb)
        result["mutation_stats"] = mut_stats
        console.print(
            f"[green]Mutaciones aplicadas: {mut_stats['mutations']} "
            f"(átomos de cadena lateral eliminados: {mut_stats['atoms_dropped']})[/green]"
        )

    protocol_dir = base_output / config.run_id / "protocol"
    if "mdps" in stages:
        console.print(f"[cyan]Escribiendo MDP del protocolo a {config.temperature} K...[/cyan]")
        mdps = write_protocol_mdps(protocol_dir, config)
        result["mdps"] = {key: str(path) for key, path in mdps.items()}

    needs_system = any(stage in stages for stage in ("prepare", "minimize", "nvt", "npt"))
    prepared = None
    if needs_system:
        existing = system_paths(base_output, config.run_id)
        if existing["topology"].exists() and existing["ions"].exists() and "prepare" not in stages:
            prepared = load_prepared_paths(base_output, config.run_id)
            result["preparation"] = {key: str(value) for key, value in prepared.items()}
            result["composition"] = parse_system_composition(prepared["topology"])
            console.print(f"[green]Reutilizando sistema en {prepared['base_dir']}[/green]")

    if "prepare" in stages:
        existing = system_paths(base_output, config.run_id)
        if existing["topology"].exists() and existing["ions"].exists():
            console.print("[yellow]Sistema ya preparado; se reutiliza (no se vuelve a solvar).[/yellow]")
            prepared = load_prepared_paths(base_output, config.run_id)
            result["preparation"] = {key: str(value) for key, value in prepared.items()}
            result["composition"] = parse_system_composition(prepared["topology"])
            console.print(f"[green]Composición: {result['composition']}[/green]")
        else:
            if not gmx_is_available(gmx_cmd):
                raise RuntimeError(
                    "GROMACS no está disponible. Instálalo o define CHA_MD_BA_GMXBIN. "
                    "Los MDP del protocolo ya pueden generarse con --stages mdps."
                )
            if config.mutations and mutated_pdb.exists():
                pdb_for_prep = mutated_pdb
            elif clean_pdb.exists():
                pdb_for_prep = clean_pdb
            else:
                pdb_for_prep = raw_pdb
            console.print(
                f"[cyan]Preparando {config.run_id}: {config.forcefield}, {config.water_model}, "
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
            prepared = load_prepared_paths(base_output, config.run_id)
        min_gro = Path(prepared["minimization_dir"]) / "minimized.gro"
        if min_gro.exists():
            console.print("[yellow]Minimización ya completa; se reutiliza minimized.gro.[/yellow]")
            result["minimization"] = {"final": str(min_gro)}
        else:
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
            prepared = load_prepared_paths(base_output, config.run_id)
        min_gro = Path(
            str(result.get("minimization", {}).get("final") or Path(prepared["minimization_dir"]) / "minimized.gro")
        )
        if not min_gro.exists():
            raise RuntimeError("No hay estructura minimizada; ejecuta la etapa 'minimize' o --resume.")
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
        if prepared is None:
            prepared = load_prepared_paths(base_output, config.run_id)
        nvt_result = result.get("nvt")
        last_fc = list(config.nvt_force_constants)[-1]
        if not nvt_result:
            if not nvt_stage_complete(prepared["nvt_dir"], last_fc):
                raise RuntimeError("NPT requiere NVT completo (POSRES 200). Ejecuta --resume.")
            last_gro = nvt_output_gro(prepared["nvt_dir"], last_fc)
            nvt_result = {last_fc: {"gro": str(last_gro)}}
            result["nvt"] = nvt_result
        nvt_last = nvt_result.get(last_fc) or nvt_result.get(str(last_fc))
        last_gro = nvt_last["gro"]
        npt_gro = Path(prepared["npt_dir"]) / "npt.gro"
        if npt_gro.exists():
            console.print("[yellow]NPT ya completo; se reutiliza npt.gro.[/yellow]")
            result["npt"] = {"final": str(npt_gro)}
        else:
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

    report_dir = base_output / config.run_id
    report_dir.mkdir(parents=True, exist_ok=True)
    report_path = report_dir / "benchmark_report.json"
    serializable = _jsonable(result)
    report_path.write_text(json.dumps(serializable, indent=2, ensure_ascii=False))
    result["report"] = str(report_path)
    console.print(f"[bold green]Benchmark {config.run_id} listo. Reporte: {report_path}[/bold green]")
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
            "con NaCl 0.5 M a 310 K. Modelos: silvestre, H41A, C145A o ambas."
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
        "--model",
        default="wt",
        help="wt | 1 (H41A) | 2 (C145A) | 3 (H41A+C145A) | all (los tres mutantes)",
    )
    parser.add_argument(
        "--mutations",
        default=None,
        help="Mutaciones extra (p. ej. H41A o H41A,C145A). Anula --model si se indica.",
    )
    parser.add_argument(
        "--stages",
        default="download,clean,prepare,mdps,minimize",
        help="Etapas: download,clean,mutate,prepare,mdps,minimize,nvt,npt",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Continúa la corrida local en work/ desde la etapa que falte (NVT/NPT inclusive)",
    )
    parser.add_argument(
        "--gpu-ids",
        default="auto",
        help='IDs de GPU ("0", "01"), "auto" (GPU 0 si hay NVIDIA) o "none" para CPU',
    )
    parser.add_argument("--gmx", default=None, help="Ejecutable de GROMACS (gmx / gmx_mpi)")
    parser.add_argument("--em-nsteps", type=int, default=None, help="Pasos máximos de minimización")
    parser.add_argument("--box-distance", type=float, default=None, help="Distancia proteína-caja en nm")
    return parser


def _apply_cli_overrides(config: Benchmark6M03Config, args) -> Benchmark6M03Config:
    if args.em_nsteps:
        config.minimization_nsteps = args.em_nsteps
    if args.box_distance:
        config.box_distance = args.box_distance
    if args.mutations:
        parsed = parse_mutations(args.mutations)
        config.mutations = tuple(mut.label for mut in parsed)
        config.run_id = run_id_from_mutations(config.pdb_id, config.mutations)
        labels = ", ".join(config.mutation_labels)
        config.description = (
            f"{labels}. SARS-CoV-2 Mpro (6M03) en agua TIP3P con NaCl 0.5 M a 310 K"
        )
    return config


def run_6m03_benchmark(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    gpu_ids = detect_gpu_ids(args.gpu_ids)
    model_keys = ["custom"] if args.mutations else resolve_model_keys(args.model)

    for model_key in model_keys:
        if args.mutations:
            config = _apply_cli_overrides(Benchmark6M03Config(), args)
        else:
            config = _apply_cli_overrides(config_for_model(model_key), args)
        if args.resume:
            stages = resume_stages(Path(args.output_dir), config)
            if not stages:
                console.print(
                    f"[bold green]{config.run_id} ya está completo (NVT + NPT).[/bold green]"
                )
                continue
            console.print(f"[cyan]{config.run_id}: reanudando {', '.join(stages)}[/cyan]")
        else:
            stages = [item.strip() for item in args.stages.split(",") if item.strip()]
            if config.mutations and "mutate" not in stages:
                if "clean" in stages:
                    idx = stages.index("clean") + 1
                    stages.insert(idx, "mutate")
                elif "prepare" in stages:
                    stages.insert(stages.index("prepare"), "mutate")
            console.print(f"[cyan]{config.run_id}: {', '.join(config.mutation_labels) or 'silvestre'}[/cyan]")
        if gpu_ids:
            console.print(f"[cyan]GPU: {gpu_ids}[/cyan]")
        else:
            console.print("[yellow]Sin GPU: mdrun en CPU (más lento).[/yellow]")
        run_benchmark(
            config=config,
            output_dir=args.output_dir,
            data_dir=args.data_dir,
            stages=stages,
            gpu_ids=gpu_ids,
            gmx=args.gmx,
        )
    return 0
