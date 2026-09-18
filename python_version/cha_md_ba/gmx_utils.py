"""Utilidades para localizar y ejecutar GROMACS."""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Optional


def find_gmx(preferred: Optional[str] = None) -> str:
    """Devuelve el ejecutable de GROMACS disponible.

    Orden de resolución:
    1. ``preferred`` si es un ejecutable existente
    2. ``CHA_MD_BA_GMXBIN`` (archivo o directorio)
    3. ``gmx`` o ``gmx_mpi`` en PATH
    4. ``gmx_mpi`` como valor por defecto (compatibilidad con tests)
    """
    if preferred:
        if os.path.isfile(preferred) or shutil.which(preferred):
            return preferred

    env_bin = os.environ.get("CHA_MD_BA_GMXBIN")
    if env_bin:
        if os.path.isdir(env_bin):
            for name in ("gmx", "gmx_mpi"):
                candidate = os.path.join(env_bin, name)
                if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                    return candidate
        elif os.path.isfile(env_bin) or shutil.which(env_bin):
            return env_bin

    for name in ("gmx", "gmx_mpi"):
        found = shutil.which(name)
        if found:
            return found

    return preferred or "gmx_mpi"


def gmx_is_available(gmx: Optional[str] = None) -> bool:
    """Indica si GROMACS está disponible en el entorno."""
    resolved = find_gmx(gmx)
    return bool(shutil.which(resolved) or (os.path.isfile(resolved) and os.access(resolved, os.X_OK)))


def detect_gpu_ids(value: Optional[str] = "auto") -> Optional[str]:
    """Resuelve IDs de GPU: ``auto`` usa la GPU 0 si hay NVIDIA; ``none``/``cpu`` fuerza CPU."""
    if value is None:
        value = "auto"
    normalized = str(value).strip().lower()
    if normalized in {"none", "cpu", ""}:
        return None
    if normalized != "auto":
        return str(value).strip()

    nvidia = shutil.which("nvidia-smi")
    if not nvidia:
        return None
    try:
        completed = subprocess.run(
            [nvidia, "-L"],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode == 0 and "GPU" in completed.stdout:
        return "0"
    return None
