"""Utilidades para localizar y ejecutar GROMACS."""

from __future__ import annotations

import os
import shutil
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
