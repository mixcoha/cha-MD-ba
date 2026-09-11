#!/usr/bin/env python3
"""Lanza el benchmark de 6M03 (silvestre o mutantes H41A / C145A) a 310 K."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_version"))

from cha_md_ba.benchmark import run_6m03_benchmark  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(run_6m03_benchmark())
