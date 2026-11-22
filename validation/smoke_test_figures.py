#!/usr/bin/env python3
"""Render Figures 2-4 and verify output formats and PNG resolution."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"

for number in (2, 3, 4):
    script = FIG / f"make_figure{number}.py"
    subprocess.run([sys.executable, str(script)], cwd=FIG, check=True)
    for ext in ("pdf", "svg", "png"):
        path = FIG / f"Figure{number}.{ext}"
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
    with Image.open(FIG / f"Figure{number}.png") as image:
        dpi = image.info.get("dpi", (0, 0))
        if min(dpi) < 599:
            raise RuntimeError(f"Figure{number}.png DPI is {dpi}, expected approximately 600")
    print(f"PASS Figure {number}: PDF, SVG and 600-dpi PNG")

