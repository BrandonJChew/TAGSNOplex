#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

import random
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


def build_highlight_string(
    regions: Sequence[Tuple[int, int]],
    colors: Optional[Sequence[str]] = None,
    seed: int = 1,
) -> str:
    """
    Build a VARNA highlightRegion string.

    Example:
        14-27:fill=#decc95;114-123:fill=#aacce3
    """
    if not regions:
        return ""

    if colors is None:
        rng = random.Random(seed)
        colors = [
            f"#{rng.randint(0,255):02X}{rng.randint(0,255):02X}{rng.randint(0,255):02X}"
            for _ in regions
        ]

    if len(colors) != len(regions):
        raise ValueError("Length of colors must match length of regions.")

    return ";".join(
        f"{start}-{end}:fill={color}"
        for (start, end), color in zip(regions, colors)
    )


def render_varna_structure(
    sequence: str,
    dbn: str,
    output_png: Path,
    varna_jar: str,
    title: str = "",
    rotation: int = 0,
    resolution: int = 3,
    highlight_regions: Optional[Sequence[Tuple[int, int]]] = None,
    highlight_colors: Optional[Sequence[str]] = None,
) -> Path:
    """
    Render an RNA secondary structure using VARNA.
    """
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    highlight_string = build_highlight_string(
        regions=highlight_regions or [],
        colors=highlight_colors,
    )

    cmd = [
        "java",
        "-cp",
        varna_jar,
        "fr.orsay.lri.varna.applications.VARNAcmd",
        "-sequenceDBN",
        sequence,
        "-structureDBN",
        dbn,
        "-algorithm",
        "radiate",
        "-baseInner",
        "#334455",
        "-baseName",
        "#334455",
        "-baseOutline",
        "#334455",
        "-bp",
        "#334455",
        "-bpStyle",
        "line",
        "-flat",
        "true",
        "-resolution",
        str(resolution),
        "-rotation",
        str(rotation),
        "-title",
        title,
        "-titleSize",
        "10",
        "-o",
        str(output_png),
    ]

    if highlight_string:
        cmd.extend(["-highlightRegion", highlight_string])

    result = subprocess.run(cmd, capture_output=True, text=True, check=False)

    if result.returncode != 0:
        raise RuntimeError(
            f"VARNA failed with exit code {result.returncode}:\n{result.stderr}"
        )

    if not output_png.exists():
        raise FileNotFoundError(f"Expected VARNA output was not created: {output_png}")

    return output_png

