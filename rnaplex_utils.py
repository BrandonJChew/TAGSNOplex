#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional


@dataclass
class RNAplexHit:
    target_start: int
    target_end: int
    query_start: int
    query_end: int
    duplex_dbn: str
    target_region: str
    query_region: str
    delta_g: float

    @property
    def duplex_length(self) -> int:
        return len(self.target_region)


@dataclass
class SnoRNA:
    name: str
    sequence: str
    dbn: Optional[str] = None


def normalize_rna_sequence(seq: str) -> str:
    """Convert DNA-style sequence to RNA and strip whitespace."""
    return seq.upper().replace("T", "U").replace(" ", "").replace("\n", "")


def write_fasta(path: Path, name: str, sequence: str) -> None:
    """Write a single-sequence FASTA file."""
    path.write_text(f">{name}\n{sequence}\n")


def parse_rnaplex_line(line: str, target_seq: str, query_seq: str) -> Optional[RNAplexHit]:
    """
    Parse a single RNAplex output line into an RNAplexHit.

    Expected information includes:
    - duplex dot-bracket string
    - two coordinate pairs
    - free energy in parentheses or brackets
    """
    line = line.strip()
    if not line or line.startswith(">"):
        return None

    parts = line.split()
    if not parts:
        return None

    duplex_dbn = parts[0]

    coord_pairs = re.findall(r"(\d+),(\d+)", line)
    if len(coord_pairs) < 2:
        return None

    target_start, target_end = map(int, coord_pairs[0])
    query_start, query_end = map(int, coord_pairs[1])

    energy_match = re.findall(r"[\(\[]\s*([-+]?\d+(?:\.\d+)?)\s*[\)\]]", line)
    delta_g = float(energy_match[-1]) if energy_match else 0.0

    target_region = (
        target_seq[target_start - 1 : target_end]
        if target_start > 0 and target_end > 0
        else ""
    )
    query_region = (
        query_seq[query_start - 1 : query_end]
        if query_start > 0 and query_end > 0
        else ""
    )

    return RNAplexHit(
        target_start=target_start,
        target_end=target_end,
        query_start=query_start,
        query_end=query_end,
        duplex_dbn=duplex_dbn,
        target_region=target_region,
        query_region=query_region,
        delta_g=delta_g,
    )


def run_rnaplex(
    target_seq: str,
    query_seq: str,
    rnaplex_executable: str = "RNAplex",
    max_interaction_length: int = 40,
    accessibility_length: int = 40,
    energy_cutoff: Optional[float] = None,
) -> List[RNAplexHit]:
    """
    Run RNAplex for one target/query pair and return all parsed hits.

    Parameters
    ----------
    target_seq : str
        Target RNA sequence.
    query_seq : str
        Query RNA sequence.
    rnaplex_executable : str
        Path to RNAplex binary.
    max_interaction_length : int
        RNAplex -l value.
    accessibility_length : int
        RNAplex -c value.
    energy_cutoff : Optional[float]
        RNAplex -e value. If set, returns all hits passing that threshold.
    """
    target_seq = normalize_rna_sequence(target_seq)
    query_seq = normalize_rna_sequence(query_seq)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        target_fa = tmpdir_path / "target.fa"
        query_fa = tmpdir_path / "query.fa"

        write_fasta(target_fa, "target", target_seq)
        write_fasta(query_fa, "query", query_seq)

        cmd = [
            rnaplex_executable,
            "-q",
            str(query_fa),
            "-t",
            str(target_fa),
            "-l",
            str(max_interaction_length),
            "-c",
            str(accessibility_length),
        ]

        if energy_cutoff is not None:
            cmd.extend(["-e", str(energy_cutoff)])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"RNAplex failed with exit code {result.returncode}:\n{result.stderr}"
            )

        hits: List[RNAplexHit] = []
        for line in result.stdout.splitlines():
            parsed = parse_rnaplex_line(line, target_seq=target_seq, query_seq=query_seq)
            if parsed is not None:
                hits.append(parsed)

        hits.sort(key=lambda hit: hit.delta_g)  # most negative first
        return hits


def best_rnaplex_hit(
    target_seq: str,
    query_seq: str,
    **kwargs,
) -> Optional[RNAplexHit]:
    """Return the single strongest RNAplex hit, or None if no hit is found."""
    hits = run_rnaplex(target_seq=target_seq, query_seq=query_seq, **kwargs)
    return hits[0] if hits else None

