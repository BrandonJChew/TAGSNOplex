#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def normalize_transcript_id(tx_id: str) -> str:
    """Remove transcript version suffix, e.g. ENST00000335137.4 -> ENST00000335137."""
    return str(tx_id).strip().split(".")[0]


def load_transcript_ids(csv_path: str | Path, transcript_id_column: str) -> set[str]:
    """Load transcript IDs from Supplementary_Table_1.csv."""
    df = pd.read_csv(csv_path)

    if transcript_id_column not in df.columns:
        raise ValueError(
            f"Column '{transcript_id_column}' not found in {csv_path}. "
            f"Available columns: {list(df.columns)}"
        )

    ids = {
        normalize_transcript_id(tx_id)
        for tx_id in df[transcript_id_column].dropna().astype(str)
        if tx_id.strip()
    }

    if not ids:
        raise ValueError("No transcript IDs were found in the input table.")

    return ids


def get_header_transcript_id(header_line: str) -> str:
    """
    Extract transcript ID from a FASTA header.

    Supports:
    - GENCODE-style:
      >ENST00000335137.4|ENSG...
    - MANE-style:
      >ENST:ENST00000372285.8 | Gene:CD40 | ...
    """
    raw = header_line[1:].strip()

    if raw.startswith("ENST:"):
        tx = raw.split()[0].replace("ENST:", "")
        return normalize_transcript_id(tx)

    first_field = raw.split("|")[0]
    return normalize_transcript_id(first_field)


def extract_transcripts_from_fasta(
    transcript_ids: set[str],
    input_fasta: str | Path,
    output_fasta: str | Path,
) -> tuple[int, int]:
    """
    Extract matching transcript records from a transcript FASTA.

    Returns:
        (number_found, number_requested)
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    found_ids: set[str] = set()
    write_current = False

    with input_fasta.open("r") as infile, output_fasta.open("w") as outfile:
        for line in infile:
            if line.startswith(">"):
                tx_id = get_header_transcript_id(line)
                write_current = tx_id in transcript_ids
                if write_current:
                    found_ids.add(tx_id)
                    outfile.write(line)
            else:
                if write_current:
                    outfile.write(line)

    return len(found_ids), len(transcript_ids)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract transcript FASTA entries listed in Supplementary_Table_1.csv."
    )
    parser.add_argument(
        "--table",
        required=True,
        help="Path to Supplementary_Table_1.csv",
    )
    parser.add_argument(
        "--transcript-column",
        required=True,
        help="Column name in Supplementary_Table_1.csv containing transcript IDs",
    )
    parser.add_argument(
        "--input-fasta",
        required=True,
        help="Path to full transcript FASTA downloaded by the user",
    )
    parser.add_argument(
        "--output-fasta",
        required=True,
        help="Path to output subset FASTA",
    )

    args = parser.parse_args()

    transcript_ids = load_transcript_ids(
        csv_path=args.table,
        transcript_id_column=args.transcript_column,
    )

    found, requested = extract_transcripts_from_fasta(
        transcript_ids=transcript_ids,
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
    )

    print(f"[INFO] Requested transcripts: {requested}")
    print(f"[INFO] Found in FASTA:        {found}")
    print(f"[INFO] Missing from FASTA:   {requested - found}")
    print(f"[INFO] Wrote subset FASTA:   {args.output_fasta}")


if __name__ == "__main__":
    main()

