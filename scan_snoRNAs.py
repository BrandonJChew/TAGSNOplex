#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import pandas as pd

from rnaplex_utils import SnoRNA, normalize_rna_sequence, run_rnaplex


def load_snornas_from_csv(csv_path: str | Path) -> List[SnoRNA]:
    """
    Load snoRNAs from a CSV file with columns:
    - gene_name
    - sequence
    - Dbn (optional but expected in your current files)
    """
    df = pd.read_csv(csv_path)

    required_columns = {"gene_name", "sequence"}
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in input CSV: {sorted(missing)}")

    snornas: List[SnoRNA] = []
    for _, row in df.dropna(subset=["gene_name", "sequence"]).iterrows():
        snornas.append(
            SnoRNA(
                name=str(row["gene_name"]),
                sequence=normalize_rna_sequence(str(row["sequence"])),
                dbn=str(row["Dbn"]) if "Dbn" in df.columns and pd.notna(row.get("Dbn")) else None,
            )
        )
    return snornas


def scan_snornas_against_target(
    snornas: List[SnoRNA],
    target_rna: str,
    rnaplex_executable: str = "RNAplex",
    max_interaction_length: int = 40,
    accessibility_length: int = 40,
    energy_cutoff: float = -18.0,
    max_hits_per_snorna: int = 5,
) -> pd.DataFrame:
    """
    Scan a list of snoRNAs against one target RNA and return a tidy DataFrame.
    """
    target_rna = normalize_rna_sequence(target_rna)
    records = []

    for sno in snornas:
        try:
            hits = run_rnaplex(
                target_seq=target_rna,
                query_seq=sno.sequence,
                rnaplex_executable=rnaplex_executable,
                max_interaction_length=max_interaction_length,
                accessibility_length=accessibility_length,
                energy_cutoff=energy_cutoff,
            )

            seen_spans = set()
            kept_hits = 0

            for hit in hits:
                span_key = (
                    sno.name,
                    hit.target_start,
                    hit.target_end,
                    hit.query_start,
                    hit.query_end,
                )
                if span_key in seen_spans:
                    continue

                seen_spans.add(span_key)
                kept_hits += 1

                records.append(
                    {
                        "snoRNA": sno.name,
                        "deltaG": hit.delta_g,
                        "duplex_dbn": hit.duplex_dbn,
                        "target_region_seq": hit.target_region,
                        "snoRNA_region_seq": hit.query_region,
                        "target_start": hit.target_start,
                        "target_end": hit.target_end,
                        "snoRNA_start": hit.query_start,
                        "snoRNA_end": hit.query_end,
                        "target_pos": f"{hit.target_start}-{hit.target_end}",
                        "snoRNA_pos": f"{hit.query_start}-{hit.query_end}",
                        "duplex_length": hit.duplex_length,
                    }
                )

                if kept_hits >= max_hits_per_snorna:
                    break

        except Exception as exc:
            print(f"[WARNING] Failed to process {sno.name}: {exc}")

    df = pd.DataFrame(records)
    if df.empty:
        return df

    df = df.sort_values(["snoRNA", "deltaG"], ascending=[True, True]).reset_index(drop=True)
    df["hit_rank"] = df.groupby("snoRNA").cumcount() + 1
    return df


def filter_ranked_hits(
    df: pd.DataFrame,
    min_delta_g: float = -18.0,
    max_duplex_length: int = 150,
    strongest_only: bool = False,
) -> pd.DataFrame:
    """
    Filter and optionally retain only the strongest hit per snoRNA.
    """
    if df.empty:
        return df.copy()

    filtered = df[
        (df["deltaG"] < min_delta_g) &
        (df["duplex_length"] <= max_duplex_length)
    ].copy()

    if filtered.empty:
        return filtered

    filtered = filtered.sort_values("deltaG", ascending=True)

    if strongest_only:
        filtered = (
            filtered.groupby("snoRNA", as_index=False)
            .first()
            .reset_index(drop=True)
        )

    return filtered


if __name__ == "__main__":
    csv_path = Path("data/snoRNA_input.csv")
    output_path = Path("results/tagsnoplex_hits.csv")

    target_rna = (
        "GCCGGGCGCGGUGGCGCGUGCCUGUAGUCCCAGCUACUCGGGAGGCUGAGGCUGGAGGAUCGCUUGAGUCCAGGAG"
        "UUCGGGGCUGUAGUGCGCUAUGCCGAUCGGGUGUCCGCACUAAGUUCGGCAUCAAUAUGGUGACCUCCCGGGAGCGG"
        "GGGACCACCAGGUUGCCUAAGGAGGGGUGAACCGGCCCAGGUCGGAAACGGAGCAGGUCAAAACUCCCGUGCUGAUC"
        "AGUAGUGGGAUCGCGCCUGUGAAUAGCCACUGCACUCCAGCCUGGGCAACAUAGCGAGACCCCGUCUCU"
    )

    snornas = load_snornas_from_csv(csv_path)

    raw_df = scan_snornas_against_target(
        snornas=snornas,
        target_rna=target_rna,
        rnaplex_executable="RNAplex",
        max_interaction_length=40,
        accessibility_length=40,
        energy_cutoff=-18.0,
        max_hits_per_snorna=5,
    )

    final_df = filter_ranked_hits(
        raw_df,
        min_delta_g=-18.0,
        max_duplex_length=150,
        strongest_only=False,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_path, index=False)
    print(f"Saved {len(final_df)} hits to {output_path}")

