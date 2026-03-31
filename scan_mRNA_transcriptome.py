#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from rnaplex_utils import run_rnaplex
from transcriptome_utils import (
    audit_transcript,
    classify_duplex_region,
    load_sno_list_df,
    load_transcript_sequence,
    parse_gtf,
    parse_subset_fasta_headers,
)


def screen_sno_vs_targets(
    sno_df: pd.DataFrame,
    fasta_path: str | Path,
    gtf_df: pd.DataFrame,
    gene_to_ensts: dict[str, set[str]],
    output_dir: str | Path,
    n_transcripts_limit: int | None = None,
    max_interaction_length: int = 40,
    accessibility_length: int = 55,
    energy_max: float = -12.0,
    duplex_max_length: int = 150,
    require_utrs: bool = True,
    tol_frac: float = 0.25,
    rnaplex_executable: str = "RNAplex",
) -> pd.DataFrame:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pairs = [(gene, enst) for gene, ensts in sorted(gene_to_ensts.items()) for enst in sorted(ensts)]
    if n_transcripts_limit is not None:
        pairs = pairs[:n_transcripts_limit]

    print(f"[INFO] Screening {len(pairs)} transcript(s).")

    results = []
    audit_skips = []
    seq_cache: dict[str, str] = {}
    fmap_cache: dict[str, list[dict]] = {}

    for gene, enst in pairs:
        if enst not in seq_cache:
            try:
                seq_cache[enst] = load_transcript_sequence(fasta_path, enst)
            except ValueError:
                audit_skips.append({"gene": gene, "enst": enst, "reason": "no_fasta_seq"})
                continue

        tx_seq = seq_cache[enst]

        audit = audit_transcript(
            gtf_df=gtf_df,
            transcript_id=enst,
            seq_len=len(tx_seq),
            require_utrs=require_utrs,
            tol_frac=tol_frac,
        )

        if not audit["ok"]:
            audit_skips.append(
                {
                    "gene": gene,
                    "enst": enst,
                    "reason": audit["reason"],
                    "len_seq": audit.get("len_seq"),
                    "len_feat": audit.get("len_feat"),
                    "len_5": audit.get("len_5"),
                    "len_cds": audit.get("len_cds"),
                    "len_3": audit.get("len_3"),
                }
            )
            continue

        fmap_cache[enst] = audit["fmap"]

        for _, sno in sno_df.iterrows():
            sno_name = sno["name"]
            sno_seq = sno["sequence"]

            try:
                hits = run_rnaplex(
                    target_seq=tx_seq,
                    query_seq=sno_seq,
                    rnaplex_executable=rnaplex_executable,
                    max_interaction_length=max_interaction_length,
                    accessibility_length=accessibility_length,
                    energy_cutoff=energy_max,
                )

                if not hits:
                    results.append(
                        {
                            "snoRNA": sno_name,
                            "target_gene": gene,
                            "target_enst": enst,
                            "deltaG": 0.0,
                            "duplex_length": 0,
                            "target_region_label": "no_hit",
                            "target_start": 0,
                            "target_end": 0,
                            "snoRNA_start": 0,
                            "snoRNA_end": 0,
                            "duplex_dbn": "",
                            "duplex_seq": "",
                            "pass_filters": False,
                        }
                    )
                    continue

                best_hit = hits[0]
                duplex_seq = f"{best_hit.target_region}&{best_hit.query_region}"
                region = classify_duplex_region(
                    best_hit.target_start,
                    best_hit.target_end,
                    fmap_cache[enst],
                )

                keep = (
                    best_hit.delta_g < energy_max
                    and best_hit.duplex_length <= duplex_max_length
                )

                results.append(
                    {
                        "snoRNA": sno_name,
                        "target_gene": gene,
                        "target_enst": enst,
                        "deltaG": best_hit.delta_g,
                        "duplex_length": best_hit.duplex_length,
                        "target_region_label": region,
                        "target_start": best_hit.target_start,
                        "target_end": best_hit.target_end,
                        "snoRNA_start": best_hit.query_start,
                        "snoRNA_end": best_hit.query_end,
                        "duplex_dbn": best_hit.duplex_dbn,
                        "duplex_seq": duplex_seq,
                        "pass_filters": bool(keep),
                        "len_seq": audit["len_seq"],
                        "len_feat": audit["len_feat"],
                        "len_5": audit["len_5"],
                        "len_cds": audit["len_cds"],
                        "len_3": audit["len_3"],
                    }
                )

            except Exception as exc:
                results.append(
                    {
                        "snoRNA": sno_name,
                        "target_gene": gene,
                        "target_enst": enst,
                        "error": str(exc),
                        "pass_filters": False,
                        "duplex_dbn": "",
                        "duplex_seq": "",
                    }
                )

    if audit_skips:
        pd.DataFrame(audit_skips).to_csv(output_dir / "audit_skips.csv", index=False)

    return pd.DataFrame(results)


def write_region_summaries(
    df_all: pd.DataFrame,
    output_dir: str | Path,
    do_site_pie: bool = True,
    do_top_pair_pie: bool = True,
    dedup_exact_sites: bool = False,
) -> None:
    output_dir = Path(output_dir)
    allowed_regions = ["5'UTR", "CDS", "3'UTR"]

    if df_all.empty or "pass_filters" not in df_all.columns or "target_region_label" not in df_all.columns:
        print("[INFO] No results available for region summaries.")
        return

    df_hits = df_all[
        df_all["pass_filters"] & df_all["target_region_label"].isin(allowed_regions)
    ].copy()

    if df_hits.empty:
        print("[INFO] No passing, region-labeled hits available.")
        return

    df_hits["deltaG"] = pd.to_numeric(df_hits["deltaG"], errors="coerce")

    if dedup_exact_sites:
        before = len(df_hits)
        df_hits = df_hits.drop_duplicates(
            subset=["snoRNA", "target_enst", "target_start", "target_end", "target_region_label"]
        )
        print(f"[INFO] Exact-site deduplication: {before} -> {len(df_hits)}")

    if do_site_pie:
        site_summary = (
            df_hits["target_region_label"]
            .value_counts()
            .rename_axis("region")
            .reset_index(name="count")
        )
        site_summary["fraction"] = site_summary["count"] / site_summary["count"].sum()
        site_summary.to_csv(output_dir / "site_region_summary.csv", index=False)
        df_hits.to_csv(output_dir / "site_hits.csv", index=False)

        fig, ax = plt.subplots()
        ax.pie(site_summary["count"], labels=site_summary["region"], autopct="%1.1f%%", startangle=90)
        ax.set_title(f"Predicted binding by region — site mode (n={int(site_summary['count'].sum())})")
        fig.savefig(output_dir / "site_region_pie.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    if do_top_pair_pie:
        df_top_pair = (
            df_hits.sort_values("deltaG")
            .groupby(["snoRNA", "target_enst"], as_index=False)
            .first()
        )

        top_summary = (
            df_top_pair["target_region_label"]
            .value_counts()
            .rename_axis("region")
            .reset_index(name="count")
        )
        top_summary["fraction"] = top_summary["count"] / top_summary["count"].sum()
        df_top_pair.to_csv(output_dir / "top_pair_hits.csv", index=False)
        top_summary.to_csv(output_dir / "top_pair_region_summary.csv", index=False)

        fig, ax = plt.subplots()
        ax.pie(top_summary["count"], labels=top_summary["region"], autopct="%1.1f%%", startangle=90)
        ax.set_title(
            f"Predicted binding by region — top per (snoRNA, transcript) (n={int(top_summary['count'].sum())})"
        )
        fig.savefig(output_dir / "top_pair_region_pie.png", dpi=200, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Screen snoRNAs against a transcriptome using RNAplex.")
    parser.add_argument("--fasta", required=True, help="Transcript FASTA file")
    parser.add_argument("--sno-list", required=True, help="CSV or Excel file of snoRNAs")
    parser.add_argument("--gtf", required=True, help="GTF annotation file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--rnaplex", default="RNAplex", help="Path to RNAplex executable")
    parser.add_argument("--l", type=int, default=40, help="RNAplex -l value")
    parser.add_argument("--c", type=int, default=55, help="RNAplex -c value")
    parser.add_argument("--energy-max", type=float, default=-12.0, help="Maximum retained energy")
    parser.add_argument("--duplex-maxlen", type=int, default=150, help="Maximum duplex length")
    parser.add_argument("--limit", type=int, default=None, help="Optional transcript limit for debugging")
    parser.add_argument("--allow-missing-utrs", action="store_true", help="Do not require both 5' and 3' UTR")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gtf_df = parse_gtf(args.gtf)
    sno_df = load_sno_list_df(args.sno_list)
    gene_to_ensts = parse_subset_fasta_headers(args.fasta)

    tx_with_cds = set(gtf_df.loc[gtf_df["feature_norm"] == "cds", "transcript_id"].unique())
    gene_to_ensts = {g: {e for e in es if e in tx_with_cds} for g, es in gene_to_ensts.items()}

    df_all = screen_sno_vs_targets(
        sno_df=sno_df,
        fasta_path=args.fasta,
        gtf_df=gtf_df,
        gene_to_ensts=gene_to_ensts,
        output_dir=outdir,
        n_transcripts_limit=args.limit,
        max_interaction_length=args.l,
        accessibility_length=args.c,
        energy_max=args.energy_max,
        duplex_max_length=args.duplex_maxlen,
        require_utrs=not args.allow_missing_utrs,
        rnaplex_executable=args.rnaplex,
    )

    all_csv = outdir / "transcriptome_screen_all_pairs.csv"
    df_all.to_csv(all_csv, index=False)
    print(f"[INFO] Wrote {len(df_all)} rows to {all_csv}")

    write_region_summaries(df_all=df_all, output_dir=outdir)


if __name__ == "__main__":
    main()

