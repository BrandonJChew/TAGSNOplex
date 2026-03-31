#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd


FeatureMap = List[dict]


def parse_gtf(gtf_file: str | Path) -> pd.DataFrame:
    """
    Parse a GTF file and retain exon/CDS/UTR-related rows.
    Transcript IDs are stored versionless for matching against FASTA records.
    """
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )

    gtf["feature_norm"] = gtf["feature"].astype(str).str.lower()
    keep_features = {"exon", "cds", "five_prime_utr", "three_prime_utr", "utr"}
    gtf = gtf[gtf["feature_norm"].isin(keep_features)].copy()

    def extract_attr(field: str, attr_str: str) -> Optional[str]:
        match = re.search(fr'{field} "([^"]+)"', str(attr_str))
        return match.group(1) if match else None

    gtf["gene_id"] = gtf["attribute"].apply(lambda s: extract_attr("gene_id", s))
    gtf["transcript_id_full"] = gtf["attribute"].apply(
        lambda s: extract_attr("transcript_id", s)
    )
    gtf["transcript_id"] = gtf["transcript_id_full"].apply(
        lambda x: x.split(".")[0] if isinstance(x, str) else x
    )
    gtf["length"] = (gtf["end"] - gtf["start"] + 1).astype(int)

    gtf["feature_label"] = gtf["feature_norm"].map(
        {
            "exon": "EXON",
            "cds": "CDS",
            "five_prime_utr": "5'UTR",
            "three_prime_utr": "3'UTR",
            "utr": "UTR",
        }
    )

    return gtf[
        [
            "seqname",
            "start",
            "end",
            "strand",
            "feature_norm",
            "feature_label",
            "transcript_id",
            "length",
        ]
    ].copy()


def load_sno_list_df(
    path: str | Path,
    name_col: str = "gene_name",
    seq_col: str = "Sequence",
) -> pd.DataFrame:
    """
    Load snoRNAs from CSV or Excel and normalize sequences to RNA alphabet.
    """
    path = Path(path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)

    df = df[[name_col, seq_col]].dropna().copy()
    df["name"] = df[name_col].astype(str).str.strip()
    df["sequence"] = (
        df[seq_col].astype(str).str.strip().str.upper().str.replace("T", "U", regex=False)
    )
    return df[["name", "sequence"]]


def parse_subset_fasta_headers(fasta_path: str | Path) -> Dict[str, set[str]]:
    """
    Return {gene_symbol -> set(versionless_ENST)} from either:
    - MANE-style headers
    - GENCODE-style headers
    """
    gene_to_ensts: Dict[str, set[str]] = {}

    with open(fasta_path, "r") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            raw = line[1:].strip()

            mane_enst = re.search(r"ENST:([A-Z0-9]+(?:\.\d+)?)", raw)
            mane_gene = re.search(r"Gene:([A-Za-z0-9_.-]+)", raw)

            if mane_enst:
                enst_base = mane_enst.group(1).split(".")[0]
                gene_symbol = mane_gene.group(1) if mane_gene else enst_base
            else:
                parts = raw.split("|")
                if not parts:
                    continue
                enst_base = parts[0].split(".", 1)[0]
                gene_symbol = parts[5].strip() if len(parts) >= 6 and parts[5] else enst_base

            gene_to_ensts.setdefault(gene_symbol, set()).add(enst_base)

    return gene_to_ensts


def load_transcript_sequence(fasta_path: str | Path, transcript_id: str) -> str:
    """
    Load one transcript sequence from FASTA using versionless ENST matching.
    Sequence is returned as RNA (T->U).
    """
    want = transcript_id.split(".")[0]
    take = False
    seq_lines: List[str] = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if take and seq_lines:
                    break

                raw = line[1:].strip()
                mane_enst = re.search(r"ENST:([A-Z0-9]+(?:\.\d+)?)", raw)

                if mane_enst:
                    enst_base = mane_enst.group(1).split(".")[0]
                else:
                    parts = raw.split("|")
                    enst_base = parts[0].split(".", 1)[0] if parts else ""

                take = enst_base == want
            else:
                if take:
                    seq_lines.append(line.strip())

    seq = "".join(seq_lines)
    if not seq:
        raise ValueError(f"Transcript {transcript_id} not found in FASTA.")

    return seq.upper().replace("T", "U")


def merge_intervals(intervals: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping or adjacent intervals.
    """
    if not intervals:
        return []

    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]

    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return merged


def split_exon_by_cds(
    exon_interval: Tuple[int, int],
    cds_intervals: Sequence[Tuple[int, int]],
) -> List[Tuple[int, int, str]]:
    """
    Split one exon into CDS and non-CDS subsegments.
    """
    exon_start, exon_end = exon_interval
    pieces: List[Tuple[int, int, str]] = []
    position = exon_start

    for cds_start, cds_end in cds_intervals:
        if cds_end < exon_start or cds_start > exon_end:
            continue

        if position < cds_start:
            pieces.append((position, min(exon_end, cds_start - 1), "UTR"))

        pieces.append((max(position, cds_start), min(exon_end, cds_end), "CDS"))
        position = cds_end + 1

        if position > exon_end:
            break

    if position <= exon_end:
        pieces.append((position, exon_end, "UTR"))

    return [(a, b, label) for a, b, label in pieces if a <= b]


def build_cdna_map(gtf_df: pd.DataFrame, transcript_id: str) -> Optional[FeatureMap]:
    """
    Build a transcript-space map of 5'UTR / CDS / 3'UTR coordinates.
    """
    tx = transcript_id.split(".")[0]
    sub = gtf_df[gtf_df["transcript_id"] == tx].copy()
    if sub.empty:
        return None

    strand = sub["strand"].iloc[0]

    has_5 = (sub["feature_norm"] == "five_prime_utr").any()
    has_3 = (sub["feature_norm"] == "three_prime_utr").any()

    if has_5 or has_3:
        sub2 = sub[sub["feature_norm"].isin({"five_prime_utr", "cds", "three_prime_utr"})].copy()
        if strand == "+":
            sub2 = sub2.sort_values(["start", "end"])
        else:
            sub2 = sub2.sort_values(["start", "end"], ascending=False)

        cdna_pos = 1
        fmap: FeatureMap = []
        for _, row in sub2.iterrows():
            start_cdna = cdna_pos
            end_cdna = cdna_pos + int(row["length"]) - 1
            fmap.append(
                {
                    "start_cdna": start_cdna,
                    "end_cdna": end_cdna,
                    "feature": row["feature_label"],
                }
            )
            cdna_pos = end_cdna + 1
        return fmap

    exons = sub[sub["feature_norm"] == "exon"].copy()
    cds = sub[sub["feature_norm"] == "cds"].copy()

    if exons.empty or cds.empty:
        return None

    if strand == "+":
        exons = exons.sort_values(["start", "end"])
        cds = cds.sort_values(["start", "end"])
    else:
        exons = exons.sort_values(["start", "end"], ascending=False)
        cds = cds.sort_values(["start", "end"], ascending=False)

    cds_intervals = merge_intervals([(int(r.start), int(r.end)) for _, r in cds.iterrows()])

    fmap: FeatureMap = []
    cdna_pos = 1
    seen_cds = False

    for _, exon in exons.iterrows():
        pieces = split_exon_by_cds((int(exon.start), int(exon.end)), cds_intervals)

        for genomic_start, genomic_end, label in pieces:
            seg_len = genomic_end - genomic_start + 1

            if label == "CDS":
                feature = "CDS"
                seen_cds = True
            else:
                feature = "5'UTR" if not seen_cds else "3'UTR"

            fmap.append(
                {
                    "start_cdna": cdna_pos,
                    "end_cdna": cdna_pos + seg_len - 1,
                    "feature": feature,
                }
            )
            cdna_pos += seg_len

    if not fmap:
        return None

    merged_fmap = [fmap[0]]
    for seg in fmap[1:]:
        prev = merged_fmap[-1]
        if (
            seg["feature"] == prev["feature"]
            and seg["start_cdna"] == prev["end_cdna"] + 1
        ):
            prev["end_cdna"] = seg["end_cdna"]
        else:
            merged_fmap.append(seg)

    return merged_fmap


def classify_duplex_region(target_start: int, target_end: int, fmap: Optional[FeatureMap]) -> str:
    """
    Assign a hit to the transcript feature with the greatest overlap.
    """
    if not fmap or target_start == 0 or target_end == 0:
        return "unknown"

    best_feature = "unknown"
    best_overlap = 0

    for seg in fmap:
        overlap_start = max(target_start, seg["start_cdna"])
        overlap_end = min(target_end, seg["end_cdna"])
        if overlap_start <= overlap_end:
            overlap = overlap_end - overlap_start + 1
            if overlap > best_overlap:
                best_feature = seg["feature"]
                best_overlap = overlap

    return best_feature


def audit_transcript(
    gtf_df: pd.DataFrame,
    transcript_id: str,
    seq_len: int,
    require_utrs: bool = True,
    tol_frac: float = 0.25,
) -> dict:
    """
    Check whether transcript feature lengths are consistent with transcript sequence length.
    """
    fmap = build_cdna_map(gtf_df, transcript_id)
    if not fmap:
        return {"ok": False, "reason": "no_map"}

    len_5 = sum(seg["end_cdna"] - seg["start_cdna"] + 1 for seg in fmap if seg["feature"] == "5'UTR")
    len_cds = sum(seg["end_cdna"] - seg["start_cdna"] + 1 for seg in fmap if seg["feature"] == "CDS")
    len_3 = sum(seg["end_cdna"] - seg["start_cdna"] + 1 for seg in fmap if seg["feature"] == "3'UTR")
    feature_sum = len_5 + len_cds + len_3

    have_5 = len_5 > 0
    have_3 = len_3 > 0

    ok_len = (seq_len > 0) and (abs(feature_sum - seq_len) <= max(20, int(tol_frac * seq_len)))
    ok_utrs = (not require_utrs) or (have_5 and have_3)

    return {
        "ok": bool(ok_len and ok_utrs),
        "reason": None if (ok_len and ok_utrs) else ("len_mismatch" if not ok_len else "missing_utr"),
        "len_seq": int(seq_len),
        "len_feat": int(feature_sum),
        "len_5": int(len_5),
        "len_cds": int(len_cds),
        "len_3": int(len_3),
        "fmap": fmap,
    }

