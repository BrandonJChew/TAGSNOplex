# TAGSNOplex

TAGSNOplex is a simple workflow for predicting **TAG-snoRNA interactions** with mRNA transcripts.

It uses:

- **RNAplex** for RNA–RNA interaction prediction
- a **reference transcript FASTA** for transcript sequences
- a **GTF annotation file** to label hits as:
  - 5′UTR
  - CDS
  - 3′UTR

---

## What you need

### Files included in this repository

You should provide:

- `Supplementary_Table_1.csv`  
  List of target transcripts to screen

- `TAGsno_DBN.csv`  
  List of TAG-snoRNAs with sequence and dot-bracket structure

- Python scripts:
  - `extract_transcript_fasta.py`
  - `rnaplex_utils.py`
  - `varna_utils.py`
  - `transcriptome_utils.py`
  - `scan_snornas.py`
  - `scan_mrna_transcriptome.py`

---

## Files the user must download separately

The user must also download:

- a **reference transcript FASTA**
- a **matching GTF annotation file**

For example:
- `gencode.v49.transcripts.fa`
- `gencode_v49.gtf`

Or a newer GENCODE release.

**Important:** the FASTA and GTF should come from the **same annotation release**.

---

## Required input files

Place the following files in the `input/` folder:

### 1. `Supplementary_Table_1.csv`
This file contains the target transcript IDs.

It must include a column containing transcript IDs, for example:

```csv
transcript_id
ENST00000372285.8
ENST00000412345.6
