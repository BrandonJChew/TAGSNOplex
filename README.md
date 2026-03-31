# TAGSNOplex

TAGSNOplex is a simple workflow for predicting TAG-snoRNA interactions with mRNA transcripts.

It uses:

- **RNAplex** for RNA–RNA interaction prediction
- a reference **transcript FASTA** for transcript sequences
- a matching **GTF annotation** file to label hits as:
  - 5′UTR
  - CDS
  - 3′UTR

## What you need

### Files included in this repository

#### Input files

- **Supplementary_Table_1.csv**  
  List of target transcripts to screen

- **TAGsno_DBN.csv**  
  List of TAG-snoRNAs with sequence and dot-bracket structure

#### Python scripts

- **extract_transcript_fasta.py**  
  Extracts transcript sequences of interest from the reference transcript FASTA using the transcript IDs in `Supplementary_Table_1.csv`

- **rnaplex_utils.py**  
  Helper functions for running RNAplex and parsing RNA–RNA interaction results

- **varna_utils.py**  
  Helper functions for generating structure visualizations of predicted interactions

- **transcriptome_utils.py**  
  Helper functions for working with transcript annotations and assigning predicted hits to 5′UTR, CDS, or 3′UTR regions using the GTF file

- **scan_snornas.py**  
  Scans one or more TAG-snoRNAs against selected target transcripts and records predicted interactions

- **scan_mrna_transcriptome.py**  
  Main transcriptome-scale workflow for screening TAG-snoRNAs against the transcript set and summarizing annotated interaction hits

### Files the user must download separately

The user must also download:

- a reference **transcript FASTA**
- a matching **GTF annotation file**

For example:

- `gencode.v49.transcripts.fa`
- `gencode.v49.annotation.gtf`

Or a newer GENCODE release.

**Important:** the FASTA and GTF must come from the same annotation release.

---

## Required input files

Place the following files in the `input/` folder.

### 1. `Supplementary_Table_1.csv`

This file contains the target transcript IDs.

It must include a column containing transcript IDs, for example:

```csv
transcript_id
ENST00000372285.8
ENST00000412345.6
