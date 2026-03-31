# TAGSNOplex

TAGSNOplex is a simple workflow for predicting TAG-snoRNA interactions with mRNA transcripts.

It uses:

- **RNAplex** for RNA–RNA interaction prediction
- a reference **transcript FASTA** for transcript sequences
- a matching **GTF annotation** file to label hits as:
  - 5′UTR
  - CDS
  - 3′UTR

## What is included

This package contains:

- `extract_transcript_fasta.py`
- `rnaplex_utils.py`
- `varna_utils.py`
- `transcriptome_utils.py`
- `scan_snornas.py`
- `scan_mrna_transcriptome.py`
- `Supplementary_Table_1.csv`
- `TAGsno_DBN.csv`

## What the user must download separately

The user must also download:

- a reference **transcript FASTA**
- a matching **GTF annotation** file

For example:

- `gencode.v49.transcripts.fa`
- `gencode.v49.annotation.gtf`

A newer GENCODE release can also be used.

**Important:** the FASTA and GTF must come from the same annotation release.

---

## Before you start

This workflow is intended to be run on your own computer.

These scripts are best run from:

- a **terminal / command prompt**
- or a **Jupyter Notebook** using commands such as `!python ...`

A terminal is usually the easiest option.

---

## Suggested folder setup

Create a folder on your computer called `TAGSNOplex` and place the files like this:

```text
TAGSNOplex/
│
├── extract_transcript_fasta.py
├── rnaplex_utils.py
├── varna_utils.py
├── transcriptome_utils.py
├── scan_snornas.py
├── scan_mrna_transcriptome.py
├── README.md
│
├── input/
│   ├── Supplementary_Table_1.csv
│   ├── TAGsno_DBN.csv
│   ├── gencode.v49.transcripts.fa
│   └── gencode.v49.annotation.gtf
│
└── output/
```

## What each script does

extract_transcript_fasta.py:
Extracts the transcript sequences of interest from the full reference transcript FASTA using the transcript IDs listed in Supplementary_Table_1.csv.
rnaplex_utils.py:
Contains helper functions for running RNAplex and parsing RNA–RNA interaction results.
varna_utils.py:
Contains helper functions for generating structure visualizations of predicted interactions.
transcriptome_utils.py:
Contains helper functions for working with transcript annotation and assigning hits to 5′UTR, CDS, or 3′UTR regions using the GTF file.
scan_snornas.py:
Runs TAG-snoRNA interaction scans for selected inputs.
scan_mrna_transcriptome.py:
Main workflow for screening TAG-snoRNAs against the extracted transcript set and saving annotated results.

## Step-by-step workflow
Step 1. Download all required files
Download the following to your computer:
extract_transcript_fasta.py
rnaplex_utils.py
varna_utils.py
transcriptome_utils.py
scan_snornas.py
scan_mrna_transcriptome.py
Supplementary_Table_1.csv
TAGsno_DBN.csv
a reference transcript FASTA
a matching GTF annotation file

Step 2. Create a working folder
Create a folder called TAGSNOplex on your computer.
Inside that folder, place:
the six .py files
README.md
Then create two subfolders:
input/
output/

Step 3. Place the input files in the input/ folder
Put the following files into input/:
Supplementary_Table_1.csv
TAGsno_DBN.csv
gencode.v49.transcripts.fa
gencode.v49.annotation.gtf
Your folder should now look like this:

```
TAGSNOplex/
│
├── extract_transcript_fasta.py
├── rnaplex_utils.py
├── varna_utils.py
├── transcriptome_utils.py
├── scan_snornas.py
├── scan_mrna_transcriptome.py
├── README.md
│
├── input/
│   ├── Supplementary_Table_1.csv
│   ├── TAGsno_DBN.csv
│   ├── gencode.v49.transcripts.fa
│   └── gencode.v49.annotation.gtf
│
└── output/

```
Step 4. Open a terminal or Jupyter Notebook
Move into the TAGSNOplex folder.
In a terminal, this means changing into that folder before running commands.
Example:
cd /path/to/TAGSNOplex
If using Jupyter Notebook, first move into the folder with:
%cd /path/to/TAGSNOplex

Step 5. Extract the target transcript sequences
Run the following command:
```
python extract_transcript_fasta.py \
  --table input/Supplementary_Table_1.csv \
  --transcript-column transcript_id \
  --input-fasta input/gencode.v49.transcripts.fa \
  --output-fasta output/target_transcripts.fasta
```
This step reads the transcript IDs from Supplementary_Table_1.csv, looks them up in the reference transcript FASTA, and writes the matching transcript sequences to:
output/target_transcripts.fasta

Step 6. Run the TAG-snoRNA transcriptome scan
After output/target_transcripts.fasta has been created, run:
```
python scan_mrna_transcriptome.py \
  --fasta output/target_transcripts.fasta \
  --sno-list input/TAGsno_DBN.csv \
  --gtf input/gencode.v49.annotation.gtf \
  --outdir output/mrna_screen
```
  
This step uses:
the extracted transcript FASTA
the TAG-snoRNA list
the matching GTF annotation
RNAplex
to predict TAG-snoRNA interactions and save the results in:
output/mrna_screen/

Step 7. Review the output
After the scan finishes, check the output/ folder and the output/mrna_screen/ folder for the generated results.
Running from Jupyter Notebook
These scripts are command-line scripts, but they can also be run from a Jupyter Notebook.
Example:
%cd /path/to/TAGSNOplex
Then run:
```
!python extract_transcript_fasta.py \
  --table input/Supplementary_Table_1.csv \
  --transcript-column transcript_id \
  --input-fasta input/gencode.v49.transcripts.fa \
  --output-fasta output/target_transcripts.fasta
```
Then run:
```
!python scan_mrna_transcriptome.py \
  --fasta output/target_transcripts.fasta \
  --sno-list input/TAGsno_DBN.csv \
  --gtf input/gencode.v49.annotation.gtf \
  --outdir output/mrna_screen
 ``` 
## Summary of the order of operations

Download the TAGSNOplex files to your computer.
Download a matching reference transcript FASTA and GTF annotation file.
Create the input/ and output/ folders.
Put Supplementary_Table_1.csv, TAGsno_DBN.csv, the FASTA, and the GTF into input/.
Run extract_transcript_fasta.py.
Run scan_mrna_transcriptome.py.
