# TAGSNOplex

**TAGSNOplex** is a structure-informed computational pipeline for predicting **TAG-snoRNA interactions** across a curated transcriptome. 
The workflow combines **R2DT** for query secondary-structure generation, **RNAplex** for thermodynamic RNA–RNA interaction prediction, and **ViennaRNA utilities** for visualization of prioritized duplexes.

---

## Overview

TAGSNOplex was designed to identify candidate interactions between query TAG-snoRNAs and a transcript set enriched for genes involved in **ER translocation** and the **secretory pathway**. 
For each query RNA, dot-bracket secondary structures are generated in **R2DT template-based mode**, and duplex interactions with target transcripts are then scanned using **RNAplex**. 
Predicted interactions are parsed, annotated by transcript region, collapsed to remove redundant coordinate-identical hits, ranked by binding energy, and visualized for downstream interpretation.

---

## Pipeline summary

The workflow consists of five major stages:

1. **Construct the TAG-snoRNA targetome**
2. **Generate query RNA secondary structures**
3. **Predict RNA–RNA interactions with RNAplex**
4. **Annotate, filter, and rank interaction calls**
5. **Visualize prioritized duplexes**

---

## Inputs

TAGSNOplex requires the following inputs:

### 1. Query RNAs
For each TAG-snoRNA:
- query name
- RNA sequence
- dot-bracket secondary structure

### 2. Curated target transcriptome
For each transcript:
- transcript ID
- gene name
- transcript sequence
- transcript-region coordinates:
  - 5′UTR
  - CDS
  - 3′UTR

### 3. Optional annotations
- UniProt functional annotations
- signal peptide annotations
- MANE Select transcript mappings

---

## Outputs

The pipeline produces:

- a table of predicted RNAplex interactions
- transcript-region annotations for each interaction
- a nonredundant ranked hit table
- optional per-query or per-target summary tables
- secondary-structure plots for prioritized interactions

Typical output fields include:

- `query_id`
- `target_transcript_id`
- `query_start`
- `query_end`
- `target_start`
- `target_end`
- `duplex_dbn`
- `deltaG`
- `interaction_length`
- `target_region`

---

## Workflow

## 1. Construction of the TAG-snoRNA targetome

Candidate target genes are compiled from **UniProt** annotations and filtered to retain proteins annotated as:

- secreted
- integral membrane
- cell-surface

Where applicable, proteins may also be required to contain an annotated **N-terminal signal peptide**, consistent with ER entry and secretory-pathway trafficking.

To minimize isoform ambiguity and standardize coordinate reporting, **one representative transcript per gene** is selected from the **MANE Select** set. Only transcripts with complete:

- 5′UTR
- CDS
- 3′UTR

are retained, enabling consistent assignment of predicted interaction sites to transcript regions.

---

## 2. Query secondary-structure generation

For each query TAG-snoRNA, secondary structure is generated using **R2DT** in **template-based mode**.

R2DT outputs:
- a sequence in FASTA-compatible format
- dot-bracket secondary-structure notation
- standardized 2D orientation for visualization

The resulting dot-bracket structures are used as structural inputs for downstream interaction prediction and plotting.

---

## 3. RNA–RNA interaction prediction and transcript scanning

RNA–RNA interactions between each query RNA and each target transcript are predicted using **RNAplex** in query/target mode.

### RNAplex settings
- `-q`: query RNA
- `-t`: target transcript
- `-l 40`: maximum interaction length of 40 nt
- `-e -18`: retain only duplexes with predicted free energy below -18 kcal/mol

For each reported interaction, RNAplex returns:
- duplex dot-bracket notation
- interacting coordinate ranges in the query and target
- predicted duplex free energy

These outputs are parsed and stored for downstream filtering and annotation.

---

## 4. Post-processing, annotation, and ranking

Predicted RNAplex hits are aggregated by query RNA and annotated by transcript region using retained transcript feature definitions.

Each interaction is assigned to:
- `5UTR`
- `CDS`
- `3UTR`

or marked as unannotated if it falls outside retained transcript features.

### Ranking and filtering
Hits are ranked primarily by **predicted interaction energy** (`deltaG`, most negative first). Additional attributes retained for prioritization include:

- interaction length
- transcript region
- query identity
- target transcript identity

### Redundancy handling
Redundant hits with identical:
- query coordinate span
- target coordinate span

are collapsed into a single record, retaining the **lowest-energy instance**.

---

## 5. Structure visualization

For prioritized interactions, duplex visualizations are generated using **ViennaRNA utilities**.

Where supported by the RNAplex output mode, plots are produced directly from RNAplex-derived duplex structures and may include:
- query and target labels
- interaction coordinates
- predicted free energy
- highlighted interaction spans

These plots are intended for downstream inspection and figure generation.

---

## Pseudocode

```text
INPUT:
    query RNAs with sequences and dot-bracket structures
    curated target transcriptome with 5′UTR/CDS/3′UTR annotations

STEP 1:
    Build a targetome enriched for secreted, membrane, and cell-surface genes
    Optionally require signal peptide annotation
    Select one MANE Select transcript per gene
    Retain only transcripts with complete 5′UTR, CDS, and 3′UTR

STEP 2:
    Generate or import dot-bracket secondary structures for each query RNA using R2DT

STEP 3:
    For each query RNA:
        For each target transcript:
            Run RNAplex with:
                -q query
                -t target
                -l 40
               

STEP 4:
    Parse RNAplex output to extract:
        query coordinates
        target coordinates
        duplex dot-bracket structure
        free energy
        interaction length

STEP 5:
    Annotate each interaction by transcript region:
        5′UTR, CDS, or 3′UTR

STEP 6:
    Collapse redundant hits with identical query–target coordinate spans
    Retain the lowest-energy instance

STEP 7:
    Rank interactions by predicted free energy
    Retain metadata for downstream filtering and prioritization

STEP 8:
    Generate duplex structure plots for prioritized interactions
