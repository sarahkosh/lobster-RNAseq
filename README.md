[README.md](https://github.com/user-attachments/files/34106680/README.md)
# Lobster RNA-seq: gill and muscle differential expression

An RNA-seq analysis I led during my PhD, comparing gene expression in American lobster (*Homarus americanus*) gill and muscle tissue between a control and warming, acidification, and combined (multistressor) treatments.

**Status:** unpublished, and not part of my dissertation. This repository contains the analysis workflow only. Raw data and results are not included.

## Overview

| Item | Details |
|---|---|
| Data | Paired-end Illumina RNA-seq, sequenced by Novogene |
| Samples | 8 gill and 9 muscle libraries |
| Reference | *H. americanus* genome, assembly Hamer_2.0 (NCBI) |

## Workflow

1. **Quality control:** FastQC on raw reads
2. **Trimming:** Trimmomatic (adapter removal, quality trimming; parameters in the script)
3. **Alignment:** HISAT2 to the reference genome (about 94-95% average overall alignment; per-sample summaries written to `logs/`)
4. **Sorting:** samtools
5. **Quantification:** featureCounts against the NCBI annotation (gene-level counts)
6. **Differential expression:** DESeq2 in R, control vs. each treatment, separately for each tissue

Steps 1-5 are in `rnaseq_pipeline.sh`. Differential expression is in `deseq2_analysis.R`, which reads the per-sample count files and a `metadata.csv` (columns: `sample`, `treatment`, `tissue_type`).

## Usage

```bash
# from the project directory, with raw reads in fastq/ and the genome and annotation in HISAT2/
./rnaseq_pipeline.sh all      # or one step: qc | trim | index | align | count
```

Alignment uses trimmed reads by default. Set `USE_TRIMMED=0` in the script to align raw reads instead.

```r
# then, in R, from the project directory
source("deseq2_analysis.R")
```

## Notes

- The NCBI annotation had records with empty `gene_id` values that caused featureCounts to fail. The script filters them out (`prep_gtf`).
- Differential expression: genes were called significant at adjusted p < 0.05 and |log2 fold change| > 1. Low-count genes were removed before testing (total count of at least 10 in gill, at least 5 in muscle). Contrasts are control vs. treatment, so a positive fold change means higher expression in the control.
- PCA plots with some muscle samples removed are for visualization only; all samples are included in the differential expression models.

## Contact

Sarah Koshak, PhD
