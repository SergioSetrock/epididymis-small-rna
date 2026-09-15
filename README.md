# Epididymal Small RNA-Seq Pipeline

This repository contains the integrated bioinformatic pipeline for the processing, alignment, and downstream statistical analysis of the **Small RNA fraction (≤ 200 nt)** from mouse epididymal tissues. 

---

## 🧬 Custom Annotation Strategy
To prevent mapping bias against diverse small non-coding RNA species inherently caused by sequential cascade-alignment approaches, a comprehensive, custom GTF annotation file was constructed. 

*   **Base Annotation:** GENCODE vM25.
*   **Integrated Databases:** The base GTF was programmatically merged with the **mm10 tRNAscan-SE database** and the **piRNAdb (v1.7.6)** mm10 database.
*   **Implementation:** The merge was executed using the `rtracklayer` and `dplyr` packages in R (`scripts/upstream/01_custom_gtf_build.R`).

---

## 🚀 Upstream Processing Workflow

### 1. Quality Control & Trimming
Raw `.fastq.gz` sequences undergo initial QC via `FastQC`. Adapter trimming and quality filtering are performed using `fastp`, optimized for short-read single-end libraries.

### 2. Alignment (Bowtie)
Trimmed single-end reads from the small RNA fraction were aligned to the GRCm38 reference genome using **Bowtie (v1.3.1)**. 
Alignment parameters were optimized for short sequences:
*   `-v 1`: Allowing a maximum of 1 mismatch.
*   `-k 1 --best`: Restricting to the best valid alignment.
*   Coordinate-sorted BAM files were output via `samtools`.

### 3. Quantification (featureCounts)
Transcript abundance was quantified against the custom integrated GTF using `featureCounts`. 
Specific parameters include:
*   `-M`: Counting multi-mapping reads.
*   `-t exon`: Summarizing at the exon level.
*   `-g gene_id`: Grouping by gene ID.

---

## 📈 Downstream Analysis
The `scripts/downstream/` directory contains R scripts for downstream evaluation of the count matrices:
*   **Biotype Distribution:** `analyze_small_rna_biotypes.R` calculates global shifts across rRNA, snRNA, snoRNA, miRNA, and tRNAs.
*   **Mitochondrial Stress:** `analyze_mt_tRNAs.R` aggregates and evaluates global Mt_tRNA abundance.
*   **Differential Expression:** `analyze_sncRNAs.R` identifies statistically significant intergenerational sRNA shifts.
*   **miRNA Target Prediction:** `predict_mirna_targets.R` and `enrich_mirna_targets.R` utilize `multiMiR` to computationally predict and map mRNA targets for altered epididymal miRNAs.
