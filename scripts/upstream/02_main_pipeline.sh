#!/bin/bash
set -e

# --- run_pipeline_genomic.sh (v4.6 - Final Version) ---
#
# Logic:
# 1. fastp command on ONE LINE to avoid '\' continuation errors.
# 2. Used --length_limit 200 (for fastp v1.0.1).
# 3. Removed --strata from bowtie (to avoid warnings).
# 4. Used --trim_poly_x for aggressive PolyA trimming.
#
# -----------------------------------------------------------------

# --- CRITICAL CONFIGURATION (Your files) ---
REF_FASTA="04_references/GCF_000001635.20_GRCm38_genomic.fna"
REF_GFF="04_references/mm10.RNA_full.custom.gtf"
REF_BOWTIE_INDEX="04_references/GRCm38_genomic"

DOCKER_IMAGE="smarter_pipeline:v1"
DOCKER_CMD="sudo docker run --rm --user $(id -u):$(id -g) -v $(pwd):/data ${DOCKER_IMAGE}"

# -----------------------------------------------------------------

echo "--- STARTING PIPELINE (v4.6 - Final Version) ---"

# --- STEP 1: INITIAL QC ---
echo "--- STEP 1: Starting raw data QC... ---"
for fq1 in 00_data_raw/*_R1_*.fastq.gz; do
    fq2=$(echo ${fq1} | sed 's/_R1_/_R2_/')
    SAMPLE_BASENAME=$(basename ${fq1})
    SAMPLE=$(echo ${SAMPLE_BASENAME} | sed 's/_R1_.*.fastq.gz//')
    echo "Processing FastQC for: ${SAMPLE}"
    $DOCKER_CMD fastqc -o 01_qc_raw -t 4 ${fq1} ${fq2}
done

# --- STEP 2: TRIMMING (One-line command) ---
echo "--- STEP 2: Starting trimming (One-line command)... ---"
mkdir -p 99_reports/fastp

for fq1 in 00_data_raw/*_R1_*.fastq.gz; do
    fq2=$(echo ${fq1} | sed 's/_R1_/_R2_/')
    SAMPLE_BASENAME=$(basename ${fq1})
    SAMPLE=$(echo ${SAMPLE_BASENAME} | sed 's/_R1_.*.fastq.gz//')
    echo "Processing fastp for: ${SAMPLE}"

    # --- CHANGE: Complete fastp command on ONE LINE to avoid formatting errors ---
    $DOCKER_CMD fastp -i ${fq1} -I ${fq2} -o 02_trimmed/${SAMPLE}.trimmed_R1.fastq.gz -O 02_trimmed/${SAMPLE}.trimmed_R2.fastq.gz -h 99_reports/fastp/${SAMPLE}.fastp.html -j 99_reports/fastp/${SAMPLE}.fastp.json --trim_front1=3 -a AAAAAAAAAA --trim_poly_x --length_required 15 -q 20 --length_limit 200

done

# --- STEP 3: POST-TRIMMING QC ---
echo "--- STEP 3: Starting clean data QC... ---"
$DOCKER_CMD fastqc -o 03_qc_trimmed -t 4 02_trimmed/*.trimmed_R1.fastq.gz

# --- STEP 4.1: BUILD REFERENCE INDEX (Genome) ---
echo "--- STEP 4.1: Building Bowtie genome index (GRCm38)... ---"
$DOCKER_CMD bowtie-build ${REF_FASTA} ${REF_BOWTIE_INDEX}

# --- STEP 4.2: ALIGNMENT AND SORTING ---
echo "--- STEP 4.2: Starting genome alignment... ---"
for fq1_trimmed in 02_trimmed/*.trimmed_R1.fastq.gz; do
    SAMPLE=$(basename ${fq1_trimmed} .trimmed_R1.fastq.gz)
    echo "Aligning: ${SAMPLE}"

    # 1. Alignment (CORRECTED: removed --strata)
    $DOCKER_CMD bowtie \
        -x ${REF_BOWTIE_INDEX} \
        -q ${fq1_trimmed} \
        -S 05_mapping/${SAMPLE}.sam \
        -v 1 -k 1 --best \
        --un 05_mapping/${SAMPLE}.unmapped.fq

    # 2. Convert SAM to BAM
    $DOCKER_CMD samtools view -bS \
        05_mapping/${SAMPLE}.sam \
        -o 05_mapping/${SAMPLE}.bam

    # 3. Sort the BAM
    $DOCKER_CMD samtools sort \
        05_mapping/${SAMPLE}.bam \
        -o 05_mapping/${SAMPLE}.sorted.bam

    # 4. Index the BAM
    $DOCKER_CMD samtools index 05_mapping/${SAMPLE}.sorted.bam

    rm 05_mapping/${SAMPLE}.sam 05_mapping/${SAMPLE}.bam
done

# --- STEP 5: QUANTIFICATION (With custom GTF) ---
echo "--- STEP 5: Starting quantification (Custom GTF)... ---"
BAM_LIST=$(ls 05_mapping/*.sorted.bam | tr '\n' ' ')

$DOCKER_CMD featureCounts \
    -a ${REF_GFF} \
    -o 06_counts/counts_matrix.txt \
    -t exon \
    -g gene_id \
    -M \
    -T 4 \
    ${BAM_LIST}

# --- STEP 6: AGGREGATED REPORT ---
echo "--- STEP 6: Generating MultiQC report... ---"
$DOCKER_CMD multiqc . -o 99_reports --force

echo "--- PIPELINE COMPLETED (v4.6 - Final Version) ---"
echo "Check the count matrix at: 06_counts/counts_matrix.txt"
echo "Check the final report at: 99_reports/multiqc_report.html"
