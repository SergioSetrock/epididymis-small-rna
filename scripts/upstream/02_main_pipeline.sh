#!/bin/bash
set -e

# --- run_pipeline_genomic.sh (v4.6 - Versión Definitiva) ---
#
# Lógica:
# 1. Comando fastp en UNA LÍNEA para evitar errores de '\'.
# 2. Se usa --length_limit 200 (para fastp v1.0.1).
# 3. Se quita --strata de bowtie (para evitar warnings).
# 4. Se usa --trim_poly_x para recorte agresivo de PolyA.
#
# -----------------------------------------------------------------

# --- CONFIGURACIÓN CRÍTICA (Tus archivos) ---
REF_FASTA="04_references/GCF_000001635.20_GRCm38_genomic.fna"
REF_GFF="04_references/mm10.RNA_full.custom.gtf"
REF_BOWTIE_INDEX="04_references/GRCm38_genomic"

DOCKER_IMAGE="smarter_pipeline:v1"
DOCKER_CMD="sudo docker run --rm --user $(id -u):$(id -g) -v $(pwd):/data ${DOCKER_IMAGE}"

# -----------------------------------------------------------------

echo "--- INICIO DEL PIPELINE (v4.6 - Versión Definitiva) ---"

# --- PASO 1: QC INICIAL ---
echo "--- PASO 1: Iniciando QC de datos crudos... ---"
for fq1 in 00_data_raw/*_R1_*.fastq.gz; do
    fq2=$(echo ${fq1} | sed 's/_R1_/_R2_/')
    SAMPLE_BASENAME=$(basename ${fq1})
    SAMPLE=$(echo ${SAMPLE_BASENAME} | sed 's/_R1_.*.fastq.gz//')
    echo "Procesando FastQC para: ${SAMPLE}"
    $DOCKER_CMD fastqc -o 01_qc_raw -t 4 ${fq1} ${fq2}
done

# --- PASO 2: TRIMMING (Comando en una sola línea) ---
echo "--- PASO 2: Iniciando trimming (Comando en una línea)... ---"
mkdir -p 99_reports/fastp

for fq1 in 00_data_raw/*_R1_*.fastq.gz; do
    fq2=$(echo ${fq1} | sed 's/_R1_/_R2_/')
    SAMPLE_BASENAME=$(basename ${fq1})
    SAMPLE=$(echo ${SAMPLE_BASENAME} | sed 's/_R1_.*.fastq.gz//')
    echo "Procesando fastp para: ${SAMPLE}"

    # --- CAMBIO: Comando fastp completo en UNA LÍNEA para evitar errores de formato ---
    $DOCKER_CMD fastp -i ${fq1} -I ${fq2} -o 02_trimmed/${SAMPLE}.trimmed_R1.fastq.gz -O 02_trimmed/${SAMPLE}.trimmed_R2.fastq.gz -h 99_reports/fastp/${SAMPLE}.fastp.html -j 99_reports/fastp/${SAMPLE}.fastp.json --trim_front1=3 -a AAAAAAAAAA --trim_poly_x --length_required 15 -q 20 --length_limit 200

done

# --- PASO 3: QC POST-TRIMMING ---
echo "--- PASO 3: Iniciando QC de datos limpios... ---"
$DOCKER_CMD fastqc -o 03_qc_trimmed -t 4 02_trimmed/*.trimmed_R1.fastq.gz

# --- PASO 4.1: CONSTRUIR ÍNDICE DE REFERENCIA (Genoma) ---
echo "--- PASO 4.1: Construyendo índice Bowtie del genoma (GRCm38)... ---"
$DOCKER_CMD bowtie-build ${REF_FASTA} ${REF_BOWTIE_INDEX}

# --- PASO 4.2: ALINEAMIENTO Y ORDENADO (Corregido) ---
echo "--- PASO 4.2: Iniciando alineamiento al genoma... ---"
for fq1_trimmed in 02_trimmed/*.trimmed_R1.fastq.gz; do
    SAMPLE=$(basename ${fq1_trimmed} .trimmed_R1.fastq.gz)
    echo "Alineando: ${SAMPLE}"

    # 1. Alineamiento (CORREGIDO: se quitó --strata)
    $DOCKER_CMD bowtie \
        -x ${REF_BOWTIE_INDEX} \
        -q ${fq1_trimmed} \
        -S 05_mapping/${SAMPLE}.sam \
        -v 1 -k 1 --best \
        --un 05_mapping/${SAMPLE}.unmapped.fq

    # 2. Convertir SAM a BAM
    $DOCKER_CMD samtools view -bS \
        05_mapping/${SAMPLE}.sam \
        -o 05_mapping/${SAMPLE}.bam

    # 3. Ordenar el BAM
    $DOCKER_CMD samtools sort \
        05_mapping/${SAMPLE}.bam \
        -o 05_mapping/${SAMPLE}.sorted.bam

    # 4. Indexar el BAM
    $DOCKER_CMD samtools index 05_mapping/${SAMPLE}.sorted.bam

    rm 05_mapping/${SAMPLE}.sam 05_mapping/${SAMPLE}.bam
done

# --- PASO 5: CUANTIFICACIÓN (Con GTF personalizado) ---
echo "--- PASO 5: Iniciando cuantificación (GTF personalizado)... ---"
BAM_LIST=$(ls 05_mapping/*.sorted.bam | tr '\n' ' ')

$DOCKER_CMD featureCounts \
    -a ${REF_GFF} \
    -o 06_counts/counts_matrix.txt \
    -t exon \
    -g gene_id \
    -M \
    -T 4 \
    ${BAM_LIST}

# --- PASO 6: REPORTE AGREGADO ---
echo "--- PASO 6: Generando reporte MultiQC... ---"
$DOCKER_CMD multiqc . -o 99_reports --force

echo "--- PIPELINE COMPLETADO (v4.6 - Versión Definitiva) ---"
echo "Revisa la matriz de conteos en: 06_counts/counts_matrix.txt"
echo "Revisa el reporte final en: 99_reports/multiqc_report.html"
