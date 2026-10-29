#!/usr/bin/env bash
#
# RNA-seq pipeline: FastQC -> (Trimmomatic) -> HISAT2 -> samtools -> featureCounts
#
# Paired-end Illumina reads from lobster gill (g_) and muscle (m_) tissue,
# aligned to the Hamer_2.0 reference genome (NCBI).
#
# Usage:  ./rnaseq_pipeline.sh [qc | trim | index | align | count | all]
# Run from the project directory. Expected layout:
#   fastq/    raw paired-end reads (<prefix>_1.fq[.gz], <prefix>_2.fq[.gz])
#   HISAT2/   Hamer_2.0_genomic.fna and genomic.gtf from NCBI
#
# Requires: fastqc, trimmomatic, hisat2, samtools, subread (featureCounts)

set -euo pipefail

# ---------- Settings ----------
THREADS=8
GENOME_FASTA="HISAT2/Hamer_2.0_genomic.fna"
GTF_RAW="HISAT2/genomic.gtf"
GTF="HISAT2/new_genomic.gtf"        # cleaned annotation (see prep_gtf)
INDEX="HISAT2/lobgenome/genome"     # HISAT2 index base name
ADAPTERS="TruSeq2-PE.fa"            # give the full path if Trimmomatic can't find it

# 1 = align trimmed reads (used for the final analysis); 0 = align raw reads
USE_TRIMMED=1
QUANT_DIR="quants_trimmed"   # per-sample featureCounts output

# Sample prefixes (file name minus _1/_2 and extension). The sample ID is the
# first two underscore-separated fields (e.g. g_003_12_12 -> g_003).
PREFIXES=(
  g_003_12_12 g_005_12_12 g_007_12_13 g_009_11_3
  g_014_12_15 g_020_12_15 g_022_11_3  g_032_12_15
  m_003_12_18 m_005_12_15 m_007_12_18 m_009_12_19 m_014_12_18
  m_015_12_19 m_020_12_19 m_022_12_19 m_032_12_19
)

mkdir -p fastqc trimmed logs bam "$QUANT_DIR" HISAT2/lobgenome

sample_id() { echo "$1" | cut -d_ -f1-2; }

# Sets R1 and R2 for a sample prefix (raw or trimmed, per USE_TRIMMED)
get_reads() {
  local prefix=$1
  if [[ $USE_TRIMMED -eq 1 ]]; then
    R1="trimmed/${prefix}_1_paired.fq.gz"
    R2="trimmed/${prefix}_2_paired.fq.gz"
  else
    R1=$(ls fastq/"${prefix}"_1.fq* | head -n 1)
    R2=$(ls fastq/"${prefix}"_2.fq* | head -n 1)
  fi
}

# ---------- Steps ----------
run_qc() {
  for p in "${PREFIXES[@]}"; do
    fastqc -t "$THREADS" -o fastqc/ fastq/"${p}"_1.fq* fastq/"${p}"_2.fq*
  done
  echo "FastQC (raw reads) finished."
}

run_trim() {
  # Removes adapters, trims low-quality ends, then a 4-base sliding window (avg quality 15)
  for p in "${PREFIXES[@]}"; do
    local r1 r2
    r1=$(ls fastq/"${p}"_1.fq* | head -n 1)
    r2=$(ls fastq/"${p}"_2.fq* | head -n 1)
    trimmomatic PE -threads "$THREADS" -phred33 "$r1" "$r2" \
      "trimmed/${p}_1_paired.fq.gz" "trimmed/${p}_1_unpaired.fq.gz" \
      "trimmed/${p}_2_paired.fq.gz" "trimmed/${p}_2_unpaired.fq.gz" \
      ILLUMINACLIP:"${ADAPTERS}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
    fastqc -t "$THREADS" -o fastqc/ "trimmed/${p}_1_paired.fq.gz" "trimmed/${p}_2_paired.fq.gz"
  done
  echo "Trimmomatic and FastQC (trimmed reads) finished."
}

run_index() {
  hisat2-build -f -p "$THREADS" "$GENOME_FASTA" "$INDEX"
  echo "HISAT2 index built."
}

run_align() {
  # Align and sort in one pass; HISAT2 alignment summaries go to logs/
  for p in "${PREFIXES[@]}"; do
    local id
    id=$(sample_id "$p")
    get_reads "$p"
    hisat2 -q --threads "$THREADS" -x "$INDEX" -1 "$R1" -2 "$R2" \
      --summary-file "logs/${id}_hisat2_summary.txt" \
      | samtools sort -@ 4 -O bam -o "bam/${id}_sorted.bam" -
  done
  echo "HISAT2 alignment and sorting finished."
}

prep_gtf() {
  # The NCBI GTF has records with empty gene_id values that make featureCounts fail.
  # Step 1: keep only the gene_id attribute in column 9.
  # Step 2: drop records whose gene_id is empty.
  if [[ ! -f "$GTF" ]]; then
    sed 's/\(gene_id "[^"]*\).*/\1"/' "$GTF_RAW" | grep -v 'gene_id ""' > "$GTF"
  fi
}

run_count() {
  prep_gtf
  # One count file per sample ($QUANT_DIR/featurecounts_<id>.txt); the R script merges them.
  # -p: paired-end. With Subread >= 2.0.2, add --countReadPairs to count fragments.
  for p in "${PREFIXES[@]}"; do
    local id
    id=$(sample_id "$p")
    featureCounts -p -T "$THREADS" -a "$GTF" \
      -o "${QUANT_DIR}/featurecounts_${id}.txt" "bam/${id}_sorted.bam"
  done
  echo "featureCounts finished."
}

# ---------- Main ----------
SECONDS=0
case "${1:-all}" in
  qc)    run_qc ;;
  trim)  run_trim ;;
  index) run_index ;;
  align) run_align ;;
  count) run_count ;;
  all)   run_qc; [[ $USE_TRIMMED -eq 1 ]] && run_trim; run_index; run_align; run_count ;;
  *)     echo "Usage: $0 [qc|trim|index|align|count|all]"; exit 1 ;;
esac
echo "$((SECONDS / 60)) minutes and $((SECONDS % 60)) seconds elapsed."
