#!/usr/bin/env bash
# Generate minimal synthetic FASTQ files for pipeline testing.
# These contain 100 random paired-end reads (150bp) — enough to exercise
# every process without requiring real reference genomes.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

NREADS=100
READ_LEN=150

generate_fastq() {
    local out="$1"
    local suffix="$2"
    local seed="$3"

    # Use awk to generate deterministic random reads
    awk -v n="$NREADS" -v len="$READ_LEN" -v seed="$seed" -v suf="$suffix" '
    BEGIN {
        srand(seed)
        bases = "ACGT"
        for (i = 1; i <= n; i++) {
            # Header
            printf "@read%d%s\n", i, suf
            # Sequence
            seq = ""
            for (j = 1; j <= len; j++) {
                seq = seq substr(bases, int(rand()*4)+1, 1)
            }
            print seq
            # Plus
            print "+"
            # Quality (all I = Q40)
            qual = ""
            for (j = 1; j <= len; j++) qual = qual "I"
            print qual
        }
    }' | gzip > "$out"
}

echo "Generating test FASTQ files..."
generate_fastq "test_R1.fastq.gz" "/1" 42
generate_fastq "test_R2.fastq.gz" "/2" 43

echo "Generating test reference..."
# Create a tiny genome (1 chromosome, 10kb)
awk 'BEGIN {
    print ">chr1"
    bases = "ACGT"
    srand(99)
    for (i = 1; i <= 100; i++) {
        line = ""
        for (j = 1; j <= 100; j++) {
            line = line substr(bases, int(rand()*4)+1, 1)
        }
        print line
    }
}' | gzip > test_genome.fa.gz

# Create a minimal GTF
printf 'chr1\ttest\tgene\t1\t5000\t.\t+\t.\tgene_id "gene1"; gene_name "TestGene";\nchr1\ttest\texon\t1\t5000\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n' | gzip > test_annotation.gtf.gz

echo "Test data generated in: $SCRIPT_DIR"
ls -lh test_*.fastq.gz test_genome.fa.gz test_annotation.gtf.gz
