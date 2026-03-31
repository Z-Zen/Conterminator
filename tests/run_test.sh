#!/usr/bin/env bash
# Run a minimal pipeline test using synthetic data.
# This validates that:
#  1. The pipeline parses inputs without errors
#  2. Sample sheet validation works
#  3. Processes that don't need real tools (subsetting, copying) execute correctly
#  4. Parameter validation catches missing required params
#
# Usage:
#   bash tests/run_test.sh              # run from repo root
#   bash tests/run_test.sh --full       # also run alignment (requires tools installed)

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_DATA="$SCRIPT_DIR/data"
WORK_DIR=$(mktemp -d)

cleanup() {
    rm -rf "$WORK_DIR"
}
trap cleanup EXIT

echo "=========================================="
echo "  Conterminator Pipeline Test Suite"
echo "=========================================="
echo "Test data: $TEST_DATA"
echo "Work dir:  $WORK_DIR"
echo ""

# Setup: create test reference directory structure
REF_DIR="$WORK_DIR/references/TestRef"
mkdir -p "$REF_DIR"
cp "$TEST_DATA/test_genome.fa.gz" "$REF_DIR/"
cp "$TEST_DATA/test_annotation.gtf.gz" "$REF_DIR/"

# Create sample sheet with absolute paths
sed "s|TESTDIR|$TEST_DATA|g" "$TEST_DATA/test_sample_sheet.tsv" > "$WORK_DIR/sample_sheet.tsv"

PASS=0
FAIL=0

run_test() {
    local name="$1"
    local expect_fail="${2:-false}"
    shift 2
    echo -n "TEST: $name ... "

    local outdir="$WORK_DIR/results_${name// /_}"

    if "$@" > "$WORK_DIR/log_${name// /_}.txt" 2>&1; then
        if [ "$expect_fail" = "true" ]; then
            echo "FAIL"
            echo "  Expected an error but pipeline succeeded"
            FAIL=$((FAIL + 1))
        else
            echo "PASS"
            PASS=$((PASS + 1))
        fi
    else
        if [ "$expect_fail" = "true" ]; then
            echo "PASS"
            PASS=$((PASS + 1))
        else
            echo "FAIL"
            echo "  Log: $WORK_DIR/log_${name// /_}.txt"
            FAIL=$((FAIL + 1))
        fi
    fi
}

# -------------------------------------------
# Test 1: Missing required params should fail
# -------------------------------------------
run_test "missing_outdir" true \
    ~/nextflow run "$REPO_DIR/main.nf" \
    --sample_sheet "$WORK_DIR/sample_sheet.tsv" \
    -profile local \
    -work-dir "$WORK_DIR/work1"

# -------------------------------------------
# Test 2: Missing sample_sheet should fail
# -------------------------------------------
run_test "missing_sample_sheet" true \
    ~/nextflow run "$REPO_DIR/main.nf" \
    --outdir "$WORK_DIR/results_missing_ss" \
    -profile local \
    -work-dir "$WORK_DIR/work2"

# -------------------------------------------
# Test 3: Help flag should succeed
# -------------------------------------------
echo -n "TEST: help_flag ... "
timeout 30 ~/nextflow run "$REPO_DIR/main.nf" --help > "$WORK_DIR/log_help_flag.txt" 2>&1 || true
if grep -q "Help Documentation" "$WORK_DIR/log_help_flag.txt" 2>/dev/null; then
    echo "PASS"
    PASS=$((PASS + 1))
else
    echo "FAIL"
    echo "  Log: $WORK_DIR/log_help_flag.txt"
    FAIL=$((FAIL + 1))
fi

# -------------------------------------------
# Test 4: Pipeline runs with all tools disabled (dry-run style)
# -------------------------------------------
echo -n "TEST: minimal_run_no_tools ... "
~/nextflow run "$REPO_DIR/main.nf" \
    --sample_sheet "$WORK_DIR/sample_sheet.tsv" \
    --outdir "$WORK_DIR/results_minimal" \
    --strains_base_dir "$WORK_DIR/references" \
    --standard_references_dir "$WORK_DIR/references" \
    --star_index_dir "$WORK_DIR/star_idx" \
    --contamination_blast_dbs "$WORK_DIR/blast_dbs" \
    --fastq_screen_conf "$WORK_DIR/dummy_fqs.conf" \
    --run_star_alignment false \
    --run_fastqc false \
    --run_fastq_screen false \
    --run_deeptools false \
    --run_picard_gc false \
    --run_bedtools_gc false \
    --run_mapinsights false \
    --run_qualimap false \
    --run_contamination_check false \
    --run_multiqc false \
    -profile local \
    -work-dir "$WORK_DIR/work4" > "$WORK_DIR/log_minimal_run_no_tools.txt" 2>&1
RC=$?
if [ $RC -eq 0 ] || grep -q "Pipeline Completed" "$WORK_DIR/log_minimal_run_no_tools.txt" 2>/dev/null; then
    echo "PASS"
    PASS=$((PASS + 1))
else
    echo "FAIL (exit code: $RC)"
    echo "  Log: $WORK_DIR/log_minimal_run_no_tools.txt"
    FAIL=$((FAIL + 1))
fi

# -------------------------------------------
# Test 5: Invalid sample sheet format should fail
# -------------------------------------------
echo -e "wrong_header\tcol2" > "$WORK_DIR/bad_sample_sheet.tsv"
run_test "invalid_sample_sheet" true \
    ~/nextflow run "$REPO_DIR/main.nf" \
    --sample_sheet "$WORK_DIR/bad_sample_sheet.tsv" \
    --outdir "$WORK_DIR/results_bad_ss" \
    --strains_base_dir "$WORK_DIR/references" \
    --standard_references_dir "$WORK_DIR/references" \
    --star_index_dir "$WORK_DIR/star_idx" \
    --contamination_blast_dbs "$WORK_DIR/blast_dbs" \
    --fastq_screen_conf "$WORK_DIR/dummy_fqs.conf" \
    --run_star_alignment false \
    --run_contamination_check false \
    -profile local \
    -work-dir "$WORK_DIR/work5"

# -------------------------------------------
# Summary
# -------------------------------------------
echo ""
echo "=========================================="
echo "  Results: $PASS passed, $FAIL failed"
echo "=========================================="

if [ "$FAIL" -gt 0 ]; then
    echo "Some tests failed. Check logs in: $WORK_DIR"
    # Don't cleanup on failure so logs can be inspected
    trap - EXIT
    exit 1
fi
