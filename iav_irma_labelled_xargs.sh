#!/bin/bash
source config.sh

# Export all variables needed for parallel processing
export SAMPLE_PREFIX BATCH THREADS MIN_LENGTH TARGET_BASES KEEP_PERCENT

# Function to display step completion
step_complete() {
    echo ""
    echo "========================================"
    echo "Step $1 complete: $2"
    echo "========================================"
    echo ""
}

START_TIME=$(date +%s)

########################################
# Step 0: Setup and dependency checks
########################################
echo "Starting fox IAV pipeline on ${BATCH}..."
echo "Sample prefix is ${SAMPLE_PREFIX}..."
echo "Step 0: Setting up directory structure, checking input files, and checking dependencies..."

# Check for required tools
command -v samtools >/dev/null 2>&1 || { echo >&2 "Error: samtools not found."; exit 1; }
command -v filtlong >/dev/null 2>&1 || { echo >&2 "Error: filtlong not found."; exit 1; }
command -v IRMA >/dev/null 2>&1 || { echo >&2 "Error: IRMA not found."; exit 1; }

# Convert BAM to FASTQ (if needed)
if ls *.bam >/dev/null 2>&1; then
    echo "Found BAM files, converting to FASTQ in parallel..."
    for i in {1..24}; do
        BARCODE_PADDED=$(printf "%02d" "$i")
        BAM_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.bam"
        FASTQ_OUT="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
        
        if [ -f "$BAM_IN" ]; then
            samtools fastq "$BAM_IN" | gzip > "$FASTQ_OUT" &
        fi
        
        # Limit concurrent jobs to $THREADS (from config.sh)
        if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
            wait -n
        fi
    done
    wait  # Ensure all jobs finish
elif ls *.fastq >/dev/null 2>&1; then
    echo "Found FASTQ files, compressing in parallel..."
    find . -maxdepth 1 -name "*.fastq" -print0 | xargs -0 -P $THREADS -I {} sh -c 'echo "Compressing {}..."; gzip {}'
else
    echo "Error: No input files found (.bam or .fastq)"
fi

step_complete "0" "Setup and dependency checks"

########################################
# Step 1: Quality filtering with filtlong (parallelized)
########################################
echo "Step 1: Running filtlong in parallel (xargs)..."

# Make directory
mkdir -p qc_reads/qc_logs

seq 1 24 | xargs -P $THREADS -I {} bash -c '
    BARCODE_PADDED=$(printf "%02d" "$1")
    FASTQ_IN="${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.fastq.gz"
    FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
    
    if [ -f "$FASTQ_IN" ]; then
        echo "Processing $FASTQ_IN..."
        filtlong --min_length $MIN_LENGTH \
                 --keep_percent $KEEP_PERCENT \
                 --target_bases $TARGET_BASES \
                 "$FASTQ_IN" 2> "qc_reads/qc_logs/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtlong.log" | \
        gzip > "$FASTQ_FILTERED"
    else
        echo "Warning: $FASTQ_IN not found" >&2
    fi
' _ {}

step_complete "1" "Quality filtering complete"

echo ""
echo "========================================"
echo "Output files: qc_reads/${SAMPLE_PREFIX}_barcode01-24.filtered.fastq.gz"
echo "========================================"

########################################
# Step 2: Run IRMA in parallel (background jobs)
########################################
echo "Step 2: Running IRMA with $THREADS parallel jobs..."

# Make directory
mkdir -p irma_results

for i in {1..24}; do
    (
        BARCODE_PADDED=$(printf "%02d" "$i")
        FASTQ_FILTERED="qc_reads/${SAMPLE_PREFIX}_barcode${BARCODE_PADDED}.filtered.fastq.gz"
        IRMA_OUTDIR="irma_results/barcode${BARCODE_PADDED}"
        
        if [ -f "$FASTQ_FILTERED" ]; then
            echo "Running IRMA on $FASTQ_FILTERED..."
            IRMA FLU_ont "$FASTQ_FILTERED" "$IRMA_OUTDIR"
        fi
    ) &
    
    # Limit concurrent IRMA jobs (RAM/CPU heavy)
    if [[ $(jobs -r -p | wc -l) -ge $THREADS ]]; then
        wait -n
    fi
done
wait  # Wait for all IRMA jobs

step_complete "2" "IRMA analysis complete"

echo ""
echo "========================================"
echo "Output files: irma_results/barcode01-24/amended_consensus/barcode01-24.fa"
echo "========================================"

########################################
# Step 3: Pool consensus sequences by segment into BatchID/BarcodeXX/SampleID labelled consensus
########################################
echo "Step 3: Pool consensus sequences by segment into BatchID/BarcodeXX/SampleID labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Segment names corresponding to the numbers
SEGMENTS=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")

# Initialize output files for each segment
for seg in "${SEGMENTS[@]}"; do
    > "irma_consensus/${seg}_consensus_${BATCH}.fasta"
done

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}/amended_consensus"
    
    echo "Processing ${BARCODE_ID}..."

    # Get sample ID for batch file
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")

    # Process each segment (1-8)
    for seg_num in {1..8}; do
        CONSENSUS_SOURCE="${IRMA_OUTDIR}/${BARCODE_ID}_${seg_num}.fa"
        SEGMENT_NAME="${SEGMENTS[$((seg_num-1))]}"
        OUTPUT_FILE="irma_consensus/${SEGMENT_NAME}_consensus_${BATCH}.fasta"
        
        # Check if consensus file exists
        if [ ! -f "$CONSENSUS_SOURCE" ]; then
            echo "WARNING: Consensus file not found for ${BARCODE_ID} segment ${SEGMENT_NAME}" >&2
            continue
        fi

        # Add to pooled consensus with barcode prefix + sample ID (Step 3 format)
        sed "s/^>/>${BATCH}_barcode${BARCODE_PADDED}|${SAMPLE_ID}|/" "$CONSENSUS_SOURCE" >> "$OUTPUT_FILE"
    done
done

# Count sequences in each segment file and remove empty ones
for seg in "${SEGMENTS[@]}"; do
    OUTPUT_FILE="irma_consensus/${seg}_consensus_${BATCH}.fasta"
    if [ ! -s "$OUTPUT_FILE" ]; then
        echo "Removing empty file: ${OUTPUT_FILE}"
        rm "$OUTPUT_FILE"
    else
        COUNT=$(grep -c "^>" "$OUTPUT_FILE")
        echo "${seg}: $COUNT sequences"
    fi
done

step_complete "3" "Pooled!! BatchID/BarcodeXX/SampleID segmented IAV sequences"

echo ""
echo "========================================"
echo "Step 3 output files:"
for seg in "${SEGMENTS[@]}"; do
    echo "- ${seg}: irma_consensus/${seg}_consensus_${BATCH}.fasta"
done
echo "========================================"

########################################
# Step 4: Pool consensus sequences by segment into SampleID labelled consensus
########################################
echo "Step 4: Pool consensus sequences by segment into SampleID only labelled consensus for ${BATCH}..."

# Make directory
mkdir -p irma_consensus

# Initialize output files for each segment
for seg in "${SEGMENTS[@]}"; do
    > "irma_consensus/${seg}_${BATCH}.fasta"
done

for i in {1..24}; do
    BARCODE_PADDED=$(printf "%02d" "$i")
    BARCODE_ID="barcode${BARCODE_PADDED}"
    IRMA_OUTDIR="irma_results/${BARCODE_ID}/amended_consensus"
    
    # Get sample ID
    SAMPLE_ID=$(get_sample_id "$BARCODE_ID")

    # Process each segment (1-8)
    for seg_num in {1..8}; do
        CONSENSUS_SOURCE="${IRMA_OUTDIR}/${BARCODE_ID}_${seg_num}.fa"
        SEGMENT_NAME="${SEGMENTS[$((seg_num-1))]}"
        OUTPUT_FILE="irma_consensus/${SEGMENT_NAME}_${BATCH}.fasta"
        
        # Check if consensus file exists
        if [ ! -f "$CONSENSUS_SOURCE" ]; then
            continue
        fi

        # Add to pooled consensus with sample ID only (Step 4 format)
        sed "s/^>.*/>${SAMPLE_ID}/" "$CONSENSUS_SOURCE" >> "$OUTPUT_FILE"
    done
done

# Count sequences in each segment file and remove empty ones
for seg in "${SEGMENTS[@]}"; do
    OUTPUT_FILE="irma_consensus/${seg}_${BATCH}.fasta"
    if [ ! -s "$OUTPUT_FILE" ]; then
        echo "Removing empty file: ${OUTPUT_FILE}"
        rm "$OUTPUT_FILE"
    else
        COUNT=$(grep -c "^>" "$OUTPUT_FILE")
        echo "${seg}: $COUNT sequences"
    fi
done

step_complete "4" "Pooled!! SampleID-labeled segmented IAV sequences"

echo ""
echo "========================================"
echo "Step 4 output files:"
for seg in "${SEGMENTS[@]}"; do
    echo "- ${seg}: irma_consensus/${seg}_${BATCH}.fasta"
done
echo "========================================"

########################################
# Step 5: Nextclade Analysis for all IAV segments
########################################
echo "Step 5: Running Nextclade for all IAV segments..."
echo "Using batch ID: ${BATCH}"

# Make directories for Nextclade results
mkdir -p nextclade_results/even_segments nextclade_results/odd_segments

# Mac-compatible segment definitions (using arrays instead of associative arrays)
SEGMENTS=("HA" "NA" "PB1" "NS" "PB2" "PA" "NP" "MP")
REF_NAMES=("A_HA_H9" "A_NA_N2" "A_PB1" "A_NS" "A_PB2" "A_PA" "A_NP" "A_MP")
SEGMENT_TYPES=("even" "even" "even" "even" "odd" "odd" "odd" "odd")

# Use explicit path to IRMA's consensus.fasta
REF_FILE="${HOME}/miniconda3/bin/IRMA_RES/modules/FLU_ont/reference/consensus.fasta"

if [[ ! -f "$REF_FILE" ]]; then
    echo "ERROR: Could not find IRMA reference file at $REF_FILE" >&2
    exit 1
fi

echo "Using IRMA reference file: $REF_FILE"

# First extract each reference from consensus.fasta
TEMP_DIR=$(mktemp -d)

for i in "${!SEGMENTS[@]}"; do
    segment="${SEGMENTS[$i]}"
    ref_name="${REF_NAMES[$i]}"
    
    # Mac-compatible reference extraction
    awk -v RS=">" -v ref="$ref_name" '$1 == ref {print ">" $0}' "$REF_FILE" | grep -v '^$' > "${TEMP_DIR}/${segment}_ref.fasta"
    
    if [[ ! -s "${TEMP_DIR}/${segment}_ref.fasta" ]]; then
        echo "ERROR: Could not extract reference $ref_name from $REF_FILE" >&2
        exit 1
    fi
done

# Process each segment with appropriate output folder
for i in "${!SEGMENTS[@]}"; do
    segment="${SEGMENTS[$i]}"
    INPUT_FILE="irma_consensus/${segment}_consensus_${BATCH}.fasta"
    
    # Determine output directory
    if [[ "${SEGMENT_TYPES[$i]}" == "even" ]]; then
        OUTPUT_DIR="nextclade_results/even_segments"
    else
        OUTPUT_DIR="nextclade_results/odd_segments"
    fi
    
    OUTPUT_PREFIX="${OUTPUT_DIR}/${segment}_consensus_${BATCH}"
    SEGMENT_REF="${TEMP_DIR}/${segment}_ref.fasta"
    
    if [[ -f "$INPUT_FILE" ]]; then
        echo "Processing $segment segment (${REF_NAMES[$i]}) → $OUTPUT_DIR"
        echo "Input file: $INPUT_FILE"
        
        nextclade run \
            "$INPUT_FILE" \
            --output-csv "${OUTPUT_PREFIX}_nextclade.csv" \
            --output-json "${OUTPUT_PREFIX}_nextclade.json" \
            --output-fasta "${OUTPUT_PREFIX}_aligned.fasta" \
            --input-ref "$SEGMENT_REF" \
            --output-tree "${OUTPUT_PREFIX}_tree.json" \
            --output-ndjson "${OUTPUT_PREFIX}_nextclade.ndjson" \
            --jobs $THREADS
        
        if [ $? -eq 0 ]; then
            echo "Successfully processed $segment segment"
        else
            echo "Error processing $segment segment with Nextclade" >&2
        fi
    else
        echo "$segment input file $INPUT_FILE not found, skipping analysis"
    fi
done

# Clean up temporary files
rm -rf "$TEMP_DIR"

step_complete "5" "Nextclade analysis complete for all IAV segments"


END_TIME=$(date +%s)
echo "Total pipeline runtime for ${BATCH}: $((END_TIME - START_TIME)) seconds"
