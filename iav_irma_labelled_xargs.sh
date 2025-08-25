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
# Step 5: Nextclade CLI Analysis for all IAV segments (macOS compatible)
########################################
echo "Step 5: Running Nextclade CLI for all IAV segments..."
echo "Using batch ID: ${BATCH}"

# Make directories for Nextclade results
mkdir -p nextclade_results/even_segments nextclade_results/odd_segments

# Define segments and their reference names in consensus.fasta
declare -A segments=(
    # Even-numbered segments (HA=4, NA=6, PB1=2, NS=8)
    ["HA"]="A_HA_H9"
    ["NA"]="A_NA_N2" 
    ["PB1"]="A_PB1"
    ["NS"]="A_NS"
    
    # Odd-numbered segments (PB2=1, PA=3, NP=5, MP=7)
    ["PB2"]="A_PB2"
    ["PA"]="A_PA"
    ["NP"]="A_NP"
    ["MP"]="A_MP"
)

# First extract each reference from consensus.fasta
REF_FILE="consensus.fasta"
TEMP_DIR=$(mktemp -d)

# macOS compatible reference extraction (using awk instead of head -n -1)
for segment in "${!segments[@]}"; do
    ref_name="${segments[$segment]}"
    # macOS compatible way to extract reference sequences
    awk -v ref="$ref_name" '
    BEGIN {RS=">"; ORS=""}
    $1 == ref {
        print ">" $0
        exit
    }' "$REF_FILE" > "${TEMP_DIR}/${segment}_ref.fasta"
    
    if [[ ! -s "${TEMP_DIR}/${segment}_ref.fasta" ]]; then
        echo "ERROR: Could not extract reference $ref_name from $REF_FILE" >&2
        exit 1
    fi
done

# Process each segment with appropriate output folder
for segment in "${!segments[@]}"; do
    INPUT_FILE="irma_consensus/${segment}_consensus_${BATCH}.fasta"
    
    # Determine output directory based on segment type
    case $segment in
        HA|NA|PB1|NS)
            OUTPUT_DIR="nextclade_results/even_segments"
            ;;
        PB2|PA|NP|MP)
            OUTPUT_DIR="nextclade_results/odd_segments"
            ;;
        *)
            echo "Unknown segment: $segment" >&2
            continue
            ;;
    esac
    
    OUTPUT_PREFIX="${OUTPUT_DIR}/${segment}_consensus_${BATCH}"
    SEGMENT_REF="${TEMP_DIR}/${segment}_ref.fasta"
    
    if [[ -f "$INPUT_FILE" ]]; then
        echo "Processing $segment segment (${segments[$segment]}) → $OUTPUT_DIR"
        echo "Input file: $INPUT_FILE"
        
        # Check if nextclade is available
        if ! command -v nextclade &> /dev/null; then
            echo "Error: nextclade command not found. Please install Nextclade CLI for macOS." >&2
            echo "Installation instructions: https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html" >&2
            exit 1
        fi
        
        # Run Nextclade with macOS compatible options
        nextclade run \
            --input-dataset "$SEGMENT_REF" \
            --input-fasta "$INPUT_FILE" \
            --output-csv "${OUTPUT_PREFIX}_nextclade.csv" \
            --output-json "${OUTPUT_PREFIX}_nextclade.json" \
            --output-fasta "${OUTPUT_PREFIX}_aligned.fasta" \
            --output-tree "${OUTPUT_PREFIX}_tree.json" \
            --output-ndjson "${OUTPUT_PREFIX}_nextclade.ndjson"
        
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

step_complete "5" "Nextclade CLI analysis complete for all IAV segments"


END_TIME=$(date +%s)
echo "Total pipeline runtime for ${BATCH}: $((END_TIME - START_TIME)) seconds"
