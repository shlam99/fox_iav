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

END_TIME=$(date +%s)
echo "Total pipeline runtime for ${BATCH}: $((END_TIME - START_TIME)) seconds"