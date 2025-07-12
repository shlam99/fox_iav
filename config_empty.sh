#!/bin/bash

# Sample Prefix
SAMPLE_PREFIX=""

# Batch Name
BATCH=""

# Sample Name
get_sample_id() {
    local barcode_id="$1"
    case "$barcode_id" in
        "barcode01") echo "" ;;
        "barcode02") echo "" ;;
        "barcode03") echo "" ;;
        "barcode04") echo "" ;;
        "barcode05") echo "" ;;
        "barcode06") echo "" ;;
        "barcode07") echo "" ;;
        "barcode08") echo "" ;;
        "barcode09") echo "" ;;
        "barcode10") echo "" ;;
        "barcode11") echo "" ;;
        "barcode12") echo "" ;;
        "barcode13") echo "" ;;
        "barcode14") echo "" ;;
        "barcode15") echo "" ;;
        "barcode16") echo "" ;;
        "barcode17") echo "" ;;
        "barcode18") echo "" ;;
        "barcode19") echo "" ;;
        "barcode20") echo "" ;;
        "barcode21") echo "" ;;
        "barcode22") echo "" ;;
        "barcode23") echo "" ;;
        "barcode24") echo "" ;;
        *) echo "$barcode_id" ;;  # Fallback to barcode ID if no match
    esac
}

# Reference Genome
REFERENCE_GENOME=""

# Filtering Parameters
THREADS=8
MIN_QUALITY=7
MIN_LENGTH=500
MAX_LENGTH=50000
TARGET_BASES=100000000
KEEP_PERCENT=90
