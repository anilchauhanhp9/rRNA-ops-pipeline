#!/bin/bash

#################################################################
# rRNA-ops_pipeline - Full-Length Bacterial rRNA Operon Analysis
# 
# Usage: ./rRNA-ops_pipeline.sh [STEP] -c config.yaml -t threads
#################################################################

set -euo pipefail

# Pipeline version
VERSION="1.0"

# Default values
THREADS=4
CONFIG_FILE=""
VERBOSE=false
STEP="all"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
ORANGE='\033[0;33m'
NC='\033[0m' # No Color

# Function to print colored messages
log_info() {
    echo -e "${GREEN}[INFO]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_step() {
    echo -e "${CYAN}[STEP]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_success() {
    echo -e "${MAGENTA}[SUCCESS]${NC} $(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Function to display help
show_help() {
    cat << EOF
${CYAN}rRNA-ops_pipeline v${VERSION}${NC} - Full-Length Bacterial rRNA Operon Analysis Pipeline

${YELLOW}Usage:${NC} $(basename "$0") [STEP] [OPTIONS]

${YELLOW}Available Steps:${NC}
    ${GREEN}all${NC}                  Run complete pipeline
    ${GREEN}basecall${NC}             Basecalling with Dorado (optional)
    ${GREEN}demultiplex${NC}          Demultiplexing with Porechop (optional)
    ${GREEN}qc_before_filter${NC}     Quality assessment BEFORE filtering with NanoPlot
    ${GREEN}filter${NC}               Read filtering with Chopper
    ${GREEN}qc_after_filter${NC}      Quality assessment AFTER filtering with NanoPlot
    ${GREEN}blast_reads${NC}          BLAST all reads combined (pre-clustering)
    ${GREEN}cluster${NC}              Sequence clustering (vsearch or swarm)
    ${GREEN}polish${NC}               Polish with Medaka (optional)
    ${GREEN}taxonomy${NC}             Taxonomic annotation with BLAST
    ${GREEN}metadata${NC}             Generate sample metadata file

    ${YELLOW}Multiple Steps:${NC} Use comma-separated list (no spaces)
    Example: demultiplex,qc_before_filter,filter,qc_after_filter,cluster,taxonomy

${YELLOW}Required Arguments:${NC}
    -c, --config FILE       Path to configuration YAML file

${YELLOW}Optional Arguments:${NC}
    -t, --threads INT       Number of threads (default: 4)
    -v, --verbose           Enable verbose output
    -h, --help              Show this help message
    --version               Show version information

${YELLOW}Examples:${NC}
    # Run complete pipeline
    $(basename "$0") all -c config.yaml -t 60

    # Run full workflow step by step
    $(basename "$0") demultiplex,qc_before_filter,filter,qc_after_filter,cluster,taxonomy -c config.yaml -t 60

    # Run from filtering onwards (demux already done)
    $(basename "$0") qc_before_filter,filter,qc_after_filter,cluster,taxonomy -c config.yaml -t 60

    # Run clustering and taxonomy only
    $(basename "$0") cluster,taxonomy -c config.yaml -t 60

    # Generate metadata only
    $(basename "$0") metadata -c config.yaml

${YELLOW}Step Dependencies:${NC}
    basecall -> demultiplex -> qc_before_filter -> filter -> qc_after_filter -> cluster -> taxonomy

${YELLOW}Note:${NC} Optional steps can be skipped if input files are already available.

EOF
    exit 0
}

# Function to show version
show_version() {
    echo "rRNA-ops_pipeline version ${VERSION}"
    exit 0
}

# Parse command line arguments
if [[ $# -eq 0 ]]; then
    show_help
fi

# First argument is the step
STEP="$1"
shift

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--config)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -h|--help)
            show_help
            ;;
        --version)
            show_version
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            ;;
    esac
done

# Parse multiple steps (comma-separated)
IFS=',' read -ra STEPS <<< "$STEP"

# Validate each step
VALID_STEPS=("all" "basecall" "demultiplex" "qc_before_filter" "filter" "qc_after_filter" "qc" "blast_reads" "cluster" "polish" "taxonomy" "metadata")
while IFS= read -r step; do
    if [[ ! " ${VALID_STEPS[*]} " =~ " ${step} " ]]; then
        log_error "Invalid step: $step"
        show_help
    fi
done < <(printf '%s\n' "${STEPS[@]}")

# Check if config file is provided
if [[ -z "$CONFIG_FILE" ]]; then
    log_error "Configuration file is required. Use -c or --config option."
    exit 1
fi

if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "Configuration file not found: $CONFIG_FILE"
    exit 1
fi

# Display banner
echo -e "${CYAN}"
cat << "EOF"
          _____  _   _    _                           
  _ __   |  __ \| \ | |  / \            ___   _ __  ___
 | '__|  | |__) |  \| | / _ \   ____   / _ \ | '_ \/ __|
 | |     |  _  /| |\  |/ ___ \ |____| | (_) || |_) \__ \
 |_|     |_| \_\| |   /_/   \_\        \___/ | .__/|___/
                                             | |        
                                             |_|        
EOF
echo -e "${NC}"
echo ""
echo "Pro-> ${BLUE}[16S]${NC} -- tRNAs -- ${GREEN}[23S]${NC} -- ${ORANGE}[5S]${NC} <-Rev"
echo ""
echo -e "${NC}"
echo -e "${MAGENTA}Version: ${VERSION}${NC}"
echo -e "${MAGENTA}Step: ${STEP}${NC}"
echo -e "${MAGENTA}Threads: ${THREADS}${NC}"
echo ""

# Source the Python helper script to parse YAML
parse_config() {
    python3 << EOF
import yaml
import sys

try:
    with open('$CONFIG_FILE', 'r') as f:
        config = yaml.safe_load(f)

    # Output as shell variables
    print(f"POD5_DIR='{config['input'].get('pod5_dir', '')}'")
    print(f"FASTQ_DIR='{config['input'].get('fastq_dir', '')}'")
    print(f"OUTPUT_DIR='{config['output']['output_dir']}'")
    # Basecalling parameters
    print(f"BASECALL_MODEL='{config['basecalling'].get('model', 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0')}'")
    print(f"BASECALL_DEVICE='{config['basecalling'].get('device', 'cuda:all')}'")
    print(f"BASECALL_NO_TRIM='{str(config['basecalling'].get('no_trim', False)).lower()}'")
    print(f"BASECALL_BATCH_SIZE={config['basecalling'].get('batch_size', 0)}")
    print(f"BASECALL_RECURSIVE='{str(config['basecalling'].get('recursive', True)).lower()}'")
    print(f"BARCODE_KIT='{config['demultiplexing'].get('barcode_kit', 'SQK-RBK114-24')}'")
    print(f"BARCODE_THRESHOLD={config['demultiplexing'].get('barcode_threshold', 85)}")
    print(f"DISCARD_UNTRIMMED='{config['demultiplexing'].get('discard_untrimmed', 'false')}'")
    print(f"MIN_LENGTH={config['filtering']['min_length']}")
    print(f"MAX_LENGTH={config['filtering']['max_length']}")
    print(f"MIN_QUALITY={config['filtering']['min_quality']}")
    print(f"CLUSTER_TOOL='{config['clustering']['tool']}'")
    print(f"CLUSTER_IDENTITY={config['clustering'].get('cluster_identity', 0.97)}")
    print(f"SWARM_D={config['clustering'].get('swarm_d', 1)}")
    print(f"SWARM_FASTIDIOUS='{config['clustering'].get('swarm_fastidious', 'true')}'")
    boundary = config['clustering'].get('swarm_boundary', 3)
    if boundary is None or boundary == '' or str(boundary).lower() == 'none':
        boundary = 3
    print(f"SWARM_BOUNDARY={boundary}")
    print(f"SWARM_NO_OTU_BREAKING='{config['clustering'].get('swarm_no_otu_breaking', 'false')}'")
    print(f"CHIMERA_CHECK='{config['clustering'].get('remove_chimeras', 'true')}'")
    print(f"USE_MEDAKA='{config.get('polishing', {}).get('use_medaka', 'false')}'")
    print(f"MEDAKA_MODEL='{config.get('polishing', {}).get('medaka_model', 'r941_min_high_g360')}'")
    print(f"POLISH_JOBS={config.get('polishing', {}).get('polish_jobs', 15)}")
    print(f"POLISH_THREADS_PER_JOB={config.get('polishing', {}).get('threads_per_job', 4)}")
    print(f"POLISH_MEDAKA_BATCH={config.get('polishing', {}).get('medaka_batch_size', 200)}")
    print(f"BLAST_DB='{config['taxonomy']['blast_db']}'")
    print(f"BLAST_EVALUE={config['taxonomy']['evalue']}")
    print(f"BLAST_MAX_TARGETS={config['taxonomy']['max_target_seqs']}")
    print(f"BLAST_MIN_IDENTITY={config['taxonomy'].get('min_identity', 80)}")
except Exception as e:
    print(f"Error parsing config: {e}", file=sys.stderr)
    sys.exit(1)
EOF
}

# Load configuration
eval "$(parse_config)"

# Create output directory structure
mkdir -p "$OUTPUT_DIR"/{basecalls,demultiplexed,filtered,qc_before_filter,qc_after_filter,clustered,polished,taxonomy,logs,metadata}

# Export paths
BASECALL_DIR="$OUTPUT_DIR/basecalls"
DEMUX_DIR="$OUTPUT_DIR/demultiplexed"
FILTERED_DIR="$OUTPUT_DIR/filtered"
QC_BEFORE_FILTER_DIR="$OUTPUT_DIR/qc_before_filter"
QC_AFTER_FILTER_DIR="$OUTPUT_DIR/qc_after_filter"
CLUSTER_DIR="$OUTPUT_DIR/clustered"
POLISH_DIR="$OUTPUT_DIR/polished"
TAXONOMY_DIR="$OUTPUT_DIR/taxonomy"
LOG_DIR="$OUTPUT_DIR/logs"
METADATA_DIR="$OUTPUT_DIR/metadata"

# Function to truncate FASTA headers
truncate_fasta_headers() {
    local input_dir="$1"
    local logfile="${input_dir}/truncate_fasta_headers.log"
    
    log_info "Truncating FASTA headers in: $input_dir"
    echo "=== FASTA Header Truncation Log ($(date)) ===" > "$logfile"
    
    # Function to process individual FASTA file
    process_fasta() {
        local infile="$1"
        local logfile="$2"
        local tmpfile=$(mktemp)
        
        if [[ "$infile" == *.gz ]]; then
            gunzip -c "$infile" | awk -v logfile="$logfile" '
                /^>/ {
                    orig = substr($0, 2)
                    short = substr(orig, 1, 36)
                    printf(">%s\n", short)
                    printf("  Header truncated: \"%s\" → \"%s\"\n", orig, short) >> logfile
                    next
                }
                { print }
            ' > "$tmpfile"
            gzip -c "$tmpfile" > "${infile}.tmp"
            mv "${infile}.tmp" "$infile"
        else
            awk -v logfile="$logfile" '
                /^>/ {
                    orig = substr($0, 2)
                    short = substr(orig, 1, 36)
                    printf(">%s\n", short)
                    printf("  Header truncated: \"%s\" → \"%s\"\n", orig, short) >> logfile
                    next
                }
                { print }
            ' "$infile" > "${infile}.tmp"
            mv "${infile}.tmp" "$infile"
        fi
        rm -f "$tmpfile"
    }
    
    export -f process_fasta
    
    # Find and count FASTA files
    local fasta_count=0
    fasta_count=$(find "$input_dir" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fasta.gz" -o -iname "*.fa.gz" \) | wc -l)

    if [[ $fasta_count -eq 0 ]]; then
        log_warn "No FASTA files found in $input_dir"
        return 0
    fi

    # Process files using while read loop
    while IFS= read -r fasta_file; do
        process_fasta "$fasta_file" "$logfile"
    done < <(find "$input_dir" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fasta.gz" -o -iname "*.fa.gz" \))

    log_info "FASTA header truncation completed. Log: $logfile"
}


#################################################################
# FUNCTION: Generate Sample Metadata
#################################################################
generate_metadata() {
    log_step "Generating sample metadata file..."
    
    local input_dir=""
    
    # Determine which directory to scan for samples
    # Check filtered directory first (FASTA or FASTQ files)
    if [[ -d "$FILTERED_DIR/fasta" ]] && [[ -n "$(find "$FILTERED_DIR/fasta" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" \) 2>/dev/null)" ]]; then
        input_dir="$FILTERED_DIR/fasta"
        log_info "Using FASTA files from: $input_dir"
    elif [[ -d "$FILTERED_DIR" ]] && [[ -n "$(find "$FILTERED_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
        input_dir="$FILTERED_DIR"
        log_info "Using FASTQ files from: $input_dir"
    elif [[ -d "$DEMUX_DIR/barcodes" ]] && [[ -n "$(find "$DEMUX_DIR/barcodes" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
        input_dir="$DEMUX_DIR/barcodes"
        log_info "Using demultiplexed files from: $input_dir"
    elif [[ -n "$FASTQ_DIR" ]] && [[ -d "$FASTQ_DIR" ]] && [[ -n "$(find "$FASTQ_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) 2>/dev/null)" ]]; then
        input_dir="$FASTQ_DIR"
        log_info "Using input FASTQ files from: $input_dir"
    else
        log_error "No input directory with FASTQ/FASTA files found for metadata generation"
        log_error "Checked directories:"
        log_error "  - $FILTERED_DIR/fasta"
        log_error "  - $FILTERED_DIR"
        log_error "  - $DEMUX_DIR/barcodes"
        log_error "  - $FASTQ_DIR"
        return 1
    fi
    
    log_info "Scanning directory: $input_dir"
    
    # Check if input_dir is empty (safety check)
    if [[ -z "$input_dir" ]]; then
        log_error "Input directory is empty - this should not happen"
        return 1
    fi
    
    python3 << PYEOF
import os
import sys
from pathlib import Path
import re

input_dir = '$input_dir'
output_file = '$METADATA_DIR/sample_metadata.tsv'

if not os.path.exists(input_dir):
    print(f'Error: Input directory not found: {input_dir}', file=sys.stderr)
    sys.exit(1)

# Find all FASTQ and FASTA files
sample_files = []
for ext in ['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz', '*.fasta', '*.fa']:
    sample_files.extend(Path(input_dir).glob(ext))

if not sample_files:
    print(f'Error: No FASTQ/FASTA files found in {input_dir}', file=sys.stderr)
    sys.exit(1)

# Create metadata
metadata = []
metadata.append(['sample-id', 'barcode', 'description', 'group'])

for sample_file in sorted(sample_files):
    # Remove all possible extensions
    filename = sample_file.stem
    for ext in ['.fastq', '.fq', '.fasta', '.fa', '.gz', '_filtered']:
        filename = filename.replace(ext, '')
    
    # Extract barcode information
    barcode_match = re.search(r'(barcode|BC)(\d+)', filename, re.IGNORECASE)
    if barcode_match:
        barcode = f'{barcode_match.group(1)}{barcode_match.group(2)}'
        sample_id = filename
    else:
        barcode = 'unclassified'
        sample_id = filename
    
    # Basic description
    description = f'Sample from {filename}'
    group = 'default'
    
    metadata.append([sample_id, barcode, description, group])

# Write metadata file
os.makedirs(os.path.dirname(output_file), exist_ok=True)
with open(output_file, 'w') as f:
    for row in metadata:
        f.write('\t'.join(row) + '\n')

print(f'Generated metadata for {len(metadata)-1} samples')
PYEOF
    
    if [[ $? -eq 0 ]]; then
        log_success "Metadata file created: $METADATA_DIR/sample_metadata.tsv"
    else
        log_error "Failed to generate metadata file"
        return 1
    fi
}

#################################################################
# STEP 1: BASECALLING WITH DORADO
#################################################################
step_basecall() {
    log_step "STEP 1: Basecalling with Dorado..."
    
    # Check if dorado is installed
    if ! command -v dorado &> /dev/null; then
        log_error "dorado is not installed or not in PATH"
        log_error "Install dorado: https://github.com/nanoporetech/dorado"
        log_error "Or download models: dorado download --model MODEL_NAME"
        return 1
    fi
    
    log_info "Dorado version: $(dorado --version 2>&1 | head -n1 || echo 'unknown')"
    
    # Validate POD5 directory
    if [[ -z "$POD5_DIR" ]]; then
        log_error "POD5 directory not set in config.yaml"
        log_error "Please set input.pod5_dir: /path/to/pod5/files"
        return 1
    fi
    
    if [[ ! -d "$POD5_DIR" ]]; then
        log_error "POD5 directory not found: $POD5_DIR"
        return 1
    fi
    
    # Count POD5 files
    local pod5_count=$(find "$POD5_DIR" -type f -name "*.pod5" | wc -l)
    if [[ $pod5_count -eq 0 ]]; then
        log_error "No POD5 files found in: $POD5_DIR"
        log_error "Please check the directory contains .pod5 files"
        return 1
    fi
    
    log_info "Found $pod5_count POD5 files to basecall"
    
    # Check if output already exists
    if [[ -f "$BASECALL_DIR/calls.fastq" ]] && [[ -s "$BASECALL_DIR/calls.fastq" ]]; then
        local existing_reads=$(grep -c "^@" "$BASECALL_DIR/calls.fastq" 2>/dev/null || true); existing_reads=${existing_reads:-0}
        if [[ $existing_reads -gt 0 ]]; then
            log_warn "Basecalling output already exists: $BASECALL_DIR/calls.fastq"
            log_warn "  Existing reads: $existing_reads"
            log_info "Skipping basecalling step..."
            return 0
        fi
    fi
    
    log_info "Basecalling parameters:"
    log_info "  Model: $BASECALL_MODEL"
    log_info "  Device: $BASECALL_DEVICE"
    log_info "  Trim adapters: $([ "$BASECALL_NO_TRIM" == "true" ] && echo "No" || echo "Yes")"
    log_info "  Output: $BASECALL_DIR/calls.fastq"
    
    # Build Dorado command
    local dorado_cmd="dorado basecaller"
    
    # Add model
    dorado_cmd="$dorado_cmd $BASECALL_MODEL"
    
    # Add POD5 directory
    dorado_cmd="$dorado_cmd $POD5_DIR"
    
    # Add device
    dorado_cmd="$dorado_cmd --device $BASECALL_DEVICE"
    
    # Add trim option
    if [[ "$BASECALL_NO_TRIM" == "true" ]]; then
        dorado_cmd="$dorado_cmd --no-trim"
    fi
    
    # Add recursive search
    dorado_cmd="$dorado_cmd --recursive"
    
    # Add emit-fastq
    dorado_cmd="$dorado_cmd --emit-fastq"
    
    log_info "Running command: $dorado_cmd"
    log_info "This may take several hours depending on data size..."
    
    # Run dorado
    eval "$dorado_cmd" > "$BASECALL_DIR/calls.fastq" 2> "$LOG_DIR/dorado.log"
    
    local exit_code=$?
    
    if [[ $exit_code -ne 0 ]]; then
        log_error "Dorado basecalling failed with exit code: $exit_code"
        log_error "Check log: $LOG_DIR/dorado.log"
        
        # Show last 20 lines of log
        if [[ -f "$LOG_DIR/dorado.log" ]]; then
            log_error "Last 20 lines of log:"
            tail -n 20 "$LOG_DIR/dorado.log" >&2
        fi
        return 1
    fi
    
    # Validate output
    if [[ ! -f "$BASECALL_DIR/calls.fastq" ]]; then
        log_error "Output file not created: $BASECALL_DIR/calls.fastq"
        return 1
    fi
    
    if [[ ! -s "$BASECALL_DIR/calls.fastq" ]]; then
        log_error "Output file is empty: $BASECALL_DIR/calls.fastq"
        return 1
    fi
    
    # Count reads
    local read_count=$(grep -c "^@" "$BASECALL_DIR/calls.fastq" 2>/dev/null || true); read_count=${read_count:-0}
    local file_size=$(du -h "$BASECALL_DIR/calls.fastq" | cut -f1)
    
    log_success "Basecalling completed successfully!"
    log_info "  Output: $BASECALL_DIR/calls.fastq"
    log_info "  Reads: $read_count"
    log_info "  Size: $file_size"
    log_info "  Log: $LOG_DIR/dorado.log"
}

#################################################################
# STEP 2: DEMULTIPLEXING WITH PORECHOP
#################################################################
step_demultiplex() {
    log_step "STEP 2: Demultiplexing with Porechop..."
    
    local input_fastq=""
    
    # First check for calls.fastq in basecall directory
    if [[ -f "$BASECALL_DIR/calls.fastq" ]]; then
        input_fastq="$BASECALL_DIR/calls.fastq"
        log_info "Using basecalled FASTQ: $input_fastq"
    # Then check if FASTQ_DIR is specified
    elif [[ -n "$FASTQ_DIR" ]]; then
        # Check if FASTQ_DIR is a file
        if [[ -f "$FASTQ_DIR" ]]; then
            input_fastq="$FASTQ_DIR"
            log_info "Using FASTQ file: $input_fastq"
        # Check if FASTQ_DIR is a directory
        elif [[ -d "$FASTQ_DIR" ]]; then
            # Find any FASTQ file in the directory
            local found_fastq=$(find "$FASTQ_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz" \) | head -1)
            if [[ -n "$found_fastq" ]]; then
                input_fastq="$found_fastq"
                log_info "Using FASTQ file from directory: $input_fastq"
            else
                log_error "No FASTQ files found in directory: $FASTQ_DIR"
                return 1
            fi
        else
            log_error "FASTQ_DIR does not exist: $FASTQ_DIR"
            return 1
        fi
    else
        log_error "No input FASTQ specified. Set fastq_dir in config or run basecalling first."
        return 1
    fi
    
    if [[ -d "$DEMUX_DIR/barcodes" ]] && [[ $(ls -A "$DEMUX_DIR/barcodes" 2>/dev/null) ]]; then
        log_warn "Demultiplexing output already exists. Skipping..."
        return 0
    fi
    
    mkdir -p "$DEMUX_DIR/barcodes"
    
    log_info "Demultiplexing with kit: $BARCODE_KIT"
    log_info "Barcode threshold: $BARCODE_THRESHOLD"
    log_info "Using $THREADS threads"

    
    porechop \
        -i "$input_fastq" \
        -b "$DEMUX_DIR/barcodes" \
        --barcode_threshold "$BARCODE_THRESHOLD" \
        --format fastq.gz \
        -t "$THREADS" \
        --discard_unassigned \
        --discard_middle \
        > "$LOG_DIR/porechop.log" 2>&1
    
    # Count barcodes
    local barcode_count=$(ls "$DEMUX_DIR/barcodes"/*.fastq.gz 2>/dev/null | wc -l)
    log_success "Demultiplexing completed: $barcode_count barcodes found"
}

#################################################################
# STEP 3: FILTERING WITH CHOPPER
#################################################################
#################################################################
# STEP 2.5: QC BEFORE FILTERING
#################################################################
step_qc_before_filter() {
    log_step "STEP 2.5: Quality assessment BEFORE filtering with NanoPlot..."

    if [[ -d "$QC_BEFORE_FILTER_DIR" ]] && [[ -f "$QC_BEFORE_FILTER_DIR/NanoStats.txt" ]]; then
        log_warn "Pre-filter QC outputs already exist in: $QC_BEFORE_FILTER_DIR"
        log_info "Skipping pre-filter QC step..."
        return 0
    fi

    local input_dir=""
    if [[ -d "$DEMUX_DIR/barcodes" ]] && [[ -n "$(ls -A "$DEMUX_DIR/barcodes" 2>/dev/null)" ]]; then
        input_dir="$DEMUX_DIR/barcodes"
    elif [[ -n "$FASTQ_DIR" ]] && [[ -d "$FASTQ_DIR" ]]; then
        input_dir="$FASTQ_DIR"
    else
        log_error "No input directory found for pre-filter QC"
        return 1
    fi

    log_info "Running pre-filter QC on: $input_dir"
    log_info "Using $THREADS threads"

    find "$input_dir" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) | \
    parallel -j "$THREADS" '
        raw_file="{}"
        barcode_name=$(basename "$raw_file" | sed "s/.fastq.gz//g" | sed "s/.fastq//g" | sed "s/.fq.gz//g" | sed "s/.fq//g")
        qc_dir="'"$QC_BEFORE_FILTER_DIR"'/$barcode_name"

        if [[ -d "$qc_dir" ]] && [[ -f "$qc_dir/NanoStats.txt" ]]; then
            echo -e "'"${YELLOW}"'[SKIP]'"${NC}"' Pre-filter QC for $barcode_name already exists"
            exit 0
        fi

        echo -e "'"${GREEN}"'[QC]'"${NC}"' Pre-filter QC for $barcode_name..."
        mkdir -p "$qc_dir"

        NanoPlot \
            -t 2 \
            --fastq "$raw_file" \
            -o "$qc_dir" \
            --tsv_stats \
            --info_in_report \
            --prefix "${barcode_name}_" \
            > "'"$LOG_DIR"'/nanoplot_prefilter_${barcode_name}.log" 2>&1

        echo -e "'"${GREEN}"'[DONE]'"${NC}"' Pre-filter QC completed for $barcode_name"
    '

    multiqc "$QC_BEFORE_FILTER_DIR" -o "$QC_BEFORE_FILTER_DIR/multiqc_report" --force \
        > "$LOG_DIR/multiqc_prefilter.log" 2>&1 || log_warn "MultiQC failed for pre-filter QC"

    log_success "Pre-filter QC completed: $QC_BEFORE_FILTER_DIR"
}

step_filter() {
    log_step "STEP 3: Filtering reads with Chopper..."
    
    local input_dir=""
    if [[ -d "$DEMUX_DIR/barcodes" ]]; then
        input_dir="$DEMUX_DIR/barcodes"
    elif [[ -n "$FASTQ_DIR" ]] && [[ -d "$FASTQ_DIR" ]]; then
        input_dir="$FASTQ_DIR"
    else
        log_error "No input directory found for filtering"
        return 1
    fi
    
    log_info "Filtering parameters:"
    log_info "  Min length: $MIN_LENGTH bp"
    log_info "  Max length: $MAX_LENGTH bp"
    log_info "  Min quality: $MIN_QUALITY"
    log_info "Using $THREADS threads"
    
    # Count input files
    local input_count=$(find "$input_dir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) | wc -l)
    log_info "Found $input_count FASTQ files to process"
    
    if [[ $input_count -eq 0 ]]; then
        log_error "No FASTQ files found in $input_dir"
        return 1
    fi

    local filtered_count=0
    
    # Process files in parallel
    find "$input_dir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" -o -name "*.fq" -o -name "*.fq.gz" \) | \
    parallel -j "$THREADS" '
        FASTQ="{}"
        barcode_name=$(basename "$FASTQ" | sed "s/.fastq.gz//g" | sed "s/.fastq//g" | sed "s/.fq.gz//g" | sed "s/.fq//g")
        OUTPUT_FILE="'"$FILTERED_DIR"'/${barcode_name}_filtered.fastq"
        
        if [[ -f "$OUTPUT_FILE" ]]; then
            echo -e "'"${YELLOW}"'[SKIP]'"${NC}"' $barcode_name already filtered"
            exit 0
        fi
        
        echo -e "'"${GREEN}"'[FILTERING]'"${NC}"' $barcode_name..."
        
        # Run Chopper with correct syntax
        chopper \
            -l '"$MIN_LENGTH"' \
            --maxlength '"$MAX_LENGTH"' \
            -q '"$MIN_QUALITY"' \
            -t 1 \
            -i "$FASTQ" \
            > "$OUTPUT_FILE" \
            2>> "'"$LOG_DIR"'/chopper_${barcode_name}.log"
        
        read_count=$(grep -c "^@" "$OUTPUT_FILE" 2>/dev/null || true); read_count=${read_count:-0}
        echo -e "'"${GREEN}"'[DONE]'"${NC}"' $barcode_name: $read_count reads"
    '
    
    # Count total filtered reads
    local total_reads=0
    while IFS= read -r filtered_file; do
        reads=$(grep -c "^@" "$filtered_file" 2>/dev/null || true); reads=${reads:-0}
        total_reads=$((total_reads + reads))
        filtered_count=$((filtered_count + 1))
    done < <(find "$FILTERED_DIR" -maxdepth 1 -name "*_filtered.fastq" -type f)
    
    log_success "Filtering completed: $filtered_count samples, $total_reads total reads"
    
# Convert FASTQ to FASTA using seqtk
    log_info "Converting filtered FASTQ files to FASTA format..."
    
    # Create fasta directory
    mkdir -p "$FILTERED_DIR/fasta"
    
    # Verify directory was created
    if [[ ! -d "$FILTERED_DIR/fasta" ]]; then
        log_error "Failed to create directory: $FILTERED_DIR/fasta"
        log_error "Check write permissions on: $FILTERED_DIR"
        return 1
    fi
    
    log_info "FASTA output directory: $FILTERED_DIR/fasta"
    
    # Check if seqtk is available
    if ! command -v seqtk &> /dev/null; then
        log_error "seqtk is not installed or not in PATH"
        log_error "Please install seqtk: conda install -c bioconda seqtk"
        return 1
    fi
    
    log_info "Using seqtk: $(which seqtk)"
    
    # Count FASTQ files using find
    local total_files=$(find "$FILTERED_DIR" -maxdepth 1 -name "*_filtered.fastq" -type f | wc -l)
    
    if [[ $total_files -eq 0 ]]; then
        log_error "No *_filtered.fastq files found in $FILTERED_DIR"
        log_error "Available files:"
        ls -lh "$FILTERED_DIR"/ | head -10
        return 1
    fi
    
    log_info "Found $total_files filtered FASTQ files to convert"
    
    # Convert each file
    local fasta_count=0
    local failed_count=0
    local file_num=0
    
    while IFS= read -r FASTQ; do
        file_num=$((file_num + 1))
        
        # Extract filename without path and extension
        full_name=$(basename "$FASTQ")
        filename="${full_name%.fastq}"
        FASTA="$FILTERED_DIR/fasta/${filename}.fasta"
        
        log_info "[$file_num/$total_files] Converting $filename to FASTA..."
        
        # Run seqtk conversion
        if seqtk seq -a "$FASTQ" > "$FASTA" 2>> "$LOG_DIR/seqtk_convert.log"; then
            # Check if file was created and has content
            if [[ -s "$FASTA" ]]; then
                seq_count=$(grep -c "^>" "$FASTA" 2>/dev/null || true); seq_count=${seq_count:-0}
                if [[ $seq_count -gt 0 ]]; then
                    log_success "  → Created with $seq_count sequences"
                    fasta_count=$((fasta_count + 1))
                else
                    log_error "  → Empty or has no sequences"
                    failed_count=$((failed_count + 1))
                fi
            else
                log_error "  → Failed to create file"
                failed_count=$((failed_count + 1))
            fi
        else
            log_error "  → seqtk conversion failed"
            failed_count=$((failed_count + 1))
        fi
    done < <(find "$FILTERED_DIR" -maxdepth 1 -name "*_filtered.fastq" -type f)
    
    log_info "Conversion summary: $fasta_count succeeded, $failed_count failed"
    
    # Final verification
    if [[ $fasta_count -eq 0 ]]; then
        log_error "No FASTA files were generated successfully"
        log_error "Check conversion log: $LOG_DIR/seqtk_convert.log"
        return 1
    fi
    
    # Truncate FASTA headers to 36 characters
    log_info "Truncating FASTA headers to 36 characters..."
    
    local truncated=0
    
    while IFS= read -r FASTA; do
        temp_fasta="${FASTA}.truncating"
        
        # Truncate headers using awk
        awk '
        BEGIN { FS="" }
        /^>/ {
            header = substr($0, 2)
            if (length(header) > 36) {
                header = substr(header, 1, 36)
            }
            print ">" header
            next
        }
        { print }
        ' "$FASTA" > "$temp_fasta"
        
        # Replace original with truncated version
        if [[ -s "$temp_fasta" ]]; then
            mv "$temp_fasta" "$FASTA"
            truncated=$((truncated + 1))
        else
            rm -f "$temp_fasta"
        fi
    done < <(find "$FILTERED_DIR/fasta" -name "*.fasta" -type f)
    
    log_info "Truncated headers in $truncated FASTA files"
    
    # List generated files (first 5)
    log_info "Generated FASTA files:"
    find "$FILTERED_DIR/fasta" -name "*.fasta" -type f | head -5 | while IFS= read -r file; do
        size=$(du -h "$file" | cut -f1)
        log_info "  $(basename "$file") ($size)"
    done
    
    if [[ $fasta_count -gt 5 ]]; then
        log_info "  ... and $(($fasta_count - 5)) more files"
    fi
    
    log_success "Filtering and conversion completed"
    log_info "  - Filtered FASTQ files: $filtered_count samples, $total_reads total reads"
    log_info "  - FASTA files generated: $fasta_count files"
    log_info "  - FASTA location: $FILTERED_DIR/fasta/"
}


#################################################################
# STEP 4: QUALITY ASSESSMENT WITH NANOPLOT
#################################################################
step_qc() {
    log_step "STEP 4: Quality assessment AFTER filtering with NanoPlot..."
    
    # Check if QC outputs already exist
    if [[ -d "$QC_AFTER_FILTER_DIR" ]] && [[ -f "$QC_AFTER_FILTER_DIR/NanoStats.txt" ]]; then
        log_warn "QC outputs already exist in: $QC_AFTER_FILTER_DIR"
        log_info "Skipping QC step..."
        return 0
    fi
    
    local input_dir="$FILTERED_DIR"
    if [[ ! -d "$input_dir" ]] || [[ ! $(ls -A "$input_dir" 2>/dev/null) ]]; then
        log_error "No filtered reads found for QC"
        return 1
    fi
    
    log_info "Using $THREADS threads"
    
    # Process QC in parallel
    find "$input_dir" -maxdepth 1 -type f \( -name "*_filtered.fastq" -o -name "*_filtered.fastq.gz" \) | \
    parallel -j "$THREADS" '
        filtered_file="{}"
        barcode_name=$(basename "$filtered_file" | sed "s/_filtered.fastq.gz//g" | sed "s/_filtered.fastq//g")
        qc_dir="'"$QC_AFTER_FILTER_DIR"'/$barcode_name"
        
        if [[ -d "$qc_dir" ]] && [[ -f "$qc_dir/NanoStats.txt" ]]; then
            echo -e "'"${YELLOW}"'[SKIP]'"${NC}"' QC for $barcode_name already exists"
            exit 0
        fi
        
        echo -e "'"${GREEN}"'[QC]'"${NC}"' Generating plots for $barcode_name..."
        
        mkdir -p "$qc_dir"
        
        NanoPlot \
            -t 2 \
            --fastq "$filtered_file" \
            -o "$qc_dir" \
            --tsv_stats \
            --info_in_report \
            --prefix "${barcode_name}_" \
            > "'"$LOG_DIR"'/nanoplot_${barcode_name}.log" 2>&1
        
        echo -e "'"${GREEN}"'[DONE]'"${NC}"' QC completed for $barcode_name"
    '
    
    # Generate MultiQC report
    log_info "Generating MultiQC report..."
    multiqc "$QC_AFTER_FILTER_DIR" -o "$QC_AFTER_FILTER_DIR/multiqc_report" --force > "$LOG_DIR/multiqc.log" 2>&1 || log_warn "MultiQC failed"
    
    log_success "QC assessment completed"
}

#################################################################
# STEP 5: BLAST ON INDIVIDUAL READS PER BARCODE
#################################################################
step_blast_reads() {
    log_step "STEP 5: BLAST ON ALL READS (COMBINED)..."

    mkdir -p "$TAXONOMY_DIR/reads"

    # Skip if all output files already exist and are non-empty
    local all_files_exist=true
    local required_files=(
        "$TAXONOMY_DIR/reads/blast_results_reads.txt"
        "$TAXONOMY_DIR/reads/taxonomy_reads.tsv"
        "$TAXONOMY_DIR/reads/taxonomy_reads_hits_only.tsv"
        "$TAXONOMY_DIR/reads/taxonomy_reads_expanded.tsv"
        "$TAXONOMY_DIR/reads/taxonomy_reads_for_microeco.tsv"
        "$TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt"
        "$FILTERED_DIR/reads_per_sample.txt"
    )

    while IFS= read -r req_file; do
        if [[ ! -f "$req_file" ]] || [[ ! -s "$req_file" ]]; then
            all_files_exist=false
            break
        fi
    done < <(printf '%s\n' "${required_files[@]}")

    if [[ "$all_files_exist" == "true" ]]; then
        log_warn "All read-level taxonomy files already exist — skipping"
        log_success "Read-level taxonomy assignment completed (using existing files)"
        return 0
    fi

    # Validate prerequisites
    if [[ ! -d "$FILTERED_DIR" ]]; then
        log_error "Filtered reads directory not found: $FILTERED_DIR"
        log_error "Run the filter step first."
        return 1
    fi

    if [[ ! -f "${BLAST_DB}.nhr" ]] && [[ ! -f "${BLAST_DB}.ndb" ]]; then
        log_error "BLAST database not found: $BLAST_DB"
        return 1
    fi

    # Locate FASTA directory
    local search_dir="$FILTERED_DIR"
    if [[ -d "$FILTERED_DIR/fasta" ]]; then
        search_dir="$FILTERED_DIR/fasta"
    fi
    log_info "Using FASTA directory: $search_dir"

    # Step 5.1: Concatenate all barcode FASTAs into one file
    log_info "Step 5.1: Combining all barcode reads into a single FASTA..."

    local combined_fasta="$TAXONOMY_DIR/reads/all_reads_combined.fasta"
    local barcode_count=0
    local total_read_count=0
    > "$combined_fasta"

    while IFS= read -r fasta_file; do
        local barcode
        barcode=$(basename "$fasta_file" | sed 's/_filtered\.fasta$//;s/_filtered\.fa$//;s/\.fasta$//;s/\.fa$//')
        local reads
        reads=$(grep -c "^>" "$fasta_file" 2>/dev/null || true); reads=${reads:-0}
        if [[ $reads -gt 0 ]]; then
            cat "$fasta_file" >> "$combined_fasta"
            total_read_count=$((total_read_count + reads))
            barcode_count=$((barcode_count + 1))
            log_info "  Added $barcode: $reads reads"
        else
            log_warn "  Skipping $barcode: no reads"
        fi
    done < <(find "$search_dir" -maxdepth 1 -type f \( -name "*_filtered.fasta" -o -name "*_filtered.fa" -o -name "*.fasta" -o -name "*.fa" \) | sort)

    if [[ $barcode_count -eq 0 ]] || [[ ! -s "$combined_fasta" ]]; then
        log_error "No reads found to combine in $search_dir"
        return 1
    fi

    log_success "Combined $total_read_count reads from $barcode_count barcodes → $combined_fasta"

    # Step 5.2: Run single BLAST on combined FASTA
    log_info "Step 5.2: Running BLAST on combined reads..."

    blastn \
        -query "$combined_fasta" \
        -db "$BLAST_DB" \
        -task megablast \
        -reward 2 -penalty -3 -word_size 65 -gapopen 0 -gapextend 4 \
        -out "$TAXONOMY_DIR/reads/blast_results_reads.txt" \
        -evalue "$BLAST_EVALUE" \
        -max_target_seqs "$BLAST_MAX_TARGETS" \
        -num_threads "$THREADS" \
        -max_hsps 1 \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        > "$LOG_DIR/blast_reads.log" 2>&1

    if [[ $? -ne 0 ]]; then
        log_error "BLAST failed — see $LOG_DIR/blast_reads.log"
        return 1
    fi

    local total_hits
    total_hits=$(wc -l < "$TAXONOMY_DIR/reads/blast_results_reads.txt" 2>/dev/null || true)
    total_hits=${total_hits:-0}
    log_success "BLAST complete: $total_hits hits across $total_read_count reads"
    
    # Step 5.3: Create reads per sample table
    log_info "Step 5.3: Creating reads per sample table..."
    
    if [[ -f "$FILTERED_DIR/reads_per_sample.txt" ]]; then
        log_info "Reads per sample table already exists: $FILTERED_DIR/reads_per_sample.txt"
    else
        export FILTERED_DIR
        
        python3 << 'PYEOF7'
import pandas as pd
import os
from pathlib import Path

filtered_dir = os.environ.get('FILTERED_DIR', '')

try:
    # Check for fasta subdirectory (common structure: filtered/fasta/)
    fasta_dir = Path(filtered_dir) / "fasta"
    if fasta_dir.exists():
        search_dir = fasta_dir
        print(f"Looking for filtered files in: {search_dir}")
    else:
        search_dir = Path(filtered_dir)
        print(f"Looking for filtered files in: {search_dir}")
    
    # Find all filtered fasta files - try multiple patterns
    filtered_files = []
    
    # Try different file patterns
    for pattern in ["*_filtered.fasta", "*_filtered.fa", "*.fasta", "*.fa"]:
        files = sorted(search_dir.glob(pattern))
        if files:
            print(f"Found {len(files)} files matching pattern: {pattern}")
            filtered_files = files
            break
    
    if not filtered_files:
        print(f"ERROR: No FASTA files found in {filtered_dir}")
        print(f"Checked patterns: *_filtered.fasta, *_filtered.fa, *.fasta, *.fa")
        # List what's actually in the directory
        if Path(filtered_dir).exists():
            all_files = list(Path(filtered_dir).glob("*"))
            print(f"Files in directory ({len(all_files)} total):")
            for f in all_files[:10]:  # Show first 10
                print(f"  {f.name}")
            if len(all_files) > 10:
                print(f"  ... and {len(all_files) - 10} more")
        else:
            print(f"Directory does not exist: {filtered_dir}")
        import sys
        sys.exit(1)
    
    print(f"Processing {len(filtered_files)} FASTA files...")
    
    # Create a dictionary to store read counts per barcode
    read_data = {}
    
    for fasta_file in filtered_files:
        barcode = fasta_file.stem.replace("_filtered", "")
        
        # Read the fasta file and extract read IDs
        read_ids = []
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    read_id = line.strip()[1:].split()[0]
                    read_ids.append(read_id)
        
        # Add read IDs to dictionary
        for read_id in read_ids:
            if read_id not in read_data:
                read_data[read_id] = {}
            read_data[read_id][barcode] = 1
    
    # Create DataFrame
    reads_df = pd.DataFrame(read_data).T
    reads_df = reads_df.fillna(0).astype(int)
    
    # Sort by total read count (descending)
    reads_df['total'] = reads_df.sum(axis=1)
    reads_df = reads_df.sort_values('total', ascending=False)
    reads_df = reads_df.drop('total', axis=1)
    
    # Sort columns (barcodes) alphabetically
    reads_df = reads_df[sorted(reads_df.columns)]
    
    # Save to file
    output_file = os.path.join(filtered_dir, "reads_per_sample.txt")
    reads_df.to_csv(output_file, sep='\t')
    
    print(f"\nReads per sample table created: {output_file}")
    print(f"  Total unique reads: {len(reads_df)}")
    print(f"  Total samples (barcodes): {len(reads_df.columns)}")
    print(f"  Total read occurrences: {reads_df.sum().sum()}")
    
except Exception as e:
    print(f"Error creating reads per sample table: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)
PYEOF7
        
        if [[ $? -eq 0 ]]; then
            log_success "Reads per sample table created: $FILTERED_DIR/reads_per_sample.txt"
        else
            log_warn "Failed to create reads per sample table"
        fi
    fi
    
    # Step 5.4: Parse combined BLAST results
    log_info "Step 5.4: Parsing BLAST results with embedded taxonomy..."
    
    # Export variables for Python
    export TAXONOMY_DIR
    export FILTERED_DIR
    export BLAST_MIN_IDENTITY
    
    python3 << 'PYEOF8'
import sys
import pandas as pd
import os
from pathlib import Path

taxonomy_dir = os.environ.get("TAXONOMY_DIR", "")
filtered_dir = os.environ.get("FILTERED_DIR", "")
blast_file = os.path.join(taxonomy_dir, "reads", "blast_results_reads.txt")
output_file = os.path.join(taxonomy_dir, "reads", "taxonomy_reads.tsv")

# Column names for TSV format (format 6)
blast_columns = ["Read_ID", "Hit_with_Taxonomy", "Identity_%", "Alignment_Length", "E-value", "Bitscore"]

try:
    min_identity = float(os.environ.get("BLAST_MIN_IDENTITY", 80))
    # Read ALL Read IDs from the filtered fasta files
    # Check for fasta subdirectory (common structure: filtered/fasta/)
    fasta_dir = Path(filtered_dir) / "fasta"
    if fasta_dir.exists():
        search_dir = fasta_dir
        print(f"Reading all Read IDs from filtered files in: {search_dir}")
    else:
        search_dir = Path(filtered_dir)
        print(f"Reading all Read IDs from filtered files in: {search_dir}")
    
    all_reads = []
    
    # Try different file patterns
    filtered_files = []
    for pattern in ["*_filtered.fasta", "*_filtered.fa", "*.fasta", "*.fa"]:
        files = sorted(search_dir.glob(pattern))
        if files:
            print(f"Found {len(files)} files matching pattern: {pattern}")
            filtered_files = files
            break
    
    if not filtered_files:
        print(f"ERROR: No FASTA files found in {filtered_dir}")
        import sys
        sys.exit(1)
    
    for fasta_file in filtered_files:
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    read_id = line.strip()[1:].split()[0]
                    all_reads.append(read_id)
    
    print(f"Total reads from filtered files: {len(all_reads)}")
    
    # Check if BLAST results exist
    if not os.path.exists(blast_file):
        print(f"Error: BLAST results file not found: {blast_file}", file=sys.stderr)
        sys.exit(1)
    
    # Load BLAST results (TSV format - tab-separated)
    print(f"Loading BLAST results from: {blast_file}")
    blast_df = pd.read_csv(blast_file, sep="\t", names=blast_columns, header=None)
    
    # Convert numeric columns to proper types
    if not blast_df.empty:
        blast_df["Identity_%"] = pd.to_numeric(blast_df["Identity_%"], errors="coerce")
        blast_df["Alignment_Length"] = pd.to_numeric(blast_df["Alignment_Length"], errors="coerce")
        blast_df["E-value"] = pd.to_numeric(blast_df["E-value"], errors="coerce")
        blast_df["Bitscore"] = pd.to_numeric(blast_df["Bitscore"], errors="coerce")
        
        # Drop rows with invalid numeric values
        blast_df = blast_df.dropna(subset=["Identity_%", "Alignment_Length"])
        
        # Apply minimum identity threshold
        before = len(blast_df)
        blast_df = blast_df[blast_df["Identity_%"] >= min_identity]
        print(f"Hits after min_identity filter ({min_identity}%): {len(blast_df)} (removed {before - len(blast_df)})")
    
    if blast_df.empty:
        print("Warning: No BLAST hits found", file=sys.stderr)
        if all_reads:
            taxonomy = pd.DataFrame({
                "Read_ID": all_reads,
                "Best_Hit": "No hit",
                "Identity_%": 0.0,
                "Alignment_Length": 0,
                "E-value": 999.0,
                "Bitscore": 0.0,
                "Taxonomy": "Unassigned"
            })
            taxonomy.to_csv(output_file, sep="\t", index=False)
        sys.exit(0)
    
    print(f"Found {len(blast_df)} total BLAST hits")
    read_count = blast_df["Read_ID"].nunique()
    print(f"Reads with BLAST hits: {read_count}")
    
    # Extract taxonomy from sseqid column
    # Expected format: accession|taxonomy (e.g., ACC123|Bacteria;Proteobacteria;...)
    def extract_taxonomy(hit_string):
        """Extract taxonomy from BLAST subject sequence ID"""
        if pd.isna(hit_string) or hit_string == "":
            return "Unknown", "Unassigned"
        
        s = str(hit_string).strip()
        # Split by pipe to separate accession and taxonomy
        parts = s.split("|")
        
        if len(parts) >= 2:
            accession = parts[0].strip()
            taxonomy = "|".join(parts[1:]).strip()
            return accession, taxonomy
        # Fallback: no pipe but contains semicolons → treat whole string as taxonomy
        elif ";" in s:
            return "Unknown", s
        else:
            return s, "Unassigned"
    
    blast_df[["Best_Hit", "Taxonomy"]] = blast_df["Hit_with_Taxonomy"].apply(
        lambda x: pd.Series(extract_taxonomy(x))
    )
    
    # Keep only first hit per read (as it appears in BLAST output - no sorting/filtering)
    blast_df = blast_df.drop_duplicates(subset=["Read_ID"], keep="first")
    
    print(f"\nFirst hits kept for {len(blast_df)} reads")
    
    # Create final taxonomy table
    taxonomy_df = blast_df[["Read_ID", "Best_Hit", "Identity_%", "Alignment_Length", "E-value", "Bitscore", "Taxonomy"]].copy()
    
    # Add reads without hits
    reads_with_hits = set(taxonomy_df["Read_ID"])
    all_reads_set = set([str(read).strip() for read in all_reads])
    reads_without_hits = sorted(list(all_reads_set - reads_with_hits))
    
    print(f"\nRead Summary:")
    print(f"  Total reads in filtered files: {len(all_reads)}")
    print(f"  Unique reads: {len(all_reads_set)}")
    print(f"  Reads with BLAST hits: {len(reads_with_hits)}")
    print(f"  Reads without hits (to be marked as Unassigned): {len(reads_without_hits)}")
    
    if reads_without_hits:
        no_hits_df = pd.DataFrame({
            "Read_ID": reads_without_hits,
            "Best_Hit": "No hit",
            "Identity_%": 0.0,
            "Alignment_Length": 0,
            "E-value": 999.0,
            "Bitscore": 0.0,
            "Taxonomy": "Unassigned"
        })
        taxonomy_df = pd.concat([taxonomy_df, no_hits_df], ignore_index=True)
        print(f"  Added {len(reads_without_hits)} reads as unassigned")
    
    # Save taxonomy table (all reads - with and without hits)
    taxonomy_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nTaxonomy table saved: {output_file}")
    print(f"Total reads: {len(taxonomy_df)}")
    with_tax = (taxonomy_df["Taxonomy"] != "Unassigned").sum()
    print(f"  With taxonomy: {with_tax}")
    unassigned = (taxonomy_df["Taxonomy"] == "Unassigned").sum()
    print(f"  Unassigned: {unassigned}")
    
    # Create taxonomy file with ONLY reads that have BLAST hits (excluding Unassigned)
    output_file_hits_only = os.path.join(taxonomy_dir, "reads", "taxonomy_reads_hits_only.tsv")
    taxonomy_hits_only = taxonomy_df[taxonomy_df["Taxonomy"] != "Unassigned"].copy()
    taxonomy_hits_only.to_csv(output_file_hits_only, sep="\t", index=False)
    print(f"\nTaxonomy hits-only table saved: {output_file_hits_only}")
    print(f"  Reads with taxonomy hits: {len(taxonomy_hits_only)}")
    
    # Create expanded taxonomy file with separate rank columns
    output_file_expanded = os.path.join(taxonomy_dir, "reads", "taxonomy_reads_expanded.tsv")
    
    def expand_taxonomy(tax_string):
        """Expand taxonomy string into separate rank columns"""
        # Expected format: Kingdom;Phylum;Class;Order;Family;Genus;Species;Strain
        ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        
        if pd.isna(tax_string) or tax_string == "Unassigned":
            return pd.Series(["Unassigned"] * 8, index=ranks)
        
        # Split by semicolon
        parts = str(tax_string).split(";")
        
        # Pad with "Unclassified" if not enough ranks
        while len(parts) < 8:
            parts.append("Unclassified")
        
        # Take first 8 ranks
        parts = parts[:8]
        
        return pd.Series(parts, index=ranks)
    
    # Apply expansion only to hits (exclude Unassigned)
    expanded_df = taxonomy_hits_only.copy()
    taxonomy_ranks = expanded_df["Taxonomy"].apply(expand_taxonomy)
    
    # Create final dataframe with Read_ID and taxonomic ranks
    final_expanded = pd.DataFrame({
        "Read_ID": expanded_df["Read_ID"]
    })
    
    # Add each rank as a column
    for rank in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]:
        final_expanded[rank] = taxonomy_ranks[rank].values
    
    # Save expanded taxonomy
    final_expanded.to_csv(output_file_expanded, sep="\t", index=False)
    print(f"\nExpanded taxonomy table saved: {output_file_expanded}")
    print(f"  Columns: Read_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain")
    
    # Create microeco-compatible file (only Read_ID + 8 ranks, no other columns)
    output_file_microeco = os.path.join(taxonomy_dir, "reads", "taxonomy_reads_for_microeco.tsv")
    final_expanded.to_csv(output_file_microeco, sep="\t", index=False)
    print(f"\nMicroeco-compatible taxonomy table saved: {output_file_microeco}")
    print(f"  Format: Read_ID + 8 taxonomic rank columns only")
    print(f"  Total reads: {len(final_expanded)}")
    
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF8
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to parse BLAST results"
        return 1
    fi
    
    # Step 5.5: Create filtered read table
    log_info "Step 5.5: Creating filtered read table (only reads with taxonomy)..."
    
    if [[ -f "$TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt" ]]; then
        log_info "Filtered read table already exists: $TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt"
    else
        if [[ -f "$FILTERED_DIR/reads_per_sample.txt" ]]; then
            export FILTERED_DIR
            export TAXONOMY_DIR
            
            python3 << 'PYEOF9'
import pandas as pd
import os
import sys

filtered_dir = os.environ.get('FILTERED_DIR', '')
taxonomy_dir = os.environ.get('TAXONOMY_DIR', '')

read_table_file = os.path.join(filtered_dir, "reads_per_sample.txt")
taxonomy_file = os.path.join(taxonomy_dir, "reads", "taxonomy_reads.tsv")
output_file = os.path.join(taxonomy_dir, "reads", "reads_per_sample_with_taxonomy.txt")

try:
    # Read taxonomy file
    print(f"Reading taxonomy file: {taxonomy_file}")
    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")
    
    # Get reads with BLAST hits (not unassigned)
    reads_with_hits = taxonomy_df[taxonomy_df["Taxonomy"] != "Unassigned"]["Read_ID"].tolist()
    print(f"Found {len(reads_with_hits)} reads with BLAST hits")
    
    # Read read table
    print(f"Reading read table: {read_table_file}")
    read_table = pd.read_csv(read_table_file, sep="\t", index_col=0)
    
    total_reads = len(read_table)
    print(f"Total reads in table: {total_reads}")
    
    # Check if all reads already have hits
    reads_in_table = set(read_table.index)
    reads_without_hits = reads_in_table - set(reads_with_hits)
    
    if len(reads_without_hits) == 0:
        print("All reads in the table already have BLAST hits. No filtering needed.")
        sys.exit(0)
    
    print(f"Reads without BLAST hits: {len(reads_without_hits)}")
    
    # Filter table to keep only reads with hits
    filtered_table = read_table.loc[read_table.index.isin(reads_with_hits)]
    
    # Save filtered table
    filtered_table.to_csv(output_file, sep="\t")
    
    print(f"\nFiltered read table saved: {output_file}")
    print(f"  Original reads: {total_reads}")
    print(f"  With taxonomy: {len(filtered_table)}")
    print(f"  Removed (no hits): {total_reads - len(filtered_table)}")
    
except Exception as e:
    print(f"Error creating filtered read table: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF9
            
            if [[ $? -eq 0 ]]; then
                if [[ -f "$TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt" ]]; then
                    log_success "Filtered read table created: $TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt"
                else
                    log_info "All reads already have BLAST hits - no filtering needed"
                fi
            fi
        else
            log_warn "No read table found to filter"
        fi
    fi
    
    # Count reads processed
    local total_reads=$(wc -l < "$FILTERED_DIR/reads_per_sample.txt" 2>/dev/null | awk '{print $1-1}' || echo 0)
    local reads_with_tax=$(wc -l < "$TAXONOMY_DIR/reads/taxonomy_reads_hits_only.tsv" 2>/dev/null | awk '{print $1-1}' || echo 0)
    local reads_unassigned=$((total_reads - reads_with_tax))
    
    log_success "Read-level taxonomy assignment completed"
    log_info "Output files:"
    log_info "  Combined BLAST:      $TAXONOMY_DIR/reads/blast_results_reads.txt"
    log_info "  Reads per sample:    $FILTERED_DIR/reads_per_sample.txt"
    log_info "  Taxonomy (all):      $TAXONOMY_DIR/reads/taxonomy_reads.tsv"
    log_info "  Taxonomy (hits):     $TAXONOMY_DIR/reads/taxonomy_reads_hits_only.tsv"
    log_info "  Expanded taxonomy:   $TAXONOMY_DIR/reads/taxonomy_reads_expanded.tsv"
    log_info "  Microeco taxonomy:   $TAXONOMY_DIR/reads/taxonomy_reads_for_microeco.tsv"
    log_info "  Filtered per sample: $TAXONOMY_DIR/reads/reads_per_sample_with_taxonomy.txt"
    log_info "Summary: total=$total_reads  assigned=$reads_with_tax  unassigned=$reads_unassigned"
}

#################################################################
# STEP 6: CLUSTERING
#################################################################
#################################################################
# STEP 6.1: CLUSTERING WITH SWARM
#################################################################
step_cluster_swarm() {
    log_step "STEP 6: Clustering with Swarm 3.0..."
    
    # Check if clustering output already exists
    if [[ -f "$CLUSTER_DIR/otus.fasta" ]] && [[ -f "$CLUSTER_DIR/otu_per_sample.txt" ]]; then
        local otu_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); otu_count=${otu_count:-0}
        if [[ $otu_count -gt 0 ]]; then
            log_warn "Clustering outputs already exist:"
            log_warn "  - OTUs: $CLUSTER_DIR/otus.fasta ($otu_count OTUs)"
            log_warn "  - OTU table: $CLUSTER_DIR/otu_per_sample.txt"
            log_info "Skipping clustering step..."
            return 0
        fi
    fi
    
    local fasta_dir="$FILTERED_DIR/fasta"
    if [[ ! -d "$fasta_dir" ]] || [[ ! $(ls -A "$fasta_dir" 2>/dev/null) ]]; then
        log_error "No FASTA files found for clustering. Run filter step first."
        return 1
    fi
    
    # Check if swarm is available
    if ! command -v swarm &> /dev/null; then
        log_error "Swarm not found. Please install Swarm 3.0"
        log_error "Installation: conda install -c bioconda swarm"
        return 1
    fi
    
    local swarm_version=$(swarm --version 2>&1 | head -1)
    log_info "Swarm version: $swarm_version"
    
    # Combine all FASTQ files for polishing
    log_info "Combining all FASTQ files..."
    if ls "$FILTERED_DIR"/*_filtered.fastq 1> /dev/null 2>&1; then
        cat "$FILTERED_DIR"/*_filtered.fastq > "$CLUSTER_DIR/all_reads.fastq"
        local fastq_count=$(grep -c "^@" "$CLUSTER_DIR/all_reads.fastq" 2>/dev/null || true); fastq_count=${fastq_count:-0}
        log_info "Combined FASTQ reads: $fastq_count"
    else
        log_warn "No FASTQ files found in $FILTERED_DIR"
        log_warn "Polishing step will not be available without FASTQ files"
    fi
    
    log_info "Combining all FASTA files..."
    cat "$fasta_dir"/*.fasta > "$CLUSTER_DIR/all_reads.fasta"
    
    local total_seqs=$(grep -c "^>" "$CLUSTER_DIR/all_reads.fasta")
    log_info "Total sequences: $total_seqs"
    
    # Step 6.1: Dereplication
    log_info "Step 6.1: Dereplicating sequences..."
    
    vsearch \
        --derep_fulllength "$CLUSTER_DIR/all_reads.fasta" \
        --output "$CLUSTER_DIR/dereplicated.fasta" \
        --sizeout \
        --threads "$THREADS" \
        --fasta_width 0 \
        --minuniquesize 1 \
        2>&1 | tee "$LOG_DIR/vsearch_derep.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/dereplicated.fasta" ]]; then
        log_error "Dereplication failed. Check log: $LOG_DIR/vsearch_derep.log"
        return 1
    fi
    
    local derep_count=$(grep -c "^>" "$CLUSTER_DIR/dereplicated.fasta" 2>/dev/null || true); derep_count=${derep_count:-0}
    if [[ $derep_count -eq 0 ]]; then
        log_error "No sequences after dereplication"
        return 1
    fi
    log_success "Dereplication complete: $derep_count unique sequences"
    
    # Step 6.2: Chimera removal (before clustering)
    if [[ "$CHIMERA_CHECK" == "true" ]]; then
        log_info "Step 6.2: Removing chimeric sequences..."
        
        vsearch \
            --uchime_denovo "$CLUSTER_DIR/dereplicated.fasta" \
            --nonchimeras "$CLUSTER_DIR/dereplicated_nochimera.fasta" \
            --chimeras "$CLUSTER_DIR/chimeras.fasta" \
            --threads "$THREADS" \
            2>&1 | tee "$LOG_DIR/vsearch_chimera.log"
        
        if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/dereplicated_nochimera.fasta" ]]; then
            log_error "Chimera removal failed. Check log: $LOG_DIR/vsearch_chimera.log"
            return 1
        fi
        
        local chimera_count=$(grep -c "^>" "$CLUSTER_DIR/chimeras.fasta" 2>/dev/null || true); chimera_count=${chimera_count:-0}
        local clean_count=$(grep -c "^>" "$CLUSTER_DIR/dereplicated_nochimera.fasta" 2>/dev/null || true); clean_count=${clean_count:-0}
        log_success "Chimera removal complete:"
        log_info "  Chimeras detected: $chimera_count"
        log_info "  Clean sequences: $clean_count"
        
        local input_for_clustering="$CLUSTER_DIR/dereplicated_nochimera.fasta"
    else
        log_info "Step 6.2: Chimera checking disabled"
        local input_for_clustering="$CLUSTER_DIR/dereplicated.fasta"
    fi
    
    # Step 6.3: Sort by abundance
    log_info "Step 6.3: Sorting by abundance..."
    vsearch \
        --sortbysize "$input_for_clustering" \
        --output "$CLUSTER_DIR/sorted.fasta" \
        --threads "$THREADS" \
        2>&1 | tee "$LOG_DIR/vsearch_sort.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/sorted.fasta" ]]; then
        log_error "Sorting failed. Check log: $LOG_DIR/vsearch_sort.log"
        return 1
    fi
    
    local sorted_count=$(grep -c "^>" "$CLUSTER_DIR/sorted.fasta" 2>/dev/null || true); sorted_count=${sorted_count:-0}
    if [[ $sorted_count -eq 0 ]]; then
        log_error "No sequences after sorting"
        return 1
    fi
    log_success "Sorting complete: $sorted_count sequences"
    
    # Step 6.4: Swarm clustering
    log_info "Step 6.4: Clustering with Swarm..."
    log_info "Parameters: d=$SWARM_D, fastidious=$SWARM_FASTIDIOUS, boundary=$SWARM_BOUNDARY"
    
    local swarm_cmd="swarm \
        -d $SWARM_D \
        -t $THREADS \
        -z \
        -o $CLUSTER_DIR/swarm_clusters.txt \
        -s $CLUSTER_DIR/swarm_statistics.txt \
        -w $CLUSTER_DIR/otus.fasta \
        -u $CLUSTER_DIR/swarm_otus.uc"
    
    # Add fastidious mode if enabled
    if [[ "$SWARM_FASTIDIOUS" == "true" ]]; then
        swarm_cmd="$swarm_cmd -f"
        # Add boundary (only valid with fastidious mode)
        if [[ -n "$SWARM_BOUNDARY" ]] && [[ "$SWARM_BOUNDARY" -gt 0 ]]; then
            swarm_cmd="$swarm_cmd -b $SWARM_BOUNDARY"
        fi
    fi
    
    # Add no OTU breaking if specified
    if [[ "$SWARM_NO_OTU_BREAKING" == "true" ]]; then
        swarm_cmd="$swarm_cmd -n"
    fi
    
    swarm_cmd="$swarm_cmd $CLUSTER_DIR/sorted.fasta"
    
    log_info "Running: $swarm_cmd"
    eval $swarm_cmd 2>&1 | tee "$LOG_DIR/swarm_cluster.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/otus.fasta" ]]; then
        log_error "Swarm clustering failed. Check log: $LOG_DIR/swarm_cluster.log"
        return 1
    fi
    
    local cluster_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); cluster_count=${cluster_count:-0}
    if [[ $cluster_count -eq 0 ]]; then
        log_error "No clusters generated"
        return 1
    fi
    
    log_success "Swarm clustering complete: $cluster_count OTUs generated"
    
    # Display Swarm statistics
    if [[ -f "$CLUSTER_DIR/swarm_statistics.txt" ]]; then
        log_info "Swarm statistics (top 5 OTUs):"
        head -6 "$CLUSTER_DIR/swarm_statistics.txt" | tail -5 | while read line; do
            log_info "  $line"
        done
    fi

    # Remove size annotations and normalize headers to lowercase, sequences to uppercase
    log_info "Removing size annotations, normalizing OTU IDs to lowercase, sequences to uppercase..."
    # Remove ;size=N then strip any leftover trailing semicolons on header lines
    sed -i '/^>/ { s/;size=[0-9]*//g; s/;[[:space:]]*$//; }' "$CLUSTER_DIR/otus.fasta"
    
    # Convert headers to lowercase, sequences to uppercase
    awk '/^>/ {print tolower($0); next} {print toupper($0)}' "$CLUSTER_DIR/otus.fasta" > "$CLUSTER_DIR/otus_temp.fasta"
    mv "$CLUSTER_DIR/otus_temp.fasta" "$CLUSTER_DIR/otus.fasta"

    
    # Create OTU table per barcode
    log_info "Creating OTU table per barcode..."
    mkdir -p "$CLUSTER_DIR/per_barcode"
    
    # Check if per-barcode tables already exist
    local existing_tables=$(find "$CLUSTER_DIR/per_barcode" -name "*_otu_table.txt" -type f 2>/dev/null | wc -l)
    local total_barcodes=$(find "$fasta_dir" -name "*.fasta" -type f | wc -l)
    
    if [[ $existing_tables -gt 0 ]]; then
        log_warn "Found $existing_tables existing per-barcode OTU tables in $CLUSTER_DIR/per_barcode"
        log_info "Checking which barcodes need processing..."
    fi
    
    log_info "Total barcodes to process: $total_barcodes"
    
    # Process each barcode
    local success_count=0
    local skipped_count=0
    local barcode_num=0
    
    while IFS= read -r fasta_file; do
        barcode_num=$((barcode_num + 1))
        barcode=$(basename "$fasta_file" _filtered.fasta)
        otu_table="$CLUSTER_DIR/per_barcode/${barcode}_otu_table.txt"
        
        # Skip if OTU table already exists
        if [[ -f "$otu_table" ]] && [[ -s "$otu_table" ]]; then
            otus=$(tail -n +2 "$otu_table" | wc -l 2>/dev/null || true)
            if [[ $otus -gt 0 ]]; then
                log_info "[$barcode_num/$total_barcodes] $barcode - already processed ($otus OTUs) - skipping"
                skipped_count=$((skipped_count + 1))
                continue
            fi
        fi
        
        reads=$(grep -c "^>" "$fasta_file" 2>/dev/null || true); reads=${reads:-0}
        log_info "[$barcode_num/$total_barcodes] Processing $barcode ($reads reads)..."
        
        [[ $reads -eq 0 ]] && { log_warn "  Skipping (empty)"; continue; }
        
        # Relabel reads
        vsearch --fastx_filter "$fasta_file" --relabel "${barcode}." --fastaout "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta" --quiet
        
        # Map to OTUs
        vsearch --usearch_global "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta" --db "$CLUSTER_DIR/otus.fasta" --id "$CLUSTER_IDENTITY" --otutabout "$otu_table" --threads "$THREADS" --quiet
        
        # Check result
        if [[ -s "$otu_table" ]]; then
            otus=$(tail -n +2 "$otu_table" | wc -l)
            log_success "  → $otus OTUs detected"
            success_count=$((success_count + 1))
        else
            log_error "  → Failed to generate OTU table"
        fi
        
        # Cleanup
        rm -f "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta"
        
    done < <(find "$fasta_dir" -name "*.fasta" -type f)
    
    log_success "Per-barcode processing complete:"
    log_info "  - Newly processed: $success_count"
    log_info "  - Skipped (existing): $skipped_count"
    log_info "  - Total: $((success_count + skipped_count))/$total_barcodes"
    
# Merge per-barcode tables
    log_info "Merging per-barcode OTU tables..."
    python3 << PYEOF
import os
import pandas as pd
from pathlib import Path

cluster_dir = os.environ.get('CLUSTER_DIR', '$CLUSTER_DIR')
if not cluster_dir or cluster_dir == '\$CLUSTER_DIR':
    raise RuntimeError("CLUSTER_DIR environment variable is not set.")

per_barcode_dir = Path(cluster_dir) / "per_barcode"
if not per_barcode_dir.exists():
    raise RuntimeError(f"per_barcode directory not found: {per_barcode_dir}")

otu_tables = []
bad_files = []

for otu_file in sorted(per_barcode_dir.glob("*_otu_table.txt")):
    try:
        df = pd.read_csv(otu_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"ERROR reading {otu_file}: {e}")
        bad_files.append(otu_file)
        continue
    
    if df.shape[1] == 0:
        print(f"Skipping {otu_file} — no columns found")
        bad_files.append(otu_file)
        continue
    
    barcode = otu_file.stem.replace("_otu_table", "")
    df = df.copy()
    df.columns = [barcode]
    otu_tables.append(df)

print(f"Found {len(otu_tables)} OTU tables, {len(bad_files)} skipped due to errors or emptiness.")

if otu_tables:
    merged = pd.concat(otu_tables, axis=1, join='outer', sort=True).fillna(0).astype(int)
    merged.index.name = '#OTU ID'
    output_file = Path(cluster_dir) / "otu_per_sample.txt"
    merged.to_csv(output_file, sep='\t')
    print(f"Merged OTU table written: {output_file}")
    print(f"Dimensions: {merged.shape[0]} OTUs × {merged.shape[1]} samples")
else:
    print("No valid OTU tables found — nothing to merge.")
PYEOF
        
    log_success "Swarm clustering completed"
    log_info "Output files:"
    log_info "  - OTUs: $CLUSTER_DIR/otus.fasta"
    log_info "  - Per-barcode: $CLUSTER_DIR/per_barcode/*_otu_table.txt"
    log_info "  - Combined: $CLUSTER_DIR/otu_table.txt"
    log_info "  - Merged: $CLUSTER_DIR/otu_per_sample.txt"    
}

step_cluster() {
    # Route to appropriate clustering tool based on config
    if [[ "$CLUSTER_TOOL" == "swarm" ]]; then
        step_cluster_swarm
    else
        # Default to vsearch
        step_cluster_size
    fi
}

#################################################################
# STEP 6: CLUSTERING WITH VSEARCH
#################################################################
step_cluster_size() {
    log_step "STEP 6: Sequence clustering (--cluster_size)..."
    
    # Check if clustering output already exists
    if [[ -f "$CLUSTER_DIR/otus.fasta" ]] && [[ -f "$CLUSTER_DIR/otu_per_sample.txt" ]]; then
        local otu_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); otu_count=${otu_count:-0}
        if [[ $otu_count -gt 0 ]]; then
            log_warn "Clustering outputs already exist:"
            log_warn "  - OTUs: $CLUSTER_DIR/otus.fasta ($otu_count OTUs)"
            log_warn "  - OTU table: $CLUSTER_DIR/otu_per_sample.txt"
            log_info "Skipping clustering step..."
            return 0
        fi
    fi
    
    local fasta_dir="$FILTERED_DIR/fasta"
    if [[ ! -d "$fasta_dir" ]] || [[ ! $(ls -A "$fasta_dir" 2>/dev/null) ]]; then
        log_error "No FASTA files found for clustering. Run filter step first."
        return 1
    fi
    
    # Combine all FASTQ files for polishing
    log_info "Combining all FASTQ files..."
    if ls "$FILTERED_DIR"/*_filtered.fastq 1> /dev/null 2>&1; then
        cat "$FILTERED_DIR"/*_filtered.fastq > "$CLUSTER_DIR/all_reads.fastq"
        local fastq_count=$(grep -c "^@" "$CLUSTER_DIR/all_reads.fastq" 2>/dev/null || true); fastq_count=${fastq_count:-0}
        log_info "Combined FASTQ reads: $fastq_count"
    else
        log_warn "No FASTQ files found in $FILTERED_DIR"
        log_warn "Polishing step will not be available without FASTQ files"
    fi
    
    log_info "Combining all FASTA files..."
    cat "$fasta_dir"/*.fasta > "$CLUSTER_DIR/all_reads.fasta"
    
    local total_seqs=$(grep -c "^>" "$CLUSTER_DIR/all_reads.fasta")
    log_info "Total sequences: $total_seqs"
    
    # Step 6.1: Dereplication
    log_info "Step 6.1: Dereplicating sequences..."
    log_info "This may take a few minutes depending on dataset size..."
    
    vsearch \
        --derep_fulllength "$CLUSTER_DIR/all_reads.fasta" \
        --output "$CLUSTER_DIR/dereplicated.fasta" \
        --sizeout \
        --threads "$THREADS" \
        --fasta_width 0 \
        --minuniquesize 1 \
        2>&1 | tee "$LOG_DIR/vsearch_derep.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/dereplicated.fasta" ]]; then
        log_error "Dereplication failed. Check log: $LOG_DIR/vsearch_derep.log"
        return 1
    fi
    
    local derep_count=$(grep -c "^>" "$CLUSTER_DIR/dereplicated.fasta" 2>/dev/null || true); derep_count=${derep_count:-0}
    if [[ $derep_count -eq 0 ]]; then
        log_error "No sequences after dereplication"
        return 1
    fi
    log_success "Dereplication complete: $derep_count unique sequences"
    
    # Step 6.2: Chimera removal (before clustering)
    if [[ "$CHIMERA_CHECK" == "true" ]]; then
        log_info "Step 6.2: Removing chimeric sequences..."
        log_info "Checking for chimeras using de novo method..."
        
        vsearch \
            --uchime_denovo "$CLUSTER_DIR/dereplicated.fasta" \
            --nonchimeras "$CLUSTER_DIR/dereplicated_nochimera.fasta" \
            --chimeras "$CLUSTER_DIR/chimeras.fasta" \
            --threads "$THREADS" \
            2>&1 | tee "$LOG_DIR/vsearch_chimera.log"
        
        if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/dereplicated_nochimera.fasta" ]]; then
            log_error "Chimera removal failed. Check log: $LOG_DIR/vsearch_chimera.log"
            return 1
        fi
        
        local chimera_count=$(grep -c "^>" "$CLUSTER_DIR/chimeras.fasta" 2>/dev/null || true); chimera_count=${chimera_count:-0}
        local clean_count=$(grep -c "^>" "$CLUSTER_DIR/dereplicated_nochimera.fasta" 2>/dev/null || true); clean_count=${clean_count:-0}
        log_success "Chimera removal complete:"
        log_info "  Chimeras detected: $chimera_count"
        log_info "  Clean sequences: $clean_count"
        
        # Use cleaned sequences for next steps
        local input_for_sorting="$CLUSTER_DIR/dereplicated_nochimera.fasta"
    else
        log_info "Step 6.2: Chimera checking disabled"
        local input_for_sorting="$CLUSTER_DIR/dereplicated.fasta"
    fi
    
    
    # Step 6.3: Sort by abundance
    log_info "Step 6.3: Sorting by abundance..."
    vsearch \
        --sortbysize "$input_for_sorting" \
        --output "$CLUSTER_DIR/sorted.fasta" \
        --threads "$THREADS" \
        2>&1 | tee "$LOG_DIR/vsearch_sort.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/sorted.fasta" ]]; then
        log_error "Sorting failed. Check log: $LOG_DIR/vsearch_sort.log"
        return 1
    fi
    
    local sorted_count=$(grep -c "^>" "$CLUSTER_DIR/sorted.fasta" 2>/dev/null || true); sorted_count=${sorted_count:-0}
    
    if [[ $sorted_count -eq 0 ]]; then
        log_error "No sequences after sorting step"
        log_error "Check log: $LOG_DIR/vsearch_sort.log"
        return 1
    fi
    
    log_success "Sorting complete: $sorted_count sequences"
    
    
    # Step 6.4: Clustering with vsearch...
    log_info "Step 6.4: Clustering with vsearch..."
    
    vsearch \
        --cluster_size "$CLUSTER_DIR/sorted.fasta" \
        --id "$CLUSTER_IDENTITY" \
        --centroids "$CLUSTER_DIR/otus.fasta" \
        --uc "$CLUSTER_DIR/clusters.uc" \
        --threads "$THREADS" \
        2>&1 | tee "$LOG_DIR/vsearch_cluster.log"
    
    if [[ $? -ne 0 ]] || [[ ! -f "$CLUSTER_DIR/otus.fasta" ]]; then
        log_error "Clustering failed. Check log: $LOG_DIR/vsearch_cluster.log"
        return 1
    fi
    
    local cluster_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); cluster_count=${cluster_count:-0}
    
    if [[ $cluster_count -eq 0 ]]; then
        log_error "No clusters generated"
        log_error "Check log: $LOG_DIR/vsearch_cluster.log"
        return 1
    fi
    
    log_success "Clustering complete: $cluster_count OTUs generated"

    # Remove size annotations and normalize headers to lowercase, sequences to uppercase
    log_info "Removing size annotations, normalizing OTU IDs to lowercase, sequences to uppercase..."
    # Remove ;size=N then strip any leftover trailing semicolons on header lines
    sed -i '/^>/ { s/;size=[0-9]*//g; s/;[[:space:]]*$//; }' "$CLUSTER_DIR/otus.fasta"
    
    # Convert headers to lowercase, sequences to uppercase
    awk '/^>/ {print tolower($0); next} {print toupper($0)}' "$CLUSTER_DIR/otus.fasta" > "$CLUSTER_DIR/otus_temp.fasta"
    mv "$CLUSTER_DIR/otus_temp.fasta" "$CLUSTER_DIR/otus.fasta"
    
    # Create OTU table per barcode
    log_info "Creating OTU table per barcode..."
    mkdir -p "$CLUSTER_DIR/per_barcode"
    
    # Check if per-barcode tables already exist
    local existing_tables=$(find "$CLUSTER_DIR/per_barcode" -name "*_otu_table.txt" -type f 2>/dev/null | wc -l)
    local total_barcodes=$(find "$fasta_dir" -name "*.fasta" -type f | wc -l)
    
    if [[ $existing_tables -gt 0 ]]; then
        log_warn "Found $existing_tables existing per-barcode OTU tables in $CLUSTER_DIR/per_barcode"
        log_info "Checking which barcodes need processing..."
    fi
    
    log_info "Total barcodes to process: $total_barcodes"
    
    # Process each barcode
    local success_count=0
    local skipped_count=0
    local barcode_num=0
    
    while IFS= read -r fasta_file; do
        barcode_num=$((barcode_num + 1))
        barcode=$(basename "$fasta_file" _filtered.fasta)
        otu_table="$CLUSTER_DIR/per_barcode/${barcode}_otu_table.txt"
        
        # Skip if OTU table already exists
        if [[ -f "$otu_table" ]] && [[ -s "$otu_table" ]]; then
            otus=$(tail -n +2 "$otu_table" | wc -l 2>/dev/null || true)
            if [[ $otus -gt 0 ]]; then
                log_info "[$barcode_num/$total_barcodes] $barcode - already processed ($otus OTUs) - skipping"
                skipped_count=$((skipped_count + 1))
                continue
            fi
        fi
        
        reads=$(grep -c "^>" "$fasta_file" 2>/dev/null || true); reads=${reads:-0}
        log_info "[$barcode_num/$total_barcodes] Processing $barcode ($reads reads)..."
        
        [[ $reads -eq 0 ]] && { log_warn "  Skipping (empty)"; continue; }
        
        # Relabel reads
        vsearch --fastx_filter "$fasta_file" --relabel "${barcode}." --fastaout "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta" --quiet
        
        # Map to OTUs
        vsearch --usearch_global "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta" --db "$CLUSTER_DIR/otus.fasta" --id "$CLUSTER_IDENTITY" --otutabout "$otu_table" --threads "$THREADS" --quiet
        
        # Check result
        if [[ -s "$otu_table" ]]; then
            otus=$(tail -n +2 "$otu_table" | wc -l)
            log_success "  → $otus OTUs detected"
            success_count=$((success_count + 1))
        else
            log_error "  → Failed to generate OTU table"
        fi
        
        # Cleanup
        rm -f "$CLUSTER_DIR/per_barcode/${barcode}_relabeled.fasta"
        
    done < <(find "$fasta_dir" -name "*.fasta" -type f)
    
    log_success "Per-barcode processing complete:"
    log_info "  - Newly processed: $success_count"
    log_info "  - Skipped (existing): $skipped_count"
    log_info "  - Total: $((success_count + skipped_count))/$total_barcodes"
    
# Merge per-barcode tables
    log_info "Merging per-barcode OTU tables..."
    python3 << PYEOF
import os
import pandas as pd
from pathlib import Path

cluster_dir = os.environ.get('CLUSTER_DIR', '$CLUSTER_DIR')
if not cluster_dir or cluster_dir == '\$CLUSTER_DIR':
    raise RuntimeError("CLUSTER_DIR environment variable is not set.")

per_barcode_dir = Path(cluster_dir) / "per_barcode"
if not per_barcode_dir.exists():
    raise RuntimeError(f"per_barcode directory not found: {per_barcode_dir}")

otu_tables = []
bad_files = []

for otu_file in sorted(per_barcode_dir.glob("*_otu_table.txt")):
    try:
        df = pd.read_csv(otu_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"ERROR reading {otu_file}: {e}")
        bad_files.append(otu_file)
        continue
    
    if df.shape[1] == 0:
        print(f"Skipping {otu_file} — no columns found")
        bad_files.append(otu_file)
        continue
    
    barcode = otu_file.stem.replace("_otu_table", "")
    df = df.copy()
    df.columns = [barcode]
    otu_tables.append(df)

print(f"Found {len(otu_tables)} OTU tables, {len(bad_files)} skipped due to errors or emptiness.")

if otu_tables:
    merged = pd.concat(otu_tables, axis=1, join='outer', sort=True).fillna(0).astype(int)
    merged.index.name = '#OTU ID'
    output_file = Path(cluster_dir) / "otu_per_sample.txt"
    merged.to_csv(output_file, sep='\t')
    print(f"Merged OTU table written: {output_file}")
    print(f"Dimensions: {merged.shape[0]} OTUs × {merged.shape[1]} samples")
else:
    print("No valid OTU tables found — nothing to merge.")
PYEOF
        
    log_success "vsearch clustering completed (--cluster_size)"
    log_info "Output files:"
    log_info "  - OTUs: $CLUSTER_DIR/otus.fasta"
    log_info "  - Per-barcode: $CLUSTER_DIR/per_barcode/*_otu_table.txt"
    log_info "  - Combined: $CLUSTER_DIR/otu_table.txt"  
    log_info "  - Merged: $CLUSTER_DIR/otu_per_sample.txt"  
}

#################################################################
# STEP 7: POLISH WITH MEDAKA (OPTIONAL)
#################################################################
step_polish() {
    log_step "STEP 7: Polishing with Medaka (controlled batch mode)..."

    if [[ "$USE_MEDAKA" != "true" ]]; then
        log_warn "Polishing is disabled in config. Skipping..."
        return 0
    fi

    # Check if polished OTUs already exist
    if [[ -f "$CLUSTER_DIR/otus_unpolished.fasta" ]] && [[ -f "$CLUSTER_DIR/otus.fasta" ]]; then
        local polished_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); polished_count=${polished_count:-0}
        local unpolished_count=$(grep -c "^>" "$CLUSTER_DIR/otus_unpolished.fasta" 2>/dev/null || true); unpolished_count=${unpolished_count:-0}
        if [[ $polished_count -gt 0 ]] && [[ $unpolished_count -gt 0 ]]; then
            log_warn "Polished OTUs already exist:"
            log_warn "  - Polished: $CLUSTER_DIR/otus.fasta ($polished_count OTUs)"
            log_warn "  - Unpolished backup: $CLUSTER_DIR/otus_unpolished.fasta ($unpolished_count OTUs)"
            log_info "Skipping polishing step..."
            return 0
        fi
    fi

    # Check for required tools
    if ! command -v minimap2 &> /dev/null; then
        log_error "minimap2 not found. Please install minimap2."
        return 1
    fi
    if ! command -v samtools &> /dev/null; then
        log_error "samtools not found. Please install samtools."
        return 1
    fi
    if ! command -v medaka_consensus &> /dev/null; then
        log_error "medaka not found. Please install medaka."
        return 1
    fi
    if ! command -v seqkit &> /dev/null; then
        log_error "seqkit not found. Please install seqkit."
        return 1
    fi

    # Validate batch config — jobs × threads/job must not exceed THREADS
    local total_threads=$(( POLISH_JOBS * POLISH_THREADS_PER_JOB ))
    if [[ $total_threads -gt $THREADS ]]; then
        log_warn "polish_jobs ($POLISH_JOBS) × threads_per_job ($POLISH_THREADS_PER_JOB) = $total_threads exceeds allocated threads ($THREADS)"
        log_warn "Reducing polish_jobs to $(( THREADS / POLISH_THREADS_PER_JOB )) to stay within allocation"
        POLISH_JOBS=$(( THREADS / POLISH_THREADS_PER_JOB ))
        total_threads=$(( POLISH_JOBS * POLISH_THREADS_PER_JOB ))
    fi

    log_info "Batch config: $POLISH_JOBS concurrent jobs × $POLISH_THREADS_PER_JOB threads/job = $total_threads total threads"
    log_info "Medaka internal batch size: $POLISH_MEDAKA_BATCH"
    log_info "Medaka model: $MEDAKA_MODEL"

    # Combine all FASTQ files for polishing
    log_info "Combining all FASTQ files..."
    if ls "$FILTERED_DIR"/*_filtered.fastq 1> /dev/null 2>&1; then
        cat "$FILTERED_DIR"/*_filtered.fastq > "$CLUSTER_DIR/all_reads.fastq"
        local fastq_count=$(grep -c "^@" "$CLUSTER_DIR/all_reads.fastq" 2>/dev/null || true); fastq_count=${fastq_count:-0}
        log_info "Combined FASTQ reads: $fastq_count"
    else
        log_warn "No FASTQ files found in $FILTERED_DIR"
        log_warn "Polishing step will not be available without FASTQ files"
    fi

    if [[ ! -f "$CLUSTER_DIR/otus.fasta" ]]; then
        log_error "OTU file not found. Run clustering step first."
        return 1
    fi
    if [[ ! -f "$CLUSTER_DIR/all_reads.fastq" ]]; then
        log_error "all_reads.fastq not found in $CLUSTER_DIR"
        return 1
    fi

    local otu_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta")
    log_info "Total OTUs to consider: $otu_count"

    # Backup original OTUs
    cp "$CLUSTER_DIR/otus.fasta" "$CLUSTER_DIR/otus_unpolished.fasta"
    log_info "Original OTUs backed up to: $CLUSTER_DIR/otus_unpolished.fasta"

    cd "$POLISH_DIR" || return 1
    cp "$CLUSTER_DIR/all_reads.fastq" "$POLISH_DIR/all_reads.fastq"
    cp "$CLUSTER_DIR/otus.fasta"      "$POLISH_DIR/otus.fasta"

    #-----------------------------------------------------------------
    # Step 7.1: Sanitize OTU headers
    #-----------------------------------------------------------------
    log_info "Step 7.1: Sanitizing OTU headers (removing ;size=N)..."
    sed -i '/^>/ { s/;size=[0-9]*//g; s/;[[:space:]]*$//; }' "$POLISH_DIR/otus.fasta"

    #-----------------------------------------------------------------
    # Step 7.2: Map all reads to OTUs with minimap2
    #-----------------------------------------------------------------
    if [[ -f "$POLISH_DIR/reads_to_otus.bam" ]] && [[ -f "$POLISH_DIR/reads_to_otus.bam.bai" ]]; then
        log_warn "BAM file already exists — skipping minimap2 mapping"
    else
        log_info "Step 7.2: Mapping reads to OTUs with minimap2 (-t $THREADS)..."
        minimap2 -ax map-ont -t "$THREADS" "$POLISH_DIR/otus.fasta" "$POLISH_DIR/all_reads.fastq" \
            2> "$LOG_DIR/minimap2.log" \
            | samtools view -bS - \
            | samtools sort -@ "$THREADS" -o "$POLISH_DIR/reads_to_otus.bam" \
            2> "$LOG_DIR/samtools_sort.log"
    fi

    if [[ $? -ne 0 ]] || [[ ! -s "$POLISH_DIR/reads_to_otus.bam" ]]; then
        log_error "Minimap2 mapping failed. Check: $LOG_DIR/minimap2.log"
        cd - > /dev/null; return 1
    fi

    samtools index "$POLISH_DIR/reads_to_otus.bam" 2> "$LOG_DIR/samtools_index.log"
    if [[ $? -ne 0 ]]; then
        log_error "BAM indexing failed. Check: $LOG_DIR/samtools_index.log"
        cd - > /dev/null; return 1
    fi
    log_success "Mapping completed and BAM indexed"

    #-----------------------------------------------------------------
    # Step 7.3: Extract OTU names
    #-----------------------------------------------------------------
    log_info "Step 7.3: Extracting OTU names..."
    grep '>' "$POLISH_DIR/otus.fasta" | sed 's/>//' > "$POLISH_DIR/otu_list.txt"
    local otu_list_count=$(wc -l < "$POLISH_DIR/otu_list.txt")
    log_info "Found $otu_list_count OTUs"

    #-----------------------------------------------------------------
    # Step 7.4: Controlled batch polishing via FIFO semaphore
    #           POLISH_JOBS concurrent jobs × POLISH_THREADS_PER_JOB each
    #           Total threads never exceeds POLISH_JOBS × POLISH_THREADS_PER_JOB
    #-----------------------------------------------------------------
    log_info "Step 7.4: Running controlled batch polishing..."
    log_info "  $POLISH_JOBS concurrent OTUs × $POLISH_THREADS_PER_JOB threads each × batch_size=$POLISH_MEDAKA_BATCH"

    # Counter files — each background job appends one line on completion
    local _cnt_done="$POLISH_DIR/.cnt_done"
    local _cnt_skip="$POLISH_DIR/.cnt_skip"
    local _cnt_fail="$POLISH_DIR/.cnt_fail"
    > "$_cnt_done"; > "$_cnt_skip"; > "$_cnt_fail"

    # Progress report interval — roughly every 10% of total OTUs
    local _prog_interval=$(( otu_list_count / 10 ))
    [[ $_prog_interval -lt 1 ]] && _prog_interval=1

    # Create FIFO semaphore
    local _sem="$POLISH_DIR/.polish_sem_$$"
    mkfifo "$_sem"
    exec 200<>"$_sem"
    rm -f "$_sem"

    # Pre-load N tokens (one per allowed concurrent job)
    local _i
    for (( _i=0; _i<POLISH_JOBS; _i++ )); do
        printf '\n' >&200
    done

    local _job_num=0
    while IFS= read -r _otu; do
        _job_num=$(( _job_num + 1 ))

        # Acquire one token — blocks here when all POLISH_JOBS slots are busy
        read -u 200 _tok

        # Launch background job (subshell inherits fd 200 and all local vars)
        {
            local _log="$LOG_DIR/medaka_${_otu}.log"

            # Extract reference sequence for this OTU
            seqkit grep -n -p "$_otu" "$POLISH_DIR/otus.fasta" \
                > "$POLISH_DIR/${_otu}.fa" 2>"$_log"

            # Extract and convert mapped reads
            samtools view -b "$POLISH_DIR/reads_to_otus.bam" "$_otu" \
                > "$POLISH_DIR/${_otu}.bam" 2>>"$_log"
            samtools index "$POLISH_DIR/${_otu}.bam"           2>>"$_log"
            samtools fastq "$POLISH_DIR/${_otu}.bam"           \
                > "$POLISH_DIR/${_otu}.fastq"                   2>>"$_log"

            # Count mapped reads for this OTU
            local _reads=0
            if [[ -s "$POLISH_DIR/${_otu}.fastq" ]]; then
                _reads=$(grep -c "^@" "$POLISH_DIR/${_otu}.fastq" 2>/dev/null || true)
                _reads=${_reads:-0}
            fi

            if [[ $_reads -gt 0 ]]; then
                # Run Medaka
                medaka_consensus \
                    -i "$POLISH_DIR/${_otu}.fastq" \
                    -d "$POLISH_DIR/${_otu}.fa" \
                    -o "$POLISH_DIR/medaka_${_otu}" \
                    -m "$MEDAKA_MODEL" \
                    -t "$POLISH_THREADS_PER_JOB" \
                    -b "$POLISH_MEDAKA_BATCH" \
                    >> "$_log" 2>&1

                if [[ $? -eq 0 ]] && [[ -f "$POLISH_DIR/medaka_${_otu}/consensus.fasta" ]]; then
                    echo "1" >> "$_cnt_done"
                    log_success "  [$_job_num/$otu_list_count] Polished: $_otu ($_reads reads)"
                else
                    echo "1" >> "$_cnt_fail"
                    log_warn "  [$_job_num/$otu_list_count] Medaka failed: $_otu ($_reads reads) — see $LOG_DIR/medaka_${_otu}.log"
                fi
            else
                echo "1" >> "$_cnt_skip"
                log_info "  [$_job_num/$otu_list_count] Skipped: $_otu (no mapped reads)"
            fi

            # Periodic progress summary (approx every 10% of total OTUs)
            local _p=$(wc -l < "$_cnt_done" 2>/dev/null || echo 0)
            local _s=$(wc -l < "$_cnt_skip" 2>/dev/null || echo 0)
            local _f=$(wc -l < "$_cnt_fail" 2>/dev/null || echo 0)
            local _tot=$(( _p + _s + _f ))
            if (( _tot % _prog_interval == 0 )) && (( _tot > 0 )); then
                log_info "--- Progress: $_tot/$otu_list_count OTUs | polished: $_p | skipped: $_s | failed: $_f ---"
            fi

            # Clean up per-OTU intermediates immediately to save disk space
            rm -f "$POLISH_DIR/${_otu}.fa"      \
                  "$POLISH_DIR/${_otu}.bam"     \
                  "$POLISH_DIR/${_otu}.bam.bai" \
                  "$POLISH_DIR/${_otu}.fastq"

            # Release semaphore token — allows the next waiting OTU to start
            printf '\n' >&200
        } &

    done < "$POLISH_DIR/otu_list.txt"

    # Wait for all background jobs to finish
    wait
    exec 200>&-

    # Tally results from counter files
    local _final_done=$(wc -l < "$_cnt_done" 2>/dev/null || echo 0)
    local _final_skip=$(wc -l < "$_cnt_skip" 2>/dev/null || echo 0)
    local _final_fail=$(wc -l < "$_cnt_fail" 2>/dev/null || echo 0)
    rm -f "$_cnt_done" "$_cnt_skip" "$_cnt_fail"

    log_success "Batch polishing completed"
    log_info "  Polished:              $_final_done / $otu_list_count OTUs"
    log_info "  Skipped (no reads):    $_final_skip / $otu_list_count OTUs"
    log_info "  Failed:                $_final_fail / $otu_list_count OTUs"

    #-----------------------------------------------------------------
    # Step 7.6: Merge all polished consensus sequences
    #-----------------------------------------------------------------
    log_info "Step 7.6: Merging polished OTUs..."

    local polished_count=$(find "$POLISH_DIR" -name "consensus.fasta" \
        -path "*/medaka_*/consensus.fasta" 2>/dev/null | wc -l)

    if [[ $polished_count -eq 0 ]]; then
        log_error "No OTUs were successfully polished"
        cd - > /dev/null; return 1
    fi

    log_info "Successfully polished $polished_count out of $otu_list_count OTUs"

    find "$POLISH_DIR" -maxdepth 2 -type f -name "consensus.fasta" \
        -exec cat {} + > "$POLISH_DIR/polished_otus.fasta" 2>/dev/null

    if [[ ! -s "$POLISH_DIR/polished_otus.fasta" ]]; then
        log_error "Failed to merge polished OTUs"
        cd - > /dev/null; return 1
    fi

    local merged_count=$(grep -c "^>" "$POLISH_DIR/polished_otus.fasta")
    log_success "Merged $merged_count polished OTUs"
    cp "$POLISH_DIR/polished_otus.fasta" "$CLUSTER_DIR/otus.fasta"

    #-----------------------------------------------------------------
    # Step 7.7: Clean up temporary files
    #-----------------------------------------------------------------
    log_info "Step 7.7: Cleaning up temporary files..."
    find "$POLISH_DIR" -maxdepth 1 -type d -name "medaka_*" -exec rm -rf {} +
    find "$POLISH_DIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.bam" \
         -o -name "*.bam.bai" -o -name "*.fastq" -o -name "*.mmi" -o -name "*.fai" \) -delete
    rm -f "$POLISH_DIR/otu_list.txt"
    rm -f "$POLISH_DIR/reads_to_otus.bam" "$POLISH_DIR/reads_to_otus.bam.bai"
    rm -f "$POLISH_DIR/otus.fasta" "$POLISH_DIR/all_reads.fastq"
    log_success "Cleanup completed"

    cd - > /dev/null

    local final_count=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta")
    log_success "Polishing completed successfully!"
    log_info "Polishing summary:"
    log_info "  - Original OTUs:        $otu_count"
    log_info "  - Polished:              $_final_done"
    log_info "  - Skipped (no reads):    $_final_skip"
    log_info "  - Failed:                $_final_fail"
    log_info "  - Final OTUs:            $final_count"
    log_info "Output files:"
    log_info "  - Unpolished OTUs: $CLUSTER_DIR/otus_unpolished.fasta"
    log_info "  - Final polished:  $CLUSTER_DIR/otus.fasta"
}

#################################################################
# STEP 8: TAXONOMY ASSIGNMENT WITH BLAST
#################################################################
step_taxonomy() {
    log_step "STEP 8: TAXONOMY ASSIGNMENT WITH BLAST..."
    
    # Check if ALL taxonomy output files already exist
    local all_files_exist=true
    local required_files=(
        "$TAXONOMY_DIR/blast_results.txt"
        "$TAXONOMY_DIR/taxonomy.tsv"
        "$TAXONOMY_DIR/taxonomy_otus_hits_only.tsv"
        "$TAXONOMY_DIR/taxonomy_otus_expanded.tsv"
        "$TAXONOMY_DIR/taxonomy_otus_expanded_hits_only.tsv"
        "$TAXONOMY_DIR/taxonomy_otus_for_microeco.tsv"
        "$TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
    )
    
    for file in "${required_files[@]}"; do
        if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
            all_files_exist=false
            break
        fi
    done
    
    if [[ "$all_files_exist" == "true" ]]; then
        log_warn "All taxonomy output files already exist"
        log_info "Skipping BLAST and parsing substeps..."
        
        if [[ -f "$TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt" ]]; then
            log_info "Filtered OTU table already exists: $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
        else
            # Check if OTU table exists
            local otu_table=""
            if [[ -f "$CLUSTER_DIR/otu_per_sample.txt" ]]; then
                otu_table="$CLUSTER_DIR/otu_per_sample.txt"
            elif [[ -f "$CLUSTER_DIR/otu_table.txt" ]]; then
                otu_table="$CLUSTER_DIR/otu_table.txt"
            fi
            
            if [[ -n "$otu_table" ]]; then
                export CLUSTER_DIR
                export TAXONOMY_DIR
                
                python3 << PYEOF3
import pandas as pd
import os
import sys

cluster_dir = os.environ.get('CLUSTER_DIR', '')
taxonomy_dir = os.environ.get('TAXONOMY_DIR', '')

otu_table_file = ""
if os.path.exists(os.path.join(cluster_dir, "otu_per_sample.txt")):
    otu_table_file = os.path.join(cluster_dir, "otu_per_sample.txt")
elif os.path.exists(os.path.join(cluster_dir, "otu_table.txt")):
    otu_table_file = os.path.join(cluster_dir, "otu_table.txt")

taxonomy_file = os.path.join(taxonomy_dir, "taxonomy.tsv")
output_file = os.path.join(taxonomy_dir, "otu_per_sample_with_taxonomy.txt")

try:
    # Read taxonomy file
    print(f"Reading taxonomy file: {taxonomy_file}")
    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")
    
    # Get OTUs with BLAST hits (not unassigned)
    otus_with_hits = taxonomy_df[taxonomy_df["Taxonomy"] != "Unassigned"]["OTU_ID"].str.lower().str.rstrip(';').tolist()
    print(f"Found {len(otus_with_hits)} OTUs with BLAST hits")
    if otus_with_hits:
        print(f"  Sample taxonomy OTU IDs: {otus_with_hits[:3]}")
    
    # Read OTU table
    print(f"Reading OTU table: {otu_table_file}")
    otu_table = pd.read_csv(otu_table_file, sep="\t", index_col=0)
    
    # Normalize index to lowercase for comparison
    otu_table.index = otu_table.index.str.lower()
    
    total_otus = len(otu_table)
    print(f"Total OTUs in table: {total_otus}")
    if total_otus > 0:
        print(f"  Sample OTU table IDs: {list(otu_table.index[:3])}")
    
    # Check if all OTUs already have hits
    otus_in_table = set(otu_table.index)
    otus_without_hits = otus_in_table - set(otus_with_hits)
    
    if len(otus_without_hits) == 0:
        print("All OTUs in the table already have BLAST hits. No filtering needed.")
        # Copy original table so the file exists
        otu_table.to_csv(output_file, sep="\t")
        sys.exit(0)
    
    print(f"OTUs without BLAST hits: {len(otus_without_hits)}")
    
    # Filter table to keep only OTUs with hits
    filtered_table = otu_table.loc[otu_table.index.isin(otus_with_hits)]
    
    if len(filtered_table) == 0 and len(otus_with_hits) > 0:
        print(f"WARNING: Zero OTU IDs matched between taxonomy and OTU table!", file=sys.stderr)
        print(f"  Taxonomy OTU IDs (first 5): {otus_with_hits[:5]}", file=sys.stderr)
        print(f"  OTU table IDs (first 5): {list(otu_table.index[:5])}", file=sys.stderr)
    
    # Save filtered table
    filtered_table.to_csv(output_file, sep="\t")
    
    print(f"\nFiltered OTU table saved: {output_file}")
    print(f"  Original OTUs: {total_otus}")
    print(f"  With taxonomy: {len(filtered_table)}")
    print(f"  Removed (no hits): {total_otus - len(filtered_table)}")
    
except Exception as e:
    print(f"Error creating filtered OTU table: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF3
                
                if [[ $? -eq 0 ]]; then
                    if [[ -f "$TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt" ]]; then
                        log_success "Filtered OTU table created: $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
                    else
                        log_info "All OTUs already have BLAST hits - no filtering needed"
                    fi
                fi
            else
                log_warn "No OTU table found to filter"
            fi
        fi
        
        log_success "OTU-level taxonomy assignment completed (using existing files)"
        log_info ""
        log_info "Output files (existing):"
        log_info "  1. BLAST results: $TAXONOMY_DIR/blast_results.txt"
        log_info "  2. OTU per sample (all): $CLUSTER_DIR/otu_per_sample.txt"
        log_info "  3. Taxonomy (all OTUs): $TAXONOMY_DIR/taxonomy.tsv"
        log_info "  4. Taxonomy (hits only): $TAXONOMY_DIR/taxonomy_otus_hits_only.tsv"
        log_info "  5. Expanded taxonomy (all): $TAXONOMY_DIR/taxonomy_otus_expanded.tsv"
        log_info "  6. Expanded taxonomy (hits only): $TAXONOMY_DIR/taxonomy_otus_expanded_hits_only.tsv"
        log_info "  7. Microeco-compatible: $TAXONOMY_DIR/taxonomy_otus_for_microeco.tsv"
        log_info "  8. Filtered OTU table: $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
        return 0
    fi
    
    if [[ ! -f "$CLUSTER_DIR/otus.fasta" ]]; then
        log_error "OTU file not found. Run clustering step first."
        return 1
    fi
    
    # Check if BLAST database exists (supports both BLAST+ v4 .nhr and v5 .ndb)
    if [[ ! -f "${BLAST_DB}.nhr" ]] && [[ ! -f "${BLAST_DB}.ndb" ]]; then
        log_error "BLAST database not found: $BLAST_DB"
        log_error "Please create the database using: makeblastdb -in your_database.fasta -dbtype nucl -out $BLAST_DB"
        return 1
    fi
    
    # Run BLAST with taxonomy in sseqid - with skip logic
    if [[ -f "$TAXONOMY_DIR/blast_results.txt" ]] && [[ -s "$TAXONOMY_DIR/blast_results.txt" ]]; then
        local blast_hits=$(wc -l < "$TAXONOMY_DIR/blast_results.txt" || echo 0)
        if [[ $blast_hits -gt 0 ]]; then
            log_warn "BLAST results already exist: $TAXONOMY_DIR/blast_results.txt ($blast_hits hits)"
            log_info "Skipping BLAST step..."
        else
            # File exists but empty, need to run BLAST
            log_info "BLAST results file exists but empty, running BLAST..."
        fi
    fi
    
    # Run BLAST if results don't exist or are empty
    if [[ ! -f "$TAXONOMY_DIR/blast_results.txt" ]] || [[ ! -s "$TAXONOMY_DIR/blast_results.txt" ]] || [[ $(wc -l < "$TAXONOMY_DIR/blast_results.txt" 2>/dev/null || true) -eq 0 ]]; then
        log_info "Running BLAST against database with embedded taxonomy..."
        log_info "Database: $BLAST_DB"
        log_info "E-value: $BLAST_EVALUE"
        log_info "Max targets: $BLAST_MAX_TARGETS"
        log_info "Using $THREADS threads"
        
        blastn \
            -query "$CLUSTER_DIR/otus.fasta" \
            -db "$BLAST_DB" \
            -task megablast \
            -reward 2 -penalty -3 -word_size 65 -gapopen 0 -gapextend 4 \
            -out "$TAXONOMY_DIR/blast_results.txt" \
            -evalue "$BLAST_EVALUE" \
            -max_target_seqs "$BLAST_MAX_TARGETS" \
            -num_threads "$THREADS" \
            -max_hsps 1 \
            -outfmt "6 qseqid sseqid pident length evalue bitscore" \
           > "$LOG_DIR/blast.log" 2>&1
    fi
    
    log_info "Step 8.1: Parsing BLAST results with embedded taxonomy..."
    
    # Export variables for Python
    export TAXONOMY_DIR
    export CLUSTER_DIR
    export BLAST_MIN_IDENTITY
    
    python3 << PYEOF
import sys
import pandas as pd
import os
import re

blast_file = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "blast_results.txt")
output_file = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "taxonomy.tsv")
cluster_dir = os.environ.get("CLUSTER_DIR", "")
otus_file = os.path.join(cluster_dir, "otus.fasta")

# Column names for TSV format (format 6)
blast_columns = ["OTU_ID", "Hit_with_Taxonomy", "Identity_%", "Alignment_Length", "E-value", "Bitscore"]

try:
    min_identity = float(os.environ.get("BLAST_MIN_IDENTITY", 80))
    # Read ALL OTU IDs from the FASTA file
    print(f"Reading all OTU IDs from: {otus_file}")
    all_otus = []
    if os.path.exists(otus_file):
        with open(otus_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    otu_id = line.strip()[1:].split()[0].lower().rstrip(';')
                    all_otus.append(otu_id)
        print(f"Total OTUs in file: {len(all_otus)}")
    else:
        print(f"Warning: OTUs file not found: {otus_file}", file=sys.stderr)
        all_otus = []
    
    # Check if BLAST results exist
    if not os.path.exists(blast_file):
        print(f"Error: BLAST results file not found: {blast_file}", file=sys.stderr)
        sys.exit(1)
    
    # Load BLAST results (TSV format - tab-separated)
    print(f"Loading BLAST results from: {blast_file}")
    blast_df = pd.read_csv(blast_file, sep="\t", names=blast_columns, header=None)
    
    # Normalize OTU IDs to lowercase
    if not blast_df.empty:
        blast_df["OTU_ID"] = blast_df["OTU_ID"].str.lower().str.rstrip(';')
        blast_df["Identity_%"] = pd.to_numeric(blast_df["Identity_%"], errors="coerce")
        blast_df["Alignment_Length"] = pd.to_numeric(blast_df["Alignment_Length"], errors="coerce")
        blast_df["E-value"] = pd.to_numeric(blast_df["E-value"], errors="coerce")
        blast_df["Bitscore"] = pd.to_numeric(blast_df["Bitscore"], errors="coerce")
        
        # Drop rows with invalid numeric values (matches blast_reads behaviour)
        blast_df = blast_df.dropna(subset=["Identity_%", "Alignment_Length"])
        
        # Apply minimum identity threshold
        before = len(blast_df)
        blast_df = blast_df[blast_df["Identity_%"] >= min_identity]
        print(f"Hits after min_identity filter ({min_identity}%): {len(blast_df)} (removed {before - len(blast_df)})")
    
    if blast_df.empty:
        print("Warning: No BLAST hits found", file=sys.stderr)
        if all_otus:
            taxonomy = pd.DataFrame({
                "OTU_ID": all_otus,
                "Best_Hit": "No hit",
                "Identity_%": 0.0,
                "Alignment_Length": 0,
                "E-value": 999.0,
                "Bitscore": 0.0,
                "Taxonomy": "Unassigned"
            })
            taxonomy.to_csv(output_file, sep="\t", index=False)
        sys.exit(0)
    
    print(f"Found {len(blast_df)} total BLAST hits")
    otu_count = blast_df["OTU_ID"].nunique()
    print(f"OTUs with BLAST hits: {otu_count}")
    
    # Extract taxonomy from sseqid column
    # Expected format: accession|taxonomy (e.g., ACC123|Bacteria;Proteobacteria;...)
    def extract_taxonomy(hit_string):
        """Extract taxonomy from BLAST subject sequence ID"""
        if pd.isna(hit_string) or hit_string == "":
            return "Unknown", "Unassigned"
        
        s = str(hit_string).strip()
        # Split by pipe to separate accession and taxonomy
        parts = s.split("|")
        
        if len(parts) >= 2:
            accession = parts[0].strip()
            taxonomy = "|".join(parts[1:]).strip()
            return accession, taxonomy
        # Fallback: no pipe but contains semicolons → treat whole string as taxonomy
        elif ";" in s:
            return "Unknown", s
        else:
            return s, "Unassigned"
    
    blast_df[["Best_Hit", "Taxonomy"]] = blast_df["Hit_with_Taxonomy"].apply(
        lambda x: pd.Series(extract_taxonomy(x))
    )
    
    # Keep only first hit per OTU (as it appears in BLAST output - no sorting/filtering)
    blast_df = blast_df.drop_duplicates(subset=["OTU_ID"], keep="first")
    
    print(f"First hits kept for {len(blast_df)} OTUs")
    
    # Create final taxonomy table
    taxonomy_df = blast_df[["OTU_ID", "Best_Hit", "Identity_%", "Alignment_Length", "E-value", "Bitscore", "Taxonomy"]].copy()
    
    # Add OTUs without hits
    otus_with_hits = set(taxonomy_df["OTU_ID"])
    otus_without_hits = [otu for otu in all_otus if otu not in otus_with_hits]
    
    if otus_without_hits:
        no_hits_df = pd.DataFrame({
            "OTU_ID": otus_without_hits,
            "Best_Hit": "No hit",
            "Identity_%": 0.0,
            "Alignment_Length": 0,
            "E-value": 999.0,
            "Bitscore": 0.0,
            "Taxonomy": "Unassigned"
        })
        taxonomy_df = pd.concat([taxonomy_df, no_hits_df], ignore_index=True)
        print(f"Added {len(otus_without_hits)} OTUs without BLAST hits")
    
    # Save complete taxonomy table (ALL OTUs - with hits + unassigned)
    taxonomy_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nTaxonomy table saved: {output_file}")
    print(f"Total OTUs: {len(taxonomy_df)}")
    with_tax = (taxonomy_df["Taxonomy"] != "Unassigned").sum()
    print(f"  With taxonomy: {with_tax}")
    unassigned = (taxonomy_df["Taxonomy"] == "Unassigned").sum()
    print(f"  Unassigned: {unassigned}")
    
    # Create taxonomy file with ONLY OTUs that have BLAST hits (excluding Unassigned)
    output_file_hits_only = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "taxonomy_otus_hits_only.tsv")
    taxonomy_hits_only = taxonomy_df[taxonomy_df["Taxonomy"] != "Unassigned"].copy()
    taxonomy_hits_only.to_csv(output_file_hits_only, sep="\t", index=False)
    print(f"\nTaxonomy hits-only table saved: {output_file_hits_only}")
    print(f"  OTUs with taxonomy hits: {len(taxonomy_hits_only)}")
    
    # Create expanded taxonomy file with separate rank columns (ALL OTUs)
    output_file_expanded = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "taxonomy_otus_expanded.tsv")
    
    def expand_taxonomy(tax_string):
        """Expand taxonomy string into separate rank columns"""
        # Expected format: Kingdom;Phylum;Class;Order;Family;Genus;Species;Strain
        ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
        
        if pd.isna(tax_string) or tax_string == "Unassigned":
            return pd.Series(["Unassigned"] * 8, index=ranks)
        
        # Split by semicolon
        parts = str(tax_string).split(";")
        
        # Pad with "Unclassified" if not enough ranks
        while len(parts) < 8:
            parts.append("Unclassified")
        
        # Take first 8 ranks
        parts = parts[:8]
        
        return pd.Series(parts, index=ranks)
    
    # Apply expansion to all OTUs (including Unassigned)
    taxonomy_ranks = taxonomy_df["Taxonomy"].apply(expand_taxonomy)
    
    # Create final dataframe with OTU_ID and taxonomic ranks (ALL OTUs)
    final_expanded = pd.DataFrame({
        "OTU_ID": taxonomy_df["OTU_ID"]
    })
    
    # Add each rank as a column
    for rank in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]:
        final_expanded[rank] = taxonomy_ranks[rank].values
    
    # Save expanded taxonomy (ALL OTUs)
    final_expanded.to_csv(output_file_expanded, sep="\t", index=False)
    print(f"\nExpanded taxonomy table saved (all OTUs): {output_file_expanded}")
    print(f"  Columns: OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, Strain")
    print(f"  Total OTUs: {len(final_expanded)}")
    
    # Create expanded taxonomy file with ONLY OTUs that have hits
    output_file_expanded_hits = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "taxonomy_otus_expanded_hits_only.tsv")
    final_expanded_hits = final_expanded[taxonomy_df["Taxonomy"] != "Unassigned"].copy()
    final_expanded_hits.to_csv(output_file_expanded_hits, sep="\t", index=False)
    print(f"\nExpanded taxonomy table saved (hits only): {output_file_expanded_hits}")
    print(f"  OTUs with taxonomy hits: {len(final_expanded_hits)}")
    
    # Create microeco-compatible file (only OTU_ID + 8 ranks, no other columns)
    output_file_microeco = os.path.join(os.environ.get("TAXONOMY_DIR", ""), "taxonomy_otus_for_microeco.tsv")
    final_expanded_hits.to_csv(output_file_microeco, sep="\t", index=False)
    print(f"\nMicroeco-compatible taxonomy table saved: {output_file_microeco}")
    print(f"  Format: OTU_ID + 8 taxonomic rank columns only")
    print(f"  Total OTUs: {len(final_expanded_hits)}")
    
    
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF
    
    if [[ $? -ne 0 ]]; then
        log_error "Failed to parse BLAST results"
        return 1
    fi
    
    # Count OTUs processed
    local total_otus=$(grep -c "^>" "$CLUSTER_DIR/otus.fasta" 2>/dev/null || true); total_otus=${total_otus:-0}
    local otus_with_tax=$(awk -F"\t" 'NR>1 && $7!="Unassigned" {count++} END {print count+0}' "$TAXONOMY_DIR/taxonomy.tsv" 2>/dev/null || true)
    local otus_unassigned=$((total_otus - otus_with_tax))
    
    log_success "OTU-level taxonomy assignment completed"
    log_info ""
    log_info "Output files created:"
    log_info "  1. BLAST results:"
    log_info "     $TAXONOMY_DIR/blast_results.txt"
    log_info "  2. OTU per sample table (all OTUs):"
    log_info "     $CLUSTER_DIR/otu_per_sample.txt"
    log_info "  3. Taxonomy file (all OTUs with/without hits):"
    log_info "     $TAXONOMY_DIR/taxonomy.tsv"
    log_info "  4. Taxonomy file (only OTUs with hits):"
    log_info "     $TAXONOMY_DIR/taxonomy_otus_hits_only.tsv"
    log_info "  5. Expanded taxonomy (all OTUs, separate rank columns):"
    log_info "     $TAXONOMY_DIR/taxonomy_otus_expanded.tsv"
    log_info "  6. Expanded taxonomy (only OTUs with hits):"
    log_info "     $TAXONOMY_DIR/taxonomy_otus_expanded_hits_only.tsv"
    log_info "  7. Microeco-compatible taxonomy (OTU_ID + ranks only):"
    log_info "     $TAXONOMY_DIR/taxonomy_otus_for_microeco.tsv"
    log_info "  8. Filtered OTU table (only OTUs with taxonomy):"
    log_info "     $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
    log_info ""
    log_info "Summary:"
    log_info "  Total OTUs processed: $total_otus"
    log_info "  OTUs with taxonomy: $otus_with_tax"
    log_info "  OTUs unassigned: $otus_unassigned"
    
    # Create filtered OTU table with only OTUs that have BLAST hits
    log_info "Step 8.2: Creating filtered OTU table (only OTUs with BLAST hits)..."
    
    # Check if filtered table already exists
    if [[ -f "$TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt" ]]; then
        log_warn "Filtered OTU table already exists: $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
        log_info "Skipping filtered OTU table creation..."
    else
        # Check if OTU table exists
        local otu_table=""
        if [[ -f "$CLUSTER_DIR/otu_per_sample.txt" ]]; then
            otu_table="$CLUSTER_DIR/otu_per_sample.txt"
        elif [[ -f "$CLUSTER_DIR/otu_table.txt" ]]; then
            otu_table="$CLUSTER_DIR/otu_table.txt"
        else
            log_warn "No OTU table found to filter"
        fi
        
        if [[ -n "$otu_table" ]]; then
            export CLUSTER_DIR
            export TAXONOMY_DIR
            
            python3 << PYEOF2
import pandas as pd
import os
import sys

cluster_dir = os.environ.get('CLUSTER_DIR', '')
taxonomy_dir = os.environ.get('TAXONOMY_DIR', '')

otu_table_file = ""
if os.path.exists(os.path.join(cluster_dir, "otu_per_sample.txt")):
    otu_table_file = os.path.join(cluster_dir, "otu_per_sample.txt")
elif os.path.exists(os.path.join(cluster_dir, "otu_table.txt")):
    otu_table_file = os.path.join(cluster_dir, "otu_table.txt")

taxonomy_file = os.path.join(taxonomy_dir, "taxonomy.tsv")
output_file = os.path.join(taxonomy_dir, "otu_per_sample_with_taxonomy.txt")

try:
    # Read taxonomy file
    print(f"Reading taxonomy file: {taxonomy_file}")
    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")
    
    # Get OTUs with BLAST hits (not unassigned)
    otus_with_hits = taxonomy_df[taxonomy_df["Taxonomy"] != "Unassigned"]["OTU_ID"].str.lower().str.rstrip(';').tolist()
    print(f"Found {len(otus_with_hits)} OTUs with BLAST hits")
    if otus_with_hits:
        print(f"  Sample taxonomy OTU IDs: {otus_with_hits[:3]}")
    
    # Read OTU table
    print(f"Reading OTU table: {otu_table_file}")
    otu_table = pd.read_csv(otu_table_file, sep="\t", index_col=0)
    
    # Normalize index to lowercase for comparison
    otu_table.index = otu_table.index.str.lower()
    
    total_otus = len(otu_table)
    print(f"Total OTUs in table: {total_otus}")
    if total_otus > 0:
        print(f"  Sample OTU table IDs: {list(otu_table.index[:3])}")
    
    # Check if all OTUs already have hits
    otus_in_table = set(otu_table.index)
    otus_without_hits = otus_in_table - set(otus_with_hits)
    
    if len(otus_without_hits) == 0:
        print("All OTUs in the table already have BLAST hits. No filtering needed.")
        # Copy original table so the file exists
        otu_table.to_csv(output_file, sep="\t")
        sys.exit(0)
    
    print(f"OTUs without BLAST hits: {len(otus_without_hits)}")
    
    # Filter table to keep only OTUs with hits
    filtered_table = otu_table.loc[otu_table.index.isin(otus_with_hits)]
    
    if len(filtered_table) == 0 and len(otus_with_hits) > 0:
        print(f"WARNING: Zero OTU IDs matched between taxonomy and OTU table!", file=sys.stderr)
        print(f"  Taxonomy OTU IDs (first 5): {otus_with_hits[:5]}", file=sys.stderr)
        print(f"  OTU table IDs (first 5): {list(otu_table.index[:5])}", file=sys.stderr)
    
    # Save filtered table
    filtered_table.to_csv(output_file, sep="\t")
    
    print(f"\nFiltered OTU table saved: {output_file}")
    print(f"  Original OTUs: {total_otus}")
    print(f"  With taxonomy: {len(filtered_table)}")
    print(f"  Removed (no hits): {total_otus - len(filtered_table)}")
    
except Exception as e:
    print(f"Error creating filtered OTU table: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
PYEOF2
            
            if [[ $? -eq 0 ]]; then
                if [[ -f "$TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt" ]]; then
                    log_success "Filtered OTU table created: $TAXONOMY_DIR/otu_per_sample_with_taxonomy.txt"
                else
                    log_info "All OTUs already have BLAST hits - no filtering needed"
                fi
            else
                log_warn "Failed to create filtered OTU table"
            fi
        fi
    fi
}



#################################################################
# MAIN PIPELINE EXECUTION
#################################################################

# Function to run a single step
run_step() {
    local step="$1"
    case $step in
        basecall)
            step_basecall
            ;;
        demultiplex)
            step_demultiplex
            ;;
        qc_before_filter)
            step_qc_before_filter
            ;;
        filter)
            step_filter
            ;;
        qc_after_filter)
            step_qc
            ;;
        qc)
            step_qc
            ;;
        blast_reads)
            step_blast_reads
            ;;
        cluster)
            step_cluster
            ;;
        polish)
            step_polish
            ;;
        taxonomy)
            step_taxonomy
            ;;
        metadata)
            generate_metadata
            ;;
        *)
            log_error "Unknown step: $step"
            return 1
            ;;
    esac
}

# Handle multi-step execution
if [[ ${#STEPS[@]} -eq 1 ]] && [[ "${STEPS[0]}" == "all" ]]; then
    # Special case: "all" runs the full pipeline
    log_info "Running complete pipeline..."
    
    # Run all steps in order
    if [[ -n "$POD5_DIR" ]] && [[ -d "$POD5_DIR" ]]; then
        step_basecall || log_warn "Basecalling failed or skipped"
    fi
    
    step_demultiplex || log_warn "Demultiplexing failed or skipped"
    step_qc_before_filter || log_warn "Pre-filter QC failed or skipped"
    step_filter || { log_error "Filtering failed"; exit 1; }
    step_qc || log_warn "Post-filter QC failed or skipped"
    
    # Use vsearch for clustering
    step_cluster || { log_error "Sequence clustering failed"; exit 1; }
    
    if [[ "$USE_MEDAKA" == "true" ]]; then
        step_polish || log_warn "Polishing failed or skipped"
    fi
    
    step_taxonomy || log_warn "Taxonomy assignment failed or skipped"
    generate_metadata || log_warn "Metadata generation failed"
    
    log_success "Pipeline completed successfully!"
else
    # Multi-step or single-step execution
    if [[ ${#STEPS[@]} -gt 1 ]]; then
        log_info "Running steps: ${STEPS[*]}"
    else
        log_info "Running step: ${STEPS[0]}"
    fi
    
    while IFS= read -r step; do
        log_info "==================== STEP: $step ===================="
        run_step "$step" || {
            log_error "Step '$step' failed"
            exit 1
        }
    done < <(printf '%s\n' "${STEPS[@]}")
    
    log_success "All steps completed successfully!"
fi

exit 0
