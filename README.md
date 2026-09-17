# rRNA-ops Pipeline

**Full-Length Bacterial rRNA Operon Analysis from Oxford Nanopore Sequencing Data**

[![Version](https://img.shields.io/badge/version-1.0-blue.svg)](https://github.com/your-repo/rRNA-ops)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)](https://github.com/your-repo/rRNA-ops)

---

## Overview

The **rRNA-ops Pipeline** is a comprehensive, automated workflow for analyzing full-length bacterial ribosomal RNA (rRNA) operon sequences from Oxford Nanopore sequencing data. The pipeline handles everything from basecalling to taxonomic annotation, producing publication-ready results.

**Key Features:**
- ✅ End-to-end automated workflow
- ✅ Modular design - run any step independently
- ✅ Multiple clustering algorithms (VSEARCH, Swarm)
- ✅ Optional GPU-accelerated basecalling
- ✅ Neural network consensus polishing
- ✅ Comprehensive quality control
- ✅ Multiple output formats for downstream analysis
- ✅ Optimized for HPC environments

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Pipeline Steps](#pipeline-steps)
5. [Output Files](#output-files)
6. [Configuration](#configuration)
7. [Examples](#examples)
8. [Citation](#citation)
9. [Support](#support)

---

## Quick Start

### Prerequisites
- Linux or macOS
- Python 3.7+
- 16+ CPU cores recommended
- 32+ GB RAM recommended
- GPU (optional, for basecalling)

### Installation

```bash
# 1. Clone repository
git clone https://github.com/your-repo/rRNA-ops.git
cd rRNA-ops

# 2. Install dependencies (see INSTALLATION_GUIDE.md)
conda env create -f rrna_conda_env.yml
conda activate rRNA-ops

# 3. Download rOpDB database
# download the database from https://www.njmicrobe.org/ropdb
# make a blast database with "makeblastdb -in sequences.fasta -dbtype nucl -out rOpDB -parse_seqids"
# give the location of the created DB in the blast_db option as blast_db: "location_of_DB/rOpDB.fasta" 
wget [rOpDB_URL]/rOpDB.fasta
makeblastdb -in rOpDB.fasta -dbtype nucl -parse_seqids

# 4. Install Dorado separately (from ONT)
# See INSTALLATION_GUIDE.md for details
```

### Basic Usage

```bash
# 1. Configure pipeline
cp config.yaml my_config.yaml
nano my_config.yaml  # Edit paths and parameters

# 2. Run the complete pipeline
./rRNA-ops_pipeline.sh all -c my_config.yaml -t 16

# 3. Check results
ls -lh results/taxonomy/
```

## Usage

### Command Structure

```bash
./rRNA-ops_pipeline.sh [STEP] [OPTIONS]
```

### Available Steps

| Step | Description | Required |
|------|-------------|----------|
| `all` | Run complete pipeline | - |
| `basecall` | Dorado basecalling | Optional* |
| `demultiplex` | Porechop demultiplexing | Optional* |
| `filter` | Chopper read filtering | **Yes** |
| `qc` | NanoPlot quality control | Recommended |
| `blast_reads` | Pre-clustering BLAST | Optional |
| `cluster` | Sequence clustering | **Yes** |
| `polish` | Medaka polishing | Optional |
| `taxonomy` | BLAST taxonomy assignment | **Yes** |
| `metadata` | Generate sample metadata | Recommended |

*Skip if starting with basecalled/demultiplexed FASTQ files

### Options

```
-c, --config FILE    Path to configuration YAML file (required)
-t, --threads INT    Number of threads (default: 4)
-v, --verbose        Enable verbose output
-h, --help           Show help message
--version            Show version
```

### Run Multiple Steps

Use a comma-separated list (no spaces):

```bash
./rRNA-ops_pipeline.sh filter,qc,cluster,taxonomy -c config.yaml -t 16
```

---

## Pipeline Steps

### 1. Basecalling (Optional)
- **Tool:** Dorado
- **Input:** POD5/FAST5 files
- **Output:** FASTQ files
- **GPU Support:** Yes

### 2. Demultiplexing (Optional)
- **Tool:** Porechop
- **Input:** FASTQ files
- **Output:** Separated FASTQ per barcode

### 3. Read Filtering
- **Tool:** Chopper
- **Input:** FASTQ files
- **Output:** Filtered reads (3500-6000 bp, Q≥10)

### 4. Quality Control
- **Tool:** NanoPlot
- **Input:** Filtered reads
- **Output:** HTML reports, statistics, plots

### 5. Pre-clustering BLAST (Optional)
- **Tool:** BLAST+
- **Input:** Filtered reads
- **Output:** Per-read taxonomy

### 6. Sequence Clustering
- **Tools:** VSEARCH or Swarm
- **Input:** Filtered reads
- **Output:** OTU sequences, abundance tables
- **Features:**
  - Chimera detection & removal
  - Denoising (VSEARCH) or high-resolution clustering (Swarm)
  - OTU table generation

### 7. Consensus Polishing (Optional)
- **Tool:** Medaka
- **Input:** OTU sequences + reads
- **Output:** Polished consensus sequences

### 8. Taxonomic Assignment
- **Tool:** BLAST+
- **Input:** OTU sequences (or polished)
- **Database:** rOpDB (bacterial rRNA operons)
- **Output:** Multiple taxonomy tables

### 9. Metadata Generation
- **Output:** Sample statistics table

For detailed information, see [PIPELINE_STEPS.md](PIPELINE_STEPS.md).

---

## Output Files

### Directory Structure

```
results/
├── filtered/                         # Filtered FASTQ files
├── qc/                              # Quality control reports
├── cluster/
│   ├── otus.fasta                   # OTU sequences
│   └── otu_per_sample.txt           # Abundance table
├── taxonomy/
│   ├── taxonomy.tsv                 # Main taxonomy file
│   ├── taxonomy_otus_expanded.tsv   # Separate rank columns
│   ├── taxonomy_otus_for_microeco.tsv  # R microeco format
│   └── otu_per_sample_with_taxonomy.txt  # Filtered abundance
└── metadata.tsv                     # Sample metadata
```

## Configuration

### Minimal config.yaml

```yaml
input:
  fastq_dir: "/path/to/fastq"

output:
  output_dir: "/path/to/results"

filtering:
  min_length: 3500
  max_length: 6000
  min_quality: 10

clustering:
  tool: "swarm"
  swarm_d: 1
  swarm_fastidious: true
  remove_chimeras: true

taxonomy:
  blast_db: "/path/to/rOpDB.fasta"  # included with the pipeline in rOpDB folder, but user can supply their own database by making one with makeblastdb -in sequences.fasta -dbtype nucl -out my_nucl_db -parse_seqids command and giving the location of the created DB in  blast_db: option
  evalue: 20
  min_identity: 50
```

## Examples

### Example 1: Complete Pipeline

```bash
# Starting from POD5 files
./rRNA-ops_pipeline.sh all -c config.yaml -t 32

# Results in: results/taxonomy/taxonomy.tsv
```

### Example 2: Skip Basecalling & Demultiplexing

```bash
# Starting from filtered FASTQ files
./rRNA-ops_pipeline.sh filter,qc,cluster,taxonomy -c config.yaml -t 16
```

### Example 3: Re-run Clustering with Different Parameters

```bash
# Edit config.yaml clustering parameters
nano config.yaml

# Re-run clustering and taxonomy
./rRNA-ops_pipeline.sh cluster,taxonomy -c config.yaml -t 16
```

### Optimization Tips

1. **Use GPU for basecalling** - 10-50x speedup
2. **Skip polishing** - Saves 30-50% time with minimal accuracy loss
3. **Adjust thread count** - More threads for BLAST/clustering

---

## Troubleshooting

### Common Issues

**Problem:** Low OTU assignment rate
```bash
# Check BLAST parameters
# Decrease min_identity or increase evalue
```

**Problem:** Too many OTUs
```bash
# Increase clustering stringency
# VSEARCH: increase unoise_alpha
# Swarm: increase swarm_d
```

**Problem:** High chimera rate
```bash
# Check if remove_chimeras: true
# Review PCR conditions
# Check filtering stringency
```

**Problem:** Memory errors
```bash
# Reduce thread count
# Process samples in batches
# Use more swap space
```

## Citation

If you use this pipeline in your research, please cite:

```
[Your Citation Here]
rRNA-ops Pipeline v1.0
https://github.com/your-repo/rRNA-ops
```

**Also cite the tools used:**

- **Dorado:** Oxford Nanopore Technologies
- **Porechop:** Wick et al. (2017)
- **VSEARCH:** Rognes et al. (2016)
- **Swarm:** Mahé et al. (2014)
- **Medaka:** Oxford Nanopore Technologies
- **BLAST:** Camacho et al. (2009)

---

## Documentation

- [Installation Guide](INSTALLATION_GUIDE.md) - Detailed setup instructions
- [Pipeline Steps](PIPELINE_STEPS.md) - Detailed pipeline steps
---

## Acknowledgments

- Oxford Nanopore Technologies for sequencing technology
- rOpDB database maintainers
- Open-source bioinformatics community

---

## Version History

**v1.0** (2026-02-13)
- Initial release
- VSEARCH and Swarm clustering
- Full taxonomy assignment pipeline
- Multiple output formats

---