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
git https://github.com/anilchauhanhp9/rRNA-ops-pipeline.git
cd rRNA-ops-pipeline-main

# 2. Install dependencies (see INSTALLATION_GUIDE.md)
conda env create -f rrna_conda_env.yml
conda activate rRNAops-pipeline

# 3. Download rOpDB database
# download the database from https://www.njmicrobe.org/ropdb
# make a blast database with "makeblastdb -in sequences.fasta -dbtype nucl -out rOpDB -parse_seqids"
# give the location of the created DB in the blast_db option as blast_db: "location_of_DB/rOpDB.fasta" 

# 4. Install Dorado separately (from ONT)
# See INSTALLATION_GUIDE.md for details
```

### Basic Usage

```bash
# 1. Configure pipeline
nano config.yaml  # Edit paths and parameters

# 2. Run the complete pipeline
./rRNA-ops_pipeline.sh all -c config.yaml -t 16

# 3. Run the pipeline with the recommended steps
./rRNA-ops_pipeline.sh basecall,demultiplex,qc_before_filter,filter,qc_after_filter,blast_reads,metadata -c config.yaml -t 16

# 4. Check results
ls -lh results/taxonomy/

# 4. Import results to microecoR for downstream figures generation
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
| `demultiplex` | Porechop demultiplexing | Yes* |
| `qc_before_filter` | NanoPlot quality control before filtring reads | **Yes** |
| `filter` | Chopper read filtering | **Yes** |
| `qc_after_filter` | NanoPlot quality control after filtring reads | **Yes** |
| `blast_reads` | Pre-clustering BLAST | Optional** |
| `cluster` | Sequence clustering | Optional |
| `polish` | Medaka polishing | Optional |
| `taxonomy` | BLAST taxonomy assignment | Optional*** |
| `metadata` | Generate sample metadata | **Yes** |

*Skip if starting with basecalled/demultiplexed FASTQ files
**If not doing clustering, run `blast_reads` step
***If clustering is carried out, skip the `blast_reads` step and run the `taxonomy` step

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
./rRNA-ops_pipeline.sh basecall,demultiplex,qc_before_filter,filter,qc_after_filter,cluster,taxonomy,metadata -c config.yaml -t 16
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
#### Change the location of input pod5/fastq reads, output directory, basecalling, demultiplexing, filtering, clustering, and database paths in the config file

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
  blast_db: "/path/to/rOpDB.fasta" 
  evalue: 20
  min_identity: 50
```

## Examples

### Example 1: Complete Pipeline

```bash
# Starting from POD5 files
./rRNA-ops_pipeline.sh all -c config.yaml -t 32

```

### Example 2: Starting from basecalling (use pod5 files as input) and running clustering [choose clustering method in config file]

```bash
# Starting from POD5 files
./rRNA-ops_pipeline.sh basecall,demultiplex,qc_before_filter,filter,qc_after_filter,cluster,taxonomy,metadata -c config.yaml -t 16

```

### Example 3: Starting from basecalled fastq reads (use fastq files as input) and running clustering [choose clustering method in config file]

```bash
# Starting from basecalled FASTQ files
./rRNA-ops_pipeline.sh demultiplex,qc_before_filter,filter,qc_after_filter,cluster,taxonomy,metadata -c config.yaml -t 16

```
### Example 4: Skip Basecalling & Demultiplexing (use basecalled and demultiplexed fastq files as input) and running clustering [choose clustering method in config file]

```bash
# Starting from basecalled and demultiplexed FASTQ files
./rRNA-ops_pipeline.sh qc_before_filter,filter,qc_after_filter,cluster,taxonomy,metadata -c config.yaml -t 16 -c config.yaml -t 16
```

### Example 5: Skip clustering and using raw pod5 or basecalled fast reads (**recommended**)

```bash
# Starting from pod5 files
./rRNA-ops_pipeline.sh basecall,demultiplex,qc_before_filter,filter,qc_after_filter,blast_reads,metadata -c config.yaml -t 16

# Starting from basecalled FASTQ files
./rRNA-ops_pipeline.sh demultiplex,qc_before_filter,filter,qc_after_filter,blast_reads,metadata -c config.yaml -t 16

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

**Problem:** Memory errors
```bash
# Reduce thread count
# Process samples in batches
# Use more swap space
```

## Citation

If you use this pipeline in your research, please cite:

```
[Citation]
rRNA-ops Pipeline v1.0
https://github.com/anilchauhanhp9/rRNA-ops-pipeline.git
```

**Also cite the tools used:**

- **Dorado:** https://github.com/nanoporetech/dorado.git
- **Porechop:** https://github.com/rrwick/Porechop.git
- **Nanoplot:** https://github.com/wdecoster/nanoplot
- **MultiQC:** https://github.com/MultiQC/MultiQC; http://dx.doi.org/10.1093/bioinformatics/btw354
- **Chopper:** https://github.com/wdecoster/chopper
- **VSEARCH:** https://peerj.com/articles/2584/
- **Swarm:** https://peerj.com/articles/593/
- **Medaka:** https://github.com/nanoporetech/medaka.git
- **BLAST:** https://link.springer.com/article/10.1186/1471-2105-10-421

---

## Documentation

- [Installation Guide](INSTALLATION_GUIDE.md) - Detailed setup instructions
- [Pipeline Steps](PIPELINE_STEPS.md) - Detailed pipeline steps
---

## Acknowledgments

- Oxford Nanopore Technologies for sequencing technology
- rOpDB database maintainers
- Open-source bioinformatics community

### AI Assistance

During the development of this pipeline, **Claude (Anthropic)** was used as an
AI-assisted programming and documentation tool. Claude was used to assist with
code generation, debugging, workflow organization, and documentation drafting.
All AI-assisted output was reviewed, tested, and adapted by the author, who is
responsible for the final implementation, configuration, scientific decisions,
and release of the software.

---

## Version History

**v1.0** (2026-02-13)
- Initial release
- VSEARCH and Swarm clustering
- Full taxonomy assignment pipeline
- Multiple output formats

---