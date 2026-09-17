# rRNA-ops Pipeline - Detailed Pipeline Steps

## Table of Contents
1. [Pipeline Overview](#pipeline-overview)
2. [Step-by-Step Workflow](#step-by-step-workflow)
3. [Configuration Parameters](#configuration-parameters)
4. [Output Files](#output-files)
5. [Best Practices](#best-practices)

---

## Pipeline Overview

The rRNA-ops pipeline processes full-length bacterial rRNA operon sequences from Oxford Nanopore data, performing quality control, clustering, and taxonomic annotation.

**Pipeline Flow:**
```
POD5/FAST5 → Basecalling → Demultiplexing → Filtering → QC → 
Pre-clustering BLAST → Clustering → Polishing → Taxonomy Assignment → Results
```

---

## Step-by-Step Workflow

### Step 1: Basecalling with Dorado (Optional)

**Purpose:** Convert raw POD5/FAST5 signal data to FASTQ sequences.

**Command:**
```bash
./rRNA-ops_pipeline.sh basecall -c config.yaml -t 16
```

**What it does:**
1. Searches for POD5 files in the specified directory
2. Uses Dorado neural network basecaller to convert signals to sequences
3. Optionally trims adapters
4. Outputs FASTQ files

**Configuration:**
```yaml
basecalling:
  model: "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"  # Basecalling model
  device: "cuda:all"      # Use all GPUs (or "cpu" for CPU-only)
  no_trim: true           # Disable adapter trimming
  batch_size: 0           # Auto-detect optimal batch size
  recursive: true         # Search subdirectories
```

**Key Points:**
- GPU acceleration significantly speeds up basecalling (10-50x faster)
- Super accuracy (sup) models are slower but more accurate
- Model selection should match your flowcell and chemistry

**Output:**
- `output_dir/basecalled/` - FASTQ files
- Typically generates 1 FASTQ file per POD5 file

**Performance:**
- ~1M reads: 2-4 hours (GPU) or 24-48 hours (CPU)
- GPU memory: 8-16 GB recommended

**Skip this step if:** You already have basecalled FASTQ files.

---

### Step 2: Demultiplexing with Porechop (Optional)

**Purpose:** Split multiplexed samples by barcode and trim adapters.

**Command:**
```bash
./rRNA-ops_pipeline.sh demultiplex -c config.yaml -t 16
```

**What it does:**
1. Identifies barcode sequences at read ends
2. Assigns reads to barcode bins
3. Trims adapter sequences
4. Creates separate FASTQ files per barcode

**Configuration:**
```yaml
demultiplexing:
  barcode_kit: "SQK-RBK114-24"  # ONT barcode kit used
  barcode_threshold: 85          # Minimum barcode score (0-100)
```

**Key Points:**
- Higher thresholds (90+) = more stringent = fewer misassigned reads
- Lower thresholds (75-85) = more lenient = more total reads
- Porechop automatically detects barcode orientation

**Output:**
```
output_dir/demultiplexed/
├── barcode01.fastq.gz
├── barcode02.fastq.gz
├── ...
└── unclassified.fastq.gz  # Reads without clear barcodes
```

**Barcode Assignment Statistics:**
- Logs show reads per barcode
- Check `unclassified.fastq.gz` size - high numbers may indicate issues

**Skip this step if:** 
- Single sample (no multiplexing)
- Already demultiplexed data

---

### Step 3: Read Filtering with Chopper

**Purpose:** Filter reads by length and quality to retain only full-length, high-quality rRNA operons.

**Command:**
```bash
./rRNA-ops_pipeline.sh filter -c config.yaml -t 16
```

**What it does:**
1. Filters reads by minimum/maximum length
2. Filters by average quality score
3. Removes low-quality sequences
4. Compresses output

**Configuration:**
```yaml
filtering:
  min_length: 3500  # Minimum read length (bp)
  max_length: 6000  # Maximum read length (bp)
  min_quality: 10   # Minimum average Q-score
```

**Key Points:**
- Full bacterial rRNA operon: ~5000-5500 bp (16S-tRNAs-23S-5S)
- Adjust length based on target organisms
- Quality threshold Q10 is permissive; Q12-15 for stricter filtering

**Length Guidelines:**
| Organism Type | Min Length | Max Length | Notes |
|--------------|------------|------------|-------|
| Most bacteria | 3500 | 6000 | Standard operon |
| E. coli | 5000 | 5500 | Well-characterized |
| High GC | 3500 | 6500 | May be longer |

**Output:**
```
output_dir/filtered/
├── barcode01/
│   └── barcode01_filtered.fastq.gz
├── barcode02/
│   └── barcode02_filtered.fastq.gz
└── ...
```

**Quality Metrics:**
- Logs show reads retained vs. removed
- Typical retention: 60-80% of input reads

---

### Step 4: Quality Control with NanoPlot

**Purpose:** Generate comprehensive quality assessment reports.

**Command:**
```bash
./rRNA-ops_pipeline.sh qc -c config.yaml -t 16
```

**What it does:**
1. Analyzes read length distributions
2. Assesses quality score distributions
3. Generates summary statistics
4. Creates visualization plots

**Output:**
```
output_dir/qc/
├── barcode01/
│   ├── NanoStats.txt           # Summary statistics
│   ├── NanoPlot-report.html    # Interactive HTML report
│   ├── LengthvsQualityScatterPlot_dot.png
│   ├── LengthHistogram.png
│   └── ...
└── barcode02/
    └── ...
```

**Key Metrics to Check:**

1. **Read Length Distribution:**
   - Should peak around 5000-5500 bp
   - Minimal reads < 3500 bp or > 6000 bp

2. **Quality Score Distribution:**
   - Median Q-score should be > 10
   - High-quality runs: Q15+

3. **Read Count:**
   - Minimum 100 reads per sample for clustering
   - 1000+ reads recommended for robust analysis

4. **N50 Length:**
   - Should be close to expected operon length (~5200 bp)

**Interpreting Results:**
- Bimodal length distribution → May indicate contamination or multiple operons
- Low quality scores → Consider re-basecalling with better model
- High read count variation → Check demultiplexing efficiency

---

### Step 5: Pre-clustering BLAST (Optional)

**Purpose:** Perform initial taxonomic screening of individual reads before clustering.

**Command:**
```bash
./rRNA-ops_pipeline.sh blast_reads -c config.yaml -t 16
```

**What it does:**
1. BLASTs each filtered read against database
2. Generates per-read taxonomy assignments
3. Creates summary tables
4. Useful for contamination screening

**Output:**
```
output_dir/blast_reads/
├── barcode01_blast_results.txt
├── barcode01_blast_summary.txt
└── ...
```

**When to use:**
- Initial exploration of sample composition
- Contamination detection
- Quality verification before clustering

**Note:** This is separate from the final OTU-level taxonomy assignment in Step 8.

---

### Step 6: Sequence Clustering

**Purpose:** Group similar sequences into Operational Taxonomic Units (OTUs).

**Command:**
```bash
./rRNA-ops_pipeline.sh cluster -c config.yaml -t 16
```

**What it does:**
1. Performs chimera detection and removal
2. Clusters sequences based on similarity
3. Generates consensus sequences per OTU
4. Creates OTU abundance tables

**Two Clustering Algorithms Available:**

#### Option A: VSEARCH (Default)
```yaml
clustering:
  tool: "vsearch"
  cluster_identity: 0.97
  remove_chimeras: true
```

**VSEARCH Workflow:**
1. Dereplication (find unique sequences)
2. Chimera removal
3. Cluster generation

#### Option B: Swarm
```yaml
clustering:
  tool: "swarm"
  swarm_d: 1                    # Max differences (resolution)
  swarm_fastidious: true        # Link low-abundance clusters
  swarm_boundary: 3             # Min abundance for fastidious
  swarm_no_otu_breaking: false  # Allow OTU splitting
```

**Swarm Workflow:**
1. Chimera removal
2. Abundance-based clustering
3. Fastidious mode (if enabled)
4. OTU refinement

**Parameters:**
- `swarm_d: 1` - High resolution (strain level)
  - d=1: ~99.7% identity (strain/species level)
  - d=5: ~97% identity (species level)
  - d=13: ~95% identity (genus level)
- `swarm_fastidious: true` - Connects clusters via low-abundance intermediates
- `swarm_boundary: 3` - Prevents linking very rare sequences

**Choosing Between VSEARCH and Swarm:**

| Feature | VSEARCH | Swarm |
|---------|---------|-------|
| Speed | Fast | Very fast |
| Accuracy | High | Very high |
| Singletons | Removes (minsize) | Keeps (fastidious) |
| Best for | Standard analysis | Fine-scale diversity |
| OTU count | Moderate | Higher (more resolved) |

**Output:**
```
output_dir/cluster/
├── otus.fasta                    # OTU representative sequences
├── otu_per_sample.txt            # Abundance table
├── chimeras.fasta                # Detected chimeras
└── clustering_stats.txt          # Summary statistics
```

**OTU Table Format:**
```
OTU_ID    barcode01    barcode02    barcode03
otu_1     125          0            89
otu_2     45           230          15
otu_3     0            12           156
```

**Typical Results:**
- 1000 input reads → 10-50 OTUs (varies by sample diversity)
- Check chimera rate: 1-5% is normal; >10% indicates issues

---

### Step 7: Consensus Polishing with Medaka (Optional)

**Purpose:** Improve OTU consensus sequence accuracy using neural networks.

**Command:**
```bash
./rRNA-ops_pipeline.sh polish -c config.yaml -t 16
```

**What it does:**
1. Aligns reads to each OTU consensus
2. Uses Medaka neural network to call better consensus
3. Corrects systematic errors from basecalling
4. Improves accuracy for taxonomy assignment

**Configuration:**
```yaml
polishing:
  use_medaka: true                              # Enable polishing
  medaka_model: "r1041_e82_400bps_sup_v5.2.0"  # Must match basecalling
  parallel_polishing: false                     # Resource management
```

**Polishing Modes:**

**Sequential Mode (parallel_polishing: false - RECOMMENDED):**
- Processes OTUs one at a time
- Uses full thread count per OTU
- Memory efficient
- Example: 16 threads → each OTU uses 16 threads

**Parallel Mode (parallel_polishing: true):**
- Processes multiple OTUs simultaneously
- Uses fewer threads per OTU
- Faster but higher CPU usage
- Example: 16 threads → runs 48 parallel jobs of 1 thread

**Model Selection:**
Critical: Model must match basecalling chemistry!

| Flowcell | Chemistry | Model |
|----------|-----------|-------|
| R10.4.1 | E8.2 | r1041_e82_400bps_sup_v5.2.0 |
| R9.4.1 | various | r941_min_high_g360 |

**Output:**
```
output_dir/polished/
├── polished_otus.fasta           # Polished sequences
└── polishing_summary.txt         # Polishing statistics
```

**When to use:**
- ✓ High-accuracy taxonomy needed
- ✓ Downstream phylogenetic analysis
- ✗ Skip if time/resources limited
- ✗ Skip if clustering quality is already excellent

**Performance Impact:**
- Time: Adds 30-60 minutes for 50 OTUs
- Accuracy improvement: 0.5-2% error rate reduction

---

### Step 8: Taxonomic Assignment with BLAST

**Purpose:** Assign taxonomy to each OTU using BLAST against rOpDB.

**Command:**
```bash
./rRNA-ops_pipeline.sh taxonomy -c config.yaml -t 16
```

**What it does:**
1. BLASTs each OTU against reference database
2. Parses taxonomy from best hit
3. Generates multiple output formats
4. Creates filtered abundance tables

**Configuration:**
```yaml
taxonomy:
  blast_db: "/path/to/rOpDB/rOpDB.fasta"
  evalue: 20              # E-value threshold
  max_target_seqs: 1      # Top hit only
  min_identity: 50        # Minimum % identity
```

**BLAST Parameters:**

- **E-value (20):** 
  - Permissive threshold for distantly related sequences
  - Use 1e-10 for more stringent matches

- **Min Identity (50%):**
  - Captures distant relatives
  - Increase to 90-97% for species-level assignments

- **Database Format:**
  Headers must be: `accession|Bacteria;Proteobacteria;...;Genus;Species`

**Output Files:**

1. **blast_results.txt** - Raw BLAST output
   ```
   otu_1    accession123    95.2    500    ...
   ```

2. **taxonomy.tsv** - All OTUs with taxonomy
   ```
   OTU_ID    Accession    Identity    Length    Taxonomy
   otu_1     ACC123       95.2        4892      Bacteria;Proteobacteria;...
   otu_2     Unassigned   0           0         Unassigned
   ```

3. **taxonomy_otus_hits_only.tsv** - Only OTUs with BLAST hits
   - Excludes unassigned OTUs

4. **taxonomy_otus_expanded.tsv** - Separate rank columns
   ```
   OTU_ID    Kingdom      Phylum           Class           Order    ...
   otu_1     Bacteria     Proteobacteria   Gammaproteobac  Entero   ...
   ```

5. **taxonomy_otus_expanded_hits_only.tsv** - Expanded format, hits only

6. **taxonomy_otus_for_microeco.tsv** - Microeco R package format
   - OTU_ID + 8 taxonomy columns only
   - Ready for `tax_table()` function

7. **otu_per_sample_with_taxonomy.txt** - Filtered abundance table
   - Only includes OTUs with BLAST hits
   - Removes unassigned OTUs from abundance data

**Taxonomy Format:**
```
Domain;Phylum;Class;Order;Family;Genus;Species;Strain
```

**Typical Assignment Rates:**
- Well-studied organisms: 80-95% OTUs assigned
- Environmental samples: 50-80% OTUs assigned
- Novel organisms: 30-50% OTUs assigned

**Interpreting Results:**

**High assignment rate (>80%):**
- ✓ Good sample quality
- ✓ Representative database coverage
- ✓ Expected organisms present

**Low assignment rate (<50%):**
- ⚠ Novel/uncharacterized organisms
- ⚠ Database gaps
- ⚠ Possible contamination
- ⚠ Check BLAST parameters

---

### Step 9: Metadata Generation

**Purpose:** Create sample metadata file for downstream analysis.

**Command:**
```bash
./rRNA-ops_pipeline.sh metadata -c config.yaml
```

**What it does:**
1. Scans output directories
2. Collects sample statistics
3. Generates metadata table
4. Ready for R/Python analysis

**Output:**
```
output_dir/metadata.tsv
```

**Format:**
```
SampleID    ReadCount    FilteredReads    OTUs    AssignedOTUs
barcode01   5234         4123             23      20
barcode02   8901         7234             45      41
...
```

**Use cases:**
- Import into R for statistical analysis
- Batch effect correction
- Quality filtering decisions

---

## Configuration Parameters

### Complete config.yaml Template

```yaml
# Input/Output Paths
input:
  pod5_dir: "/path/to/pod5_files"
  fastq_dir: "/path/to/fastq_files"

output:
  output_dir: "/path/to/results"

# Basecalling (optional)
basecalling:
  model: "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
  device: "cuda:all"
  no_trim: true
  batch_size: 0
  recursive: true

# Demultiplexing (optional)
demultiplexing:
  barcode_kit: "SQK-RBK114-24"
  barcode_threshold: 85

# Read Filtering
filtering:
  min_length: 3500
  max_length: 6000
  min_quality: 10

# Clustering
clustering:
  tool: "swarm"                    # "vsearch" or "swarm"
    
  # SWARM parameters
  swarm_d: 1
  swarm_fastidious: true
  swarm_boundary: 3
  swarm_no_otu_breaking: false
  
  # Common
  remove_chimeras: true

# Polishing (optional)
polishing:
  use_medaka: false
  medaka_model: "r1041_e82_400bps_sup_v5.2.0"
  parallel_polishing: false

# Taxonomy
taxonomy:
  blast_db: "/path/to/rOpDB.fasta"
  evalue: 20
  max_target_seqs: 1
  min_identity: 50
```

---

## Output Files

### Directory Structure

```
output_dir/
├── basecalled/              # Raw FASTQ from basecalling
├── demultiplexed/           # Separated by barcode
│   ├── barcode01.fastq.gz
│   └── ...
├── filtered/                # Quality-filtered reads
│   ├── barcode01/
│   └── ...
├── qc/                      # NanoPlot reports
│   ├── barcode01/
│   │   └── NanoPlot-report.html
│   └── ...
├── blast_reads/             # Pre-clustering BLAST (optional)
├── cluster/                 # Clustering results
│   ├── otus.fasta
│   ├── otu_per_sample.txt
│   └── chimeras.fasta
├── polished/               # Polished OTUs (optional)
│   └── polished_otus.fasta
├── taxonomy/               # Taxonomy assignments
│   ├── blast_results.txt
│   ├── taxonomy.tsv
│   ├── taxonomy_otus_expanded.tsv
│   ├── taxonomy_otus_for_microeco.tsv
│   └── otu_per_sample_with_taxonomy.txt
└── metadata.tsv            # Sample metadata
```

---

## Best Practices

### 1. Data Quality
- Always run QC step to assess read quality
- Check filtering statistics - aim for >60% retention
- Review chimera detection rates

### 2. Clustering Strategy
- **For microbial ecology:** Use Swarm with fastidious mode
- **For strain typing:** Use Swarm with d=1
- **For general analysis:** VSEARCH is faster and robust

### 3. Resource Management
- Basecalling: Use GPU when available
- Polishing: Use sequential mode unless you have high CPU count
- BLAST: More threads = faster (diminishing returns >16 threads)

### 4. Database Curation
- Keep rOpDB updated
- Consider building custom database for specific environments
- Validate taxonomy strings in headers

### 5. Reproducibility
- Save config.yaml with results
- Document pipeline version
- Record tool versions used

### 6. Troubleshooting
- Low OTU count: Check filtering stringency
- High chimera rate: Review PCR conditions or filtering
- Poor BLAST hits: Check database coverage for your organisms

---

## Running Multiple Steps

**Single step:**
```bash
./rRNA-ops_pipeline.sh filter -c config.yaml -t 16
```

**Multiple steps (comma-separated, no spaces):**
```bash
./rRNA-ops_pipeline.sh filter,qc,cluster,taxonomy -c config.yaml -t 16
```

**Complete pipeline:**
```bash
./rRNA-ops_pipeline.sh all -c config.yaml -t 16
```

---

## Performance Benchmarks

| Step | 1000 Reads | 10000 Reads | Threads | Notes |
|------|-----------|-------------|---------|-------|
| Filtering | 2 min | 15 min | 4 | Fast |
| QC | 5 min | 30 min | 4 | Visualization |
| Clustering | 10 min | 60 min | 16 | CPU intensive |
| Polishing | 30 min | 180 min | 16 | Optional |
| Taxonomy | 5 min | 20 min | 16 | BLAST speed |

