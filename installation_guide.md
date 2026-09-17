# Installation Guide for rRNAops_pipeline

Complete installation instructions for all pipeline components.

## Prerequisites

- **Operating System**: Linux or macOS
- **Conda/Miniconda**: [Install here](https://docs.conda.io/en/latest/miniconda.html)
- **Disk Space**: At least 10 GB free
- **RAM**: Minimum 8 GB (16 GB recommended)

## Step 1: Clone Repository

```bash
git clone https://github.com/yourusername/rRNAops_pipeline.git
cd rRNAops_pipeline
```

## Step 2: Create Conda Environment

```bash
# Create environment from file
conda env create -f rrna_conda_env.yml

# Activate environment
conda activate rRNAops-pipeline

# Verify installation
conda list
```

**Note**: This installs all tools except Dorado

## Step 3: Install Dorado (for Basecalling)

**Note**: Dorado is only needed if you're starting from POD5 files. Skip this if you already have basecalled FASTQ files.
**Note**: Download the latest version from the Dorado GitHub webpage

#### For Linux (x86_64)

```bash
# Download latest release
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64.tar.gz

# Extract
tar -xzf dorado-0.5.3-linux-x64.tar.gz

# Move to installation directory
sudo mv dorado-0.5.3-linux-x64/bin/dorado /usr/local/bin/

# OR add to your PATH in ~/.bashrc
echo 'export PATH="/path/to/dorado-0.5.3-linux-x64/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Verify installation
dorado --version
```

#### For macOS (Apple Silicon)

```bash
# Download for Apple Silicon
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-osx-arm64.tar.gz

# Extract
tar -xzf dorado-0.5.3-osx-arm64.tar.gz

# Move to installation directory
sudo mv dorado-0.5.3-osx-arm64/bin/dorado /usr/local/bin/

# Verify installation
dorado --version
```

#### For macOS (Intel)

```bash
# Download for Intel Mac
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-osx-x64.tar.gz

# Extract and install
tar -xzf dorado-0.5.3-osx-x64.tar.gz
sudo mv dorado-0.5.3-osx-x64/bin/dorado /usr/local/bin/

# Verify installation
dorado --version
```
### Download Basecalling Models

```bash
# List available models
dorado download --list

# Download a model (example for R10.4.1 super accuracy) in the diectory where rRNAop-pipeline is running
dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.2.0
```

## Step 4: Verify Installation

Run the verification script:

```bash
# Make pipeline executable
chmod +x rRNA-ops_pipeline.sh

# Check version
rRNA-ops_pipeline.sh --version
```

## Step 5: Database
```bash
# download the database from https://www.njmicrobe.org/ropdb
# make a blast database with "makeblastdb -in sequences.fasta -dbtype nucl -out rOpDB -parse_seqids"
# give the location of the created DB in the blast_db option as blast_db: "location_of_DB/rOpDB.fasta" 
```
## Troubleshooting

### Issue: Conda environment creation fails

**Solution 1**: Update conda
```bash
conda update -n base -c defaults conda
```

**Solution 2**: Try mamba (faster resolver)
```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yml
```
### Issue: Dorado not found

**Check PATH**:
```bash
echo $PATH
which dorado

# Add to PATH if needed
export PATH="/path/to/dorado/bin:$PATH"
```

**Verify binary is executable**:
```bash
chmod +x /path/to/dorado
```
### Issue: Permission denied

```bash
# Make script executable
chmod +x rRNA-ops_pipeline.sh

# If installing to system directories, use sudo
sudo mv tool /usr/local/bin/
```
### Issue: GPU support for Dorado

For GPU acceleration:
```bash
# Check CUDA availability
nvidia-smi

# Download GPU version of Dorado
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.3-linux-x64-cuda.tar.gz

# Use GPU for basecalling
dorado basecaller --device cuda:0 ...
```

## Platform-Specific Notes

### Linux

Most tools install without issues. Ensure you have:
```bash
sudo apt-get update
sudo apt-get install -y build-essential git wget curl
```

### macOS

May need Xcode command line tools:
```bash
xcode-select --install
```

For M1/M2 Macs:
- Use ARM64 versions of tools
- Some bioconda packages may require Rosetta 2

### Windows (WSL2)

Install WSL2 first:
1. Install WSL2: https://docs.microsoft.com/en-us/windows/wsl/install
2. Install Ubuntu 22.04 from Microsoft Store
3. Follow Linux installation instructions inside WSL2

## System Requirements Summary

### Minimum System
- 4 CPU cores
- 8 GB RAM
- 50 GB free disk space
- Linux/macOS

### Recommended System
- 16+ CPU cores
- 32 GB RAM
- 500 GB free disk space (SSD)
- Linux with GPU (for Dorado)

### For Large Datasets (>10M reads)
- 32+ CPU cores
- 64+ GB RAM
- 1+ TB disk space (SSD)
- High-memory cluster node

## Next Steps

After successful installation:

1. **Read documentation**: `README.md`
2. **Configure pipeline**: Edit `config.yaml`
3. **Run test**: Use small dataset first

---

**Installation complete!** You're ready to analyze rRNA operons! 🧬
