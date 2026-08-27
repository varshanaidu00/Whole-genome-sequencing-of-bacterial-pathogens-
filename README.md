Bacterial Whole-Genome Sequencing Analysis Pipeline
=================================================================================================================

## Background

A Snakemake-based bacterial whole-genome sequencing (WGS) analysis pipeline for paired-end Illumina FASTQ data. The workflow performs read trimming, read quality assessment, taxonomic classification, de novo genome assembly, assembly quality assessment, genome completeness assessment, read-mapping coverage, MLST typing, antimicrobial resistance (AMR) detection, and virulence gene detection.
The pipeline is designed so that each bioinformatics tool is run in its own Conda environment, providing reproducibility and simplifying software dependency management.

Instructions on how to create and execute a basic workflow in Snakemake can be found in this [wiki](https://github.com/varshanaidu00/Whole-genome-sequencing-of-bacterial-pathogens-/wiki/).

### Overview

The pipeline takes paired-end raw FASTQ files as input and produces:

- Quality control metrics from `FastQC`
- Adapter and quality-trimmed reads using `Trimmomatic`
- Taxonomic classification using `Kraken2`
- De novo genome assemblies using `Unicycler`
- Assembly statistics using `QUAST`
- Genome completeness using `BUSCO`
- Read depth/coverage using `CoverM`
- Sequence typing using `MLST`
- Antimicrobial resistance gene detection using `AMRFinderPlus`
- Virulence gene detection using `VirulenceFinder`

A final tabulated summary containing key metrics for every sample

```
                     Raw paired-end FASTQ
                              |
                              v
                     +----------------+
                     |  Trimmomatic   |
                     +-------+--------+
                             |
                     Trimmed paired reads
                             |
             +---------------+----------------+
             |               |                |
             v               v                v
         +-------+       +---------+      +----------+
         |FastQC |       | Kraken2 |      | Unicycler|
         +-------+       +---------+      +----+-----+
                                               |
                                        assembly.fasta
                                               |
                  +------------+---------------+---------------+
                  |            |               |               |
                  v            v               v               v
                QUAST        BUSCO           CoverM           MLST
                  |            |               |               |
                  +------------+---------------+---------------+
                                               |
                                  +------------+-------------+
                                  |                          |
                                  v                          v
                             AMRFinderPlus            VirulenceFinder
                                  |                          |
                                  +------------+-------------+
                                               |
                                               v
                                      WGS_summary.tsv
```
### Repository structure
```
bacterial-wgs-snakemake/
├── Snakefile
├── config.yaml
├── samples.tsv
├── environment.yml
├── README.md
├── .gitignore
│
├── envs/
│   ├── fastqc.yaml
│   ├── trimmomatic.yaml
│   ├── kraken2.yaml
│   ├── unicycler.yaml
│   ├── quast.yaml
│   ├── busco.yaml
│   ├── coverm.yaml
│   ├── mlst.yaml
│   ├── amrfinderplus.yaml
│   ├── virulencefinder.yaml
│   └── summary.yaml
│
├── scripts/
│   └── make_summary.py
│
├── data/
├── results/
└── logs/
```
---

## Requirements

Before running the pipeline, ensure that the following are available:

* [ ] Snakemake (≥6.0)
* [ ] Conda or Mamba
* [ ] All required Conda environments (created via provided `envs/` files)
* [ ] Kraken2 database (Standard or Standard-16)
* [ ] VirulenceFinder database
* [ ] AMRFinderPlus database
* [ ] BUSCO lineage dataset (bacteria_odb12 or similar)
* [ ] Paired-end FASTQ files for every sample
* [ ] Appropriate local databases for Kraken2, BUSCO, AMRFinderPlus and VirulenceFinder
* [ ] `Snakefile`
* [ ] `config.yaml`
* [ ] `samples.tsv` (sample metadata file)

The pipeline assumes that the required software environments and databases already exist and are correctly configured in `config.yaml`.

---

## Installation
Create the Snakemake environment:
```
conda env create -f environment.yml
```
Activate it:
```
conda activate bacterial-wgs-snakemake
```
Verify:
```
snakemake --version
```
## Input data
Paired-end FASTQ files should be stored separately from the workflow repository.

Example:
```
project/
├── data/
│   ├── sample01_R1.fastq.gz
│   ├── sample01_R2.fastq.gz
│   ├── sample02_R1.fastq.gz
│   └── sample02_R2.fastq.gz
│
└── bacterial-wgs-snakemake/
```
## Sample sheet
Create samples.tsv:
```
sample	R1	R2
sample01	data/sample01_R1.fastq.gz	data/sample01_R2.fastq.gz
sample02	data/sample02_R1.fastq.gz	data/sample02_R2.fastq.gz
```
The three required columns are:

| Column | Description |
|--------|-------------|
| `Sample`| Unique sample identifier |
| `R1`| Forward FASTQ |
| `R2`| Reverse FASTQ |

Sample identifiers should avoid spaces and shell/path special characters.

Paired-end sequencing reads must be registered in a tab-separated sample metadata file specified by the `samples` parameter in `config.yaml`.

### Creating the samples.tsv File

Create a tab-separated values (TSV) file with exactly three columns: `sample`, `R1`, and `R2`.

**Column definitions:**
- `sample`: Sample identifier (used in output directory names)
- `R1`: Full or relative path to the forward-read FASTQ file
- `R2`: Full or relative path to the reverse-read FASTQ file

Files can be gzip-compressed (`.fastq.gz`) or uncompressed (`.fastq`).

**Example `samples.tsv`:**

```
sample	R1	R2
isolate01	reads/isolate01_R1.fastq.gz	reads/isolate01_R2.fastq.gz
isolate02	reads/isolate02_R1.fastq.gz	reads/isolate02_R2.fastq.gz
isolate03	/absolute/path/to/isolate03_1.fastq	/absolute/path/to/isolate03_2.fastq
```

**Important notes:**
- Column headers are required and must be exactly as shown.
- Paths can be absolute or relative to the directory where Snakemake is executed.
- All files must exist before running the pipeline.
- Sample IDs must be unique and suitable for use in directory names (alphanumeric, underscore, hyphen).

---

## 3. Configuration

The pipeline reads configuration values from `config.yaml`. This file specifies directories, databases, tool parameters, and per-tool thread allocations.

### Example config.yaml

```yaml
samples: "samples.tsv"

results: "results"
logs: "logs"

# Databases
kraken2_db: "/path/to/kraken2/database"
busco_lineage: "bacteria_odb12"
virulencefinder_db: "/path/to/virulencefinder_db"

# Optional reference genome for QUAST
# Set to null for de novo-only analysis (no genome fraction metric)
reference: null

# Adapter file for Trimmomatic
adapters: "/path/to/TruSeq3-PE.fa"

# Trimmomatic parameters
trimmomatic:
  leading: 3
  trailing: 3
  slidingwindow: "4:20"
  minlen: 36

# Kraken2 parameters
kraken2:
  confidence: 0.0

# CoverM parameters
coverm:
  min_identity: 95
  min_aligned_percent: 50

# VirulenceFinder parameters
virulencefinder:
  min_coverage: 0.60
  min_identity: 0.90

# Thread allocation per tool
threads:
  fastqc: 4
  trimmomatic: 8
  kraken2: 16
  unicycler: 8
  quast: 4
  busco: 8
  coverm: 8
  mlst: 2
  amrfinderplus: 8
  virulencefinder: 8
```

### Configuration Parameter Reference

| Parameter                    | Type             | Default         | Description                                    |
|------------------------------|------------------|-----------------|------------------------------------------------|
| `samples`                    | string (path)    | —               | Path to sample metadata TSV file               |
| `results`                    | string (path)    | `results`       | Output directory for results                   |
| `logs`                       | string (path)    | `logs`          | Directory for log files                        |
| `kraken2_db`                 | string (path)    | —               | Path to Kraken2 database                       |
| `busco_lineage`              | string           | `bacteria_odb12`| BUSCO lineage dataset                          |
| `virulencefinder_db`         | string (path)    | —               | Path to VirulenceFinder database               |
| `reference`                  | string (path) or null | `null`      | Reference genome for QUAST (optional)          |
| `adapters`                   | string (path)    | —               | Adapter FASTA file for Trimmomatic             |
| `trimmomatic.leading`        | integer          | `3`             | Remove leading bases below quality             |
| `trimmomatic.trailing`       | integer          | `3`             | Remove trailing bases below quality            |
| `trimmomatic.slidingwindow`  | string           | `"4:20"`        | Sliding window quality threshold                |
| `trimmomatic.minlen`         | integer          | `36`            | Minimum read length after trimming             |
| `kraken2.confidence`         | float (0-1)      | `0.0`           | Kraken2 confidence threshold                   |
| `coverm.min_identity`        | integer (0-100)  | `95`            | Minimum sequence identity (%)                  |
| `coverm.min_aligned_percent` | integer (0-100)  | `50`            | Minimum aligned length (%)                     |
| `virulencefinder.min_coverage`| float (0-1)     | `0.60`          | Minimum coverage threshold                     |
| `virulencefinder.min_identity`| float (0-1)     | `0.90`          | Minimum identity threshold                     |
| `threads.*`                  | integer          | Per-tool        | CPU threads for each tool                      |

---

## 4. Running the Pipeline

### Quick Start

Run the complete workflow with:

```bash
snakemake --use-conda -j 8 --configfile config.yaml
```

**Parameter explanations:**
- `--use-conda` — Activates Conda environments for each rule.
- `-j 8` — Run up to 8 jobs concurrently (adjust based on your CPU cores).
- `--configfile config.yaml` — Specifies the configuration file.

### Before Running: Dry Run

Always perform a dry run first to verify the pipeline will execute correctly:

```bash
snakemake -n --use-conda -j 8 --configfile config.yaml
```

This shows what jobs Snakemake intends to run **without executing them**.

### Verbose Execution

For detailed command execution information:

```bash
snakemake -p --use-conda -j 8 --configfile config.yaml
```

The `-p` flag prints the shell commands as they execute, useful for debugging.

### Generate a DAG Visualization

To visualize the workflow as a directed acyclic graph (DAG):

```bash
snakemake --dag --use-conda --configfile config.yaml | dot -Tpdf > pipeline_dag.pdf
```

This requires Graphviz to be installed.

---

## 5. Pipeline Workflow Overview

The workflow follows this execution order:

```
Paired-end FASTQ reads
        │
        ├─ Trimmomatic ────────► Quality trimming & adapter removal
        │
        ├─ FastQC ─────────────► Read QC statistics
        │
        ├─ Kraken2 ────────────► Taxonomic classification
        │
        └─ Unicycler ──────────► De novo genome assembly
                                  │
                                  ├─ QUAST ──────────────► Assembly quality metrics
                                  │
                                  ├─ BUSCO ──────────────► Assembly completeness
                                  │
                                  ├─ CoverM ─────────────► Sequencing depth
                                  │
                                  ├─ MLST ───────────────► Sequence typing
                                  │
                                  ├─ AMRFinderPlus ──────► AMR gene detection
                                  │
                                  └─ VirulenceFinder ────► Virulence gene screening
                                  
                    All results
                        │
                        └─ Summary rule ────────────────► WGS_summary.tsv
```

---

## 6. Step 1 — Trimmomatic

Trimmomatic removes adapter sequences and performs quality-based trimming of paired-end reads.

**Inputs:**
- Forward and reverse FASTQ files from `samples.tsv`

**Outputs:**
```
results/trimmed/{sample}_R1_paired.fastq.gz
results/trimmed/{sample}_R2_paired.fastq.gz
results/trimmed/{sample}_R1_unpaired.fastq.gz  (discarded after FastQC)
results/trimmed/{sample}_R2_unpaired.fastq.gz  (discarded after FastQC)
```

**Log:**
```
logs/trimmomatic/{sample}.log
```

**Parameters (from config.yaml):**
- `adapters` — Path to adapter sequence file (e.g., TruSeq3-PE.fa)
- `leading` — Remove leading bases below quality threshold (default: 3)
- `trailing` — Remove trailing bases below quality threshold (default: 3)
- `slidingwindow` — Window size and quality (default: 4:20, meaning 4bp windows with avg Q≥20)
- `minlen` — Discard reads shorter than this after trimming (default: 36bp)

**Conda environment:** `envs/trimmomatic.yaml`

---

## 7. Step 2 — FastQC

FastQC generates read quality statistics on the trimmed paired-end reads.

**Inputs:**
```
results/trimmed/{sample}_R1_paired.fastq.gz
results/trimmed/{sample}_R2_paired.fastq.gz
```

**Outputs:**
```
results/fastqc/{sample}_R1_paired_fastqc.html
results/fastqc/{sample}_R2_paired_fastqc.html
results/fastqc/{sample}_R1_paired_fastqc/fastqc_data.txt
results/fastqc/{sample}_R2_paired_fastqc/fastqc_data.txt
```

**Log:**
```
logs/fastqc/{sample}.log
```

The `fastqc_data.txt` files are parsed by the summary rule to extract:
- Total reads
- Duplicate percentage
- Average quality score
- Q30 percentage
- GC content
- Sequence length

**Conda environment:** `envs/fastqc.yaml`

---

## 8. Step 3 — Kraken2

Kraken2 performs taxonomic classification of trimmed reads using k-mer matching.

**Inputs:**
```
results/trimmed/{sample}_R1_paired.fastq.gz
results/trimmed/{sample}_R2_paired.fastq.gz
```

**Outputs:**
```
results/kraken2/{sample}.kraken2.txt    (classification details)
results/kraken2/{sample}.report.txt     (summary report)
```

**Log:**
```
logs/kraken2/{sample}.log
```

The report is parsed by the summary rule to extract:
- Percentage of reads classified vs. unclassified
- Top taxon at species rank (S)

**Parameters:**
- `kraken2_db` — Path to the Kraken2 database
- `kraken2.confidence` — Confidence threshold (default: 0.0)

**Conda environment:** `envs/kraken2.yaml`

---

## 9. Step 4 — Unicycler

Unicycler performs de novo hybrid genome assembly using SPAdes.

**Inputs:**
```
results/trimmed/{sample}_R1_paired.fastq.gz
results/trimmed/{sample}_R2_paired.fastq.gz
```

**Outputs:**
```
results/assembly/{sample}/assembly.fasta
results/assembly/{sample}/assembly.gfa    (graphical fragment assembly, intermediate)
```

**Log:**
```
logs/unicycler/{sample}.log
```

The assembly FASTA file is used as input for all downstream analyses (QUAST, BUSCO, CoverM, MLST, AMRFinderPlus, VirulenceFinder).

**Parameters:**
- `threads` — CPU threads (default: 8)
- Mode is fixed to `normal` (medium computational cost, good accuracy)

**Conda environment:** `envs/unicycler.yaml`

---

## 10. Step 5 — QUAST

QUAST evaluates the quality and characteristics of each genome assembly.

**Inputs:**
```
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/quast/{sample}/report.tsv
```

**Log:**
```
logs/quast/{sample}.log
```

**Metrics extracted by the summary rule:**
- Genome fraction (%) — Only if a reference genome is provided
- Assembly N50 in kilobases (converted from QUAST's base pair output)
- Total assembly size in base pairs
- Number of contigs
- GC percentage

**Parameters:**
- `reference` — Optional reference genome path (if provided, enables genome fraction metric)
- `threads` — CPU threads (default: 4)

**Conda environment:** `envs/quast.yaml`

---

## 11. Step 6 — BUSCO

BUSCO evaluates genome completeness by identifying conserved single-copy orthologs.

**Inputs:**
```
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/busco/{sample}/short_summary.txt
```

**Log:**
```
logs/busco/{sample}.log
```

The completeness percentage (C:X.X%) is extracted and added to the summary.

**Parameters:**
- `busco_lineage` — Lineage dataset (default: bacteria_odb12)
- `threads` — CPU threads (default: 8)

**Conda environment:** `envs/busco.yaml`

---

## 12. Step 7 — CoverM

CoverM calculates sequencing depth and coverage metrics by mapping trimmed reads back to the assembly.

**Inputs:**
```
results/trimmed/{sample}_R1_paired.fastq.gz
results/trimmed/{sample}_R2_paired.fastq.gz
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/coverm/{sample}.tsv
```

**Log:**
```
logs/coverm/{sample}.log
```

The mean sequencing depth is extracted and added to the summary as "Depth_of_Coverage".

**Parameters:**
- `coverm.min_identity` — Minimum sequence identity threshold (default: 95%)
- `coverm.min_aligned_percent` — Minimum aligned length threshold (default: 50%)
- `threads` — CPU threads (default: 8)

**Methods used:** `mean` and `covered_fraction`

**Conda environment:** `envs/coverm.yaml`

---

## 13. Step 8 — MLST

MLST performs sequence typing using PubMLST typing schemes.

**Inputs:**
```
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/mlst/{sample}.tsv
```

**Log:**
```
logs/mlst/{sample}.log
```

The sequence type (ST) is extracted from the output and added to the summary.

**Conda environment:** `envs/mlst.yaml`

---

## 14. Step 9 — AMRFinderPlus

AMRFinderPlus identifies antimicrobial resistance genes from the assembled genome.

**Inputs:**
```
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/amrfinderplus/{sample}.tsv
```

**Log:**
```
logs/amrfinderplus/{sample}.log
```

Gene symbols are extracted and formatted as a semicolon-separated list in the summary.

**Parameters:**
- `threads` — CPU threads (default: 8)

**Conda environment:** `envs/amrfinderplus.yaml`

---

## 15. Step 10 — VirulenceFinder

VirulenceFinder searches the assembled genome for virulence-associated genes.

**Inputs:**
```
results/assembly/{sample}/assembly.fasta
```

**Outputs:**
```
results/virulencefinder/{sample}/results_tab.tsv
```

**Log:**
```
logs/virulencefinder/{sample}.log
```

Virulence factor names are extracted and formatted as a semicolon-separated list in the summary.

**Parameters:**
- `virulencefinder_db` — Path to VirulenceFinder database
- `virulencefinder.min_coverage` — Minimum coverage threshold (default: 0.60)
- `virulencefinder.min_identity` — Minimum identity threshold (default: 0.90)
- `threads` — CPU threads (default: 8)

**Conda environment:** `envs/virulencefinder.yaml`

---

## 16. Final Summary

The `summary` rule aggregates results from all preceding steps into a single tab-separated summary table.

**Output:**
```
results/WGS_summary.tsv
```

**Processing:**
The `scripts/make_summary.py` script parses output files from each tool and combines metrics into one row per sample. If a file is missing or cannot be parsed, the corresponding column is set to `NA`.

### Summary Table Columns

| Column                      | Source          | Description                                  |
|-----------------------------|-----------------|----------------------------------------------|
| `Sample`                    | Input           | Sample identifier from samples.tsv           |
| `Total_Reads`               | FastQC R1       | Total sequence count (paired-end R1)         |
| `Unique_Reads`              | FastQC R1       | Estimated unique sequences after dedup       |
| `Duplicate_Reads`           | FastQC R1       | Estimated number of duplicates               |
| `Pct_Duplicates(%)`         | FastQC R1       | Percentage of reads that are duplicates      |
| `Avg_Q_Score`               | FastQC R1       | Average Phred quality score                  |
| `Pct_Q30(%)`                | FastQC R1       | Percentage of bases with Q≥30               |
| `GC_Content(%)`             | FastQC R1       | GC content of reads                          |
| `Sequence_Length`           | FastQC R1       | Read length                                  |
| `Kraken2_Classified(%)`     | Kraken2         | Percentage of reads assigned to a taxon      |
| `Kraken2_Unclassified(%)` | Kraken2         | Percentage of unclassified reads             |
| `Kraken2_Top_Taxon`         | Kraken2         | Most abundant species (rank S)               |
| `Genome_Fraction(%)`        | QUAST           | Fraction of reference covered (if ref provided) |
| `Assembly_N50(kbp)`         | QUAST           | Assembly N50 in kilobases                    |
| `Assembly_Size(bp)`         | QUAST           | Total assembly length in base pairs          |
| `Number_of_Contigs`         | QUAST           | Number of contigs                            |
| `Percent_GC(%)`             | QUAST           | GC content of assembly                       |
| `BUSCO_Completeness`        | BUSCO           | BUSCO complete percentage (single-copy)      |
| `Depth_of_Coverage`         | CoverM          | Mean sequencing depth across assembly        |
| `ST`                        | MLST            | Sequence type from PubMLST scheme            |
| `AMR_Genes`                 | AMRFinderPlus   | Semicolon-separated list of AMR genes found  |
| `Virulence_Genes`           | VirulenceFinder | Semicolon-separated list of virulence genes  |

**Missing values:** If a tool's output file is missing or cannot be parsed, the corresponding cell is populated with `NA`.

**Conda environment:** `envs/summary.yaml`

---

## 17. Expected Output Structure

For two example samples, the output directory should resemble:

```
results/
├── fastqc/
│   ├── isolate01_R1_paired_fastqc.html
│   ├── isolate01_R1_paired_fastqc.zip
│   ├── isolate01_R1_paired_fastqc/
│   │   └── fastqc_data.txt
│   ├── isolate01_R2_paired_fastqc.html
│   ├── isolate01_R2_paired_fastqc.zip
│   ├── isolate01_R2_paired_fastqc/
│   │   └── fastqc_data.txt
│   ├── isolate02_R1_paired_fastqc.html
│   ├── isolate02_R1_paired_fastqc.zip
│   ├── isolate02_R1_paired_fastqc/
│   │   └── fastqc_data.txt
│   ├── isolate02_R2_paired_fastqc.html
│   ├── isolate02_R2_paired_fastqc.zip
│   └── isolate02_R2_paired_fastqc/
│       └── fastqc_data.txt
│
├── trimmed/
│   ├── isolate01_R1_paired.fastq.gz
│   ├── isolate01_R1_unpaired.fastq.gz
│   ├── isolate01_R2_paired.fastq.gz
│   ├── isolate01_R2_unpaired.fastq.gz
│   ├── isolate02_R1_paired.fastq.gz
│   ├── isolate02_R1_unpaired.fastq.gz
│   ├── isolate02_R2_paired.fastq.gz
│   └── isolate02_R2_unpaired.fastq.gz
│
├── kraken2/
│   ├── isolate01.kraken2.txt
│   ├── isolate01.report.txt
│   ├── isolate02.kraken2.txt
│   └── isolate02.report.txt
│
├── assembly/
│   ├── isolate01/
│   │   ├── assembly.fasta
│   │   └── assembly.gfa
│   └── isolate02/
│       ├── assembly.fasta
│       └── assembly.gfa
│
├── quast/
│   ├── isolate01/
│   │   ├── report.tsv
│   │   └── ... (other QUAST outputs)
│   └── isolate02/
│       ├── report.tsv
│       └── ... (other QUAST outputs)
│
├── busco/
│   ├── isolate01/
│   │   └── short_summary.txt
│   └── isolate02/
│       └── short_summary.txt
│
├── coverm/
│   ├── isolate01.tsv
│   └── isolate02.tsv
│
├── mlst/
│   ├── isolate01.tsv
│   └── isolate02.tsv
│
├── amrfinderplus/
│   ├── isolate01.tsv
│   └── isolate02.tsv
│
├── virulencefinder/
│   ├── isolate01/
│   │   └── results_tab.tsv
│   └── isolate02/
│       └── results_tab.tsv
│
├── logs/
│   ├── fastqc/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── trimmomatic/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── kraken2/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── unicycler/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── quast/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── busco/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── coverm/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── mlst/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   ├── amrfinderplus/
│   │   ├── isolate01.log
│   │   └── isolate02.log
│   └── virulencefinder/
│       ├── isolate01.log
│       └── isolate02.log
│
└── WGS_summary.tsv
```

---

## 18. Re-running and Partial Execution

Snakemake automatically tracks which outputs are missing or outdated. You can use this to resume interrupted runs or re-execute specific rules.

### Run a Specific Rule

To run a specific rule, target its output file:

```bash
snakemake --use-conda -j 8 --configfile config.yaml \
    results/quast/isolate01/report.tsv
```

### Rebuild a Specific Rule

To force re-execution of a rule (even if outputs exist):

```bash
snakemake --use-conda -j 8 --configfile config.yaml \
    --forcerun quast
```

### Dry Run Before Changes

Always use `--dry-run` (or `-n`) to preview which jobs will be executed:

```bash
snakemake -n --use-conda -j 8 --configfile config.yaml --forcerun quast
```

---

## 19. Important Notes and Troubleshooting

### Environment Activation

This workflow requires pre-existing Conda environments. Ensure all `envs/*.yaml` files are present in the repository and that conda is available in your PATH:

```bash
conda --version
```

### Database Setup

Before running the pipeline, verify that all required databases exist and are correctly specified in `config.yaml`:

```bash
# Test Kraken2 database
ls -lah /path/to/kraken2/database/

# Test VirulenceFinder database
ls -lah /path/to/virulencefinder_db/
```

### Resource Planning

The pipeline's resource usage depends on:
- **Sample count** — More samples = longer overall runtime, but Snakemake can parallelize
- **Read depth** — Higher coverage = longer trimming, assembly, and mapping steps
- **Kraken2 database size** — Larger databases take longer but are more sensitive
- **Thread allocation** — Configure `threads.*` in config.yaml based on your CPU cores

For example, with 8 cores available:
- FastQC: 4 threads
- Trimmomatic: 8 threads
- Kraken2: 8 threads
- Unicycler: 8 threads
- Others: 2–8 threads

Snakemake will serialize jobs if total requested threads exceed available cores.

### Common Issues

**Issue:** `FileNotFoundError: samples.tsv`
- **Solution:** Verify the path to `samples.tsv` in `config.yaml` and that the file exists.

**Issue:** `ModuleNotFoundError: No module named 'pandas'`
- **Solution:** The summary rule uses `scripts/make_summary.py` which requires pandas. Ensure `envs/summary.yaml` includes pandas, or run within the Conda environment: `snakemake --use-conda ...`

**Issue:** Database not found errors
- **Solution:** Verify database paths in `config.yaml` are absolute paths and that files are readable: `ls -la /path/to/database/`

**Issue:** Insufficient disk space
- **Solution:** Trimmed FASTQ files and intermediate assemblies can consume significant space. Ensure at least 50 GB free for a typical 10-sample WGS analysis.

---

## 20. Recommended Validation Before Production

Run these checks before processing a large dataset:

* [ ] Verify all samples listed in `samples.tsv` match input files
* [ ] Confirm every sample has both R1 and R2 FASTQ files
* [ ] Test that all file paths in `samples.tsv` are readable: `ls -l <path>`
* [ ] Confirm Kraken2 database exists and is accessible
* [ ] Confirm VirulenceFinder database exists and is accessible
* [ ] Confirm BUSCO lineage dataset is installed: `busco --list-datasets | grep bacteria`
* [ ] Confirm all Conda environments work: `snakemake --use-conda --list-resources`
* [ ] Perform a **dry run**: `snakemake -n --use-conda -j 8 --configfile config.yaml`
* [ ] Test with **one representative sample** (subset of reads if desired)
* [ ] Inspect all generated logs before processing full dataset

---

## 21. Quick Command Reference

```bash
# Perform a dry run (recommended before first execution)
snakemake -n --use-conda -j 8 --configfile config.yaml

# Run the full pipeline
snakemake --use-conda -j 8 --configfile config.yaml

# Verbose output (prints executed shell commands)
snakemake -p --use-conda -j 8 --configfile config.yaml

# Generate DAG visualization
snakemake --dag --use-conda --configfile config.yaml | dot -Tpdf > dag.pdf

# Run a specific rule
snakemake --use-conda -j 8 --configfile config.yaml results/quast/sample01/report.tsv

# Force re-run of a specific rule
snakemake --use-conda -j 8 --configfile config.yaml --forcerun quast

# Generate summary statistics
snakemake --use-conda -j 8 --configfile config.yaml results/WGS_summary.tsv
```

The final deliverable is:

```
results/WGS_summary.tsv
```

alongside per-sample output directories under:
- `results/fastqc/`
- `results/trimmed/`
- `results/kraken2/`
- `results/assembly/`
- `results/quast/`
- `results/busco/`
- `results/coverm/`
- `results/mlst/`
- `results/amrfinderplus/`
- `results/virulencefinder/`

---

## References

<div id="refs" class="references" markdown="1">

<div id="ref-Bankevich2012-of" markdown="1">

Bankevich, Anton, Sergey Nurk, Dmitry Antipov, Alexey A Gurevich,
Mikhail Dvorkin, Alexander S Kulikov, Valery M Lesin, et al. 2012.
"SPAdes: A New Genome Assembly Algorithm and Its Applications to
Single-Cell Sequencing." *J. Comput. Biol.* 19 (5): 455–77.

</div>

<div id="ref-Bolger2014-qe" markdown="1">

Bolger, Anthony M, Marc Lohse, and Bjoern Usadel. 2014. "Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data." *Bioinformatics* 30 (15):
2114–20.

</div>

<div id="ref-Gurevich2013-cw" markdown="1">

Gurevich, Alexey, Vladislav Saveliev, Nikolay Vyahhi, and Glenn Tesler.
2013. "QUAST: Quality Assessment Tool for Genome Assemblies."
*Bioinformatics* 29 (8): 1072–5.

</div>

<div id="ref-Jolley2018-jn" markdown="1">

Jolley, Keith A, James E Bray, and Martin C J Maiden. 2018. "Open-Access
Bacterial Population Genomics: BIGSdb Software, the PubMLST.org Website
and Their Applications." *Wellcome Open Res* 3 (September): 124.

</div>

<div id="ref-Koster2012-cf" markdown="1">

Köster, Johannes, and Sven Rahmann. 2012. "Snakemake–a Scalable
Bioinformatics Workflow Engine." *Bioinformatics* 28 (19): 2520–2.

</div>

<div id="ref-Kurtzer2017-se" markdown="1">

Kurtzer, Gregory M, Vanessa Sochat, and Michael W Bauer. 2017.
"Singularity: Scientific Containers for Mobility of Compute." *PLoS One*
12 (5): e0177459.

</div>

<div id="ref-Ondov2016-gn" markdown="1">

Ondov, Brian D, Todd J Treangen, Páll Melsted, Adam B Mallonee, Nicholas
H Bergman, Sergey Koren, and Adam M Phillippy. 2016. "Mash: Fast Genome
and Metagenome Distance Estimation Using MinHash." *Genome Biol.* 17
(1): 132.

</div>

<div id="ref-Seemann2018-yj" markdown="1">

Seemann, Torsten. 2018. "mlst: Scan Contig Files Against PubMLST Typing Schemes."
<https://github.com/tseemann/mlst/>.

</div>

<div id="ref-Wood2014-we" markdown="1">

Wood, Derrick E, and Steven L Salzberg. 2014. "Kraken: Ultrafast
Metagenomic Sequence Classification Using Exact Alignments." *Genome
Biol.* 15 (3): R46.

</div>

<div id="ref-Yoshida2016-uu" markdown="1">

Yoshida, Catherine E, Peter Kruczkiewicz, Chad R Laing, Erika J Lingohr,
Victor P J Gannon, John H E Nash, and Eduardo N Taboada. 2016. "The
Salmonella in Silico Typing Resource (SISTR): An Open Web-Accessible
Tool for Rapidly Typing and Subtyping Draft Salmonella Genome
Assemblies." *PLoS One* 11 (1): e0147101.

</div>

</div>
