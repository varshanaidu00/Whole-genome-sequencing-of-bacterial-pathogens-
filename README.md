WGS of Bacterial Pathogens for the detection of AMR and virulence genes
=================================================================================================================

## Background
----------
This repository provides a reproducible, containerized Snakemake pipeline to analyse paired‑end Illumina whole‑genome sequencing (WGS) FASTQ data from bacterial isolates. The workflow performs taxonomic classification, de novo assembly, assembly quality assessment, and detection of clinically relevant antimicrobial resistance (AMR) genes and known virulence determinants, where applicable. Outputs include per‑sample Kraken2 reports, Unicycler assemblies with QUAST/BUSCO quality metrics, AMRFinderPlus AMR annotations, VirulenceFinder results, and a consolidated pipeline summary suitable for PHL or research. 

Instructions on how to create and execute a basic workflow in Snakemake can be found in this wiki <a href="https://github.com/varshanaidu00/Whole-genome-sequencing-of-bacterial-pathogens-/wiki/Instructions-for-writing-a-Snakemake-workflow-with-Conda" target="_blank">Instructions for writing a Snakemake workflow with Conda</a>.


The pipeline uses:

- `FastQC` - to perform read QC statistics on reads.
- `Trimmomatic` (Bolger, Lohse, and Usadel 2014) — adapter removal and quality trimming if needed.
- `Kraken2` (Wood and Salzberg 2014) — k‑mer based taxonomic classification (Standard-16 DB; Standard with DB capped at 16 GB).
- `Mash` (Ondov et al. 2016) — k‑mer based genome distance and genome‑size/coverage estimation.
- `Unicycler` (`SPAdes` v4.30) — hybrid/long‑read-aware assembly (for illumina and ONT workflows).
- `QUAST` (Gurevich et al. 2013) — assembly quality assessment and metrics.
- `BUSCO` — conserved single‑copy ortholog assessment for assembly completeness.
- `mlst` — sequence typing using PubMLST schemes.
- `NCBI‑AMRFinderPlus` — antimicrobial resistance gene detection.
- `VirulenceFinder` — virulence gene screening (optional).

Each major tool will run in its own pre-existing Conda environment where specified by the workflow.

---

## 1. Requirements

Before running the pipeline, ensure that the following are available:

* [ ] Snakemake
* [ ] Conda
* [ ] All required Conda environments
* [ ] Kraken2 database
* [ ] VirulenceFinder database
* [ ] AMRFinderPlus database
* [ ] Paired-end FASTQ files for every sample
* [ ] `Snakefile`
* [ ] `config.yaml`

The pipeline assumes that the required software environments and databases already exist.

---

## 2. Input FASTQ Files

Paired-end reads must be stored in the directory specified by `read_dir` in `config.yaml`.

Each sample must have two FASTQ files following this naming convention:

```text
SAMPLE_1.fastq
SAMPLE_2.fastq
```

For example:

```text
reads/
├── isolate01_1.fastq
├── isolate01_2.fastq
├── isolate02_1.fastq
└── isolate02_2.fastq
```

The sample names listed in `config.yaml` must correspond exactly to the FASTQ filename prefixes.

---

## 3. Configuration

The pipeline reads configuration values from `config.yaml`.

A compatible configuration should contain the following fields:

```yaml
samples:
  - isolate01
  - isolate02

read_dir: "reads"
out_dir: "results"

kraken2_db: "/path/to/kraken2/database"
virulencefinder_db: "/path/to/virulencefinder/database"

threads: 8
busco_lineage: "bacteria_odb10"
```

### Configuration parameters

| Parameter            | Description                                 |
| -------------------- | ------------------------------------------- |
| `samples`            | List of sample names to process             |
| `read_dir`           | Directory containing paired-end FASTQ files |
| `out_dir`            | Main output directory                       |
| `kraken2_db`         | Path to the Kraken2 database                |
| `virulencefinder_db` | Path to the VirulenceFinder database        |
| `threads`            | Number of CPU threads allocated to tools    |
| `busco_lineage`      | BUSCO lineage dataset                       |

If `threads` is not specified, the workflow defaults to **8 threads**.

If `busco_lineage` is not specified, the workflow defaults to:

```text
bacteria_odb10
```

---

## 4. Running the Pipeline

Run the complete workflow with:

```bash
snakemake --use-conda -j 8 --configfile config.yaml
```

Here:

* `--use-conda` enables Conda environments defined by the workflow.
* `-j 8` allows up to 8 jobs to run concurrently.
* `--configfile config.yaml` specifies the configuration file.

Before executing the complete workflow, it is recommended to perform a dry run:

```bash
snakemake -n --use-conda -j 8 --configfile config.yaml
```

For more detailed execution information:

```bash
snakemake -p --use-conda -j 8 --configfile config.yaml
```

---

## 5. Pipeline Workflow

The workflow follows this general structure:

```text
Paired-end FASTQ reads
        │
        ├── Kraken2 ───────────────► Taxonomic classification
        │
        └── Unicycler ─────────────► Genome assembly
                                      │
                                      ├── QUAST ─────────► Assembly quality
                                      │
                                      ├── BUSCO ─────────► Genome completeness
                                      │
                                      ├── AMRFinderPlus ─► AMR genes
                                      │
                                      └── VirulenceFinder ► Virulence genes
                                               
All results
    │
    └── Summary rule ──────────────► pipeline_summary.tsv
```

---

## 6. Step 1 — Kraken2

Kraken2 performs taxonomic classification of the paired-end reads.

The rule uses:

* Paired-end FASTQ files
* The configured Kraken2 database
* The configured number of threads

The Kraken2 classification output is written to:

```text
results/kraken2/<sample>/report.txt
```

The workflow suppresses the normal classification output by sending it to `/dev/null` and retains the Kraken2 report.

A log file is created at:

```text
results/logs/kraken2/<sample>.log
```

---

## 7. Step 2 — Unicycler

Unicycler performs de novo assembly of the paired-end reads.

For each sample, the workflow produces:

```text
results/unicycler/<sample>/assembly.fasta
results/unicycler/<sample>/assembly.gfa
```

The FASTA assembly is subsequently used as input for QUAST, BUSCO, AMRFinderPlus, and VirulenceFinder.

The Unicycler log is stored at:

```text
results/logs/unicycler/<sample>.log
```

---

## 8. Step 3 — QUAST

QUAST evaluates the quality and characteristics of each genome assembly.

The input is:

```text
results/unicycler/<sample>/assembly.fasta
```

The main report required by the workflow is:

```text
results/quast/<sample>/report.tsv
```

The workflow uses QUAST results to extract:

* Total assembly length
* Number of contigs
* N50
* GC content

These metrics are subsequently incorporated into the final summary.

---

## 9. Step 4 — BUSCO

BUSCO evaluates genome completeness using the configured bacterial lineage dataset.

The default lineage is:

```text
bacteria_odb10
```

The expected summary file is:

```text
results/busco/<sample>/short_summary.txt
```

The pipeline extracts the BUSCO complete percentage from this file for the final summary.

> **Important:** The current Snakefile contains hard-coded paths to the BUSCO executable and Python interpreter under `/home/varsha.naidu@rcpaqap.local/...`. These paths are specific to the original environment and should be replaced with portable environment-based commands if the workflow is moved to another system.

---

## 10. Step 5 — NCBI-AMRFinderPlus

AMRFinderPlus identifies antimicrobial resistance genes from the assembled genome.

Input:

```text
results/unicycler/<sample>/assembly.fasta
```

Output:

```text
results/amrfinder/<sample>.tsv
```

The final summary counts the number of AMR-gene result rows reported by AMRFinderPlus.

> **Important:** The current Snakefile uses a hard-coded AMRFinderPlus executable and database path. These paths should be checked and updated when running the workflow on a different machine or cluster.

---

## 11. Step 6 — VirulenceFinder

VirulenceFinder searches the assembled genome for virulence-associated genes.

Input:

```text
results/unicycler/<sample>/assembly.fasta
```

Output directory:

```text
results/virulencefinder/<sample>/
```

The final workflow target is:

```text
results/virulencefinder/<sample>/data.json
```

The VirulenceFinder database is supplied through the `virulencefinder_db` configuration parameter.

> **Important:** The current Snakefile also contains hard-coded paths to the VirulenceFinder Python executable and `blastn`. These should be made environment-independent before distributing the workflow.

---

## 12. Final Summary

The `summary` rule combines selected results from Kraken2, QUAST, BUSCO, and AMRFinderPlus into:

```text
results/pipeline_summary.tsv
```

The summary contains one row per sample and the following columns:

| Column              | Description                                               |
| ------------------- | --------------------------------------------------------- |
| `Sample`            | Sample identifier                                         |
| `Top_species(%)`    | First species-level Kraken2 classification and percentage |
| `Unclassified(%)`   | Percentage of reads classified as unclassified            |
| `Assembly_size(bp)` | Total assembly length                                     |
| `Num_contigs`       | Number of assembly contigs                                |
| `N50(bp)`           | Assembly N50                                              |
| `GC(%)`             | Assembly GC percentage                                    |
| `BUSCO_complete(%)` | BUSCO completeness percentage                             |
| `AMR_genes_found`   | Number of AMR genes detected by AMRFinderPlus             |

If a metric cannot be parsed from its expected output file, the pipeline records:

```text
NA
```

---

## 13. Expected Output Structure

For two example samples, the output directory should resemble:

```text
results/
├── kraken2/
│   ├── isolate01/
│   │   └── report.txt
│   └── isolate02/
│       └── report.txt
│
├── unicycler/
│   ├── isolate01/
│   │   ├── assembly.fasta
│   │   └── assembly.gfa
│   └── isolate02/
│       ├── assembly.fasta
│       └── assembly.gfa
│
├── quast/
│   ├── isolate01/
│   │   └── report.tsv
│   └── isolate02/
│       └── report.tsv
│
├── busco/
│   ├── isolate01/
│   │   └── short_summary.txt
│   └── isolate02/
│       └── short_summary.txt
│
├── amrfinder/
│   ├── isolate01.tsv
│   └── isolate02.tsv
│
├── virulencefinder/
│   ├── isolate01/
│   │   └── data.json
│   └── isolate02/
│       └── data.json
│
├── logs/
│   ├── kraken2/
│   ├── unicycler/
│   ├── quast/
│   ├── busco/
│   ├── amrfinder/
│   └── virulencefinder/
│
└── pipeline_summary.tsv
```

---

## 14. Re-running and Partial Execution

Snakemake automatically determines which outputs are missing or outdated.

To run a specific rule, for example QUAST:

```bash
snakemake --use-conda -j 8 --configfile config.yaml \
    results/quast/isolate01/report.tsv
```

To rebuild a specific output:

```bash
snakemake --use-conda -j 8 --configfile config.yaml \
    --forcerun quast
```

Use `--dry-run` before executing changes to verify the jobs Snakemake intends to run.

---

## 15. Important Portability Considerations

The Snakefile currently mixes Conda-managed environments with hard-coded absolute executable paths.

The following paths are system-specific:

```text
/home/varsha.naidu@rcpaqap.local/anaconda3/envs/busco/...
/home/varsha.naidu@rcpaqap.local/anaconda3/envs/ncbi-amrfinderplus/...
/home/varsha.naidu@rcpaqap.local/anaconda3/envs/virulencefinder/...
```

For a portable Snakemake workflow, these should ideally be replaced with commands available through the corresponding Conda environments, for example:

```bash
busco
amrfinder
python -m virulencefinder
blastn
```

This allows the workflow to run on different machines without modifying user-specific filesystem paths.

The AMRFinderPlus database path is also currently hard-coded and should preferably be moved into `config.yaml`.

---

## 16. Recommended Validation Before a Production Run

Run the following checks before processing a large dataset:

* [ ] Confirm all sample names in `config.yaml`.
* [ ] Confirm every sample has both `_1.fastq` and `_2.fastq`.
* [ ] Confirm the Kraken2 database exists and is accessible.
* [ ] Confirm the VirulenceFinder database exists and is accessible.
* [ ] Confirm the AMRFinderPlus database is available.
* [ ] Confirm the BUSCO lineage dataset is installed.
* [ ] Confirm all Conda environments are functional.
* [ ] Run a Snakemake dry run.
* [ ] Test the pipeline with one representative sample.
* [ ] Inspect the generated logs before processing the full dataset.

---

## 17. Minimal Command Reference

```bash
# Dry run
snakemake -n --use-conda -j 8 --configfile config.yaml

# Run pipeline
snakemake --use-conda -j 8 --configfile config.yaml

# Verbose command execution
snakemake -p --use-conda -j 8 --configfile config.yaml

# Generate a DAG
snakemake --dag --use-conda --configfile config.yaml | dot -Tpdf > pipeline_dag.pdf
```

The final deliverable is:

```text
<out_dir>/pipeline_summary.tsv
```
alongside the per-sample Kraken2, assembly, QUAST, BUSCO, AMRFinderPlus, and VirulenceFinder results.



References 
----------

<div id="refs" class="references" markdown="1">

<div id="ref-Bankevich2012-of" markdown="1">

Bankevich, Anton, Sergey Nurk, Dmitry Antipov, Alexey A Gurevich,
Mikhail Dvorkin, Alexander S Kulikov, Valery M Lesin, et al. 2012.
“SPAdes: A New Genome Assembly Algorithm and Its Applications to
Single-Cell Sequencing.” *J. Comput. Biol.* 19 (5): 455–77.

</div>

<div id="ref-Bolger2014-qe" markdown="1">

Bolger, Anthony M, Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics* 30 (15):
2114–20.

</div>

<div id="ref-Gurevich2013-cw" markdown="1">

Gurevich, Alexey, Vladislav Saveliev, Nikolay Vyahhi, and Glenn Tesler.
2013. “QUAST: Quality Assessment Tool for Genome Assemblies.”
*Bioinformatics* 29 (8): 1072–5.

</div>

<div id="ref-Jolley2018-jn" markdown="1">

Jolley, Keith A, James E Bray, and Martin C J Maiden. 2018. “Open-Access
Bacterial Population Genomics: BIGSdb Software, the PubMLST.org Website
and Their Applications.” *Wellcome Open Res* 3 (September): 124.

</div>

<div id="ref-Koster2012-cf" markdown="1">

Köster, Johannes, and Sven Rahmann. 2012. “Snakemake–a Scalable
Bioinformatics Workflow Engine.” *Bioinformatics* 28 (19): 2520–2.

</div>

<div id="ref-Kurtzer2017-se" markdown="1">

Kurtzer, Gregory M, Vanessa Sochat, and Michael W Bauer. 2017.
“Singularity: Scientific Containers for Mobility of Compute.” *PLoS One*
12 (5): e0177459.

</div>

<div id="ref-Li2018-ow" markdown="1">

Li, Heng. 2018. “Seqtk: Toolkit for Processing Sequences in FASTA/Q
Formats.” <https://github.com/lh3/seqtk>.

</div>

<div id="ref-Ncbi_undated-gz" markdown="1">

NCBI. n.d. “NCBI AMR Reference Gene Database.”
<https://www.ncbi.nlm.nih.gov/pathogens/isolates#/refgene/>.

</div>

<div id="ref-Ondov2016-gn" markdown="1">

Ondov, Brian D, Todd J Treangen, Páll Melsted, Adam B Mallonee, Nicholas
H Bergman, Sergey Koren, and Adam M Phillippy. 2016. “Mash: Fast Genome
and Metagenome Distance Estimation Using MinHash.” *Genome Biol.* 17
(1): 132.


</div>

<div id="ref-Seemann2018-yj" markdown="1">

———. 2018b. “Mlst: Scan Contig Files Against PubMLST Typing Schemes.”
<https://github.com/tseemann/mlst/>.


</div>

<div id="ref-Wood2014-we" markdown="1">

Wood, Derrick E, and Steven L Salzberg. 2014. “Kraken: Ultrafast
Metagenomic Sequence Classification Using Exact Alignments.” *Genome
Biol.* 15 (3): R46.

</div>

<div id="ref-Yoshida2016-uu" markdown="1">

Yoshida, Catherine E, Peter Kruczkiewicz, Chad R Laing, Erika J Lingohr,
Victor P J Gannon, John H E Nash, and Eduardo N Taboada. 2016. “The
Salmonella in Silico Typing Resource (SISTR): An Open Web-Accessible
Tool for Rapidly Typing and Subtyping Draft Salmonella Genome
Assemblies.” *PLoS One* 11 (1): e0147101.

</div>

</div>


