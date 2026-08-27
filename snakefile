configfile: "config.yaml"

import pandas as pd


# ============================================================
# SAMPLE INFORMATION
# ============================================================

samples = pd.read_csv(config["samples"], sep="\t")

SAMPLES = samples["sample"].tolist()

R1 = dict(zip(samples["sample"], samples["R1"]))
R2 = dict(zip(samples["sample"], samples["R2"]))


RESULTS = config["results"]
LOGS = config["logs"]


# ============================================================
# FINAL TARGET
# ============================================================

rule all:

    input:

        # ----------------------------------------------------
        # FastQC
        # ----------------------------------------------------

        expand(
            RESULTS + "/fastqc/{sample}_R1_paired_fastqc.html",
            sample=SAMPLES
        ),

        expand(
            RESULTS + "/fastqc/{sample}_R2_paired_fastqc.html",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # Trimmomatic
        # ----------------------------------------------------

        expand(
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",
            sample=SAMPLES
        ),

        expand(
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # Kraken2
        # ----------------------------------------------------

        expand(
            RESULTS + "/kraken2/{sample}.kraken2.txt",
            sample=SAMPLES
        ),

        expand(
            RESULTS + "/kraken2/{sample}.report.txt",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # Assembly
        # ----------------------------------------------------

        expand(
            RESULTS + "/assembly/{sample}/assembly.fasta",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # QUAST
        # ----------------------------------------------------

        expand(
            RESULTS + "/quast/{sample}/report.tsv",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # BUSCO
        # ----------------------------------------------------

        expand(
            RESULTS + "/busco/{sample}/short_summary.txt",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # CoverM
        # ----------------------------------------------------

        expand(
            RESULTS + "/coverm/{sample}.tsv",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # MLST
        # ----------------------------------------------------

        expand(
            RESULTS + "/mlst/{sample}.tsv",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # AMRFinderPlus
        # ----------------------------------------------------

        expand(
            RESULTS + "/amrfinderplus/{sample}.tsv",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # VirulenceFinder
        # ----------------------------------------------------

        expand(
            RESULTS + "/virulencefinder/{sample}/results_tab.tsv",
            sample=SAMPLES
        ),

        # ----------------------------------------------------
        # Final summary
        # ----------------------------------------------------

        RESULTS + "/WGS_summary.tsv"


# ============================================================
# TRIMMOMATIC
# ============================================================

rule trimmomatic:

    input:

        r1=lambda wc: R1[wc.sample],

        r2=lambda wc: R2[wc.sample]

    output:

        r1_paired=
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",

        r1_unpaired=
            RESULTS + "/trimmed/{sample}_R1_unpaired.fastq.gz",

        r2_paired=
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz",

        r2_unpaired=
            RESULTS + "/trimmed/{sample}_R2_unpaired.fastq.gz"

    threads:
        config["threads"]["trimmomatic"]

    params:

        adapters=config["adapters"],

        leading=config["trimmomatic"]["leading"],

        trailing=config["trimmomatic"]["trailing"],

        slidingwindow=config["trimmomatic"]["slidingwindow"],

        minlen=config["trimmomatic"]["minlen"]

    log:
        LOGS + "/trimmomatic/{sample}.log"

    conda:
        "envs/trimmomatic.yaml"

    shell:

        """
        mkdir -p {RESULTS}/trimmed
        mkdir -p {LOGS}/trimmomatic

        trimmomatic PE \
            -threads {threads} \
            {input.r1} \
            {input.r2} \
            {output.r1_paired} \
            {output.r1_unpaired} \
            {output.r2_paired} \
            {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.slidingwindow} \
            MINLEN:{params.minlen} \
            > {log} 2>&1
        """


# ============================================================
# FASTQC ON TRIMMED READS
# ============================================================

rule fastqc:

    input:

        r1=
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",

        r2=
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz"

    output:

        r1_html=
            RESULTS + "/fastqc/{sample}_R1_paired_fastqc.html",

        r2_html=
            RESULTS + "/fastqc/{sample}_R2_paired_fastqc.html",

        r1_zip=
            RESULTS + "/fastqc/{sample}_R1_paired_fastqc.zip",

        r2_zip=
            RESULTS + "/fastqc/{sample}_R2_paired_fastqc.zip",

        r1_data=
            RESULTS + "/fastqc/{sample}_R1_paired_fastqc/fastqc_data.txt",

        r2_data=
            RESULTS + "/fastqc/{sample}_R2_paired_fastqc/fastqc_data.txt"

    threads:
        config["threads"]["fastqc"]

    log:
        LOGS + "/fastqc/{sample}.log"

    conda:
        "envs/fastqc.yaml"

    shell:

        """
        mkdir -p {RESULTS}/fastqc
        mkdir -p {LOGS}/fastqc

        fastqc \
            {input.r1} \
            {input.r2} \
            --outdir {RESULTS}/fastqc \
            --threads {threads} \
            --extract \
            > {log} 2>&1
        """


# ============================================================
# KRAKEN2
# ============================================================

rule kraken2:

    input:

        r1=
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",

        r2=
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz"

    output:

        kraken=
            RESULTS + "/kraken2/{sample}.kraken2.txt",

        report=
            RESULTS + "/kraken2/{sample}.report.txt"

    threads:
        config["threads"]["kraken2"]

    params:

        database=config["kraken2_db"],

        confidence=config["kraken2"]["confidence"]

    log:
        LOGS + "/kraken2/{sample}.log"

    conda:
        "envs/kraken2.yaml"

    shell:

        """
        mkdir -p {RESULTS}/kraken2
        mkdir -p {LOGS}/kraken2

        kraken2 \
            --db {params.database} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --confidence {params.confidence} \
            --output {output.kraken} \
            --report {output.report} \
            {input.r1} \
            {input.r2} \
            > {log} 2>&1
        """


# ============================================================
# UNICYCLER DE NOVO ASSEMBLY
# ============================================================

rule unicycler:

    input:

        r1=
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",

        r2=
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz"

    output:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    threads:
        config["threads"]["unicycler"]

    params:

        outdir=lambda wc:
            f"{RESULTS}/assembly/{wc.sample}"

    log:
        LOGS + "/unicycler/{sample}.log"

    conda:
        "envs/unicycler.yaml"

    shell:

        """
        mkdir -p {params.outdir}
        mkdir -p {LOGS}/unicycler

        unicycler \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            --mode normal \
            --threads {threads} \
            > {log} 2>&1
        """


# ============================================================
# QUAST
# ============================================================

rule quast:

    input:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        report=
            RESULTS + "/quast/{sample}/report.tsv"

    threads:
        config["threads"]["quast"]

    params:

        outdir=lambda wc:
            f"{RESULTS}/quast/{wc.sample}",

        reference=config["reference"]

    log:
        LOGS + "/quast/{sample}.log"

    conda:
        "envs/quast.yaml"

    shell:

        r"""
        mkdir -p {params.outdir}
        mkdir -p {LOGS}/quast

        if [ "{params.reference}" != "None" ] && [ -n "{params.reference}" ]; then

            quast.py \
                {input.assembly} \
                -r {params.reference} \
                -t {threads} \
                -o {params.outdir} \
                > {log} 2>&1

        else

            quast.py \
                {input.assembly} \
                -t {threads} \
                -o {params.outdir} \
                > {log} 2>&1

        fi
        """


# ============================================================
# BUSCO
# ============================================================

rule busco:

    input:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        summary=
            RESULTS + "/busco/{sample}/short_summary.txt"

    threads:
        config["threads"]["busco"]

    params:

        lineage=config["busco_lineage"],

        outdir=lambda wc:
            f"{RESULTS}/busco/{wc.sample}"

    log:
        LOGS + "/busco/{sample}.log"

    conda:
        "envs/busco.yaml"

    shell:

        """
        mkdir -p {params.outdir}
        mkdir -p {LOGS}/busco

        busco \
            -i {input.assembly} \
            -m genome \
            -l {params.lineage} \
            -o {wildcards.sample} \
            --out_path {params.outdir} \
            -c {threads} \
            --force \
            > {log} 2>&1

        cp {params.outdir}/{wildcards.sample}/short_summary*.txt \
            {output.summary}
        """


# ============================================================
# COVERM
# ============================================================

rule coverm:

    input:

        r1=
            RESULTS + "/trimmed/{sample}_R1_paired.fastq.gz",

        r2=
            RESULTS + "/trimmed/{sample}_R2_paired.fastq.gz",

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        RESULTS + "/coverm/{sample}.tsv"

    threads:
        config["threads"]["coverm"]

    params:

        identity=config["coverm"]["min_identity"],

        aligned=config["coverm"]["min_aligned_percent"]

    log:
        LOGS + "/coverm/{sample}.log"

    conda:
        "envs/coverm.yaml"

    shell:

        """
        mkdir -p {RESULTS}/coverm
        mkdir -p {LOGS}/coverm

        coverm genome \
            --coupled {input.r1} {input.r2} \
            --genome-fasta-files {input.assembly} \
            --threads {threads} \
            --min-read-percent-identity {params.identity} \
            --min-read-aligned-percent {params.aligned} \
            --methods mean covered_fraction \
            --output-file {output} \
            > {log} 2>&1
        """


# ============================================================
# MLST
# ============================================================

rule mlst:

    input:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        RESULTS + "/mlst/{sample}.tsv"

    log:
        LOGS + "/mlst/{sample}.log"

    conda:
        "envs/mlst.yaml"

    shell:

        """
        mkdir -p {RESULTS}/mlst
        mkdir -p {LOGS}/mlst

        mlst \
            {input.assembly} \
            > {output} \
            2> {log}
        """


# ============================================================
# AMRFINDERPLUS
# ============================================================

rule amrfinderplus:

    input:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        RESULTS + "/amrfinderplus/{sample}.tsv"

    threads:
        config["threads"]["amrfinderplus"]

    log:
        LOGS + "/amrfinderplus/{sample}.log"

    conda:
        "envs/amrfinderplus.yaml"

    shell:

        """
        mkdir -p {RESULTS}/amrfinderplus
        mkdir -p {LOGS}/amrfinderplus

        amrfinder \
            -n {input.assembly} \
            -o {output} \
            --threads {threads} \
            > {log} 2>&1
        """


# ============================================================
# VIRULENCEFINDER
# ============================================================

rule virulencefinder:

    input:

        assembly=
            RESULTS + "/assembly/{sample}/assembly.fasta"

    output:

        results=
            RESULTS +
            "/virulencefinder/{sample}/results_tab.tsv"

    threads:
        config["threads"]["virulencefinder"]

    params:

        outdir=lambda wc:
            f"{RESULTS}/virulencefinder/{wc.sample}",

        database=config["virulencefinder_db"],

        mincov=config["virulencefinder"]["min_coverage"],

        minid=config["virulencefinder"]["min_identity"]

    log:
        LOGS + "/virulencefinder/{sample}.log"

    conda:
        "envs/virulencefinder.yaml"

    shell:

        """
        mkdir -p {params.outdir}
        mkdir -p {LOGS}/virulencefinder

        virulencefinder.py \
            -i {input.assembly} \
            -o {params.outdir} \
            -p {params.database} \
            -l {params.mincov} \
            -t {params.minid} \
            -x \
            > {log} 2>&1
        """


# ============================================================
# FINAL SUMMARY
# ============================================================

rule summary:

    input:

        fastqc_r1=lambda wc:
            expand(
                RESULTS +
                "/fastqc/{sample}_R1_paired_fastqc/fastqc_data.txt",
                sample=SAMPLES
            ),

        fastqc_r2=lambda wc:
            expand(
                RESULTS +
                "/fastqc/{sample}_R2_paired_fastqc/fastqc_data.txt",
                sample=SAMPLES
            ),

        kraken=lambda wc:
            expand(
                RESULTS + "/kraken2/{sample}.kraken2.txt",
                sample=SAMPLES
            ),

        kraken_report=lambda wc:
            expand(
                RESULTS + "/kraken2/{sample}.report.txt",
                sample=SAMPLES
            ),

        quast=lambda wc:
            expand(
                RESULTS + "/quast/{sample}/report.tsv",
                sample=SAMPLES
            ),

        busco=lambda wc:
            expand(
                RESULTS + "/busco/{sample}/short_summary.txt",
                sample=SAMPLES
            ),

        coverm=lambda wc:
            expand(
                RESULTS + "/coverm/{sample}.tsv",
                sample=SAMPLES
            ),

        mlst=lambda wc:
            expand(
                RESULTS + "/mlst/{sample}.tsv",
                sample=SAMPLES
            ),

        amr=lambda wc:
            expand(
                RESULTS + "/amrfinderplus/{sample}.tsv",
                sample=SAMPLES
            ),

        virulence=lambda wc:
            expand(
                RESULTS +
                "/virulencefinder/{sample}/results_tab.tsv",
                sample=SAMPLES
            )

    output:

        RESULTS + "/WGS_summary.tsv"

    conda:
        "envs/summary.yaml"

    script:
        "scripts/make_summary.py"
