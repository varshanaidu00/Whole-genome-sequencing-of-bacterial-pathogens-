import os
import re
import pandas as pd


RESULTS = snakemake.config["results"]

SAMPLES = pd.read_csv(
    snakemake.config["samples"],
    sep="\t"
)["sample"].tolist()


# ============================================================
# FASTQC
# ============================================================

def parse_fastqc(data_file):

    metrics = {
        "Total_Reads": "NA",
        "Unique_Reads": "NA",
        "Duplicate_Reads": "NA",
        "Pct_Duplicates": "NA",
        "Avg_Q_Score": "NA",
        "Pct_Q30": "NA",
        "GC_Content": "NA",
        "Sequence_Length": "NA"
    }

    if not os.path.exists(data_file):
        return metrics

    with open(data_file) as f:
        lines = f.readlines()

    for line in lines:

        parts = line.rstrip().split("\t")

        if len(parts) < 2:
            continue

        key = parts[0]
        value = parts[1]

        if key == "Total Sequences":
            metrics["Total_Reads"] = value

        elif key == "Sequence length":
            metrics["Sequence_Length"] = value

        elif key == "%GC":
            metrics["GC_Content"] = value

        elif key == "Total Deduplicated Percentage":

            try:

                dedup = float(value)

                total = float(metrics["Total_Reads"])

                unique = round(total * dedup / 100)

                duplicate = round(total - unique)

                metrics["Unique_Reads"] = unique

                metrics["Duplicate_Reads"] = duplicate

                metrics["Pct_Duplicates"] = round(
                    100 - dedup, 2
                )

            except (ValueError, TypeError):
                pass

    # --------------------------------------------------------
    # Average Q score
    # --------------------------------------------------------

    in_quality = False

    quality_sum = 0
    quality_count = 0

    for line in lines:

        if line.startswith(
            ">>Per sequence quality scores"
        ):
            in_quality = True
            continue

        if line.startswith(">>END_MODULE"):
            in_quality = False

        if in_quality:

            parts = line.rstrip().split("\t")

            if (
                len(parts) >= 2
                and re.match(r"^[0-9]", parts[0])
            ):

                try:

                    q = float(parts[0])
                    count = float(parts[1])

                    quality_sum += q * count
                    quality_count += count

                except ValueError:
                    pass

    if quality_count > 0:

        metrics["Avg_Q_Score"] = round(
            quality_sum / quality_count,
            2
        )

    # --------------------------------------------------------
    # Q30
    #
    # Percentage of base positions with mean quality >= Q30.
    # --------------------------------------------------------

    in_quality = False

    total_positions = 0
    q30_positions = 0

    for line in lines:

        if line.startswith(
            ">>Per base sequence quality"
        ):
            in_quality = True
            continue

        if line.startswith(">>END_MODULE"):
            in_quality = False

        if in_quality:

            parts = line.rstrip().split("\t")

            if (
                len(parts) >= 2
                and re.match(r"^[0-9]", parts[0])
            ):

                try:

                    quality = float(parts[1])

                    total_positions += 1

                    if quality >= 30:
                        q30_positions += 1

                except ValueError:
                    pass

    if total_positions > 0:

        metrics["Pct_Q30"] = round(
            q30_positions /
            total_positions *
            100,
            2
        )

    return metrics


# ============================================================
# KRAKEN2
# ============================================================

def parse_kraken(report):

    result = {
        "Kraken2_Classified_Pct": "NA",
        "Kraken2_Unclassified_Pct": "NA",
        "Kraken2_Top_Taxon": "NA"
    }

    if not os.path.exists(report):
        return result

    taxa = []

    with open(report) as f:

        for line in f:

            parts = line.rstrip().split("\t")

            if len(parts) < 6:
                continue

            try:

                pct = float(parts[0])

            except ValueError:
                continue

            taxon = parts[5].strip()

            rank = parts[3]

            if rank == "U":
                result["Kraken2_Unclassified_Pct"] = pct

            elif rank == "S":

                taxa.append(
                    (pct, taxon)
                )

    if taxa:

        taxa.sort(reverse=True)

        result["Kraken2_Top_Taxon"] = taxa[0][1]

        result["Kraken2_Classified_Pct"] = round(
            100 -
            float(result["Kraken2_Unclassified_Pct"]),
            2
        )

    return result


# ============================================================
# QUAST
# ============================================================

def parse_quast(report):

    result = {
        "Genome_Fraction": "NA",
        "Assembly_N50_kbp": "NA",
        "Assembly_Size_bp": "NA",
        "Number_of_Contigs": "NA",
        "Percent_GC": "NA"
    }

    if not os.path.exists(report):
        return result

    df = pd.read_csv(
        report,
        sep="\t"
    )

    if df.empty:
        return result

    row = df.iloc[0]

    columns = {
        "Genome fraction (%)": "Genome_Fraction",
        "N50": "Assembly_N50_kbp",
        "Total length": "Assembly_Size_bp",
        "# contigs": "Number_of_Contigs",
        "GC (%)": "Percent_GC"
    }

    for quast_name, summary_name in columns.items():

        if quast_name in df.columns:

            value = row[quast_name]

            if pd.notna(value):
                result[summary_name] = value

    # QUAST N50 is reported in bp; convert to kb.
    try:

        result["Assembly_N50_kbp"] = round(
            float(result["Assembly_N50_kbp"]) / 1000,
            2
        )

    except (ValueError, TypeError):
        pass

    return result


# ============================================================
# BUSCO
# ============================================================

def parse_busco(summary):

    result = {
        "BUSCO_Completeness": "NA"
    }

    if not os.path.exists(summary):
        return result

    with open(summary) as f:

        text = f.read()

    match = re.search(
        r"C:(\d+(?:\.\d+)?)%",
        text
    )

    if match:

        result["BUSCO_Completeness"] = (
            match.group(1)
        )

    return result


# ============================================================
# COVERM
# ============================================================

def parse_coverm(file):

    result = {
        "Depth_of_Coverage": "NA"
    }

    if not os.path.exists(file):
        return result

    try:

        df = pd.read_csv(
            file,
            sep="\t"
        )

        if df.empty:
            return result

        # CoverM's Mean column is normally the mean coverage.
        for column in df.columns:

            if column.lower() == "mean":

                result["Depth_of_Coverage"] = round(
                    float(df.iloc[0][column]),
                    2
                )

                break

    except Exception:
        pass

    return result


# ============================================================
# MLST
# ============================================================

def parse_mlst(file):

    result = {
        "ST": "NA"
    }

    if not os.path.exists(file):
        return result

    try:

        with open(file) as f:

            line = f.readline().strip()

        parts = line.split("\t")

        if len(parts) >= 2:

            result["ST"] = parts[1]

    except Exception:
        pass

    return result


# ============================================================
# AMRFINDERPLUS
# ============================================================

def parse_amr(file):

    result = {
        "AMR_Genes": "None"
    }

    if not os.path.exists(file):
        return result

    genes = []

    try:

        df = pd.read_csv(
            file,
            sep="\t"
        )

        if "Element symbol" in df.columns:

            genes = sorted(
                set(
                    df["Element symbol"]
                    .dropna()
                    .astype(str)
                )
            )

        elif "Gene symbol" in df.columns:

            genes = sorted(
                set(
                    df["Gene symbol"]
                    .dropna()
                    .astype(str)
                )
            )

    except Exception:
        pass

    if genes:

        result["AMR_Genes"] = ";".join(genes)

    return result


# ============================================================
# VIRULENCEFINDER
# ============================================================

def parse_virulence(file):

    result = {
        "Virulence_Genes": "None"
    }

    if not os.path.exists(file):
        return result

    genes = []

    try:

        df = pd.read_csv(
            file,
            sep="\t"
        )

        possible_columns = [
            "Virulence factor",
            "Gene",
            "Gene name",
            "VF"
        ]

        for column in possible_columns:

            if column in df.columns:

                genes = sorted(
                    set(
                        df[column]
                        .dropna()
                        .astype(str)
                    )
                )

                break

    except Exception:
        pass

    if genes:

        result["Virulence_Genes"] = ";".join(
            genes
        )

    return result


# ============================================================
# BUILD FINAL TABLE
# ============================================================

records = []


for sample in SAMPLES:

    # --------------------------------------------------------
    # FastQC R1
    # --------------------------------------------------------

    r1_fastqc = (
        f"{RESULTS}/fastqc/"
        f"{sample}_R1_paired_fastqc/"
        f"fastqc_data.txt"
    )

    r1 = parse_fastqc(r1_fastqc)

    # --------------------------------------------------------
    # FastQC R2
    #
    # Use R1 for the primary summary metrics.
    # --------------------------------------------------------

    # --------------------------------------------------------
    # Kraken2
    # --------------------------------------------------------

    kraken_report = (
        f"{RESULTS}/kraken2/"
        f"{sample}.report.txt"
    )

    kraken = parse_kraken(
        kraken_report
    )

    # --------------------------------------------------------
    # QUAST
    # --------------------------------------------------------

    quast = parse_quast(
        f"{RESULTS}/quast/"
        f"{sample}/report.tsv"
    )

    # --------------------------------------------------------
    # BUSCO
    # --------------------------------------------------------

    busco = parse_busco(
        f"{RESULTS}/busco/"
        f"{sample}/short_summary.txt"
    )

    # --------------------------------------------------------
    # CoverM
    # --------------------------------------------------------

    coverm = parse_coverm(
        f"{RESULTS}/coverm/"
        f"{sample}.tsv"
    )

    # --------------------------------------------------------
    # MLST
    # --------------------------------------------------------

    mlst = parse_mlst(
        f"{RESULTS}/mlst/"
        f"{sample}.tsv"
    )

    # --------------------------------------------------------
    # AMRFinderPlus
    # --------------------------------------------------------

    amr = parse_amr(
        f"{RESULTS}/amrfinderplus/"
        f"{sample}.tsv"
    )

    # --------------------------------------------------------
    # VirulenceFinder
    # --------------------------------------------------------

    virulence = parse_virulence(
        f"{RESULTS}/virulencefinder/"
        f"{sample}/results_tab.tsv"
    )

    # --------------------------------------------------------
    # Combine
    # --------------------------------------------------------

    record = {
        "Sample": sample,

        # FastQC
        "Total_Reads":
            r1["Total_Reads"],

        "Unique_Reads":
            r1["Unique_Reads"],

        "Duplicate_Reads":
            r1["Duplicate_Reads"],

        "Pct_Duplicates(%)":
            r1["Pct_Duplicates"],

        "Avg_Q_Score":
            r1["Avg_Q_Score"],

        "Pct_Q30(%)":
            r1["Pct_Q30"],

        "GC_Content(%)":
            r1["GC_Content"],

        "Sequence_Length":
            r1["Sequence_Length"],

        # Kraken2
        "Kraken2_Classified(%)":
            kraken["Kraken2_Classified_Pct"],

        "Kraken2_Unclassified(%)":
            kraken["Kraken2_Unclassified_Pct"],

        "Kraken2_Top_Taxon":
            kraken["Kraken2_Top_Taxon"],

        # QUAST
        "Genome_Fraction(%)":
            quast["Genome_Fraction"],

        "Assembly_N50(kbp)":
            quast["Assembly_N50_kbp"],

        "Assembly_Size(bp)":
            quast["Assembly_Size_bp"],

        "Number_of_Contigs":
            quast["Number_of_Contigs"],

        "Percent_GC(%)":
            quast["Percent_GC"],

        # BUSCO
        "BUSCO_Completeness":
            busco["BUSCO_Completeness"],

        # CoverM
        "Depth_of_Coverage":
            coverm["Depth_of_Coverage"],

        # MLST
        "ST":
            mlst["ST"],

        # AMR
        "AMR_Genes":
            amr["AMR_Genes"],

        # Virulence
        "Virulence_Genes":
            virulence["Virulence_Genes"]
    }

    records.append(record)


# ============================================================
# WRITE SUMMARY
# ============================================================

summary = pd.DataFrame(records)

summary.to_csv(
    snakemake.output[0],
    sep="\t",
    index=False
)

print(
    f"Summary written to: "
    f"{snakemake.output[0]}"
)
