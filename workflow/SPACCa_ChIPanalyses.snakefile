#########################################
## SPACCa: Snakefile for ChIP analyses ##
#########################################

import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
from itertools import combinations


# Define general variables
genome_used = (str(config["genomic_annotations"]["genome_id"])).lower()
genome_fasta = str(config["genomic_annotations"]["genome_fasta"])
blacklist = str(config["genomic_annotations"]["blacklist"])
genomeSize = config["genomic_annotations"]["effective_genomeSize"]
ignore_for_normalization = str(config["genomic_annotations"]["ignore_for_normalization"])


if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    read_extension = "--extendReads"
else:
    read_extension = "--extendReads "+str(config["fragment_length"])


if ((eval(str(config["use_macs3"])) == True)):
    macs_version = "macs3"
else:
    macs_version = "macs2"


### working directory
home_dir = os.path.join(config["workflow_configuration"]["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir


### get the unique samples names and other variables
# loading the sample table
sample_metadata = pd.read_csv(str(config["workflow_configuration"]["sample_config_table"]),  sep='\t+', engine='python')   # target_id | input_id | broad
#sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'input_id', 'broad'], axis=1, inplace=False)
sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'input_id', 'broad'], axis=1, copy=False)
TARGETNAMES = list(numpy.unique(list(sample_metadata.target_id)))
INPUTNAMES = list(numpy.unique(list(sample_metadata.input_id)))
SAMPLENAMES = list(numpy.unique(TARGETNAMES + INPUTNAMES))


# Get bam list
if not (os.path.exists(config["workflow_configuration"]["runs_directory"])):
    os.system("printf '\033[1;31m\\n!!! *runs_directory* does not exist !!!\\n\\n\033[0m'")
else:
    BAMS = next(os.walk(config["workflow_configuration"]["runs_directory"]))[2]
    RUNNAMES = numpy.unique([re.sub(rf"{config['bam_features']['bam_suffix']}$", "", i) for i in BAMS])



### MACS2 or MACS3?
if ((eval(str(config["use_macs3"])) == True)):
    PEAKCALLER = "macs3"
else:
    PEAKCALLER = "macs2"

### other generic variable
if ((eval(str(config["bam_features"]["remove_duplicates"])) == True)):
    DUP = "dedup"
else:
    DUP = "mdup"



### Optional analysis outputs
if ((eval(str(config["perform_plotFingerprint"])) == True)):
    plotFingerprint_results = expand("05_Quality_controls_and_statistics/plotFingerprint/{target}_fingerPrinting_plot.pdf", target = TARGETNAMES)
else:
    plotFingerprint_results = []

if ((eval(str(config["perform_fragmentSizeDistribution"])) == True) & (eval(str(config["bam_features"]["paired_end"])) == True)):
    fragmentSizeDistribution_results = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt"
else:
    fragmentSizeDistribution_results = []

if (len(SAMPLENAMES) > 2):
    if (len(INPUTNAMES) > 2):
        PCA_wholeGenome_12 = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome_{ip}.pdf", ip = ["targets", "inputs"])
        PCA_wholeGenome_23 = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome_{ip}.pdf", ip = ["targets", "inputs"])
    else:
        PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome_targets.pdf"
        PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome_targets.pdf"
else:
    PCA_wholeGenome_12 = []
    PCA_wholeGenome_23 = []


if (len(TARGETNAMES) > 2):
    PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf"
    PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
else:
    PCA_atPeaks_12 = []
    PCA_atPeaks_23 = []


if (len(SAMPLENAMES) > 1):
    if (len(INPUTNAMES) > 2):
        correlation_heatmap_wholeGenome_pearson = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_{ip}.pdf", ip = ["targets", "inputs"])
        correlation_heatmap_wholeGenome_spearman = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_{ip}.pdf", ip = ["targets", "inputs"])
    else:
        correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_targets.pdf"
        correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_targets.pdf"
else:
    correlation_heatmap_wholeGenome_pearson = []
    correlation_heatmap_wholeGenome_spearman = []


if (len(TARGETNAMES) > 1):
    correlation_heatmap_atPeaks_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.pdf"
    correlation_heatmap_atPeaks_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.pdf"
else:
    correlation_heatmap_atPeaks_pearson = []
    correlation_heatmap_atPeaks_spearman = []


if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
    peak_type = "BAMPE"
else:
    peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)
    peak_type = "BAM"


### Somatic Variants & CopyWriteR
if (eval(str(config["somatic_variants"]["call_variants"])) == True):
    GATKDIR = "06a_Somatic_variants"
    COPYWRITERDIR = "06b_Copy_Number_Variation"
else:
    GATKDIR = ""
    COPYWRITERDIR = "06b_Copy_Number_Variation"


if (eval(str(config["somatic_variants"]["call_variants"])) == True):
    somatic_variants = expand(os.path.join(GATKDIR, "SV_count_plots/all.samples_{variant}_counts_plot.pdf"), variant = ["SNP", "InDel"])
else:
    somatic_variants = []


if (eval(str(config["copy_number_variation"]["call_CNV"])) == True):
    CNV_results = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized_CNA.corrected/", ''.join(["{target}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_CNA.corrected.bw"])), target = TARGETNAMES)
    CNV_coverage_plot = expand(os.path.join(GATKDIR, "coverage_plots/{target}_plotCoverage.pdf"), target = TARGETNAMES)

    if (eval(str(config["GCbias_correction"]["correct_GCbias"])) == True):
        CNV_results_GCbias = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected_CNA.corrected/", ''.join(["{target}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_GC.corrected_CNA.corrected.bw"])), target = TARGETNAMES)
    else:
        CNV_results_GCbias = []
else:
    CNV_results = []
    CNV_results_GCbias = []
    CNV_coverage_plot = []


### GC-bias correction
if (eval(str(config["GCbias_correction"]["correct_GCbias"])) == True):
    GC_corrected_bam = expand(os.path.join("01_BAM_filtered/GCbias_corrected_files", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_GC.corrected.bam"])), sample = SAMPLENAMES)
    GC_corrected_bw = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_GC.corrected.bw"])), sample = SAMPLENAMES)
else:
    GC_corrected_bam = []
    GC_corrected_bw = []




### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    TARGET = constraint_to(TARGETNAMES),
    INPUT = constraint_to(INPUTNAMES)


ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
if (set(SAMPLENAMES) <= set(RUNNAMES)):
    rule AAA_initialization:
        input:
            fastqc_bam_zip = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            plotFingerprint_results = plotFingerprint_results,
            fragmentSizeDistribution_results = fragmentSizeDistribution_results,
            normalized_bigWig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            raw_bigWig = expand(os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            correlation_heatmap_wholeGenome_pearson = correlation_heatmap_wholeGenome_pearson,
            correlation_heatmap_wholeGenome_spearman = correlation_heatmap_wholeGenome_spearman,
            PCA_wholeGenome_12 = PCA_wholeGenome_12,
            PCA_wholeGenome_23 = PCA_wholeGenome_23,
            peaks = peaks,
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html",
            correlation_heatmap_atPeaks_pearson = correlation_heatmap_atPeaks_pearson,
            correlation_heatmap_atPeaks_spearman = correlation_heatmap_atPeaks_spearman,
            PCA_atPeaks_12 = PCA_atPeaks_12,
            PCA_atPeaks_23 = PCA_atPeaks_23,
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv",
            somatic_variants = somatic_variants,
            CNV_results = CNV_results,
            CNV_coverage_plot = CNV_coverage_plot,
            CNV_results_GCbias = CNV_results_GCbias,
            GC_corrected_bam = GC_corrected_bam,
            GC_corrected_bw = GC_corrected_bw
        shell:
            """
            printf '\033[1;36mPipeline ended!\\n\033[0m'
            """
else:
    missing_samples = '\\n  - '.join(list(set(SAMPLENAMES) - set(RUNNAMES)))
    os.system("printf '\033[1;31m\\n!!! Not all bam files are avalable in the input directory !!!\\n\\nPlease provide files for:\\n  - "+missing_samples+"\\n\\n\033[0m'")

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

if (eval(str(config["bam_features"]["skip_bam_filtering"])) == False):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLE}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), ".bam"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam.bai"])))
        params:
            sample = "{SAMPLE}",
            MAPQ_threshold = config["bam_features"]["MAPQ_threshold"]
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        benchmark:
            "benchmarks/MAPQ_filter/MAPQ_filter---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ and re-indexing...\\n\033[0m'

            $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

            $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}
            """


    if ((eval(str(config["bam_features"]["paired_end"])) == True) & (eval(str(config["bam_features"]["umi_present"])) == True)):
        rule gatk4_markdups_umiAware:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                umi_metrics = "01_BAM_filtered/umi_metrics/{SAMPLE}_UMI_metrics.txt",
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["bam_features"]["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                max((workflow.cores-1), 1)
            benchmark:
                "benchmarks/gatk4_markdups_umiAware/gatk4_markdups_umiAware---{SAMPLE}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.sample}: UMI-aware gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/umi_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                $CONDA_PREFIX/bin/gatk UmiAwareMarkDuplicatesWithMateCigar \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --MAX_EDIT_DISTANCE_TO_JOIN 1 \
                --UMI_METRICS_FILE {output.umi_metrics} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --UMI_TAG_NAME RX \
                --CREATE_INDEX true \
                --VALIDATION_STRINGENCY STRICT \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
    else: # Single-end/no-UMI dedup
        rule gatk4_markdups:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["bam_features"]["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                max((workflow.cores-1), 1)
            benchmark:
                "benchmarks/gatk4_markdups/gatk4_markdups---{SAMPLE}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                $CONDA_PREFIX/bin/gatk MarkDuplicates \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --VALIDATION_STRINGENCY LENIENT \
                --CREATE_INDEX true \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLE}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
        params:
            sample = "{SAMPLE}"
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        benchmark:
            "benchmarks/bam_link__skip_filtering/bam_link__skip_filtering---{SAMPLE}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s $BAM_REAL {output.bam_mdup}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """



rule fastQC_filtered_BAM:
    input:
        bam_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        html = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.html"])),
        zip = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join("02_fastQC_on_BAM_filtered/"),
        sample = "{SAMPLE}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    benchmark:
        "benchmarks/fastQC_filtered_BAM/fastQC_filtered_BAM---{SAMPLE}_benchmark.txt"
    shell:
        """
        mkdir -p 02_fastQC_on_BAM_filtered

        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.bam_mapq}
        """

# ------------------------------------------------------------------------------

rule plotFingerprint:
    input:
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
        input_bam = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input=INPUTNAMES),
        input_bai = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])), input=INPUTNAMES)
    output:
        lorenz_curve_pdf = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/{TARGET}_fingerPrinting_plot.pdf"),
        quality_metrics = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/quality_metrics/{TARGET}_fingerPrinting_quality_metrics.txt")
    params:
        sample = "{TARGET}",
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        input_suffix = "_mapq"+str(config["bam_features"]["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
        read_extension = read_extension,
        blacklist = blacklist
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.out",
        err = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.err"
    benchmark:
        "benchmarks/plotFingerprint/plotFingerprint---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: plotting fingerprint...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/logs
        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/quality_metrics

        INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)

        $CONDA_PREFIX/bin/plotFingerprint \
        -b {input.target_bam} \
        01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        --JSDsample 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        -plot {output.lorenz_curve_pdf} \
        {params.read_extension} \
        --ignoreDuplicates \
        --outQualityMetrics {output.quality_metrics} \
        --labels {params.sample} ${{INPUT_ID}} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule fragmentSizeDistribution:
    input:
        all_bams = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), sample = SAMPLENAMES),
        all_bais = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])), sample = SAMPLENAMES)
    output:
        fragment_distribution_plot = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_plot.pdf",
        fragmentSize_metrics = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt",
        fragmentSize_RawFragmentLengths = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_RawFragmentLengths.tab"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist
    threads:
        max((workflow.cores-1), 1)
    log:
        out = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.out",
        err = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.err"
    benchmark:
        "benchmarks/fragmentSizeDistribution/fragmentSizeDistribution---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting fragment size distribution...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/fragmentSize_distribution/logs

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        --bamfiles {input.all_bams} \
        --binSize 1000000 \
        --blackListFileName {params.blacklist} \
        --samplesLabel {params.labels} \
        --histogram {output.fragment_distribution_plot} \
        --table {output.fragmentSize_metrics} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule normalized_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_bamCoverage_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_bamCoverage_log.err"
    benchmark:
        "benchmarks/normalized_bigWig/normalized_bigWig---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.normalized_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        #--ignoreDuplicates \
	--samFlagExclude 1024 \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule raw_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        raw_bigWig = os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.out",
        err = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.err"
    benchmark:
        "benchmarks/raw_bigWig/raw_bigWig---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating raw coverage bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/raw_coverage/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.raw_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing None \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        #--ignoreDuplicates \
	--samFlagExclude 1024 \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule multiBigwigSummary_wholeGenome_targets:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = TARGETNAMES),
    output:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_targets.npz"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        max((workflow.cores-1), 1)
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_targets_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_targets_log.err"
    benchmark:
        "benchmarks/multiBigwigSummary_wholeGenome_targets/multiBigwigSummary_wholeGenome_targets---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (whole genome - TARGETS)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs

        $CONDA_PREFIX/bin/multiBigwigSummary bins \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_wholeGenome_targets:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_targets.npz"
    output:
        correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_targets.pdf",
        correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_targets.pdf",
        correlation_heatmap_wholeGenome_pearson_matrix = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_targets_matrix.txt",
        correlation_heatmap_wholeGenome_spearman_matrix = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_targets_matrix.txt"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_targets_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_targets_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_targets_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_targets_log.err"
    benchmark:
        "benchmarks/correlations_wholeGenome_targets/correlations_wholeGenome_targets---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (whole genome - TARGETS)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Pearson correlation whole genome RPGC normalized coverage (TARGETS)' \
        --plotFile {output.correlation_heatmap_wholeGenome_pearson} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_pearson_matrix} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Spearman correlation whole genome RPGC normalized coverage (TARGETS)' \
        --plotFile {output.correlation_heatmap_wholeGenome_spearman} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_spearman_matrix} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_wholeGenome_targets:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_targets.npz"
    output:
        PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome_targets.pdf",
        PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome_targets.pdf"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_targets_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_targets_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_targets_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_targets_log.err"
    benchmark:
        "benchmarks/multiBigWig_matrix_wholeGenome_targets/multiBigWig_matrix_wholeGenome_targets---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting PCA (whole genome - TARGETS)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA whole genome (TARGETS): PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_12} > {log.out_12} 2> {log.err_12}

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA whole genome (TARGETS): PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_23} > {log.out_23} 2> {log.err_23}
        """

##

rule multiBigwigSummary_wholeGenome_inputs:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = INPUTNAMES),
    output:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_inputs.npz"
    params:
        labels = ' '.join(INPUTNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        max((workflow.cores-1), 1)
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_inputs_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_inputs_log.err"
    benchmark:
        "benchmarks/multiBigwigSummary_wholeGenome_inputs/multiBigwigSummary_wholeGenome_inputs---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (whole genome - INPUTS)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs

        $CONDA_PREFIX/bin/multiBigwigSummary bins \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_wholeGenome_inputs:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_inputs.npz"
    output:
        correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_inputs.pdf",
        correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_inputs.pdf",
        correlation_heatmap_wholeGenome_pearson_matrix = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome_inputs_matrix.txt",
        correlation_heatmap_wholeGenome_spearman_matrix = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome_inputs_matrix.txt"
    params:
        labels = ' '.join(INPUTNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_inputs_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_inputs_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_inputs_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_inputs_log.err"
    benchmark:
        "benchmarks/correlations_wholeGenome_inputs/correlations_wholeGenome_inputs---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (whole genome - INPUTS)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Pearson correlation whole genome RPGC normalized coverage (INPUTS)' \
        --plotFile {output.correlation_heatmap_wholeGenome_pearson} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_pearson_matrix} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Spearman correlation whole genome RPGC normalized coverage (INPUTS)' \
        --plotFile {output.correlation_heatmap_wholeGenome_spearman} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_spearman_matrix} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_wholeGenome_inputs:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome_inputs.npz"
    output:
        PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome_inputs.pdf",
        PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome_inputs.pdf"
    params:
        labels = ' '.join(INPUTNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_inputs_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_inputs_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_inputs_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_inputs_log.err"
    benchmark:
        "benchmarks/multiBigWig_matrix_wholeGenome_inputs/multiBigWig_matrix_wholeGenome_inputs---allSamples_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting PCA (whole genome - INPUTS)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA whole genome (INPUTS): PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_12} > {log.out_12} 2> {log.err_12}

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA whole genome (INPUTS): PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_23} > {log.out_23} 2> {log.err_23}
        """


# ------------------------------------------------------------------------------

if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    rule macs_callpeak_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
        output:
            peaksPE = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["workflow_configuration"]["sample_config_table"],
            input_suffix = "_mapq"+str(config["bam_features"]["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            read_extension = read_extension,
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.err"
        benchmark:
            "benchmarks/macs_callpeak_PE/macs_callpeak_PE---{TARGET}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep -w {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAMPE \
            -g {params.genomeSize} \
            -q {params.macs_qValue_cutoff} \
            --keep-dup all \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAMPE ${{BROAD}} > {log.err} 2> {log.out}
            """
else:
    rule phantom_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES)
        output:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        params:
            sample = "{TARGET}",
            sample_config_table = config["workflow_configuration"]["sample_config_table"],
            input_suffix = "_mapq"+str(config["bam_features"]["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            max((workflow.cores-1), 1)
        log:
            out = '04_Called_peaks/phantom/logs/{TARGET}.phantom.log'
        benchmark:
            "benchmarks/phantom_SE/phantom_SE---{TARGET}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: calculating phantom peak...\\n\033[0m'
            mkdir -p 04_Called_peaks/phantom/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)

            ${{CONDA_PREFIX}}/bin/Rscript ${{CONDA_PREFIX}}/bin/run_spp.R -rf -c='{input.target_bam}' -i="01_BAM_filtered/${{INPUT_ID}}{params.input_suffix}" -savp -out={output.phantom} &> {log.out}
            """



    rule fragment_length:
        input:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        output:
            fragment_length_phanthom = temp('04_Called_peaks/phantom/{TARGET}.fragment_length')
        benchmark:
            "benchmarks/fragment_length/fragment_length---{TARGET}_benchmark.txt"
        shell:
            """
            awk '{{print $3}}' < {input.phantom} | tr ',' '\\t' | awk '{{if($1!=0) print $1; else print $2}}' > {output.fragment_length_phanthom}
            """



    rule macs_callpeak_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
            phantom = '04_Called_peaks/phantom/{TARGET}.fragment_length'
        output:
            peaksSE = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["workflow_configuration"]["sample_config_table"],
            input_suffix = "_mapq"+str(config["bam_features"]["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.err"
        benchmark:
            "benchmarks/macs_callpeak_SE/macs_callpeak_SE---{TARGET}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep -w {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            EXTSIZEPHANTOM=$(cat {input.phantom}) ${{BROAD}}

            if [ "$EXTSIZEPHANTOM" -lt 1 ]; then
              EXTSIZEPHANTOM=200
            fi

            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAM \
            -g {params.genomeSize} \
            --nomodel \
            -q {params.macs_qValue_cutoff} \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAM \
            --extsize $EXTSIZEPHANTOM > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if (eval(str(config["bam_features"]["skip_bam_filtering"])) == False):
    picard_metrics_file = expand("01_BAM_filtered/MarkDuplicates_metrics/{sample}_MarkDuplicates_metrics.txt", sample = SAMPLENAMES)
    picard_metrics_dir = "01_BAM_filtered/MarkDuplicates_metrics"
else:
    picard_metrics_file = []
    picard_metrics_dir = []


if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    rule multiQC_PE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        benchmark:
            "benchmarks/multiQC_PE/multiQC_PE---allSamples_benchmark.txt"
        shell:
            """
            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """
else:
    rule multiQC_SE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES),
            phanthom = expand('04_Called_peaks/phantom/{target}.phantom.spp.out', target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        benchmark:
            "benchmarks/multiQC_SE/multiQC_SE---allSamples_benchmark.txt"
        shell:
            """
            printf '\033[1;36mGenerating multiQC report...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    rule merge_all_peaks_PE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            max((workflow.cores-1), 1)
        benchmark:
            "benchmarks/merge_all_peaks_PE/merge_all_peaks_PE---allTargets_benchmark.txt"
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*_peaks.*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            $CONDA_PREFIX/bin/bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """
else:
    rule merge_all_peaks_SE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            workflow.cores
        benchmark:
            "benchmarks/merge_all_peaks_SE/merge_all_peaks_SE---allTargets_benchmark.txt"
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*_peaks.*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            $CONDA_PREFIX/bin/bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """



rule multiBigwigSummary_atPeaks:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{target}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), target = TARGETNAMES),
        merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
    output:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        max((workflow.cores-1), 1)
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.err"
    benchmark:
        "benchmarks/all_norm_bigwig/all_norm_bigwig---allTargets_benchmark.txt"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (at peaks)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs

        $CONDA_PREFIX/bin/multiBigwigSummary BED-file \
        --BED {input.merged_peaks_sorted} \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        correlation_heatmap_atPeaks_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.pdf",
        correlation_heatmap_atPeaks_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.pdf",
        correlation_heatmap_atPeaks_pearson_matrix = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks_matrix.txt",
        correlation_heatmap_atPeaks_spearman_matrix = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks_matrix.txt"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.err"
    benchmark:
        "benchmarks/correlations_atPeaks/correlations_atPeaks---allTargets_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (at peaks)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Pearson correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_pearson} \
        --outFileCorMatrix {output.correlation_heatmap_atPeaks_pearson_matrix} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --zMin 0 \
        --zMax 1 \
        --plotTitle 'Spearman correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_spearman} \
        --outFileCorMatrix {output.correlation_heatmap_atPeaks_spearman_matrix} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf",
        PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.err"
    benchmark:
        "benchmarks/PCA_atPeaks/PCA_atPeaks---allTargets_benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting PCA (at peaks)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA at peaks: PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_12} > {log.out_12} 2> {log.err_12}

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA at peaks: PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_23} > {log.out_23} 2> {log.err_23}
        """

# ------------------------------------------------------------------------------

if ((eval(str(config["bam_features"]["paired_end"])) == True)):
    rule MACS_peak_QC_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc.txt")
        params:
            sample_config_table = config["workflow_configuration"]["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        benchmark:
            "benchmarks/MACS_peak_QC_PE/MACS_peak_QC_PE---{TARGET}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep -w {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_PE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        benchmark:
            "benchmarks/aggregate_FRiP_PE/aggregate_FRiP_PE---allTargets_benchmark.txt"
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """

#*******************************************************************
else:
    rule MACS_peak_QC_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAM_peaks.qc.txt")
        params:
            sample_config_table = config["workflow_configuration"]["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAM_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        benchmark:
            "benchmarks/MACS_peak_QC_SE/MACS_peak_QC_SE---{TARGET}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep -w {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_SE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        benchmark:
            "benchmarks/aggregate_FRiP_SE/aggregate_FRiP_SE---allTargets_benchmark.txt"
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# GC bias correction
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Create 2bit genome
rule convert_genome_to_2bit:
    input:
        genome = ancient(genome_fasta)
    output:
        genome_2bit = ''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.2bit'])
    benchmark:
        "benchmarks/convert_genome_to_2bit/convert_genome_to_2bit---benchmark.txt"
    shell:
        """
        printf '\033[1;36mConverting genome fasta to .2bit format...\\n\033[0m'

        $CONDA_PREFIX/bin/faToTwoBit {input.genome} {output.genome_2bit}
        """

# Compute GC bias
rule compute_GC_bias:
    input:
        genome_2bit = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.2bit'])),
        sample_bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        sample_bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
    output:
        GCbiasFrequenciesFile = "01_BAM_filtered/GCbias_corrected_files/GCbias_frequencies_files/{SAMPLE}_GCbiasFrequencies.txt",
        biasPlot = "01_BAM_filtered/GCbias_corrected_files/bias_plots/{SAMPLE}_biasPlot.pdf"
    params:
        sample = "{SAMPLE}",
        genome = genome_fasta,
        genomeSize = genomeSize,
        blacklist = blacklist,
        fragmentLength = str(config["GCbias_correction"]["GCbias_fragment_length"])
    log:
        out = "01_BAM_filtered/GCbias_corrected_files/GCbias_frequencies_files/logs/{SAMPLE}_GCbiasFrequencies.out",
        err = "01_BAM_filtered/GCbias_corrected_files/GCbias_frequencies_files/logs/{SAMPLE}_GCbiasFrequencies.err"
    threads:
        max(math.floor(workflow.cores/2), 1)
    benchmark:
        "benchmarks/compute_GC_bias/compute_GC_bias---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Compute GCbias frequencies...\\n\033[0m'

        $CONDA_PREFIX/bin/computeGCBias \
        --bamfile {input.sample_bam} \
        --effectiveGenomeSize {params.genomeSize} \
        --genome {input.genome_2bit} \
        --fragmentLength {params.fragmentLength} \
        --blackListFileName {params.blacklist} \
        --numberOfProcessors {threads} \
        --biasPlot {output.biasPlot} \
        --GCbiasFrequenciesFile {output.GCbiasFrequenciesFile} > {log.out} 2> {log.err}
        """


# GC bias correction
rule correct_GC_bias:
    input:
        genome_2bit = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.2bit'])),
        sample_bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        sample_bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
        GCbiasFrequenciesFile = "01_BAM_filtered/GCbias_corrected_files/GCbias_frequencies_files/{SAMPLE}_GCbiasFrequencies.txt"
    output:
        corrected_bam = os.path.join("01_BAM_filtered/GCbias_corrected_files", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_GC.corrected.bam"])),
        corrected_bai = os.path.join("01_BAM_filtered/GCbias_corrected_files", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_GC.corrected.bai"]))
    params:
        sample = "{SAMPLE}",
        genome = genome_fasta,
        genomeSize = genomeSize
    log:
        out = "01_BAM_filtered/GCbias_corrected_files/logs/{SAMPLE}_GCbias_correction.out",
        err = "01_BAM_filtered/GCbias_corrected_files/logs/{SAMPLE}_GCbias_correction.err"
    threads:
        max(math.floor(workflow.cores/2), 1)
    benchmark:
        "benchmarks/correct_GC_bias/correct_GC_bias---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Compute GCbias frequencies...\\n\033[0m'

        $CONDA_PREFIX/bin/correctGCBias \
        --bamfile {input.sample_bam} \
        --effectiveGenomeSize {params.genomeSize} \
        --genome {input.genome_2bit} \
        --GCbiasFrequenciesFile {input.GCbiasFrequenciesFile} \
        --correctedFile {output.corrected_bam} \
        --numberOfProcessors {threads} > {log.out} 2> {log.err}

        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.corrected_bam} {output.corrected_bai}
        """


# Normalize GC corrected bam files
rule GCcorrected_normalized_bigWig:
    input:
        corrected_bam = os.path.join("01_BAM_filtered/GCbias_corrected_files", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_GC.corrected.bam"])),
        corrected_bai = os.path.join("01_BAM_filtered/GCbias_corrected_files", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_GC.corrected.bai"]))
    output:
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected/", ''.join(["{SAMPLE}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_GC.corrected.bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized_GC.corrected/logs/{SAMPLE}_bamCoverage_GCcorrection_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized_GC.corrected/logs/{SAMPLE}_bamCoverage_GCcorrection_log.err"
    benchmark:
        "benchmarks/GCcorrected_normalized_bigWig/GCcorrected_normalized_bigWig---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.corrected_bam} \
        -o {output.normalized_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        #--ignoreDuplicates \
	--samFlagExclude 1024 \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  VARIANT CALLING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ----------------------------------------------------------------------------------------
# Create reference genome dictionary
rule make_reference_genome_dictionary:
    input:
        genome = ancient(genome_fasta)
    output:
        genome_dict = ''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.dict'])
    benchmark:
        "benchmarks/make_reference_genome_dictionary/make_reference_genome_dictionary---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating genome dictionary...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk CreateSequenceDictionary REFERENCE={input.genome} OUTPUT={output.genome_dict}
        """
# ----------------------------------------------------------------------------------------

# run base score recalibration (BSQR) of the bams
rule GATK_bam_base_quality_score_recalibration:
    input:
        genome_dict = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",genome_fasta),'.dict'])),
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
    output:
        target_bam_withRG = temp(os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted_withRG.bam"]))),
        bsqr_table = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["bsqr_tables/{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.table"])),
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        target_bam_bsqr_index = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bai"]))
    params:
        sample = "{TARGET}",
        gatk_directory = GATKDIR,
        genome = genome_fasta,
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    log:
        readGroup_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{TARGET}_addReadGroup.log"),
        BaseRecalibrator_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{TARGET}_BaseRecalibrator.log"),
        ApplyBQSR_log = os.path.join(GATKDIR, "recalibrated_bams/logs/{TARGET}_ApplyBQSR.log")
    threads:
        max(math.floor(workflow.cores/2), 1)
    benchmark:
        "benchmarks/GATK_bam_base_quality_score_recalibration/GATK_bam_base_quality_score_recalibration---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Adding Read Groups to filtered bams...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk AddOrReplaceReadGroups \
        -I {input.target_bam} \
        -O {output.target_bam_withRG} \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM {params.sample} &> {log.readGroup_log}


        printf '\033[1;36m{params.sample}: Base Quality Score Recalibration of the deduplicated bam...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk BaseRecalibrator \
        --input {output.target_bam_withRG} \
        --known-sites {params.dbsnp} \
        --output {output.bsqr_table} \
        --reference {params.genome} &> {log.BaseRecalibrator_log}

        $CONDA_PREFIX/bin/gatk ApplyBQSR \
        -R {params.genome} \
        -I {output.target_bam_withRG} \
        --bqsr-recal-file {output.bsqr_table} \
        -O {output.target_bam_bsqr} &> {log.ApplyBQSR_log}

        printf '\033[1;36m{params.sample}: Indexing recalibrated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.target_bam_bsqr} {output.target_bam_bsqr_index}
        """

# ----------------------------------------------------------------------------------------


# compute the coverage
rule plotCoverage_merged_peaks:
    input:
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        concatenation_bed_collapsed_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
    output:
        coverage_plot = os.path.join(GATKDIR, "coverage_plots/{TARGET}_plotCoverage.pdf")
    params:
        label = "{TARGET}",
        blacklist = blacklist
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "coverage_plots/logs/{TARGET}_plotCoverage.log")
    benchmark:
        "benchmarks/plotCoverage_merged_peaks/plotCoverage_merged_peaks---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.label}: Plotting coverage at ChIP (merged) peaks...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCoverage \
        --bamfiles {input.target_bam_bsqr} \
        --labels {params.label} \
        --BED {input.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --plotFile {output.coverage_plot} \
        --plotTitle {params.label} \
        --numberOfProcessors {threads} &> {log.out}
        """

# ----------------------------------------------------------------------------------------

# run gatk haplotype caller
rule GATK_haplotype_calling:
    input:
        concatenation_bed_collapsed_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed",
        target_bam_bsqr = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bam"])),
        target_bam_bsqr_index = os.path.join(GATKDIR, "recalibrated_bams", ''.join(["{TARGET}_mapQ", str(config["bam_features"]["MAPQ_threshold"]), "_sorted_woMT_", DUP, "_bsqr.bai"]))
    params:
        genome = genome_fasta,
        sample = "{TARGET}",
        to_copy_bed = os.path.join(GATKDIR, "all_samples_peak_concatenation_collapsed_sorted.bed")
    output:
        gvcf = temp(os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.g.vcf.gz"]))),
        gvcf_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.g.vcf.gz.tbi"])))
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/logs/{TARGET}_HaplotypeCaller.log")
    benchmark:
        "benchmarks/GATK_haplotype_calling/GATK_haplotype_calling---{TARGET}_benchmark.txt"
    shell:
        """
        cp {input.concatenation_bed_collapsed_sorted} {params.to_copy_bed}

        printf '\033[1;36m{params.sample}: GATK Haplotype calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk HaplotypeCaller \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -I {input.target_bam_bsqr} \
        -O {output.gvcf} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation &> {log.out}
        """


# correct the genotypes that come out of haplotype caller
rule GATK_haplotype_calling_correction:
    input:
        concatenation_bed_collapsed_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed",
        gvcf = os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.g.vcf.gz"])),
        gvcf_idx = os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.g.vcf.gz.tbi"]))
    params:
        sample = "{TARGET}",
        genome = genome_fasta
    output:
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.vcf.gz"]))
    threads:
        max(math.floor(workflow.cores/4), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/logs/{TARGET}_GenotypeGVCFs.log")
    benchmark:
        "benchmarks/GATK_haplotype_calling_correction/GATK_haplotype_calling_correction---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Haplotype call correction...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk GenotypeGVCFs \
        --include-non-variant-sites \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.gvcf} \
        -O {output.vcf} &> {log.out}
        """

# ----------------------------------------------------------------------------------------


# Call SNPs
rule GATK_call_SNPs:
    input:
        concatenation_bed_collapsed_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed",
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.vcf.gz"]))
    params:
        sample = "{TARGET}",
        genome = genome_fasta
    output:
        snp = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp.vcf"]))),
        snp_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp.vcf.idx"])))
    threads:
        max(math.floor(workflow.cores/4), 1)
    log:
        out = os.path.join(GATKDIR, ''.join(["VCF/SNP/logs/{TARGET}_SNP_SelectVariants.log"]))
    benchmark:
        "benchmarks/GATK_call_SNPs/GATK_call_SNPs---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK SNP calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type SNP \
        --select-type NO_VARIATION \
        --select-type-to-exclude INDEL \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.snp} &> {log.out}
        """


# Filter SNPs
rule SnpSift_filter_SNPs:
    input:
        snp_vcf = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp.vcf"]))
    output:
        filtered_snp_vcf = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), ".vcf"])))
    params:
        sample = "{TARGET}",
        DP_snp_threshold = config["somatic_variants"]["DP_snp_threshold"],
        QUAL_snp_threshold = config["somatic_variants"]["QUAL_snp_threshold"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    benchmark:
        "benchmarks/SnpSift_filter_SNPs/SnpSift_filter_SNPs---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift SNP filtering...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift filter '( DP > {params.DP_snp_threshold} & ( QUAL > {params.QUAL_snp_threshold} ))' {input.snp_vcf} > {output.filtered_snp_vcf}
        """


# Annotate SNP
rule SnpSift_annotate_SNPs:
    input:
        filtered_snp_vcf = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), ".vcf"]))
    output:
        filtered_snp_vcf_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf"]))),
        filtered_snp_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz"])),
        filtered_snp_vcf_gz_annotated_idx = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz.tbi"]))
    params:
        sample = "{TARGET}",
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/SNP/logs/{TARGET}_SnpSift_annotate_SNPs.log")
    benchmark:
        "benchmarks/SnpSift_annotate_SNPs/SnpSift_annotate_SNPs---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift SNP annotation...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift annotate -a {params.dbsnp} {input.filtered_snp_vcf} > {output.filtered_snp_vcf_annotated} 2> {log.out}

        printf '\033[1;36m{params.sample}: Filtered SNP vcf bgzipping...\\n\033[0m'
        $CONDA_PREFIX/bin/bgzip -@ {threads} -c {output.filtered_snp_vcf_annotated} > {output.filtered_snp_vcf_annotated_gz}
        $CONDA_PREFIX/bin/tabix -p vcf {output.filtered_snp_vcf_annotated_gz}
        """


# Export SNPs table
rule GATK_vcf2txt_SNPs:
    input:
        filtered_snp_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.vcf.gz"]))
    output:
        filtered_snp_allGT_tb_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_allGT_annotated.txt"]))),
        filtered_snp_tb_annotated = os.path.join(GATKDIR, ''.join(["VCF/SNP/{TARGET}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    params:
        sample = "{TARGET}",
        fileds = ''.join(['"', '" "'.join(config["somatic_variants"]["SnpSift_vcf_fields_to_extract"]), '"'])
    threads:
        max(math.floor(workflow.cores/5), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/SNP/logs/{TARGET}_GATK_vcf2txt_SNPs.log")
    benchmark:
        "benchmarks/GATK_vcf2txt_SNPs/GATK_vcf2txt_SNPs---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Filtered SNP export to txt...\\n\033[0m'

        $CONDA_PREFIX/bin/SnpSift extractFields \
        -s "," \
        -e "." \
        {input.filtered_snp_vcf_annotated_gz} \
        {params.fileds} > {output.filtered_snp_allGT_tb_annotated} 2> {log.out}

        grep -v '0/0' {output.filtered_snp_allGT_tb_annotated} > {output.filtered_snp_tb_annotated}
        """


# Merge all SNPs tables
rule GATK_merge_SNP_tables:
    input:
        snp_txt = expand(os.path.join(GATKDIR, ''.join(["VCF/SNP/{target}_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"])), target=TARGETNAMES)
    params:
        sample = TARGETNAMES,
        gatk_dir = GATKDIR,
        suffix_tb = ''.join(["_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"])
    output:
        merged_snps = os.path.join(GATKDIR, ''.join(["VCF/SNP/all.samples_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    threads: 1
    benchmark:
        "benchmarks/GATK_merge_SNP_tables/GATK_merge_SNP_tables---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mMerging all SNP txt tables in one...\\n\033[0m'")

        import pandas as pd

        merged_table = pd.DataFrame()

        for s in params.sample:
            if os.path.getsize(''.join([GATKDIR, "/VCF/SNP/", s, params.suffix_tb])) > 0:
                tb = pd.read_table(''.join([GATKDIR, "/VCF/SNP/", s, params.suffix_tb]))
                tb.insert(0, "sample_ID", s, allow_duplicates=True)
                merged_table = pd.concat([merged_table, tb], ignore_index = True, sort = False)

        merged_table.to_csv(output.merged_snps, encoding="utf-8", index=False, header=True, sep="\t")



# Plot SNP occurences
rule GATK_plot_SNPs:
    input:
        merged_SNPs = os.path.join(GATKDIR, ''.join(["VCF/SNP/all.samples_", DUP, "_gatk-snp_filtered.DP", str(config["somatic_variants"]["DP_snp_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_snp_threshold"]), "_annotated.txt"]))
    params:
        sample = TARGETNAMES,
        plot_title = ''.join(["SNPs (on ", DUP, " bams): DP > ", str(config["somatic_variants"]["DP_snp_threshold"]), ", QUAL > ", str(config["somatic_variants"]["QUAL_snp_threshold"]), ", w\o 0|0"])
    output:
        snp_plot = os.path.join(GATKDIR, "SV_count_plots/all.samples_SNP_counts_plot.pdf")
    threads: 1
    benchmark:
        "benchmarks/GATK_plot_SNPs/GATK_plot_SNPs---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mPlotting SNP occurences per sample...\\n\033[0m'")

        import pandas as pd
        import numpy as np
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        tb = pd.read_table(input.merged_SNPs)

        occurences = []
        for s in np.flip(params.sample):
            occurences.append(len(tb[tb["sample_ID"] == s].index))

        occurences_tb = pd.DataFrame({'Sample':np.flip(params.sample), 'SNP.counts':occurences})
        plot = occurences_tb.plot.barh(x='Sample', y='SNP.counts', title = params.plot_title, legend=False, xlabel = "SNP count")
        plot.figure.savefig(output.snp_plot, bbox_inches='tight')


# ----------------------------------------------------------------------------------------


# Call InDels
rule GATK_call_InDels:
    input:
        concatenation_bed_collapsed_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed",
        vcf = os.path.join(GATKDIR, ''.join(["VCF/{TARGET}_", DUP, "_gatk.vcf.gz"]))
    output:
        indels = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel.vcf"]))),
        indels_idx = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel.vcf.idx"])))
    params:
        sample = "{TARGET}",
        genome = genome_fasta
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{TARGET}_InDel_SelectVariants.log")
    benchmark:
        "benchmarks/GATK_call_InDels/GATK_call_InDels---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Indel calling...\\n\033[0m'

        $CONDA_PREFIX/bin/gatk SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type INDEL \
        --select-type NO_VARIATION \
        --select-type-to-exclude SNP \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.indels} &> {log.out}
        """


# Filter InDels
rule SnpSift_filter_InDels:
    input:
        indel_vcf = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel.vcf"]))
    output:
        filtered_indel_vcf = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), ".vcf"])))
    params:
        sample = "{TARGET}",
        DP_indel_threshold = config["somatic_variants"]["DP_indel_threshold"],
        QUAL_indel_threshold = config["somatic_variants"]["QUAL_indel_threshold"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{TARGET}_filter_InDels.log")
    benchmark:
        "benchmarks/SnpSift_filter_InDels/SnpSift_filter_InDels---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift InDel filtering...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift filter '( (DP > {params.DP_indel_threshold}) & ( QUAL > {params.QUAL_indel_threshold} ))' {input.indel_vcf} > {output.filtered_indel_vcf} 2> {log.out}
        """


# Annotate InDels
rule SnpSift_annotate_InDels:
    input:
        filtered_indel_vcf = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), ".vcf"]))
    output:
        filtered_indel_vcf_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf"]))),
        filtered_indel_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz"])),
        filtered_indel_vcf_gz_annotated_idx = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz.tbi"]))
    params:
        sample = "{TARGET}",
        dbsnp = config["somatic_variants"]["dbsnp_file"]
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{TARGET}_SnpSift_annotate_InDels.log")
    benchmark:
        "benchmarks/SnpSift_annotate_InDels/SnpSift_annotate_InDels---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: SnpSift InDel annotation...\\n\033[0m'
        $CONDA_PREFIX/bin/SnpSift annotate -a {params.dbsnp} {input.filtered_indel_vcf} > {output.filtered_indel_vcf_annotated} 2> {log.out}

        printf '\033[1;36m{params.sample}: Filtered InDel vcf bgzipping...\\n\033[0m'
        $CONDA_PREFIX/bin/bgzip -@ {threads} -c {output.filtered_indel_vcf_annotated} > {output.filtered_indel_vcf_annotated_gz}
        $CONDA_PREFIX/bin/tabix -p vcf {output.filtered_indel_vcf_annotated_gz}
        """


# Export InDels table
rule GATK_vcf2txt_InDels:
    input:
        filtered_indel_vcf_annotated_gz = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.vcf.gz"]))
    output:
        filtered_indel_allGT_tb_annotated = temp(os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_allGT_annotated.txt"]))),
        filtered_indel_tb_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/{TARGET}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = "{TARGET}",
        fileds = ''.join(['"', '" "'.join(config["somatic_variants"]["SnpSift_vcf_fields_to_extract"]), '"'])
    threads:
        max(math.floor(workflow.cores/5), 1)
    log:
        out = os.path.join(GATKDIR, "VCF/InDel/logs/{TARGET}_GATK_vcf2txt_InDels.log")
    benchmark:
        "benchmarks/GATK_vcf2txt_InDels/GATK_vcf2txt_InDels---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Filtered InDel export to txt...\\n\033[0m'

        $CONDA_PREFIX/bin/SnpSift extractFields \
        -s "," \
        -e "." \
        {input.filtered_indel_vcf_annotated_gz} \
        {params.fileds} > {output.filtered_indel_allGT_tb_annotated} 2> {log.out}

        grep -v '0/0' {output.filtered_indel_allGT_tb_annotated} > {output.filtered_indel_tb_annotated}
        """



# Merge all InDels tables
rule GATK_merge_InDel_tables:
    input:
        indel_txt_annotated = expand(os.path.join(GATKDIR, ''.join(["VCF/InDel/{target}_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"])), target=TARGETNAMES)
    output:
        merged_indels_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/all.samples_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = TARGETNAMES,
        gatk_dir = GATKDIR,
        suffix_tb = ''.join(["_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"])
    threads: 1
    benchmark:
        "benchmarks/GATK_merge_InDel_tables/GATK_merge_InDel_tables---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mMerging all InDel txt tables in one...\\n\033[0m'")

        import pandas as pd

        merged_table = pd.DataFrame()

        for s in params.sample:
            if os.path.getsize(''.join([GATKDIR, "/VCF/InDel/", s, params.suffix_tb])) > 0:
                tb = pd.read_table(''.join([GATKDIR, "/VCF/InDel/", s, params.suffix_tb]))
                tb.insert(0, "sample_ID", s, allow_duplicates=True)
                merged_table = pd.concat([merged_table, tb], ignore_index = True, sort = False)

        merged_table.to_csv(output.merged_indels_annotated, encoding="utf-8", index=False, header=True, sep="\t")



# Plot InDels occurences
rule GATK_plot_InDels:
    input:
        merged_indels_annotated = os.path.join(GATKDIR, ''.join(["VCF/InDel/all.samples_", DUP, "_gatk-indel_filtered.DP", str(config["somatic_variants"]["DP_indel_threshold"]), ".QUAL", str(config["somatic_variants"]["QUAL_indel_threshold"]), "_annotated.txt"]))
    params:
        sample = TARGETNAMES,
        plot_title = ''.join(["InDels (on ", DUP, " bams): DP > ", str(config["somatic_variants"]["DP_indel_threshold"]), ", QUAL > ", str(config["somatic_variants"]["QUAL_indel_threshold"]), ", w\o 0|0"])
    output:
        indel_plot = os.path.join(GATKDIR, "SV_count_plots/all.samples_InDel_counts_plot.pdf")
    threads: 1
    benchmark:
        "benchmarks/GATK_plot_InDels/GATK_plot_InDels---allTargets_benchmark.txt"
    run:
        shell("printf '\033[1;36mPlotting indel occurences per sample...\\n\033[0m'")

        import pandas as pd
        import numpy as np
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        tb = pd.read_table(input.merged_indels_annotated)

        occurences = []
        for s in np.flip(params.sample):
            occurences.append(len(tb[tb["sample_ID"] == s].index))

        occurences_tb = pd.DataFrame({'Sample':np.flip(params.sample), 'InDel.counts':occurences})
        plot = occurences_tb.plot.barh(x='Sample', y='InDel.counts', title = params.plot_title, legend=False, xlabel = "InDel count")
        plot.figure.savefig(output.indel_plot, bbox_inches='tight')



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# generate copywriteR script
rule generate_copywriteR_Rscript:
    output:
        script = os.path.join(COPYWRITERDIR, "CopywriteR_script.R")
    threads: 1
    benchmark:
        "benchmarks/generate_copywriteR_Rscript/generate_copywriteR_Rscript---genomeMappability_benchmark.txt"
    run:
        shell("printf '\033[1;36mGenerating the Rscript used to compute and visualize the log2 fequencies (CopywriteR)...\\n\033[0m'")

        x = """# Load parameters
        args = commandArgs(trailingOnly = TRUE)

        sample_id = as.character(strsplit(grep('--SAMPLENAME*', args, value = TRUE), split = '=')[[1]][[2]])
        target_id = as.character(strsplit(grep('--BAMFILE*', args, value = TRUE), split = '=')[[1]][[2]])
        peaks = as.character(strsplit(grep('--PEAKS*', args, value = TRUE), split = '=')[[1]][[2]])
        cores = as.numeric(strsplit(grep('--CORES*', args, value = TRUE), split = '=')[[1]][[2]])
        data.folder = as.character(strsplit(grep('--COPYWRITERDIR*', args, value = TRUE), split = '=')[[1]][[2]])
        kb.resolution = as.numeric(strsplit(grep('--RESOLUTION*', args, value = TRUE), split = '=')[[1]][[2]])
        genome = as.character(strsplit(grep('--GENOME*', args, value = TRUE), split = '=')[[1]][[2]])
        CNA_threshold = as.numeric(strsplit(grep('--CNATHRESHOLD*', args, value = TRUE), split = '=')[[1]][[2]])
        CNA_color = as.character(strsplit(grep('--LINECOLOR*', args, value = TRUE), split = '=')[[1]][[2]])
        point_size = as.numeric(strsplit(grep('--POINTSIZE*', args, value = TRUE), split = '=')[[1]][[2]])
        point_alpha = as.numeric(strsplit(grep('--POINTALPHA*', args, value = TRUE), split = '=')[[1]][[2]])



        #########################
        # Load libs
        require(dplyr)
        require(ggplot2)


        ### Set-up copywriteR
        bp.param = BiocParallel::SnowParam(workers = cores, type = "SOCK")


        # Run CopyWriteR
        sample_dir = paste0(gsub("/$","",data.folder), "/", sample_id)

        CopywriteR::CopywriteR(sample.control = data.frame(target_id, target_id),
                               destination.folder = sample_dir,
                               reference.folder = paste0(data.folder,"/",genome,"_",kb.resolution,"kb"),
                               capture.regions.file = peaks,
                               bp.param = bp.param,
                               keep.intermediary.files = F)


        # Plot by copyWriteR
        CopywriteR::plotCNA(destination.folder = sample_dir)


        # Read segmentation
        loadRData = function(fileName){
          load(fileName)
          get(ls()[ls() != "fileName"])}

        segmentation = loadRData(paste0(sample_dir, "/CNAprofiles/segment.Rdata"))


        # Read log2_readCounts and plot it
        log2.tb =
          read.table(file = file.path(sample_dir, "CNAprofiles", "log2_read_counts.igv"), header = TRUE) %>%
          #dplyr::filter(Chromosome != "Y") %>%
          dplyr::mutate(mid.point = (Start+End)/2,
                        Chromosome = factor(Chromosome, levels = unique(Chromosome)))
        names(log2.tb)[5] = "log2.value"


        # Re-plot CNA
        ymax = ceiling(max(abs(boxplot(log2.tb$log2.value, na.rm = T)$stats[,1]))) + 1

        CNA.plot =
          ggplot(log2.tb,
                 aes(x = mid.point,
                     y = log2.value)) +
          geom_point(size = point_size,
                     stroke = NA,
                     alpha = point_alpha) +
          theme_classic() +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(color = "black"),
                axis.text.y = element_text(color = "black"),
                panel.spacing = unit(0,'lines'),
                strip.background = element_blank(),
                strip.placement = "outside",
                panel.grid = element_blank(),
                panel.border = element_blank()) +
          geom_hline(yintercept = c(-1,1)*log2(CNA_threshold), linetype = 2, color = CNA_color) +
          geom_vline(xintercept = +Inf, linetype = 3, color = "gray60") +
          ylim(c(-1,1)*ymax) +
          xlab("Chromosome") +
          ylab("log2(copy number)") +
          ggtitle(sample_id) +
          geom_segment(data =
                         dplyr::filter(segmentation$output, abs(seg.mean) > 0) %>%
                         dplyr::mutate(chrom = gsub("23", "X", chrom)) %>%
                         dplyr::mutate(chrom = gsub("24", "Y", chrom)) %>%
                         dplyr::mutate(chrom = factor(chrom)) %>%
                         dplyr::rename(Chromosome = chrom,
                                       log2.value = seg.mean),
                       aes(x = loc.start,
                           xend = loc.end,
                           y = log2.value,
                           yend = log2.value),
                       show.legend = F,
                       inherit.aes = F,
                       color = CNA_color,
                       linewidth = 0.75) +
          facet_grid(~ Chromosome,
                     space = "free_x",
                     scales = "free_x",
                     switch = "x")


        pdf(paste0(sample_dir, "/CNA.plot_", sample_id, "_all.chr_",kb.resolution,"kb.pdf"), width = 16, height = 5)
        print(CNA.plot)
        invisible(dev.off())"""

        with open(str(output.script),"w+") as f:
          f.writelines(x)


# Calculate genome GC_mappability
rule generate_copywriteR_genome_map:
    output:
        genome_mappability = os.path.join(COPYWRITERDIR, ''.join([str(re.sub("_.*$", "", genome_used).lower()),"_",str(config["copy_number_variation"]["kb_bin_resolution"]),"kb/GC_mappability.rda"]))
    params:
        data_folder = os.path.join(home_dir, COPYWRITERDIR),
        resolution = str(config["copy_number_variation"]["kb_bin_resolution"]),
        genome = genome_used
    log:
        out = os.path.join(COPYWRITERDIR, "logs/genome_mappability_file_generation_copywriteR.log")
    benchmark:
        "benchmarks/generate_copywriteR_genome_map/generate_copywriteR_genome_map---genomeMappability_benchmark.txt"
    shell:
        """
        ### Generate genome index
        printf '\033[1;36mGenerating genome mappability files (CopywriteR)...\\n\033[0m'
        $CONDA_PREFIX/bin/Rscript -e 'CopywriteR::preCopywriteR(output.folder = "{params.data_folder}", bin.size = {params.resolution}*1000, ref.genome = "{params.genome}", prefix = "")' &> {log.out}
        """


# Running CopyWriteR
rule CopywriteR:
    input:
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
        peaks = "04_Called_peaks/{TARGET}.filtered."+peak_type+"_peaks.xls",
        copywriteR_script = os.path.join(COPYWRITERDIR, "CopywriteR_script.R"),
        genome_mappability = os.path.join(COPYWRITERDIR, ''.join([str(re.sub("_.*$", "", genome_used).lower()),"_",str(config["copy_number_variation"]["kb_bin_resolution"]),"kb/GC_mappability.rda"]))
    output:
        #CNA_plot = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNA.plot_{TARGET}_all.chr_", str(config["copy_number_variation"]["kb_bin_resolution"]), "kb.pdf"])),
        collapsed_peak = temp(os.path.join(COPYWRITERDIR, "{TARGET}/{TARGET}_collapsed_peaks_chr.bed")),
        log2_read_counts = os.path.join(COPYWRITERDIR, "{TARGET}/CNAprofiles/log2_read_counts.igv")
    params:
        sample_config_table = config["workflow_configuration"]["sample_config_table"],
        outdir = os.path.join(home_dir, COPYWRITERDIR, "{TARGET}"),
        sample = "{TARGET}",
        target = os.path.join(home_dir, "01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        collapsed_peak = os.path.join(home_dir, COPYWRITERDIR, "{TARGET}/{TARGET}_collapsed_peaks_chr.bed"),
        data_folder = os.path.join(home_dir, COPYWRITERDIR),
        resolution = str(config["copy_number_variation"]["kb_bin_resolution"]),
        genome = re.sub("_.*$", "", genome_used).lower(),
        peak_suffix = "04_Called_peaks/{TARGET}.filtered."+peak_type+"_peaks_chr",
        CNA_threshold = abs(config["copy_number_variation"]["CNA_threshold"]),
        CNA_plot_line_colors = config["copy_number_variation"]["CNA_plot_line_colors"],
        CNA_plot_point_size = config["copy_number_variation"]["CNA_plot_point_size"],
        CNA_plot_point_alpha = config["copy_number_variation"]["CNA_plot_point_transparency"]
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        out = os.path.join(COPYWRITERDIR, "logs/{TARGET}_copywriteR.log")
    benchmark:
        "benchmarks/CopywriteR/CopywriteR---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Running CopywriteR...\\n\033[0m'

        mkdir -p {params.outdir}/CNAprofiles
        rm -r {params.outdir}/CNAprofiles
        mkdir -p {params.outdir}

        BROAD=$(grep -w {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

        if [ $BROAD == "false" ]; then
            EXT="narrowPeak"
        else
            EXT="broadPeak"
        fi

        $CONDA_PREFIX/bin/bedtools merge -i {params.peak_suffix}.${{EXT}} > {output.collapsed_peak}

        $CONDA_PREFIX/bin/Rscript {input.copywriteR_script} \
        --SAMPLENAME={params.sample} \
        --BAMFILE={params.target} \
        --PEAKS={params.collapsed_peak} \
        --CORES={threads} \
        --COPYWRITERDIR={params.data_folder} \
        --RESOLUTION={params.resolution} \
        --GENOME={params.genome} \
        --CNATHRESHOLD={params.CNA_threshold} \
        --LINECOLOR={params.CNA_plot_line_colors} \
        --POINTSIZE={params.CNA_plot_point_size} \
        --POINTALPHA={params.CNA_plot_point_alpha} &> {log.out}
        """


# Conversion of CNA counts (filtered to threshold) to bedGraph
rule convert_CNA_to_bedGraph:
    input:
        log2_read_counts = os.path.join(COPYWRITERDIR, "{TARGET}/CNAprofiles/log2_read_counts.igv")
    output:
        bedGraph_filtered = temp(os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bedGraph"])))
    params:
        sample = "{TARGET}",
        CNA_threshold = int(config["copy_number_variation"]["CNA_threshold"])
    threads: 1
    benchmark:
        "benchmarks/convert_CNA_to_bedGraph/convert_CNA_to_bedGraph---{TARGET}_benchmark.txt"
    run:
        shell("printf '\033[1;36m {params.sample}: Filter and converting CNA counts to bedGraph...\\n\033[0m'")

        counts = pd.read_csv(input.log2_read_counts,  sep='\t+', engine='python', skiprows=1)
        counts.iloc[:,4] = 2**counts.iloc[:,4]

        bdg = counts.iloc[:,[0,1,2,4]]
        bdg_filt = bdg[abs(bdg.iloc[:,3]) >= params.CNA_threshold]
        bdg_filt.to_csv(output.bedGraph_filtered, header=False, index=False, sep='\t')



# Create reference genome index if not present
if not any([os.path.exists(''.join([re.sub(".gz", "", config["genomic_annotations"]["genome_fasta"], count=0, flags=0),".fai"]))]):
    rule make_reference_genome_index:
        input:
            genome = ancient(genome_fasta)
        output:
            genome_fai = ''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"])
        benchmark:
            "benchmarks/make_reference_genome_index/make_reference_genome_index---benchmark.txt"
        shell:
            """
            printf '\033[1;36mGenerating genome index (.fai)...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools faidx {input.genome_fasta} -o {output.genome_fai}
            """


# Extract chromosome sizes fro converstion to bigWig
rule extract_chromosome_sizes:
    input:
        genome_fai = ancient(''.join([re.sub(".gz", "", genome_fasta, count=0, flags=0),".fai"]))
    output:
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt"
    benchmark:
        "benchmarks/extract_chromosome_sizes/extract_chromosome_sizes---benchmark.txt"
    shell:
        """
        printf '\033[1;36mExtracting chromosome sizes from genome index (.fai)...\\n\033[0m'
        cut -f1,2 {input.genome_fai} > {output.chrSizes}
        """


# Correct CNV on RPGC normalized bigWigs
rule correct_CNV_normalized_bigWigs:
    input:
        bedGraph_filtered = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bedGraph"])),
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt",
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"]))
    output:
        sorted_bedGraph = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts_sorted.bedGraph"])),
        CNA_corrected_bw = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_CNA.corrected/", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_CNA.corrected.bw"]))
    params:
        sample = "{TARGET}",
        bigWig_filtered = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bw"])),
        bw_binSize = config["bigWig_binSize"]
    threads:
        max((workflow.cores-1), 1)
    log:
        bedGraphToBigWig_out = os.path.join(COPYWRITERDIR, "logs/{TARGET}_bedGraphToBigWig.log"),
        bigwigCompare_out = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_CNA.corrected/log/{TARGET}_CNA.correction_bigwig.out")
    benchmark:
        "benchmarks/correct_CNV_normalized_bigWigs/correct_CNV_normalized_bigWigs---{TARGET}_benchmark.txt"
    shell:
        """
        sort -k1,1 -k2,2n {input.bedGraph_filtered} > {output.sorted_bedGraph}
        NROWSBDG=$(wc -l {output.sorted_bedGraph} | head -n 1 | cut -f 1 -d ' ')

        if [[ $NROWSBDG -gt 0 ]]
        then
            printf '\033[1;36m{params.sample}: signal correction for CNA...\\n\033[0m'

            $CONDA_PREFIX/bin/bedGraphToBigWig {output.sorted_bedGraph} {input.chrSizes} {params.bigWig_filtered} &> {log.bedGraphToBigWig_out}

            $CONDA_PREFIX/bin/bigwigCompare \
            --bigwig1 {input.normalized_bigWig} \
            --bigwig2 {params.bigWig_filtered} \
            -o {output.CNA_corrected_bw} \
            --operation ratio \
            --binSize {params.bw_binSize} \
            -p {threads} \
            -of bigwig \
            --pseudocount 0 1 &> {log.bigwigCompare_out}
        else
            cp {input.normalized_bigWig} {output.CNA_corrected_bw}
        fi
        """


# Correct GCbiased bigWigs
rule correct_CNV_GCbias_corrected_bigWigs:
    input:
        bedGraph_filtered = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bedGraph"])),
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt",
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected/", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_GC.corrected.bw"])),
        sorted_bedGraph = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts_sorted.bedGraph"]))
    output:
        CNA_corrected_bw = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected_CNA.corrected/", ''.join(["{TARGET}_mapq", str(config["bam_features"]["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), "_GC.corrected_CNA.corrected.bw"]))
    params:
        sample = "{TARGET}",
        bigWig_filtered = os.path.join(COPYWRITERDIR, ''.join(["{TARGET}/CNAprofiles/{TARGET}_filtered.abs.", str(config["copy_number_variation"]["CNA_threshold"]), "_linear_CNAcounts.bw"])),
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/2), 1)
    log:
        bedGraphToBigWig_out = os.path.join(COPYWRITERDIR, "logs/{TARGET}_bedGraphToBigWig.log"),
        bigwigCompare_out = os.path.join("03_bigWig_bamCoverage/RPGC_normalized_GC.corrected_CNA.corrected/log/{TARGET}_CNA.correction_GC.corrected.bigwig.out")
    benchmark:
        "benchmarks/correct_CNV_GCbias_corrected_bigWigs/correct_CNV_GCbias_corrected_bigWigs---{TARGET}_benchmark.txt"
    shell:
        """
        NROWSBDG=$(wc -l {input.sorted_bedGraph} | head -n 1 | cut -f 1 -d ' ')

        if [[ $NROWSBDG -gt 0 ]]
        then
            printf '\033[1;36m{params.sample}: signal correction for CNA of GCbias corrected bigWigs...\\n\033[0m'

            #$CONDA_PREFIX/bin/bedGraphToBigWig {input.sorted_bedGraph} {input.chrSizes} {params.bigWig_filtered} &> {log.bedGraphToBigWig_out}

            $CONDA_PREFIX/bin/bigwigCompare \
            --bigwig1 {input.normalized_bigWig} \
            --bigwig2 {params.bigWig_filtered} \
            -o {output.CNA_corrected_bw} \
            --operation ratio \
            --binSize {params.bw_binSize} \
            -p {threads} \
            -of bigwig \
            --pseudocount 0 1 &> {log.bigwigCompare_out}
        else
            cp {input.normalized_bigWig} {output.CNA_corrected_bw}
        fi
        """



# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
