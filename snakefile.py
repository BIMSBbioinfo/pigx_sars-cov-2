# PiGx SARS-CoV-2 wastewater sequencing pipeline
#
# Copyright © 2021 Akalin lab.
#
# This file is part of the PiGx SARS-CoV-2 wastewater sequencing pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Snakefile for PiGx SARS-CoV-2 wastewater sequencing pipeline
"""

import os
import csv
import yaml
from itertools import chain
import re
import inspect
from pathlib import Path


# Helper functions


def toolArgs(name):
    """
    Helper function to retrieve all preset arguments for a given tool.
    """
    if "args" in config["tools"][name]:
        return config["tools"][name]["args"]
    else:
        return ""


def tool(name):
    """
    Helper function to bundle a tool with its preset arguments.
    """
    cmd = config["tools"][name]["executable"]
    args = toolArgs(name)

    if args:
        return cmd + " " + args
    else:
        return cmd


# Convenience function to access fields of sample sheet columns that
# match the predicate.  The predicate may be a string.
def lookup(column, predicate, fields=[]):
    """
    Convenience function to access fields of sample sheet columns that
    match the predicate.  The predicate may be a string.
    """
    if inspect.isfunction(predicate):
        records = [line for line in SAMPLE_SHEET if predicate(line[column])]
    else:
        records = [line for line in SAMPLE_SHEET if line[column] == predicate]
    return [record[field] for record in records for field in fields]


def fastq_ext(fastq_file):
    "Function to determine the fastq file extension"
    root, ext = os.path.splitext(fastq_file)
    if ext == ".gz":
        root_root, root_ext = os.path.splitext(root)
        ext = "".join([root_ext, ext])
    return ext


# WIP create a dummy entry if no variant is found - use this as long as the input-function solution doesn't work
def no_variant_vep(sample, lofreq_output):
    """
    Work-around to create dummy entries in lofreq output to ensure smooth VEP
    running in case no variants are found by lofreq (?)
    """
    content = open(lofreq_output.format(sample=sample), "r").read()
    if re.findall("^NC", content, re.MULTILINE):  # regex ok or not?
        # trigger vep path
        logger.info("File can be used for downstream processing")
    else:
        # write smth so that vep does not crash - deal with everything later in the variant_report
        logger.info("adding dummy entry to vcf file, because no variants were found")
        open(lofreq_output.format(sample=sample), "a").write(
            "NC_000000.0\t00\t.\tA\tA\t00\tPASS\tDP=0;AF=0;SB=0;DP4=0,0,0,0"
        )


# Input functions

def samtools_sort_preprimertrim_input(wildcards):
    sample = wildcards[0]

    if START_POINT == "bam":
        # take bam files directly from the reads dir
        input_file = os.path.join(INPUT_DIR, f"{sample}.bam")

    else:
        input_file = os.path.join(MAPPED_READS_DIR, f"{sample}_aligned.bam")

    return input_file


# function to pass read files to trim/filter/qc improvement
def trim_reads_input(args):
    """
    Get a list of all files related to a sample from the sample sheet. Helps in
    working with both single and paired end data.
    """
    sample = args[0]
    return [
        os.path.join(INPUT_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]


def bwa_align_input(args):
    """
    Input function to return the trimmed read files belonging to a sample, 
    independent of their end mode (single vs paired).
    """
    sample = args[0]
    reads_files = [
        os.path.join(INPUT_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]
    if len(reads_files) > 1:
        return [
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq.gz".format(sample=sample)
            ),
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq.gz".format(sample=sample)
            ),
        ]
    elif len(reads_files) == 1:
        return [
            os.path.join(
                TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz".format(sample=sample)
            )
        ]


def lofreq_input(wildcards):
    sample = wildcards[0]

    if RUN_IVAR_PRIMER_TRIMMING and not START_POINT == "bam":
        file_descript = "_aligned_sorted_primer-trimmed_sorted"

    else:
        file_descript = "_aligned_sorted"        

    files = {
        "aligned_bam": os.path.join(
            MAPPED_READS_DIR,
            f"{sample}{file_descript}.bam"),
        "aligned_bai":  os.path.join(
            MAPPED_READS_DIR,
            f"{sample}{file_descript}.bai"),
        "ref": os.path.join(
            INDEX_DIR,
            f"{os.path.basename(REFERENCE_FASTA)}"
            )
        }

    return files


def vcf2csv_input(wildcards):
    sample = wildcards[0]

    if START_POINT == "vcf":
        # take vcf files directly from the reads dir
        input_file = os.path.join(INPUT_DIR, f"{sample}.vcf")

    else:
        input_file = os.path.join(VARIANTS_DIR, f"{sample}_snv.vcf")

    return input_file


def vep_input(wildcards):
    sample = wildcards[0]

    input={"db_dir": VEP_DB}

    if START_POINT == "vcf":
        # take vcf files directly from the reads dir
        input_file = os.path.join(INPUT_DIR, f"{sample}.vcf")

    else:
        input_file = os.path.join(VARIANTS_DIR, f"{sample}_snv.vcf")

    input["input_file"] = input_file

    return input


# dynamically define the multiqc input files created by FastQC and fastp
# TODO add kraken reports per sample
def multiqc_input(args):
    """
    Dynamically define the multiqc input files created by FastQC and fastp for
    each sample.
    """
    sample = args[0]
    reads_files = [
        os.path.join(INPUT_DIR, f)
        for f in lookup("name", sample, ["reads", "reads2"])
        if f
    ]
    # read_num is either ["_R1", "_R2"] or [""] depending on number of read files
    read_num = [
        "_R" + str(f) if len(reads_files) > 1 else ""
        for f in range(1, len(reads_files) + 1)
    ]
    se_or_pe = ["pe" if len(reads_files) > 1 else "se"]
    files = [
        # fastp on raw files
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}_{end}_fastp.html"),
            sample=sample,
            end=se_or_pe,
        ),
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}_{end}_fastp.json"),
            sample=sample,
            end=se_or_pe,
        ),
        # fastqc on raw files
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}{read_num}_fastqc.html"),
            sample=sample,
            read_num=read_num,
        ),
        expand(
            os.path.join(FASTQC_DIR, "{sample}", "{sample}{read_num}_fastqc.zip"),
            sample=sample,
            read_num=read_num,
        ),
        # fastqc after trimming
        expand(
            os.path.join(
                FASTQC_DIR, "{sample}", "{sample}_trimmed{read_num}_fastqc.html"
            ),
            sample=sample,
            read_num=read_num,
        ),
        expand(
            os.path.join(
                FASTQC_DIR, "{sample}", "{sample}_trimmed{read_num}_fastqc.zip"
            ),
            sample=sample,
            read_num=read_num,
        ),
        # fastqc after primer trimming
        expand(
            os.path.join(
                FASTQC_DIR,
                "{sample}",
                "{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.html",
            ),
            sample=sample,
        ),
        expand(
            os.path.join(
                FASTQC_DIR,
                "{sample}",
                "{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.zip",
            ),
            sample=sample,
        ),
    ]
    return list(chain.from_iterable(files))


def render_qc_report_input(wildcards):
    sample = wildcards[0]

    input = {
        "script": os.path.join(SCRIPTS_DIR, "renderReport.R"),
        "report": os.path.join(SCRIPTS_DIR, "report_scripts", "qc_report_per_sample.Rmd"),
        "header": os.path.join(REPORT_DIR, "_navbar.html"),
        "coverage": os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
        "logo": LOGO
    }

    if START_POINT != "bam":
        input["multiqc"] = os.path.join(MULTIQC_DIR, "{sample}", "multiqc_report.html")

    return input


def render_qc_report_params(wildcards, input, output = None, threads = None, resources = None):
    params = {"rscript_exec": RSCRIPT_EXEC}

    if "multiqc" in input.keys():
        params["multiqc_ran"]      = True
        params["multiqc_rel_path"] = input.multiqc[len(REPORT_DIR) + 1 :]

    else:
        params["multiqc_ran"]      = False

    return params


SAMPLE_SHEET_CSV    = config["locations"]["sample-sheet"]
MUTATION_SHEET_CSV  = config["locations"]["mutation-sheet"]
INPUT_DIR           = config["locations"]["input-dir"]
REFERENCE_FASTA     = config["locations"]["reference-fasta"]
AMPLICONS_BED       = config["locations"]["amplicons-bed"]
MUTATIONS_BED       = config["locations"]["mutations-bed"]
KRAKEN_DB           = config["locations"]["kraken-db-dir"]
KRONA_DB            = config["locations"]["krona-db-dir"]
VEP_DB              = config["locations"]["vep-db-dir"]
OUTPUT_DIR          = config["locations"]["output-dir"]

kraken_dl_params = config["databases"]["kraken2"]
krona_dl_params  = config["databases"]["krona"]
vep_dl_params    = config["databases"]["vep"]

KRAKEN_DB_URL        = kraken_dl_params["archive-url"]
KRAKEN_DB_DOWNSAMPLE = kraken_dl_params["downsample-db"]
KRAKEN_DB_MAX_SIZE   = kraken_dl_params["max-db-size-bytes"]

KRONA_DB_USE_PREBUILT = krona_dl_params["use-prebuilt"]
KRONA_DB_URL          = krona_dl_params["archive-url"]

VEP_DB_URL = vep_dl_params["archive-url"]

# TODO: get default read length from multiqc
parameters = config["parameters"]

# vep parameters
VEP_BUFFER_SIZE         = parameters["vep"]["buffer-size"]
SPECIES                 = parameters["vep"]["species"]
VEP_TRANSCRIPT_DISTANCE = parameters["vep"]["transcript-distance"]
VEP_DB_VERSION          = parameters["vep"]["db-version"]

# ivar parameters
IVAR_QUALITY_CUTOFF = parameters["ivar_trimming"]["quality-cutoff"]
IVAR_LENGTH_CUTOFF  = parameters["ivar_trimming"]["length-cutoff"]
IVAR_WINDOW_WIDTH   = parameters["ivar_trimming"]["window-width"]

# mutation regression parameters
MUTATION_COVERAGE_THRESHOLD = parameters['reporting']['mutation-coverage-threshold']

# deconvolution parameters
MUTATION_DEPTH_THRESHOLD    = parameters["reporting"]["mutation-depth-threshold"]
DECONVOLUTION_METHOD = config["deconvolution"]["method"]

START_POINT = config["control"]["start"].lower()
TARGETS     = config["control"]["targets"]

RUN_IVAR_PRIMER_TRIMMING = config["control"]["run-ivar-primer-trimming"]

INDEX_DIR         = os.path.join(OUTPUT_DIR, 'index')
TRIMMED_READS_DIR = os.path.join(OUTPUT_DIR, 'trimmed_reads')
LOG_DIR           = os.path.join(OUTPUT_DIR, 'logs')
MAPPED_READS_DIR  = os.path.join(OUTPUT_DIR, 'mapped_reads')
VARIANTS_DIR      = os.path.join(OUTPUT_DIR, 'variants')
MUTATIONS_DIR     = os.path.join(OUTPUT_DIR, 'mutations')
KRAKEN_DIR        = os.path.join(OUTPUT_DIR, 'kraken')
COVERAGE_DIR      = os.path.join(OUTPUT_DIR, 'coverage')
REPORT_DIR        = os.path.join(OUTPUT_DIR, 'report')
FASTQC_DIR        = os.path.join(REPORT_DIR, 'fastqc')
MULTIQC_DIR       = os.path.join(REPORT_DIR, 'multiqc')
SCRIPTS_DIR       = os.path.join(config['locations']['pkglibexecdir'], 'scripts/')
TMP_DIR           = os.path.join(config['locations']['output-dir'], 'pigx_work')

if os.getenv("PIGX_UNINSTALLED"):
    LOGO = os.path.join(config['locations']['pkgdatadir'], "images/Logo_PiGx.png")
else:
    LOGO = os.path.join(config['locations']['pkgdatadir'], "Logo_PiGx.png")

BWA_EXEC             = tool("bwa")
FASTP_EXEC           = tool("fastp")
FASTQC_EXEC          = tool("fastqc")
GUNZIP_EXEC          = tool("gunzip")
GZIP_EXEC            = tool("gzip")
MULTIQC_EXEC         = tool("multiqc")
IMPORT_TAXONOMY_EXEC = tool("import_taxonomy")
KRAKEN2_EXEC         = tool("kraken2")
KRAKEN2_BUILD_EXEC   = tool("kraken2_build")
LOFREQ_EXEC          = tool("lofreq")
PYTHON_EXEC          = tool("python")
RSCRIPT_EXEC         = tool("Rscript")
SAMTOOLS_EXEC        = tool("samtools")
VEP_EXEC             = tool("vep")
IVAR_EXEC            = tool("ivar")

KRONA_TAXUPDATE = os.path.join(
    os.path.dirname(IMPORT_TAXONOMY_EXEC),
    "../share/krona-tools/updateTaxonomy.sh")

# start file types and the rules that will be skipped:
#   fastq.gz: 
#     None
#   bam:
#     * fastp
#     * fastp_se
#     * bwa_align
#     * samtools_filter_aligned
#     * samtools_filter_unaligned
#     * ivar_primer_trim
#     * samtools_sort_postprimertrim
#     * samtools_index_postprimertrim
#     * fastqc_raw_se
#     * fastqc_raw
#     * fastqc_trimmed_se
#     * fastqc_trimmed_pe
#     * fastqc_primer_trimmed
#   vcf:
#     * bam2fastq
#     * bwa_align
#     * bwa_index
#     * create_mutations_summary
#     * create_overviewQC_table
#     * create_variants_summary
#     * fastp
#     * fastp_se
#     * fastqc_primer_trimmed
#     * fastqc_raw
#     * fastqc_raw_se
#     * fastqc_trimmed_pe
#     * fastqc_trimmed_se
#     * get_qc_table
#     * ivar_primer_trim
#     * kraken
#     * krona_report
#     * lofreq
#     * multiqc
#     * render_index
#     * render_kraken2_report
#     * render_qc_report
#     * run_mutation_regression
#     * samtools_bedcov
#     * samtools_coverage
#     * samtools_filter_aligned
#     * samtools_filter_unaligned
#     * samtools_index_postprimertrim
#     * samtools_index_preprimertrim
#     * samtools_sort_postprimertrim
#     * samtools_sort_preprimertrim


## Load sample sheet
with open(SAMPLE_SHEET_CSV, 'r') as fp:
  rows =  [row for row in csv.reader(fp, delimiter=',')]
  header = rows[0]; rows = rows[1:]
  SAMPLE_SHEET = [dict(zip(header, row)) for row in rows]

SAMPLES = [line['name'] for line in SAMPLE_SHEET]

# predefine files for targets
final_report_files = (
    expand(os.path.join(REPORT_DIR, '{sample}.variantreport_p_sample.html'), sample=SAMPLES)
)

if START_POINT != "vcf":
    final_report_files = (
        final_report_files +
        expand(os.path.join(REPORT_DIR, '{sample}.qc_report_per_sample.html'), sample=SAMPLES) +    
        [os.path.join(REPORT_DIR, 'index.html')]
    )

if START_POINT not in ["bam", "vcf"]:
    final_report_files = (
        final_report_files +
        expand(os.path.join(REPORT_DIR, '{sample}.taxonomic_classification.html'), sample=SAMPLES) +
        expand(os.path.join(REPORT_DIR, '{sample}.Krona_report.html'), sample=SAMPLES)
    )

targets = {
    'help': {
        'description': "Print all rules and their descriptions.",
        'files': []
    },
    'final_reports': {
        'description': "Produce a comprehensive report. This is the default target.",
        'files': final_report_files
    },
    "deconvolution": {
        "description": (
            "Run deconvolution for all provided samples and create a summary "
            "table containing abundances from all samples."),
        "files": [os.path.join(VARIANTS_DIR, "data_variant_plot.csv")]
    },
    'lofreq': {
        'description': "Call variants and produce .vcf file and overview .csv file.",
        'files': (
            expand(os.path.join(VARIANTS_DIR, '{sample}_snv.csv'), sample=SAMPLES)
        )
    },
    'multiqc': {
        'description': "Create MultiQC reports for including raw and trimmed reads.",
        'files': (
            expand(os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html'), sample=SAMPLES)
        )
    }
}

selected_targets = config["control"]["targets"]
OUTPUT_FILES = list(chain.from_iterable([targets[name]['files'] for name in selected_targets]))

run_params_info = (
    f"Run parameters:\n"
    f"\tStart point: {START_POINT}\n"
    f"\tTargets: {TARGETS}\n"
)

logger.info(run_params_info)

rule all:
    input: OUTPUT_FILES

# Record any existing output files, so that we can detect if they have
# changed.
expected_files = {}
onstart:
    if OUTPUT_FILES:
        for name in OUTPUT_FILES:
            if os.path.exists(name):
                expected_files[name] = os.path.getmtime(name)

# Print generated target files.
onsuccess:
    if OUTPUT_FILES:
        # check if any existing files have been modified
        generated = []
        for name in OUTPUT_FILES:
            if name not in expected_files or os.path.getmtime(name) != expected_files[name]:
                generated.append(name)
        if generated:
            print("The following files have been generated:")
            for name in generated:
                logger.info("  - {}".format(name))


# ANNOT: If the any of the databases is not present, download it to the location
# specified in the settings file.
rule download_kraken_db:
    output: directory(KRAKEN_DB)
    params:
        kraken2_build=KRAKEN2_BUILD_EXEC,
        dl_url=KRAKEN_DB_URL,
        downsample_db=KRAKEN_DB_DOWNSAMPLE,
        max_db_size=KRAKEN_DB_MAX_SIZE
    log:
       os.path.join(LOG_DIR, "database_downloads", "download_kraken_db.log")
    script: "snakefile_scripts/rule_download_kraken_db.py"
    
rule download_krona_db:
    output: directory(KRONA_DB)
    params:
        krona_update_tax_script=KRONA_TAXUPDATE,
        dl_url=KRONA_DB_URL,
        use_prebuilt=KRONA_DB_USE_PREBUILT
    log:
       os.path.join(LOG_DIR, "database_downloads", "download_krona_db.log")
    script: "snakefile_scripts/rule_download_krona_db.py"
    
rule download_vep_db:
    output: directory(VEP_DB)
    params:
        dl_url=VEP_DB_URL
    log:
       os.path.join(LOG_DIR, "database_downloads", "download_vep_db.log")
    script: "snakefile_scripts/rule_download_vep_db.py"

# Trimming in three steps: general by qual and cutoff, get remaining adapters out, get remaining primers out

# TODO the output suffix should be dynamic depending on the input
# TODO with the use of fastp the use of fastqc becomes partly reduntant, fastqc should be removed or adjusted
rule fastp:
    input: trim_reads_input
    output:
        r1 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R1.fastq.gz"),
        r2 = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R2.fastq.gz"),
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_pe_fastp.html'),
        json = os.path.join(FASTQC_DIR, '{sample}', '{sample}_pe_fastp.json')
    log: os.path.join(LOG_DIR, 'fastp_{sample}.log')
    shell: """
         {FASTP_EXEC} -i {input[0]} -I {input[1]} -o {output.r1} -O {output.r2} --html {output.html} --json {output.json} >> {log}t 2>&1
     """

rule fastp_se:
    input: trim_reads_input
    output:
        r = os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz"),
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_se_fastp.html'),
        json = os.path.join(FASTQC_DIR, '{sample}', '{sample}_se_fastp.json')
    log: os.path.join(LOG_DIR, 'fastp_{sample}.log')
    shell: """
        {FASTP_EXEC} -i {input[0]} -o {output.r} --html {output.html} --json {output.json} >> {log}t 2>&1
    """

rule bwa_index:
    input: REFERENCE_FASTA
    output:
      ref=os.path.join(INDEX_DIR, os.path.basename(REFERENCE_FASTA)),
      index=os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    log: os.path.join(LOG_DIR, 'bwa_index.log')
    shell: """
        mkdir -p {INDEX_DIR};
        ln -sf {input} {INDEX_DIR};
        cd {INDEX_DIR};
        {BWA_EXEC} index {output.ref} >> {log} 2>&1
        """

# alignment works with both single and paired-end files
rule bwa_align:
    input:
        fastq = bwa_align_input,
        ref = os.path.join(INDEX_DIR, "{}".format(os.path.basename(REFERENCE_FASTA))),
        index = os.path.join(INDEX_DIR, "{}.bwt".format(os.path.basename(REFERENCE_FASTA)))
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    params:
        threads = 4
    log: os.path.join(LOG_DIR, 'bwa_align_{sample}.log')
    shell: "{BWA_EXEC} mem -t {params.threads} {input.ref} {input.fastq} > {output} 2>> {log} 3>&2"

# TODO verify that subsequent tools do not require filtering for proper pairs
# NOTE verification of flags can be done with "samtools view -h -f 4 <file.sam|file.bam> | samtools flagstat -
rule samtools_filter_aligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned.bam')
    params:
        # add 'proper-pair' filter (-f 2) if sample is paired-end
        proper_pair = lambda wc: "-f 2" if len(trim_reads_input(wc))>1 else ""
    log: os.path.join(LOG_DIR, 'samtools_filter_aligned_{sample}.log')
    shell: # exclude (F) reads that are not mapped (4) and supplementary (2048)
        "{SAMTOOLS_EXEC} view -bh {params.proper_pair} -F 4 -F 2048 {input} > {output} 2>> {log} 3>&2"

rule samtools_filter_unaligned:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_tmp.sam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    log: os.path.join(LOG_DIR, 'samtools_filter_unaligned_{sample}.log')
    shell: # keep (-f) reads that are unmapped (4)
        "{SAMTOOLS_EXEC} view -bh -f 4 {input} > {output} 2>> {log} 3>&2"

rule samtools_sort_preprimertrim:
    input: samtools_sort_preprimertrim_input
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"

rule samtools_index_preprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"

rule ivar_primer_trim:
    # FIXME Currently this rule may be skipped by giving an extra option.
    # I am wondering whether you could also just set the parameters so that
    # the rule leaves everything as is. But that might be more trouble than it
    # is worth.
    input:
        primers=AMPLICONS_BED,
        aligned_bam=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bam"),
        aligned_bai=os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted.bai"),
    output:
        os.path.join(MAPPED_READS_DIR, "{sample}_aligned_sorted_primer-trimmed.bam"),
    params:
        output=lambda wildcards, output: os.path.splitext(f"{output}")[0],
        quality_cutoff=IVAR_QUALITY_CUTOFF,
        length_cutoff=IVAR_LENGTH_CUTOFF,
        window_width=IVAR_WINDOW_WIDTH
    log:
        os.path.join(LOG_DIR, "ivar_{sample}.log"),
    shell:
        """
        {IVAR_EXEC} trim \
        -b {input.primers} \
        -p {params.output} \
        -i {input.aligned_bam} \
        -q {params.quality_cutoff} \
        -m {params.length_cutoff} \
        -s {params.window_width} \
        -e \
        >> {log} 2>&1
        """


rule samtools_sort_postprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    log: os.path.join(LOG_DIR, 'samtools_sort_{sample}.log')
    shell: "{SAMTOOLS_EXEC} sort -o {output} {input} >> {log} 2>&1"

rule samtools_index_postprimertrim:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bai')
    log: os.path.join(LOG_DIR, 'samtools_index_{sample}.log')
    shell: "{SAMTOOLS_EXEC} index {input} {output} >> {log} 2>&1"


# FIXME: single-end version needed
# NOTE: fastqc does not process reads in pairs. files are processed as single units.
# ANNOT: Do quality control on single end sample fastq read files.
rule fastqc_raw_se:
    input: trim_reads_input
    output:
        # all outputs are provided to ensure atomicity
        rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_fastqc.zip'),
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_raw.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    run:
        # renaming the ".fastq.gz" suffix to "_fastqc.html"
        tmp_output = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.html')
        tmp_zip = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.zip')
        shell("""{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1;
                if [[ {tmp_output} != {wildcards.sample}_fastqc.html ]]; then
                    mv {params.output_dir}/{tmp_output} {output.rep} &&\
                    mv {params.output_dir}/{tmp_zip} {output.zip}
                fi """)


# FIXME: or discard completely and change multiqc to use fastp --> fastp rule would have to be adjusted to create reasonable outputs
# ANNOT: Do quality control on paired end sample fastq read files.
# FIXME: Should this rule be called `fastqc_raw_pe`?
rule fastqc_raw:
    input: trim_reads_input
    output:
        r1_rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R1_fastqc.html'),
        r1_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R1_fastqc.zip'),
        r2_rep = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R2_fastqc.html'),
        r2_zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_R2_fastqc.zip') # all outputs are provided to ensure atomicity
    log: [os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R1.log'), os.path.join(LOG_DIR, 'fastqc_{sample}_raw_R2.log')]
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    run:
        # renaming the ".fastq.gz" suffix to "_fastqc.html"
        tmp_R1_output = os.path.basename(input[0]).replace(fastq_ext(input[0]), '_fastqc.html')
        tmp_R1_zip = os.path.basename(input[0]).replace(fastq_ext(input[0]),  '_fastqc.zip')
        tmp_R2_output = os.path.basename(input[1]).replace(fastq_ext(input[0]),'_fastqc.html')
        tmp_R2_zip = os.path.basename(input[1]).replace(fastq_ext(input[0]),'_fastqc.zip')
        shell("""{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1;
                if [[ {tmp_R1_output} != {wildcards.sample}_R1_fastqc.html ]]; then
                    mv {params.output_dir}/{tmp_R1_output} {output.r1_rep} &&\
                    mv {params.output_dir}/{tmp_R1_zip} {output.r1_zip} &&\
                    mv {params.output_dir}/{tmp_R2_output} {output.r2_rep} &&\
                    mv {params.output_dir}/{tmp_R2_zip} {output.r2_zip}
                fi """)

# TODO: can probably be done by using bwa_align_input, no seperate functions neccessary?
rule fastqc_trimmed_se:
    input: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed.fastq.gz")
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_trimmed.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

rule fastqc_trimmed_pe:
    input: os.path.join(TRIMMED_READS_DIR, "{sample}_trimmed_R{read_num}.fastq.gz")
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_trimmed_R{read_num}_fastqc.zip')
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_trimmed_R{read_num}.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

rule fastqc_primer_trimmed:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted_primer-trimmed_sorted.bam')
    output:
        html = os.path.join(FASTQC_DIR, '{sample}', '{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.html'),
        zip = os.path.join(FASTQC_DIR, '{sample}', '{sample}_aligned_sorted_primer-trimmed_sorted_fastqc.zip'),
    log: os.path.join(LOG_DIR, 'fastqc_{sample}_aligned_primer-trimmed.log')
    params:
        output_dir = os.path.join(FASTQC_DIR, '{sample}')
    shell: "{FASTQC_EXEC} -o {params.output_dir} {input} >> {log} 2>&1"

# TODO think about adding a global version to include all samples
# ANNOT: Generate overall QC report from all the per sample QC reports already
# generated in previous rules. Generates a subdir per sample, with a static
# structure determined by the multiqc executable.
rule multiqc:
  input: multiqc_input
  output: os.path.join(MULTIQC_DIR, '{sample}', 'multiqc_report.html')
  params:
    output_dir = os.path.join(MULTIQC_DIR, '{sample}')
  log: os.path.join(LOG_DIR, 'multiqc_{sample}.log')
  shell: "{MULTIQC_EXEC} -f -o {params.output_dir} {input} >> {log} 2>&1"


# TODO it should be possible to add customized parameter
rule lofreq:
    input:
        unpack(lofreq_input)
    output: vcf = os.path.join(VARIANTS_DIR, '{sample}_snv.vcf')
    log: os.path.join(LOG_DIR, 'lofreq_{sample}.log')
    run:
        call = (f"{LOFREQ_EXEC} call "
                f"-f {input.ref} "
                f"-o {output} "
                f"--verbose " 
                f"{input.aligned_bam} "
                f">> {log} 2>&1")

        shell(f"echo {call} > {log}")
        shell(call)
        
        # WIP create a dummy entry if no variant is found - use this as long as
        # the input-function solution doesn't work
        no_variant_vep(wildcards.sample, output.vcf)


rule vcf2csv:
    input: vcf2csv_input
    output: os.path.join(VARIANTS_DIR, '{sample}_snv.csv')
    params:
        script = os.path.join(SCRIPTS_DIR, 'vcfTocsv.py')
    log: os.path.join(LOG_DIR, 'vcf2csv_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} {output} >> {log} 2>&1"

rule vep:
    input:
        unpack(vep_input),
    output:
        os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2.txt"),
    params:
        buffer_size=VEP_BUFFER_SIZE,
        species=SPECIES,
        transcript_distance=VEP_TRANSCRIPT_DISTANCE,
        db_version=VEP_DB_VERSION
    log:
        os.path.join(LOG_DIR, "vep_{sample}.log"),
    shell:
        """
        {VEP_EXEC} --verbose --offline \
        --dir_cache {input.db_dir} \
        --DB_VERSION {params.db_version} \
        --buffer_size {params.buffer_size} \
        --species {params.species} \
        --check_existing \
        --distance {params.transcript_distance} \
        --biotype \
        --protein \
        --symbol \
        --transcript_version \
        --input_file {input.input_file} \
        --output_file {output} \
        --force_overwrite \
        >> {log} 2>&1
        """

rule parse_vep:
    input: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2.txt')
    output: os.path.join(VARIANTS_DIR, '{sample}_vep_sarscov2_parsed.txt')
    params:
        script = os.path.join(SCRIPTS_DIR, 'parse_vep.py')
    log: os.path.join(LOG_DIR, 'parse_vep_{sample}.log')
    shell: "{PYTHON_EXEC} {params.script} {input} {output} >> {log} 2>&1"


rule bam2fastq:
    input: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.bam')
    output: os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq')
    log: os.path.join(LOG_DIR, 'bam2fastq_{sample}.log')
    shell: "{SAMTOOLS_EXEC} fastq {input} > {output} 2>> {log} 3>&2"


rule kraken:
    input:
        unaligned_fastq = os.path.join(MAPPED_READS_DIR, '{sample}_unaligned.fastq'),
        database = KRAKEN_DB
    output: os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt')
    log: os.path.join(LOG_DIR, 'kraken_{sample}.log')
    shell: "{KRAKEN2_EXEC} --report {output} --db {input.database} {input.unaligned_fastq} >> {log} 2>&1"


rule krona_report:
    input:
        kraken_output = os.path.join(KRAKEN_DIR, '{sample}_classified_unaligned_reads.txt'),
        database = KRONA_DB
    output: os.path.join(REPORT_DIR, '{sample}.Krona_report.html')
    log: os.path.join(LOG_DIR, 'krona_report_{sample}.log')
    shell: "{IMPORT_TAXONOMY_EXEC} -m 3 -t 5 {input.kraken_output} -tax {input.database} -o {output} >> {log} 2>&1"


rule samtools_bedcov:
    input:
        mutations_bed = MUTATIONS_BED,
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_mut_cov.tsv')
    log: os.path.join(LOG_DIR, 'samtools_bedcov_{sample}.log')
    shell: "{SAMTOOLS_EXEC} bedcov {input.mutations_bed} {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule samtools_coverage:
    input:
        aligned_bam = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bam'),
        aligned_bai = os.path.join(MAPPED_READS_DIR, '{sample}_aligned_sorted.bai')
    output: os.path.join(COVERAGE_DIR, '{sample}_genome_cov.tsv')
    log: os.path.join(LOG_DIR, 'samtools_coverage_{sample}.log')
    shell: "{SAMTOOLS_EXEC} coverage {input.aligned_bam} > {output} 2>> {log} 3>&2"


rule create_sample_quality_table:
    input:
        genome_cov_file=os.path.join(COVERAGE_DIR, "{sample}_genome_cov.tsv"),
        mut_cov_file=os.path.join(COVERAGE_DIR, "{sample}_mut_cov.tsv"),
    output:
        os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
    params:
        script=os.path.join(SCRIPTS_DIR, "create_sample_quality_table.R"),
    log:
        os.path.join(LOG_DIR, "create_sample_quality_table_{sample}.log"),
    shell:
        """
        {RSCRIPT_EXEC} {params.script} \
            {wildcards.sample} \
            {input.genome_cov_file} \
            {input.mut_cov_file} \
            {output} > {log} 2>&1
        """


rule create_sample_quality_summary:
    input:
        script=os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files=expand(
            os.path.join(COVERAGE_DIR, "{sample}_quality.csv"),
            sample=SAMPLES,
        ),
    output:
        os.path.join(COVERAGE_DIR, "sample_quality_summary.csv"),
    log:
        os.path.join(LOG_DIR, "create_sample_quality_summary.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """


rule generate_navbar:
    input:
      script = os.path.join(SCRIPTS_DIR, "generateNavigation.R")
    output:
      os.path.join(REPORT_DIR, "_navbar.html")
    params:
      report_scripts_dir = os.path.join(SCRIPTS_DIR, "report_scripts")
    log: os.path.join(LOG_DIR, "generate_navigation.log")
    shell: "{RSCRIPT_EXEC} {input.script} \
{params.report_scripts_dir} {SAMPLE_SHEET_CSV} {output} > {log} 2>&1"


rule render_kraken2_report:
    input:
      script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
      report=os.path.join(SCRIPTS_DIR, "report_scripts", "taxonomic_classification.Rmd"),
      header=os.path.join(REPORT_DIR, "_navbar.html"),
      kraken=os.path.join(KRAKEN_DIR, "{sample}_classified_unaligned_reads.txt"),
      krona=os.path.join(REPORT_DIR, "{sample}.Krona_report.html")
    output: os.path.join(REPORT_DIR, "{sample}.taxonomic_classification.html")
    log: os.path.join(LOG_DIR, "reports", "{sample}_taxonomic_classification.log")
    shell: """{RSCRIPT_EXEC} {input.script} \
{input.report} {output} {input.header} \
'{{\
  "sample_name": "{wildcards.sample}",  \
  "site_dir":    "{REPORT_DIR}",        \
  "krona_file":  "{input.krona}",       \
  "kraken_file": "{input.kraken}",      \
  "logo": "{LOGO}" \
}}' > {log} 2>&1"""


rule run_deconvolution:
    input:
        script=os.path.join(SCRIPTS_DIR, "deconvolution.R"),
        deconvolution_functions=os.path.join(
            SCRIPTS_DIR, "deconvolution_funs.R"
        ),
        vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
        snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
    output:
        sigmut_df=os.path.join(MUTATIONS_DIR, "{sample}_sigmuts.csv"),
        non_sigmut_df=os.path.join(MUTATIONS_DIR, "{sample}_non_sigmuts.csv"),
        variant_proportions=os.path.join(
            VARIANTS_DIR, "{sample}_variant_abundance.csv"
        ),
        variants_with_meta=os.path.join(VARIANTS_DIR,
         "{sample}_variants_with_meta.csv"),
        mutations=os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"),
    log:
        os.path.join(LOG_DIR, "{sample}_deconvolution.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
        "{input.deconvolution_functions}" \
        "{wildcards.sample}" \
        "{MUTATION_SHEET_CSV}" \
        "{SAMPLE_SHEET_CSV}" \
        "{input.vep}" \
        "{input.snv}" \
        "{MUTATION_DEPTH_THRESHOLD}" \
        "{output.sigmut_df}" \
        "{output.non_sigmut_df}" \
        "{output.variant_proportions}" \
        "{output.variants_with_meta}" \
        "{output.mutations}" \
        "{DECONVOLUTION_METHOD}" \
        > {log} 2>&1
        """


rule render_variant_report:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "variantreport_p_sample.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        sigmut_file=os.path.join(MUTATIONS_DIR, "{sample}_sigmuts.csv"),
        non_sigmut_file=os.path.join(MUTATIONS_DIR, "{sample}_non_sigmuts.csv"),
        variant_abundance_file=os.path.join(
            VARIANTS_DIR, "{sample}_variant_abundance.csv"
        ),
        mutations=os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"),
        vep_raw=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2.txt"),
        vep=os.path.join(VARIANTS_DIR, "{sample}_vep_sarscov2_parsed.txt"),
        snv=os.path.join(VARIANTS_DIR, "{sample}_snv.csv"),
    output:
        varreport=os.path.join(REPORT_DIR, "{sample}.variantreport_p_sample.html"),
    params:
        vep_transcript_distance=VEP_TRANSCRIPT_DISTANCE,
        vep_species=SPECIES
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_variant_report.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
        {input.report} {output.varreport} {input.header} \
        '{{ \
          "sample_name": "{wildcards.sample}", \
          "sigmut_file": "{input.sigmut_file}", \
          "non_sigmut_file": "{input.non_sigmut_file}", \
          "variant_abundance_file": "{input.variant_abundance_file}", \
          "snv_file": "{input.snv}", \
          "vep_file": "{input.vep}", \
          "vep_file_raw": "{input.vep_raw}", \
          "vep_transcript_distance": "{params.vep_transcript_distance}", \
          "vep_species": "{params.vep_species}", \
          "logo": "{LOGO}" \
        }}' > {log} 2>&1
        """


# ANNOT: Render quality control report, summarizing coverage, and linking to the
# per sample fastq reports for the different stages of processing (see rule 
# multiqc)
rule render_qc_report:
    input:
        unpack(render_qc_report_input)
    output:
        html_report=os.path.join(REPORT_DIR, "{sample}.qc_report_per_sample.html"),
    params:
        render_qc_report_params,
    log:
        os.path.join(LOG_DIR, "reports", "{sample}_qc_report.log"),
    script:
        "snakefile_scripts/rule_render_qc_report.py"


rule create_variants_summary:
    input:
        script = os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files = expand(os.path.join(VARIANTS_DIR, "{sample}_variants_with_meta.csv"), sample = SAMPLES)
    output: os.path.join(VARIANTS_DIR, 'data_variant_plot.csv')
    log: os.path.join(LOG_DIR, "create_variants_summary.log")
    shell: """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """

rule create_mutations_summary:
    input:
        script = os.path.join(SCRIPTS_DIR, "create_summary_table.R"),
        files = expand(os.path.join(MUTATIONS_DIR, "{sample}_mutations.csv"), sample = SAMPLES)
    output: os.path.join(MUTATIONS_DIR, 'data_mutation_plot.csv')
    log: os.path.join(LOG_DIR, "create_mutations_summary.log")
    shell: """
        {RSCRIPT_EXEC} {input.script} {output} {input.files} > {log} 2>&1
        """


# TODO integrate the output of fastp.json to get the number of raw and trimmed reads
rule create_overviewQC_table:
    input:
        script=os.path.join(SCRIPTS_DIR, "overview_QC_table.R"),
        quality_table_file=os.path.join(COVERAGE_DIR, "sample_quality_summary.csv"),
    output:
        os.path.join(OUTPUT_DIR, "overview_QC.csv"),
    log:
        os.path.join(LOG_DIR, "create_overviewQC_table.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
            {SAMPLE_SHEET_CSV} \
            {output} \
            {INPUT_DIR} \
            {TRIMMED_READS_DIR} \
            {MAPPED_READS_DIR} \
            {input.quality_table_file} > {log} 2>&1
        """


rule run_mutation_regression:
    input:
        script=os.path.join(SCRIPTS_DIR, "mutation_regression.R"),
        mutations_csv=os.path.join(MUTATIONS_DIR, "data_mutation_plot.csv"),
        overviewQC=os.path.join(OUTPUT_DIR, "overview_QC.csv"),
        fun_lm=os.path.join(SCRIPTS_DIR, "pred_mutation_increase.R"),
        fun_tbls=os.path.join(SCRIPTS_DIR, "table_extraction.R"),
        mutation_sheet=MUTATION_SHEET_CSV
    output:
        mut_count_outfile=os.path.join(OUTPUT_DIR, "mutations_counts.csv"),
        unfilt_mutation_sig_outfile=os.path.join(
            OUTPUT_DIR, "unfiltered_mutations_sig.csv"
        ),
    log:
        os.path.join(LOG_DIR, "mutation_regression.log"),
    shell:
        """
        {RSCRIPT_EXEC} {input.script} \
            {input.mutations_csv} \
            {input.mutation_sheet} \
            {input.fun_lm} \
            {input.fun_tbls} \
            {MUTATION_COVERAGE_THRESHOLD} \
            {input.overviewQC} \
            {output.mut_count_outfile} \
            {output.unfilt_mutation_sig_outfile} \
            > {log} 2>&1
        """


rule render_index:
    input:
        script=os.path.join(SCRIPTS_DIR, "renderReport.R"),
        report=os.path.join(SCRIPTS_DIR, "report_scripts", "index.Rmd"),
        header=os.path.join(REPORT_DIR, "_navbar.html"),
        variants=os.path.join(VARIANTS_DIR, "data_variant_plot.csv"),
        mutations=os.path.join(MUTATIONS_DIR, "data_mutation_plot.csv"),
        overviewQC=os.path.join(OUTPUT_DIR, "overview_QC.csv"),
        mut_count_file=os.path.join(OUTPUT_DIR, "mutations_counts.csv"),
        unfiltered_mutation_sig_file=os.path.join(
            OUTPUT_DIR, "unfiltered_mutations_sig.csv"
        ),
    params:
        fun_tbls=os.path.join(SCRIPTS_DIR, "table_extraction.R"),
        fun_pool=os.path.join(SCRIPTS_DIR, "pooling.R"),
        fun_index=os.path.join(SCRIPTS_DIR, "fun_index.R"),
    output:
        report=os.path.join(REPORT_DIR, "index.html"),
    log:
        os.path.join(LOG_DIR, "reports", "index.log"),
    shell:
        """{RSCRIPT_EXEC} {input.script} \
        {input.report} {output.report} {input.header}   \
        '{{ \
          "variants_csv": "{input.variants}", \
          "mutations_csv": "{input.mutations}", \
          "sample_sheet": "{SAMPLE_SHEET_CSV}", \
          "mutation_sheet": "{MUTATION_SHEET_CSV}", \
          "mutation_coverage_threshold": "{MUTATION_COVERAGE_THRESHOLD}", \
          "mut_count_file": "{input.mut_count_file}", \
          "unfiltered_mutation_sig_file": "{input.unfiltered_mutation_sig_file}", \
          "logo": "{LOGO}", \
          "fun_tbls": "{params.fun_tbls}", \
          "fun_pool": "{params.fun_pool}", \
          "fun_index": "{params.fun_index}", \
          "overviewQC": "{input.overviewQC}", \
          "output_dir": "{OUTPUT_DIR}" \
        }}' > {log} 2>&1"""
