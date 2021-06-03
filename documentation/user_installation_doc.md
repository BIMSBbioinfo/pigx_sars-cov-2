# Introduction 
PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater samples and identifying given variants-of-concern of SARS-CoV-2. Currently wastewater samples are used, which are enriched for SARS-CoV-2. The pipeline can be used for continuous sampling. The output of the PiGx SARS-CoV-2 pipeline is summarized in a report which provides an intuitive visual overview about the development of variant abundance over time and location. Additionally there will be more detailed reports per sample, which cover the quality control of the samples, the detected variants and a taxonomic classification of all reads which are not aligned to SARS-CoV-2.

## Workflow

First the raw reads are trimmed by using [Prinseq](http://prinseq.sourceforge.net/) to improve alignment rates and mutation calling. Next the trimmed reads are aligned to the reference genome of SARS-CoV-2 using [BWA](https://github.com/lh3/bwa), the results are *SAM*/*BAM* files of **aligned** and **unaligned** reads. Following the alignment a quality check on raw and processed reads is performed by using [MultiQC](https://multiqc.info/). The next step is calling the variants and inferring SNVs (single nucleotide polymorphisms) on all **aligned** reads with [LoFreg](https://csb5.github.io/lofreq/). [ ... mutation step ?]
To check the wastewater samples also for the abundance of other species the **unaligned** reads will be taxonomicly classified with [Kraken2](https://github.com/DerrickWood/kraken2). The Kraken2 requires a database of all genomes the reads are getting aligned against, therefore keep in mind that you can only find those species which are included in the chosen database, for documentation how to set this up, see: [Prepare databases](#Prepare databases). For a better and interactive visualization of all species present in the wastewater [Krona](https://github.com/marbl/Krona/wiki) is used. Also here a small step of setting up a database is needed before running the pipeline. 

## Output
* Merged Report including:
* Overview of development of variant and mutation abundance over time and locations
* Quality Control report of raw and processed (trimmed) reads
* Variant report
* Taxonomic classification
* SAM/ BAM files of the aligned and unaligned reads against SARS-CoV-2


# Installation

## Install via guix

This pipeline is going to be packaged in the [GNU Guix](https://gnu.org/s/guix) package manager soon. Right now in this preliminary state it has to be installed from source manually.

## Install from source

First, please clone this repository and change directory accordingly:

```sh
git clone https://github.com/BIMSBbioinfo/pigx_sarscov2_ww.git
cd pigx_sarscov2_ww
```

To fetch code that is common to all PiGx pipelines run this:

```sh
git submodule update --init
```

Before setting everything up, though, make sure all dependencies are met. For this either install the following software in a directory listed in the `PATH` environment variable, or enter the provided reproducible Guix environment. If you are using Guix we definitely recommend the latter. This command spawns a sub-shell in which all dependencies are available:

```sh
USE_GUIX_INFERIOR=t guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

To use your current Guix channels instead of the fixed set of
channels, just omit the `USE_GUIX_INFERIOR` shell variable:

```sh
guix environment --pure -m manifest.scm --preserve=GUIX_LOCPATH
```

<details>
<summary>Software dependencies</summary>

- R
    - minimal
    - base64url
    - dplyr
    - dt
    - ggplot2
    - magrittr
    - plotly
    - qpcr
    - rmarkdown
    - stringr
    - tidyr
    - reshape2
- python
    - wrapper
    - pyyaml
- bash-minimal
- bwa
- ensembl-vep
- fastqc
- multiqc
- kraken2
- krona-tools
- lofreq
- prinseq
- samtools
- snakemake
- autoconf
- automake
- coreutils
- gawk
- grep
- make
- sed


</details>
</br>

Inside the environment you can then perform the usual build steps:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

# Getting started

At this point you are able to run the pipeline from within the current directory `pigx_sarscov2_ww`. Use `--help` to see the available options.

```sh
$ ./pigx-sars-cov2-ww --help
```
<details>
<summary>toggle output</summary>

```sh
usage: pigx-sars-cov2-ww [-h] [-v] (--init [{settings,sample-sheet,both}] | -s SETTINGS) [-c CONFIGFILE]
                         [--target TARGET] [-n] [--graph GRAPH] [--force] [--reason] [--unlock] [--verbose]
                         [--printshellcmds]
                         [sample_sheet]

PiGx SARS-CoV-2 wastewater sequencing Pipeline.

This is a pipeline for analyzing data from sequenced wastewater
samples and identifying given variants-of-concern of SARS-CoV-2.  The
pipeline can be used for continuous sampling.  The output report
provides an intuitive visual overview about the development of variant
abundance over time and location.

positional arguments:
  sample_sheet                            The sample sheet containing sample data in CSV format.

optional arguments:
  -h, --help                              show this help message and exit
  -v, --version                           show program's version number and exit
  --init [{settings,sample-sheet,both}]   Generate a template SETTINGS file, a SAMPLE-SHEET.  Leave empty for both.
  -s SETTINGS, --settings SETTINGS        A YAML file for settings that deviate from the defaults.
  -c CONFIGFILE, --configfile CONFIGFILE  The config file used for calling the underlying snakemake process.  By
                                          default the file 'config.json' is dynamically created from the sample
                                          sheet and the settings file.
  --target TARGET                         Stop when the named target is completed instead of running the whole
                                          pipeline.  The default target is "final-report".  Pass "--target=help"
                                          to describe all available targets.
  -n, --dry-run                           Only show what work would be performed.  Do not actually run the
                                          pipeline.
  --graph GRAPH                           Output a graph in PDF format showing the relations between rules of
                                          this pipeline.  You must specify a graph file name such as
                                          "graph.pdf".
  --force                                 Force the execution of rules, even though the outputs are considered
                                          fresh.
  --reason                                Print the reason why a rule is executed.
  --unlock                                Recover after a snakemake crash.
  --verbose                               Print supplementary info on job execution.
  --printshellcmds                        Print commands being executed by snakemake.

This pipeline was developed by the Akalin group at MDC in Berlin in 2017-2021.
```
</details>
</br>

Though, to actually use it on your experimental data still more setup is required. Please follow [the steps to prepare the required databases](##Prepare-databases) first. Then afterwards, you can run the pipeline from the source directory.

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww [options] sample_sheet.csv
```

For example, with this command the pipeline is used for the included test data: 

```sh
PIGX_UNINSTALLED=t ./pigx-sars-cov2-ww -s tests/settings.yaml tests/sample_sheet.csv
```

These links guide you to more information on the [settings file](##settings-file) and the [sample sheet](##sample-sheet).

## Structure overview

This guide will lead to the following structure of directories and files.

```
tests/
|
|__sample_sheet.csv
|
|__settings.yaml
|
|__sample_data/__reads/
|
|__databases/
   │
   │__kraken_db/
   │
   │__krona_db/
   │
   │__vep_db/
   |
   |__sigmut_db/
``` 


## Settings file and sample sheet

A draft for these two files are given in the `tests/` directory. The pre-filled rows should work but feel free to modify either locations in the settings file or given samples in the sample sheet.

## Dummy sample fastq file

For testing purposes we created small dummy reads `Test_R1.fastq` and `Test_R2.fastq` which can be found in `tests/sample_data/reads/`.

## Prepare databases

Before the pipeline can be run, 3 databases must be downloaded and their location will need to be provided in the settings file. Depending on the size of the databses this can take some time.
Be sure that the pigx-sarscov2-ww pipeline is downloaded and the tools are installed or used via the provided and suggested guix environment, bevore preparing the databases. One database (signature mutations, `sigmut_db`) is already provided via the repository. The folder structure is suggested like the following and pre-filled accordingly in the settings file.

### Kraken2 database

There are several libraries of genomes that can be used to classify the (unaligned) reads. It is up to you which one to use, but be sure that they fulfill the necessities stated by Kraken2 [Kraken2 manual](https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases). We recommend to use the Plus-PFP library provided [here](https://benlangmead.github.io/aws-indexes/k2).
It is also possible to have multiple Kraken2 databases, just be sure to provide the wanted one to the settings file.

First download and unpack the database in the `tests/databases/kraken_db/`:

```
DIR=tests/databases/kraken_db/
mkdir -p $DIR

# NOTE: This command will download a very large file, after unpacking this will
# require about 100GB of disc space. If this is not feasible use another
# database instead. For this please see link above or commented lines below.
wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20210127.tar.gz | tar -C $DIR -xzv

# Use the following two lines to use smaller 8GB version instead.
# wget -qO- https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_8gb_20210127.tar.gz | tar -C $DIR -xzv
```

Next go to the `tests/databases/` directory and build the Kraken database. This might take a while (depends on the size of the downloaded database):

```
cd tests/databases
DBNAME=kraken_db/
kraken2-build --use-ftp --download-taxonomy --db $DBNAME # if this fails, you might want to try it without the --use-ftp flag
kraken2-build --build --db $DBNAME
```

Kraken might tell you that it can't find a library subdirectory. If that's the case it should be fine though. 

### Krona database

Krona Tools needs a two files, which have to be installed in the `tests/databases/krona_db/` folder. Also this might take a while:

```
DBNAME=tests/databases/krona_db/
mkdir -p $DBNAME
KRONA=$(dirname $(which ktImportTaxonomy))/../share/krona-tools/ # this is just a workaround until tool paths are declared
$KRONA/updateTaxonomy.sh $DBNAME # the scripts are stored a priori in that folder
$KRONA/updateAccessions.sh $DBNAME
```

### VEP database

Just download the `sars_cov_2` database for VEP and unpack it in the `tests/databases/vep_db/` directory.

```
DBNAME=tests/databases/vep_db/
mkdir -p $DBNAME
wget -qO- ftp://ftp.ensemblgenomes.org/pub/viruses/variation/indexed_vep_cache/sars_cov_2_vep_101_ASM985889v3.tar.gz | tar -C $DBNAME -xzv
```

### sigmut database

Nothing to be done here. Necessary files are provided in `tests/databases/sigmut_db/`.
