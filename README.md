<a name="logo"/>
<div align="center">
<img src="images/Logo_PiGx.png" alt="PiGx Logo"  width="30%" height="30%" ></img>
</a>
</div>

# PiGx SARS-CoV-2 Wastewater Sequencing Pipeline

**Copyright 2021-2022: Vic-Fabienne Schumann, Ricardo Wurmus, Miriam Faxel, Jan
Dohmen, Rafael Cuadrat, Bora Uyar, Vedran Franke, Alexander Blume, Jonas
Freimuth, and Altuna Akalin.**  
**This work is distributed under the terms of the GNU General Public
License, version 3 or later. It is free to use for all purposes.**

-----------

PiGx SARS-CoV-2 is a pipeline for analysing data from sequenced wastewater
samples and identifying given lineages of SARS-CoV-2. It was developed for
wastewater samples, which are enriched for SARS-CoV-2. The pipeline can be used
for continuous sampling. The output of the PiGx SARS-CoV-2 pipeline is
summarized in a report which provides an intuitive visual overview about the
development of variant abundance over time and location (see our
[example report](https://github.com/BIMSBbioinfo/pigx_sars-cov-2#sample-reports)).
Additionally there will be more detailed reports per sample, which cover the
quality control of the samples, the detected variants and a taxonomic
classification of all unaligned reads. This version of the pipeline was designed
to work with single- and paired-end amplicon sequencing data generated using the
ARTIC nCoV-2019 primers. But it works with any set of primer.

Features to enable e.g. single-end input are currently under development. We are
happy to hear about any other feature request! Please consider using the Issue
Tracker for this.

# Publication

The publication that describes and uses this pipeline is as follows:

_SARS-CoV-2 infection dynamics revealed by wastewater sequencing analysis and
deconvolution_

Vic-Fabienne Schumann, Rafael Ricardo de Castro Cuadrat, Emanuel Wyler, Ricardo
Wurmus, Aylina Deter, Claudia Quedenau, Jan Dohmen, Miriam Faxel, Tatiana
Borodina, Alexander Blume, Jonas Freimuth, Martin Meixner, José Horacio Grau,
Karsten Liere, Thomas Hackenbeck, Frederik Zietzschmann, Regina Gnirss, Uta
Böckelmann, Bora Uyar, Vedran Franke, Niclas Barke, Janine Altmüller, Nikolaus
Rajewsky, Markus Landthaler, Altuna Akalin

Science of The Total Environment, 2022; doi:
<https://doi.org/10.1016/j.scitotenv.2022.158931>

## Reproducing the analysis

The presented analysis results in the publication were produced using
three different development versions of PiGx SARS-CoV-2.  The commits
are as follows:

- dataset-Berlin250, dataset-NYC(RBD) (MiSeq data and all samples merged)
  commit: 524ed4832a6972fd695c0eeec25264188710a143

- dataset-Berlin35, dataset-NYC(RBD) (iSeq data), insilico-simulation
  commit: 0a150c4bec58a5a8296c870586e225e49ee2b6f8

- UCSD-spike in
  commit: bd87e7f2d83317e9d83f6fd81abb631af95476f6

To enter a faithful reproduction of the development environment used
for the analysis at these three commits execute this command:

```sh
guix time-machine -C reproducibility/channels.scm -- shell --pure -m reproducibility/manifest.scm
```

This will build the exact version of Guix described in
`reproducibility/channels.scm` and then use that version of Guix to
launch a shell session where the environment specified in
`reproducibility/manifest.scm` is activated.  All other environment
variables are cleared, so consider augmenting this environment with
`git` and other tools as needed.

Proceed to check out the repository commit in question (see above) and
build the pipeline in the usual way:

```sh
/bin/git checkout $commit
git submodule update --init
./bootstrap.sh
./configure --prefix=$HOME/pigx-sars-cov-2--reproduce
make
make install
```

You can then launch the pipeline with

```sh
$HOME/pigx-sars-cov-2--reproduce/bin/pigx-sars-cov-2 ...
```


# Installation

## Databases

Some tools require some databases to be installed locally. For this please
either see the Documentation below or run
[download_databases.sh](https://github.com/BIMSBbioinfo/pigx_sars-cov-2/blob/main/scripts/download_databases.sh.in)
after the pipeline is installed.

## Installation via Guix (recommended)

You can install this pipeline with all its dependencies using GNU Guix:

```sh
guix install pigx-sars-cov-2
```

Using GNU Guix has many advantages in terms of reproducibility of your projects
(see e.g.
[here](https://academic.oup.com/gigascience/article/7/12/giy123/5114263)) If you
don't have GNU Guix on your system you can also consider using
[GNU Guix in a VM](https://guix.gnu.org/manual/en/html_node/Running-Guix-in-a-VM.html)

## Installation from source

You can also install PiGx-SARS-CoV-2 from source manually. Make sure that all
the required dependencies (for this see e.g
[configure.ac](https://github.com/BIMSBbioinfo/pigx_sars-cov-2/blob/main/configure.ac))
are installed e.g. by installing them through a package manager like Guix or
Conda. We recommend using Guix to install the complete pipeline (see above).

<details>
    <summary> The following tools must be available (latest version): </summary>

- snakemake
- samtools
- bwa
- bedtools
- fastp
- fastqc
- R
- Rscript
- kraken2
- kraken2-build
- ktImportKrona
- ktImportTaxonomy
- ivar
- lofreq
- vep
- multiqc
- pandoc

And the R-packages:  

- DT
- base64url
- deconvr
- dplyr
- ggplot2
- htmltools
- jsonlite
- magrittr
- mass
- plotly
- qpcR
- reshape2
- R.utils
- rmarkdown
- stringr
- tidyr
- viridis

All of these dependencies must be present in the environment at
configuration time.
</details>

You will also need to fetch code that is common to all PiGx pipelines by running
this:

```sh
git submodule update --init
```

Then follow the usual build steps of the GNU build system:

```sh
./bootstrap.sh # to generate the "configure" script
./configure
make
make check
```

## Use the Docker Image

Will be available soon

# Documentation

Please see our documentation in order to find out about more details e.g. about
the required structure of the input files:
<https://bioinformatics.mdc-berlin.de/pigx_docs/pigx-sars-cov-2.html>

# Sample Reports

See as example the HTML report produced by PiGx SARS-CoV-2 used for
["SARS-CoV-2 infection dynamics revealed by wastewater sequencing analysis and deconvolution"](https://doi.org/10.1016/j.scitotenv.2022.158931):
<https://bimsbstatic.mdc-berlin.de/akalin/AAkalin_pathogenomics/sarscov2_ww_reports/211104_pub_version/index.html>

# Getting help

If you have any questions please e-mail: pigx@googlegroups.com or use the web
form to ask questions <https://groups.google.com/forum/#!forum/pigx/>.

If you run into any bugs, please open an issue here:
<https://github.com/BIMSBbioinfo/pigx_sars-cov-2/issues>.

# Links

- [The Bioinformatics and Omics Data Science Platform](https://bioinformatics.mdc-berlin.de)
