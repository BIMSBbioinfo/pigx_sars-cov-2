;;; PiGx SARS-CoV2 wastewater sequencing pipeline
;;; Copyright © 2021, 2022 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx SARS-CoV2 wastewater sequencing pipeline
;;;
;;; This is free software; see LICENSE file for details.

(define %packages
  (list "bash-minimal"
        "bwa"
        "bedtools"
        "ensembl-vep"
        "fastp"
        "fastqc"
        "multiqc"
        "ivar"
        "kraken2"
        "krona-tools"
        "lofreq"
        "samtools"
        "snakemake"
        "r-minimal"
        "r-base64url"
        "r-dplyr"
        "r-dt"
        "r-ggplot2"
        "r-magrittr"
        "r-plotly"
        "r-qpcr"
        "r-rmarkdown"
        "r-stringr"
        "r-tidyr"
        "r-tidyverse"
        "r-reshape2"
        "r-r-utils"
        "r-viridislite"
        "r-viridis"
        "python-wrapper"
        "python-pyyaml"))

(define %native-packages
  (list "autoconf"
        "automake"
        "coreutils"
        "diffutils" ; for tests
        "findutils" ; for make distcheck
        "gawk"
        "grep"
        "gzip" ; for make distcheck
        "make"
        "sed"
        "tar"  ; for make distcheck and db downloads
        "wget")) ; for db downloads

(specifications->manifest
 (append %packages %native-packages))
