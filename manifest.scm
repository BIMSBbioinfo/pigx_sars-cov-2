;;; PiGx SARS-CoV-2 wastewater sequencing pipeline
;;; Copyright © 2021, 2022 Ricardo Wurmus <rekado@elephly.net>
;;;
;;; This file is part of PiGx SARS-CoV-2 wastewater sequencing pipeline
;;;
;;; This is free software; see LICENSE file for details.
;;;
;;; Run the following command to enter a development environment:
;;;
;;;  $ guix environment -m manifest.scm
;;;
;;; When the shell variable USE_GUIX_INFERIOR is set (to any value), a
;;; specific set of Guix channels will be used to build the
;;; environment.

(use-modules
 (guix profiles)
 (guix channels)
 (guix inferior)
 (ice-9 match))

(define channels
  (list (channel
          (name 'guix)
          (url "https://git.savannah.gnu.org/git/guix.git")
          (branch "master")
          (commit
            "e3ed1d09f9d490eff6becd6e9cb85a4d36c48e85")
          (introduction
            (make-channel-introduction
              "9edb3f66fd807b096b48283debdcddccfea34bad"
              (openpgp-fingerprint
                "BBB0 2DDF 2CEA F6A8 0D1D  E643 A2A0 6DF2 A33A 54FA"))))))

(define (lookup name)
  (specification->package name))

(define (lookup-inferior name)
  (define inferior
    (inferior-for-channels channels))
  (match (lookup-inferior-packages inferior name)
    ((first . rest) first)
    (_ (error
        (format #false "Could not find package `~a'.~%" name)))))

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
        "r-data-table"
        "r-deconvr"
        "r-dplyr"
        "r-dt"
        "r-ggplot2"
        "r-htmltools"
        "r-jsonlite"
        "r-knitr"
        "r-magrittr"
        "r-mass"
        "r-plotly"
        "r-qpcr"
        "r-reshape2"
        "r-rmarkdown"
        "r-r-utils"
        "r-stringr"
        "r-tidyr"
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

(packages->manifest
 (let ((how (if (getenv "USE_GUIX_INFERIOR")
                lookup-inferior lookup)))
   (map how (append %packages %native-packages))))
