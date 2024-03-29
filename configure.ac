dnl -*- Autoconf -*-

AC_INIT([PiGx_sars-cov-2], [m4_translit(m4_esyscmd([cat VERSION]),m4_newline)])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([pigx-common/common/m4])
AM_INIT_AUTOMAKE([gnu color-tests tar-pax -Wall -Wno-portability foreign])
AM_SILENT_RULES([yes])

dnl We need sed for version checks
AC_PROG_SED

dnl Load Python support for Automake and ensure we have a Python
dnl interpreter for at least version 3.5.  This sets the "PYTHON"
dnl variable.
dnl https://www.gnu.org/software/automake/manual/html_node/Python.html
AM_PATH_PYTHON([3.5])

AX_PYTHON_MODULE([yaml], "required")

dnl Check for required programmes and store their full path in the
dnl given variables.  The variables are used to substitute
dnl placeholders in the scripts.

dnl The variables only need to be set manually when they cannot be
dnl found via the PATH environment variable.
AC_DEFUN([find_or_override_prog],
[AC_ARG_VAR($1, override location of $2 executable)dnl
AS_IF([test -z "$$1"],dnl
      [AC_PATH_PROG([$1], [$2], no)dnl
      if test x"$$1" == x"no" ; then
        AC_MSG_ERROR([Could not find $2.])
      fi
      ],dnl
      [AC_MSG_NOTICE([Using $$1 as $2 executable.])])])

find_or_override_prog([GUNZIP],          [gunzip])
find_or_override_prog([GZIP],            [gzip])
find_or_override_prog([GNUBASH],         [bash])
find_or_override_prog([SNAKEMAKE],       [snakemake])
find_or_override_prog([SAMTOOLS],        [samtools])
find_or_override_prog([BWA],             [bwa])
find_or_override_prog([BEDTOOLS],        [bedtools])
find_or_override_prog([FASTP],           [fastp])
find_or_override_prog([FASTQC],          [fastqc])
find_or_override_prog([R],               [R])
find_or_override_prog([RSCRIPT],         [Rscript])
find_or_override_prog([KRAKEN2],         [kraken2])
find_or_override_prog([KRAKEN2_BUILD],   [kraken2-build])
find_or_override_prog([IMPORT_KRONA],    [ktImportKrona])
find_or_override_prog([IMPORT_TAXONOMY], [ktImportTaxonomy])
find_or_override_prog([IVAR],            [ivar trim])
find_or_override_prog([LOFREQ],          [lofreq])
find_or_override_prog([VEP],             [vep])
find_or_override_prog([MULTIQC],         [multiqc])
find_or_override_prog([WGET],            [wget])
find_or_override_prog([PANDOC],          [pandoc])

AC_ARG_ENABLE([r-packages-check],
  AS_HELP_STRING([--disable-r-packages-check], [Do not check any R packages.]),
  [AC_MSG_NOTICE([Skipping R packages check.  Be careful!])],
  [dnl

AX_R_PACKAGE([base64url])
AX_R_PACKAGE([data.table])
AX_R_PACKAGE([deconvR])
AX_R_PACKAGE([dplyr])
AX_R_PACKAGE([DT])
AX_R_PACKAGE([ggplot2])
AX_R_PACKAGE([htmltools])
AX_R_PACKAGE([jsonlite])
AX_R_PACKAGE([knitr])
AX_R_PACKAGE([magrittr])
AX_R_PACKAGE([MASS])
AX_R_PACKAGE([plotly])
AX_R_PACKAGE([qpcR])
AX_R_PACKAGE([reshape2])
AX_R_PACKAGE([rmarkdown])
AX_R_PACKAGE([R.utils])
AX_R_PACKAGE([stringr])
AX_R_PACKAGE([tidyr])
AX_R_PACKAGE([viridis])
  ])

AC_ARG_WITH([pigx-db-dir],
  AS_HELP_STRING([--with-pigx-db-dir],
    [Default database dir used by the pipeline. [default: "~/.local/share/pigx/databases]]),
  [PIGX_DB_DIR=$withval], [PIGX_DB_DIR="~/.local/share/pigx/databases"])

AC_ARG_VAR(PIGX_DB_DIR, Default database dir used by the pipeline.)

AC_ARG_ENABLE([environment-capture],
  AS_HELP_STRING([--disable-environment-capture], [Do not capture environment variables.]),
  [AC_MSG_NOTICE([Environment variables will not be captured.])],
  [dnl
AC_MSG_NOTICE([Environment variables will be captured.])
capture_environment=yes
AC_ARG_VAR(R_LIBS_SITE, The search path for all R packages)
AC_ARG_VAR(PYTHONPATH, The search path for all Python packages)
  ])

AC_SUBST([capture_environment])
AC_SUBST_FILE(PIPELINE_META)
PIPELINE_META=./META
AC_CONFIG_FILES([Makefile META etc/settings.yaml tests/setup_test_settings.yaml])
AC_CONFIG_FILES([tests/tests/vep.sh], [chmod +x tests/tests/vep.sh])
AC_CONFIG_FILES([test.sh], [chmod +x test.sh])
AC_CONFIG_FILES([scripts/download_databases.sh], [chmod +x scripts/download_databases.sh])
AC_CONFIG_FILES([qsub-template.sh:pigx-common/common/qsub-template.sh.in])
AC_CONFIG_FILES([pigx-sars-cov-2:pigx-common/common/pigx-runner.in], [chmod +x pigx-sars-cov-2])
AC_OUTPUT
