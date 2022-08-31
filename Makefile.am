ACLOCAL_AMFLAGS = -I pigx-common/common/m4

nodist_bin_SCRIPTS = \
  pigx-sars-cov-2

dist_pkglibexec_SCRIPTS = \
  snakefile.py

pkglibexec_scriptsdir = $(pkglibexecdir)/scripts
dist_pkglibexec_scripts_SCRIPTS = \
  scripts/create_summary_table.R \
  scripts/deconvolution.R \
  scripts/deconvolution_funs.R \
  scripts/download_databases.sh \
  scripts/fun_index.R \
  scripts/generateNavigation.R \
  scripts/get_qc_table.py \
  scripts/mutation_regression.R \
  scripts/overview_QC_table.R \
  scripts/parse_vep.py \
  scripts/pooling.R \
  scripts/pred_mutation_increase.R \
  scripts/renderReport.R \
  scripts/sample_coverage_score.R \
  scripts/table_extraction.R \
  scripts/vcfTocsv.py \
  scripts/vep_rule.py

pkglibexec_reportsdir = $(pkglibexecdir)/scripts/report_scripts/
dist_pkglibexec_reports_SCRIPTS = \
  scripts/report_scripts/index.Rmd \
  scripts/report_scripts/qc_report_per_sample.Rmd \
  scripts/report_scripts/taxonomic_classification.Rmd \
  scripts/report_scripts/variantreport_p_sample.Rmd

dist_pkgdata_DATA =            \
  etc/pretty.txt \
  etc/settings.yaml \
  etc/sample_sheet.csv.example \
  images/Logo_PiGx.png

dist_doc_DATA = 									\
  README.md

EXTRA_DIST = \
  VERSION \
  META.in \
  tests/sample_data/Test_vep_sarscov2_parsed.txt \
  tests/sample_data/test_vep_parsing_input.txt \
  tests/sample_sheet.csv \
  tests/settings.yaml \
  $(wildcard tests/reads/*) \
  tests/sample_data/NC_045512.2.fasta \
  tests/sample_data/nCoV-2019_NCref.bed \
  tests/sample_data/Covid_CG_NTmutation_t07.bed \
  tests/sample_data/mutation_sheet_211006_covidCG_NT_location.csv

AM_TESTS_ENVIRONMENT = srcdir="$(abs_top_srcdir)" builddir="$(abs_top_builddir)" PIGX_UNINSTALLED=1 PIGX_UGLY=1

TESTS = \
  tests/tests/vep.sh \
  test.sh

clean-local:
	-rm -rf .snakemake
	-rm -rf tests/output

CLEANFILES = $(nodist_bin_SCRIPTS) config.json

integration:
	$(AM_TESTS_ENVIRONMENT) $(abs_top_srcdir)/pigx-sars-cov-2 --printshellcmds -s tests/settings.yaml tests/sample_sheet.csv

# Build docker image with Guix
docker: dist
	guix pack -C none -e '(load "guix.scm")' --with-source=pigx_sars-cov-2-$(VERSION).tar.gz -f docker -S /bin=bin -S /lib=lib -S /share=share glibc-utf8-locales tzdata coreutils bash
