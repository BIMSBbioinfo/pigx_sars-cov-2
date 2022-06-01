
# argparsing -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# give default parameters
# FIXME: Order these args sensibly, needs also to be adjusted in snakefile
if (length(args) == 0) {
  args <- c(
    mutations_csv = "",
    coverage_dir = "",
    mutation_sheet = "",
    fun_cvrg_scr = "",
    fun_lm = "",
    fun_pool = "",
    mutation_coverage_threshold = "",
    overviewQC = "",
    mut_count_outfile = "",
    unfilt_mutation_sig_outfile = ""
  )
}

params <- list(
  mutations_csv = args[[1]],
  coverage_dir = args[[2]],
  mutation_sheet = args[[3]],
  fun_cvrg_scr = args[[4]],
  fun_lm = args[[5]],
  fun_pool = args[[6]],
  mutation_coverage_threshold = args[[7]],
  overviewQC = args[[8]],
  mut_count_outfile = args[[9]],
  unfilt_mutation_sig_outfile = args[[10]]
)


## ----libraries----------------------------------------------------------------
library(dplyr)

## ----input--------------------------------------------------------------------
df_mut <- read.csv(params$mutations_csv,
  header = TRUE,
  check.names = FALSE
)


## ----filter_plot_frames_samplescore, warning=TRUE-----------------------------
source(params$fun_pool)
source(params$fun_cvrg_scr)

# FIXME: Check if all this computation is necessary for the tasks below
good_samples_df <- merge(get_genome_cov(params$coverage_dir),
  get_mutation_cov(params$coverage_dir),
  by = "samplename"
) %>%
  mutate(proport_muts_covered = round(
    (as.numeric(total_muts_cvrd) * 100) / as.numeric(total_num_muts), 1
  )) %>%
  filter(as.numeric(proport_muts_covered) >= params$mutation_coverage_threshold)

approved_mut_plot <- df_mut %>%
  filter(samplename %in% good_samples_df$samplename)

# pool the samples per day, discard locations
weights <- read.csv(params$overviewQC, header = TRUE, check.names = FALSE) %>%
  dplyr::select(c(samplename, total_reads))


## ----linear_regression, eval=RUN_MUTATION_REGRESSION--------------------------
# only run regression if there are values after filtering and at least two dates
# are left over
if (nrow(approved_mut_plot) > 0 &&
  length(unique(approved_mut_plot$dates)) > 1) {
  source(params$fun_lm)

  mutation_sheet <- params$mutation_sheet

  sigmuts_df <- read.csv(mutation_sheet, header = TRUE) %>%
    na_if("") %>%
    # split gene name of for easier matching
    mutate_all(funs(str_replace(., "^[^:]*:", "")))


  changing_muts <- parsing_mutation_plot_data(approved_mut_plot)
  mutations_sig_unfiltered <- refined_lm_model(changing_muts)
  mutations_sig <- filter_lm_res_top20(mutations_sig_unfiltered, 0.05)

  ## ----mutation_counts--------------------------------------------------------
  # get functions for counting and writing
  # TODO: check where the fun_tbls script is needed
  source(params$fun_tbls)
  count_frame <- write_mutations_count(df_mut, sigmuts_df, mutations_sig)

  # write to file
  write.csv(count_frame,
    file.path(
      params$mut_count_outfile
    ),
    na = "NA", row.names = FALSE, quote = FALSE
  )

  # mutations_sig
  write.csv(mutations_sig_unfiltered,
    file.path(
      params$unfilt_mutation_sig_outfile
    ),
    na = "NA", row.names = FALSE, quote = FALSE
  )

} else {
  # TODO: This will generate an error when read as csv
  # write empty files
  err_msg <- "Not significantly increasing mutations found..."
  writeLines(err_msg,  params$mut_count_outfile)
  writeLines(err_msg,  params$unfilt_mutation_sig_outfile)
}