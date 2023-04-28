
# argparsing -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# give default parameters
# FIXME: Order these args sensibly, needs also to be adjusted in snakefile
if (length(args) == 0) {
  args <- c(
    mutations_csv = "",
    mutation_sheet = "",
    fun_lm = "",
    fun_tbls = "",
    mutation_coverage_threshold = "",
    overviewQC = "",
    mut_count_outfile = "",
    unfilt_mutation_sig_outfile = ""
  )
}

# names must match order in snakefile and defaults
arg_names <- c(
  "mutations_csv",
  "mutation_sheet",
  "fun_lm",
  "fun_tbls",
  "mutation_coverage_threshold",
  "overviewQC",
  "mut_count_outfile",
  "unfilt_mutation_sig_outfile"
  )

params <- lapply(args, function(x) x)
names(params) <- arg_names

# pretty print parameters for easier debugging
cat("Script running with parameters:\n\n")
par_vec <- c()
for (i in seq_along(params)) {
  par_vec[i] <- paste0(names(params)[i], " = \"", params[[i]], "\"")
}
cat(paste(par_vec, collapse = ",\n"))
cat("\n\n")


## ----libraries----------------------------------------------------------------
library(dplyr)
library(data.table)

## ----input--------------------------------------------------------------------
df_mut <- fread(params$mutations_csv)

quality_df <- fread(params$overviewQC)


## ----filter_plot_frames_samplescore, warning=TRUE-----------------------------

mutation_coverage_threshold <- params$mutation_coverage_threshold %>%
  # check if value is given as fraction [0,1] or percentage [0,100]
  ifelse(. >= 0 & . <= 1,
         . * 100,
         .
  ) %>%
  as.numeric()

# FIXME: Check if all this computation is necessary for the tasks below
good_samples_df <- quality_df %>%
  filter(as.numeric(perc_muts_covered) >= mutation_coverage_threshold)

approved_mut_plot <- df_mut %>%
  filter(samplename %in% good_samples_df$samplename)

# pool the samples per day, discard locations
weights <- quality_df %>%
  dplyr::select(c(samplename, total_reads))


## ----linear_regression, eval=RUN_MUTATION_REGRESSION--------------------------
# only run regression if there are values after filtering and at least two dates
# are left over
if (nrow(approved_mut_plot) > 0 &&
  length(unique(approved_mut_plot$dates)) > 1) {
  source(params$fun_lm)

  mutation_sheet <- params$mutation_sheet

  sigmuts_df <- fread(mutation_sheet) %>%
    mutate(across(everything(), ~dplyr::na_if(.x, ""))) %>%
    # split gene name of for easier matching
    mutate_all(funs(str_replace(., "^[^:]*:", "")))

  # create vector of metadata col names to be excluded
  meta_cols_excl <- c(
    "samplename",
    "location_name",
    "coordinates_lat",
    "coordinates_long"
  )

  changing_muts <- approved_mut_plot %>%

    # enforcing date type for column dates ...it will get rid of the time
    mutate(dates = as.Date(dates)) %>%

    select(- all_of(meta_cols_excl)) %>%

    # remove mutations with NA for all rows and create a new dataframe
    select(where(function(x) any(!is.na(x))))


  mutations_sig_unfiltered <- refined_lm_model(changing_muts)
  mutations_sig <- mutations_sig_unfiltered %>%
    filter(pvalues >= 0.05) %>%
    arrange(desc(coefficients))

  ## ----mutation_counts--------------------------------------------------------
  # get functions for counting and writing
  # TODO: check where the fun_tbls script is needed
  source(params$fun_tbls)
  count_frame <- write_mutations_count(df_mut, sigmuts_df, mutations_sig)

  # write to file
  fwrite(count_frame,
    file.path(
      params$mut_count_outfile
    )
  )

  # mutations_sig
  fwrite(mutations_sig_unfiltered,
    file.path(
      params$unfilt_mutation_sig_outfile
    )
  )
} else {
  # write empty files
  err_msg <- "No significantly increasing mutations found..."
  writeLines(err_msg, params$mut_count_outfile)
  writeLines(err_msg, params$unfilt_mutation_sig_outfile)
}
