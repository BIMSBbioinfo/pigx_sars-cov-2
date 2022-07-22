# argparsing -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) <= 0) {
  args <- c(
    "sample_name" = "",
    "genome_cov_file" = "",
    "mut_cov_file" = "",
    "quality_table_outfile" = ""
  )
}

arg_names <- c(
  "sample_name",
  "genome_cov_file",
  "mut_cov_file",
  "quality_table_outfile"
)

params <- lapply(args[seq_along(arg_names)], function(x) x)
names(params) <- arg_names

# pretty print parameters for easier debugging
cat("Script running with parameters:\n\n")
par_vec <- c()
for (i in seq_along(params)) {
  par_vec[i] <- paste0(names(params)[i], " = \"", params[[i]], "\"")
}
cat(paste(par_vec, collapse = ",\n"))
cat("\n\n")


library("stringr")
library("dplyr")
library("magrittr")
library("data.table")

genome_cov_df <- fread(params$genome_cov_file)

mut_cov_names <- c(
  "chrom",
  "start",
  "end",
  "name",
  "total_coverage"
)

mut_cov_df <- fread(params$mut_cov_file, col.names = mut_cov_names) %>%
  mutate(
    mut_loc = as.numeric(str_extract(name, "(?<=_)[0-9]+(?=_)")),
    region_length = end - start,
    avg_depth = total_coverage / region_length
  )


# summarize mut_cov_df
mut_cov_sum_df <- mut_cov_df %>%
  summarise(
    n_loc_total   = n(),
    n_loc_covered = sum(avg_depth > 0),
    n_loc_no_cov  = sum(avg_depth <= 0),
    avg_loc_cvrg  = mean(avg_depth)
  )

# join dfs
sample_qual_df <- bind_cols(genome_cov_df, mut_cov_sum_df) %>% 
  mutate(
    samplename = params$sample_name,
    perc_muts_covered = (n_loc_covered / n_loc_total) * 100
  )

# save df
fwrite(sample_qual_df, params$quality_table_outfile)
