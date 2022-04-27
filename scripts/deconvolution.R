## ----setup, include = FALSE, warning = FALSE----------------------------------

library(knitr)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)
library(qpcR)
library(stringr)
library(magrittr)
library(base64url)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

# give default parameters
if (length(args) == 0) {
  args <- c(
    sample_name = "",
    output_dir = "",
    vep_file = "",
    snv_file = "",
    sample_sheet = "",
    mutation_sheet = "",
    deconvolution_functions = ""
  )
}

params <- list(
  sample_name = args[[1]],
  output_dir = args[[2]],
  vep_file = args[[3]],
  snv_file = args[[4]],
  sample_sheet = args[[5]],
  mutation_sheet = args[[6]],
  deconvolution_functions = args[[7]]
)


## function loading
source(params$deconvolution_functions)


## ----print_input_settings, echo = FALSE---------------------------------------
sample_name         <- params$sample_name
sample_sheet        <- data.table::fread(params$sample_sheet)
mutation_sheet      <- params$mutation_sheet

variants_output_dir <- file.path(params$output_dir, "variants")
mutation_output_dir <- file.path(params$output_dir, "mutations")


sample_entry <- sample_sheet %>%
  filter(name == sample_name)

date                <- as.character(sample_entry$date)
location_name       <- as.character(sample_entry$location_name)
coordinates_lat     <- as.character(sample_entry$coordinates_lat)
coordinates_long    <- as.character(sample_entry$coordinates_long)


sigmut_output_file <- file.path(
  mutation_output_dir,
  paste0(
    sample_name,
    "_sigmuts.csv"
  )
)

non_sigmut_output_file <- file.path(
  mutation_output_dir,
  paste0(
    sample_name,
    "_non_sigmuts.csv"
  )
)

mutation_output_file <- file.path(
  mutation_output_dir,
  paste0(sample_name,
    "_mutations.csv")
)

variant_abundance_file <- file.path(
  variants_output_dir,
  paste0(
    sample_name,
    "_variant_abundance.csv"
  )
)

variants_with_meta_file <- file.path(
  variants_output_dir,
  paste0(sample_name, "_variants_with_meta.csv")
)


## ----process_signature_mutations, include = FALSE-----------------------------
# Read signature data
sigmut_df <- read.csv(mutation_sheet, header = TRUE) %>%
  dplyr::select(-matches("source")) %>%
  dplyr::na_if("") %>%
  tidyr::pivot_longer(everything(), values_drop_na = TRUE) %>%
  dplyr::select(variant = name, mutation = value)

vep_output_df <- read.table(params$vep_file, sep = ",", header = TRUE) %>%
  dplyr::na_if("-")

sigmuts_deduped <- sigmut_df %>%
  group_by(mutation) %>%
  summarise(variant = paste(variant, collapse = ","))

# remove mutation type info
sigmuts_deduped_no_gene <- sigmuts_deduped %>%
  mutate(mutation = str_extract(mutation, "(?<=:)[[:alnum:]]+"))

## ----match_snvs_to_signature_mutations, include = FALSE-----------------------
variant_protein_mut <- get_protein_mut(params$vep_file)
# match the variant names according to the signature mutations, if there is no
# signature mutation no name will be given
# variant characterizing is done by NT mutations
variant_protein_mut <- dplyr::left_join(variant_protein_mut,
  sigmuts_deduped_no_gene,
  by = c("gene_mut" = "mutation")
)


## ----merge_vep_with_lofreq_info, include = FALSE------------------------------
# get the SNV frequency values and coverage information for the mutations from
# the LoFreq output
lofreq.info <- as_data_frame(parse_snv_csv(params$snv_file))
vep.info <- variant_protein_mut

complete_df <- dplyr::left_join(lofreq.info, vep.info,
  by = c("gene_mut" = "gene_mut"), copy = TRUE
) %>%
  rowwise() %>%
  mutate(gene_mut_collapsed = paste(genes, gene_mut, sep = ":"))

# TODO: let the read coverage filter be set dynamically over the setting file
complete_cov_filtered.df <- complete_df %>% filter(!(as.numeric(cov) < 100))

# filter for mutations which are signature mutations
match.df <- complete_cov_filtered.df %>%
  filter(!is.na(variant))

# filter for everything that is not a signature mutation
nomatch.df <- complete_cov_filtered.df %>%
  filter(is.na(variant))

cat("Writing signature mutation file to ", sigmut_output_file, "...\n")
write.csv(
  match.df,
  sigmut_output_file
)

cat("Writing non signature mutation file to ", non_sigmut_output_file, "...\n")
write.csv(
  nomatch.df,
  non_sigmut_output_file
)

# Tables are displayed here in report


## ----getting_unique_muts_bulk, include = FALSE--------------------------------
# get  NT mutations only, input for the signature matrix
muations_vec <- match.df$gene_mut_collapsed

# only execute the deconvolution when at least one signature mutation was found
execute_deconvolution <- length(muations_vec) > 0


if (execute_deconvolution) {
  ## ----creating_signature_matrix, include = FALSE-----------------------------
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  # for the deconvolution to work we need the "wild type" frequencies too.
  # The matrix from above got mirrored, wild type mutations are simulated the
  # following: e.g. T210I (mutation) -> T210T ("wild type")
  msig_simple <- createSigMatrix(muations_vec, mutation_sheet) %>%

    # When multiple columns look like the same, the deconvolution will not work,
    # because the function can't distinguish between those columns. The
    # workaround for now is to identify those equal columns and merge them into
    # one, returning also a vector with the information about which of the
    # columns were merged.
    cbind(muts = muations_vec, .)

  msig_transposed <- dedupeDF(msig_simple)
  msig_stable_transposed <- msig_transposed[[1]]
  msig_dedupe_transposed <- msig_transposed[[2]]

  dropped_variants <- c()

  # for every variant update the rownames with the group they are in
  # FIXME: Shorten this and similar constructs
  for (variant in rownames(
    msig_stable_transposed[-(rownames(msig_stable_transposed) %in% "muts"), ])
  ) {
    grouping_res <- dedupeVariants(
      variant,
      msig_stable_transposed,
      msig_dedupe_transposed
    )

    msig_dedupe_transposed <- grouping_res[[1]]
    dropped_variants <- c(dropped_variants, grouping_res[[2]])
  }

  # transpose the data frame back to column format for additional processing
  if (length(msig_dedupe_transposed) >= 1) {
    # the 1 get's rid of the additional first row which is an transposing
    # artifact
    msig_simple_unique <- as.data.frame(t(msig_dedupe_transposed[, -1])) %>%
      mutate(across(!c("muts"), as.numeric))
  }

  # clean the vector to know which variants has to be add with value 0 after
  # deconvolution
  dropped_variants <- unique(dropped_variants)
  dropped_variants <- dropped_variants[!is.na(dropped_variants)]


  ## ----calculate_sigmat_weigths, include = FALSE------------------------------
  deconv_lineages <- colnames(
    msig_simple_unique[, -which(names(msig_simple_unique) %in% c("muts", "WT"))]
  )

  # create list of proportion values that will be used as weigths
  sigmut_proportion_weights <- list()
  for (lineage in deconv_lineages) {
    if (lineage == "Others") {
      # !! 17/02/2022 It's not yet tested how robust this behaves when one would
      # mindlessly clutter the mutationsheet
      # with lineages that are very unlikely to detect or not detected
      value <- nrow(msig_simple_unique) / nrow(sigmuts_deduped)
    } else if (grepl(",", lineage)) {
      group <- unlist(str_split(lineage, ","))
      avrg <- sum(sigmut_df$variant %in% group) / length(group)
      value <- sum(msig_simple_unique[lineage]) / avrg
    } else {
      value <- sum(msig_simple_unique[lineage]) /
        sum(sigmut_df$variant == lineage)
    }
    sigmut_proportion_weights[lineage] <- value
  }

  sigmut_proportion_weights <- as_tibble(sigmut_proportion_weights)

  # applying weights on signature matrix
  # FIXME: there should be a way to do this vectorized
  msig_simple_unique_weighted <- msig_simple_unique
  for (lineage in deconv_lineages) {
    msig_simple_unique_weighted[lineage] <- msig_simple_unique_weighted[lineage] / as.numeric(sigmut_proportion_weights[lineage])
  }


  ## ----simulating_WT_mutations, include = FALSE-------------------------------
  # construct additional WT mutations that are not weighted

  # get bulk frequency values, will be input for the deconvolution function
  bulk_freq_vec <- as.numeric(match.df$freq)

  # construct additional WT mutations that are not weighted
  Others_weight <- as.numeric(sigmut_proportion_weights["Others"])
  msig_stable_all <- simulateOthers(
    muations_vec, bulk_freq_vec,
    msig_simple_unique_weighted[, -which(names(msig_simple_unique_weighted) == "muts")],
    match.df$cov,
    Others_weight
  )

  msig_stable_unique <- msig_stable_all[[1]]


  ## ----deconvolution, include = FALSE-----------------------------------------
  # this hack is necessary because otherwise the deconvolution will throw:
  # Error in x * wts: non-numeric argument to binary operator
  # also see: https://stackoverflow.com/questions/37707060/converting-data-frame-column-from-character-to-numeric/37707117
  sig <- apply(
    msig_stable_unique[, -which(names(msig_stable_unique) %in% "muts")],
    2,
    function(x) {
      as.numeric(as.character(x))
    }
  )

  bulk_all <- as.numeric(msig_stable_all[[2]])

  # central deconvolution step
  variant_abundance <- deconv(bulk_all, sig)


  ## ----plot, echo = FALSE-----------------------------------------------------
  # work in progress...only to show how it theoretically can look like in the
  # report
  variants <- colnames(msig_stable_unique[, -1])
  df <- data.frame(rbind(variant_abundance))

  colnames(df) <- variants
  df <- df %>%
    tidyr::pivot_longer(everything()) %>%
    dplyr::select(variant = name, abundance = value)

  # Handling of ambiguous cases and grouped variants

  # case 1: add dropped variants again with value 0 in case all of the other
  # variants add up to 1
  if (round(sum(df$abundance), 1) == 1) {
    for (variant in dropped_variants) {
      df <- rbind(df, c(variant, 0))
    }
  }

  # case 2: in case "others" == 0, both variants can be split up again and being
  while (any(str_detect(df$variant, ","))) {
    grouped_rows <- which(str_detect(df$variant, ","))
    # fixme: this loop might be unneccessary, since only the first row should
    # been picked, everything else will be handled by the while loop
    for (row in grouped_rows) {
      if (df[row, "abundance"] == 0) {
        grouped_variants <- unlist(str_split(df[row, "variant"], ","))
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          df <- rbind(df, c(variant, 0))
        }
      } else if (df[row, "abundance"] != 0) {
        grouped_variants <- unlist(str_split(df[row, "variant"], ","))
        # normal distribution, devide deconv value by number of grouped variants
        distributed_freq_value <-
          as.numeric(as.numeric(df[row, "abundance"]) /
            length(grouped_variants))
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          df <- rbind(df, c(variant, distributed_freq_value))
        }
      }
      # drop grouped row
      df <- df[-row, ]
    }
  }

  df <- transform(df, abundance = as.numeric(abundance))


  cat("Writing variant abundance file to ", variant_abundance_file, "...\n")
  write.csv(df, variant_abundance_file)

  # plot comes here in report
} else {
  # write dummy variants file
  # TODO: do this as a proper emty table with the correct col names
  file.create(variant_abundance_file)
}

# TODO: check if the else of the above if is handled correctly


## ----csv_output_variant_plot, include = F-------------------------------------
# prepare processed variant values to output them as a csv which will be used for the plots in index.rmd
# those outputs are not offically declared as outputs which can lead to issues - that part should be handled by a seperate
# file (and maybe rule)
output_variant_plot <- data.frame(
  samplename = character(),
  dates = character(),
  location_name = character(),
  coordinates_lat = character(),
  coordinates_long = character()
)
if (!execute_deconvolution) {
  # if no signatur mutation found write empty output file
  # TODO: sombody should check whether this empty file with header is enough, or a more sensible default is required
  write.table(output_variant_plot, variants_with_meta_file,
    sep = "\t",
    na = "NA", row.names = FALSE, quote = FALSE
  )
} else {
  # get all possible variants
  all_variants <- colnames(msig_simple[, -which(names(msig_simple) %in% "muts")])
  # add columns for all possible variants to the dataframe
  for (variant in all_variants) {
    output_variant_plot[, variant] <- numeric()
  }
  meta_data <- c(
    samplename = sample_name,
    dates = date,
    location_name = location_name,
    coordinates_lat = coordinates_lat,
    coordinates_long = coordinates_long
  )

  output_variant_plot <- bind_rows(output_variant_plot, meta_data)

  # get rownumber for current sample
  sample_row <- which(grepl(sample_name, output_variant_plot$samplename))

  # write mutation frequency values to df
  for (i in all_variants) {
    if (i %in% df$variant) {
      # check if variant already has a column
      if (i %in% colnames(output_variant_plot)) {
        output_variant_plot[sample_row, ][i] <- df$abundance[df$variant == i]
        output_variant_plot <- output_variant_plot %>% mutate(others = 1 - rowSums(across(all_of(all_variants)), na.rm = TRUE))
      }
    }
  }

  ## # TODO: This chunk hast to go into a seperate rule
  ## # 2. check if file exists already
  ## if (file.exists (variants_with_meta_file)) {
  ##   previous_df <- read.csv (variants_with_meta_file,
  ##                              header = TRUE, colClasses = "character", check.names = FALSE)
  ##   # convert numeric values to character
  ##   output_variant_plot <- as.data.frame(lapply(output_variant_plot, as.character), check.names = FALSE)
  ##   # merge with adding cols and rows
  ##   output_variant_plot <- full_join(previous_df, output_variant_plot, by = colnames(previous_df), copy = TRUE)
  ## }

  # 3. write to output file
  write.table(output_variant_plot, variants_with_meta_file,
    sep = "\t",
    na = "NA", row.names = FALSE, quote = FALSE
  )
}

## ----csv_output_mutation_plot, include = FALSE------------------------------
# prepare processed mutation values to output them as a csv which will be used
# for the plots in index.rmd those outputs are not officially declared as
# outputs which can lead to issues - that part should be handled by a seperate
# file (and maybe rule)
# get all possible mutations
# one aa mutation can have different codon mutations reported with
# different freqs- for the summary table they have to be summed up
# (process see line 1872 of documentation)
complete_df <- complete_df %>%
  group_by(across(c(-freq, -gene_mut, -gene_mut_collapsed, AA_mut))) %>%
  summarise(
    freq = sum(as.numeric(freq)),
    gene_mut = paste(gene_mut, collapse = ",")
  ) %>%
  rowwise() %>%
  mutate(AA_mut = replace(AA_mut, is.na(AA_mut), "\\:\\")) %>%
  # 211006 this exclusion is necessary because this mutation has a wrong entry
  # in VEP which gives two AA_muts instead of probably 1 deletion
  filter(!(gene_mut %in% "G13477A")) %>%
  ungroup()
# report the gene, translated_AA_mut and NT mut accordingly
# easier to spot translation inconsitentcies that way
all_mutations <- paste(complete_df$AA_mut[!is.na(complete_df$AA_mut)],
  complete_df$gene_mut,
  sep = "::"
)
# 1. write dataframe with this information here
output_mutation_frame <- data.frame(
  samplename = character(),
  dates = character(),
  location_name = character(),
  coordinates_lat = character(),
  coordinates_long = character()
)
# add columns for all possible mutations to the dataframe
for (mutation in all_mutations) {
  output_mutation_frame[, mutation] <- numeric()
}
meta_data <- c(
  samplename = sample_name,
  dates = date,
  location_name = location_name,
  coordinates_lat = coordinates_lat,
  coordinates_long = coordinates_long
)
output_mutation_frame <- bind_rows(output_mutation_frame, meta_data)
# write mutation frequency values to df
for (i in all_mutations) {
  i_nt <- str_split(i, "::")[[1]][2]
  if (i_nt %in% complete_df$gene_mut) { # split gene name to match with AA mut
    # check if variant already has a column
    if (i %in% colnames(output_mutation_frame)) {
      output_mutation_frame[, i] <- complete_df$freq[which(
        complete_df$gene_mut == i_nt
      )]
    }
  }
}
colnames(output_mutation_frame) <- as.character(colnames(
  output_mutation_frame
))
output_mutation_frame <- output_mutation_frame %>%
  dplyr::select(-contains("NA", ignore.case = FALSE))
# convert numeric values to character
output_mutation_frame <- as.data.frame(lapply(
  output_mutation_frame,
  as.character
),
check.names = FALSE
)
# 3. write to output file
cat("Writing mutation file to ", mutation_output_file, "...\n")
write.table(output_mutation_frame, mutation_output_file,
  sep = "\t",
  row.names = FALSE, quote = FALSE
)