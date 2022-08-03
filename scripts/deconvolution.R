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
    deconvolution_functions = "",
    mutation_depth_threshold = ""
  )
}

params <- list(
  sample_name = args[[1]],
  output_dir = args[[2]],
  vep_file = args[[3]],
  snv_file = args[[4]],
  sample_sheet = args[[5]],
  mutation_sheet = args[[6]],
  deconvolution_functions = args[[7]],
  mutation_depth_threshold = args[[8]]
)

# pretty print parameters for easier debugging
cat("Script running with parameters:\n\n")
par_vec <- c()
for (i in seq_along(params)) {
  par_vec[i] <- paste0(names(params)[i], " = \"", params[[i]], "\"")
}
cat(paste(par_vec, collapse = ",\n"))
cat("\n\n")

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
  paste0(
    sample_name,
    "_mutations.csv"
  )
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
  by = c("mut_str" = "mutation")
)


## ----merge_vep_with_lofreq_info, include = FALSE------------------------------
# get the SNV frequency values and read depth information for the mutations from
# the LoFreq output
lofreq_info <- parse_snv_csv(params$snv_file)

complete_df <- dplyr::right_join(
  variant_protein_mut,
  lofreq_info,
  by = c("mut_str" = "gene_mut"), copy = TRUE
) %>%
  mutate(mut_str_collapsed = paste(gene_name, mut_str, sep = ":"))

complete_dep_filtered_df <- complete_df %>%
  filter(as.numeric(dep) > as.numeric(params$mutation_depth_threshold))

# filter for mutations which are signature mutations
match_df <- complete_dep_filtered_df %>%
  filter(!is.na(variant))

# filter for everything that is not a signature mutation
nomatch_df <- complete_dep_filtered_df %>%
  filter(is.na(variant))

cat("Writing signature mutation file to ", sigmut_output_file, "...\n")
write.csv(
  match_df,
  sigmut_output_file
)

cat("Writing non signature mutation file to ", non_sigmut_output_file, "...\n")
write.csv(
  nomatch_df,
  non_sigmut_output_file
)

# Tables are displayed here in report


## ----getting_unique_muts_bulk, include = FALSE--------------------------------
# get  NT mutations only, input for the signature matrix
mutations_vec <- match_df$mut_str_collapsed

# only execute the deconvolution when at least one signature mutation was found
execute_deconvolution <- length(mutations_vec) > 0


if (execute_deconvolution) {
  ## ----creating_signature_matrix, include = FALSE-----------------------------
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  # for the deconvolution to work we need the "wild type" frequencies too.
  # The matrix from above got mirrored, wild type mutations are simulated the
  # following: e.g. T210I (mutation) -> T210T ("wild type")
  msig_simple <- create_sig_matrix(mutations_vec, mutation_sheet)

  # deduplicate matrix

  variant_names <- colnames(msig_simple)

  is_dupe <- duplicated(msig_simple, MARGIN = 2)
  dupe_variants <- variant_names[is_dupe]

  # coerce back to dataframe for easier processing
  msig_simple_df <- as.data.frame(msig_simple)

  # find out of which variant a dupe variant is a dupe of, generate groups
  # of variants which are duplicates of each other
  dupe_group_list <- list()

  for (dupe_var in dupe_variants) {
    if (!dupe_var %in% unique(unlist(dupe_group_list))) {
      dupe_var_col <- msig_simple_df[[dupe_var]]

      dupe_group_logi <- apply(
        msig_simple,
        2,
        function(col, dupe_col) {
          all(col == dupe_col)
        }, dupe_var_col
      )

      dupe_group_vec <- variant_names[dupe_group_logi]

      dupe_group_list[[dupe_group_vec[1]]] <- dupe_group_vec
    }
  }

  # concat dupe groups to form a new composite name for the now unique col
  dupe_group_names <- lapply(dupe_group_list, paste, collapse = ",") %>%
    unlist()

  # juggle names to get a named vector with names and values flipped
  # needed by dplyr::rename()
  old_names <- names(dupe_group_names)
  new_names <- dupe_group_names

  dupe_group_names <- old_names %>%
    set_names(new_names)

  # generate deduped signature matrix
  # is a col was duplicated this contains only the first col of each dupe group
  msig_deduped_df <- msig_simple_df[, !is_dupe] %>%
    rename(!!dupe_group_names)

  ## ----calculate_sigmat_weigths, include = FALSE------------------------------
  deconv_lineages <- colnames(msig_deduped_df)

  # create list of proportion values that will be used as weigths
  sigmut_proportion_weights <- list()
  for (lineage in deconv_lineages) {
    if (lineage == "Others") {
      # !! 17/02/2022 It's not yet tested how robust this behaves when one would
      # mindlessly clutter the mutationsheet
      # with lineages that are very unlikely to detect or not detected

      # n all detected mutations / n all known mutations
      # TODO Find out why
      value <- nrow(msig_deduped_df) / nrow(sigmuts_deduped)
    } else if (grepl(",", lineage)) {
      # as we can not weight the variants separately by their n detected sig
      # muts / n known sig muts (for this group), we use n detected sigmuts
      # (that is still accurate) / group average n known signature mutations

      # TODO Currently, if "Others" is within a group, it biases that groups
      # weight upwards as no "Others" variant is present in the sigmut_df
      # variant column.
      group <- unlist(str_split(lineage, ","))
      avrg <- sum(sigmut_df$variant %in% group) / length(group)
      value <- sum(msig_deduped_df[lineage]) / avrg
    } else {

      # n lineage signature mutations detected in sample /
      # n known lineage signature mutations (provided in the mutation
      # sheet)
      value <- sum(msig_deduped_df[lineage]) /
        sum(sigmut_df$variant == lineage)
    }
    sigmut_proportion_weights[lineage] <- value
  }

  # apply weights to signature matrix
  msig_deduped_df_weighted <- msig_deduped_df %>%
    mutate(across(
      everything(), ~ .x / sigmut_proportion_weights[[cur_column()]]
    )) %>%
    replace(is.na(.), 0)

  ## ----simulating_WT_mutations, include = FALSE-------------------------------
  # construct additional WT mutations that are not weighted
  # for the deconvolution to work we need the "wild type" frequencies too. The
  # matrix from above got (elementwise) inverted, wild type mutations are
  # simulated the following: e.g. T210I (mutation) -> T210T ("wild type")
  #
  # for each mutation, generate a dummy mutation that results in no change

  # get bulk frequency values, will be input for the deconvolution function
  bulk_freq_vec <- as.numeric(match_df$freq)

  # generate frequency values for the dummy mutations. As they represent
  # all the variants without the respective mutation, they are the remainder
  # to one of the original mutations frequencies
  others_freq_vec <- 1 - bulk_freq_vec

  others_select_vec <- str_detect(names(sigmut_proportion_weights), "Others")

  others_prop_wt <- sigmut_proportion_weights[others_select_vec]

  # generate vector with all mutation frequencies: dummy and real
  bulk_all <- c(others_freq_vec, bulk_freq_vec)

  # make matrix with Others mutations and inverse the values and wild type
  # freqs
  msig_inverse <- msig_deduped_df_weighted %>%
    mutate(across(everything(), ~ as.numeric(!as.logical(.x)))) %>%

    # apply weights right away
    mutate(across(
      everything(),
      ~ .x / as.numeric(others_prop_wt)
    ))


  # generate combined signature matrix for variants, dummy and real
  msig_all <- rbind(msig_inverse, msig_deduped_df_weighted) %>%
    as.matrix()


  ## central deconvolution step ------------------------------------------------
  variant_abundance <- deconv(bulk_all, msig_all)


  ## ----plot, echo = FALSE-----------------------------------------------------
  variant_abundance_df <- data.frame(
    variant = deconv_lineages,
    abundance = variant_abundance
  ) %>%
    separate_rows(variant, sep = ",")

  # go trough all groups and assign each group member the group abundance
  # divided by the number of group members
  for (group in dupe_group_list) {
    group_ind <- variant_abundance_df$variant %in% group
    group_abundance <- variant_abundance_df$abundance[group_ind][1]

    variant_abundance_df$abundance[group_ind] <- group_abundance / length(group)
  }

  cat("Writing variant abundance file to ", variant_abundance_file, "...\n")
  write.csv(variant_abundance_df, variant_abundance_file)

  ## ----csv_output_variant_plot, include = F-------------------------------------
  # prepare processed variant values to output them as a csv which will be
  # concatenated across samples and used for the plots in index.rmd.

  # NOTE: previously an additional column called lowercase "others" was
  # calculated as 1-sum(all other variants) due to the idea behind the uppercase
  # "Others" col, this was always 0 and supposed to be already fixed.
  output_variant_plot <- variant_abundance_df %>%
    pivot_wider(names_from = variant, values_from = abundance) %>%
    mutate(
      samplename = sample_name,
      dates = date,
      location_name = location_name,
      coordinates_lat = coordinates_lat,
      coordinates_long = coordinates_long
    ) %>%

    # ensure metadata cols are first
    dplyr::select(all_of(c(
      "samplename",
      "dates",
      "location_name",
      "coordinates_lat",
      "coordinates_long"
    )), everything())

  write.table(output_variant_plot, variants_with_meta_file,
    sep = "\t",
    na = "NA", row.names = FALSE, quote = FALSE
  )
  
} else {
  cat("Writing dummy variants file to ", variant_abundance_file, "...\n")

  # write dummy variants file
  writeLines(
    "Deconvolution not run, this is a dummy file.",
    variant_abundance_file
  )

  cat(
    "Writing dummy variants file with metadata to ",
    variants_with_meta_file,
    "...\n"
  )

  # write dummy variants file
  writeLines(
    "Deconvolution not run, this is a dummy file.",
    variants_with_meta_file
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
  group_by(across(c(-freq, -mut_str, -mut_str_collapsed, aa_str))) %>%
  summarise(
    freq = sum(as.numeric(freq)),
    mut_str = paste(mut_str, collapse = ",")
  ) %>%
  rowwise() %>%
  mutate(aa_str = replace(aa_str, is.na(aa_str), "\\:\\")) %>%
  # 211006 this exclusion is necessary because this mutation has a wrong entry
  # in VEP which gives two AA_muts instead of probably 1 deletion
  filter(!(mut_str %in% "G13477A")) %>%
  ungroup()
# report the gene, translated_AA_mut and NT mut accordingly
# easier to spot translation inconsitentcies that way
all_mutations <- paste(complete_df$aa_str[!is.na(complete_df$aa_str)],
  complete_df$mut_str,
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
  if (i_nt %in% complete_df$mut_str) { # split gene name to match with AA mut
    # check if variant already has a column
    if (i %in% colnames(output_mutation_frame)) {
      output_mutation_frame[, i] <- complete_df$freq[which(
        complete_df$mut_str == i_nt
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