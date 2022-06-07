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
    deconvolution_functions = "",
    sample_name = "",
    mutation_sheet = "",
    sample_sheet = "",
    vep_file = "",
    snv_file = "",
    mutation_depth_threshold = "",
    sigmut_output_file = "",
    non_sigmut_output_file = "",
    variants_output_file = "",
    variants_with_meta_output_file = "",
    mutation_output_file = ""
  )
}

params <- list(
    deconvolution_functions = args[[1]],
    sample_name = args[[2]],
    mutation_sheet = args[[3]],
    sample_sheet = args[[4]],
    vep_file = args[[5]],
    snv_file = args[[6]],
    mutation_depth_threshold = args[[7]],
    sigmut_output_file = args[[8]],
    non_sigmut_output_file = args[[9]],
    variants_output_file = args[[10]],
    variants_with_meta_output_file = args[[11]],
    mutation_output_file = args[[12]]
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

# set script wide separator used for variants that contain the same signature
# muations
var_sep <- ","

# separator for same mutations
mut_sep <- ","

# separator for the string containing mutation infos
mut_str_sep  <- ":"

# character used for na entries in the mut info string
mut_str_na_char <- "\\"

## ----printInputSettings, echo = FALSE-----------------------------------------
sample_name         <- params$sample_name
sample_sheet        <- data.table::fread(params$sample_sheet)
mutation_sheet      <- params$mutation_sheet

sample_entry <- sample_sheet %>%
  filter(name == sample_name)

date                <- as.character(sample_entry$date)
location_name       <- as.character(sample_entry$location_name)
coordinates_lat     <- as.character(sample_entry$coordinates_lat)
coordinates_long    <- as.character(sample_entry$coordinates_long)


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
  summarise(variant = paste(variant, collapse = var_sep))

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
  mutate(mut_str_collapsed = paste(gene_name, mut_str, sep = mut_str_sep))

complete_dep_filtered_df <- complete_df %>%
  filter(as.numeric(dep) > as.numeric(params$mutation_depth_threshold))

# filter for mutations which are signature mutations
match_df <- complete_dep_filtered_df %>%
  filter(!is.na(variant))

# filter for everything that is not a signature mutation
nomatch_df <- complete_dep_filtered_df %>%
  filter(is.na(variant))

# write sigmuts to file
cat(
  "Writing signature mutation file to ",
  params$sigmut_output_file,
  "...\n")

write.csv(
  match_df,
  params$sigmut_output_file
)

# write non sigmuts to file
cat(
  "Writing non signature mutation file to ",
  params$non_sigmut_output_file,
  "...\n")

write.csv(
  nomatch_df,
  params$non_sigmut_output_file
)

# Tables are displayed here in report

## ----echo = FALSE-------------------------------------------------------------
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
  dupe_group_names <- lapply(dupe_group_list, paste, collapse = var_sep) %>%
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
  # FIXME Rename this to something like deconv_groups
  deconv_lineages <- colnames(msig_deduped_df)

  # calculate each variants weighting value
  n_found_vec <- lapply(
    variant_names,
    function(variant) {
      if (variant == "Others") {
        n_mut_found <- nrow(msig_simple_df)
      } else {
        n_mut_found <- sum(msig_simple_df[variant])
      }
      names(n_mut_found) <- variant

      return(n_mut_found)
    }
  ) %>%
    unlist()

  # average weights per group
  # TODO Currently, if "Others" is within a group, it biases that groups
  # weight upwards as no "Others" variant is present in the sigmut_df
  # variant column.
  group_weights_vec <- lapply(
    deconv_lineages,
    function(group) {
      group_vec <- str_split(group, var_sep) %>%
        unlist()

      sel_vec <- names(n_found_vec) %in% group_vec
      mean_n_found <- sum(n_found_vec[sel_vec]) / length(group_vec)

      # Special case: Group consists only of "Others"
      # TODO For now it stays like this to produce the original output,
      # but it is actually not what it is supposed to do, as it ignores
      # the case when a group of more than one variant includes "Others"
      if (group == "Others") {
        # Note: When the mutation sheet contains a large number of variants
        # that are not being detected, that will strongly upweight the XXX

        mean_n_total <- nrow(sigmuts_deduped)
      } else {
        mean_n_total <- sum(sigmut_df$variant %in% group_vec) /
         length(group_vec)
      }

      weight <- mean_n_found / mean_n_total
      names(weight) <- group

      return(weight)
    }
  ) %>%
    unlist()

  # apply weights to signature matrix
  msig_deduped_df_weighted <- msig_deduped_df %>%
    mutate(across(
      everything(), ~ .x / group_weights_vec[[cur_column()]]
    )) %>%
    replace(is.na(.), 0)

  # TODO For downstream compatability only, remove once no longer needed
  sigmut_proportion_weights <- group_weights_vec

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
    separate_rows(variant, sep = var_sep)

  # go trough all groups and assign each group member the group abundance
  # divided by the number of group members
  for (group in dupe_group_list) {
    group_ind <- variant_abundance_df$variant %in% group
    group_abundance <- variant_abundance_df$abundance[group_ind][1]

    variant_abundance_df$abundance[group_ind] <- group_abundance / length(group)
  }

  cat(
    "Writing variant abundance file to ",
    params$variants_output_file,
    "...\n"
  )

  write.csv(variant_abundance_df, params$variants_output_file)

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

  write.table(output_variant_plot, params$variants_with_meta_output_file,
    sep = "\t",
    na = "NA", row.names = FALSE, quote = FALSE
  )

} else {

  # write dummy variants file
  cat("Writing dummy variants file to ", params$variants_output_file, "...\n")
  writeLines(
    "Deconvolution not run, this is a dummy file.",
    params$variants_output_file
  )

  # write dummy variants with meta file
  cat(
    "Writing dummy variants file with metadata to ",
    params$variants_with_meta_output_file,
    "...\n"
  )

  writeLines(
    "Deconvolution not run, this is a dummy file.",
    params$variants_with_meta_output_file
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

# TODO These muations here are not filtered for coverage. Is this intended?
# NOTE Previously NA aa mutations were not deduplicated, this may lead to
# differences in the output for that case.
output_mutation_frame <- complete_df %>%
  group_by(aa_str) %>%
  summarise(
    freq = sum(as.numeric(freq)),
    mut_str = paste(mut_str, collapse = mut_sep)
  ) %>%
  mutate(aa_str = replace(aa_str, is.na(aa_str), paste0(
    mut_str_na_char, mut_str_sep, mut_str_na_char
  ))) %>%
  # 211006 this exclusion is necessary because this mutation has a wrong entry
  # in VEP which gives two AA_muts instead of probably 1 deletion
  filter(!(mut_str %in% "G13477A")) %>%
  ungroup() %>%

  # report the gene, translated_AA_mut and NT mut accordingly
  # easier to spot translation inconsitentcies that way
  mutate(nuc_aa_mut = paste(
    aa_str, mut_str,
    sep = str_glue("{mut_str_sep}{mut_str_sep}")
  )) %>%
  dplyr::select(nuc_aa_mut, freq) %>%

  # Filter aa muts that contain NA
  # TODO Why do some contain NA? I suspect this is not necessary, as all NAs
  # should be converted above.
  filter(!str_detect(nuc_aa_mut, "NA")) %>%
  pivot_wider(names_from = nuc_aa_mut, values_from = freq) %>%

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


# 3. write to output file
cat("Writing mutation file to ", params$mutation_output_file, "...\n")
write.table(output_mutation_frame, params$mutation_output_file,
  sep = "\t",
  row.names = FALSE, quote = FALSE
)
