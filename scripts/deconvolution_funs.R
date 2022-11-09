library("data.table")
library("MASS")
library("stringr")
library("dplyr")
library("tidyr")

parse_snv_csv <- function(snvfile) { # allele frequency from v-pipe vcf
  #' input: csv file derived from vpipe vcf when using LoFreq,
  #' parsing snv-csv file for coverage, frequency and genomic mutation information

  snvtable <- fread(snvfile)

  ref <- snvtable$Ref
  pos <- snvtable$Pos
  var <- snvtable$Var

  snv_info_df <- data.frame(
    gene_pos = pos,
    gene_mut = paste0(ref, pos, var),
    freq     = snvtable$AF,
    dep      = snvtable$DP
  )

  return(snv_info_df)
}

get_protein_mut <- function(vepfile) {
  # parse together from vep output "Protein_position" and "Amino_acid"
  #' input: []_sarscov2_parsed.txt, parsed vcf from VEP CLI, comma separated
  #' this function is for parsing the information about mutation position in the
  #' amino acid sequence, the reference aa and the alternative mutation into the
  #' aa-mutation-notation which is later on comparable to the lists of signature
  #' mutations

  # reading in whole vep txt output
  # you should include in the vep script to parse out the #
  # in the beginning of the line or include that step here.
  vepfile_df <- fread(vepfile)
  # parsing snv and protmut location


  locations <- vepfile_df %>%

    rename(
      prot_pos  = Protein_position,
      conseq    = Consequence,
      gene_name = SYMBOL
    ) %>%
    # get general info on mutation and its position
    mutate(
      mut_nucs   = str_extract(Uploaded_variation, "[A-Z*]/[A-Z*]$"),
      mut_chrom  = str_extract(Location, "^[^:]+"),
      mut_start  = str_extract(Location, "(?<=:)[0-9]+[^-]"),
      mut_end    = str_extract(Location, "(?<=-)[0-9]+")
    ) %>%
    # generate unique string describing mutation
    mutate(
      mut_ref = str_extract(mut_nucs, "^[A-Z*-]"),
      mut_var = str_extract(mut_nucs,  "[A-Z*-]$"),

      mut_str = paste0(mut_ref, mut_start, mut_var)
    ) %>%

    # get infos on mutation consequences for protein
    mutate(
      # Note: This may behave unexpectedly; in the case of no change aa_ref
      # and aa_var are the same, aa_var will not be empty / NA
      aa_ref = str_extract(Amino_acids, "^[A-Z*-]+"),
      aa_var = str_extract(Amino_acids,  "[A-Z*-]+$"),

      aa_str = paste0(aa_ref, prot_pos, aa_var)
    ) %>%

    # delete all rows with no protein position value
    # TODO Why do we delete these mutations? This just means they do not
    # affect the proteins, the mutations may still be informative...
    filter(!str_detect(prot_pos, "^-")) %>%
    # TODO What is the relevance of this information?
    # specific B117 mutations: 21990-21993, 21764-21770, maybe also 3675-3677,
    # 69-70 - all there

    # remove unneeded cols
    dplyr::select(-matches(names(vepfile_df), ignore.case = FALSE)) %>%

    # This also removes variants that essentially look the same, but differ in
    # their VEP determined Gene
    # e.g. Gene:ENSSASG00005000002 <> Gene:ENSSASG00005000003
    distinct()

  if (any(!is.na(locations$mut_end))) {
    # TODO: look into whether any other consequence might lead to more than one
    # AA ref / var
    # (https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences)
    indel_locations_idcs <- with(
      locations,
      str_detect(conseq, "(deletion)|(insertion)")
    )

    locations_indel_df <- apply(
      locations[indel_locations_idcs, ],
      1,
      function(row) {

        # FIXME Probably more complicated than it needs to be. This just
        # converts the named vector to a dataframe, but with one row instead of
        # one col.
        row <- row %>%
          as.data.frame() %>%
          transpose() %>%
          set_names(names(row))

        if (str_detect(row["conseq"], "insertion")) {
          sep_col <- "aa_var"
          rpl_col <- "aa_ref"
        } else {
          sep_col <- "aa_ref"
          rpl_col <- "aa_var"
        }

        row_sepd <- row %>%
          separate_rows(all_of(sep_col), sep = "(?<=[[:alpha:]])(?!$)")

        row_sepd[2:nrow(row_sepd), rpl_col] <- "-"

        return(row_sepd)
      }
    ) %>%
      bind_rows()

    locations <- locations[- which(indel_locations_idcs), ]

    locations <- locations %>%
      bind_rows(locations_indel_df)
  }

  return(locations)
}

create_sig_matrix <- function(mutations_vector, mutation_sheet_file) {
  #' for making the signature matrix based on the signature mutations found in
  #' the sample (given as input as a vector)for it self
  #' returns simple signature matrix as data.frame without frequency values

  # read in provided mutation sheet
  mutations_df <- fread(mutation_sheet_file) %>%

    # remove source col if there. TODO What problems would this col cause?
    dplyr::select(-matches("source"))

  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant
  # to see how many characterising mutations where found by variant

  sig_mat <- lapply(colnames(mutations_df), function(variant) {
    return(as.numeric(mutations_vector %in% mutations_df[[variant]]))
  }) %>%

    magrittr::set_names(names(mutations_df)) %>%

    dplyr::bind_cols() %>%

    # add "Others" col of all "0"s to indicate possible other variants not
    # possessing any of the mutations
    # TODO ensure this is how this column is supposed to work
    mutate(Others = rep(0, length(mutations_vector))) %>%

    as.matrix() %>%

    magrittr::set_rownames(mutations_vector)

  return(sig_mat)
}


dedupe_sigmut_mat <- function(sigmut_mat, var_sep = "_") {
  # input:
  #   * sigmut_mat:
  #     A matrix with rows corresponing to singnature mutations and cols
  #     corresponding to variants, each dimension named. Entries are 1 if the
  #     mutation of the row is a signature mutation of the variant in the col,
  #     and 0 otherwise.
  #   * var_sep:
  #     A string separating the names of variables with identical columns.
  # output:
  #   * The input matrix with identical columns and their colnames merged.
  #   * A list of character vectors, giving names of variants which have equal
  #     cols.
  variant_names <- colnames(sigmut_mat)

  is_dupe <- duplicated(sigmut_mat, MARGIN = 2)
  dupe_variants <- variant_names[is_dupe]

  # coerce back to dataframe for easier processing
  sigmut_mat_df <- as.data.frame(sigmut_mat)

  # find out of which variant a dupe variant is a dupe of, generate groups
  # of variants which are duplicates of each other
  dupe_group_list <- list()

  for (dupe_var in dupe_variants) {
    if (!dupe_var %in% unique(unlist(dupe_group_list))) {
      dupe_var_col <- sigmut_mat_df[[dupe_var]]

      dupe_group_logi <- apply(
        sigmut_mat,
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
  msig_deduped_df <- sigmut_mat_df[, !is_dupe, drop = FALSE] %>%
    rename(!!dupe_group_names) %>%

    # Simple replace may complain about too long variable names.
    apply(
      2,
      function(col) {
        replace(col, is.na(col), 0)
      },
      simplify = FALSE
    ) %>%
    bind_cols() %>%
    as.data.frame()

  return(msig_deduped_df)
}