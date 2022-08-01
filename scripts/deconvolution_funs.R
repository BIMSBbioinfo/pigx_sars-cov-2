library("MASS")
library("stringr")
library("dplyr")
library("tidyr")

parse_snv_csv <- function(snvfile) { # allele frequency from v-pipe vcf
  #' input: csv file derived from vpipe vcf when using LoFreq,
  #' parsing snv-csv file for coverage, frequency and genomic mutation information

  snvtable <- read.table(snvfile, sep = ",", header = TRUE)

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
  vepfile_df <- read.table(vepfile, sep = ",", header = TRUE)
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
  mutations_df <- read.csv(mutation_sheet_file) %>%

    # remove source col if there. TODO What problems would this col cause?
    dplyr::select(-matches("source"))

  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant
  # to see how many characterising mutations where found by variant

  sig_mat <- lapply(colnames(mutations_df), function(variant) {
    return(as.numeric(mutations_vector %in% mutations_df[[variant]]))
  }) %>%

    dplyr::bind_cols() %>%

    magrittr::set_names(names(mutations_df)) %>%

    # add "Others" col of all "0"s to indicate possible other variants not
    # possessing any of the mutations
    # TODO ensure this is how this column is supposed to work
    mutate(Others = rep(0, length(mutations_vector))) %>%

    magrittr::set_rownames(mutations_vector)


  return(sig_mat)
}

simulate_others <- function ( mutations.vector, bulk_freq.vector, simple_sigmat.dataframe, coverage.vector, Others_weight) {
  #' for the deconvolution to work we need the "wild type" frequencies too. The matrix from above got mirrored, 
  #' wild type mutations are simulated the following: e.g. T210I (mutation) -> T210T ("wild type")
  
  # 1. make "Others mutations" 
  muts_Others <- lapply(mutations.vector,function(x) str_replace(x,regex(".$"), 
                                                             str_sub(str_split(x,":")[[1]][2], 1,1)))
  muts_Others.df <- data.frame(muts = unlist(muts_Others))
  # 2. make frequency values, subtract the observed freqs for the real mutations from 1
  bulk_Others <- lapply(bulk_freq.vector, function (x) {1-x})
  
  # 3. make matrix with Others mutations and inverse the values and wild type freqs
  msig_inverse <- bind_cols(muts_Others.df, as.data.frame(+(!simple_sigmat.dataframe)))
    
  # 4. apply Others weight 
  # fixme: it could be this can be implemented in the step above already
  msig_inverse[ msig_inverse == 1] <- 1/Others_weight
  
  # fixme: not sure if this really is a nice way to concat those things...
    # no it's not you could use dplyr and mutate
  muts_all <- c(muts_Others,mutations.vector)
  muts_all.df <- data.frame(muts = unlist(muts_all))
  
  bulk_all <- c(bulk_Others, bulk_freq.vector)
  bulk_all.df <- data.frame(freq = unlist(bulk_all))
    
  coverage_all <- c(coverage.vector,coverage.vector)
  coverage_all.df <- data.frame(cov = unlist(coverage_all))
  
  msig_all <- rbind(msig_inverse[,-which(names(msig_inverse) %in% 'muts')],simple_sigmat.dataframe)
  
  # 4. concat the data frames
  # without bulk freq for building the signature matrix
  msig_stable <- bind_cols(muts_all.df,msig_all)
  
  # with bulk freq for export and overview
  msig_stable_complete <- bind_cols(muts_all.df,msig_all,bulk_all.df,coverage_all.df)
  
  return ( list(msig_stable, bulk_all, msig_stable_complete) )
}

# When multiple columns look like the same, the deconvolution will not work, because the function can't distinguish 
# between those columns. The workaround for now is to identify those equal columns and merge them into one, returning also
# a vector with the information about which of the columns were merged. 
# deduplicate dataframe
dedupe_df <- function( msig_stable ){
  # transpose and add mutations as first column
  msig_stable_transposed <- as.data.frame(cbind(
    variants = colnames(msig_stable),
    t(msig_stable)
  ))

  # mark duplicated columns, forward and backwards to get ALL the duplicates,
  # otherwise the first one would missing
  dupes_variants <- duplicated(
    msig_stable_transposed[, -which(names(msig_stable_transposed) %in% "variants")],
    fromLast = TRUE
  )

  msig_dedupe_transposed <- msig_stable_transposed[!dupes_variants, ]

  return(list(msig_stable_transposed, msig_dedupe_transposed))
}

dedupe_variants <- function (variant, variants.df, dedup_variants.df) {
        # get variant group per mutation pattern
        # duped_variants <- grep (variant, variants.df$variants)
        duped_variants <- c()
        row_number_variant <- which( grepl( variant, variants.df$variants ))
        for (row in 1:nrow( variants.df )) { 
            if (all ( variants.df[row_number_variant,-1] == variants.df[row,-1] )) { # TODO: what are those magic numbers?
              duped_variants <- c(duped_variants, variants.df[row,"variants"])
            }
        }
       # grouped_variants <- variants.df$variants[duped_variants]
        groupName_variants <- paste( duped_variants, collapse = "," )
        for ( row in dedup_variants.df$variants ){
          if ( grepl( row,groupName_variants )) {
            
            # if variants are getting pooled with Others they are just Others and nothing else
            if (str_detect(groupName_variants, "Others")){
              rownames ( dedup_variants.df )[rownames(dedup_variants.df) == row] <- "Others"
              variants_to_drop <- duped_variants[!grepl("Others",duped_variants)]
            } else{
              rownames ( dedup_variants.df )[rownames(dedup_variants.df) == row] <- groupName_variants
              variants_to_drop <- NA
              # TODO you can stop after this ( I think)
            }
        }
        }
        # clean the vector to know which variants has to be add with value 0 after deconvolution
        variants_to_drop <- unique(variants_to_drop)[!is.na(variants_to_drop)]
        return ( list(dedup_variants.df, variants_to_drop) )
}  


deconv <- function(bulk, sig) {
  #' This function performs the deconvolution using a signature matrix for the
  #' mutations found in the sample and bulk frequency values derived by the SNV
  #' caller
  #' it was build by Altuna

  rlm_model <- suppressWarnings(MASS::rlm(sig, bulk, maxit = 100, method = "M"))


  rlm_coefficients <- rlm_model$coefficients

  rlm_coefficients <- ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

  sum_of_cof <- sum(rlm_coefficients)

  # normalize so coefficients add to 1
  rlm_coefficients <- rlm_coefficients / sum_of_cof

  return(rlm_coefficients)
}

deconv_debug <- function(bulk, sig) {
  #' This function performs the deconvolution using a signature matrix for the
  #' mutations found in the sample and bulk frequency values derived by the SNV
  #' caller
  #' it was build by Altuna

  rlm_model <- suppressWarnings(MASS::rlm(sig, bulk, maxit = 100, method = "M"))

  rlm_coefficients <- rlm_model$coefficients

  rlm_coefficients <- ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

  sum_of_cof <- sum(rlm_coefficients)

  # normalize so coefficients add to 1
  rlm_coefficients <- rlm_coefficients / sum_of_cof

  return(list(as.vector(rlm_coefficients), rlm_model$fitted.values))
}
