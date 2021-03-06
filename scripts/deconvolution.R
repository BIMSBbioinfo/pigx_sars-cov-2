# requires dplyr and stringr

createSigMatrix <- function ( mutations.vector, mutation_sheet ) {
  #' for making the signature matrix based on the signature mutations found in the sample (given as input as a vector)
  #' for it self
  #' returns simple signature matrix as data.frame without frequency values
  
  # read in provided mutation sheet
  mutations.df <- read.csv(mutation_sheet)
  if ("source" %in% colnames(mutations.df)){
      mutations.df <- mutations.df[,-(which(names(mutations.df) %in% "source"))]
  }
  # create an empty data frame add a column for the Wildtype
  # "Others" means that the particular mutation is not found and the mutation site could be mutated otherwise or not at all
  
  msig <- setNames( data.frame( matrix( ncol = ncol(mutations.df)+1, nrow = 0 )), c("Others", colnames(mutations.df)))
  msig <- bind_rows(tibble(muts=mutations.vector), msig)
  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant to see how many characterising mutations
  # where found by variant
  
  for ( variant in names(mutations.df) ){
    msig[[ variant ]] <- msig$muts %in% mutations.df[[ variant ]]
  }
  msig[is.na(msig)] <- 0
  
  return( msig[,-match('muts', names(msig))]*1 ) # use the *1 to turn true/false to 0/1
}

simulateOthers <- function ( mutations.vector, bulk_freq.vector, simple_sigmat.dataframe, coverage.vector, Others_weight) {
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
dedupeDF <- function( msig_stable ){
  # transpose and add mutations as first column
  msig_stable_transposed <- as.data.frame( cbind( variants = colnames(msig_stable), t( msig_stable ) ))
  # mark duplicated columns, forward and backwards to get ALL the duplicates, otherwise the first one would missing
  dupes_variants <- duplicated ( msig_stable_transposed[,-which(names(msig_stable_transposed) %in% 'variants')], 
                                 fromLast=TRUE)
  msig_dedupe_transposed <- msig_stable_transposed[!dupes_variants,]

  return( list( msig_stable_transposed, msig_dedupe_transposed) )
}

dedupeVariants <- function (variant, variants.df, dedup_variants.df) {
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

deconv <- function (bulk,sig){
  #' This function performs the deconvolution using a signature matrix for the mutations found in the sample 
  #' and bulk frequency values derived by the SNV caller
  #' it was build by Altuna
    
  rlm_model = suppressWarnings(MASS::rlm(sig,bulk, maxit = 100, method = "M"))
      
      
  rlm_coefficients = rlm_model$coefficients
  
  rlm_coefficients = ifelse(rlm_coefficients < 0, 0, rlm_coefficients)
  
  sumOfCof = sum(rlm_coefficients)
  
  rlm_coefficients = rlm_coefficients / sumOfCof  #normalize so coefficients add to 1
  
  as.vector(rlm_coefficients)
}

deconv_debug <- function (bulk,sig){
  #' This function performs the deconvolution using a signature matrix for the mutations found in the sample
  #' and bulk frequency values derived by the SNV caller
  #' it was build by Altuna

  rlm_model = suppressWarnings(MASS::rlm(sig,bulk, maxit = 100, method = "M"))
    
  rlm_coefficients = rlm_model$coefficients

  rlm_coefficients = ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

  sumOfCof = sum(rlm_coefficients)

  rlm_coefficients = rlm_coefficients / sumOfCof  #normalize so coefficients add to 1

  return(list(as.vector(rlm_coefficients), rlm_model$fitted.values))
}