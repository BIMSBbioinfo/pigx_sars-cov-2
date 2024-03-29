library(dplyr)

count_muts <- function (x, mutation_sheet.df) { # x = sample row
  #' function used in rowwise apply() call
  #' takes row as input, calculates mutation counts and returns a dataframe
  #'
  # transform mutation_sheet to one comparable vector
  mutation_sheet.v <- unique( unlist( mutation_sheet.df, use.names = FALSE))
  mutation_sheet.v <- mutation_sheet.v[!is.na(mutation_sheet.v)] 
  
  # transform char. vector into dataframe
  x <- as_tibble(t(as.matrix(x)))
  mutations_ps <- x[ (which( names(x) %in% "coordinates_long")+1) : length( names( x ))]
  count_frame <- data.frame(
     # mutations only
      sample =  as.character(x["samplename"]),
      # count all mutations which are not NA
      total_muts = as.numeric(rowSums(!is.na(mutations_ps))),
      # count all mutations that are signature mutations
      total_sigmuts = as.numeric(rowSums(!is.na(mutations_ps %>% dplyr::select( dplyr::contains(mutation_sheet.v))))),
      # get num of muts with significant increase over time
      sig_incr_muts = as.numeric(rowSums(!is.na( mutations_ps %>% dplyr::select( dplyr::contains(mutations_sig$mutation)))))
  )
  # get number of mutations which aren't signature mutations
  count_frame <- count_frame %>% mutate( non_sigmuts = total_muts - total_sigmuts)
  # get number of siganture mutation per variant
  for ( var in colnames(mutation_sheet.df)){
    count_frame[,paste0("sigmuts_",var)] <- as.numeric(rowSums(!is.na(mutations_ps %>% dplyr::select( dplyr::contains( na.omit(mutation_sheet.df[[var]]))))))
  }
  return(count_frame)
}

write_mutations_count <- function ( mutation_plot_data, mutation_sheet.df, mutations_sig ){
  #' takes data_mut_plot.csv df, mutation_sheet.df with NAs at empty cells, mutations_sig.df as input
  #' counts mutations and return them as a dataframe

  # transform mutation_sheet to one comparable vector
  mutation_sheet.v <- unique( unlist( mutation_sheet.df, use.names = FALSE))
  mutation_sheet.v <- mutation_sheet.v[!is.na(mutation_sheet.v)] # fixme this feels like there should be a more elegent way to remove na right away
  # get names of mutations without meta data
  mutations <- names ( mutation_plot_data[ (which( names(mutation_plot_data) %in% "coordinates_long")+1) : length( names( mutation_plot_data ))])
  # signature mutations found across samples
  sigmuts_found.df <- mutation_plot_data %>% dplyr::select( dplyr::contains(mutation_sheet.v))
  
  # total counts across all samples, without duplicated counts (compared to what I'd get when summing up all the columns)
  if (length(mutations_sig) >0){
    count_frame <- data.frame( sample = "Total",
                               total_muts = length(mutations),
                               total_sigmuts = length(sigmuts_found.df),
                               # get how many of all found mutations will be tracked because of significant increase over time
                               sig_incr_muts = length(mutation_plot_data %>% dplyr::select( dplyr::contains(mutations_sig$mutation)))
    )
    # get number of mutations which aren't signature mutations
    count_frame <- count_frame %>% mutate( non_sigmuts = total_muts - total_sigmuts)
    # get number of siganture mutation per variant
    for ( var in colnames(mutation_sheet.df)){
      sigmut_pv <- mutation_plot_data %>% dplyr::select( dplyr::contains( na.omit(mutation_sheet.df[[var]])))
      count_frame[,paste0("sigmuts_",var)] <- length(sigmut_pv)
    }
    
    counts_per_sample <- do.call(bind_rows, apply(mutation_plot_data,1, count_muts, mutation_sheet.df))
    count_frame <- bind_rows(count_frame, counts_per_sample)
  } else{count_frame <- data.frame()}
  
  return(count_frame)
}