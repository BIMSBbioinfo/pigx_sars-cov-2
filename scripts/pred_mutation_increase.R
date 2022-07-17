library(dplyr)
library(stringr)

refined_lm_model <- function(mutations_df) {
  #' takes data frames with mutations, frequency values over time
  #' returns dataframe with related pvalues

  # TODO check file format assumptions

  # for mutation frequency - a missing value can be assumed as "not found", so
  # NA can be set to 0
  mutations_df <- mutations_df %>%
    replace(is.na(.), 0)

  mutations <- names(mutations_df %>% select(- dates))
  dates     <- mutations_df$dates

  # run a linear model of each mutation col (all cols except dates)
  # to get the change of each mutations abundance over all dates
  all_lm_res <- lapply(
    mutations,
    function(mutation) {
      mutation_col <- mutations_df[[mutation]]

      # only do regression if >= 20% of dates have detected mutations
      # (if not, this will return NULL)
      if (sum(mutation_col > 0) >= 0.2 * length(dates)) {
        lm      <- lm(mutation_col ~ dates)
        lm_sum  <- summary(lm)

        data.frame(
          mutation     = mutation,
          coefficients = lm_sum$coefficients["dates", "Estimate"],
          pvalues      = lm_sum$coefficients["dates", "Pr(>|t|)"]
        )
      }
    }
  )

  lm_res_df <- bind_rows(all_lm_res)

  return(lm_res_df)
}

filter_lm_res_top20 <- function ( lm_res.df, pvalue_cutoff){
  #' lm_res.df: a dataframe with variables, pvalues and regression coefficients 
  #' pvalue_cutoff: a numeric value to use as p-value cutoff for filtering
  
  if (!(all(is.nan(lm_res.df$pvalues)) & all(is.nan(lm_res.df$coefficients)))){
  
    # generate dataframe with significant results only
    lm_res_sig.df <- lm_res.df %>% 
                      # filter for increasing trends only
                      filter( coefficients > 0) %>%
                      # filter for significance
                      filter( pvalues < pvalue_cutoff) %>% 
                      # sort for decreasing coeffs
                      arrange( desc(coefficients)) %>% 
                      # only take the 20 strongest trends
                      slice_head(n = 20)
    
  } else { lm_res_sig.df <- data.frame() }
  
  return ( lm_res_sig.df )
}
