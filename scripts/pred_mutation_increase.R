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
