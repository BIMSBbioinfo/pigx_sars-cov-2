library(dplyr)
library(data.table)

group_by_day_location <- function( df ){
  df_grouped <- df %>% 
        # discard location
        dplyr::select ( -location_name, -coordinates_long, -coordinates_lat) %>%
        # discard time only keep day
        mutate(dates = as.Date( dates )) %>%
        # pool samples per date and calc. the mean for every mutation column
        group_by(dates)
  return(df_grouped)
}
group_by_day <- function( df ){
  df_grouped <- df %>% # discard time only keep day
        mutate(dates = as.Date( dates )) %>%
        # pool samples per date and calc. the mean for every mutation column
        group_by(dates, location_name, coordinates_long, coordinates_lat)
  return(df_grouped)
}

pool_by_weighted_mean <- function(df, weights, group_fun = c("day_location", "day")) {
  #' docstring missing
  #' weigths is a dataframe with minimum samplenames and total_reads as column
  #' total reads is the number of reads used for alignment of one sample, should
  #' be the sum of read1 and read2 with paired end data

  weights <- weights %>%
    # only take weights from approved samples
    filter(samplename %in% df$samplename)

  # get variants from data frame, all cols not known to be predefined metadata
  # columns
  meta_cols <- c(
    "samplename",
    "dates",
    "location_name",
    "coordinates_lat",
    "coordinates_long"
  )

  variants <- names(df)[!grepl(paste(meta_cols, collapse = "|"), names(df))]

  if (group_fun == "day_location") {
    df_grouped <- group_by_day_location(df)
  } else if (group_fun == "day") {
    df_grouped <- group_by_day(df)
  }

  df_pooled <- df_grouped %>%
    left_join(weights, by = "samplename") %>%
    relocate(total_reads) %>%

    # summarize by calc weighted mean, with the num of raw reads as weight
    summarise_at(
      vars(all_of(variants)),
      list(~ weighted.mean(., total_reads))
    ) %>%

    # rename samples to indicated that they were pooled
    mutate(samplename = ifelse(
      group_fun == "day_location",
      paste0(dates, "_pooled"),
      paste0(dates, "_", location_name, "_pooled")
    )) %>%

    # put names first again
    relocate(samplename) %>%
    ungroup()

  return(df_pooled)
}

pool_by_mean <- function(df, na_handling, group_fun = c("day_location", "day")) {
  #' docstring missing
  #' 
  if (group_fun == "day_location"){
     df_grouped <- group_by_day_location(df)
  } else if (group_fun == "day"){
    df_grouped <- group_by_day(df)
  } 
  
  df_pooled <- df_grouped %>%
                summarize(across(where(is.numeric), mean, na.rm = na_handling)) %>%
                # rename samples to indicated that they were pooled
                mutate(samplename = ifelse(group_fun == "day_location", paste0(dates, "_pooled"), paste0(dates, "_", location_name, "_pooled"))) %>%
                # put names first again
                relocate(samplename) %>%
                ungroup()
  return (df_pooled)
}