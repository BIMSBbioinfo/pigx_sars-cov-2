library(dplyr)

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

read_files <- function ( sample_sheet.df ){
  samples <- sample_sheet.df$name
  do.call( bind_rows, lapply(samples, apply_fun_lookup, sample_sheet = sample_sheet.df))
}

apply_fun_lookup <- function ( sample, sample_sheet.df ){
  sample_row <- which(sample_sheet.df$name == sample)
  # works only for paired end for now - I don't know how to make num of returned read files conditional based on sample_sheet columns
  # TODO make num of returend reads depended on paired or single end
  data.frame( samplename = sample,
              raw_reads1 = sample_sheet.df[sample_row,"reads"],
              raw_reads2 = sample_sheet.df[sample_row, "reads2"])
}

read_num_raw <- function ( raw_reads_vector, reads_dir){
  do.call( bind_rows, lapply(raw_reads_vector, apply_fun_get_read_num, reads_dir = reads_dir))
}

apply_fun_get_read_num <- function (read, reads_dir) {
  read_num <- as.numeric(
                      system2(command = "echo", 
                      args = c ("$(zcat ", file.path(reads_dir, read), "|wc -l)/4|bc"),
                      stdout = TRUE))
  data.frame( read_num )
}

get_num_raw_reads <- function (reads_dir, sample_sheet){
  
  sample_sheet.df <- read.csv(sample_sheet, header = TRUE)
  # get read files matching samples
  cat("get samples and reads from sample_sheet...\n")
  read_counts <- read_files ( sample_sheet.df )
  # fixme can do in one line with dplyr and across I think
  read_counts$reads_r1 <- read_num_raw(read_counts$raw_reads1, reads_dir)$read_num
  read_counts$reads_r2 <- read_num_raw(read_counts$raw_reads2, reads_dir)$read_num
  read_counts <- read_counts %>% mutate( total_reads = reads_r1 + reads_r2 )
  
  return(read_counts)
}

pool_by_weighted_mean <- function(df, weights, group_fun = c("day_location", "day")) {
  #' docstring missing
  #' weigths is a dataframe with minimum samplenames and total_reads as column 
  #' total reads is the number of reads used for alignment of one sample, should be the sum of read1 and read2 with 
  #' paired end data

  weights <- weights  %>% 
              # only take weights from approved samples
              semi_join(df, by = "samplename")
  # get variants from data frame, being every column after the metadata
  variants <- names (
              df[ ( which( names(df) %in% "coordinates_long")+1) : length( names( df ))]
  )
  
  if (group_fun == "day_location"){
     df_grouped <- group_by_day_location(df)
  } else if (group_fun == "day"){
    df_grouped <- group_by_day(df)
  }

  df_pooled <- df_grouped %>%
                left_join(weights, by = "samplename")  %>%
                relocate(total_reads) %>%
                # summarize by calc weighted mean, with the num of raw reads as weight
                summarise_at(vars( all_of(variants) ), list(~ weighted.mean(., total_reads))) %>%
                # rename samples to indicated that they were pooled
                mutate(samplename = ifelse(group_fun == "day_location", paste0(dates, "_pooled"), paste0(dates, "_", location_name, "_pooled"))) %>%
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