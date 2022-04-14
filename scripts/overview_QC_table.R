library(dplyr)
library(R.utils)
library(stringr)

concat_overview_table <- function ( sample_sheet, reads_dir, sample_dir ) {
  
  sample_sheet.df <- read.csv(sample_sheet, header = TRUE, stringsAsFactors = FALSE)
  # get read files matching samples
  # TODO No need for having the reads in different rows, since I'm always counting the whole alignment, or then count the read numbers for the single read files each...still...don't do this over multiple lines - and remove the read file columns afterwards for better readibility
  cat("get samples and reads from sample_sheet...\n")
  read_counts <- read_files( sample_sheet.df)
  # get read number of raw reads
  cat("get num of total raw reads...\n ")
  read_counts$reads_r1 <- read_num_raw(read_counts$file_raw_reads1, reads_dir)$read_num
  read_counts$reads_r2 <- read_num_raw(read_counts$file_raw_reads2, reads_dir)$read_num
  # NOTE Total read number should be halved for paired end reads
  read_counts <- read_counts %>% mutate( total_reads = ifelse(paired_end,
                                                              yes = reads_r1 + reads_r2,
                                                              no = reads_r1))
  # get read number after trimming
  cat("get num of reads after trimming...\n")
  read_counts$read_num_trimmed_r1 <- read_num_raw( read_counts$file_trimmed_reads_r1, file.path(sample_dir,"trimmed_reads"))$read_num
  read_counts$read_num_trimmed_r2 <- read_num_raw( read_counts$file_trimmed_reads_r2, file.path(sample_dir, "trimmed_reads"))$read_num

  cat("get num of unaligned reads ...\n")
  # get read number unaligned from files
  read_counts$unaligned_reads_from_file <- read_num_raw( read_counts$file_unaligned_reads, file.path(sample_dir,"mapped_reads"))$read_num
  
  cat("join counts together...\n")
  read_counts <- left_join(read_counts, parse_amplicons( sample_sheet.df, sample_dir), by = "samplename") 
  
  # for double check, unaligend reads by calculation
  read_counts <- read_counts %>% mutate(unaligned_reads_by_calc = total_reads - as.numeric(aligned_reads) )
  # difference between unaligned reads from file and by calc must be reads filtered by QC
  read_counts <- read_counts %>% mutate(reads_removed_by_QC = unaligned_reads_by_calc - unaligned_reads_from_file )
  # TODO check if removing file names is okay
  read_counts <- read_counts %>% select(!starts_with("file"))
  return( read_counts)
}

parse_amplicons <- function ( sample_sheet.df, sample_dir ){
  samples <- sample_sheet.df$name
  do.call(bind_rows, lapply(samples, apply_fun_parse_coverage_file, sample_dir = sample_dir))
}

read_files <- function ( sample_sheet.df ){
  samples <- sample_sheet.df$name
  do.call( bind_rows, lapply(samples, apply_fun_lookup, sample_sheet = sample_sheet.df))
}

apply_fun_lookup <- function ( sample, sample_sheet.df ){
  sample_row <- which(sample_sheet.df$name == sample)
  # works only for paired end for now - I don't know how to make num of returned read files conditional based on sample_sheet columns
  # TODO make num of returend reads depended on paired or single end
  sample_df <- data.frame(
    samplename = sample,
    date = sample_sheet.df[sample_row, "date"],
    paired_end = sample_sheet.df[sample_row, "reads2"] != "",
    file_raw_reads1 = sample_sheet.df[sample_row,"reads"],
    file_raw_reads2 = sample_sheet.df[sample_row, "reads2"],
    stringsAsFactors = FALSE
  ) %>%
    mutate(file_trimmed_reads_r1 = ifelse(paired_end,
                                     paste0(sample,"_trimmed_R1.fastq.gz"),
                                     paste0(sample, "_trimmed.fastq.gz")),
           file_trimmed_reads_r2 = ifelse(paired_end,
                                     paste0(sample,"_trimmed_R1.fastq.gz"),
                                     ""),
           file_unaligned_reads = paste0(sample,"_unaligned.fastq")
           )

  return(sample_df)
}

read_num_raw <- function ( file_reads_vector, reads_dir){
  do.call( bind_rows, lapply(file_reads_vector, apply_fun_get_read_num, reads_dir = reads_dir))
}

apply_fun_get_read_num <- function (read, reads_dir) {

  file <- file.path(reads_dir, read)
  read_num <-
    if(file.exists(file) & !isDirectory(file)) {
     countLines(file) / 4
  } else {
    0
  }

  data.frame( read_num )
}

apply_fun_parse_coverage_file <- function ( sample, sample_dir ){
  
  coverage_file <- file.path(sample_dir, "coverage", paste0(sample, "_merged_covs.csv"))
  
  coverage.df <- read.table(coverage_file, sep = "\t", header = TRUE)
  coverage.df <- dplyr::na_if(coverage.df, "[]")
  coverage.df[] <- lapply(coverage.df, function (x) gsub('(\\[|\\])','', x))
  data.frame( samplename = sample,
              aligned_reads = coverage.df$Total.number.aligned.reads,
              num_sigmuts_covered = coverage.df$Total.number.of.mutations.covered,
              percentage_refgenome_covered = coverage.df$Percentage.ref.genome.covered)
}   

args <- commandArgs(trailingOnly = TRUE)
sample_dir <- args[1]
sample_sheet <- args[2]
reads_dir <- args[3]
output <- args[4]

df <- concat_overview_table(sample_sheet, reads_dir, sample_dir)
write.csv(df, output, row.names = FALSE)
