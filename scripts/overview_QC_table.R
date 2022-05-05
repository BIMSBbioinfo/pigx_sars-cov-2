library(dplyr)
library(R.utils)
library(stringr)

concat_overview_table <- function(sample_sheet,
                                  raw_reads_dir,
                                  trimmed_reads_dir,
                                  mapped_reads_dir,
                                  coverage_dir) {
  sample_sheet.df <- read.csv(sample_sheet, header = TRUE, stringsAsFactors = FALSE)
  # get read files matching samples
  cat("get samples and reads from sample_sheet...\n")
  read_counts <- parse_sample_sheet(sample_sheet.df,
                                    raw_reads_dir,
                                    trimmed_reads_dir,
                                    mapped_reads_dir,
                                    coverage_dir)
  # get read number of raw reads
  cat("get num of total raw reads...\n")
  read_counts <- read_counts %>%
    mutate(
      reads_r1 = read_num_fastq(file_raw_reads1),
      reads_r2 = read_num_fastq(file_raw_reads2)
    ) %>%
    mutate(
      # NOTE Total read number should be halved for paired end reads
      total_reads = ifelse(paired_end,
        yes = reads_r1 + reads_r2,
        no = reads_r1
      )
    )
  # get read number after trimming
  cat("get num of reads after trimming...\n")
  read_counts <- read_counts %>%
    mutate(
      read_num_trimmed_r1 = read_num_fastq(file_trimmed_reads_r1),
      read_num_trimmed_r2 = read_num_fastq(file_trimmed_reads_r2)
    )
  cat("get num of unaligned reads ...\n")
  # get read number unaligned from files
  read_counts <- read_counts %>% mutate(
    unaligned_reads_from_file = read_num_fastq(file_unaligned_reads)
  )

  cat("join counts together...\n")
  read_counts <- left_join(read_counts,
                           parse_amplicon_coverage(sample_sheet.df$name,
                                                   coverage_dir),
                           by = "samplename")
  # for double check, unaligend reads by calculation
  read_counts <- read_counts %>%
    mutate(unaligned_reads_by_calc = total_reads - as.numeric(aligned_reads))
  # difference between unaligned reads from file and by calc must be reads filtered by QC
  read_counts <- read_counts %>%
    mutate(reads_removed_by_QC = unaligned_reads_by_calc - unaligned_reads_from_file)
  # TODO check if removing file names is okay
  read_counts <- read_counts %>%
    select(!starts_with("file"))
  return(read_counts)
}


parse_sample_sheet <- function(sample_sheet.df,
                               raw_reads_dir,
                               trimmed_reads_dir,
                               mapped_reads_dir,
                               coverage_dir) {
  # function to interpolate file paths for a given sample sheet
  # works for single and paired end samples
  sample_df <- sample_sheet.df %>%
    transmute(
      samplename = name,
      date,
      file_raw_reads1 = file.path(raw_reads_dir, reads),
      file_raw_reads2 = file.path(raw_reads_dir, reads2)
    ) %>%
    mutate(
      paired_end = file_raw_reads2 != ""
    ) %>%
    mutate(
      file_trimmed_reads_r1 = file.path(
        trimmed_reads_dir,
        ifelse(paired_end,
          paste0(samplename, "_trimmed_R1.fastq.gz"),
          paste0(samplename, "_trimmed.fastq.gz")
        )
      ),
      file_trimmed_reads_r2 = ifelse(paired_end,
        file.path(trimmed_reads_dir, paste0(samplename, "_trimmed_R1.fastq.gz")),
        ""
      ),
      file_unaligned_reads = file.path(
        mapped_reads_dir,
        paste0(samplename, "_unaligned.fastq")
      )
    )

  return(sample_df)
}

read_num_fastq <- function(file_reads_vector) {
  # vectorized function to count the number of reads in a fastq file
  sapply(
    file_reads_vector,
    FUN = function(fastq_file) {
      if (file.exists(fastq_file) & !isDirectory(fastq_file)) {
        countLines(fastq_file) / 4
      } else {
        0
      }
    }
  )
}

parse_amplicon_coverage <- function(samples, coverage_dir) {
  # vectorized function to query and aggregate coverage files
  lapply(samples,
         coverage_dir = coverage_dir,
    FUN = function(sample, coverage_dir) {
      # TODO: Convert below file to csv and pass file name from args
      coverage.df <- read.csv(
        file.path(coverage_dir, paste0(sample, "_merged_covs.csv")),
          header = TRUE) %>%
        dplyr::na_if("[]") %>%
        mutate(across(everything(), str_replace, "(\\[|\\])", "")) %>%
        transmute(
          samplename = sample,
          aligned_reads = as.numeric(Total.number.aligned.reads),
          num_sigmuts_covered = as.numeric(Total.number.of.mutations.covered),
          percentage_refgenome_covered = as.numeric(Percentage.ref.genome.covered)
        )
    }
  ) %>%
    bind_rows()
}

args <- commandArgs(trailingOnly = TRUE)
sample_sheet <- args[1]
output_file <- args[2]
raw_reads_dir     <- args[3]
trimmed_reads_dir <- args[4]
mapped_reads_dir  <- args[5]
coverage_dir      <- args[6]

df <- concat_overview_table(sample_sheet,
                            raw_reads_dir,
                            trimmed_reads_dir,
                            mapped_reads_dir,
                            coverage_dir)
write.csv(df, output_file, row.names = FALSE)
