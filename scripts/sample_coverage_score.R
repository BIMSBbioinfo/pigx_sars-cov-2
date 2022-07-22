library(data.table)

get_genome_cov <- function ( coverage_dir, samples_names ) {

  files <- list.files(path = coverage_dir,
                      pattern = "_genome_cov.csv",
                      full.names = TRUE,
                      recursive = FALSE) # why false?
  # stop if no coverage files exists in dir or if the dir does not exist
  if (length(files) == 0) {
    stop("coverage directory is empty or does not exist")
  }

  # define expected file header
  samtools_coverage_header <-
    c(
      "#rname",
      "startpos",
      "endpos",
      "numreads",
      "covbases",
      "coverage",
      "meandepth",
      "meanbaseq",
      "meanmapq"
    )

  # iterate over files and return either data.frame with correct columns
  # ("samplename","ref_genome_coverage") or NULL if column names do not
  # match expected header
  genome_cov.dfs <- lapply( files, function (x){
    rd_tbl <- data.table::fread(x)
    data <- if (all(names(rd_tbl) %in% samtools_coverage_header)) {
      data.frame(
        samplename = strsplit(basename(x), "\\_genome_cov.csv")[[1]],
        ref_genome_coverage = as.numeric(rd_tbl$coverage)
      )
    } else {
      NULL
    }

    return(data)
  })
  # combine data.frames by row, automatically dropping NULL elements
  # NOTE: rbindlist is almost thirtee times faster than bind_rows()
  genome_cov.df <- as.data.frame(data.table::rbindlist(genome_cov.dfs))

  return(genome_cov.df)
}

get_mutation_cov <- function ( coverage_dir ) {

  # TODO: convert these files to proper csv and make names dynamic
  files <- list.files(path = coverage_dir,
                      pattern = "_merged_covs.csv",
                      full.names = TRUE,
                      recursive = FALSE) # why false?
  # stop if no coverage files exists in dir or if the dir does not exist
  if (length(files) == 0) {
    stop("coverage directory is empty or does not exist")
  }

  mutation_cov.df <- data.frame(samplename = c(),
                                total_num_muts = c(),
                                total_muts_cvrd = c(),
                                drop_out_muts = c())

  mutation_cov_header <- c(
    "Total number of tracked mutations",
    "Total number of mutations covered",
    "Number of mutations not covered",
    "Total number aligned reads",
    "Percentage ref.genome covered",
    "Mean depth ref.genome coverage"
  )

  mutation_cov.dfs <- lapply( files, function (x){
    rd_tbl <- fread(x)
    # check if file has the correct header format
    data <- if (all(names(rd_tbl) %in% mutation_cov_header)) {
      data.frame(
        samplename = strsplit(basename (x), "\\_merged_covs.csv")[[1]],
        total_num_muts = rd_tbl$`Total number of tracked mutations`,
        total_muts_cvrd = rd_tbl$`Total number of mutations covered`,
        drop_out_muts = gsub(pattern = "\\[|\\]|\\s",
                             replacement = "",
                             x = rd_tbl$`Number of mutations not covered`)
      )
    } else {
      NULL
    }
    return(data)
  })

  # combine data.frames by row, automatically dropping NULL elements
  # NOTE: rbindlist is almost thirtee times faster than bind_rows()
  mutation_cov.df <- as.data.frame(data.table::rbindlist(mutation_cov.dfs))


  return(mutation_cov.df)
}
