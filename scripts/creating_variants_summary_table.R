

library(dplyr)
library(stringr)

get_files <- function ( variants_dir){
  files <- list.files(path = variants_dir,
                    pattern = "_variants.csv",
                    full.names = TRUE,
                    recursive = FALSE)
  return (files)
}

create_summary <- function ( files, output_file ){

  require(dplyr)

  cat(paste("Summarizing", length(files), "variant files."))
  # TODO: turn this while loop into a vectorized binary function with Reduce and lapply

  i <- 0
  while ( i != length(files)){

    if ( file.exists (output_file) ) {
      df1 <- read.table( output_file,
                         sep = "\t", header = TRUE, colClasses = "character", check.names = FALSE )
    } else {
      df1 <- read.table(files[1],
                        sep = "\t", header = TRUE, colClasses = "character", check.names = FALSE)
        # write to output file
      write.table(df1, output_file, sep = "\t",
          row.names = FALSE, quote = FALSE)
      i <- i + 1
    }

    i <- i + 1
    df2 <- read.table(files[i],
                      sep = "\t", header = TRUE, colClasses = "character", check.names = FALSE)
    output_df <- dplyr::full_join(df1, df2, by = intersect(colnames(df1),
                                                    colnames(df2)), copy = TRUE)
    # write to output file
    write.table(output_df, output_file, sep = "\t",
          row.names = FALSE, quote = FALSE)
  }
}

args <- commandArgs (trailingOnly=TRUE)
variants_dir <- args[1]
output_file <- args[2]

files <- get_files( variants_dir )
create_summary( files, output_file )
