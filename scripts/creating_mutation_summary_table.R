
get_files <- function ( mutations_dir){
  files <- list.files(path = mutations_dir,
                    pattern = "_mutations.csv",
                    full.names = TRUE,
                    recursive = FALSE)
  return (files)
}

create_summary <- function ( files, output_file ){

  cat( paste("Summarizing", length(files), "mutation files.") )

    # write to output file
    write.table(output_df, output_file, sep = "\t",
          row.names = FALSE, quote = FALSE)
  }
  # read files into list
  mutations_list <- lapply(X = files,
                          FUN = read.table,
                          sep = "\t",
                          header = TRUE,
                          colClasses = "character",
                          check.names = FALSE )

  # merge mutation files in pairs
  merged_mutations <- Reduce(f = function(df1,df2) merge(df1,
                                                    df2,
                                                    by = intersect(colnames(df1),colnames(df2)),
                                                    all.x = TRUE ,
                                                    all.y = TRUE),
                            x = mutations_list)
}

args <- commandArgs (trailingOnly=TRUE)
mutations_dir <- args[1]
output_file <- args[2]

files <- get_files( mutations_dir )
create_summary( files, output_file )

