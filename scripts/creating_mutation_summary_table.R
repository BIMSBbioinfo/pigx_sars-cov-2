
get_files <- function ( mutations_dir){
  files <- list.files(path = mutations_dir,
                    pattern = "_mutations.csv",
                    full.names = TRUE,
                    recursive = FALSE)
  return (files)
}

create_summary <- function ( files, output_file ){

  cat( paste("Summarizing", length(files), "mutation files.") )

  # read files into list
  mutations_list <- lapply(X = files,
                          FUN = read.csv,                          
                          header = TRUE,
                          colClasses = "character",
                          check.names = FALSE )

  # remove empty files from list
  mutations_list_has_rows <- sapply(mutations_list, nrow)
  mutations_list <- mutations_list[mutations_list_has_rows > 0]

  # merge mutation files in pairs
  merged_mutations <- Reduce(f = function(df1,df2) merge(df1,
                                                    df2,
                                                    by = intersect(colnames(df1),colnames(df2)),
                                                    all.x = TRUE ,
                                                    all.y = TRUE),
                            x = mutations_list)

  return(merged_mutations)
}

args <- commandArgs (trailingOnly=TRUE)
#args <- c("/home/jfreige/proj/pigx_sars-cov-2/tests/output/mutations", "/home/jfreige/proj/pigx_sars-cov-2/tests/output/variants/data_mutation_plot.csv")
mutations_dir <- args[1]
output_file <- args[2]

files <- get_files( mutations_dir )
output <- create_summary( files)
# write to output file
write.csv(output, output_file, row.names = FALSE, quote = FALSE)
