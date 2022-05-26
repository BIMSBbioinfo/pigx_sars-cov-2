
get_files <- function ( variants_dir){
  files <- list.files(path = variants_dir,
                    # TODO: pass file names from args
                    pattern = "_variants_with_meta.csv",
                    full.names = TRUE,
                    recursive = FALSE)
  return (files)
}

create_summary <- function ( files ){

  cat(paste("Summarizing", length(files), "variant files."))

  # read files into list
  variants_list <- lapply(X = files,
                          FUN = read.csv,
                          header = TRUE,
                          colClasses = "character",
                          check.names = FALSE )

  # remove empty files from list
  variants_list_has_rows <- sapply(variants_list, nrow)
  variants_list <- variants_list[variants_list_has_rows > 0]

  # merge variant files in pairs
  merged_variants <- Reduce(f = function(df1,df2) merge(df1,
                                                    df2,
                                                    by = intersect(colnames(df1),colnames(df2)),
                                                    all.x = TRUE ,
                                                    all.y = TRUE),
                            x = variants_list)

  return(merged_variants)

}

args <- commandArgs (trailingOnly=TRUE)
variants_dir <- args[1]
output_file <- args[2]

files <- get_files( variants_dir )
output <- create_summary( files)
# write to output file
write.csv(output, output_file, row.names = FALSE, quote = FALSE)
