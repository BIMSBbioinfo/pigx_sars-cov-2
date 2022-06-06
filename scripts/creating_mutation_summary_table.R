library("data.table")

create_summary <- function(files) {
  cat(paste("Summarizing", length(files), "mutation files.\n"))

  # read files into list
  mutations_list <- lapply(
    X = files,
    FUN = fread
  )

  # remove empty files from list
  mutations_list_has_rows <- sapply(mutations_list, nrow)
  mutations_list <- mutations_list[mutations_list_has_rows > 0]

  # merge mutation files in pairs
  merged_mutations <- Reduce(
    f = function(df1, df2) {
      merge(df1,
        df2,
        by = intersect(colnames(df1), colnames(df2)),
        all.x = TRUE,
        all.y = TRUE
      )
    },
    x = mutations_list
  )

  return(merged_mutations)
}

args <- commandArgs(trailingOnly = TRUE)

output_file <- args[1]
files       <- args[2:length(args)]

output <- create_summary(files)

# write to output file
write.csv(output, output_file, row.names = FALSE, quote = FALSE)
