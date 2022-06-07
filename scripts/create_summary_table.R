library("data.table")

create_summary <- function(files) {
  cat(paste("Summarizing", length(files), "files.\n"))

  # read files into list
  table_list <- lapply(
    X = files,
    FUN = fread
  )

  # remove empty files from list
  table_list <- table_list[sapply(table_list, nrow) > 0]

  # merge variant tables in pairs
  merged_table <- Reduce(
    f = function(df1, df2) {
      merge(df1,
        df2,
        by = intersect(colnames(df1), colnames(df2)),
        all.x = TRUE,
        all.y = TRUE
      )
    },
    x = table_list
  )

  return(merged_table)
}

args <- commandArgs(trailingOnly = TRUE)

output_file <- args[1]
files       <- args[2:length(args)]

output <- create_summary(files)

# write to output file
fwrite(output, output_file, row.names = FALSE, quote = FALSE)
