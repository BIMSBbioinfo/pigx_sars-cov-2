create_summary <- function(files) {
  cat(paste("Summarizing", length(files), "mutation files.\n"))

  # read files into list
  mutations_list <- lapply(
    X = files,
    FUN = read.table,
    sep = "\t",
    header = TRUE,
    colClasses = "character",
    check.names = FALSE
  )

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

cat("\n\"")
cat(paste(args, collapse = "\",\n\""))
cat("\"\n\n")

output_file <- args[1]
files <- args[2:length(args)]

output <- create_summary(files)

# write to output file
write.table(output, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
