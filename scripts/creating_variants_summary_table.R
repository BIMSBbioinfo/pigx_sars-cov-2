library("data.table")

create_summary <- function(files) {
  cat(paste("Summarizing", length(files), "variant files.\n"))

  # read files into list
  variants_list <- lapply(
    X = files,
    FUN = fread
  )

  # remove empty files from list
  variants_list_has_rows <- sapply(variants_list, nrow)
  variants_list <- variants_list[variants_list_has_rows > 0]

  # merge variant files in pairs
  merged_variants <- Reduce(
    f = function(df1, df2) {
      merge(df1,
        df2,
        by = intersect(colnames(df1), colnames(df2)),
        all.x = TRUE,
        all.y = TRUE
      )
    },
    x = variants_list
  )

  return(merged_variants)
}

args <- commandArgs(trailingOnly = TRUE)

output_file <- args[1]
files       <- args[2:length(args)]

output <- create_summary(files)

# write to output file
fwrite(output, output_file, row.names = FALSE, quote = FALSE)
