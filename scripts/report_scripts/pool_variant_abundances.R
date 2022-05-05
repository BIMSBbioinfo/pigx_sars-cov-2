  # prepare processed variant values to output them as a csv which will be used
  # for the plots in index.rmd those outputs are not offically declared as
  # outputs which can lead to issues - that part should be handled by a seperate
  # file (and maybe rule) get all possible variants
  all_variants <- colnames(msig_simple[, -which(
      names(msig_simple) %in% "muts"
  )])

  output_variants <- file.path(csv_output_dir, "data_variant_plot.csv")

  # 1. write dataframe with this information here
  output_variant_plot <- data.frame(
      samplename = character(),
      dates = character(),
      location_name = character(),
      coordinates_lat = character(),
      coordinates_long = character()
  )

  # add columns for all possible variants to the dataframe
  for (variant in all_variants) {
      output_variant_plot[, variant] <- numeric()
  }

  meta_data <- c(
      samplename = sample_name,
      dates = date,
      location_name = location_name,
      coordinates_lat = coordinates_lat,
      coordinates_long = coordinates_long
  )

  output_variant_plot <- bind_rows(output_variant_plot, meta_data)

  # get rownumber for current sample
  sample_row <- which(grepl(sample_name, output_variant_plot$samplename))

  # write mutation frequency values to df
  for (i in all_variants) {
      if (i %in% df$name) {
          # check if variant already has a column
          if (i %in% colnames(output_variant_plot)) {
              output_variant_plot[sample_row, ][i] <- df$value[df$name == i]
              output_variant_plot <- output_variant_plot %>%
                  mutate(others = 1 - rowSums(across(all_of(all_variants)),
                      na.rm = TRUE
                  ))
          }
      }
  }

  # 2. check if file exists already
  if (file.exists(output_variants)) {
      previous_df <- read.csv(
          output_variants,          
          header = TRUE,
          colClasses = "character",
          check.names = FALSE
      )
      # convert numeric values to character
      output_variant_plot <-
          as.data.frame(lapply(output_variant_plot, as.character),
              check.names = FALSE
          )
      # merge with adding cols and rows
      output_variant_plot <- full_join(previous_df,
          output_variant_plot,
          by = colnames(previous_df),
          copy = TRUE
      )
  }

  # 3. write to output file
  cat("Writing pooled variant file to ", output_variants, "...\n")
  write.csv(output_variant_plot, output_variants,
      na = "NA", row.names = FALSE, quote = FALSE
  )