library("DT")

## ----functions--------------------------------------------------------------------
# function for downloadable tables
create_dt <- function(x) {
  DT::datatable(x, options = list(
    fixedColumns = TRUE,
    scrollY = 180,
    scrollX = TRUE,
    scroller = TRUE,
    dom = "Slfrtip"
  ))
}

# function to dump data in csv formatted string
dump_csv <- function(data, format = "csv", rownames = FALSE) {
  if (format != "csv") stop("Data can only be dumped in CSV format.")
  capture.output(write.csv(data,
                           file = "",
                           quote = FALSE, row.names = rownames
  )) %>%
    paste0(collapse = "\n")
}

# function to create html link with named downloadable files
create_html_download_link <- function(base64enc, filename, text = NULL) {
  if (is.null(text)) text <- sprintf("Download %s", filename)
  sprintf(
    paste0(
      '<a download="%s" href="%s">%s</a>'
    ),
    filename,
    sprintf("data:text/csv;base64,%s", base64enc),
    text
  )
}
# function to create html formatted notes in panel box with optional toggling
# see https://about.gitlab.com/handbook/markdown-guide/#additional-information
create_html_note <- function(x, toggle = FALSE) {
  paste(
    '\n\n<div class="panel panel-info">',
    ifelse(toggle,
           "<details>\n<summary>**Note**</summary>",
           "&nbsp;**Note**"
    ),
    '<div class="panel-body">',
    gsub("##", "\n", x),
    "</div>",
    ifelse(toggle, "", "</details>"),
    "</div>\n\n",
    sep = "\n"
  )
}