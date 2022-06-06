library("DT")
library("data.table")

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
  capture.output(fwrite(data,
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

# function to order levels in a dataframe column such that strings matching
# a pattern are always the first / last ones (controlled via leading)
# default leading FALSE as this is needed for the plots
give_priority_in_order <- function(vec, pattern, leading = FALSE,
  decreasing = TRUE) {

  # extract all variables that look like "Others". This should probably be only
  # one, but currently there is also lowercase "others"
  # FIXME This could probably be done more elegantly
  others_str <- vec[str_detect(vec, pattern)] %>%
    as.character()

  rest_str <- vec[!str_detect(vec, pattern)] %>%
    as.character() %>%
    sort(decreasing = decreasing)

  if (leading) {
    ordering <- c(others_str, rest_str)
  } else {
    ordering <- c(rest_str, others_str)
  }

  # deduplicate ordering as factor levels need to be unique
  ordering <- unique(ordering)

  return(factor(vec, levels = ordering))
}
