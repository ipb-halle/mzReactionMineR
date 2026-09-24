#' se_to_long
#'
#' Generates a long form data.frame from a summarized experiment produced by
#'     mzmine_to_se
#'
#' @importFrom dplyr %>% inner_join
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom SummarizedExperiment assays rowData SummarizedExperiment
#' @param object a SummarizedExperiment object
#' @param assay a character string indicating the assay to be used
#' @param filename a character string indicating the name of the column in the colData that corresponds to the assay columns. Default is "filename".
#'
#' @returns a data.frame in long format
#' @export
#'
se_to_long <- function(
  object,
  assay,
  filename = "filename"
) {
  object_name <- object_argument_name(substitute(object))

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if(length(assay) != 1L || !is.character(assay) || is.na(assay) ||
      !nzchar(assay) || !assay %in% names(assays(object))) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }
  if(length(filename) != 1L || !is.character(filename) || is.na(filename) ||
     !nzchar(filename)) {
    stop("'filename' must be one non-empty column name.", call. = FALSE)
  }

  column_data <- get_col_data(object)
  if(!filename %in% names(column_data)) {
    stop("'filename' must name a column in colData(", object_name, ").", call. = FALSE)
  }

  assay_data <- assays(object)[[assay]]
  assay_columns <- colnames(assay_data)
  sample_names <- column_data[[filename]]
  if(is.null(assay_columns) || anyNA(assay_columns) ||
     any(!nzchar(assay_columns)) || anyDuplicated(assay_columns)) {
    stop("The selected assay must have unique, non-missing column names.", call. = FALSE)
  }
  if(anyNA(sample_names) || anyDuplicated(sample_names)) {
    stop("colData(", object_name, ")[['", filename,
         "']] must contain unique, non-missing sample names.", call. = FALSE)
  }
  if(!setequal(assay_columns, sample_names)) {
    stop("The selected assay column names must match colData(", object_name,
         ")[['", filename, "']].", call. = FALSE)
  }

  df <- cbind(
    get_row_data(object), assay_data
  ) %>%
    pivot_longer(
      cols = -all_of(colnames(rowData(object))),
      names_to = filename,
      values_to = "value"
    ) %>%
    inner_join(
      column_data,
      by = filename
    )
  return(df)
}
