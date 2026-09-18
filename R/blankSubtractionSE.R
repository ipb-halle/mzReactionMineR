#' blankSubtractionSE
#'
#' Use blank samples to filter features in a SummarizedExperiment object based on a specified ratio of sample signal to blank signal.
#'
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#'
#' @param object A SummarizedExperiment object containing the data to be filtered.
#' @param assay  Character. what assay to use. (i.e. "height")
#' @param sample_col Character. Column name in colData that contains sample names (i.e. "filename")
#' @param blanks Character vector. Names of blank samples in the sample_col column of colData
#' @param ratio  Numeric. Minimum ratio of sample signal to blank signal for a feature to be retained. Default is 3.
#' @param ratio_type Character. Method to calculate the blank signal for each feature. Options are "maximum", "mean", or "median". Default is "maximum".
#' @param min_detection_blank Integer. Minimum number of blank samples in which a feature must be detected (i.e. non-NA value) to be considered for filtering. Default is 1.
#' @param remove Logical. Whether to remove the blank samples from the resulting SummarizedExperiment object. Default is TRUE.
#' @param id_col Character. Column name in rowData that contains unique feature identifiers (i.e. "id")
#'
#' @returns A SummarizedExperiment object with features filtered based on the specified ratio.
#' @export
#'

blankSubtractionSE <- function(
    object = NULL,
    assay = NULL,
    blanks = NULL,
    ratio_type = "maximum",
    ratio = 3,
    min_detection_blank = 1,
    remove = TRUE,
    id_col = "id",
    sample_col = "filename"
) {
  object_name <- object_argument_name(substitute(object))
  ratio_type <- match.arg(ratio_type, c("maximum", "mean", "median"))

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }
  if(length(ratio) != 1L || !is.numeric(ratio) || is.na(ratio) ||
     !is.finite(ratio) || ratio < 0) {
    stop("'ratio' must be one non-negative finite number.", call. = FALSE)
  }
  if(length(min_detection_blank) != 1L ||
     !is.numeric(min_detection_blank) || is.na(min_detection_blank) ||
     min_detection_blank < 1 || min_detection_blank %% 1 != 0) {
    stop("'min_detection_blank' must be a positive integer.", call. = FALSE)
  }
  if(length(remove) != 1L || !is.logical(remove) || is.na(remove)) {
    stop("'remove' must be TRUE or FALSE.", call. = FALSE)
  }

  row_data <- as.data.frame(rowData(object))
  column_data <- as.data.frame(colData(object))
  for (column in c(id_col, sample_col)) {
    if(length(column) != 1L || !is.character(column) || is.na(column) ||
       !nzchar(column)) {
      stop("Column names must be non-empty character scalars.", call. = FALSE)
    }
  }
  if(!id_col %in% names(row_data)) {
    stop("'id_col' is missing from rowData(", object_name, ").", call. = FALSE)
  }
  if(!sample_col %in% names(column_data)) {
    stop("'sample_col' is missing from colData(", object_name, ").", call. = FALSE)
  }
  if(anyNA(row_data[[id_col]]) || anyDuplicated(row_data[[id_col]])) {
    stop("'id_col' must contain unique feature identifiers in rowData(", object_name, ").", call. = FALSE)
  }
  if(anyNA(column_data[[sample_col]]) || anyDuplicated(column_data[[sample_col]])) {
    stop("'sample_col' must contain unique, non-missing sample names in colData(", object_name, ").", call. = FALSE)
  }
  if(length(blanks) == 0L || !is.character(blanks) || anyNA(blanks) ||
     any(!nzchar(blanks)) || anyDuplicated(blanks)) {
    stop("'blanks' must contain unique, non-empty sample names.", call. = FALSE)
  }
  if(!all(blanks %in% column_data[[sample_col]])) {
    stop("All values in 'blanks' must be present in 'sample_col'.", call. = FALSE)
  }
  if(!any(!column_data[[sample_col]] %in% blanks)) {
    stop("'blanks' must leave at least one non-blank sample.", call. = FALSE)
  }

  id_cols <- names(row_data)
  measurement_data <- cbind(row_data, as.data.frame(assays(object)[[assay]]))
  sample_data <- column_data[, sample_col, drop = FALSE]

  summarize_signal <- function(sample_names, min_detection = NULL) {
    data <- measurement_data %>%
      pivot_longer(cols = -all_of(id_cols), names_to = sample_col,
                   values_to = "Value") %>%
      inner_join(sample_data, by = sample_col) %>%
      filter(.data[[sample_col]] %in% sample_names, !is.na(.data$Value))

    if(!is.null(min_detection)) {
      data <- data %>%
        group_by(.data[[id_col]]) %>%
        filter(n() >= min_detection)
    }

    data %>%
      group_by(.data[[id_col]]) %>%
      summarize(
        Value = switch(
          ratio_type,
          maximum = max(.data$Value),
          mean = mean(.data$Value),
          median = median(.data$Value)
        ),
        .groups = "drop"
      )
  }

  # Compare the strongest, average, or middle signal in blanks and samples.
  data_blank <- summarize_signal(blanks, min_detection_blank)
  data_samples <- summarize_signal(setdiff(column_data[[sample_col]], blanks))

  ids_to_remove <- data_blank %>%
    left_join(data_samples, by = id_col, suffix = c("_blank", "_sample")) %>%
    filter(.data$Value_sample < ratio * .data$Value_blank |
             is.na(.data$Value_sample)) %>%
    pull(.data[[id_col]])

  if(length(ids_to_remove) == 0L) {
    warning("No features were filtered based on the provided criteria. Returning original ", object_name, ".")
    result <- object
  } else {
    result <- object[!rowData(object)[[id_col]] %in% ids_to_remove, ]
  }

  if(remove) {
    result <- result[, !colData(result)[[sample_col]] %in% blanks]
  }
  return(result)
}

 utils::globalVariables(c(".",".data","Value"))
