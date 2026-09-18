#' Convert a generic feature table to a SummarizedExperiment
#'
#' The feature table must contain one column for each sample named in
#' `sample_meta_data`. All other columns are stored as row data.
#'
#' @importFrom utils read.csv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @param path_to_file Character. Path to the feature table CSV file.
#' @param sample_meta_data Data frame containing sample metadata.
#' @param assay_name Character. Name of the assay to create. Defaults to
#'   `"intensity"`.
#' @param filename_column Character. Name of the metadata column containing
#'   feature-table sample names. Defaults to `"filename"`.
#'
#' @returns A `SummarizedExperiment` object.
#' @export
#'
feature_table_to_se <- function(
    path_to_file,
    sample_meta_data,
    assay_name = "intensity",
    filename_column = "filename"
) {
  if (!is.data.frame(sample_meta_data)) {
    stop("`sample_meta_data` must be a data.frame")
  }
  if (!is.character(filename_column) || length(filename_column) != 1L ||
      !filename_column %in% names(sample_meta_data)) {
    stop("`filename_column` must name a column in `sample_meta_data`")
  }
  if (!is.character(assay_name) || length(assay_name) != 1L ||
      is.na(assay_name) || !nzchar(assay_name)) {
    stop("`assay_name` must be a non-empty character string")
  }

  feature_table <- tryCatch(
    read.csv(path_to_file, check.names = TRUE),
    error = function(error) stop("Could not read feature table: ", error$message)
  )

  # Match metadata names to CSV column names in the same way as R's CSV reader.
  sample_names <- make.names(
    gsub("\\.d$", "", basename(as.character(sample_meta_data[[filename_column]])))
  )
  if (anyNA(sample_names) || any(!nzchar(sample_names))) {
    stop("`sample_meta_data` contains missing or empty sample names")
  }
  if (anyDuplicated(sample_names)) {
    stop("`sample_meta_data` contains duplicate sample names")
  }
  sample_meta_data[[filename_column]] <- sample_names

  available_samples <- intersect(sample_names, names(feature_table))
  if (!length(available_samples)) {
    stop("Feature table does not contain any samples from `sample_meta_data`")
  }
  missing_samples <- setdiff(sample_names, available_samples)
  if (length(missing_samples)) {
    warning(
      "Samples in `sample_meta_data` are not in the feature table and will be removed: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  assay_data <- as.matrix(feature_table[, available_samples, drop = FALSE])
  row_data <- feature_table[, setdiff(names(feature_table), sample_names), drop = FALSE]
  sample_meta_data <- sample_meta_data[
    match(available_samples, sample_meta_data[[filename_column]]),
    ,
    drop = FALSE
  ]

  SummarizedExperiment(
    assays = setNames(list(assay_data), assay_name),
    rowData = row_data,
    colData = sample_meta_data
  )
}