#' imputeMinFrac
#'
#' Function to impute missing values in a SummarizedExperiment assay with a
#' fraction of the minimum value for each row.
#'
#' @param object A SummarizedExperiment object
#' @param assay The assay to impute
#' @param fraction The fraction of the minimum value to impute (default is 5)
#' @param new_assay_name The name of the new assay to store the imputed values
#'        (default is "imputed")
#'
#' @returns a SummarizedExperiment object with a new assay containing the imputed values
#' @export
#'
impute_min_frac <- function(
    object,
    assay,
    fraction = 5,
    new_assay_name = "imputed"
) {
  object_name <- object_argument_name(substitute(object))

  if (!inherits(object, "SummarizedExperiment")) {
    stop(sprintf("`%s` must be a SummarizedExperiment object.", object_name))
  }
  if (!assay %in% names(assays(object))) {
    stop(sprintf("Assay '%s' was not found in `%s`.", assay, object_name))
  }
  if (!is.numeric(fraction) || length(fraction) != 1L || is.na(fraction) || fraction <= 0) {
    stop("`fraction` must be a single positive numeric value.")
  }
  if (!is.character(new_assay_name) || length(new_assay_name) != 1L || !nzchar(new_assay_name)) {
    stop("`new_assay_name` must be a non-empty character string.")
  }
  if (new_assay_name %in% names(assays(object))) {
    warning(sprintf("Assay '%s' already exists in `%s`. It will be overwritten.", new_assay_name, object_name))
  }

  mat <- assays(object)[[assay]]
  if (!is.numeric(mat)) {
    stop(sprintf("Assay '%s' must contain numeric values.", assay))
  }

  new_mat <- t(apply(mat, 1, function(row) {
    if (all(is.na(row))) {
      stop(sprintf("Cannot impute missing values in assay '%s' because at least one row contains only NA values.", assay))
    }
    replace(row, is.na(row), min(row, na.rm = TRUE) / fraction)
  }))

  assays(object)[[new_assay_name]] <- new_mat
  object
}

#' @export
imputeMinFrac <- function(...) impute_min_frac(...)
