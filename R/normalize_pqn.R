#' Probabilistic Quotient Normalization (PQN)
#'
#' This function performs probabilistic quotient normalization on a
#' SummarizedExperiment assay.
#'
#' Probabilistsic quotient normalization (PQN) is a method to normalize data by
#'    dividing each value within a sample by a normalization factor.
#'
#'    \deqn{normalization\_factor = median \left( \frac{value_{i,j}}{median(feature_i)} \right)}
#'
#'    Where \eqn{value_{i,j}} is the value of the i-th feature in the j-th sample and
#'    \eqn{median(feature_i)} is the median of the i-th feature across all samples.
#'
#'    If measure = "mean", the normalization factor is calculated as:
#'
#'    \deqn{normalization\_factor = median \left( \frac{value_{i,j}}{mean(feature_i)} \right)}
#'
#'    The normalized matrix is then simply the quotient of a \eqn{value_{i,j}} by
#'    the \eqn{normalization_factor_{j}} of a given sample j.
#'
#'    \deqn{normalized\_data = \frac{value_{i,j}}{normalization\_factor_j}}
#'
#' @importFrom MatrixGenerics rowMedians
#' @importFrom stats median
#'
#' @param object A SummarizedExperiment object.
#' @param assay Character. The name of the assay to be normalized.
#' @param type What measure to be used for the normalization factor.
#'     Can be "median" or "mean".
#' @param reference_samples Samples used to calculate the feature reference.
#'     Character vector of unique sample names. Defaults to all samples.
#' @param new_assay_name Character. The name of the new assay containing the
#'     normalized values. Defaults to "pqn_normalized".
#'
#' @returns A SummarizedExperiment object with a new assay containing the
#'     normalized values.
#' @export
#

normalize_pqn <- function(
    object,
    assay,
    type = c("median", "mean"),
    reference_samples = NULL,
    new_assay_name = "pqn_normalized"
) {
  object_name <- object_argument_name(substitute(object))

  if (!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if (length(assay) != 1L || !is.character(assay) ||
      !assay %in% names(assays(object))) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }

  type <- match.arg(type)
  raw_data <- as.matrix(assays(object)[[assay]])

  if (!is.numeric(raw_data) || length(raw_data) == 0L ||
      any(dim(raw_data) == 0L)) {
    stop("The selected assay must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(raw_data), na.rm = TRUE)) {
    stop("The selected assay may contain NA values but no infinite values.",
         call. = FALSE)
  }

  # Resolve the named reference subset to assay column positions.
  if (is.null(reference_samples)) {
    reference_samples <- seq_len(ncol(raw_data))
  } else {
    if (length(reference_samples) == 0L || !is.character(reference_samples) ||
        anyNA(reference_samples) || any(!nzchar(reference_samples)) ||
        anyDuplicated(reference_samples)) {
      stop("`reference_samples` must contain unique, non-empty sample names.",
           call. = FALSE)
    }
    if (is.null(colnames(raw_data))) {
      stop("Sample names require column names in the selected assay.",
           call. = FALSE)
    }
    if (!all(reference_samples %in% colnames(raw_data))) {
      stop("All values in `reference_samples` must be present in the assay.",
           call. = FALSE)
    }
    reference_samples <- match(reference_samples, colnames(raw_data))
  }

  # Build feature references from the selected samples, then compare every
  # sample against those references to obtain its normalization factor.
  reference_data <- raw_data[, reference_samples, drop = FALSE]

  reference <- switch(
    type,
    median = rowMedians(reference_data, na.rm = TRUE),
    mean = rowMeans(reference_data, na.rm = TRUE)
  )
  if (any(!is.finite(reference) | reference == 0)) {
    stop("Each feature must have a finite, non-zero reference value.", call. = FALSE)
  }

  normalization_factor <- apply(
    sweep(raw_data, 1L, reference, "/"), 2L, median, na.rm = TRUE
  )
  if (any(!is.finite(normalization_factor) | normalization_factor == 0)) {
    stop("Each sample must have a finite, non-zero normalization factor.", call. = FALSE)
  }

  divide_by_feature(
    object = object,
    assay = assay,
    vector = normalization_factor,
    new_assay_name = new_assay_name
  )

}

