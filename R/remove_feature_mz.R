#' removeFeatureMz
#'
#' Removes features from a SummarizedExperiment object based on m/z values.
#' Useful for removing background signals (i.e. column bleed).
#'
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#'
#' @param object a SummarizedExperiment object
#' @param mz a vector of m/z values to remove
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerances for matching features.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#' @param id_col Character. The name of the column in rowData that contains the ids.
#'
#' @returns a SummarizedExperiment object with the features removed
#' @export
#'
remove_feature_mz <- function(
    object,
    mz,
    mz_tolerance = c(0.005, 5),
    mz_col = "mz",
    id_col = "id"
) {
  object_name <- object_argument_name(substitute(object))

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if(!is.numeric(mz) || anyNA(mz) || any(!is.finite(mz))) {
    stop("'mz' must contain only finite numeric values.", call. = FALSE)
  }
  if(length(mz_tolerance) != 2L || !is.numeric(mz_tolerance) ||
     anyNA(mz_tolerance) || any(!is.finite(mz_tolerance)) ||
     any(mz_tolerance < 0)) {
    stop("'mz_tolerance' must contain two non-negative finite numbers.", call. = FALSE)
  }
  if(any(vapply(c(mz_col, id_col), function(column) {
    length(column) != 1L || !is.character(column) ||
      is.na(column) || !nzchar(column)
  }, logical(1)))) {
    stop("'mz_col' and 'id_col' must be single non-empty strings.", call. = FALSE)
  }

  row_data <- rowData(object)
  missing_columns <- setdiff(c(mz_col, id_col), names(row_data))
  if(length(missing_columns) > 0L) {
    stop(
      "Required rowData(", object_name, ") columns are missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  # Resolve each requested m/z to the closest feature within its tolerance.
  ids_to_remove <- unlist(lapply(mz, function(value) {
    get_id(
      object = object,
      mz = value,
      mz_range = calc_mz_range(value, mz_tolerance),
      id_col = id_col,
      mz_col = mz_col
    )
  }), use.names = FALSE)

  # Keep rows whose feature id was not matched by any requested m/z value.
  object[!row_data[[id_col]] %in% ids_to_remove, ]

}

#' @export
removeFeatureMz <- function(...) remove_feature_mz(...)
