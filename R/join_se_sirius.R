#' join_se_sirius
#'
#' A convenience function that joins a SIRIUS output file to a
#' SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment rowData rowData<- SummarizedExperiment
#' @importFrom dplyr distinct left_join
#' @importFrom utils read.delim
#' @importFrom methods as is
#' @importFrom stats setNames
#'
#' @param object a SummarizedExperiment object
#' @param path_to_sirius path to the SIRIUS output file that should be joined
#' @param id_col character. The name of the column in the rowData that contains the feature ids. Default is "id".
#'
#' @returns a SummarizedExperiment object with the SIRIUS data joined to the
#' rowData
#' @export
#'
join_se_sirius <- function(
    object = NULL,
    path_to_sirius = NULL,
    id_col = "id"
) {
  if (!is(object, "SummarizedExperiment")) {
    stop("object must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!is.character(path_to_sirius) || length(path_to_sirius) != 1L ||
      is.na(path_to_sirius) || !file.exists(path_to_sirius) ||
      isTRUE(file.info(path_to_sirius)$isdir)) {
    stop("path_to_sirius must be an existing file path.", call. = FALSE)
  }
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col)) {
    stop("id_col must be a single, non-missing column name.", call. = FALSE)
  }

  row_data <- as.data.frame(rowData(object))
  if (!id_col %in% names(row_data)) {
    stop("id_col must be a column in rowData(object).", call. = FALSE)
  }

  sirius_data <- read.delim(path_to_sirius, check.names = FALSE)
  if (ncol(sirius_data) == 0L && nrow(sirius_data) == 0L) {
    rowData(object) <- row_data
    return(object)
  }
  if (!"mappingFeatureId" %in% names(sirius_data)) {
    stop(
      "The SIRIUS output must contain a 'mappingFeatureId' column.",
      call. = FALSE
    )
  }

  # Match SIRIUS identifiers to the type used by the feature row data.
  sirius_data$mappingFeatureId <- as(
    sirius_data$mappingFeatureId,
    class(row_data[[id_col]])
  )
  duplicate_ids <- duplicated(sirius_data$mappingFeatureId)
  if (any(duplicate_ids, na.rm = TRUE)) {
    warning(
      "Duplicate mappingFeatureId values found; keeping the first row for each ID.",
      call. = FALSE
    )
    sirius_data <- distinct(sirius_data, mappingFeatureId, .keep_all = TRUE)
  }

  row_data <- left_join(
    row_data,
    sirius_data,
    by = setNames("mappingFeatureId", id_col)
  )
  rowData(object) <- row_data
  object
}

utils::globalVariables(c("mappingFeatureId"))
