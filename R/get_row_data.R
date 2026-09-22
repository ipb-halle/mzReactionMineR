#' get_rowData
#'
#' Convenience function that extracts the rowData from a SummarizedExperiment
#'     object as a data.frame
#'
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#' @param object a SummarizedExperiment object
#'
#' @returns rowData as data.frame
#' @export
#'
get_row_data <- function(
  object = NULL
) {
  data <- as.data.frame(
    rowData(object)
  )
  return(data)
}

#' @export
get_rowData <- function(...) get_row_data(...)
