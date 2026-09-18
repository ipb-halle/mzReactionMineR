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
  data <- base::as.data.frame(
    SummarizedExperiment::rowData(object)
  )
  return(data)
}

get_rowData <- function(...) get_row_data(...)
