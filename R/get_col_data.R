#' get_colData
#'
#' Convenience function that extracts the colData from a SummarizedExperiment
#'     object as a data.frame
#'
#' @importFrom SummarizedExperiment colData SummarizedExperiment
#' @param object a SummarizedExperiment object
#'
#' @returns colData as data.frame
#' @export
#'
get_col_data <- function(
    object = NULL
) {
  data <- as.data.frame(
    SummarizedExperiment::colData(object)
  )
  return(data)
}

#' @export
get_colData <- function(...) get_col_data(...)
