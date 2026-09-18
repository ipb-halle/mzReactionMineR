#' normalizeIS
#'
#' A function to normalize the intensities of an assay in a SummarizedExperiment
#'  object based on the intensities of an internal standard. The internal
#'  standard can be specified by its id, or by its rt and mz values. The
#'  function will find the closest feature to the provided rt and mz values, and
#'  use its intensities for normalization.
#'
#'  Normalization is defined as:
#'
#'    \deqn{\hat{Y_ij} = Y_ij\frac{\overline{S}{S_j}}}
#'
#'  Where \eqn{\hat{Y_ij}} is the normalized value of peak i in sample j,
#'  \eqn{Y_ij} is the original value  of peak i in sample j, \eqn{S_j} is the
#'  internal standard intensity in sample j and \eqn{\overline{S}} is the mean
#'  or median of the internal standard across samples.
#'
#' @importFrom dplyr %>% select mutate filter slice_min pull between across
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData assays<-
#' @importFrom stats median
#'
#' @param object A SummarizedExperiment object
#' @param assay Character. The name of the assay to be normalized.
#' @param rt Numeric. The theoretical retention time of the internal standard.
#' @param mz Numeric. The theoretical m/z of the internal standard.
#' @param id Character. The id of the internal standard. If "none", the function
#'        will find the closest feature to the provided rt and mz values.
#'        Default is "none".
#' @param type Either "mean" oder "median". How to calculate the reference IS
#'        value for normalization. Default is "median".
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerance for finding the internal standard.
#' @param rt_tolerance Numeric. Absolute retention time tolerance for finding
#'        the internal standard.
#' @param new_assay_name Character. The name of the new assay that will contain
#'        the normalized intensities. Default is "is_normalized".
#' @param remove Logical. Whether to remove the internal standard from the
#'        object after normalization. Default is TRUE.
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns A SummarizedExperiment object with a new assay containing the
#'        normalized intensities.
#' @export
normalize_is <- function(
    object,
    assay = NULL,
    rt = NULL,
    mz = NULL,
    id = "none",
    type = "median",
    mz_tolerance = c(0.005, 10),
    rt_tolerance = 0.1,
    new_assay_name = "is_normalized",
    remove = TRUE,
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {
  object_name <- object_argument_name(substitute(object))
  type <- match.arg(type, c("mean", "median"))
  is_data <- get_is_data(
    object, assay, rt, mz, id, mz_tolerance, rt_tolerance,
    id_col, rt_col, mz_col, object_name = object_name
  )

  correction_factor <- switch(
    type,
    mean = is_data$intensities/mean(is_data$intensities),
    median = is_data$intensities/median(is_data$intensities)
  )

  new_object <- divide_by_feature(
    object = object,
    assay = assay,
    vector = correction_factor,
    new_assay_name = new_assay_name,
    object_name = object_name
    )

  if(remove) {
    new_object <- new_object[-which(rowData(new_object)[[id_col]] == is_data$id), ]
  }

  return(new_object)

}

normalizeIS <- function(...) normalize_is(...)

utils::globalVariables(c(".data", ".", "mz_diff", "rt_diff"))
