#' Align peaks across multiple samples using ion mobility.
#'
#' This is a convenience wrapper around [join_aligner()] with ion mobility
#' alignment enabled.
#'
#' @param input A SummarizedExperiment object or list of SummarizedExperiment
#'     objects containing peak information in rowData columns for id, rt, mz,
#'     and ion mobility.
#' @param mz_tolerance Numeric vector of length 2,
#'     the tolerance for m/z alignment window (absolute, relative).
#' @param rt_tolerance Numeric,
#'     the tolerance for retention time alignment window.
#' @param ion_mobility_tolerance Numeric,
#'     the tolerance for ion mobility alignment window.
#' @param mz_weight Numeric,
#'     the weight for m/z in the alignment score calculation.
#' @param rt_weight Numeric,
#'     the weight for retention time in the alignment score calculation.
#' @param ion_mobility_weight Numeric,
#'     the weight for retention time in the alignment score calculation.
#' @param id_col Character,
#'     the name of the column containing peak IDs.
#' @param rt_col Character,
#'     the name of the column containing retention times.
#' @param mz_col Character,
#'     the name of the column containing m/z values.
#' @param ion_mobility_col Character,
#'     the name of the column containing ion mobility values.
#' @return A SummarizedExperiment object with aligned rowData, joined colData,
#'     and the assays shared by all input objects.
#' @export
#'

join_aligner_ims <- function(
  input,
  mz_tolerance = c(0.05, 20),
  rt_tolerance = 0.2,
  ion_mobility_tolerance = 0.1,
  mz_weight = 3,
  rt_weight = 1,
  ion_mobility_weight = 1,
  id_col = "id",
  rt_col = "rt",
  mz_col = "mz",
  ion_mobility_col = "ion_mobility"
) {
  join_aligner(
    input = input,
    mz_tolerance = mz_tolerance,
    rt_tolerance = rt_tolerance,
    mz_weight = mz_weight,
    rt_weight = rt_weight,
    id_col = id_col,
    rt_col = rt_col,
    mz_col = mz_col,
    ims = TRUE,
    ion_mobility_tolerance = ion_mobility_tolerance,
    ion_mobility_weight = ion_mobility_weight,
    ion_mobility_col = ion_mobility_col
  )
}
