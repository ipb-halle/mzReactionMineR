#' Convert MetaboScape output to a SummarizedExperiment
#'
#' This legacy wrapper delegates to [feature_table_to_se()].
#'
#' @param path_to_file chr. Path to the csv file containing the metaboscape
#'        output.
#' @param sample_meta_data data.frame. Sample meta data used as colData.
#' @param assay_name chr. Name of the assay to be created from the data. Default is "intensity".
#' @param filename_column chr. The name of the column in the sample_meta_data
#'        that contains the column names of the metaboscape output. Default is "filename".
#'
#' @returns a SummarizedExperiment Object
#' @export
#'
metaboscape_to_se <- function(
    path_to_file,
    sample_meta_data,
    assay_name = "intensity",
    filename_column = "filename"
) {
  feature_table_to_se(
    path_to_file = path_to_file,
    sample_meta_data = sample_meta_data,
    assay_name = assay_name,
    filename_column = filename_column
  )
}
