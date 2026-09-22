#' filterSE
#'
#' A function that filters a SummarizedExperiment object based on various
#'     measures.
#'
#' @importFrom dplyr %>% mutate filter group_by summarize ungroup pull sym n inner_join
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assays colData
#' @importFrom tidyselect all_of
#'
#' @param object A SummarizedExperiment object
#' @param assay Assay name to filter on
#' @param sample_col Character. Name of the column in the colData that
#'     that corresponds to the assay columns
#' @param group_col Character. Either "none" or a columnname in the colData
#' @param not_in Character. Either "none" or a value in the group_col.
#'     Specify, if an id should not be present in a group.
#' @param min_abundance Numeric. Minimum abundance for a feature to be present
#'     in a sample
#' @param min_pct Numeric. Minimum percentage of samples a feature must be
#'     present in. If group_col is specified, the percentage is calculated for
#'     for each group.
#' @param min_n Integer. Minimum number of samples a feature must be present in.
#' @param min_pct_total Numeric. Minimum ratio of all samples a feature must be
#'     present in. Defaults to min_pct.
#' @param min_n_total Integer. Minimum number of all samples a feature must be
#'     present in. Defaults to min_n.
#' @param min_pct_group Numeric. Minimum ratio of samples per group a feature
#'     must be present in. Defaults to min_pct.
#' @param min_n_group Integer. Minimum number of samples per group a feature
#'     must be present in. Defaults to min_n.
#' @param min_count_rule Character. Either "or" or "and". Controls whether
#'     features must pass the min_n/min_pct threshold per group and total, or
#'     per group or total. If group_col is "none", only the total threshold is
#'     used.
#' @param mz_range Numeric. A vector of length 2. The range of kept m/z values.
#' @param rt_range Numeric. A vector of length 2. The range of kept rt values.
#' @param mobility_range Numeric. A vector of length 2.
#'     The range of kept ion mobility values.
#' @param id_col Character. The respective column name in the rowData.
#' @param rt_col Character. The respective column name in the rowData.
#' @param mz_col Character. The respective column name in the rowData.
#' @param ion_mobility_col Character. The respective column name in the rowData.
#'
#' @returns A filtered summarizeExperiment object
#' @export
filter_se_ims <- function(
    object = NULL,
    assay = NULL,
    sample_col = "filename",
    group_col = "none",
    not_in = "none",
    min_abundance = 0,
    min_pct = 0.8,
    min_n = 0L,
    min_pct_total = min_pct,
    min_n_total = min_n,
    min_pct_group = min_pct,
    min_n_group = min_n,
    min_count_rule = c("or", "and"),
    mz_range = c(0,Inf),
    rt_range = c(0,Inf),
    mobility_range = c(0,Inf),
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz",
    ion_mobility_col = "ion_mobility"
) {
  filter_se(
    object = object,
    assay = assay,
    sample_col = sample_col,
    group_col = group_col,
    not_in = not_in,
    min_abundance = min_abundance,
    min_pct = min_pct,
    min_n = min_n,
    min_pct_total = min_pct_total,
    min_n_total = min_n_total,
    min_pct_group = min_pct_group,
    min_n_group = min_n_group,
    min_count_rule = min_count_rule,
    mz_range = mz_range,
    rt_range = rt_range,
    mobility_range = mobility_range,
    id_col = id_col,
    rt_col = rt_col,
    mz_col = mz_col,
    ion_mobility_col = ion_mobility_col
  )

}

#' @export
filterSe_ims <- function(...) filter_se_ims(...)

utils::globalVariables(c(".",".data","Value"))
