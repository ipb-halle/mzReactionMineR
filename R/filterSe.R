#' filterSE
#'
#' A function that filters a SummarizedExperiment object based on various
#'     measures.
#'
#' @importFrom dplyr %>% mutate filter group_by summarize ungroup pull n inner_join
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
#' @param min_abundance Numeric \[0,1\]. Minimum value in assay for a feature. Values
#'     below are treated as absent. Default is 0.
#' @param min_pct Numeric. Minimum ratio of samples a feature must be
#'     present in. If group_col is specified, the ratio is calculated for
#'     for each group. If min_n is specified, the larger value will be used.
#'     This is applied after filtering for min_abundance.
#' @param min_n Integer. Minimum number of samples a feature must be
#'     present in. If group_col is specified, the percentage is calculated for
#'     for each group. If min_pct is specified, this is ignored, the larger value
#'     will be used.
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
#' @param mobility_range Numeric or NULL. A vector of length 2 with the range
#'     of kept ion mobility values. If NULL, ion mobility is not used for
#'     filtering.
#' @param specific_mz Numeric vector or "none". remove rows that do not have a
#'     specific m/z value. This is applied after filtering for min_abundance.
#' @param id_col Character. The respective column name in the rowData.
#' @param rt_col Character. The respective column name in the rowData.
#' @param mz_col Character. The respective column name in the rowData.
#' @param ion_mobility_col Character. The respective column name in the rowData.
#'
#' @returns A filtered summarizeExperiment object
#' @export
filterSe <- function(
  object = NULL,
  assay = NULL,
  sample_col = "filename",
  group_col = "none",
  not_in = "none",
  min_abundance = 0,
  min_pct = 0.8,
  min_n = 1L,
  min_pct_total = min_pct,
  min_n_total = min_n,
  min_pct_group = min_pct,
  min_n_group = min_n,
  min_count_rule = c("or", "and"),
  mz_range = c(0,Inf),
  rt_range = c(0,Inf),
  mobility_range = NULL,
  specific_mz = "none",
  id_col = "id",
  rt_col = "rt",
  mz_col = "mz",
  ion_mobility_col = "ion_mobility"
) {
  object_name <- object_argument_name(substitute(object))

  min_count_rule <- match.arg(min_count_rule)

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }

  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }

  character_arguments <- list(
    sample_col = sample_col,
    id_col = id_col,
    rt_col = rt_col,
    mz_col = mz_col,
    ion_mobility_col = ion_mobility_col
  )

  invalid_character_arguments <- names(character_arguments)[
    !vapply(character_arguments, function(x) {
      length(x) == 1L && is.character(x) && !is.na(x) && nzchar(x)
    }, logical(1))
  ]

  if(length(invalid_character_arguments) > 0L) {
    stop(
      "The following arguments must be single non-empty strings: ",
      paste(invalid_character_arguments, collapse = ", "),
      call. = FALSE
    )
  }

  if(length(group_col) != 1L || !is.character(group_col) ||
     is.na(group_col) || !nzchar(group_col)) {
    stop("'group_col' must be a single non-empty string or 'none'.", call. = FALSE)
  }

  if(!is.character(not_in) || anyNA(not_in) || any(!nzchar(not_in))) {
    stop("'not_in' must contain non-empty strings.", call. = FALSE)
  }

  if(group_col == "none" && all(not_in != "none")) {
    stop("'not_in' can only be used when 'group_col' names a colData(", object_name, ") column.", call. = FALSE)
  }

  if(length(min_abundance) != 1L || !is.numeric(min_abundance) ||
     is.na(min_abundance) || !is.finite(min_abundance) || min_abundance < 0) {
    stop("'min_abundance' must be a single non-negative finite number.", call. = FALSE)
  }

  pct_arguments <- list(
    min_pct = min_pct,
    min_pct_total = min_pct_total,
    min_pct_group = min_pct_group
  )

  invalid_pct_arguments <- names(pct_arguments)[
    !vapply(pct_arguments, function(x) {
      length(x) == 1L && is.numeric(x) && !is.na(x) && x >= 0 && x <= 1
    }, logical(1))
  ]

  if(length(invalid_pct_arguments) > 0L) {
    stop(
      paste(invalid_pct_arguments, collapse = ", "),
      " must be single numbers between 0 and 1.",
      call. = FALSE
    )
  }

  n_arguments <- list(
    min_n = min_n,
    min_n_total = min_n_total,
    min_n_group = min_n_group
  )

  invalid_n_arguments <- names(n_arguments)[
    !vapply(n_arguments, function(x) {
      length(x) == 1L && is.numeric(x) && !is.na(x) && is.finite(x) &&
        x >= 0 && x == floor(x)
    }, logical(1))
  ]

  if(length(invalid_n_arguments) > 0L) {
    stop(
      paste(invalid_n_arguments, collapse = ", "),
      " must be single non-negative integers.",
      call. = FALSE
    )
  }

  range_arguments <- list(mz_range = mz_range, rt_range = rt_range)
  if(!is.null(mobility_range)) {
    range_arguments$mobility_range <- mobility_range
  }

  invalid_range_arguments <- names(range_arguments)[
    !vapply(range_arguments, function(x) {
      length(x) == 2L && is.numeric(x) && !anyNA(x) && x[1] < x[2]
    }, logical(1))
  ]

  if(length(invalid_range_arguments) > 0L) {
    stop(
      "The following arguments must be numeric length-2 ranges with the lower value first: ",
      paste(invalid_range_arguments, collapse = ", "),
      call. = FALSE
    )
  }

  if(!identical(specific_mz, "none") &&
     (!is.numeric(specific_mz) || anyNA(specific_mz))) {
    stop("'specific_mz' must be numeric values or 'none'.", call. = FALSE)
  }

  row_data <- as.data.frame(rowData(object))
  sample_data <- get_colData(object)

  missing_row_columns <- setdiff(
    c(id_col, rt_col, mz_col, if(!is.null(mobility_range)) ion_mobility_col),
    colnames(row_data)
  )
  if(length(missing_row_columns) > 0L) {
    stop(
      "The following columns are missing from rowData(", object_name, "): ",
      paste(missing_row_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_col_columns <- setdiff(
    c(sample_col, if(group_col != "none") group_col),
    colnames(sample_data)
  )
  if(length(missing_col_columns) > 0L) {
    stop(
      "The following columns are missing from colData(", object_name, "): ",
      paste(missing_col_columns, collapse = ", "),
      call. = FALSE
    )
  }

  missing_assay_samples <- setdiff(colnames(assays(object)[[assay]]), sample_data[[sample_col]])
  if(length(missing_assay_samples) > 0L) {
    stop(
      "The following assay column names are missing from colData(", object_name, ")[[\'", sample_col, "\']]: ",
      paste(missing_assay_samples, collapse = ", "),
      call. = FALSE
    )
  }

  if(group_col != "none" && all(not_in != "none")) {
    missing_not_in <- setdiff(not_in, sample_data[[group_col]])
    if(length(missing_not_in) > 0L) {
      stop(
        "The following 'not_in' values are missing from colData(", object_name, ")[[\'", group_col, "\']]: ",
        paste(missing_not_in, collapse = ", "),
        call. = FALSE
      )
    }
  }

  data_long <- cbind(
    row_data,
    as.data.frame(assays(object)[[assay]]) %>%
      replace(is.na(.data),0)
  ) %>%
    pivot_longer(
      cols = -all_of(names(row_data)),
      names_to = sample_col,
      values_to = "Value"
    ) %>%
    mutate(
      Value = as.integer(.data$Value > min_abundance)
    ) %>%
    inner_join(
      sample_data,
      by = sample_col
    ) %>%
    filter(
      .data[[mz_col]] > mz_range[1],
      .data[[mz_col]] < mz_range[2],
      .data[[rt_col]] > rt_range[1],
      .data[[rt_col]] < rt_range[2]
    )

  if(!is.null(mobility_range)) {
    data_long <- data_long %>%
      filter(
        .data[[ion_mobility_col]] > mobility_range[1],
        .data[[ion_mobility_col]] < mobility_range[2]
      )
  }

  if(group_col != "none" && all(not_in != "none")) {
    kept_ids <- data_long %>%
      filter(.data[[group_col]] %in% not_in) %>%
      group_by(.data[[id_col]]) %>%
      summarize(Value = sum(.data$Value, na.rm = TRUE)) %>%
      filter(.data$Value == 0) %>%
      pull(.data[[id_col]]) %>%
      unique()

    data_long <- data_long %>%
      filter(.data[[id_col]] %in% kept_ids)

  }

  total_ids <- data_long %>%
    group_by(.data[[id_col]]) %>%
    summarize(
      Value = sum(.data$Value, na.rm = TRUE),
      cutoff = max(floor(n()*min_pct_total), min_n_total)
    ) %>%
    filter(.data$Value >= .data$cutoff) %>%
    ungroup() %>%
    pull(.data[[id_col]]) %>%
    unique()

  if(group_col == "none") {
    ids <- total_ids
  } else {
    group_ids <- data_long %>%
      group_by(.data[[id_col]], .data[[group_col]]) %>%
      summarize(
        Value = sum(.data$Value, na.rm = TRUE),
        cutoff = max(floor(n()*min_pct_group), min_n_group)
      ) %>%
      filter(.data$Value >= .data$cutoff) %>%
      ungroup() %>%
      pull(.data[[id_col]]) %>%
      unique()

    if(min_count_rule == "and") {
      ids <- intersect(group_ids, total_ids)
    } else {
      ids <- union(group_ids, total_ids)
    }
  }

  object[rowData(object)[[id_col]] %in% ids,]

}

utils::globalVariables(c(".",".data","Value"))
