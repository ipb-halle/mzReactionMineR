#' Align peaks across multiple samples
#'
#' This function aligns peaks across multiple samples based on m/z and retention time (rt) tolerances.
#' It is a simplified implementation of mzMine's Join Aligner algorithm.
#' mzmine.github.io/mzmine_documentation/module_docs/align_join_aligner/join_aligner.html
#'
#' @importFrom dplyr %>% select mutate bind_rows sym row_number filter arrange across between
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom rlang :=
#' @param input A SummarizedExperiment object or list of SummarizedExperiment
#' objects containing peak information in rowData columns for id, rt, and mz.
#' @param mz_tolerance Numeric of length 2,
#' the tolerance for m/z alignment window in absolute and relative (ppm) units.
#' @param rt_tolerance Numeric,
#' the tolerance for retention time alignment window.
#' @param mz_weight Numeric,
#' the weight for m/z in the alignment score calculation.
#' @param rt_weight Numeric,
#' the weight for retention time in the alignment score calculation.
#' @param ims Logical, whether ion mobility should be included in the alignment.
#' @param ion_mobility_tolerance Numeric, the tolerance for ion mobility alignment.
#' @param ion_mobility_weight Numeric, the weight for ion mobility in the alignment score.
#' @param id_col Character,
#' the name of the column containing peak IDs.
#' @param rt_col Character,
#' the name of the column containing retention times.
#' @param mz_col Character,
#' the name of the column containing m/z values.
#' @param ion_mobility_col Character,
#' the name of the column containing ion mobility values.
#' @return A SummarizedExperiment object with aligned rowData, joined colData,
#' and the assays shared by all input objects.
#' @export
#'

join_aligner <- function(
  input,
  mz_tolerance = c(0.01, 20),
  rt_tolerance = 0.5,
  mz_weight = 3,
  rt_weight = 1,
  id_col = "id",
  rt_col = "rt",
  mz_col = "mz",
  ims = FALSE,
  ion_mobility_tolerance = 0.1,
  ion_mobility_weight = 1,
  ion_mobility_col = "ion_mobility"
) {
  if (!is.logical(ims) || length(ims) != 1L || is.na(ims)) {
    stop("ims must be a single TRUE or FALSE value")
  }
  if (!is.character(c(id_col, rt_col, mz_col, ion_mobility_col)) ||
      anyNA(c(id_col, rt_col, mz_col, ion_mobility_col)) ||
      any(!nzchar(c(id_col, rt_col, mz_col, ion_mobility_col)))) {
    stop("Column names must be non-empty character strings")
  }
  if (anyDuplicated(c(id_col, rt_col, mz_col, if (ims) ion_mobility_col))) {
    stop("Column names must be unique")
  }
  if (!is.numeric(mz_tolerance) || length(mz_tolerance) != 2L ||
      any(!is.finite(mz_tolerance)) || any(mz_tolerance < 0)) {
    stop("mz_tolerance must be a numeric vector of two non-negative finite values")
  }
  tolerances <- c(rt_tolerance, if (ims) ion_mobility_tolerance)
  weights <- c(mz_weight, rt_weight, if (ims) ion_mobility_weight)
  if (any(!is.numeric(tolerances)) || any(lengths(lapply(tolerances, length)) != 1L) ||
      any(!is.finite(tolerances)) || any(tolerances <= 0)) {
    stop("Alignment tolerances must be positive finite numeric values")
  }
  if (!is.numeric(weights) || any(vapply(weights, length, integer(1)) != 1L) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("Alignment weights must be non-negative finite numeric values")
  }

  required <- c(id_col, rt_col, mz_col, if (ims) ion_mobility_col)
  se_list <- if (methods::is(input, "SummarizedExperiment")) {
    list(input)
  } else if (is.list(input) && !is.data.frame(input)) {
    input
  } else {
    stop("Input must be a SummarizedExperiment object or a list of SummarizedExperiment objects")
  }

  if (!length(se_list) || any(!vapply(se_list, methods::is, logical(1), class2 = "SummarizedExperiment"))) {
    stop("Input must be a SummarizedExperiment object or a list of SummarizedExperiment objects")
  }

  assay_names <- lapply(se_list, function(object) names(SummarizedExperiment::assays(object)))
  shared_assays <- Reduce(intersect, assay_names)
  if (!length(shared_assays)) {
    stop("Input SummarizedExperiment objects must share at least one assay")
  }

  samples <- unlist(
    Map(se_to_samples, se_list, seq_along(se_list), MoreArgs = list(required = required)),
    recursive = FALSE
  )

  if (length(samples) < 2L) {
    stop("At least two samples must be supplied")
  }
  if (anyNA(names(samples)) || any(!nzchar(names(samples))) || anyDuplicated(names(samples))) {
    stop("Input sample names must be unique and non-empty")
  }

  sample_sources <- rep(seq_along(se_list), vapply(se_list, ncol, integer(1)))
  sample_columns <- unlist(lapply(se_list, function(object) seq_len(ncol(object))), use.names = FALSE)
  names(sample_sources) <- names(samples)
  names(sample_columns) <- names(samples)

  output_col_data <- lapply(names(samples), function(sample_name) {
    source_id <- sample_sources[[sample_name]]
    sample_col <- sample_columns[[sample_name]]
    sample_col_data <- as.data.frame(SummarizedExperiment::colData(se_list[[source_id]]))[
      sample_col,
      ,
      drop = FALSE
    ]
    rownames(sample_col_data) <- sample_name
    sample_col_data
  })
  names(output_col_data) <- names(samples)
  col_data_names <- unique(unlist(lapply(output_col_data, names), use.names = FALSE))
  output_col_data <- lapply(output_col_data, function(sample_col_data) {
    missing <- setdiff(col_data_names, names(sample_col_data))
    sample_col_data[missing] <- NA
    sample_col_data[col_data_names]
  })

  measurement_cols <- c("rt", "mz", if (ims) "ion_mobility")
  samples <- lapply(samples, function(sample) {
    sample <- sample[, c(required, ".source", ".row"), drop = FALSE]
    names(sample) <- c("id", "rt", "mz", if (ims) "ion_mobility", ".source", ".row")
    sample[measurement_cols] <- suppressWarnings(lapply(sample[measurement_cols], as.numeric))
    if (anyNA(sample[, measurement_cols, drop = FALSE])) {
      stop("Alignment measurements must be numeric and cannot contain missing values")
    }
    sample
  })

  if (any(vapply(samples, nrow, integer(1)) == 0L)) {
    stop("Input data frames must contain at least one peak")
  }

  samples <- samples[order(vapply(samples, nrow, integer(1)), decreasing = TRUE)]
  output_col_data <- if (length(col_data_names)) {
    do.call(rbind, output_col_data[names(samples)])
  } else {
    data.frame(row.names = names(samples))
  }
  master_list <- samples[[1L]]

  for (sample_id in seq_along(samples)) {
    sample <- samples[[sample_id]]
    sample_name <- names(samples)[sample_id]
    source_col <- paste0(sample_name, ".source")
    row_col <- paste0(sample_name, ".row")

    master_list[[sample_name]] <- NA
    master_list[[source_col]] <- NA_integer_
    master_list[[row_col]] <- NA_integer_
    master_list$score <- NA

    for (i in seq_len(nrow(master_list))) {
      master_mz <- master_list$mz[i]
      master_rt <- master_list$rt[i]
      mz_window <- max(mz_tolerance[1L], master_mz * mz_tolerance[2L] / 1e6)
      candidate_peaks <- sample %>% filter(
        between(mz, master_mz - mz_window, master_mz + mz_window),
        between(rt, master_rt - rt_tolerance, master_rt + rt_tolerance)
      )
      if (ims) {
        master_ims <- master_list$ion_mobility[i]
        candidate_peaks <- candidate_peaks %>% filter(
          between(ion_mobility,
                  master_ims - ion_mobility_tolerance,
                  master_ims + ion_mobility_tolerance)
        )
      }

      if (nrow(candidate_peaks) > 0) {
        score_terms <- list(
          (1 - abs(candidate_peaks$mz - master_mz) / mz_window) * mz_weight,
          (1 - abs(candidate_peaks$rt - master_rt) / rt_tolerance) * rt_weight
        )
        if (ims) {
          score_terms[[3L]] <- (1 - abs(candidate_peaks$ion_mobility - master_ims) /
            ion_mobility_tolerance) * ion_mobility_weight
        }
        scores <- rowSums(do.call(cbind, score_terms))
        best_match <- candidate_peaks[which.max(scores), ]
        master_list[i, "mz"] <- mean(c(master_mz, best_match$mz))
        master_list[i, "rt"] <- mean(c(master_rt, best_match$rt))
        if (ims) {
          master_list[i, "ion_mobility"] <- mean(c(master_ims, best_match$ion_mobility))
        }
        master_list[i, sample_name] <- best_match$id
        master_list[i, source_col] <- best_match$.source
        master_list[i, row_col] <- best_match$.row
        master_list[i, "score"] <- max(scores)
      }
    }

    matched <- !is.na(master_list[[sample_name]])
    for (peak_id in unique(master_list[[sample_name]][matched])) {
      matches <- which(master_list[[sample_name]] == peak_id)
      remove_matches <- setdiff(matches, matches[which.max(master_list$score[matches])])
      master_list[[sample_name]][remove_matches] <- NA
      master_list[[source_col]][remove_matches] <- NA_integer_
      master_list[[row_col]][remove_matches] <- NA_integer_
    }
    master_list$score <- NULL

    missing_features <- sample %>%
      filter(!id %in% master_list[[sample_name]])

    if (nrow(missing_features) > 0) {
      missing_features[[sample_name]] <- missing_features$id
      missing_features[[source_col]] <- missing_features$.source
      missing_features[[row_col]] <- missing_features$.row
      master_list <- bind_rows(master_list, missing_features) %>%
        arrange(rt, mz, !!!if (ims) list(ion_mobility) else list()) %>%
        mutate(id = row_number())
    }
  }

  feature_cols <- c("id", "rt", "mz", if (ims) "ion_mobility")
  output_row_data <- master_list[, feature_cols, drop = FALSE]
  names(output_row_data) <- c(id_col, rt_col, mz_col, if (ims) ion_mobility_col)

  output_assays <- lapply(shared_assays, function(assay_name) {
    aligned_assay <- matrix(
      NA,
      nrow = nrow(master_list),
      ncol = length(samples),
      dimnames = list(NULL, names(samples))
    )
    for (sample_name in names(samples)) {
      source_col <- paste0(sample_name, ".source")
      row_col <- paste0(sample_name, ".row")
      matched <- !is.na(master_list[[source_col]]) & !is.na(master_list[[row_col]])
      if (any(matched)) {
        source_id <- sample_sources[[sample_name]]
        sample_col <- sample_columns[[sample_name]]
        aligned_assay[matched, sample_name] <- SummarizedExperiment::assays(se_list[[source_id]])[[assay_name]][
          master_list[[row_col]][matched],
          sample_col
        ]
      }
    }
    aligned_assay
  })
  names(output_assays) <- shared_assays

  SummarizedExperiment::SummarizedExperiment(
    rowData = output_row_data,
    assays = output_assays,
    colData = output_col_data
  )
}

utils::globalVariables(c(".data", ".", "id", "mz", "rt", "ion_mobility"))
