#' mzmine_to_se
#'
#' Convert mzMine feature table to a SummarizedExperiment object. Assumes the
#'     "Export to CSV (modular)" was used.
#'
#' @importFrom utils read.csv
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @param path_to_file path to the mzmine feature table
#' @param sample_meta_data data.frame. sample meta data to become colData
#'     the first columnn has to be the sample names
#' @param filenames character. name of the column in the meta data that
#'     corresponds to the sample names in the feature table
#' @param assays character vector. what result type to use. (i.e. "height")
#'
#' @returns a summarizedExperiment Object
#' @export
#'

mzmine_to_se <- function(
    path_to_file,
    sample_meta_data = NULL,
    assays = c("area", "height", "mz", "rt"),
    filenames = "filename"
) {

  if (!is.data.frame(sample_meta_data)) {
    stop("`sample_meta_data` must be a data.frame")
  }
  if (!is.character(filenames) || length(filenames) != 1L ||
      !filenames %in% names(sample_meta_data)) {
    stop("`filenames` must name a column in `sample_meta_data`")
  }

  features <- tryCatch(
    read.csv(path_to_file, check.names = TRUE),
    error = function(error) stop("Could not read feature table: ", error$message)
  )
  if (!"id" %in% names(features)) {
    stop("Feature table must contain an `id` column")
  }
  features$id <- as.character(features$id)

  datafile_columns <- grep("^datafile[.]", names(features), value = TRUE)
  if (!length(datafile_columns)) {
    stop("Feature table does not contain datafile columns")
  }

  available_assays <- unique(sub(".*[.]", "", datafile_columns))
  available_assays <- setdiff(available_assays, c("min", "max"))
  if (is.null(assays)) {
    assays <- intersect(c("area", "height", "mz", "rt"), available_assays)
  } else if (!is.character(assays) || !length(assays) ||
             anyNA(assays) || any(!assays %in% available_assays)) {
    stop("`assays` must contain only assay names found in the feature table")
  }
  if (!length(assays)) {
    stop("Feature table does not contain any usable assays")
  }

  assays_list <- lapply(assays, function(assay) {
    suffix <- paste0(".", assay)
    columns <- datafile_columns[endsWith(datafile_columns, suffix)]
    values <- features[columns]
    names(values) <- make.names(sub("^datafile[.]", "", substr(
      columns,
      1L,
      nchar(columns) - nchar(suffix)
    )))
    as.matrix(values)
  })
  names(assays_list) <- assays

  samples <- colnames(assays_list[[1]])
  if (anyDuplicated(samples)) {
    stop("Feature table contains duplicate sample names")
  }
  if (any(vapply(assays_list, function(x) !identical(colnames(x), samples), logical(1)))) {
    stop("Sample names in the feature table are not equal across assays")
  }

  sample_meta_data[[filenames]] <- make.names(as.character(sample_meta_data[[filenames]]))
  if (anyDuplicated(sample_meta_data[[filenames]])) {
    stop("`sample_meta_data` contains duplicate sample names")
  }
  missing_metadata <- setdiff(samples, sample_meta_data[[filenames]])
  if (length(missing_metadata)) {
    warning(
      "Samples in the feature table are not in `sample_meta_data` and will be removed: ",
      paste(missing_metadata, collapse = ", ")
    )
    assays_list <- lapply(
      assays_list,
      function(x) x[, setdiff(samples, missing_metadata), drop = FALSE]
    )
    samples <- colnames(assays_list[[1]])
  }
  removed_metadata <- setdiff(sample_meta_data[[filenames]], samples)
  if (length(removed_metadata)) {
    warning(
      "Samples in `sample_meta_data` are not in the feature table and will be removed: ",
      paste(removed_metadata, collapse = ", ")
    )
    sample_meta_data <- sample_meta_data[
      sample_meta_data[[filenames]] %in% samples,
      ,
      drop = FALSE
    ]
  }
  if (!length(samples)) {
    stop("No feature-table samples are present in `sample_meta_data`")
  }

  sample_meta_data <- sample_meta_data[match(samples, sample_meta_data[[filenames]]), , drop = FALSE]
  rowData_cols <- setdiff(names(features), datafile_columns)

  SummarizedExperiment(
    rowData = features[,rowData_cols],
    assays = assays_list,
    colData = sample_meta_data
  )

}


