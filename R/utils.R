# A set of helper functions used in the main exported functions.

object_argument_name <- function(object_expression) {
  paste(deparse(object_expression), collapse = "")
}


# function to find the closest id of a given rt/mz pair ------------------------

#' get_id
#'
#' Finds the id of the closest feature to the provided rt and mz values, within
#' the specified ranges. If rt is not provided, it will find the closest feature
#'  based on mz values only.
#'
#' @param object a SummarizedExperiment object
#' @param rt Numeric. Retention time value
#' @param mz Numeric. mz value
#' @param mz_range Numeric vector of length 2. The range of mz values to consider when finding the closest feature.
#' @param rt_range Numeric vector of length 2. The range of rt values to consider when finding the closest feature.umer
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns the id of the closest feature to the provided rt and mz values
#'
get_id <- function(
    object,
    rt = NULL,
    mz = NULL,
    mz_range = c(0, Inf),
    rt_range = c(0, Inf),
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {
  object_name <- object_argument_name(substitute(object))

  if(length(mz) != 1L || !is.numeric(mz) || is.na(mz) || !is.finite(mz)) {
    stop("'mz' must be one finite number.", call. = FALSE)
  }
  if(!is.null(rt) &&
     (length(rt) != 1L || !is.numeric(rt) || is.na(rt) || !is.finite(rt))) {
    stop("'rt' must be NULL or one finite number.", call. = FALSE)
  }
  if(length(mz_range) != 2L || !is.numeric(mz_range) || anyNA(mz_range) ||
     mz_range[1] > mz_range[2]) {
    stop("'mz_range' must contain two ordered numbers.", call. = FALSE)
  }
  if(length(rt_range) != 2L || !is.numeric(rt_range) || anyNA(rt_range) ||
     rt_range[1] > rt_range[2]) {
    stop("'rt_range' must contain two ordered numbers.", call. = FALSE)
  }

  row_data <- get_rowData(object)
  required_columns <- c(id_col, mz_col, if(!is.null(rt)) rt_col)
  if(any(!required_columns %in% names(row_data))) {
    stop(
      "Required rowData(", object_name, ") columns are missing: ",
      paste(setdiff(required_columns, names(row_data)), collapse = ", "),
      call. = FALSE
    )
  }

  if(is.null(rt)) {

    # find id
    result_id <- row_data %>%
      filter(
        between(.data[[mz_col]], mz_range[1], mz_range[2])
      ) %>%
      mutate(
        mz_diff = abs(.data[[mz_col]]-mz)
      ) %>%
      slice_min(
        across(c(mz_diff)), n=1
      ) %>%
      pull(.data[[id_col]])

  } else {

    # find id
    result_id <- row_data %>%
      filter(
        between(.data[[mz_col]], mz_range[1], mz_range[2]),
        between(.data[[rt_col]], rt_range[1], rt_range[2])
      ) %>%
      mutate(
        mz_diff = abs(.data[[mz_col]]-mz),
        rt_diff = abs(.data[[rt_col]]-rt)
      ) %>%
      slice_min(
        across(c(rt_diff,mz_diff)), n=1
      ) %>%
      pull(.data[[id_col]])

  }

  return(result_id)

}

# function to get an intensity vector from a specified assay and id ------------

get_intensities_id <- function(
    object,
    assay = NULL,
    id,
    id_col = "id",
    object_name = object_argument_name(substitute(object))
) {

  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }
  row_data <- rowData(object)
  if(!id_col %in% names(row_data)) {
    stop("'id_col' is missing from rowData(", object_name, ").", call. = FALSE)
  }
  matching_rows <- which(row_data[[id_col]] == id)
  if(length(matching_rows) != 1L) {
    stop("'id' must identify exactly one feature in rowData(", object_name, ").", call. = FALSE)
  }

  values <- assays(object[matching_rows, ])[[assay]]
  return(as.numeric(values))

}

# Resolve and validate the internal standard used by normalizeIS and plotIS.
resolve_is_id <- function(
    object,
    id,
    rt,
    mz,
    mz_tolerance,
    rt_tolerance,
    id_col,
    rt_col,
    mz_col,
    object_name = object_argument_name(substitute(object))
) {
  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if(length(id) != 1L || !is.character(id) || is.na(id) || !nzchar(id)) {
    stop("'id' must be one non-empty string or 'none'.", call. = FALSE)
  }
  columns <- c(id_col, rt_col, mz_col)
  if(any(vapply(columns, function(x) {
    length(x) != 1L || !is.character(x) || is.na(x) || !nzchar(x)
  }, logical(1)))) {
    stop("Feature column names must be single non-empty strings.", call. = FALSE)
  }
  if(!is.numeric(mz_tolerance) || length(mz_tolerance) != 2L ||
     anyNA(mz_tolerance) || any(!is.finite(mz_tolerance)) ||
     any(mz_tolerance < 0)) {
    stop("'mz_tolerance' must contain two non-negative finite numbers.", call. = FALSE)
  }
  if(length(rt_tolerance) != 1L || !is.numeric(rt_tolerance) ||
     is.na(rt_tolerance) || !is.finite(rt_tolerance) || rt_tolerance < 0) {
    stop("'rt_tolerance' must be one non-negative finite number.", call. = FALSE)
  }

  row_data <- get_rowData(object)
  if(!id_col %in% names(row_data)) {
    stop("'id_col' is missing from rowData(", object_name, ").", call. = FALSE)
  }
  if(anyNA(row_data[[id_col]]) || anyDuplicated(row_data[[id_col]])) {
    stop("'id_col' must contain unique, non-missing feature identifiers in rowData(", object_name, ").", call. = FALSE)
  }

  if(id == "none") {
    if(length(rt) != 1L || !is.numeric(rt) || is.na(rt) || !is.finite(rt) ||
       length(mz) != 1L || !is.numeric(mz) || is.na(mz) || !is.finite(mz)) {
      stop("'rt' and 'mz' must each be one finite number when 'id' is 'none'.", call. = FALSE)
    }
    result_id <- get_id(
      object = object,
      rt = rt,
      mz = mz,
      mz_range = calc_mz_range(mz, mz_tolerance),
      rt_range = calc_rt_range(rt, rt_tolerance),
      id_col = id_col,
      rt_col = rt_col,
      mz_col = mz_col
    )
    if(length(result_id) != 1L || is.na(result_id)) {
      stop("No feature matching the supplied internal-standard coordinates was found.", call. = FALSE)
    }
    return(result_id)
  }

  if(!id %in% row_data[[id_col]]) {
    stop("The internal-standard 'id' was not found in rowData(", object_name, ").", call. = FALSE)
  }
  id
}

# Resolve the internal standard and return its id and intensities together.
get_is_data <- function(
    object,
    assay,
    rt,
    mz,
    id,
    mz_tolerance,
    rt_tolerance,
    id_col,
    rt_col,
    mz_col,
    object_name = object_argument_name(substitute(object))
) {
  is_id <- resolve_is_id(
    object, id, rt, mz, mz_tolerance, rt_tolerance,
    id_col, rt_col, mz_col, object_name = object_name
  )
  is_intensities <- get_intensities_id(object, assay, is_id, id_col, object_name = object_name)

  if(anyNA(is_intensities)) {
    stop("Internal-standard intensities contain missing values.", call. = FALSE)
  }
  if(any(!is.finite(is_intensities)) || any(is_intensities <= 0)) {
    stop("Internal-standard intensities must be positive and finite.", call. = FALSE)
  }

  list(id = is_id, intensities = is_intensities)
}

# functions to divide an assay by a vector of intensities ----------------------

divide_by_feature <- function(
    object,
    assay = NULL,
    vector = NULL,
    new_assay_name = NULL,
    object_name = object_argument_name(substitute(object))
) {

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.", call. = FALSE)
  }
  if(length(new_assay_name) != 1L || !is.character(new_assay_name) ||
     is.na(new_assay_name) || !nzchar(new_assay_name)) {
    stop("'new_assay_name' must be one non-empty string.", call. = FALSE)
  }
  if(length(vector) != ncol(object) || anyNA(vector) ||
     any(!is.finite(vector)) || any(vector == 0)) {
    stop("'vector' must contain one finite, non-zero value per sample.", call. = FALSE)
  }
  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.", call. = FALSE)
  }

  new_object <- object
  assays(new_object, withDimnames = FALSE)[[new_assay_name]] <- t(apply(
    t(assays(new_object)[[assay]]),2,function(col) {
      col/vector
    }))

  return(new_object)

}

# function to calculate mz range based on absolute tolerance and ppm -----------

calc_mz_range <- function(
    mz,
    tolerance = c(0.005, 5)
) {

  mz_tolerance <- max(tolerance[1], (tolerance[2] * (mz/1e6)))
  mz_range <- c(mz - mz_tolerance, mz + mz_tolerance)
  return(mz_range)

}

# function to calculate rt range based on absolute tolerance and ppm -----------

calc_rt_range <- function(
    rt,
    tolerance = 0.1
) {

  rt_range <- c(rt - tolerance, rt + tolerance)
  return(rt_range)

}

# function to calculate cosine distance between rows of a matrix ---------------

calc_cosine_distance <- function(mat) {

  similarity <- (mat %*% t(mat)) / # calculate pairwise dotproduct
    # divide by the product of the length of the vectors
    outer(sqrt(rowSums(mat**2)),sqrt(rowSums(mat**2)))

  distance <- 1 - similarity

  return(distance)

}

# function to make a knn graph from a distance matrix --------------------------

#' make_knn_graph
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @param distance_mat distance matrix
#' @param k number of neighbors
#'
#' @returns an igraph object representing the knn graph
#'
make_knn_graph <- function(distance_mat, k) {

  # make 0 matrix to store knn graph
  knn_graph <- matrix(0, nrow = nrow(distance_mat), ncol = ncol(distance_mat))

  # loop over each row and find the k nearest neighbors of each point, excluding itself
  for(i in 1:nrow(knn_graph)) {

    knn_graph[i, order(distance_mat[i, ], decreasing = FALSE)[2:(k+1)]] <- 1

  }

  # make symmetric
  knn_graph <- pmax(knn_graph, t(knn_graph))

  # make igraph object from adjacency matrix
  knn_graph <- graph_from_adjacency_matrix(knn_graph, mode = "undirected" )

  return(knn_graph)

}
