#' knn_clustering_samples
#'
#' A function that performs clustering (community detection) on samples in a
#' SummarizedExperiment object based on a k-nearest neighbor graph.
#'
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom MatrixGenerics rowVars rowMaxs
#' @importFrom igraph cluster_louvain membership
#' @importFrom stats prcomp dist
#'
#' @param object a SummarizedExperiment object
#' @param assay Character. The name of the assay to be used for clustering.
#' @param filter_type Character. The method to select features for clustering.
#'        Must be either "intsnity" or "variance". Default is "area".
#' @param n_top Integer. The number of features to select based on
#'        filter_type.
#' @param calc_PCA Logical. Wether to perform PCA before clustering.
#'        Default is TRUE.
#' @param min_PC Integer. The minimum number of PCs to retain. Default is 5.
#' @param PC_var Numeric. The minimum cumulative variance that the retained PCs should explain.
#' @param k Integer. The number of neighbors to use for the k-nearest neighbor graph. Default is 5.
#' @param resolution Numeric. Resolution parameter for the Louvain algorithm. Default is 1.
#' @param metric Either "euclidean", "manhattan" or "cosine". The distance
#'        metric to use for calculating the distance matrix. Default is "cosine".
#' @param scale Logical. Whether to scale the data before calculating the
#'        distance matrix. Default is TRUE.
#' @param id_col Character. The name of the column in rowData(object) that
#'        contains the feature ids. Default is "id".
#'
#' @returns a data.frame with two columns: "sample" and "cluster".
#' @export
#'
knn_clustering_samples <- function(
  object,
  assay,
  filter_type = "intensity",
  n_top = nrow(object),
  calc_PCA = TRUE,
  min_PC = 1,
  PC_var = 0.8,
  k = 5,
  resolution = 1,
  metric = "cosine",
  scale = TRUE,
  id_col = "id"
) {

  # basic validation only; preserve original behavior for valid inputs
  if (!inherits(object, "SummarizedExperiment")) {
    stop("'object' must be a SummarizedExperiment.", call. = FALSE)
  }
  if (length(assay) != 1L || !is.character(assay) ||
      !assay %in% names(assays(object))) {
    stop("'assay' must name an assay in 'object'.", call. = FALSE)
  }
  if (length(id_col) != 1L || !is.character(id_col) ||
      !id_col %in% names(rowData(object))) {
    stop("'id_col' must name a column in rowData(object).", call. = FALSE)
  }
  if (length(filter_type) != 1L || !filter_type %in% c("intensity", "variance")) {
    stop("Invalid value for 'filter_type'. Must be either 'intensity' or 'variance'.", call. = FALSE)
  }
  if (length(n_top) != 1L || !is.numeric(n_top) || n_top < 1 ||
      n_top != as.integer(n_top) || n_top > nrow(object)) {
    stop("'n_top' must be an integer between 1 and the number of features.", call. = FALSE)
  }
  if (!is.logical(calc_PCA) || length(calc_PCA) != 1L || is.na(calc_PCA) ||
      !is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'calc_PCA' and 'scale' must each be TRUE or FALSE.", call. = FALSE)
  }
  if (length(min_PC) != 1L || !is.numeric(min_PC) || min_PC < 1 ||
      min_PC != as.integer(min_PC)) {
    stop("'min_PC' must be a positive integer.", call. = FALSE)
  }
  if (length(PC_var) != 1L || !is.numeric(PC_var) || !is.finite(PC_var) ||
      PC_var <= 0 || PC_var > 1) {
    stop("'PC_var' must be greater than 0 and at most 1.", call. = FALSE)
  }
  if (length(k) != 1L || !is.numeric(k) || k < 1 || k != as.integer(k) ||
      k >= ncol(object)) {
    stop("'k' must be an integer between 1 and one less than the number of samples.", call. = FALSE)
  }
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      !is.finite(resolution) || resolution <= 0) {
    stop("'resolution' must be one positive finite number.", call. = FALSE)
  }
  if (length(metric) != 1L || !metric %in% c("euclidean", "manhattan", "cosine")) {
    stop("Invalid value for 'metric'. Must be either 'euclidean', 'manhattan' or 'cosine'.", call. = FALSE)
  }

  print(paste0("Select top, ", n_top, " features based on", filter_type))

  if (filter_type == "intensity") {
    top_ids <- get_row_data(object)[
      order(rowMaxs(assays(object)[[assay]]), decreasing = TRUE)[1:n_top],
      , drop = FALSE
    ][[id_col]]
  } else {
    top_ids <- get_row_data(object)[
      order(rowVars(assays(object)[[assay]]), decreasing = TRUE)[1:n_top],
      , drop = FALSE
    ][[id_col]]
  }

  input_data <- assays(object[rowData(object)[[id_col]] %in% top_ids, ])[[assay]]

  if (scale) {
    print("Scaling data.")
    input_data <- scale(t(input_data))
  } else {
    input_data <- t(input_data)
  }

  if (calc_PCA) {
    print("Calculating PCA.")
    pca_res <- prcomp(input_data)
    var_explained <- cumsum((pca_res$sdev^2) / sum(pca_res$sdev^2))
    num_PCs <- max(which.max(var_explained > PC_var), min_PC)
    print(paste0(
      "Retaining ", num_PCs, " PCs, explaining ",
      round(var_explained[num_PCs] * 100, 2), "% of the total variance."
    ))
    input_data <- pca_res$x[, 1:num_PCs]
  }

  print("Calculating distance matrix.")

  if (metric %in% c("euclidean", "manhattan")) {
    distance_mat <- as.matrix(dist(input_data, method = metric))
  } else {
    distance_mat <- calc_cosine_distance(input_data)
  }

  print(paste0("Generating knn (k=", k, ") graph."))
  knn_graph <- make_knn_graph(distance_mat, k)

  print(paste0("Execute clustering using Louvain algorithm with resolution ", resolution, "."))
  knn_communities <- cluster_louvain(knn_graph, resolution = resolution)

  data.frame(
    sample = row.names(distance_mat),
    cluster = as.character(membership(knn_communities))
  )

}

