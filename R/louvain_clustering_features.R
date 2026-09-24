#' louvain_clustering_features
#'
#' A function that performs clustering (community detection) on features in a
#' SummarizedExperiment object based on correlation across samples.
#'
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom MatrixGenerics rowVars rowMaxs
#' @importFrom igraph cluster_louvain membership vertex_attr vertex_attr<- V degree delete_vertices subgraph
#' @importFrom stats prcomp cor
#' @importFrom dplyr right_join
#'
#' @param object a SummarizedExperiment object
#' @param assay Character. The name of the assay to be used for clustering.
#' @param R_trsh Numeric. Correlation threshold for an edge to be drawn.
#' @param resolution Numeric. Resolution parameter for the Louvain algorithm.
#'        Default is 1.
#' @param min_degree Integer. minimum number of edges for a feature to be kept.
#' @param min_cluster_size Integer. Minimum cluster size to be kept after clustering.
#' @param metric Either "pearson", "kendall", or "spearman". The correlation
#'        metric to use for calculating the correlation matrix. Default is "pearson".
#' @param type Either "pos", "neg", or "both". What kinds of correlation will be
#'        used for the graph: only positive, only negative or both.
#' @param filter_type Character. The method to select features for clustering.
#'        Must be either "intsnity" or "variance". Default is "area".
#' @param n_top Integer. The number of features to select based on
#'        filter_type.
#' @param return_removed Logical. Whether to return features removed during the clustering process.
#'        removed features will be added as cluster -1.
#' @param calc_PCA Logical. Wether to perform PCA before clustering.
#'        Default is TRUE.
#' @param min_PC Integer. The minimum number of PCs to retain. Default is 5.
#' @param PC_var Numeric. The minimum cumulative variance that the retained PCs should explain.
#' @param scale Logical. Whether to scale the data before calculating the
#'        distance matrix. Default is TRUE.
#' @param id_col Character. The name of the column in rowData(object) that
#'        contains the feature ids. Default is "id".
#'
#' @returns a data.frame with feature ids, cluster labels, and within-cluster connectivity.
#' @export
#'
louvain_clustering_features <- function(
    object,
    assay,
    R_trsh = 0.8,
    resolution = 1,
    min_degree = 20,
    min_cluster_size = 50,
    metric = "pearson",
    type = "pos",
    filter_type = "intensity",
    n_top = nrow(object),
    return_removed = TRUE,
    calc_PCA = FALSE,
    min_PC = 5,
    PC_var = 0.8,
    scale = TRUE,
    id_col = "id"
) {

  # Validate the object, selected columns, and scalar controls.
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
  if (anyNA(rowData(object)[[id_col]]) || anyDuplicated(rowData(object)[[id_col]])) {
    stop("'id_col' must contain unique, non-missing feature identifiers.", call. = FALSE)
  }
  if (length(filter_type) != 1L || !filter_type %in% c("intensity", "variance")) {
    stop("'filter_type' must be either 'intensity' or 'variance'.", call. = FALSE)
  }
  if (length(n_top) != 1L || !is.numeric(n_top) || n_top < 1 ||
      n_top != as.integer(n_top) || n_top > nrow(object)) {
    stop("'n_top' must be an integer between 1 and the number of features.", call. = FALSE)
  }
  if (!is.logical(return_removed) || length(return_removed) != 1L || is.na(return_removed) ||
      !is.logical(calc_PCA) || length(calc_PCA) != 1L || is.na(calc_PCA) ||
      !is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("'return_removed', 'calc_PCA', and 'scale' must each be TRUE or FALSE.", call. = FALSE)
  }
  if (length(R_trsh) != 1L || !is.numeric(R_trsh) || !is.finite(R_trsh) ||
      R_trsh < -1 || R_trsh > 1) {
    stop("'R_trsh' must be one finite number between -1 and 1.", call. = FALSE)
  }
  if (length(resolution) != 1L || !is.numeric(resolution) ||
      !is.finite(resolution) || resolution <= 0) {
    stop("'resolution' must be one positive finite number.", call. = FALSE)
  }
  if (length(min_degree) != 1L || !is.numeric(min_degree) || min_degree < 0 ||
      min_degree != as.integer(min_degree)) {
    stop("'min_degree' must be a non-negative integer.", call. = FALSE)
  }
  if (length(min_cluster_size) != 1L || !is.numeric(min_cluster_size) ||
      min_cluster_size < 1 || min_cluster_size != as.integer(min_cluster_size)) {
    stop("'min_cluster_size' must be a positive integer.", call. = FALSE)
  }
  if (length(metric) != 1L || !metric %in% c("pearson", "kendall", "spearman")) {
    stop("'metric' must be either 'pearson', 'kendall', or 'spearman'.", call. = FALSE)
  }
  if (length(type) != 1L || !type %in% c("pos", "neg", "both")) {
    stop("'type' must be either 'pos', 'neg', or 'both'.", call. = FALSE)
  }
  if (length(min_PC) != 1L || !is.numeric(min_PC) || min_PC < 1 ||
      min_PC != as.integer(min_PC)) {
    stop("'min_PC' must be a positive integer.", call. = FALSE)
  }
  if (length(PC_var) != 1L || !is.numeric(PC_var) || !is.finite(PC_var) ||
      PC_var <= 0 || PC_var > 1) {
    stop("'PC_var' must be greater than 0 and at most 1.", call. = FALSE)
  }

  # Select the highest-intensity or highest-variance features.
  feature_data <- assays(object)[[assay]]
  feature_score <- switch(
    filter_type,
    intensity = rowMaxs(feature_data),
    variance = rowVars(feature_data)
  )
  selected_rows <- order(feature_score, decreasing = TRUE)[seq_len(n_top)]
  feature_ids <- rowData(object)[[id_col]][selected_rows]

  # Keep the same subset-based preparation pattern used by knn_clustering_samples
  # so row ordering and filtering semantics stay comparable to the last working implementation.
  input_data <- assays(object[rowData(object)[[id_col]] %in% feature_ids, ])[[assay]]
  if (any(!is.finite(input_data))) {
    stop("Selected assay values must be finite.", call. = FALSE)
  }

  # Arrange samples as rows, optionally scale, and optionally reduce dimensions.
  input_data <- t(input_data)
  if (scale) input_data <- scale(input_data)
  if (any(!is.finite(input_data))) {
    stop("Selected assay values must remain finite after scaling.", call. = FALSE)
  }
  if (calc_PCA) {
    pca_res <- prcomp(input_data)
    var_explained <- cumsum(pca_res$sdev^2) / sum(pca_res$sdev^2)
    target_pc <- which(var_explained >= PC_var)[1L]
    num_pcs <- max(min_PC, target_pc)
    if (num_pcs > ncol(pca_res$x)) {
      stop("'min_PC' cannot exceed the number of available principal components.", call. = FALSE)
    }
    input_data <- pca_res$x[, seq_len(num_pcs), drop = FALSE]
  }

  # Calculate the feature-to-feature correlation matrix and adjacency matrix.
  cor_mat <- cor(input_data, method = metric)
  adj_mat <- switch(
    type,
    pos = (cor_mat > R_trsh) * 1,
    neg = (cor_mat < R_trsh) * 1,
    both = (abs(cor_mat) > R_trsh) * 1
  )
  diag(adj_mat) <- 0
  if (any(!is.finite(adj_mat))) {
    stop("The correlation matrix contains non-finite values.", call. = FALSE)
  }

  # Build the graph and remove features with insufficient connectivity.
  graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  vertex_attr(graph)$name <- feature_ids
  features_to_remove <- V(graph)[degree(graph) < min_degree]$name
  if (length(features_to_remove) == length(V(graph))) {
    stop(paste0("All features have degree < ", min_degree, ". Try lowering the threshold."), call. = FALSE)
  }
  graph_clean <- delete_vertices(graph, features_to_remove)
  communities <- cluster_louvain(graph_clean, resolution = resolution)
  clusters <- data.frame(
    feature_id = V(graph_clean)$name,
    cluster = as.character(membership(communities)),
    stringsAsFactors = FALSE
  )
  names(clusters)[1L] <- id_col

  # Remove small communities and renumber the remaining communities by size.
  clusters_to_remove <- clusters %>%
    group_by(cluster) %>%
    summarize(n = n(), .groups = "drop") %>%
    filter(n < min_cluster_size) %>%
    pull(cluster)
  if (length(clusters_to_remove) == length(unique(clusters$cluster))) {
    stop(paste0("No clusters survive threshold of ", min_cluster_size,
                " features. Try lowering the threshold."), call. = FALSE)
  }
  clusters_clean <- clusters %>%
    filter(!cluster %in% clusters_to_remove) %>%
    left_join(
      clusters %>%
        group_by(cluster) %>%
        summarize(cluster_size = n(), .groups = "drop") %>%
        arrange(desc(cluster_size), cluster) %>%
        mutate(cluster_new = row_number()),
      by = "cluster"
    ) %>%
    mutate(cluster = cluster_new) %>%
    select(all_of(id_col), cluster)

  # Calculate within-community degree and optionally append removed features.
  clusters_clean$cluster_connectivity <- vapply(seq_len(nrow(clusters_clean)), function(i) {
    vertices <- clusters_clean[[id_col]][clusters_clean$cluster == clusters_clean$cluster[i]]
    degree(subgraph(graph_clean, vertices))[clusters_clean[[id_col]][i]]
  }, numeric(1))
  if (return_removed) {
    removed_features <- data.frame(
      feature_id = c(features_to_remove, clusters[clusters$cluster %in% clusters_to_remove, id_col]),
      cluster = -1,
      cluster_connectivity = NA_real_,
      stringsAsFactors = FALSE
    )
    names(removed_features)[1L] <- id_col
    clusters_clean <- rbind(clusters_clean, removed_features)
  }
  row.names(clusters_clean) <- NULL
  clusters_clean

}



utils::globalVariables(c("cluster", "desc", "cluster_new", "cluster_size"))
