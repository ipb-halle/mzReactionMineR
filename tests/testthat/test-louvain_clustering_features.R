make_louvain_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      intensity = matrix(
        c(
          1, 2, 3, 4, 5,
          4, 3, 2, 1, 0,
          0, 1, 2, 2, 3,
          2, 1, 2, 1, 2,
          1, 2, 1, 2, 1
        ),
        nrow = 5,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:5), paste0("sample", 1:5))
      )
    ),
    rowData = data.frame(id = paste0("feature", 1:5), row.names = paste0("feature", 1:5)),
    colData = data.frame(row.names = paste0("sample", 1:5))
  )
}

test_that("louvain_clustering_features clusters selected features", {
  result <- louvain_clustering_features(
    make_louvain_object(),
    assay = "intensity",
    n_top = 3,
    R_trsh = 0,
    min_degree = 0,
    min_cluster_size = 1,
    return_removed = FALSE
  )

  expect_named(result, c("id", "cluster", "cluster_connectivity"))
  expect_equal(nrow(result), 3)
  expect_true(all(result$id %in% paste0("feature", 1:3)))
  expect_true(all(result$cluster >= 1))
})

test_that("louvain_clustering_features supports variance filtering", {
  result <- louvain_clustering_features(
    make_louvain_object(),
    assay = "intensity",
    filter_type = "variance",
    n_top = 4,
    R_trsh = 0.2,
    type = "both",
    min_degree = 0,
    min_cluster_size = 1,
    return_removed = FALSE
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$cluster > 0))
})

test_that("louvain_clustering_features validates inputs", {
  object <- make_louvain_object()

  expect_error(
    louvain_clustering_features(object, assay = "missing"),
    "assay.*name.*assay"
  )
  expect_error(
    louvain_clustering_features(object, assay = "intensity", filter_type = "wrong"),
    "filter_type"
  )
  expect_error(
    louvain_clustering_features(object, assay = "intensity", type = "wrong"),
    "type"
  )
  expect_error(
    louvain_clustering_features(object, assay = "intensity", R_trsh = 2),
    "R_trsh"
  )
})