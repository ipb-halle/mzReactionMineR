make_knn_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      intensity = matrix(
        c(
          1, 2, 3, 4,
          4, 3, 2, 1,
          1, 1, 2, 2,
          2, 4, 1, 3
        ),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:4), paste0("sample", 1:4))
      )
    ),
    rowData = data.frame(id = paste0("feature", 1:4), row.names = paste0("feature", 1:4)),
    colData = data.frame(row.names = paste0("sample", 1:4))
  )
}

test_that("knn_clustering_samples returns one cluster per sample", {
  result <- knn_clustering_samples(
    make_knn_object(),
    assay = "intensity",
    n_top = 3,
    calc_PCA = FALSE,
    k = 2
  )

  expect_named(result, c("sample", "cluster"))
  expect_identical(result$sample, paste0("sample", 1:4))
  expect_equal(nrow(result), 4)
  expect_true(all(nzchar(result$cluster)))
})

test_that("knn_clustering_samples supports variance selection and euclidean distance", {
  result <- knn_clustering_samples(
    make_knn_object(),
    assay = "intensity",
    filter_type = "variance",
    metric = "euclidean",
    calc_PCA = TRUE,
    min_PC = 1,
    PC_var = 0.5,
    k = 1
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$sample == paste0("sample", 1:4)))
})

test_that("knn_clustering_samples validates inputs", {
  object <- make_knn_object()

  expect_error(
    knn_clustering_samples(object, assay = "missing", k = 1),
    "assay.*name.*assay"
  )
  expect_error(
    knn_clustering_samples(object, assay = "intensity", filter_type = "wrong", k = 1),
    "filter_type"
  )
  expect_error(
    knn_clustering_samples(object, assay = "intensity", k = 4),
    "k.*between 1"
  )
  expect_error(
    knn_clustering_samples(object, assay = "intensity", PC_var = 1.5, k = 1),
    "PC_var"
  )
})