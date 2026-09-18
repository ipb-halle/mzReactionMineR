make_remove_feature_mz_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = matrix(1:6, nrow = 3)),
    rowData = S4Vectors::DataFrame(
      id = c("feature1", "feature2", "feature3"),
      mz = c(100, 200, 300)
    )
  )
}

test_that("remove_feature_mz removes features within the m/z tolerance", {
  result <- remove_feature_mz(
    make_remove_feature_mz_object(),
    mz = c(100.001, 300),
    mz_tolerance = c(0.01, 5)
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    c("feature2")
  )
})

test_that("remove_feature_mz validates its inputs", {
  object <- make_remove_feature_mz_object()

  expect_error(remove_feature_mz(list(), mz = 100), "SummarizedExperiment")
  expect_error(remove_feature_mz(object, mz = NA_real_), "finite numeric")
  expect_error(
    remove_feature_mz(object, mz = 100, mz_tolerance = c(0.01, -1)),
    "non-negative finite"
  )
  expect_error(
    remove_feature_mz(object, mz = 100, mz_col = "missing"),
    "columns are missing"
  )
})