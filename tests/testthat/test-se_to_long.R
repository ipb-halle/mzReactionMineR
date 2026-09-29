make_se_to_long_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      raw = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(NULL, c("sample1", "sample2"))
      )
    ),
    rowData = S4Vectors::DataFrame(id = c("feature1", "feature2")),
    colData = S4Vectors::DataFrame(
      filename = c("sample1", "sample2"),
      group = c("control", "case")
    )
  )
}

test_that("se_to_long returns assay values with feature and sample metadata", {
  result <- se_to_long(make_se_to_long_object(), assay = "raw")

  expect_s3_class(result, "data.frame")
  expect_named(result, c("id", "filename", "value", "group"))
  expect_equal(nrow(result), 4L)
  expect_equal(result$value, c(1, 3, 2, 4))
})

test_that("se_to_long reports a missing filename column", {
  object <- make_se_to_long_object()
  names(SummarizedExperiment::colData(object))[1] <- "sample"

  expect_error(
    se_to_long(object, assay = "raw"),
    "'filename' must name a column in colData\\(object\\)"
  )
})

test_that("se_to_long validates object, assay, and filename inputs", {
  object <- make_se_to_long_object()
  input_object <- matrix(1:4, nrow = 2)

  expect_error(se_to_long(input_object, assay = "raw"), "input_object.*SummarizedExperiment")
  expect_error(se_to_long(object, assay = "missing"), "assay.*object")
  expect_error(se_to_long(object, assay = "raw", filename = NA_character_), "filename.*non-empty")
})

test_that("se_to_long requires assay and metadata sample names to match", {
  object <- make_se_to_long_object()
  SummarizedExperiment::colData(object)$filename[2] <- "other"

  expect_error(
    se_to_long(object, assay = "raw"),
    "assay column names must match"
  )
})