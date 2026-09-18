make_impute_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      raw = matrix(
        c(
          1, NA, 3,
          4, NA, 6
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          c("feature1", "feature2"),
          c("sample1", "sample2", "sample3")
        )
      )
    )
  )
}

test_that("imputeMinFrac imputes missing values using row minima", {
  object <- make_impute_object()

  result <- imputeMinFrac(
    object = object,
    assay = "raw",
    fraction = 5,
    new_assay_name = "imputed"
  )

  expect_identical(
    names(SummarizedExperiment::assays(result)),
    c("raw", "imputed")
  )
  expect_equal(
    SummarizedExperiment::assay(result, "imputed"),
    matrix(
      c(
        1, 0.2, 3,
        4, 0.8, 6
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("feature1", "feature2"),
        c("sample1", "sample2", "sample3")
      )
    )
  )
})

test_that("imputeMinFrac errors when a row is entirely missing", {
  object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      raw = matrix(
        c(NA_real_, NA_real_),
        nrow = 1,
        dimnames = list("feature1", c("sample1", "sample2"))
      )
    )
  )

  expect_error(
    imputeMinFrac(object = object, assay = "raw"),
    "only NA values"
  )
})

test_that("imputeMinFrac reports the supplied object name", {
  input_se <- make_impute_object()

  expect_error(
    imputeMinFrac(input_se, assay = "missing"),
    "missing.*`input_se`"
  )
  expect_warning(
    imputeMinFrac(input_se, assay = "raw", new_assay_name = "raw"),
    "already exists in `input_se`"
  )
})
