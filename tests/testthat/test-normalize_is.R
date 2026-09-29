make_normalize_is_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      raw = matrix(
        c(
          10, 20,
          2, 4,
          30, 60
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:3), c("sample1", "sample2"))
      )
    ),
    rowData = data.frame(
      id = paste0("feature", 1:3),
      rt = c(1, 2, 3),
      mz = c(100, 200, 300),
      row.names = paste0("feature", 1:3)
    )
  )
}

test_that("normalize_is normalizes by an internal-standard id", {
  result <- normalize_is(
    object = make_normalize_is_object(),
    assay = "raw",
    id = "feature2",
    remove = FALSE
  )

  expect_equal(
    SummarizedExperiment::assay(result, "is_normalized"),
    matrix(
      c(
        15, 15,
        3, 3,
        45, 45
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(paste0("feature", 1:3), c("sample1", "sample2"))
    )
  )
})

test_that("normalize_is resolves an internal standard by coordinates", {
  result <- normalize_is(
    object = make_normalize_is_object(),
    assay = "raw",
    rt = 2,
    mz = 200,
    remove = TRUE
  )

  expect_identical(SummarizedExperiment::rowData(result)$id,
                   c("feature1", "feature3"))
})

test_that("internal-standard functions report invalid inputs", {
  object <- make_normalize_is_object()

  expect_error(
    normalize_is(object, assay = "missing", id = "feature2"),
    "assay.*name.*assay"
  )
  expect_error(
    normalize_is(object, assay = "raw", id = "missing"),
    "id.*not found"
  )
  expect_error(
    normalize_is(object, assay = "raw", rt = 20, mz = 200),
    "No feature matching"
  )
  expect_error(
    normalize_is(object, assay = "raw", id = "feature2", mz_tolerance = -1),
    "mz_tolerance.*non-negative"
  )

  invalid_intensities <- object
  SummarizedExperiment::assay(invalid_intensities, "raw")[2, 1] <- 0
  expect_error(
    normalize_is(invalid_intensities, assay = "raw", id = "feature2"),
    "positive and finite"
  )
})
