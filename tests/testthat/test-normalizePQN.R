make_normalize_pqn_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      raw = matrix(
        c(1, 2, 10, 20),
        nrow = 2,
        dimnames = list(c("feature1", "feature2"), c("sample1", "sample2"))
      )
    )
  )
}

test_that("normalizePQN returns a normalized assay", {
  object <- make_normalize_pqn_object()

  result <- normalizePQN(object, assay = "raw")

  expect_s4_class(result, "SummarizedExperiment")
  expect_true("pqn_normalized" %in% names(SummarizedExperiment::assays(result)))
  expect_identical(
    SummarizedExperiment::assay(result, "pqn_normalized"),
    matrix(
      c(5.5, 11, 5.5, 11),
      nrow = 2,
      dimnames = dimnames(SummarizedExperiment::assay(object, "raw"))
    )
  )
})

test_that("normalizePQN calculates references from named samples", {
  object <- make_normalize_pqn_object()
  expected <- matrix(
    c(1, 2, 1, 2),
    nrow = 2,
    dimnames = dimnames(SummarizedExperiment::assay(object, "raw"))
  )

  expect_identical(
    SummarizedExperiment::assay(
      normalizePQN(object, "raw", reference_samples = "sample1"),
      "pqn_normalized"
    ),
    expected
  )
})

test_that("normalizePQN supports custom assay names and mean references", {
  result <- normalizePQN(
    make_normalize_pqn_object(),
    assay = "raw",
    type = "mean",
    new_assay_name = "normalized"
  )

  expect_true("normalized" %in% names(SummarizedExperiment::assays(result)))
})

test_that("normalizePQN reports invalid object and assay inputs", {
  object <- make_normalize_pqn_object()
  input_object <- matrix(1:4, nrow = 2)

  expect_error(normalizePQN(input_object), "input_object.*SummarizedExperiment")
  expect_error(normalizePQN(object, assay = "missing"), "assay.*object")
})

test_that("normalizePQN reports invalid reference sample selections", {
  object <- make_normalize_pqn_object()
  unnamed_object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = unname(SummarizedExperiment::assay(object, "raw")))
  )

  expect_error(
    normalizePQN(object, "raw", reference_samples = 1),
    "sample names"
  )
  expect_error(
    normalizePQN(object, "raw", reference_samples = c(TRUE, FALSE)),
    "sample names"
  )
  expect_error(
    normalizePQN(object, "raw", reference_samples = character()),
    "unique"
  )
  expect_error(
    normalizePQN(object, "raw", reference_samples = c("sample1", "sample1")),
    "unique"
  )
  expect_error(
    normalizePQN(object, "raw", reference_samples = "missing"),
    "present in the assay"
  )
  expect_error(
    normalizePQN(unnamed_object, "raw", reference_samples = "sample1"),
    "column names"
  )
})
