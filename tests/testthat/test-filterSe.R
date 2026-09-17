make_filter_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      intensity = matrix(
        c(
          10, 10, 0,
          10, 0, 0,
          10, 10, 10,
          10, 0, 10
        ),
        nrow = 3,
        byrow = FALSE,
        dimnames = list(paste0("feature", 1:3), paste0("sample", 1:4))
      )
    ),
    rowData = data.frame(
      id = paste0("feature", 1:3),
      rt = c(1, 2, 3),
      mz = c(100, 101, 102),
      ion_mobility = c(0.8, 1.2, 1.6),
      row.names = paste0("feature", 1:3)
    ),
    colData = data.frame(
      filename = paste0("sample", 1:4),
      group = c("control", "control", "treated", "treated"),
      row.names = paste0("sample", 1:4)
    )
  )
}

test_that("filterSe filters by abundance without ion mobility", {
  object <- make_filter_object()

  result <- filterSe(
    object = object,
    assay = "intensity",
    min_pct = 0.75,
    min_n = 1L
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    "feature1"
  )
})

test_that("filterSe combines grouped and total abundance filters with or", {
  object <- make_filter_object()

  result <- filterSe(
    object = object,
    assay = "intensity",
    group_col = "group",
    min_pct_total = 1,
    min_n_total = 4L,
    min_pct_group = 1,
    min_n_group = 2L,
    min_count_rule = "or"
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    c("feature1", "feature3")
  )

  result <- filterSe(
    object = object,
    assay = "intensity",
    group_col = "group",
    min_pct_total = 1,
    min_n_total = 4L,
    min_pct_group = 0.5,
    min_n_group = 1L,
    min_count_rule = "or"
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    paste0("feature", 1:3)
  )
})

test_that("filterSe combines grouped and total abundance filters with and", {
  object <- make_filter_object()

  result <- filterSe(
    object = object,
    assay = "intensity",
    group_col = "group",
    min_pct_total = 0.5,
    min_n_total = 2L,
    min_pct_group = 1,
    min_n_group = 2L,
    min_count_rule = "and"
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    c("feature1", "feature3")
  )
})

test_that("filterSe optionally filters by ion mobility", {
  object <- make_filter_object()

  result <- filterSe(
    object = object,
    assay = "intensity",
    min_pct = 0,
    min_n = 0L,
    mobility_range = c(1, 2)
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    c("feature2", "feature3")
  )
})

test_that("filterSe_ims delegates to filterSe with ion mobility", {
  object <- make_filter_object()

  result <- filterSe_ims(
    object = object,
    assay = "intensity",
    min_pct = 0,
    mobility_range = c(1, 1.5)
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    "feature2"
  )
})

test_that("filterSe_ims forwards grouped and total abundance filters", {
  object <- make_filter_object()

  result <- filterSe_ims(
    object = object,
    assay = "intensity",
    group_col = "group",
    min_pct_total = 1,
    min_n_total = 4L,
    min_pct_group = 1,
    min_n_group = 2L,
    min_count_rule = "or"
  )

  expect_identical(
    SummarizedExperiment::rowData(result)$id,
    c("feature1", "feature3")
  )
})

test_that("filterSe validates object and assay inputs", {
  object <- make_filter_object()

  expect_error(
    filterSe(object = data.frame(), assay = "intensity"),
    "object.*SummarizedExperiment"
  )
  expect_error(
    filterSe(object = object, assay = "missing"),
    "assay.*name.*assay"
  )
})

test_that("filterSe validates metadata columns", {
  object <- make_filter_object()

  expect_error(
    filterSe(object = object, assay = "intensity", mz_col = "missing"),
    "missing from rowData"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", sample_col = "missing"),
    "missing from colData"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", group_col = "missing"),
    "missing from colData"
  )
  expect_error(
    filterSe(
      object = object,
      assay = "intensity",
      mobility_range = c(1, 2),
      ion_mobility_col = "missing"
    ),
    "missing from rowData"
  )
})

test_that("filterSe validates numeric filter arguments", {
  object <- make_filter_object()

  expect_error(
    filterSe(object = object, assay = "intensity", min_pct = 2),
    "min_pct.*between 0 and 1"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_n = 1.5),
    "min_n.*integer"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_pct_total = 2),
    "min_pct_total.*between 0 and 1"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_n_total = 1.5),
    "min_n_total.*integer"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_pct_group = 2),
    "min_pct_group.*between 0 and 1"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_n_group = 1.5),
    "min_n_group.*integer"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", min_count_rule = "missing"),
    "'arg' should be one of"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", mz_range = c(200, 100)),
    "mz_range"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", mobility_range = c(2, 1)),
    "mobility_range"
  )
})

test_that("filterSe validates sample and not_in values", {
  object <- make_filter_object()
  object_missing_sample <- object
  SummarizedExperiment::colData(object_missing_sample)$filename[1] <- "missing"

  expect_error(
    filterSe(object = object_missing_sample, assay = "intensity"),
    "assay column names are missing"
  )
  expect_error(
    filterSe(object = object, assay = "intensity", not_in = "control"),
    "not_in.*group_col"
  )
  expect_error(
    filterSe(
      object = object,
      assay = "intensity",
      group_col = "group",
      not_in = "missing"
    ),
    "not_in.*missing from colData"
  )
})