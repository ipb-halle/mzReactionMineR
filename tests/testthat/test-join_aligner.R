make_join_aligner_se <- function(
  sample_names,
  feature_id,
  rt,
  mz,
  assays,
  mobility = NULL,
  col_data = data.frame(row.names = sample_names)
) {
  row_data <- data.frame(id = feature_id, rt = rt, mz = mz)
  if (!is.null(mobility)) {
    row_data$mobility <- mobility
  }

  SummarizedExperiment::SummarizedExperiment(
    rowData = row_data,
    assays = assays,
    colData = col_data
  )
}

test_that("join_aligner aligns a multi-sample SummarizedExperiment", {
  object <- make_join_aligner_se(
    sample_names = c("first", "second"),
    feature_id = c("a", "b", "c"),
    rt = c(1, 5, 8),
    mz = c(100, 200, 300),
    assays = list(area = matrix(
      c(10, 20, 30, 11, 21, 31),
      nrow = 3,
      dimnames = list(NULL, c("first", "second"))
    )),
    col_data = data.frame(
      batch = c("a", "b"),
      row.names = c("first", "second")
    )
  )

  result <- join_aligner(object, mz_tolerance = c(0.05, 20), rt_tolerance = 0.2)

  expect_s4_class(result, "SummarizedExperiment")
  expect_equal(nrow(result), 3)
  expect_identical(colnames(result), c("first", "second"))
  expect_identical(SummarizedExperiment::colData(result)$batch, c("a", "b"))
  expect_identical(names(SummarizedExperiment::assays(result)), "area")
  expect_equal(SummarizedExperiment::assay(result, "area"), SummarizedExperiment::assay(object, "area"))
})

test_that("join_aligner returns only assays shared by all input objects", {
  first <- make_join_aligner_se(
    sample_names = "first",
    feature_id = c("a", "c"),
    rt = c(1, 8),
    mz = c(100, 300),
    assays = list(
      area = matrix(c(10, 30), ncol = 1, dimnames = list(NULL, "first")),
      height = matrix(c(100, 300), ncol = 1, dimnames = list(NULL, "first"))
    ),
    col_data = data.frame(batch = "a", row.names = "first")
  )
  second <- make_join_aligner_se(
    sample_names = "second",
    feature_id = c("b", "d"),
    rt = c(1.1, 8.1),
    mz = c(100.01, 300.01),
    assays = list(
      area = matrix(c(11, 31), ncol = 1, dimnames = list(NULL, "second")),
      intensity = matrix(c(1000, 3000), ncol = 1, dimnames = list(NULL, "second"))
    ),
    col_data = data.frame(group = "treated", row.names = "second")
  )

  result <- join_aligner(
    list(first, second),
    mz_tolerance = c(0.05, 20),
    rt_tolerance = 0.2
  )

  expect_equal(nrow(result), 2)
  expect_identical(colnames(result), c("first", "second"))
  expect_identical(as.data.frame(SummarizedExperiment::colData(result))$batch, c("a", NA))
  expect_identical(as.data.frame(SummarizedExperiment::colData(result))$group, c(NA, "treated"))
  expect_identical(names(SummarizedExperiment::assays(result)), "area")
  expect_equal(
    SummarizedExperiment::assay(result, "area"),
    matrix(c(10, 30, 11, 31), nrow = 2, dimnames = list(NULL, c("first", "second")))
  )
})

test_that("join_aligner returns colData in output sample order", {
  first <- make_join_aligner_se(
    sample_names = "first",
    feature_id = "a",
    rt = 1,
    mz = 100,
    assays = list(area = matrix(10, ncol = 1, dimnames = list(NULL, "first"))),
    col_data = data.frame(batch = "a", row.names = "first")
  )
  second <- make_join_aligner_se(
    sample_names = "second",
    feature_id = c("b", "c"),
    rt = c(1.1, 8),
    mz = c(100.01, 300),
    assays = list(area = matrix(c(11, 31), ncol = 1, dimnames = list(NULL, "second"))),
    col_data = data.frame(batch = "b", row.names = "second")
  )

  result <- join_aligner(
    list(first, second),
    mz_tolerance = c(0.05, 20),
    rt_tolerance = 0.2
  )

  expect_identical(rownames(SummarizedExperiment::colData(result)), colnames(result))
  expect_identical(as.data.frame(SummarizedExperiment::colData(result))$batch, c("b", "a"))
})

test_that("join_aligner_ims returns a SummarizedExperiment", {
  object <- make_join_aligner_se(
    sample_names = c("first", "second"),
    feature_id = c("a", "b"),
    rt = c(1, 5),
    mz = c(100, 200),
    mobility = c(1, 2),
    assays = list(area = matrix(
      c(10, 20, 11, 21),
      nrow = 2,
      dimnames = list(NULL, c("first", "second"))
    ))
  )

  result <- join_aligner_ims(
    object,
    ion_mobility_col = "mobility",
    mz_tolerance = c(0.05, 20),
    rt_tolerance = 0.2
  )

  expect_s4_class(result, "SummarizedExperiment")
  expect_true("mobility" %in% names(SummarizedExperiment::rowData(result)))
  expect_identical(names(SummarizedExperiment::assays(result)), "area")
})

test_that("join_aligner validates SummarizedExperiment inputs", {
  one_sample <- make_join_aligner_se(
    sample_names = "sample",
    feature_id = "a",
    rt = 1,
    mz = 100,
    assays = list(area = matrix(10, ncol = 1, dimnames = list(NULL, "sample")))
  )
  missing_mz <- SummarizedExperiment::SummarizedExperiment(
    rowData = data.frame(id = "a", rt = 1),
    assays = list(area = matrix(10, ncol = 1, dimnames = list(NULL, "sample"))),
    colData = data.frame(row.names = "sample")
  )
  no_shared_assays <- list(
    make_join_aligner_se(
      sample_names = "first",
      feature_id = "a",
      rt = 1,
      mz = 100,
      assays = list(area = matrix(10, ncol = 1, dimnames = list(NULL, "first")))
    ),
    make_join_aligner_se(
      sample_names = "second",
      feature_id = "b",
      rt = 1.1,
      mz = 100.01,
      assays = list(height = matrix(11, ncol = 1, dimnames = list(NULL, "second")))
    )
  )

  expect_error(join_aligner(data.frame(id = 1, rt = 1, mz = 1)), "SummarizedExperiment")
  expect_error(join_aligner(one_sample), "At least two samples")
  expect_error(join_aligner(missing_mz), "rowData.*missing column")
  expect_error(join_aligner(no_shared_assays), "share at least one assay")
  expect_error(join_aligner(one_sample, mz_tolerance = 0), "mz_tolerance")
  expect_error(join_aligner(one_sample, rt_tolerance = 0), "tolerances")
})
