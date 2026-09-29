test_that("feature_table_to_se converts a generic feature table", {
  feature_table_path <- tempfile(fileext = ".csv")
  on.exit(unlink(feature_table_path), add = TRUE)
  utils::write.csv(
    data.frame(
      feature_id = c("feature-1", "feature-2"),
      mz = c(100, 200),
      sample_a = c(1, 2),
      sample_b = c(3, 4),
      check.names = FALSE
    ),
    feature_table_path,
    row.names = FALSE
  )
  metadata <- data.frame(
    filename = c("sample_b.d", "sample_a.d"),
    sample_type = c("sample", "blank")
  )

  result <- feature_table_to_se(
    feature_table_path,
    metadata,
    assay_name = "area"
  )

  expect_s4_class(result, "SummarizedExperiment")
  expect_identical(colnames(result), c("sample_b", "sample_a"))
  expect_identical(names(assays(result)), "area")
  expect_identical(rowData(result)$feature_id, c("feature-1", "feature-2"))
  expect_equal(rowData(result)$mz, c(100, 200))
  expect_identical(colData(result)$sample_type, c("sample", "blank"))
  expect_equal(
    unname(SummarizedExperiment::assay(result, "area")),
    matrix(c(3, 4, 1, 2), nrow = 2)
  )
})

test_that("metaboscape_to_se delegates to feature_table_to_se", {
  feature_table_path <- tempfile(fileext = ".csv")
  on.exit(unlink(feature_table_path), add = TRUE)
  utils::write.csv(
    data.frame(id = 1:2, sample_a = c(5, 6)),
    feature_table_path,
    row.names = FALSE
  )
  metadata <- data.frame(filename = "sample_a.d", group = "sample")

  expect_equal(
    metaboscape_to_se(feature_table_path, metadata),
    feature_table_to_se(feature_table_path, metadata)
  )
})

test_that("feature_table_to_se validates metadata and missing samples", {
  feature_table_path <- tempfile(fileext = ".csv")
  on.exit(unlink(feature_table_path), add = TRUE)
  utils::write.csv(
    data.frame(id = 1, sample_a = 2, row.names = FALSE),
    feature_table_path
  )

  expect_error(
    feature_table_to_se(feature_table_path, NULL),
    "sample_meta_data.*data.frame"
  )
  expect_error(
    feature_table_to_se(
      feature_table_path,
      data.frame(sample = "sample_a"),
      filename_column = "filename"
    ),
    "filename_column.*column"
  )
  expect_error(
    feature_table_to_se(
      feature_table_path,
      data.frame(filename = c("sample_a", "sample_a"))
    ),
    "duplicate sample names"
  )

  expect_warning(
    result <- feature_table_to_se(
      feature_table_path,
      data.frame(filename = c("sample_a", "missing.d"))
    ),
    "sample_meta_data.*not in the feature table"
  )
  expect_identical(colnames(result), "sample_a")
})