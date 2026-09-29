feature_table_path <- testthat::test_path(
  "test-data",
  "test_mzmine_full_feature_table.csv"
)
metadata_path <- testthat::test_path(
  "test-data",
  "test_mzmine_metadata.csv"
)

test_that("mzmine_to_se imports the supplied MZmine fixtures", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  feature_table <- utils::read.csv(feature_table_path)

  result <- mzmine_to_se(feature_table_path, metadata)

  expect_s4_class(result, "SummarizedExperiment")
  expect_identical(
    colnames(result),
    make.names(metadata$filename)
  )
  expect_identical(
    rowData(result)$id,
    as.character(feature_table$id)
  )
  expect_true(all(c("area", "height", "mz", "rt") %in% names(rowData(result))))
  expect_identical(rowData(result)$area, feature_table$area)
  expect_identical(rowData(result)$height, feature_table$height)
  expect_identical(rowData(result)$mz, feature_table$mz)
  expect_identical(rowData(result)$rt, feature_table$rt)
  expect_identical(names(assays(result)), c("area", "height", "mz", "rt"))
  expect_equal(nrow(result), nrow(feature_table))
  expect_equal(ncol(result), nrow(metadata))
})

test_that("mzmine_to_se imports a selected assay", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)

  result <- mzmine_to_se(
    feature_table_path,
    sample_meta_data = metadata,
    assays = "height"
  )

  expect_identical(names(assays(result)), "height")
  expect_equal(
    dim(SummarizedExperiment::assay(result, "height")),
    c(nrow(utils::read.csv(feature_table_path)), nrow(metadata))
  )
})

test_that("mzmine_to_se validates metadata and assay inputs", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)

  expect_error(
    mzmine_to_se(feature_table_path),
    "sample_meta_data.*data.frame"
  )
  expect_error(
    mzmine_to_se(feature_table_path, metadata, filenames = "missing"),
    "filenames.*column"
  )
  expect_error(
    mzmine_to_se(feature_table_path, metadata, assays = "missing"),
    "assays.*assay names"
  )
  expect_error(
    mzmine_to_se(feature_table_path, metadata, assays = character()),
    "assays.*assay names"
  )

  duplicate_metadata <- rbind(metadata, metadata[1, ])
  expect_error(
    mzmine_to_se(feature_table_path, duplicate_metadata, assays = "height"),
    "duplicate sample names"
  )
})

test_that("mzmine_to_se warns about metadata-only samples", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  metadata <- rbind(
    metadata,
    data.frame(
      filename = "metadata-only.d",
      sample_type = "sample",
      treatment = "extra",
      stringsAsFactors = FALSE
    )
  )

  expect_warning(
    result <- mzmine_to_se(feature_table_path, metadata, assays = "height"),
    "metadata.only.d"
  )
  expect_identical(
    colnames(result),
    make.names(metadata$filename[metadata$filename != "metadata-only.d"])
  )
})

test_that("mzmine_to_se warns about feature-only samples", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  removed_sample <- metadata$filename[1]
  metadata <- metadata[-1, , drop = FALSE]

  expect_warning(
    result <- mzmine_to_se(feature_table_path, metadata, assays = "height"),
    make.names(removed_sample)
  )
  expect_identical(
    colnames(result),
    make.names(metadata$filename)
  )
})

test_that("mzmine_to_se reports tables without MZmine datafile columns", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  invalid_table <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(id = 1), invalid_table, row.names = FALSE)
  on.exit(unlink(invalid_table), add = TRUE)

  expect_error(
    mzmine_to_se(invalid_table, metadata),
    "datafile columns"
  )
})

test_that("mzmine_to_se reports unusable assays", {
  metadata <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
  invalid_table <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      id = 1,
      `datafile.sample.min` = 1,
      check.names = FALSE
    ),
    invalid_table,
    row.names = FALSE
  )
  on.exit(unlink(invalid_table), add = TRUE)

    expect_error(
      mzmine_to_se(invalid_table, metadata),
      "assay names"
  )
})