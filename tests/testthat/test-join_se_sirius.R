make_sirius_test_object <- function(ids = as.character(1:3)) {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = matrix(seq_along(ids), nrow = length(ids))),
    rowData = data.frame(id = ids, row.names = paste0("feature", ids))
  )
}

test_that("join_se_sirius joins every populated SIRIUS fixture", {
  sirius_dir <- testthat::test_path("test-data", "sirius")
  sirius_files <- list.files(sirius_dir, pattern = "\\.tsv$", full.names = TRUE)
  populated_files <- sirius_files[vapply(
    sirius_files,
    function(path) nrow(utils::read.delim(path)) > 0,
    logical(1)
  )]

  for (path in populated_files) {
    result <- suppressWarnings(
      join_se_sirius(make_sirius_test_object(), path)
    )
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), 3L, info = basename(path))
    expect_true("molecularFormula" %in% names(SummarizedExperiment::rowData(result)))
  }
})

test_that("join_se_sirius accepts empty SIRIUS output", {
  path <- testthat::test_path(
    "test-data", "sirius", "denovo_structure_identifications.tsv"
  )

  result <- join_se_sirius(make_sirius_test_object(), path)

  expect_equal(nrow(result), 3L)
})

test_that("join_se_sirius validates inputs and duplicate IDs", {
  object <- make_sirius_test_object()
  invalid_file <- tempfile(fileext = ".tsv")
  duplicate_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(invalid_file, duplicate_file)))
  utils::write.table(data.frame(other = 1), invalid_file, sep = "\t", row.names = FALSE)
  utils::write.table(
    data.frame(mappingFeatureId = c(1, 1), value = c("first", "second")),
    duplicate_file,
    sep = "\t",
    row.names = FALSE
  )

  expect_error(join_se_sirius(data.frame(), duplicate_file), "SummarizedExperiment")
  expect_error(join_se_sirius(object, "missing.tsv"), "existing file path")
  expect_error(join_se_sirius(object, invalid_file), "mappingFeatureId")
  expect_error(join_se_sirius(object, duplicate_file, id_col = "missing"), "rowData")
  expect_warning(
    result <- join_se_sirius(object, duplicate_file),
    "Duplicate mappingFeatureId"
  )
  expect_identical(SummarizedExperiment::rowData(result)$value[1], "first")
})