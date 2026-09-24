make_blank_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      intensity = matrix(
        c(
          1, 1, 10, 10,
          5, NA, 10, NA,
          1, 1, 2, 2
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:3),
                        c("blank1", "blank2", "sample1", "sample2"))
      )
    ),
    rowData = data.frame(
      id = paste0("feature", 1:3),
      row.names = paste0("feature", 1:3)
    ),
    colData = data.frame(
      filename = c("blank1", "blank2", "sample1", "sample2"),
      row.names = c("blank1", "blank2", "sample1", "sample2")
    )
  )
}

test_that("blank_subtraction_se filters features and removes blanks", {
  result <- blank_subtraction_se(
    object = make_blank_object(),
    assay = "intensity",
    blanks = c("blank1", "blank2")
  )

  expect_identical(SummarizedExperiment::rowData(result)$id, "feature1")
  expect_identical(SummarizedExperiment::colData(result)$filename,
                   c("sample1", "sample2"))
})

test_that("blank_subtraction_se validates inputs", {
  object <- make_blank_object()
  input_object <- data.frame()

  expect_error(
    blank_subtraction_se(object = input_object, assay = "intensity",
                       blanks = "blank1"),
    "input_object.*SummarizedExperiment"
  )
  expect_error(
    blank_subtraction_se(object = object, assay = "missing", blanks = "blank1"),
    "assay.*name.*assay"
  )
  expect_error(
    blank_subtraction_se(object = object, assay = "intensity", blanks = "unknown"),
    "blanks.*present"
  )
  expect_error(
    blank_subtraction_se(object = object, assay = "intensity", blanks = "blank1",
                       ratio = -1),
    "ratio.*non-negative"
  )
  expect_error(
    blank_subtraction_se(object = object, assay = "intensity", blanks = "blank1",
                       ratio_type = "sum"),
    "arg.*should be one of"
  )
  duplicate_samples <- object
  SummarizedExperiment::colData(duplicate_samples)$filename[2] <- "blank1"
  expect_error(
    blank_subtraction_se(object = duplicate_samples, assay = "intensity",
                       blanks = "blank1"),
    "sample_col.*unique"
  )
  expect_error(
    blank_subtraction_se(object = object, assay = "intensity",
                       blanks = colData(object)$filename),
    "at least one non-blank"
  )
})
