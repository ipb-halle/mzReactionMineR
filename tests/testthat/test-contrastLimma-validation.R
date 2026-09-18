make_contrast_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = matrix(1:12, nrow = 3, ncol = 4)),
    rowData = data.frame(id = paste0("feature", 1:3)),
    colData = data.frame(
      group = factor(c("control", "control", "treated", "treated")),
      batch = factor(c("one", "two", "one", "two")),
      row.names = paste0("sample", 1:4)
    )
  )
}

test_that("contrastLimma validates its object and assay", {
  input_object <- data.frame()

  expect_error(
    contrastLimma(object = input_object, assay = "intensity", contrast_variable = "group"),
    "input_object.*SummarizedExperiment"
  )

  object <- make_contrast_object()
  expect_error(
    contrastLimma(object = object, assay = "missing", contrast_variable = "group"),
    "assay.*name.*assay.*object"
  )
})

test_that("contrastLimma validates design and contrast levels", {
  object <- make_contrast_object()

  expect_error(
    contrastLimma(object = object, assay = "intensity", contrast_variable = "missing"),
    "design columns are missing"
  )
  expect_error(
    contrastLimma(object = object, assay = "intensity", contrast_variable = "group", blocking_variables = "group"),
    "Contrast and blocking variables"
  )
  expect_error(
    contrastLimma(object = object, assay = "intensity", contrast_variable = "group", controls = "missing"),
    "levels from"
  )
})

test_that("contrastLimma validates its adjustment method", {
  expect_error(
    contrastLimma(
      object = make_contrast_object(),
      assay = "intensity",
      contrast_variable = "group",
      padj_method = "invalid"
    ),
    "padj_method.*valid"
  )
})