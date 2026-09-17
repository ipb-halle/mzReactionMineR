make_anova_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = matrix(1:12, nrow = 3, ncol = 4)),
    rowData = data.frame(
      id = paste0("feature", 1:3),
      rt = 1:3,
      mz = 100:102,
      row.names = paste0("feature", 1:3)
    ),
    colData = data.frame(
      group = factor(c("control", "control", "treated", "treated")),
      row.names = paste0("sample", 1:4)
    )
  )
}

test_that("anovaLimma validates its object and assay", {
  expect_error(
    anovaLimma(object = data.frame(), assay = "intensity", test_variables = "group"),
    "object.*SummarizedExperiment"
  )

  object <- make_anova_object()
  expect_error(
    anovaLimma(object = object, assay = "missing", test_variables = "group"),
    "assay.*name.*assay"
  )
})

test_that("anovaLimma validates design columns", {
  object <- make_anova_object()

  expect_error(
    anovaLimma(object = object, assay = "intensity", test_variables = "missing"),
    "design columns are missing"
  )
  expect_error(
    anovaLimma(object = object, assay = "intensity", test_variables = "group", blocking_variables = "group"),
    "Test and blocking variables"
  )

  object_with_one_level <- make_anova_object()
  SummarizedExperiment::colData(object_with_one_level)$group <- factor("control")
  expect_error(
    anovaLimma(object = object_with_one_level, assay = "intensity", test_variables = "group"),
    "at least two unique values"
  )

  object_with_na <- make_anova_object()
  SummarizedExperiment::colData(object_with_na)$group[1] <- NA
  expect_error(
    anovaLimma(object = object_with_na, assay = "intensity", test_variables = "group"),
    "must not contain missing values"
  )
})

test_that("anovaLimma validates result and adjustment arguments", {
  object <- make_anova_object()

  expect_error(
    anovaLimma(object = object, assay = "intensity", test_variables = "group", padj_method = "invalid"),
    "padj_method.*valid"
  )
  expect_error(
    anovaLimma(object = object, assay = "intensity", test_variables = "group", return_colums = "missing"),
    "return columns are missing"
  )
})
