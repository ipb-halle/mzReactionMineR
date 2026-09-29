make_qc_object <- function() {
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      rt = matrix(
        c(
          0.5, 1.5,
          1.5, 2.5,
          2.5, 3.5
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:3), c("sample1", "sample2"))
      ),
      mz = matrix(
        c(
          95, 105,
          195, 205,
          295, 305
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(paste0("feature", 1:3), c("sample1", "sample2"))
      ),
      ims = matrix(
        c(
          5, 15,
          15, 25,
          25, 35
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
      ims = c(10, 20, 30),
      row.names = paste0("feature", 1:3)
    ),
    colData = data.frame(
      filename = c("sample1", "sample2")
    )
  )
}

test_that("qc_plots calculates correct mz_diff, rt_diff and ims_diff assays", {
  qc_object <- make_qc_object()
  output <- qc_plots(
    qc_object,
    rt_assay = "rt",
    mz_assay = "mz",
    ims_assay = "ims",
    rt_col = "rt",
    mz_col = "mz",
    ims_col = "ims",
    filename = "filename",
    what = c(
      "rt_diff",
      "mz_diff",
      "ims_diff"
    ),
    return = TRUE,
    print = FALSE
  )

  expect_equal(c(assays(output$se)$mz_diff), c(-5, -5, -5, 5, 5, 5))
  expect_equal(c(assays(output$se)$rt_diff), c(-.5, -.5, -.5, .5, .5, .5))
  expect_equal(c(assays(output$se)$ims_diff), c(-5, -5, -5, 5, 5, 5))
})

test_that("qc_plots returns list of ggplot objects", {

  qc_object <- make_qc_object()
  output <- qc_plots(
    qc_object,
    rt_assay = "rt",
    mz_assay = "mz",
    ims_assay = "ims",
    rt_col = "rt",
    mz_col = "mz",
    ims_col = "ims",
    what = c(
      "rt_diff",
      "mz_diff",
      "ims_diff",
      "rt_diff_sample",
      "mz_diff_sample",
      "ims_diff_sample"
    ),
    return = FALSE,
    print = FALSE
  )

  expect_s3_class(output[["rt_diff"]], "ggplot")
  expect_s3_class(output[["mz_diff"]], "ggplot")
  expect_s3_class(output[["ims_diff"]], "ggplot")
  expect_s3_class(output[["rt_diff_sample"]], "ggplot")
  expect_s3_class(output[["mz_diff_sample"]], "ggplot")
  expect_s3_class(output[["ims_diff_sample"]], "ggplot")

})

