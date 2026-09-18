test_that("camelCase APIs wrap snake_case implementations", {
  snake_case_names <- c(
    "anova_limma",
    "blank_subtraction_se",
    "contrast_limma",
    "filter_se",
    "filter_se_ims",
    "impute_min_frac",
    "knn_clustering_samples",
    "louvain_clustering_features",
    "normalize_is",
    "normalize_pqn",
    "plot_is",
    "remove_feature_mz"
    ,"qc_plots",
    "get_col_data",
    "get_row_data",
    "mzmine_alignment_qc"
  )

  expect_true(all(snake_case_names %in% getNamespaceExports("mzReactionMineR")))
  camel_case_names <- c(
    "anovaLimma",
    "blankSubtractionSE",
    "contrastLimma",
    "filterSe",
    "filterSe_ims",
    "imputeMinFrac",
    "knnClusteringSamples",
    "louvainClusteringFeatures",
    "normalizeIS",
    "normalizePQN",
    "plotIS",
    "removeFeatureMz"
    ,"QC_plots",
    "get_colData",
    "get_rowData",
    "mzmine_alignment_QC"
  )

  expect_true(all(camel_case_names %in% getNamespaceExports("mzReactionMineR")))
  expect_true(all(vapply(
    c(snake_case_names, camel_case_names),
    function(name) is.function(get(name, asNamespace("mzReactionMineR"))),
    logical(1)
  )))
})
