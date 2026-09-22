#' anovaLimma
#'
#' A wrapper function that performs an ANOVA test on a summarized experiment using limma.
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats contr.sum as.formula model.matrix contrasts<-
#' @param object A summarized experiment object.
#' @param assay Character. What is the name of the assay the ANOVA should be
#'     executed on.
#' @param blocking_variables Character vector. The names of the columns in
#'     colData(object) that should be used as blocking factors
#' @param test_variables Character vector. The names of the columns in
#'     colData(object) that should be used as test variables.
#' @param padj_method Character. The method used to adjust the p-values.
#'     Default is "fdr".
#' @param return_colums Character vector. The columns from rowData(object) that are reported in results. Default is c("id", "rt", "mz").
#' @return A data.frame containing the requested feature metadata and the
#'     limma F-statistic, raw p-value, and adjusted p-value.
#' @export

anova_limma <- function(object = NULL,
                       assay = NULL,
                       blocking_variables = NULL,
                       test_variables = NULL,
                       padj_method = "fdr",
                       return_colums = c("id", "rt", "mz")) {
  object_name <- object_argument_name(substitute(object))

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.")
  }

  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.")
  }

  if(length(test_variables) == 0L || !is.character(test_variables)) {
    stop("'test_variables' must contain at least one column name.")
  }

  if(anyNA(test_variables) || any(!nzchar(test_variables))) {
    stop("'test_variables' must contain non-empty column names.")
  }

  if(anyDuplicated(test_variables)) {
    stop("'test_variables' must not contain duplicate column names.")
  }

  if(!is.null(blocking_variables) &&
     (length(blocking_variables) == 0L || !is.character(blocking_variables) ||
      anyNA(blocking_variables) || any(!nzchar(blocking_variables)))) {
    stop("'blocking_variables' must contain non-empty column names or be NULL.")
  }

  if(!is.null(blocking_variables) && anyDuplicated(blocking_variables)) {
    stop("'blocking_variables' must not contain duplicate column names.")
  }

  if(!is.null(blocking_variables) &&
     any(test_variables %in% blocking_variables)) {
    stop("Test and blocking variables must be different columns.")
  }

  if(!is.character(return_colums) || anyNA(return_colums) ||
     any(!nzchar(return_colums))) {
    stop("'return_colums' must contain non-empty column names.")
  }

  sample_data <- as.data.frame(SummarizedExperiment::colData(object))
  design_variables <- c(test_variables, blocking_variables)
  missing_variables <- setdiff(design_variables, colnames(sample_data))

  if(length(missing_variables) > 0L) {
    stop(
      "The following design columns are missing from colData(", object_name, "): ",
      paste(missing_variables, collapse = ", ")
    )
  }

  if(any(vapply(sample_data[design_variables], anyNA, logical(1)))) {
    stop("Design columns must not contain missing values.")
  }

  if(any(vapply(sample_data[design_variables], function(x) {
    length(unique(x)) < 2L
  }, logical(1)))) {
    stop("Each design column must contain at least two unique values.")
  }

  if(length(padj_method) != 1L || !is.character(padj_method) ||
     !padj_method %in% stats::p.adjust.methods) {
    stop("'padj_method' must be a valid stats::p.adjust method.")
  }

  row_data <- as.data.frame(SummarizedExperiment::rowData(object))
  missing_return_columns <- setdiff(return_colums, colnames(row_data))
  if(length(missing_return_columns) > 0L) {
    stop(
      "The following return columns are missing from rowData(", object_name, "): ",
      paste(missing_return_columns, collapse = ", ")
    )
  }

  # Use one coefficient per test level so limma can test the full factor.
  formula_string <- paste("~ 0 +", paste(test_variables, collapse = " + "))

  X <- model.matrix(
    data = sample_data,
    as.formula(formula_string)
  )

  # Keep coefficient names aligned with the original factor levels.
  for(name in test_variables){
    colnames(X) <- gsub(name, "", colnames(X))
  }

  n_blocking_coefficients <- 0L
  if(!is.null(blocking_variables)) {
    # Add sum-coded blocking factors before the variables being tested.
    blocking_design <- lapply(blocking_variables, function(block_name) {
      block <- factor(sample_data[[block_name]])
      contrasts(block) <- contr.sum(nlevels(block))
      block_design <- model.matrix(~ block)[, -1, drop = FALSE]
      colnames(block_design) <- levels(block)[-1]
      block_design
    })

    n_blocking_coefficients <- sum(vapply(blocking_design, ncol, integer(1)))
    X <- cbind(do.call(cbind, blocking_design), X)
  }

  # limma cannot fit missing assay values; match the package's existing
  # behavior by replacing them with zero before fitting.
  fit <- lmFit(
    replace(SummarizedExperiment::assays(object)[[assay]], is.na(SummarizedExperiment::assays(object)[[assay]]), 0),
    design = X
  )

  # Empirical Bayes moderation stabilizes the feature-wise variance estimates.
  fit2 <- eBayes(fit)

  anova_res <- cbind(
    as.data.frame(SummarizedExperiment::rowData(object)),
    # Test only the coefficients for test variables, excluding blocking terms.
    topTable(
      fit2,
      number=Inf,
      sort.by="none",
      adjust.method = padj_method,
      coef = seq.int(n_blocking_coefficients + 1L, ncol(X))
    )
  )

  return(anova_res[, c(return_colums, "F", "P.Value", "adj.P.Val"), drop = FALSE])
}

#' @export
anovaLimma <- function(...) anova_limma(...)