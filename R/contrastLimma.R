#' contrastLimma
#'
#' A wrapper function that performs pairwise contrasts on a summarized experiment using limma.
#'
#' @importFrom dplyr %>% mutate bind_rows
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom stats contr.sum as.formula model.matrix contrasts<- contrasts
#' @param object A summarized experiment object.
#' @param assay Character. The name of the assay on which contrasts should be executed.
#' @param blocking_variables Character vector. The names of columns in colData(object) to use as blocking factors.
#' @param contrast_variable Character. The name of the column in colData(object) containing the contrast levels.
#' @param controls Character vector. The levels in contrast_variable to use as controls.
#' @param padj_method Character. The method used to adjust p-values. Default is "fdr".
#' @return A data.frame containing the feature metadata and pairwise limma results.
#'
#' @export

contrastLimma <- function(object = NULL,
                       assay = NULL,
                       blocking_variables = NULL,
                       contrast_variable = NULL,
                       controls = NULL,
                       padj_method = "fdr") {
  object_name <- object_argument_name(substitute(object))

  if(!inherits(object, "SummarizedExperiment")) {
    stop("'", object_name, "' must be a SummarizedExperiment.")
  }

  if(length(assay) != 1L || !is.character(assay) ||
     !assay %in% SummarizedExperiment::assayNames(object)) {
    stop("'assay' must name an assay in '", object_name, "'.")
  }

  if(length(contrast_variable) != 1L || !is.character(contrast_variable) ||
     is.na(contrast_variable) || !nzchar(contrast_variable)) {
    stop("'contrast_variable' must be one non-empty column name.")
  }

  sample_data <- as.data.frame(SummarizedExperiment::colData(object))

  if(!is.null(blocking_variables) &&
     (length(blocking_variables) == 0L || !is.character(blocking_variables) ||
      anyNA(blocking_variables) || any(!nzchar(blocking_variables)))) {
    stop("'blocking_variables' must contain non-empty column names or be NULL.")
  }

  if(!is.null(blocking_variables) && anyDuplicated(blocking_variables)) {
    stop("'blocking_variables' must not contain duplicate column names.")
  }

  if(!is.null(blocking_variables) && contrast_variable %in% blocking_variables) {
    stop("Contrast and blocking variables must be different columns.")
  }

  if(!is.null(controls) &&
     (!is.character(controls) || length(controls) == 0L || anyNA(controls) ||
      any(!nzchar(controls)) || anyDuplicated(controls))) {
    stop("'controls' must contain unique, non-empty level names or be NULL.")
  }

  if(length(padj_method) != 1L || !is.character(padj_method) ||
     !padj_method %in% stats::p.adjust.methods) {
    stop("'padj_method' must be a valid stats::p.adjust method.")
  }

  design_variables <- c(contrast_variable, blocking_variables)
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

  sample_data[[contrast_variable]] <- make.names(sample_data[[contrast_variable]])
  contrast_levels <- unique(sample_data[[contrast_variable]])
  if(!is.null(controls)) {
    controls <- make.names(controls)
    if(!all(controls %in% contrast_levels)) {
      stop("'controls' must contain levels from 'contrast_variable'.")
    }
  }

  formula_string <- paste("~ 0 +", contrast_variable)

  X <- model.matrix(
    data = sample_data,
    as.formula(formula_string)
  )

  colnames(X) <- gsub(contrast_variable, "", colnames(X))

  if(!is.null(blocking_variables)) {
    blocking_design <- lapply(blocking_variables, function(block_name) {
      block <- factor(sample_data[[block_name]])
      contrasts(block) <- contr.sum(nlevels(block))
      block_design <- model.matrix(~ block)[, -1, drop = FALSE]
      colnames(block_design) <- levels(block)[-1]
      block_design
    })
    X <- cbind(do.call(cbind, blocking_design), X)
  }

  if(!is.null(controls)) {
    treatments <- unique(
      sample_data[,contrast_variable][
        !sample_data[,contrast_variable] %in% controls
      ]
    )

    contrast_pairs <- expand.grid(
      controls,
      treatments
    )

  } else {

    treatments <- unique(
      sample_data[,contrast_variable]
    )

    contrast_pairs <- expand.grid(
      treatments,
      treatments
    )

  }

  contrast_names <- contrast_pairs %>%
    mutate(
      contrast = paste0(.data$Var2, "-", .data$Var1)
    ) %>%
    `[[`("contrast")

  contrast_matrix <- makeContrasts(
    contrasts = contrast_names,
    levels = colnames(X)
  )

  # Replace missing intensities because limma cannot fit NA values.
  fit <- lmFit(
    replace(SummarizedExperiment::assays(object)[[assay]],
            is.na(SummarizedExperiment::assays(object)[[assay]]), 0),
    design = X
  )

  fit_contrast <- contrasts.fit(fit, contrasts = contrast_matrix)

  fit2 <- eBayes(fit_contrast)

  pairwise_res <- list()

  for(i in seq_along(contrast_names)) {

    pairwise_res[[contrast_names[i]]] <- data.frame(
      as.data.frame(SummarizedExperiment::rowData(object)),
      topTable(
        fit2,
        number = Inf,
        sort.by = "none",
        adjust.method = padj_method,
        coef = i
      )
    )

  }

  pairwise_res <- pairwise_res %>%
    bind_rows(.id = "contrast")

  return(pairwise_res)

}

utils::globalVariables(c(".data", "."))

