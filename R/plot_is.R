#' Title
#'
#' @param object A SummarizedExperiment object
#' @param assay Character. The name of the assay to be normalized.
#' @param rt Numeric. The theoretical retention time of the internal standard.
#' @param mz Numeric. The theoretical m/z of the internal standard.
#' @param id Character. The id of the internal standard. If "none", the function
#'        will find the closest feature to the provided rt and mz values.
#'        Default is "none".
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerance for finding the internal standard.
#' @param rt_tolerance Numeric. Absolute retention time tolerance for finding
#'        the internal standard.
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns a ggplot object
#' @export
#'
plot_is <- function(
    object,
    assay,
    rt = NULL,
    mz = NULL,
    id = "none",
    mz_tolerance = c(0.005, 10),
    rt_tolerance = 0.1,
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {
  object_name <- object_argument_name(substitute(object))
  is_data <- get_is_data(
    object, assay, rt, mz, id, mz_tolerance, rt_tolerance,
    id_col, rt_col, mz_col, object_name = object_name
  )

  plot_data <- data.frame(
    x = 1:ncol(object),
    y = is_data$intensities
  )

  p <- ggplot(
    data = plot_data,
    aes(x = x, y = y, fill = y)
  ) +
    geom_point(shape = 21, size = 4, color = "black") +
    scale_y_continuous(limits = c(0,NA)) +
    scale_fill_viridis_c(option = "plasma") +
    theme_classic() +
    geom_smooth() +
    labs(x = "index", y = "Peak Area [a.u.]") +
    theme(
      legend.position = "none",
      legend.justification = "center",
      legend.background = element_rect(fill = NA, color = "black"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 7, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      strip.background = element_rect(color = NA, fill = NA),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.background = element_rect(fill = NA, color = NA),
      panel.background = element_rect(fill = NA, color = NA)
    )

  print(p)

  return(p)

}

#' @export
plotIS <- function(...) plot_is(...)

utils::globalVariables(c("x", "y"))
