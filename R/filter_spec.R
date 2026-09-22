#' filter_spec
#' A function that filters a spectrum object based on the number of peaks and
#'     intensity threshold.
#'
#' @importFrom Spectra filterIntensity applyProcessing mz intensity filterMzValues filterPrecursorPeaks
#' @param sps spectrum object. Input data
#' @param intensity_threshold numeric. minimum intensity threshold. Defaults to 0
#'     only meaningful for non-normalized data.
#' @param intensity_threshold_pct numeric. minimum intensity in % of max peak
#'     Peaks below will be removed.
#' @param min_peaks Integer. Minimum number of peaks after thresholding.
#' @param max_peaks Integer. Maximum number of peaks after thresholding.
#'     Choses the most intense peaks if more than max_peaks are present.
#' @param max_pct numeric. Maximum percentage of peaks to keep, if n_peaks > max_peaks.
#' @param remove_precursor logical. If TRUE, precursor peaks (+/- 10 mz) will be removed.
#' @param tolerance numeric. Tolerance for precursor peak removal in Da.
#' @param ppm numeric. ppm for precursor peak removal.
#' @returns A spectrum object with filtered peaks.
#' @export
filter_spec <- function(sps,
                        intensity_threshold = 0,
                        intensity_threshold_pct = 1,
                        min_peaks = 4,
                        max_peaks = 50,
                        max_pct = 0.95,
                        remove_precursor = TRUE,
                        tolerance = 0.005,
                        ppm = 20) {
  if (!methods::is(sps, "Spectra")) {
    stop("sps must be a Spectra object.", call. = FALSE)
  }
  numeric_args <- list(
    intensity_threshold = intensity_threshold,
    intensity_threshold_pct = intensity_threshold_pct,
    min_peaks = min_peaks,
    max_peaks = max_peaks,
    max_pct = max_pct,
    tolerance = tolerance,
    ppm = ppm
  )
  if (any(vapply(numeric_args, function(x) {
    length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x)
  }, logical(1)))) {
    stop("numeric arguments must each be one finite numeric value.", call. = FALSE)
  }
  if (min_peaks < 0 || max_peaks < min_peaks ||
      min_peaks != as.integer(min_peaks) || max_peaks != as.integer(max_peaks)) {
    stop("min_peaks and max_peaks must be non-negative integers with max_peaks >= min_peaks.", call. = FALSE)
  }
  if (intensity_threshold < 0 || intensity_threshold_pct < 0 ||
      intensity_threshold_pct > 100 || max_pct < 0 || max_pct > 1 ||
      tolerance < 0 || ppm < 0) {
    stop("thresholds, max_pct, tolerance, and ppm must be within their valid ranges.", call. = FALSE)
  }
  if (length(remove_precursor) != 1L || is.na(remove_precursor) ||
      !is.logical(remove_precursor)) {
    stop("remove_precursor must be a single logical value.", call. = FALSE)
  }

  filter_one <- function(spectrum) {
    if(remove_precursor) {
      spectrum <- filterPrecursorPeaks(spectrum,
                                       tolerance = tolerance,
                                       ppm = ppm,
                                       mz = ">=")
    }
    intensities <- unlist(intensity(spectrum), use.names = FALSE)
    if (!length(intensities)) {
      return(NULL)
    }
    min_intensity <- max(
      intensity_threshold_pct/100 * max(intensities),
      intensity_threshold
    )
    spectrum <- filterIntensity(
      spectrum,
      intensity = c(min_intensity,Inf)
    )
    n_peaks <- length(unlist(mz(spectrum), use.names = FALSE))
    if (n_peaks < min_peaks) {
      return(NULL)
    }
    if (n_peaks > max_peaks) {
      intensities <- unlist(intensity(spectrum), use.names = FALSE)
      idx <- order(intensities, decreasing = TRUE)[1:max_peaks]
      cumulative_pct <- cumsum(intensities[idx] / sum(intensities[idx]))
      keep <- which(cumulative_pct >= max_pct)[1L]
      idx_new <- idx[seq_len(ifelse(is.na(keep), max_peaks, keep))]
      spectrum <- filterMzValues(
        spectrum,
        mz = unlist(mz(spectrum), use.names = FALSE)[idx_new],
        ppm = 0
      )
    }
    applyProcessing(spectrum)
  }

  filtered <- Filter(Negate(is.null), lapply(seq_along(sps), function(i) {
    filter_one(sps[i])
  }))
  if (!length(filtered)) {
    return(sps[integer(0)])
  }
  do.call(c, filtered)
}
