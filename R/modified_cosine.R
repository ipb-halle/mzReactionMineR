#' modified cosine.
#'
#' A wrapper function to calculate all pairwise modified cosine scores based on
#'     the Spectra package.
#'
#'@importFrom Spectra precursorMz peaksData acquisitionNum
#'@importFrom MsCoreUtils gnps join_gnps
#'
#' @param sps Spectra object. Must contain multiple MS/MS spectra
#' @param tolerance numeric. Absolute mz tolerance
#' @param ppm numeric. PPM tolerance for matching precursor m/z values
#'
#' @returns An upper triangular matrix of modified cosine scores
#' @export
#'
modified_cosine <- function(sps, tolerance = 0, ppm = 10) {
  # Validate the core inputs before accessing spectrum metadata.
  if (!methods::is(sps, "Spectra")) {
    stop("sps must be a Spectra object.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.numeric(tolerance) || is.na(tolerance) ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a single non-negative numeric value.", call. = FALSE)
  }
  if (length(ppm) != 1L || !is.numeric(ppm) || is.na(ppm) ||
      !is.finite(ppm) || ppm < 0) {
    stop("ppm must be a single non-negative numeric value.", call. = FALSE)
  }

  n <- length(sps)
  if (!n) {
    return(matrix(NA_real_, nrow = 0L, ncol = 0L))
  }

  # Fill only the upper triangle; diagonal and lower triangle remain NA.
  scores <- matrix(NA_real_, nrow = n, ncol = n)
  rownames(scores) <- colnames(scores) <- acquisitionNum(sps)

  for (i in seq_len(n)) {
    spectrum_i <- peaksData(sps[i])[[1L]]
    precursor_mz_i <- precursorMz(sps[i])

    for (j in seq_len(n)) {
      if (i >= j) {
        next
      }
      spectrum_j <- peaksData(sps[j])[[1L]]
      precursor_mz_j <- precursorMz(sps[j])

      if (!nrow(spectrum_i) || !nrow(spectrum_j)) {
        next
      }

      # Match peaks using precursor m/z and the requested tolerance window.
      map <- join_gnps(
        x = spectrum_i[, 1L],
        y = spectrum_j[, 1L],
        xPrecursorMz = precursor_mz_i,
        yPrecursorMz = precursor_mz_j,
        tolerance = tolerance,
        ppm = ppm
      )

      if (!length(map[[1L]]) || !length(map[[2L]])) {
        next
      }

      scores[i, j] <- gnps(spectrum_i[map[[1L]], ], spectrum_j[map[[2L]], ])
    }
  }

  scores
}
