test_that("filter_spec validates its arguments", {
  spectra <- Spectra::Spectra(backend = Spectra::MsBackendMemory(), data = base::data.frame(
    spectrum_id = 1L,
    mz = base::I(list(c(100, 101, 102, 103))),
    intensity = base::I(list(c(1, 2, 3, 4)))
  ))

  expect_error(filter_spec(data.frame()), "Spectra object")
  expect_error(filter_spec(spectra, min_peaks = 2.5), "non-negative integers")
  expect_error(filter_spec(spectra, max_pct = 2), "valid ranges")
  expect_error(filter_spec(spectra, remove_precursor = NA), "single logical")
})

test_that("filter_spec removes spectra below the minimum peak count", {
  spectra <- Spectra::Spectra(backend = Spectra::MsBackendMemory(), data = base::data.frame(
    spectrum_id = 1:2,
    mz = base::I(list(c(100, 101), c(200, 201, 202, 203))),
    intensity = base::I(list(c(1, 2), c(1, 2, 3, 4)))
  ))

  result <- filter_spec(spectra, min_peaks = 4, remove_precursor = FALSE)

  expect_s4_class(result, "Spectra")
  expect_length(result, 1)
  expect_equal(length(Spectra::mz(result)[[1]]), 4)
})

test_that("filter_spec returns an empty Spectra object when nothing passes", {
  spectra <- Spectra::Spectra(backend = Spectra::MsBackendMemory(), data = base::data.frame(
    spectrum_id = 1L,
    mz = base::I(list(c(100, 101))),
    intensity = base::I(list(c(1, 2)))
  ))

  result <- filter_spec(spectra, min_peaks = 4, remove_precursor = FALSE)

  expect_s4_class(result, "Spectra")
  expect_length(result, 0)
})