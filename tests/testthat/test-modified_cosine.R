test_that("modified_cosine validates inputs", {
  spectra <- Spectra::Spectra(
    backend = Spectra::MsBackendMemory(),
    data = data.frame(
      spectrum_id = 1:2,
      precursorMz = c(100, 100),
      mz = I(list(c(50, 60), c(50, 60))),
      intensity = I(list(c(1, 2), c(2, 1)))
    )
  )

  expect_error(modified_cosine(data.frame()), "Spectra object")
  expect_error(modified_cosine(spectra, tolerance = c(0, 1)), "single non-negative numeric")
  expect_error(modified_cosine(spectra, ppm = -1), "single non-negative numeric")
})
