path_to_data <- system.file("extdata", package = "flowHelpers")

test_that("gs_create creates a valid gatingset",{
  expect_s4_class(gs_create(path_to_data,
                            ncdf = FALSE,
                            linearChannels = 6),
                  "GatingSet")
})
