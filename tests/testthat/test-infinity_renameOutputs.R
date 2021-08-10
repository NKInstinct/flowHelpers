# Make sure that this sucker can properly detect old filepaths and then rename
# them. Obviously, don't run these tests with trialRun = FALSE or else it'll
# change all the test files and then will fail on a subsequent test!
dir <- system.file("extdata", "sampleInfinity", package = "flowHelpers")
paths <- infinity_renameOutputs(dir, pattern = "\\.fcs", concat_tag_to_remove = "/FCS/concat",trialRun = TRUE)

test_that("Old paths contains the FCS/concat string",{
  expect_equal(str_detect(paths$`Old paths`, "infinityOutput/.*/FCS/concat/"), c(TRUE, TRUE))
})

test_that("New paths do not contain the FCS/concat string",{
  expect_equal(str_detect(paths$`New paths`, "/FCS/concat"), c(FALSE, FALSE))
})

test_that("Files are renamed appropriately", {
  expect_equal(str_detect(paths$`New paths`, "/Sample.\\.fcs$"), c(TRUE, TRUE))
})