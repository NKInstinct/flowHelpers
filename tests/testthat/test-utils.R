path_to_data <- system.file("extdata", package = "flowHelpers")

test_that("fcs files return correct flowSet classes", {
  expect_s4_class(readFS(FCSDirectory = path_to_data,
                         pattern = ".fcs$",
                         recursive = FALSE,
                         ncdf = FALSE),
                  "flowSet")
  expect_s4_class(readFS(FCSDirectory = path_to_data,
                         pattern = ".fcs$",
                         recursive = FALSE,
                         ncdf = TRUE),
                  "ncdfFlowSet")
})

compCSV <- system.file("extdata", "comp.csv", package = "flowHelpers")

compStandard <- readr::read_csv(compCSV)
compStandard <- dplyr::select(compStandard, -X1)
compStandard <- as.matrix(compStandard)

fs <- readFS(FCSDirectory = path_to_data,
             pattern = ".fcs$",
             recursive = FALSE,
             ncdf = FALSE)

test_that("comp matrices are correctly applied from specified source", {
  expect_equal(applyComp(fs, "acquisition"), compStandard)
  expect_s4_class(applyComp(fs, compCSV), "compensation")
})

gs <- flowWorkspace::load_gs(system.file("extdata", "gatingSet.gs", package = "flowHelpers"))

nonLin <- c("FITC-A", "BV421-A", "APC-A", "PE-Cy5-A", "PE-CF594-A")

test_that("each way of specifying channels gives back a length-5 list
          (the number of non-linear channels in the set)", {
            expect_equal(getLinearChannels(gs, 6), nonLin)
            expect_equal(getLinearChannels(gs, c("FSC-A",
                                                 "FSC-H",
                                                 "FSC-W",
                                                 "SSC-A",
                                                 "SSC-H",
                                                 "SSC-W",
                                                 "Time")), nonLin)
          })

fs <- readFS(path_to_data,
             pattern = ".fcs$",
             recursive = FALSE,
             ncdf = FALSE)

gs <- flowWorkspace::GatingSet(fs)

test_that("applying a transform returns a length-5 list", {
  expect_equal(length(applyTransform(gs, 6, "biex", 150)), 5)
  expect_equal(length(applyTransform(gs, 6, "arcsinh", 150)), 5)
})
