test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

# Tests for readFS: if given a directory with FCS files, do you return a flowSet?

# Tests for applyComp: if given a comp matrix in .csv format, do you return a comp object?
# If given a flowset with an acquisition-defined matrix, can you extract it?

# Tests for getLinearChannels: if given a length-1 number and a gatingset, do you return the correct vector of channels?
# If given a vector of linear channels, do you return the correct vector of non-linear channels?

# Tests for applyTransform: if given a gs and linear channels, do you return a transformed gs?
# note - this one is a little tough, but if you save the output from
# flowWorkspace::gh_get_transformations to a variable you can at least test that
# the length of the resulting list is the number of expected non-linear
# channels.
