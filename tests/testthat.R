library("testthat")
library("guidedPLS")

options(testthat.use_colours = FALSE)

# Basic usage
test_file("testthat/test_PLSSVD.R")
test_file("testthat/test_sPLSDA.R")
test_file("testthat/test_guidedPLS.R")
test_file("testthat/test_toyModel.R")
