# This file is part of the standard testthat setup.
# It ensures tests are run when R CMD check is executed.

library(testthat)
library(neurothresh)

test_check("neurothresh")
