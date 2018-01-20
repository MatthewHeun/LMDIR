# Contains tests for the LMDIR package.

# Need to put dplyr before testthat.
# If not, the "matches" function in dplyr overrides the "matches" function in testthat,
# and tests containing the string "(" don't work as expected.

library(dplyr)
library(tidyr)
library(rlang)
library(magrittr)
library(mosaic)
library(testthat)


###########################################################
context("additive")
###########################################################

test_that("simple additive LMDI works as expected", {
  # Create a simple example
  DF <- data.frame(Year = c(1950, 2000, 2013),
                   x11 = c(1, 4, 8),
                   x21 = c(10, 5, 2),
                   x31 = c(2, 3, 4),
                   x12 = c(4, 5, 5),
                   x22 = c(5, 6, 7),
                   x32 = c(3, 4, 5))
  result <- lmdi(DF, agg_structure = list(V1 = c("x11", "x21", "x31"), V2 = c("x12", "x22", "x32")))

})