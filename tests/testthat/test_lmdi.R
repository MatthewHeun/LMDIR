# Contains tests for the LMDIR package.

# Need to put dplyr before testthat.
# If not, the "matches" function in dplyr overrides the "matches" function in testthat,
# and tests containing the string "(" don't work as expected.

library(dplyr)
library(tidyr)
library(rlang)
library(magrittr)
library(mosaic)
library(byname)
library(matsindf)
library(testthat)


###########################################################
context("Additive")
###########################################################

create_simple_LMDI <- function(){
  simple.factor.names <- c("factor 1", "factor 2", "factor 3")
  simple.subcat.names <- c("subcat 1", "subcat 2")
  # Create a tidy data frame of x values
  # which will be collapsed to matrices
  AB <- data.frame(Year = rep.int(1971, 6), x = c(1, 10, 2, 4, 5, 3)) %>%
    rbind(data.frame(Year = rep.int(1972, 6), x = c(4, 5, 3, 5, 6, 4))) %>%
    rbind(data.frame(Year = rep.int(1973, 6), x = c(8, 2, 4, 5, 7, 5))) %>%
    rbind(data.frame(Year = rep.int(1974, 6), x = c(10, 1, 5, 6, 8, 6))) %>%
    mutate(
      matnames = "X",
      rowtypes = "subcat",
      coltypes = "factor",
      rownames = rep.int(c(rep.int(simple.subcat.names[[1]], 3), rep.int(simple.subcat.names[[2]], 3)), 4),
      colnames = rep.int(c(rep.int(simple.factor.names, 2)), 4),
      Country = "AB"
    ) %>%
    group_by(Country, Year) %>%
    collapse_to_matrices(matnames = "matnames", values = "x",
                         rownames = "rownames", colnames = "colnames",
                         rowtypes = "rowtypes", coltypes = "coltypes") %>%
    rename(X = x)
  rbind(AB, AB %>% mutate(Country = "YZ"))
}


test_that("simple additive LMDI works as expected", {
  # Verify that grouping on time_colname fails.
  expect_error(create_simple_LMDI() %>% group_by(Country, Year) %>% lmdi(),
               "'Year' is a grouping variable, but you can't group on time_colname in argument .DF of collapse_to_matrices.")
  create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
})
