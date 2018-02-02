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


###########################################################
context("Utilities")
###########################################################

test_that("Z_byname works as expected", {
  simple <- create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
  X_0 <- simple$X_0[[1]]
  X_T <- simple$X_T[[1]]
  Z_1 <- matrix(c(50.47438029, -25.23719014, 14.76280986,
                  19.31568569, 15.78206435, 24.90224996), byrow = TRUE, nrow = 2, ncol = 3,
                dimnames = list(c("subcat 1", "subcat 2"), c("factor 1", "factor 2", "factor 3"))) %>%
    setrowtype("subcat") %>% setcoltype("factor")
  Z_2 <- matrix(c(42.96021467, -56.79031473, 17.83010006,
                  0, 22.47128816, 32.52871184), byrow = TRUE, nrow = 2, ncol = 3,
                dimnames = list(c("subcat 1", "subcat 2"), c("factor 1", "factor 2", "factor 3"))) %>%
    setrowtype("subcat") %>% setcoltype("factor")
  Z_3 <- matrix(c(12.6549815, -39.30996299, 12.6549815,
                  41.35566084, 30.28867832, 41.35566084), byrow = TRUE, nrow = 2, ncol = 3,
                dimnames = list(c("subcat 1", "subcat 2"), c("factor 1", "factor 2", "factor 3"))) %>%
    setrowtype("subcat") %>% setcoltype("factor")

  expect_equal(Z_byname(X_0 = X_0, X_T = X_T), Z_1)
  expect_equal(Z_byname(X_0 = simple$X_0, X_T = simple$X_T),
               list(Z_1, Z_2, Z_3, Z_1, Z_2, Z_3))
  # Now try in the context of a data frame.
  simple2 <- simple %>%
    mutate(
      Z = Z_byname(X_0 = X_0, X_T = X_T)
    )
  expect_equal(simple2$Z, list(Z_1, Z_2, Z_3, Z_1, Z_2, Z_3))
})


test_that("Zij works as expected", {
  # Test degenerate cases.
  PN <- 9999 # Positive Number
  ANS <- 42 # the answer (modulo the sign)
  # Case 1
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = ANS, X_0ij = 0, X_Tij = PN), ANS)
  # Case 2
  expect_equal(LMDIR:::Zij(v_0i1 = ANS, v_Ti1 = 0, X_0ij = PN, X_Tij = 0), -ANS)
  # Case 3
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = PN, X_0ij = PN, X_Tij = PN), 0)
  # Case 4
  expect_equal(LMDIR:::Zij(v_0i1 = PN, v_Ti1 = 0, X_0ij = PN, X_Tij = PN), 0)
  # Case 5
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = 0, X_0ij = PN, X_Tij = PN), 0)
  # Case 6
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = 0, X_0ij = 0, X_Tij = 0), 0)
  # Case 7
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = 0, X_0ij = PN, X_Tij = 0), 0)
  # Case 8
  expect_equal(LMDIR:::Zij(v_0i1 = 0, v_Ti1 = 0, X_0ij = 0, X_Tij = PN), 0)

  simple <- create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
  X_T <- simple$X_T[[1]]
  X_0 <- simple$X_0[[1]]
  expect_equal(LMDIR:::Zij(1, 1, X_0 = X_0, X_T = X_T), 50.47438029)
  expect_equal(LMDIR:::Zij(1, 2, X_0 = X_0, X_T = X_T), -25.23719014)
  expect_equal(LMDIR:::Zij(1, 3, X_0 = X_0, X_T = X_T), 14.76280986)
  expect_equal(LMDIR:::Zij(2, 1, X_0 = X_0, X_T = X_T), 19.31568569)
  expect_equal(LMDIR:::Zij(2, 2, X_0 = X_0, X_T = X_T), 15.78206435)
  expect_equal(LMDIR:::Zij(2, 3, X_0 = X_0, X_T = X_T), 24.90224996)
})

###########################################################
context("Additive")
###########################################################



test_that("simple additive LMDI works as expected", {
  # Verify that grouping on time_colname fails.
  expect_error(create_simple_LMDI() %>% group_by(Country, Year) %>% lmdi(),
               "'Year' is a grouping variable, but you can't group on time_colname in argument .DF of collapse_to_matrices.")
  create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
})
