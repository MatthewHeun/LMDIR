# Contains tests for the LMDIR package.

# Need to put dplyr before testthat.
# If not, the "matches" function in dplyr overrides the "matches" function in testthat,
# and tests containing the string "(" don't work as expected.

library(dplyr)
library(tidyr)
library(rlang)
library(magrittr)
library(matsbyname)
library(matsindf)
library(testthat)

###########################################################
context("Linear LMDI")
###########################################################

test_that("linear LMDI works as expected", {
  res <- create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
  expect_equal(res$dV_agg, list(0, 100, 59, 99, 0, 100, 59, 99))
  expect_equal(res$dV_agg_cum, list(0, 100, 159, 258, 0, 100, 159, 258))
  expect_equal(res$dV[[1]], matrix(c(0, 0, 0),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$dV[[2]], matrix(c(69.79006598, -9.455125793, 39.66505981),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$dV[[3]], matrix(c(42.96021467, -34.31902656, 50.3588119),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$dV[[4]], matrix(c(54.01064234, -9.021284671, 54.01064234),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))

  expect_equal(res$dV[[1]], res$dV[[5]])
  expect_equal(res$dV[[2]], res$dV[[6]])
  expect_equal(res$dV[[3]], res$dV[[7]])
  expect_equal(res$dV[[4]], res$dV[[8]])

  expect_equal(res$dV_cum[[1]], res$dV[[1]])
  expect_equal(res$dV_cum[[2]], sum_byname(res$dV[[1]], res$dV[[2]]))
  expect_equal(res$dV_cum[[3]], sum_byname(res$dV_cum[[2]], res$dV[[3]]))
  expect_equal(res$dV_cum[[4]], sum_byname(res$dV_cum[[3]], res$dV[[4]]))
})


###########################################################
context("Multiplicative LMDI")
###########################################################

test_that("multiplicative LMDI works as expected", {
  res <- create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
  expect_equal(res$D_agg, list(1, 2.25, 1.327777778, 1.414225941, 1, 2.25, 1.327777778, 1.414225941))
  expect_equal(res$D_agg_cum, list(1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941,
                                   1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941))
  expect_equal(res$D[[1]], matrix(c(1, 1, 1),
                                  nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                      c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$D[[2]], matrix(c(1.761117821, 0.926191306, 1.379410116),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$D[[3]], matrix(c(1.229284572, 0.847970248, 1.273773912),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))
  expect_equal(res$D[[4]], matrix(c(1.208140223, 0.968911503, 1.208140223),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subsubcat"))

  expect_equal(res$D[[1]], res$D[[5]])
  expect_equal(res$D[[2]], res$D[[6]])
  expect_equal(res$D[[3]], res$D[[7]])
  expect_equal(res$D[[4]], res$D[[8]])

  expect_equal(res$D_cum[[1]], res$D[[1]])
  expect_equal(res$D_cum[[2]], hadamardproduct_byname(res$D[[1]], res$D[[2]]))
  expect_equal(res$D_cum[[3]], hadamardproduct_byname(res$D_cum[[2]], res$D[[3]]))
  expect_equal(res$D_cum[[4]], hadamardproduct_byname(res$D_cum[[3]], res$D[[4]]))
})