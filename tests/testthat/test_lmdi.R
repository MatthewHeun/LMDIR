
test_that("linear LMDI works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi()
  expect_equal(res$dV_agg, list(0, 100, 59, 99, 0, 100, 59, 99))
  expect_equal(res$dV_agg_cum, list(0, 100, 159, 258, 0, 100, 159, 258))
  expect_equal(res$dV[[1]], matrix(c(0, 0, 0),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[2]], matrix(c(69.79006598, -9.455125793, 39.66505981),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[3]], matrix(c(42.96021467, -34.31902656, 50.3588119),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[4]], matrix(c(54.01064234, -9.021284671, 54.01064234),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))

  expect_equal(res$dV[[1]], res$dV[[5]])
  expect_equal(res$dV[[2]], res$dV[[6]])
  expect_equal(res$dV[[3]], res$dV[[7]])
  expect_equal(res$dV[[4]], res$dV[[8]])

  expect_equal(res$dV_cum[[1]], res$dV[[1]])
  expect_equal(res$dV_cum[[2]], matsbyname::sum_byname(res$dV[[1]], res$dV[[2]]))
  expect_equal(res$dV_cum[[3]], matsbyname::sum_byname(res$dV_cum[[2]], res$dV[[3]]))
  expect_equal(res$dV_cum[[4]], matsbyname::sum_byname(res$dV_cum[[3]], res$dV[[4]]))
})


test_that("multiplicative LMDI works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi()
  expect_equal(res$D_agg, list(1, 2.25, 1.327777778, 1.414225941, 1, 2.25, 1.327777778, 1.414225941))
  expect_equal(res$D_agg_cum, list(1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941,
                                   1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941))
  expect_equal(res$D[[1]], matrix(c(1, 1, 1),
                                  nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                      c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[2]], matrix(c(1.761117821, 0.926191306, 1.379410116),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[3]], matrix(c(1.229284572, 0.847970248, 1.273773912),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[4]], matrix(c(1.208140223, 0.968911503, 1.208140223),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))

  expect_equal(res$D[[1]], res$D[[5]])
  expect_equal(res$D[[2]], res$D[[6]])
  expect_equal(res$D[[3]], res$D[[7]])
  expect_equal(res$D[[4]], res$D[[8]])

  expect_equal(res$D_cum[[1]], res$D[[1]])
  expect_equal(res$D_cum[[2]], matsbyname::hadamardproduct_byname(res$D[[1]], res$D[[2]]))
  expect_equal(res$D_cum[[3]], matsbyname::hadamardproduct_byname(res$D_cum[[2]], res$D[[3]]))
  expect_equal(res$D_cum[[4]], matsbyname::hadamardproduct_byname(res$D_cum[[3]], res$D[[4]]))
})


test_that("LMDI works as expected with a custom fillrow", {
  simple_lmdi <- create_simple_LMDI() %>%
    dplyr::group_by(Country)
  fillrow <- matrix(1e-10, nrow = 1, ncol = 3, dimnames = list("row", c("factor 1", "factor 2", "factor 3")))

  res <- lmdi(simple_lmdi, fillrow = fillrow)

  # Do a simple check to ensure that things worked correctly.
  expect_equal(res$D[[1]], res$D[[5]])
})

