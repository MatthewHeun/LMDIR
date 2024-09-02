# Contains tests for the LMDIR package.

test_that("additive LMDI-I works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi()
  # Values tested below come from data-raw/LMDI_testing.xslx
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


test_that("multiplicative LMDI-I works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi()
  # Values tested below come from data-raw/LMDI_testing.xslx
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


test_that("additive LMDI-II works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi(weights="LMDI-II")
  # Values tested below come from data-raw/LMDI_testing.xslx
  expect_equal(res$dV_agg, list(0, 100, 59, 99, 0, 100, 59, 99))
  expect_equal(res$dV_agg_cum, list(0, 100, 159, 258, 0, 100, 159, 258))
  expect_equal(res$dV[[1]], matrix(c(0, 0, 0),
                                      nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                          c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[2]], matrix(c(69.18298477, -8.87773383,	39.69474906),
                                      nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                          c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[3]], matrix(c(43.25676087,	-34.72232042,	50.46555955),
                                      nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                          c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$dV[[4]], matrix(c(54.45209029,	-9.90418057,	54.45209029),
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


test_that("multiplicative-II LMDI works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi(weights="LMDI-II")
  # Values tested below come from data-raw/LMDI_testing.xslx
  expect_equal(res$D_agg, list(1, 2.25, 1.327777778, 1.414225941, 1, 2.25, 1.327777778, 1.414225941))
  expect_equal(res$D_agg_cum, list(1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941,
                                   1, 2.25, 2.25*1.327777778, 2.25*1.327777778*1.414225941))
  expect_equal(res$D[[1]], matrix(c(1, 1, 1),
                                     nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                         c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[2]], matrix(c(1.752469135,	0.930538130,	1.379742261),
                                     nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                         c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[3]], matrix(c(1.231037506,	0.846328552,	1.274427454),
                                     nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                         c("subsubcat"))) %>%
                 matsbyname::setrowtype("factor") %>% matsbyname::setcoltype("subsubcat"))
  expect_equal(res$D[[4]], matrix(c(1.210008769,	0.965921347,	1.210008769),
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


test_that("argument for weights works as expected", {
  res <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi()

  res_I <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi(weights="LMDI-I")

  res_II <- create_simple_LMDI() %>%
    dplyr::group_by(Country) %>%
    lmdi(weights="LMDI-II")

  # res and res_I should yield the same results, because "LMDI-I" is the default value for weights
  expect_equal(res, res_I)
  expect_false(identical(res, res_II))
})


test_that("LMDI works as expected with a custom fillrow", {
  simple_lmdi <- create_simple_LMDI() %>%
    dplyr::group_by(Country)
  fillrow <- matrix(1e-10, nrow = 1, ncol = 3, dimnames = list("row", c("factor 1", "factor 2", "factor 3")))

  res <- lmdi(simple_lmdi, fillrow = fillrow)

  # Do a simple check to ensure that things worked correctly.
  expect_equal(res$D[[1]], res$D[[5]])
})


test_that("wrong naming for weights returns an error", {
  expect_error(
    create_simple_LMDI() %>%
      dplyr::group_by(Country) %>%
      lmdi(weights = "IDML-I")
  )
})



