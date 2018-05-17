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

test_that("Zij works as expected", {
  # Test degenerate cases.
  PN <- 9999 # Positive Number
  ANS <- 42 # the answer (modulo the sign)
  #
  # The cases below are taken from Table 2, p. 492 in
  # B. Ang, F. Zhang, and K.-H. Choi.
  # Factorizing changes in energy and environmental indicators through decomposition.
  # Energy, 23(6):489–495, Jun 1998.
  #
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

  simple <- create_simple_LMDI()

  X_0 <- simple$X[[1]]
  X_T <- simple$X[[2]]
  expect_equal(LMDIR:::Zij(1, 1, X_0 = X_0, X_T = X_T), 50.47438029)
  expect_equal(LMDIR:::Zij(1, 2, X_0 = X_0, X_T = X_T), -25.23719014)
  expect_equal(LMDIR:::Zij(1, 3, X_0 = X_0, X_T = X_T), 14.76280986)
  expect_equal(LMDIR:::Zij(2, 1, X_0 = X_0, X_T = X_T), 19.31568569)
  expect_equal(LMDIR:::Zij(2, 2, X_0 = X_0, X_T = X_T), 15.78206435)
  expect_equal(LMDIR:::Zij(2, 3, X_0 = X_0, X_T = X_T), 24.90224996)
})

test_that("Z_byname works as expected", {
  simple <- create_simple_LMDI()

  X_T <- simple$X[[2]]
  X_0 <- simple$X[[1]]

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
  expect_equal(Z_byname(X_0 = simple$X[1:3], X_T = simple$X[2:4]),
               list(Z_1, Z_2, Z_3))
  # Now try in the context of a data frame.
  simple2 <- simple %>%
    mutate(
      X_0 = list(simple$X[[1]], simple$X[[2]], simple$X[[3]], NULL,
                 simple$X[[5]], simple$X[[6]], simple$X[[7]], NULL),
      X_T = list(simple$X[[2]], simple$X[[3]], simple$X[[4]], NULL,
                 simple$X[[6]], simple$X[[7]], simple$X[[8]], NULL)
    ) %>%
    filter(Year != 1974) %>%
    mutate(
      Z = Z_byname(X_0 = X_0, X_T = X_T)
    )
  expect_equal(simple2$Z, list(Z_1, Z_2, Z_3, Z_1, Z_2, Z_3))

  # Try with some degenerate values in the X matrix.
  # When a 0 is present, the row product is zero, making the logarithms blow up.
  # These are the conditions under which we need to refer to
  # Table 2, p. 492 of
  # B.W. Ang and F.Q. Zhang and Ki-Hong Choi, 1998,
  # Factorizing changes in energy and environmental indicators through decomposition,
  # Energy, Volume 23, Number 6, pp. 489-495.
  # We employ the method suggested by
  # R. Wood and M. Lenzen. 2006.
  # Zero-value problems of the logarithmic mean divisia index decomposition method.
  # Energy Policy, volume 34, number 12, pp. 1326–1331.

  # First, set one of the values in X_0 to 0.
  X_0_2 <- X_0
  X_0_2[[1, 1]] <- 0
  # In this situation, the second row of Z remains same as Z_1 above.
  # However, with X_0_11 = 0, we also get v_0_11 = 0.
  # For Z_11, we have Case 0, and Z_11 = v_T_11 = 60.
  # For Z_12 and Z_13, we have Case 3, and both are 0.
  Z_expected_1 <- matrix(c(60, 0, 0,
                           19.31568569, 15.78206435, 24.90224996), byrow = TRUE, nrow = 2, ncol = 3,
                         dimnames = list(c("subcat 1", "subcat 2"), c("factor 1", "factor 2", "factor 3"))) %>%
    setrowtype("subcat") %>% setcoltype("factor")
  expect_equal(Z_byname(X_0 = X_0_2, X_T = X_T), Z_expected_1)

  X_T_2 <- X_T
  X_T_2[[2, 3]] <- 0
  # In this situation, the first row of Z remains same as Z_1 above.
  # However, with X_T_23 = 0, we also get v_T_21 = 0.
  # For Z_21 and Z_22, we have Case 4, and both are 0.
  # For Z_23, we have Case 2, and we obtain Z_23 = -v_0_21 = -60.
  Z_expected_2 <- matrix(c(50.47438029, -25.23719014, 14.76280986,
                           0, 0, -60), byrow = TRUE, nrow = 2, ncol = 3,
                         dimnames = list(c("subcat 1", "subcat 2"), c("factor 1", "factor 2", "factor 3"))) %>%
    setrowtype("subcat") %>% setcoltype("factor")
  expect_equal(Z_byname(X_0 = X_0, X_T = X_T_2), Z_expected_2)

})


###########################################################
context("Group error")
###########################################################

test_that("errors are given when grouping errors are present", {
  # Verify that grouping on time_colname fails.
  expect_error(create_simple_LMDI() %>% group_by(Country, Year) %>% lmdi(),
               "'Year' is a grouping variable, but you can't group on time_colname in argument .lmdidata of collapse_to_matrices.")
})


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
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$dV[[2]], matrix(c(69.79006598, -9.455125793, 39.66505981),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$dV[[3]], matrix(c(42.96021467, -34.31902656, 50.3588119),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$dV[[4]], matrix(c(54.01064234, -9.021284671, 54.01064234),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))

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
                                                                      c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$D[[2]], matrix(c(1.761117821, 0.926191306, 1.379410116),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$D[[3]], matrix(c(1.229284572, 0.847970248, 1.273773912),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))
  expect_equal(res$D[[4]], matrix(c(1.208140223, 0.968911503, 1.208140223),
                                   nrow = 3, ncol = 1, dimnames = list(c("factor 1", "factor 2", "factor 3"),
                                                                       c("subcat"))) %>%
                 setrowtype("factor") %>% setcoltype("subcat"))

  expect_equal(res$D[[1]], res$D[[5]])
  expect_equal(res$D[[2]], res$D[[6]])
  expect_equal(res$D[[3]], res$D[[7]])
  expect_equal(res$D[[4]], res$D[[8]])

  expect_equal(res$D_cum[[1]], res$D[[1]])
  expect_equal(res$D_cum[[2]], elementproduct_byname(res$D[[1]], res$D[[2]]))
  expect_equal(res$D_cum[[3]], elementproduct_byname(res$D_cum[[2]], res$D[[3]]))
  expect_equal(res$D_cum[[4]], elementproduct_byname(res$D_cum[[3]], res$D[[4]]))
})


###########################################################
context("Preserve grouping")
###########################################################

test_that("Preserving grouping works as expected", {
  res <- create_simple_LMDI() %>%
    mutate(
      groupingcol = paste(Country, "group")
    ) %>%
    group_by(Country, groupingcol) %>%
    lmdi()
  expect_equal(group_vars(res), c("Country", "groupingcol"))
})


###########################################################
context("First row 0s and 1s")
###########################################################

test_that("First row contains 0s and 1s", {
  res <- create_simple_LMDI() %>%
    group_by(Country) %>%
    lmdi()
  expect_equal(res$dV_agg[[1]], 0)
  expect_equal(res$D_agg[[1]], 1)
  expect_equal(res$dV_agg_cum[[1]], 0)
  expect_equal(res$D_agg_cum[[1]], 1)
  dV0 <- matrix(c(0, 0, 0), nrow = 3, ncol = 1,
                dimnames = list(c("factor 1", "factor 2", "factor 3"), c("subcat"))) %>%
    setrowtype("factor") %>% setcoltype("subcat")
  D0 <- matrix(c(1, 1, 1), nrow = 3, ncol = 1,
               dimnames = list(c("factor 1", "factor 2", "factor 3"), c("subcat"))) %>%
    setrowtype("factor") %>% setcoltype("subcat")
  expect_equal(res$dV[[1]], dV0)
  expect_equal(res$D[[1]], D0)
  expect_equal(res$dV_cum[[1]], dV0)
  expect_equal(res$D_cum[[1]], D0)
})


###########################################################
context("Fillrow")
###########################################################

test_that("fillrow option works as expected on Z_byname", {
  # Create X_0 and X_T matrices that have different rows.
  # This example comes from Ghana's energy LMDI in the years 2003 and 2004.
  dn <- list(c("HTH.600.C - Electric heaters", "KE - Fans"),
             c("E.ktoe", "eta_ij", "phi_i", "phi_ij"))
  X_0 <- matrix(c(7909.898576, 0.168054745, 0.521283202, 0.009258785,
                  7909.898576, 0.072489945, 0.078700891, 0.036152202), byrow = TRUE, nrow = 2, ncol = 4,
                dimnames = dn) %>%
    setrowtype("categories") %>% setcoltype("factors")
  X_T <- matrix(c(7962.921168, 0.101321321, 0.059104816, 0.036042366), byrow = TRUE, nrow = 1, ncol = 4,
                dimnames = list("KE - Fans", dn[[2]])) %>%
    setrowtype("categories") %>% setcoltype("factors")
  # Z1 should be a 2-row matrix formed by assuming small numbers for all of the missing values.
  Z1 <- Z_byname(X_0 = X_0, X_T = X_T)
  expect_equal(nrow(Z1), 2)
  expect_equal(rownames(Z1), c("HTH.600.C - Electric heaters", "KE - Fans"))
  expect_equal(Z1, matrix(c(-2.185092, -1.450439688, -1.527733412, -1.252513979,
                            0.011188553, 0.56076961, -0.479535332, -0.005095744),
                          byrow = TRUE, nrow = 2, ncol = 4, dimnames = dn) %>%
                 setrowtype("categories") %>% setcoltype("factors"),
               tolerance = 1e-6)
  # Now try with a fillrow argument.
  # The following fillrow value sets ONLY the allocation from subcategory to subsubcategory to 0.
  # Doing so kicks the algorithm into another state compared to the previous example.
  # In this one, we go to case 2 of Table 2, p. 492 in Ang et al. 1998.
  # These conditions are what you find when an energy type disappears at a later year.
  # In GH, HTH.600.C was present in 2003 but disappeared in 2004 due to
  # shutdown of the VALCO smelters.
  fr <- matrix(c(42, 42, 42, 0), nrow = 1, ncol = 4,
                    dimnames = list("row", c("E.ktoe", "eta_ij", "phi_i", "phi_ij"))) %>%
    setrowtype("categories") %>% setcoltype("factors")
  Z2 <- Z_byname(X_0 = X_0, X_T = X_T, fillrow = fr)
  expect_equal(Z2, matrix(c(0, 0, 0, -6.415779079,
                            0.011188553, 0.56076961, -0.479535332, -0.005095744),
                          byrow = TRUE, nrow = 2, ncol = 4, dimnames = dn) %>%
                 setrowtype("categories") %>% setcoltype("factors"),
               tolerance = 1e-6)
  # If we switch the order of X_0 and X_T, we kick to case 1 of of Table 2, p. 492 in Ang et al. 1998.
  # These conditions are what you find when an energy type turns appears in a subsequent year.
  Z3 <- Z_byname(X_0 = X_T, X_T = X_0, fillrow = fr)
  expect_equal(Z3, matrix(c(0, 0, 0, 6.415779079,
                            -0.011188553, -0.56076961, 0.479535332, 0.005095744),
                          byrow = TRUE, nrow = 2, ncol = 4, dimnames = dn) %>%
                 setrowtype("categories") %>% setcoltype("factors"),
               tolerance = 1e-6)

  # Ensure that fillrow works properly from the lmdi method.
  # First using the small values approach.
  DF1 <- data.frame(Year = c(2003, 2004))
  DF1$X <- list(X_0, X_T)
  res1 <- lmdi(DF1, time_colname = "Year", X_colname = "X")
  expect_equal(res1$dV_agg[[1]], 0)
  expect_equal(res1$dV_agg[[2]], -6.328451992, tolerance = 1e-6)
  expect_equal(res1$dV[[1]], matrix(0, nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))
  expect_equal(res1$dV[[2]], matrix(c(-2.173903446, -0.889670079, -2.007268743, -1.257609724),
                                   nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"),
               tolerance = 1e-6)

  expect_equal(res1$D_agg[[1]], 1)
  expect_equal(res1$D_agg[[2]], 0.213582281)
  expect_equal(res1$D[[1]], matrix(1, nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))
  expect_equal(res1$D[[2]], matrix(c(0.588433187, 0.804912278, 0.612844653, 0.7358158),
                                  nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))

  # Now using a fillrow.
  DF2 <- data.frame(Year = c(2003, 2004))
  DF2$X <- list(X_0, X_T)
  res2 <- lmdi(DF2, time_colname = "Year", X_colname = "X", fillrow = fr)
  expect_equal(res2$dV_agg[[1]], 0)
  expect_equal(res2$dV_agg[[2]], -6.328451992, tolerance = 1e-6)
  expect_equal(res2$dV[[1]], matrix(0, nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))
  expect_equal(res2$dV[[2]], matrix(c(0.011188553, 0.56076961, -0.479535332, -6.420874823),
                                   nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"),
               tolerance = 1e-6)

  expect_equal(res2$D_agg[[1]], 1)
  expect_equal(res2$D_agg[[2]], 0.213582281)
  expect_equal(res2$D[[1]], matrix(1, nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))
  expect_equal(res2$D[[2]], matrix(c(1.002733012, 1.146589093, 0.889606884, 0.208820902),
                                  nrow = 4, ncol = 1, dimnames = list(dn[[2]], "categories")) %>%
                 setrowtype("factors") %>% setcoltype("categories"))
})


###########################################################
# context("Interactive")
###########################################################

# test_that("interactive programming works, too", {
#   res <- create_simple_LMDI() %>%
#     group_by(Country) %>%
#     lmdi(time_colname = Year)
#
# })
