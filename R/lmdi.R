#' Log-Mean Divisia Index (LMDI) decomposition analysis
#'
#' Performs log-mean divisia index decomposition analysis on a suitably-formatted data frame.
#'
#' @param .lmdidata a grouped data frame.
#'        Group by columns of variables within which you want an LMDI analysis conducted.
#'        `time` should not be one of the grouping variables.
#' @param time the name of the column in `.lmdidata` that contains times at which
#'        data are available (a string).
#'        Default is "`Year`".
#' @param X the name (as a string) of a column in `.lmdidata` containing
#'        `X` matrices with
#'        named rows representing subcategories of the energy aggregate (`V`) and
#'        named columns representing factors contributing to changes in `V` over time.
#'        Default is "`X`".
#' @param fillrow a row vector of type `matrix` passed to [Z_byname()].
#'        (See [Z_byname()] for details.)
#' @param V the name for the `V` output column (a string).
#'        Default is "`V`".
#' @param Z the name for the `Z` output column (a string).
#'        Default is "`Z`".
#' @param deltaV the name for the `deltaV` output column (a string).
#'        Default is "`dV`".
#' @param D the name for the `D` output column (a string).
#'        Default is "`D`".
#'
#' @return a data frame containing several columns.
#'
#' @export
#'
lmdi <- function(.lmdidata, time = "Year", X = "X", fillrow = NULL,
                 # Output columns
                 V = "V", Z = "Z", deltaV = "dV", D = "D"){

  # Establish names for some intermediate columns.
  v_colname <- ".v"
  L_name <- ".L"
  w_colname <- ".w"

  # Establish names for new columns.
  zero_suffix <- "_0"
  T_suffix <- "_T"
  vec_suffix <- "_vec"
  X0_colname <- paste0(X, zero_suffix)
  v0_colname <- paste0(v_colname, zero_suffix)
  V0_colname <- paste0(V, zero_suffix)
  XT_colname <- paste0(X, T_suffix)
  vT_colname <- paste0(v_colname, T_suffix)
  VT_colname <- paste0(V, T_suffix)
  LV_colname <- paste0(L_name, "(", V, ")")

  # Ensure that all of these intermediate columns are missing.
  .lmdidata %>%
    matsindf::verify_cols_missing(c(deltaV, D, Z, v_colname, V, L_name, w_colname,
                                    X0_colname, v0_colname, V0_colname,
                                    XT_colname, vT_colname, VT_colname, LV_colname))

  # Ensure that time is NOT a grouping variable.
  if (time %in% dplyr::groups(.lmdidata)) {
    stop(paste0("'", time, "'", " is a grouping variable, but you can't group on ",
               "time in argument .lmdidata of collapse_to_matrices."))
  }

  XvV <- .lmdidata %>% dplyr::mutate(
    !!as.name(v_colname) := matsbyname::rowprods_byname(!!as.name(X)),
    !!as.name(V) := matsbyname::sumall_byname(!!as.name(v_colname))
  )

  # Create a data frame of metadata and X matrices, v column vectors, and V values
  # for time 0 and time T.
  XvV0T <- create0Tcolumns(XvV, time_colname = time,
                           X_colname = X, v_colname = v_colname, V_colname = V,
                           zero_suffix = zero_suffix, T_suffix = T_suffix)
  # Do year-by-year LMDI calcs.
  dVD <- XvV0T %>%
    dplyr::mutate(
      !!as.name(LV_colname) := matsbyname::logarithmicmean_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      !!as.name(Z) := Z_byname(X_0 = !!as.name(X0_colname), X_T = !!as.name(XT_colname),
                                       fillrow = fillrow),
      !!as.name(deltaV) := matsbyname::colsums_byname(!!as.name(Z)) %>% matsbyname::transpose_byname(),
      !!as.name(D) := matsbyname::quotient_byname(!!as.name(deltaV), !!as.name(LV_colname)) %>%
        matsbyname::exp_byname()
    )

  # Test to ensure that everything works as expected.
  # We can calculate deltaV and D in two ways.
  # The first way is called "raw" and is calculated with the
  # X data prior to any of the LMDI calculations.
  # The second way is called "dc" (decomposed) and is calculated with the
  # LMDI-decomposed results.
  # If all calculations have gone well, the "raw" and "dc" versions must be exactly the same.
  # We make these calculations here.
  raw_suffix <- "_raw"
  dc_suffix <- "_decomp"
  agg_suffix <- "_agg"
  cum_suffix <- "_cum"
  err_suffix <- "_err"
  dV_raw_colname <- paste0(deltaV, raw_suffix)
  dV_decomp_colname <- paste0(deltaV, dc_suffix)
  dV_agg_colname <- paste0(deltaV, agg_suffix)
  dV_agg_cum_colname <- paste0(dV_agg_colname, cum_suffix)
  dV_cum_colname <- paste0(deltaV, cum_suffix)
  dV_err_colname <- paste0(deltaV, err_suffix)
  D_raw_colname <- paste0(D, raw_suffix)
  D_decomp_colname <- paste0(D, dc_suffix)
  D_agg_colname <- paste0(D, agg_suffix)
  D_agg_cum_colname <- paste0(D_agg_colname, cum_suffix)
  D_cum_colname <- paste0(D, cum_suffix)
  D_err_colname <- paste0(D, err_suffix)
  # Verify that these columns are missing.
  dVD %>%
    matsindf::verify_cols_missing(c(dV_raw_colname, dV_decomp_colname, dV_agg_colname, dV_agg_cum_colname, dV_cum_colname, dV_err_colname,
                                    D_raw_colname, D_decomp_colname, D_agg_colname, D_agg_cum_colname, D_cum_colname, D_err_colname))
  chk <- dVD %>%
    dplyr::mutate(
      # The "raw" way of calaculating deltaV at each time comes from the "raw" data in X.
      !!as.name(dV_raw_colname) := matsbyname::difference_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      # The "decomp" way of calculating deltaV at each time comes after we calculate all deltaV's for all factors
      # The "raw" and "decomp" methods of calculating deltaV should be identical.
      !!as.name(dV_decomp_colname) := matsbyname::sumall_byname(!!as.name(deltaV)),
      # Calculate error column
      !!as.name(dV_err_colname) := matsbyname::difference_byname(!!as.name(dV_decomp_colname), !!as.name(dV_raw_colname)),
      # The "raw" way of calaculating D at each time comes from the "raw" data in X.
      !!as.name(D_raw_colname) := matsbyname::quotient_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      # The "decomp" way of calculating D at each time comes after we calculate all D's for all factors
      # The "raw" and "decomp" methods of calculating D should be identical.
      !!as.name(D_decomp_colname) := matsbyname::prodall_byname(!!as.name(D)),
      # Calculate D error column
      !!as.name(D_err_colname) := matsbyname::difference_byname(!!as.name(D_decomp_colname), !!as.name(D_raw_colname))
    )
  # If these assertions pass, the calculations are internally consistent.
  assertthat::assert_that(all(Map(f = all.equal, chk[[dV_raw_colname]], chk[[dV_decomp_colname]]) %>% as.logical()),
              msg = "dV_raw and dV_decomp are not all identical in lmdi()")
  assertthat::assert_that(all(Map(f = all.equal, chk[[D_raw_colname]], chk[[D_decomp_colname]]) %>% as.logical()),
              msg = "D_raw and D_decomp are not all identical in lmdi()")

  cumulatives <- chk %>%
    dplyr::select(!!!dplyr::group_vars(chk), !!as.name(time), !!as.name(Z),
                  !!as.name(dV_raw_colname), !!as.name(D_raw_colname),
                  !!as.name(deltaV), !!as.name(D)) %>%
    dplyr::rename(
      !!as.name(dV_agg_colname) := !!as.name(dV_raw_colname),
      !!as.name(D_agg_colname) := !!as.name(D_raw_colname)
    ) %>%
    dplyr::mutate(
      # These cululative sums and products are performed by group,
      # which is exactly what we want!
      !!as.name(dV_agg_cum_colname) := matsbyname::cumsum_byname(!!as.name(dV_agg_colname)),
      !!as.name(D_agg_cum_colname) := matsbyname::cumprod_byname(!!as.name(D_agg_colname)),
      !!as.name(dV_cum_colname) := matsbyname::cumsum_byname(!!as.name(deltaV)),
      !!as.name(D_cum_colname) := matsbyname::cumprod_byname(!!as.name(D))
    )

  # Now join the group_vars and Year column of .lmdidata and out by the group_vars and Year.
  XvV %>%
    dplyr::select(dplyr::group_vars(XvV), dplyr::all_of(time), X, V) %>%
    dplyr::left_join(cumulatives, by = c(dplyr::group_vars(.lmdidata), time))
}
