#' Log-Mean Divisia Index (LMDI) decomposition analysis
#'
#' Performs log-mean divisia index decomposition analysis on a suitably-formatted data frame.
#'
#' @param .lmdidata a grouped data frame.
#'        Group by columns of variables within which you want an LMDI analysis conducted.
#'        \code{time_colname} should not be one of the grouping variables.
#' @param time_colname the name of the column in \code{.lmdidata} that contains times at which
#'        data are available (a string).
#'        Default is "\code{Year}".
#' @param X_colname the name (as a string) of a column in \code{.lmdidata} containing
#'        \strong{\code{X}} matrices with
#'        named rows representing subcategories of the energy aggregate (\code{V}) and
#'        named columns representing factors contributing to changes in \code{V} over time.
#'        Default is "\code{X}".
#' @param deltaV_colname the name for the \code{deltaV} column (a string).
#'        Default is "\code{dV}".
#' @param D_colname the name for the \code{D} column (a string).
#'        Default is "\code{D}".
#'
#' @return a data frame containing several columns.
#'
#' @importFrom matsbyname elementquotient_byname
#' @importFrom matsbyname difference_byname
#' @importFrom matsbyname logarithmicmean_byname
#' @importFrom matsbyname rowprods_byname
#' @importFrom matsbyname colsums_byname
#' @importFrom dplyr groups
#' @importFrom dplyr left_join
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#'
#' @export
#'
lmdi <- function(.lmdidata, time_colname = "Year", X_colname = "X",
                 # Output columns
                 deltaV_colname = "dV", D_colname = "D"){

  # Establish names for some intermediate columns.
  Z_colname <- ".Z"
  v_colname <- ".v"
  V_colname <- ".V"
  L_name <- ".L"
  w_colname <- ".w"

  # Establish names for new columns.
  zero_suffix <- "_0"
  T_suffix <- "_T"
  vec_suffix <- "_vec"
  X0_colname <- paste0(X_colname, zero_suffix)
  v0_colname <- paste0(v_colname, zero_suffix)
  V0_colname <- paste0(V_colname, zero_suffix)
  XT_colname <- paste0(X_colname, T_suffix)
  vT_colname <- paste0(v_colname, T_suffix)
  VT_colname <- paste0(V_colname, T_suffix)
  LV_colname <- paste0(L_name, "(", V_colname, ")")

  # Ensure that time_colname is NOT a grouping variable.
  if (time_colname %in% groups(.lmdidata)) {
    stop(paste0("'", time_colname, "'", " is a grouping variable, but you can't group on ",
               "time_colname",
               " in argument .lmdidata of collapse_to_matrices."))
  }

  XvV <- .lmdidata %>% mutate(
    !!as.name(v_colname) := rowprods_byname(!!as.name(X_colname)),
    !!as.name(V_colname) := sumall_byname(!!as.name(v_colname))
  )

  # Create a data frame of metadata and X matrices, v column vectors, and V values
  # for time 0 and time T.
  XvV0T <- create0Tcolumns(XvV, time_colname = time_colname,
                           X_colname = X_colname, v_colname = v_colname, V_colname = V_colname,
                           zero_suffix = zero_suffix, T_suffix = T_suffix)
  # Do year-by-year LMDI calcs.
  dVD <- XvV0T %>%
    mutate(
      !!as.name(LV_colname) := logarithmicmean_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      !!as.name(Z_colname) := Z_byname(X_0 = !!as.name(X0_colname), X_T = !!as.name(XT_colname)),
      !!as.name(deltaV_colname) := colsums_byname(!!as.name(Z_colname)) %>% transpose_byname(),
      # !!as.name(w_colname) := elementquotient_byname(!!as.name(Lv_colname), !!as.name(LV_colname)),
      !!as.name(D_colname) := elementquotient_byname(!!as.name(Z_colname), !!as.name(LV_colname)) %>%
        colsums_byname() %>% transpose_byname() %>% elementexp_byname()
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
  dV_raw_colname <- paste0(deltaV_colname, raw_suffix)
  dV_decomp_colname <- paste0(deltaV_colname, dc_suffix)
  dV_agg_colname <- paste0(deltaV_colname, agg_suffix)
  dV_agg_cum_colname <- paste0(dV_agg_colname, cum_suffix)
  dV_cum_colname <- paste0(deltaV_colname, cum_suffix)
  dV_err_colname <- paste0(deltaV_colname, err_suffix)
  D_raw_colname <- paste0(D_colname, raw_suffix)
  D_decomp_colname <- paste0(D_colname, dc_suffix)
  D_agg_colname <- paste0(D_colname, agg_suffix)
  D_agg_cum_colname <- paste0(D_agg_colname, cum_suffix)
  D_cum_colname <- paste0(D_colname, cum_suffix)
  D_err_colname <- paste0(D_colname, err_suffix)
  chk <- dVD %>%
    mutate(
      # The "raw" way of calaculating deltaV at each time comes from the "raw" data in X.
      !!as.name(dV_raw_colname) := difference_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      # The "decomp" way of calculating deltaV at each time comes after we calculate all deltaV's for all factors
      # The "raw" and "decomp" methods of calculating deltaV should be identical.
      !!as.name(dV_decomp_colname) := sumall_byname(!!as.name(deltaV_colname)),
      # Calculate error column
      !!as.name(dV_err_colname) := difference_byname(!!as.name(dV_decomp_colname), !!as.name(dV_raw_colname)),
      # The "raw" way of calaculating D at each time comes from the "raw" data in X.
      !!as.name(D_raw_colname) := elementquotient_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      # The "decomp" way of calculating D at each time comes after we calculate all D's for all factors
      # The "raw" and "decomp" methods of calculating D should be identical.
      !!as.name(D_decomp_colname) := prodall_byname(!!as.name(D_colname)),
      # Calculate D error column
      !!as.name(D_err_colname) := difference_byname(!!as.name(D_decomp_colname), !!as.name(D_raw_colname))
    )
  # If these tests pass, the calculations are internally consistent.
  if (!all(Map(f = all.equal, chk[[dV_raw_colname]], chk[[dV_decomp_colname]]) %>% as.logical)) {
    # There is a problem. dV_raw and dV_decomp should be identical.
    warning("dV_raw and dV_decomp are not all identical in lmdi()")
  }
  if (!all(Map(f = all.equal, chk[[D_raw_colname]], chk[[D_decomp_colname]]) %>% as.logical)) {
    # There is a problem. D_raw and D_decomp should be identical.
    warning("D_raw and D_decomp are not all identical in lmdi()")
  }

  cumulatives <- chk %>%
    select(!!!group_vars(chk), !!as.name(time_colname),
           !!as.name(dV_raw_colname), !!as.name(D_raw_colname),
           !!as.name(deltaV_colname), !!as.name(D_colname)) %>%
    rename(
      !!as.name(dV_agg_colname) := !!as.name(dV_raw_colname),
      !!as.name(D_agg_colname) := !!as.name(D_raw_colname)
    ) %>%
    mutate(
      # These cululative sums and products are performed by group,
      # which is exactly what we want!
      !!as.name(dV_agg_cum_colname) := cumsum_byname(!!as.name(dV_agg_colname)),
      !!as.name(D_agg_cum_colname) := cumprod_byname(!!as.name(D_agg_colname)),
      !!as.name(dV_cum_colname) := cumsum_byname(!!as.name(deltaV_colname)),
      !!as.name(D_cum_colname) := cumprod_byname(!!as.name(D_colname))
    )

  # Now join the group_vars and Year column of .lmdidata and out by the group_vars and Year.
  .lmdidata %>%
    select(group_vars(.lmdidata), time_colname) %>%
    left_join(cumulatives, by = c(group_vars(.lmdidata), time_colname))
}