#' Log-Mean Divisia Index (LMDI) decomposition analysis
#'
#' Performs log-mean divisia index decomposition analysis on a suitably-formatted data frame.
#'
#' @param .lmdidata a grouped data frame.
#'        Group by columns of variables within which you want an LMDI analysis conducted.
#'        \code{time_colname} should not be one of the grouping variables.
#' @param time_colname the name of the column in \code{.data} that contains times at which
#'        data are available (a string).
#'        Default is "\code{Year}".
#' @param X_colname the name (as a string) of a column in \code{.DF} containing
#'        \strong{\code{X}} matrices with
#'        named rows representing subcategories of the energy aggregate (\code{V}) and
#'        named columns representing factors contributing to changes in \code{V} over time.
#' @param pad either "\code{tail}" or "\code{head}" to indicate whether the first or last row,
#'        respectively, should contain \code{pad.value}.
#' @param pad.value the value to be used in the padding row.  Default is \code{NA}.
#' @param D_colname the name for the \code{D} column (a string).
#' @param deltaV_colname the name for the \code{deltaV} column (a string).
#' @param F_colname the name for the \code{F} column (a string).
#'
#' @return a data frame containing several columns.
#'
#' @importFrom byname elementquotient_byname
#' @importFrom byname difference_byname
#' @importFrom byname logarithmicmean_byname
#' @importFrom byname rowprods_byname
#' @importFrom byname colsums_byname
#' @importFrom dplyr groups
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom magrittr %>%
#'
#' @export
#'
lmdi <- function(.lmdidata, time_colname = "Year", X_colname = "X",
                 pad = c("tail", "head"), pad.value = NA,
                 # Output columns
                 D_colname = "D", deltaV_colname = "dV", F_colname = "F"){

  pad <- match.arg(pad)

  # Establish names for some intermediate columns.
  v_colname <- ".v"
  V_colname <- ".V"
  L_name <- "L"
  w_name <- "w"

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
  Lv_colname <- paste0(L_name, "(", v_colname, ")")
  wv_colname <- paste0(w_name, "(", v_colname, ")")
  deltaV_vec_colname <- paste0(deltaV_colname, vec_suffix)
  D_vec_colname <- paste0(D_colname, vec_suffix)

  # Ensure that time_colname is NOT a grouping variable.
  if (time_colname %in% groups(.lmdidata)) {
    stop(paste0("'", time_colname, "'", " is a grouping variable, but you can't group on ",
               "time_colname",
               " in argument .DF of collapse_to_matrices."))
  }

  XvV <- .lmdidata %>% mutate(
    !!as.name(v_colname) := rowprods_byname(!!as.name(X_colname)),
    # Add as.numeric() here to get a single number, not a 1x1 matrix.
    !!as.name(V_colname) := colsums_byname(!!as.name(v_colname)) %>% as.numeric()
  )

  XvV0T <- create0Tcolumns(XvV, time_colname = time_colname,
                             X_colname = X_colname, v_colname = v_colname, V_colname = V_colname,
                             pad = pad, zero_suffix = zero_suffix, T_suffix = T_suffix)
  # Do year-by-year LMDI calcs.
  XvV0T %>%
    mutate(
      !!as.name(D_colname) := elementquotient_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      !!as.name(deltaV_colname) := difference_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      !!as.name(LV_colname) := logarithmicmean_byname(!!as.name(VT_colname), !!as.name(V0_colname)),
      !!as.name(Lv_colname) := logarithmicmean_byname(!!as.name(vT_colname), !!as.name(v0_colname)),
      !!as.name(wv_colname) := elementquotient_byname(!!as.name(Lv_colname), !!as.name(LV_colname)),
      !!as.name(F_colname) := elementquotient_byname(X_T, X_0) %>% elementlog_byname(),
      !!as.name(deltaV_vec_colname) := matrixproduct_byname(hatize_byname(!!as.name(Lv_colname)), !!as.name(F_colname)) %>%
        colsums_byname() %>% transpose_byname(),
      !!as.name(D_vec_colname) := matrixproduct_byname(hatize_byname(!!as.name(wv_colname)), !!as.name(F_colname)) %>%
        colsums_byname() %>% elementexp_byname() %>% transpose_byname()
    )
}



