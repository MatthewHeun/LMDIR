#' Z matrix from \code{X_0} and \code{X_T} matrices
#'
#' The formula for \code{Z} is
#' \deqn{Z_ij = logmean(v_Ti1, v_0i1) * log(X_Tij / X_0ij)}
#' where \code{v} is a column vector formed by row products of \code{X}.
#' The \code{0} and \code{T} subscripts on \code{X} and \code{v} indicate an initial time (\code{0})
#' and a final time (\code{T}).
#' The \code{i} and \code{j} subscripts on \code{Z}, \code{v}, and \code{X} are
#' matrix indices.
#'
#' \code{Z} and \code{X} are category-by-factor matrices, and
#' \code{v} is a category-by-1 column vector.
#'
#' When a category comes and goes through time,
#' the category row can be absent from the \code{X} matrix at one time
#' but present in the \code{X} matrix at an adjacent time.
#' If the missing category row is replaced by \code{0}s,
#' the LMDI algorithm fails due to \code{log(0)} errors.
#' The usual advice is to insert a row consisting not of \code{0}s
#' but rather filled with small numbers (e.g., \code{1e-100}).
#'
#' But the usual advice doesn't represent reality for some LMDI decomposition analyses.
#' To provide greater flexibility, this function provides the \code{emptyrow} argument.
#' Callers can supply their own \code{emptyrow},
#' a single-row \code{matrix} with column names identical
#' to \code{X_0} and \code{X_T}.
#' If \code{emptyrow} is not specified,
#' the usual advice will be followed, and
#' a row vector consisting of very small values (\code{1e-100}) will be inserted into
#' \code{X_0} or \code{X_T}, as appropriate.
#'
#' Note that the \code{lmdi} function passes its \code{emptyrow} argument, if present,
#' to \code{Z_byname}.
#'
#' The nomenclature for this function comes from
#' \href{https://doi.org/10.1016/s0360-5442(98)00016-4}{Ang, Zhang, and Choi (1998)}.
#' This function fully accounts for the degenerate cases
#' found in Table 2, p. 492 of
#' \href{https://doi.org/10.1016/s0360-5442(98)00016-4}{Ang, Zhang, and Choi (1998)}.
#' When \code{0}s are encountered in \code{emptyrow},
#' this function employs the method for dealing with the \code{log(0)} problem suggested by
#' \href{https://doi.org/10.1016/j.enpol.2004.11.010}{Wood and Lenzen (2006)}.
#'
#' @param X_0 an \code{X} matrix for initial time \code{0}
#' @param X_T an \code{X} matrix for final time \code{T}
#' @param fillrow a row vector of type \code{matrix} with column names identical
#'                to \code{X_0} and \code{X_T} (see details)
#'
#' @return a \code{Z} matrix
#'
#' @importFrom matsbyname samestructure_byname
#' @importFrom matsbyname setrownames_byname
#' @importFrom matsbyname setcolnames_byname
#' @importFrom matsbyname binaryapply_byname
#' @importFrom matsbyname logmean
#' @importFrom matsbyname sumall_byname
#' @importFrom matsbyname transpose_byname
#' @importFrom matsbyname elementexp_byname
#' @importFrom matsbyname prodall_byname
#' @importFrom matsbyname cumsum_byname
#' @importFrom matsbyname cumprod_byname
#' @importFrom dplyr everything
#' @importFrom dplyr group_by
#'
#' @export
Z_byname <- function(X_0, X_T, fillrow = NULL){
  Z.func <- function(X_0, X_T, fillrow = NULL){
    # At this point, X_0 and X_T are single matrices.
    # We need to take control of completing and sorting X_0 and X_T matrices here, because
    # we have a more-complex situation than simply filling the missing rows with 0s.
    if (is.null(fillrow)) {
      emptyrow <- matrix(1e-100, nrow = 1, ncol = ncol(X_0),
                         dimnames = list("row", colnames(X_0)))
    }
    X_0_comp <- complete_rows_cols(X_0, X_T, fillrow = fillrow, margin = 1)
    X_T_comp <- complete_rows_cols(X_T, X_0, fillrow = fillrow, margin = 1)
    X_0_comp_sort <- sort_rows_cols(X_0_comp)
    X_T_comp_sort <- sort_rows_cols(X_T_comp)
    # At this point, X_0_comp_sort and X_T_comp_sort should be the same type of matrices.
    # I.e., they have the same row and column names.
    # And they have the same row and column types.
    # Ensure that this is so!
    stopifnot(samestructure_byname(X_0_comp_sort, X_T_comp_sort))

    # Create an empty Z matrix.  Z will be filled with default entries (NA).
    Z <- matrix(nrow = nrow(X_0_comp_sort), ncol = ncol(X_0_comp_sort)) %>%
      setrownames_byname(rownames(X_0_comp_sort)) %>% setcolnames_byname(colnames(X_0_comp_sort)) %>%
      setrowtype(rowtype(X_0_comp_sort)) %>% setcoltype(coltype(X_0_comp_sort))
    # Use an old-fashioned for loop to fill all elements ofthe Z matrix
    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        Z[[i, j]] <- Zij(i = i, j = j, X_0 = X_0_comp_sort, X_T = X_T_comp_sort)
      }
    }
    return(Z)
  }

  # If we are missing a row in X_0 or X_T compared to the other,
  # it is because that particular type subsubcategory of useful exergy is present in one year
  # but absent in the other.
  # In this situation, we want to fill values in the missing row with non-zero numbers (42) for
  # * primary exergy (E.ktoe),
  # * allocation from primary exergy to subcategory (phi_i), and
  # * primary-to-useful efficiency (eta_ij).
  # Then, we set the allocation from subcategory to subsubcategory (phi_ij) to 0.
  # This approach correctly models the fact that
  # despite the fact that we have no useful exergy of this type being produced in
  # one of the years, we still have
  # total primary en/xergy (E.ktoe),
  # there is still an allocation of primary exergy to the subcategory (phi_i), and
  # if there were machines making this subsubcategory of useful exergy
  # in this time period, it would have a certain primary-to-useful efficiency (eta_ij).
  # It turns out that we don't need to know the exact values of
  # primary exergy (E.ktoe),
  # allocation to subcategory (phi_i), or
  # primary-to-useful efficiency (eta_ij).
  # These values must simply be non-zero so long as
  # allocation from subcategory to subsubcategory (phi_ij) is zero.
  #
  #
  # fr <- matrix(c(42, 42, 42, 0), nrow = 1, ncol = 4,
  #              dimnames = list("row", c("E.ktoe", "eta_ij", "phi_i", "phi_ij"))) %>%
  #   setrowtype("category") %>% setcoltype("factor")

  binaryapply_byname(Z.func, a = X_0, b = X_T,
                     .FUNdots = list(fillrow = fillrow), match_type = "all", .organize = FALSE)
  # binaryapply_byname(Z.func, a = X_0, b = X_T, match_type = "all")
}

#' Calculate element \code{Z_i,j}
#'
#' There are many special cases for calculating the \code{i,j}th term in \code{Z}.
#' These special cases must be handled on a term-by term basis.
#' This function handles all those cases.
#' See \code{\link{Z_byname}} for details.
#' The \strong{\code{X}} matrices must be same size, must have same row and column names, and
#' must have same row and column types.
#' Note that \code{i} and \code{j} must be positive and less than
#' \code{nrow(X)} and \code{ncol(X)}, respectively.
#'
#' @param i row index
#' @param j column index
#' @param X_0 the sub-sector by factor matrix for time 0
#' @param X_T the sub-sector by factor matrix for time T
#' @param v_0i1 the i,1th element of the \code{v_0} column vector. (v_0 is formed from the row products of the \code{X_0} matrix.)
#' @param v_Ti1 the i,1th element of the \code{v_T} column vector. (v_T is formed from the row products of the \code{X_T} matrix.)
#' @param X_0ij the i,jth element of the \code{X_0} matrix
#' @param X_Tij the i,jth element of the \code{X_T} matrix
#'
#' @return the Z value corresponding to the \code{v_0i1}, \code{v_Ti1}, \code{X_0ij}, and \code{X_Tij} values
#'
#' @export
Zij <- function(i, j, X_0, X_T,
                v_0i1 = rowprods_byname(X_0)[i, 1],
                v_Ti1 = rowprods_byname(X_T)[i, 1],
                X_0ij = X_0[i, j],
                X_Tij = X_T[i, j]){

  # Check the conditions, found in Table 2, p. 492 of
  # B.W. Ang and F.Q. Zhang and Ki-Hong Choi, 1998,
  # Factorizing changes in energy and environmental indicators through decomposition,
  # Energy, Volume 23, Number 6, pp. 489-495.
  if (v_0i1 == 0 & v_Ti1 > 0 & X_0ij == 0 & X_Tij > 0) {
    # Case 1
    return(v_Ti1)

  } else if (v_0i1 > 0 & v_Ti1 == 0 & X_0ij > 0 & X_Tij == 0) {
    # Case 2
    return(-v_0i1)

  } else if (v_0i1 == 0 & v_Ti1 > 0 & X_0ij > 0 & X_Tij > 0) {
    # Case 3
    return(0)

  } else if (v_0i1 > 0 & v_Ti1 == 0 & X_0ij > 0 & X_Tij > 0) {
    # Case 4
    return(0)

  } else if (v_0i1 == 0 & v_Ti1 == 0 & X_0ij > 0 & X_Tij > 0) {
    # Case 5
    return(0)

  } else if (v_0i1 == 0 & v_Ti1 == 0 & X_0ij == 0 & X_Tij == 0) {
    # Case 6
    return(0)

  } else if (v_0i1 == 0 & v_Ti1 == 0 & X_0ij > 0 & X_Tij == 0) {
    # Case 7
    return(0)

  } else if (v_0i1 == 0 & v_Ti1 == 0 & X_0ij == 0 & X_Tij > 0) {
    # Case 8
    return(0)

  } else if (v_0i1 > 0 & v_Ti1 > 0 & X_0ij > 0 & X_Tij > 0) {
    # This is the non-degenerate case
    return(logmean(v_Ti1, v_0i1) * log(X_Tij / X_0ij))
  }

  stop("Unknown conditions for v_0i1, v_Ti1, X_0ij, and X_Tij in Zij")
}

#' Make columns for time 0 and time T
#'
#' @param XvV a data frame containing \code{X}, \code{v}, and \code{V} columns.
#' @param time_colname the name of the time column (a string)
#' @param X_colname the name of the \code{X} column (default is "\code{X}")
#' @param v_colname the name of the \code{v} column (default is "\code{v}")
#' @param V_colname the name of the \code{V} column (default is "\code{V}")
#' @param zero_suffix suffix for "\code{0}" variables (default is "\code{_0}")
#' @param T_suffix suffix for "\code{T}" variables (default is "\code{_T}")
#'
#' @importFrom dplyr do
#' @importFrom dplyr group_vars
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @return a data frame containing metadata (\code{time_colname} and grouping variables)
#'         and columns for X_0, v_0, V_0, X_T, v_T, and V_T.
create0Tcolumns <- function(XvV,
                            time_colname,
                            X_colname = "X", v_colname = "v", V_colname = "V",
                            zero_suffix = "_0",
                            T_suffix = "_T"){

  # Establish names for new columns.
  X0_colname <- paste0(X_colname, zero_suffix)
  v0_colname <- paste0(v_colname, zero_suffix)
  V0_colname <- paste0(V_colname, zero_suffix)
  XT_colname <- paste0(X_colname, T_suffix)
  vT_colname <- paste0(v_colname, T_suffix)
  VT_colname <- paste0(V_colname, T_suffix)

  # Add a repeat of the first row at the top.
  # Doing so ensures that
  # the calculation of the first row deltaV values gives 0 and
  # the calculation of the first row D values gives 1, as it should.
  # Performing this action with "do" ensures that each group has a repeated first row.
  XvV_repeat1strow <- XvV %>%
    do(
      rbind(.data[1, ], .data)
    )

  # Set up for aligning the rows for further calculations.
  .DF0 <- XvV_repeat1strow %>%
    do(
      # do works in groups, which is what we want.
      head(.data, -1)
    ) %>%
    rename(
      !!as.name(X0_colname) := !!as.name(X_colname),
      !!as.name(v0_colname) := !!as.name(v_colname),
      !!as.name(V0_colname) := !!as.name(V_colname)
    ) %>%
    # Don't need grouping variables or time variable here.
    # We'll pick them up from .DFT.
    ungroup() %>%
    select(!!as.name(X0_colname), !!as.name(v0_colname), !!as.name(V0_colname))

  .DFT <- XvV_repeat1strow %>%
    do(
      # do works in groups, which is what we want.
      tail(.data, -1)
    ) %>%
    rename(
      !!as.name(XT_colname) := !!as.name(X_colname),
      !!as.name(vT_colname) := !!as.name(v_colname),
      !!as.name(VT_colname) := !!as.name(V_colname))

  # Bind everything together and return it
  cbind(.DF0 %>% ungroup(), .DFT %>% ungroup()) %>%
    select(group_vars(XvV), time_colname, everything()) %>%
    group_by(!!!groups(XvV))
}

