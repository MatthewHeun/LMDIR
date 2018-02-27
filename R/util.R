#' Z matrix from \code{X_0} and \code{X_T} matrices
#'
#' The formula for Z is \deqn{Z_ij = logmean(v_Ti1, v_0i1) * log(X_Tij / X_0ij)}
#' where \code{v} is a column vector formed by row products of \code{X}.
#' The \code{T} and \code{0} subscripts on \code{X} and \code{v} indicate an initial time (\code{0})
#' and a final time (\code{T}).
#'
#' The nomenclature for this function comes from
#' \href{https://doi.org/10.1016/s0360-5442(98)00016-4}{Ang, Zhang, and Choi (1998)}.
#' This function fully accounts for the degenerate cases
#' found in Table 2, p. 492 of
#' \href{https://doi.org/10.1016/s0360-5442(98)00016-4}{Ang, Zhang, and Choi (1998)}.
#' We're employing the method suggested by
#' \href{https://doi.org/10.1016/j.enpol.2004.11.010}{Wood and Lenzen (2006)}.
#'
#' @param X_0 an \code{X} matrix for initial time \code{0}
#' @param X_T an \code{X} matrix for final time \code{T}
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
Z_byname <- function(X_0, X_T){
  Z.func <- function(X_0, X_T){
    # Ensure that X_0 and X_T are same type of matrices.
    # I.e., they have the same row and column names.
    # And they have the same row and column types.
    stopifnot(samestructure_byname(X_0, X_T))
    # Create an empty z matrix.  This will be filled with default value (NA) entires.
    Z <- matrix(nrow = nrow(X_0), ncol = ncol(X_0)) %>%
      setrownames_byname(rownames(X_0)) %>% setcolnames_byname(colnames(X_0))
    # Use an old-fashioned for loop to fill all elements ofthe Z matrix
    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        Z[[i, j]] <- Zij(i = i, j = j, X_0 = X_0, X_T = X_T)
      }
    }
    return(Z)
  }
  # Take control of completing and sorting matrices here, because
  # we have a more-complex situation than simply filling the missing rows with 0s.
  # temp <- complete_and_sort(X_0, X_T, rowfill = c(E.ktoe = 1, eta_ij = 1, phi_i = 1, phi_ij = 0))
  # X_0_completed <- temp$m1
  # X_T_completed <- temp$m2
  # binaryapply_byname(Z.func, a = X_0_completed, b = X_T_completed, match_type = "all")
  binaryapply_byname(Z.func, a = X_0, b = X_T, match_type = "all")
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

