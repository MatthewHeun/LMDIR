#' Z matrix from `X_0` and `X_T` matrices
#'
#' The formula for `Z` is
#' \deqn{Z_ij = logmean(v_Ti1, v_0i1) * log(X_Tij / X_0ij)}
#' where `v` is a column vector formed by row products of `X`.
#' The `0` and `T` subscripts on `X` and `v` indicate an initial time (`0`)
#' and a final time (`T`).
#' The `i` and `j` subscripts on `Z`, `v`, and `X` are
#' matrix indices.
#'
#' `Z` and `X` are category-by-factor matrices, and
#' `v` is a category-by-1 column vector.
#'
#' When a category comes and goes through time,
#' the category row can be absent from the `X` matrix at one time
#' but present in the `X` matrix at an adjacent time.
#' If the missing category row is replaced by `0`s,
#' the LMDI algorithm fails due to `log(0)` errors.
#' The usual advice is to insert a row consisting not of `0`s
#' but rather filled with small numbers (e.g., `1e-10`).
#'
#' But the usual advice doesn't represent reality for some LMDI decomposition analyses.
#' For example, the missing row may be caused by only one of the many factors being zero,
#' the other factors being non-zero.
#' To provide greater flexibility, this function provides the `fillrow` argument.
#' Callers can supply their own `fillrow`,
#' a single-row `matrix` with column names identical
#' to `X_0` and `X_T`.
#' If `fillrow` is not specified,
#' the usual advice will be followed, and
#' a row vector consisting of very small values (`1e-10`)
#' for each factor will be inserted into
#' `X_0` or `X_T`, as appropriate.
#'
#' Note that the [lmdi()] function passes its `fillrow` argument, if present,
#' to `Z_byname`.
#'
#' The nomenclature for this function comes from
#' [Ang, Zhang, and Choi (1998)](\doi{10.1016/s0360-5442(98)00016-4}).
#' This function fully accounts for the degenerate cases
#' found in Table 2, p. 492 of
#' [Ang, Zhang, and Choi (1998)](\doi{10.1016/s0360-5442(98)00016-4}).
#' When `0`s are encountered in `emptyrow`,
#' this function employs the method for dealing with the `log(0)` problem suggested by
#' [Wood and Lenzen (2006)](\doi{10.1016/j.enpol.2004.11.010}).
#'
#' @param X_0 an `X` matrix for initial time `0`
#' @param X_T an `X` matrix for final time `T`
#' @param fillrow a row vector of type `matrix` with column names identical
#'                to `X_0` and `X_T`. (See details.)
#' @param weights the name of the weights type chosen for the LMDI calculations.
#'        `weights` should be one of "LMDI-I" or "LMDI-II". Default is "LMDI-I".
#'
#' @return A `Z` list of matrices for LMDI-I and LMDI-II weights.
#'
#' @export
#'
Z_byname <- function(X_0, X_T, fillrow = NULL, weights = c("LMDI-I", "LMDI-II")){
  fillrow <- matsbyname::prep_vector_arg(X_0, vector_arg = fillrow)
  weights <- match.arg(weights, choices = c("LMDI-I", "LMDI-II"))
  Z_func <- function(X_0, X_T, fillrow = NULL){
    # At this point, X_0 and X_T are single matrices.
    # We need to take control of completing and sorting X_0 and X_T matrices here, because
    # we have a more-complex situation than simply filling the missing rows with 0s.
    if (is.null(fillrow)) {
      # If the caller didn't supply a value for fillrow, we default to the
      # "small value" approximation suggested by
      # B.W. Ang and F.Q. Zhang and Ki-Hong Choi, 1998,
      # Factorizing changes in energy and environmental indicators through decomposition,
      # Energy, Volume 23, Number 6, pp. 489-495.
      fillrow <- matrix(1e-10, nrow = 1, ncol = ncol(X_0),
                        dimnames = list("row", colnames(X_0))) %>%
        matsbyname::setrowtype(matsbyname::rowtype(X_0)) %>% matsbyname::setcoltype(matsbyname::coltype(X_0))
    }
    # Complete the matrices relative to one another, using fillrow.
    X_0_comp <- matsbyname::complete_rows_cols(X_0, X_T, fillrow = fillrow, margin = 1)
    X_T_comp <- matsbyname::complete_rows_cols(X_T, X_0, fillrow = fillrow, margin = 1)
    # Sort the matrices relative to one another.
    X_0_comp_sort <- matsbyname::sort_rows_cols(X_0_comp)
    X_T_comp_sort <- matsbyname::sort_rows_cols(X_T_comp)
    # At this point, X_0_comp_sort and X_T_comp_sort should be the same type of matrices.
    # I.e., they have the same row and column names.
    # And they have the same row and column types.
    # Ensure that this is so!
    stopifnot(matsbyname::samestructure_byname(X_0_comp_sort, X_T_comp_sort))

    # Create an empty Z matrix.
    # The empty Z is filled with default entries (NA).
    Z <- matrix(nrow = nrow(X_0_comp_sort), ncol = ncol(X_0_comp_sort)) %>%
      matsbyname::setrownames_byname(rownames(X_0_comp_sort)) %>% matsbyname::setcolnames_byname(colnames(X_0_comp_sort)) %>%
      matsbyname::setrowtype(matsbyname::rowtype(X_0_comp_sort)) %>% matsbyname::setcoltype(matsbyname::coltype(X_0_comp_sort))

    # Use an old-fashioned for loop to fill all elements of the Z matrix
    for (i in 1:nrow(Z)) {
      for (j in 1:ncol(Z)) {
        Z[i, j] <- Zij(i = i, j = j, X_0 = X_0_comp_sort, X_T = X_T_comp_sort, weights = weights)
      }
    }
    return(Z)
  }

  matsbyname::binaryapply_byname(Z_func, a = X_0, b = X_T,
                                 .FUNdots = list(fillrow = fillrow), match_type = "all", .organize = FALSE)
}

#' Calculate element `Z_i,j`
#'
#' There are many special cases for calculating the `i,j`th term in `Z`.
#' These special cases must be handled on a term-by term basis.
#' This function handles all those cases.
#' See [Z_byname()] for details.
#' The `X` matrices must be same size, must have same row and column names, and
#' must have same row and column types.
#' Note that `i` and `j` must be positive and less than or equal to
#' `nrow(X)` and `ncol(X)`, respectively.
#'
#' Arguments `i`, `j`, `X_0`, and `X_T` are optional.
#' If they are not specified,
#' arguments `v_0`, `v_T`, `V_0`, `V_T`, `v_0i1`, `v_Ti1`, `X_0ij`, and `X_Tij`
#' must be given.
#'
#' @param i optional row index for `X_0` and `X_T`.
#' @param j optional column index for `X_0` and `X_T`.
#' @param X_0 optional sub-sector by factor matrix for time 0.
#' @param X_T optional sub-sector by factor matrix for time T.
#' @param weights the name of the weights type chosen for the LMDI calculations.
#'        `weights` should be one of "LMDI-I" or "LMDI-II". Default is "LMDI-I".
#' @param v_0 the column vector formed from the row products of the `X_0` matrix (only necessary for LMDI-II weights).
#' @param v_T the column vector formed from the row products of the `X_T` matrix (only necessary for LMDI-II weights).
#' @param V_0 the column sum of the `v_0` column vector (only necessary for LMDI-II weights).
#' @param V_T the column sum of the `v_T` column vector (only necessary for LMDI-II weights).
#' @param v_0i1 the i,1th element of the `v_0` column vector.
#' @param v_Ti1 the i,1th element of the `v_T` column vector.
#' @param X_0ij the i,jth element of the `X_0` matrix
#' @param X_Tij the i,jth element of the `X_T` matrix
#'
#' @return the `Z_I_ij` the LMDI-I value corresponding to the `v_0i1`, `v_Ti1`, `X_0ij`, and `X_Tij` values and the
#' `Z_II_ij` the LMDI-II Z value corresponding to the `v_0`, `v_T`, `V_0`, `V_T`, `v_0i1`, `v_Ti1`, `X_0ij`, and
#' `V_Tij` values for row `i` and column `j`.
#'
#' @export
Zij <- function(i = NULL, j = NULL, X_0 = NULL, X_T = NULL,
                weights = c("LMDI-I", "LMDI-II"),
                v_0 = matsbyname::rowprods_byname(X_0)[, 1],
                v_T = matsbyname::rowprods_byname(X_T)[, 1],
                V_0 = matsbyname::colsums_byname(matsbyname::rowprods_byname(X_0))[, 1],
                V_T = matsbyname::colsums_byname(matsbyname::rowprods_byname(X_T))[, 1],
                v_0i1 = matsbyname::rowprods_byname(X_0)[i, 1],
                v_Ti1 = matsbyname::rowprods_byname(X_T)[i, 1],
                X_0ij = X_0[i, j],
                X_Tij = X_T[i, j]){

  weights <- match.arg(weights, choices = c("LMDI-I", "LMDI-II"))

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
    if (weights == "LMDI-I") {
      return(matsbyname::logmean(v_Ti1, v_0i1) * log(X_Tij / X_0ij))
    } else {
      return(
        matsbyname::logmean((v_Ti1/V_T), (v_0i1/V_0)) *
          matsbyname::logmean(V_T, V_0) /
          sum(
            matsbyname::logarithmicmean_byname((v_T/as.numeric(V_T)), (v_0/as.numeric(V_0)))) *
          log(X_Tij / X_0ij)
      )}
  }
  # We should never get here.
  stop("Unknown conditions for v_0i1, v_Ti1, X_0ij, and X_Tij in Zij")
}

#' Make columns for time 0 and time T
#'
#' @param XvV a data frame containing `X`, `v`, and `V` columns.
#' @param time_colname the name of the time column (a string)
#' @param X_colname the name of the `X` column (default is "`X`")
#' @param v_colname the name of the `v` column (default is "`v`")
#' @param V_colname the name of the `V` column (default is "`V`")
#' @param zero_suffix suffix for "`0`" variables (default is "`_0`")
#' @param T_suffix suffix for "`T`" variables (default is "`_T`")
#'
#' @importFrom dplyr do
#' @importFrom dplyr group_vars
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom rlang .data
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @return a data frame containing metadata (`time_colname` and grouping variables)
#'         and columns for X_0, v_0, V_0, X_T, v_T, and V_T.
create0Tcolumns <- function(XvV,
                            time_colname,
                            X_colname = "X", v_colname = "v", V_colname = "V",
                            zero_suffix = "_0",
                            T_suffix = "_T"){

  # Eliminate the NOTE in R CMD chk about no visible global binding for "."
  . <- NULL

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
  # . refers to the current group.
  # (See https://dplyr.tidyverse.org/reference/do.html for details.)
  XvV_repeat1strow <- XvV %>%
    dplyr::do(
      rbind(head(., n = 1), .)
    )

  # Set up for aligning the rows for further calculations.
  .DF0 <- XvV_repeat1strow %>%
    dplyr::do(
      # do works in groups, which is what we want.
      # . refers to the current group.
      head(., -1)
    ) %>%
    rename(
      !!as.name(X0_colname) := !!as.name(X_colname),
      !!as.name(v0_colname) := !!as.name(v_colname),
      !!as.name(V0_colname) := !!as.name(V_colname)
    ) %>%
    # Don't need grouping variables or time variable here.
    # We'll pick them up from .DFT.
    dplyr::ungroup() %>%
    dplyr::select(!!as.name(X0_colname), !!as.name(v0_colname), !!as.name(V0_colname))

  .DFT <- XvV_repeat1strow %>%
    dplyr::do(
      # do works in groups, which is what we want.
      # . refers to the current group.
      tail(., -1)
    ) %>%
    dplyr::rename(
      !!as.name(XT_colname) := !!as.name(X_colname),
      !!as.name(vT_colname) := !!as.name(v_colname),
      !!as.name(VT_colname) := !!as.name(V_colname))

  # Bind everything together and return it
  cbind(.DF0 %>% dplyr::ungroup(), .DFT %>% dplyr::ungroup()) %>%
    dplyr::select(dplyr::group_vars(XvV), dplyr::all_of(time_colname), dplyr::everything()) %>%
    dplyr::group_by(!!!dplyr::groups(XvV))
}


#' Create a simple LMDI data frame
#'
#' This function creates a simple LMDI data frame that can be used for testing or examples.
#' It contains `X` matrices for two fictitious countries (`AB` and `XY`)
#' and four years (`1971-1974`).
#' Countries `AB` and `YZ` are identical for the purposes of illustration.
#'
#' @return a data frame of example LMDI data
#'
#' @importFrom matsindf collapse_to_matrices
#'
#' @export
#'
#' @examples
#' library(matsindf)
#' DF <- create_simple_LMDI()
#' DF
#' DF$X[[1]]
create_simple_LMDI <- function(){
  simple.factor.names <- c("factor 1", "factor 2", "factor 3")
  simple.subsubcat.names <- c("subsubcat 1", "subsubcat 2")
  # Eliminate the NOTE in R CMD chk about no visible global binding for these variables
  Country <- NULL
  Year <- NULL
  x <- NULL
  # Create a tidy data frame of x values
  # which will be collapsed to matrices
  AB <- data.frame(Year = rep.int(1971, 6), x = c(1, 10, 2, 4, 5, 3)) %>%
    rbind(data.frame(Year = rep.int(1972, 6), x = c(4, 5, 3, 5, 6, 4))) %>%
    rbind(data.frame(Year = rep.int(1973, 6), x = c(8, 2, 4, 5, 7, 5))) %>%
    rbind(data.frame(Year = rep.int(1974, 6), x = c(10, 1, 5, 6, 8, 6))) %>%
    dplyr::mutate(
      matnames = "X",
      rowtypes = "subsubcat",
      coltypes = "factor",
      rownames = rep.int(c(rep.int(simple.subsubcat.names[[1]], 3), rep.int(simple.subsubcat.names[[2]], 3)), 4),
      colnames = rep.int(c(rep.int(simple.factor.names, 2)), 4),
      Country = "AB"
    ) %>%
    dplyr::group_by(Country, Year) %>%
    matsindf::collapse_to_matrices(matnames = "matnames", matvals = "x",
                         rownames = "rownames", colnames = "colnames",
                         rowtypes = "rowtypes", coltypes = "coltypes") %>%
    dplyr::rename(X = x)
  rbind(AB, AB %>% dplyr::mutate(Country = "YZ"))
}

