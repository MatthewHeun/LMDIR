#' Log-Mean Divisia Index (LMDI) decomposition analysis
#'
#' @param .data a grouped data frame
#' @param time_colname the name of the column in \code{.data} that contains times at which
#'        data are available (a string).
#'        Default is "\code{Year}".
#' @param X the name (as a string) of a column in \code{.data} containing
#'        \strong{\code{X}} matrices with
#'        named rows representing subcategories of the energy aggregate (\code{V}) and
#'        named columns representing factors contributing to changes in \code{V} over time.
#'
#' @return a data frame containing several columns.
#' @export
#'
#' @examples
lmdi <- function(.data, time_colname = "Year", X_colname = "X",
                 # Output columns
                 v_colname = "v",
                 V_colname = "V"){
  .data %>% mutate(
    !!as.name(v_colname) := rowprods_byname(!!as.name(X_colname)),
    !!as.name(V_colname) := colsums_byname(!!as.name(v_colname))
  )
}



#' #' Log-Mean Divisia Index (LMDI) decomposition analysis
#' #'
#' #' @param .data a grouped data frame
#' #' @param time_colname the name of the column in \code{.data} that contains times at which
#' #'        data are available (a string).
#' #'        Default is "\code{Year}".
#' #' @param agg_structure a named list of vectors of column names (strings) in \code{.data}
#' #'        that gives the aggregation structure of \code{.data}.
#' #'        The names of \code{agg_structure} must be names of columns in \code{.data},
#' #'        and they represent the names of sub-categories \code{V_i} of the aggregate \code{V}.
#' #'        Each item in \code{agg_structure} must be a list or vector of strings.
#' #'        Each item in \code{agg_structure} is interpreted as names of columns in \code{.data}
#' #'        that comprise sub-aggregates.
#' #'        The columns \code{.data} named in the items of \code{agg_structure}
#' #'        must be names for each \code{x} in the LMDI decomposition analysis.
#' #'
#' #' @return a data frame containing several columns.
#' #' @export
#' #'
#' #' @examples
#' lmdi <- function(.data, time_colname = "Year", agg_structure,
#'                  # Output columns
#'                  V_colname = "V",
#'                  L_colname = "L"){
#'   groups_orig <- group_vars(.data)
#'   if (length(groups_orig > 0)) {
#'     # Ensure that .data is not grouped by the Vx_colnames.
#'     if (grepl(paste(agg_structure, collapse = "|"), groups_orig)) {
#'       # This is an error
#'       stop("One of the agg_structure columns was also a grouping variable.")
#'     }
#'   }
#'   # Create sub-aggregates
#'   prodsumcol <- ".prodsumcol"
#'   V_i <- DF %>%
#'     # Be sure we are grouped by year (at least). .data may come in with other groups, too.
#'     group_by(!!as.name(time_colname), add = TRUE) %>%
#'     gather(key = "var", value = "val", -!!!as.name(c(time_colname, groups_orig))) %>%
#'     # Add a column of variable names which are the V_i's
#'     # associated with each x_ki.
#'     mutate(
#'       agg_var = lapply(1:length(agg_structure), FUN = function(i){
#'         rep.int(names(agg_structure)[[i]], length(agg_structure[[i]]))
#'       }) %>% unlist
#'     ) %>%
#'     group_by(agg_var, add = TRUE) %>%
#'     summarise(!!as.name(prodsumcol) := prod(val))
#'   # Create aggregate
#'   V <- V_i %>%
#'     group_by(!!!c(as.name(time_colname), groups_orig)) %>%
#'     summarise(!!as.name(prodsumcol) := sum(!!as.name(prodsumcol))) %>%
#'     rename(!!as.name(V_colname) := !!as.name(prodsumcol))
#'
#'   Aggregates <- full_join(V, V_i %>% spread(key = agg_var, value = !!as.name(prodsumcol)),
#'     by = c(time_colname, groups_orig)) %>%
#'     full_join(.data, by = c(time_colname, groups_orig)) %>%
#'     ungroup()
#'
#'   out <- Aggregates %>%
#'     # Calculate all the "L" values
#'     do(Louterfun(., c(V_colname, names(agg_structure)))) %>%
#'     # Calculate all w_i values where w_i = L(V_i)/L(V) values
#'     do(wouterfun(., names(agg_structure), V_colname)) %>%
#'     do(fouterfun(., unlist(agg_structure)))
#' }
#'
#' wfouterfun <- function(.data, x_names)
#'
#' wfinnerfun <- function(.data, x_name){
#'   w_col <- as.name(paste0("w(", x_name, ")"))
#'   f_col <- as.name(paste0("f(", x_name, ")"))
#'   .data <- mutate(
#'
#'   )
#' }
#'
#' fouterfun <- function(.data, x_names){
#'   lapply(x_names, FUN = function(xn){
#'     finnerfun(.data, xn)
#'   }) %>%
#'     cbind(.data, .)
#' }
#'
#' finnerfun <- function(.data, x_name, log.string = "ln"){
#'   x_col <- as.name(x_name)
#'   lnx_col <- as.name(paste0(log.string, "(", x_name, ")"))
#'   f_col <- as.name(paste0("f(", x_name, ")"))
#'   .data %>%
#'     mutate(
#'       # f is defined as ln(xT/x0).
#'       # But that's same as ln(xT) - ln(x0).
#'       # Choose second approach, as that is compatible with ediff
#'       !!lnx_col := log(!!x_col),
#'       !!f_col := ediff(!!lnx_col, pad = "tail")
#'     ) %>%
#'     select(!!f_col)
#' }
#'
#' wouterfun <- function(.data, Vi_names, V_colname){
#'   lapply(Vi_names, FUN = function(cn){
#'     winnerfun(.data, Vi_colname = cn, V_colname = V_colname)
#'   }) %>%
#'     cbind(.data, .)
#' }
#'
#' winnerfun <- function(.data, Vi_colname, V_colname){
#'   w_col <- as.name(paste0("w(", Vi_colname, ")"))
#'   LVi_col <- as.name(paste0("L(", Vi_colname, ")"))
#'   LV_col <- as.name(paste0("L(", V_colname, ")"))
#'   .data %>%
#'     mutate(
#'       !!w_col := (!!LVi_col) / !!LV_col
#'     ) %>%
#'     select(!!w_col)
#' }
#'
#'
#' Louterfun <- function(.data, V_names){
#'   lapply(V_names, FUN = function(cn){
#'     Linnerfun(.data, Vi_colname = cn)
#'   }) %>%
#'     cbind(.data, .)
#' }
#'
#' Linnerfun <- function(.data, Vi_colname, delta.string = "∆", log.string = "ln"){
#'   Vi_col <- as.name(Vi_colname)
#'   L_col <- as.name(paste0("L(", Vi_colname, ")"))
#'   deltaVi_col <- as.name(paste0(delta.string, Vi_colname))
#'   logVi_col <- as.name(paste0(ln.string, "(", Vi_colname, ")"))
#'   deltalogVi_col <- as.name(paste0(delta.string, logVi_col))
#'
#'   .data %>%
#'     mutate(
#'       !!logVi_col := log(!!Vi_col),
#'       !!deltaVi_col := ediff(!!Vi_col, pad = "tail"),
#'       !!deltalogVi_col := ediff(!!logVi_col, pad = "tail"),
#'       !!L_col := (!!deltaVi_col) / !!deltalogVi_col
#'     ) %>%
#'     select(!!L_col)
#' }
