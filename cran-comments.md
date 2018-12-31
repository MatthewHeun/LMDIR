## Context
`LMDIR` is a new package that performs log-mean divisia index analysis on a data frame.

## Test environments
* local macOS X install 10.14.2 (Mojave), R3.5.2
* ubuntu 14.04.5 (on Travis CI), R3.5.1
* windows (on win-builder)
    * `devtools::check_win_oldrelease()`, R3.4.4
    * `devtools::check_win_release()`, R3.5.2
    * `devtools::check_win_devel()`, r75909
* rhub
    * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
    * Ubuntu Linux 16.04 LTS, R-release, GCC
    * Fedora Linux, R-devel, clang, gfortran

## R CMD check results
* ERRORs: 0
* WARNINGs: 0
* NOTEs: 3
    * The first states that `LMDIR` is a new submission to CRAN.
    * The second states that the Author field in DESCRIPTION
      differs from that derived from Authors@R.
      I have include my ORCID for `pgkdown` per the documentation found at
      <https://pkgdown.r-lib.org/articles/pkgdown.html#home-page>.
    * The third states that there is 
      no visible binding for global variable `.`.
      This note arises from the need to refer to the current group 
      within a call to `dplyr::do()`.
      This approach is appropriate via the documentation at
      <https://dplyr.tidyverse.org/reference/do.html>.

## Downstream dependencies
There are currently no downstream dependencies for `LMDIR`.