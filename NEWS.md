---
title: "Release notes for `LMDIR`"
output: html_document
---


# LMDIR 0.1.13 (2020-04-12)

* Transition from Travis to GitHub actions for continuous integration testing.
* No new tests.
    * Still at 80 tests, all passing.
    * Test coverage remains at 100 %.


# LMDIR 0.1.12 (2020-05-01)

* Update to automatically wrap argument `fillrow` in a list
  in `Z_byname()`.


# LMDIR 0.1.11 (2020-04-23)

* Updates to keep up with changes to the `matsindf` package.


# LMDIR 0.1.10 (2020-02-12)

* Added badges on README file.
* Added automatic code coverage calculation.
* All function calls are now fully-qualified.
* Now at 100% testing coverage.


# LMDIR 0.1.9 (2019-01-07)

* Compatibility with removal of `element` prefixes to function names in `matsbyname`.


# LMDIR 0.1.8 (2019-01-02)

* Vignette no longer needs correct order for fillrow columns.


# LMDIR 0.1.7 (2019-01-02)

* Improved the vignette with details of replacing missing values.


# LMDIR 0.1.6 (2019-01-01)

* Fixed a bug in the website. Vignette wasn't showing up due to capitalization error.


# LMDIR 0.1.5 (2019-01-01)

* Now includes a detailed vignette.
* Now includes a new dataset `XGH` for the vignette.
* Output from `lmdi()` now includes `X`, `V`, and `Z` columns.


# LMDIR 0.1.4 (2018-12-29)

* Now including pkgdown website.
* Removed `_colname` from argument names.
* When calling `matsindf::collapse_to_matrices`, 
  now using `matvals` argument name instead of `values`.


# LMDIR 0.1.3

* Now calculating the D vector by D_j = exp(deltaV_j / L(V)).
  Doing so reduces the possibility that we'll hit an overflow error.
  Previous approach was D_j = exp(deltaV_j)^(1/L(V)).
  The two approaches are equivalent, mathematically.
  But the new approach avoids taking exp(deltaV_j), 
  which could overflow when deltaV_j was greater than about 725.


# LMDIR 0.1.2

* `lmdi` now defaults to supplying very small values for missing categories (1e-10)
  when calculating `Z` matrices.
* Both `lmdi` and `Z_byname` functions now accept `fillrow` arguments.
  `fillrow` provides an alternative to assuming small values for missing categories.
  The `fillrow` approach is especially useful when a row in X_0 or X_T 
  is missing compared to the other,
  due to a particular category being present in one year but absent in the other.
  In this situation, we want to fill values in the missing row with non-zero numbers for
  all categories except for allocation to the final subsubcategory.
  Then, we set only the allocation from subcategory to subsubcategory (phi_ij) to 0.
  This approach correctly models the fact that
  despite the fact that we have no useful exergy of this type being produced in
  one of the years, we still have
  total primary en/xergy (E.ktoe),
  there is still an allocation of primary exergy to the subcategory (phi_i), and
  if there were machines making this subsubcategory of useful exergy
  in this time period, they would have a certain primary-to-useful efficiency (eta_ij).
  It turns out that we don't need to know the exact values of
  primary exergy (E.ktoe),
  allocation to subcategory (phi_i), or
  primary-to-useful efficiency (eta_ij).
  These values must simply be non-zero so long as
  allocation from subcategory to subsubcategory (phi_ij) is zero.


# LMDIR 0.1.1

* `lmdi()` now creates first row with 0s (for $\Delta V$ terms) and 1s (for $D$ terms).
  This change means that the outgoing data frame has same number of rows as the incoming data frame,
  eliminating the need for head or tail padding.


# LMDIR 0.1.0

Initial version.