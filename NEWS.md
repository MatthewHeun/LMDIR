# News for `LMDIR`

## LMDIR 0.1.3

* Now calculating the D vector by D_j = exp(deltaV_j / L(V)).
  Doing so reduces the possibility that we'll hit an overflow error.
  Previous approach was D_j = exp(deltaV_j)^(1/L(V)).
  The two approaches are equivalent, mathematically.
  But the new approach avoids taking exp(deltaV_j), 
  which could overflow when deltaV_j was greater than about 725.


## LMDIR 0.1.2

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


## LMDIR 0.1.1

* `lmdi` function now creates first row with 0s (for $\Delta V$ terms) and 1s (for $D$ terms).
  This change means that the outgoing data frame has same number of rows as the incoming data frame,
  eliminating the need for head or tail padding.


## LMDIR 0.1.0

Initial version.