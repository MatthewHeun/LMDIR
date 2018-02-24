# News for `LMDIR`

## LMDIR 0.1.1

* `lmdi` function now creates first row with 0s (for $\Delta V$ terms) and 1s (for $D$ terms).
This change means that the outgoing data frame has same number of rows as the incoming data frame,
eliminating the need for head or tail padding.


## LMDIR 0.1.0

Initial version.