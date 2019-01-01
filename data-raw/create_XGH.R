# This script documents the process of producing the XGH object used in the LMDIR.Rmd vignette
#
#   -- Matthew Kuperus Heun, 1 January 2019
#
# This script requires the availability of a data frame of X matrices for LMDI analysis.
# I first entered the Useful-Work-Analysis project and said

ProjectTemplate::load.project().

# load.project() creates an LMDI object that contains all of the results.

# Next, I ran the following code in the Useful-Work-Analysis project:

XGH <- LMDI %>%
  ungroup() %>%
  filter(Method == "PCM", Country == "GH", Last.stage == "useful",
         Energy.type == "X.ktoe", Nonenergy == "w/o non-energy") %>%
  select(Country, Year, X)

# At this point XGH contains only Country, Year, and X columns.
# To make XGH available to the LMDI package,
# I said the following in the Useful-Work-Analysis project:

saveRDS(XGH, "~/Desktop/XGH.Rds")

# Next, I switched to the LMDIR project and said:

XGH <- readRDS("~/Desktop/XGH.Rds")

# Finally, to create the object that will be a part of the package, I said:

use_data(XGH)

# Don't forget to document the data frame in R/data.R.