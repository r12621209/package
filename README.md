This package includes three functions:

Simulates 3 environments with 3-fold and 4-fold cross-validation (CV-O, CV-NO)

Simulates 5 environments with 5-fold and 8-fold cross-validation (CV-O, CV-NO)

Performs 3-fold and 4-fold cross-validation with CV-O and CV-NO

## Installation

You can install from GitHub using the `devtools` package:

```r
# First, make sure devtools is installed
install.packages("devtools")

# Then install EHPGS from GitHub
library(devtools)
install_github("r12621209/package", dependencies = TRUE, force = TRUE)
