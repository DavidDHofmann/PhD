################################################################################
#### Installing My Own Package
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(devtools)
library(roxygen2)

################################################################################
#### davidoff
################################################################################
setwd("/home/david/ownCloud/University/15. PhD/General/R-Packages")

# Create package
# create("davidoff")

# Once you added some files containing functions, process the documentation
setwd("davidoff")
use_rcpp()
document()

# Install the package
setwd("..")
install("davidoff")

################################################################################
#### floodmapr
################################################################################
setwd("/home/david/ownCloud/University/15. PhD/General/R-Packages")

# Create package
# create("floodmapr")

# Once you added some files containing functions, process the documentation
setwd("floodmapr")
document()

# Install the package
setwd("..")
install("floodmapr")
