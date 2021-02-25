################################################################################
#### Update Step Length and Turning Angle Distributions
################################################################################
# Here, I'll follow Fieberg et al. to update the tentative parameters for our
# step length and turning angle distributions. We will use our back-transformed
# movement model for this.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling

# Load the backtransformed coefficients and our distributions
mod <- read_rds("03_Data/03_Results/99_MovementModelBacktransformed.rds")
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")
