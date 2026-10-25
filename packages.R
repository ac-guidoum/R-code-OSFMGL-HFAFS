## =============================================================================
##  packages.R -- package calls and global settings
##
##  Supplementary R code for:
##   Guidoum, A.C.; Bassoudi, M.; Boukhetala, K.
##   "One-Sided Fractional Moments of Gaussian Laws via Hermite Functions
##    with Applications to Fractional SDEs" 
## =============================================================================

pkgs <- c("stats", "graphics", "grDevices", "ggplot2", "scales", "patchwork",
          "gsl", "yuima", "parallel", "foreach", "doParallel", "doRNG",
          "dplyr", "ggsci", "MASS", "cowplot")
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

options(digits = 10)
n_cores <- max(1, detectCores() - 1)
