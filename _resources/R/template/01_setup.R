############################################################
# 01 - SETUP
############################################################

# Load package and metabarcoding helper functions
source("./R_functions/load_packages.R")
source("./R_functions/metabarcoding_tools.R")

# Check, optionally install, and load required packages
load_metabarcoding_packages()


# ==========================================================
# OUTPUT DIRECTORIES

dir.create("plots", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)


# ==========================================================
# PLOT COLOUR SCHEMES

## Discrete categories
# Publication-friendly colours. Uses ggsci for up to 10 groups and automatically
# switches to a broader qualitative palette when many categories are required.

palette_discrete <- function(n) {

  if (!is.numeric(n) || length(n) != 1L || n < 1) {
    stop("n must be a positive integer.", call. = FALSE)
  }

  n <- as.integer(n)

  if (n <= 10) {
    ggsci::pal_npg("nrc")(10)[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

## Alternative discrete palette
# Useful when a second independent grouping variable should look visually distinct.
palette_discrete_alt <- function(n) {
  if (n <= 9) {
    ggsci::pal_lancet("lanonc")(9)[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

## Continuous data
palette_continuous <- function(n = 100) {
  grDevices::hcl.colors(n, palette = "Viridis")
}

## Diverging data
# Useful for effects centred around zero, correlations, deviations, etc.
palette_diverging <- function(n = 100) {
  grDevices::hcl.colors(n, palette = "Blue-Red 3")
}


