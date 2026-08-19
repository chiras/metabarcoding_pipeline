############################################################
# Marker Setup
############################################################

if (!exists("marker") || is.null(marker) || !nzchar(marker)) {
  stop(
    "Marker is not defined. Expected 'ITS2', 'COI', or '16S'.",
    call. = FALSE
  )
}

marker <- toupper(marker)

if (!marker %in% c("ITS2", "COI", "16S", "fITS")) {
  stop(
    "Unsupported marker: ", marker,
    ". Expected 'ITS2', 'COI', 'fITS' or '16S'.",
    call. = FALSE
  )
}

analysis_rank <- switch(
  marker,
  ITS2 = "species",
  COI  = "species",
  fITS  = "species",
  `16S` = "genus"
)

############################################################
# Package setup for metabarcoding downstream analyses
############################################################

load_metabarcoding_packages <- function(ask_install = interactive()) {

  cran_packages <- c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "stringr",
    "reshape2",
    "ggraph",
    "iNEXT",
    "ggsci",
    "bipartite",
    "ape"
  )

  bioc_packages <- c(
    "phyloseq",
    "speedyseq"
  )

  all_packages <- c(cran_packages, bioc_packages)

  missing <- all_packages[
    !vapply(
      all_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(missing) > 0L) {

    cat("\nMissing R packages:\n")
    cat("  ", paste(missing, collapse = ", "), "\n\n", sep = "")

    install_missing <- FALSE

    if (isTRUE(ask_install)) {
      answer <- readline(
        "Install missing packages now? [y/N]: "
      )

      install_missing <- tolower(trimws(answer)) %in% c(
        "y", "yes"
      )
    }

    if (!install_missing) {
      stop(
        "Required packages are missing: ",
        paste(missing, collapse = ", "),
        "\nInstall them and rerun the script.",
        call. = FALSE
      )
    }

    # CRAN packages
    missing_cran <- intersect(missing, cran_packages)

    if (length(missing_cran) > 0L) {
      message(
        "Installing CRAN packages: ",
        paste(missing_cran, collapse = ", ")
      )

      install.packages(missing_cran)
    }

    # Bioconductor packages
    missing_bioc <- intersect(missing, bioc_packages)

    if (length(missing_bioc) > 0L) {

      if (!requireNamespace("BiocManager", quietly = TRUE)) {

        message("Installing BiocManager...")
        install.packages("BiocManager")
      }

      message(
        "Installing Bioconductor packages: ",
        paste(missing_bioc, collapse = ", ")
      )

      BiocManager::install(
        missing_bioc,
        ask = FALSE,
        update = FALSE
      )
    }
  }

  # Verify installation before attaching packages
  still_missing <- all_packages[
    !vapply(
      all_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(still_missing) > 0L) {
    stop(
      "Package installation failed or packages are still unavailable: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(phyloseq)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(speedyseq)
    library(ape)
    library(stringr)
    library(reshape2)
    library(ggsci)
    library(bipartite)
  })

  message(
    "Loaded ",
    length(all_packages),
    " required R packages."
  )

  invisible(TRUE)
}