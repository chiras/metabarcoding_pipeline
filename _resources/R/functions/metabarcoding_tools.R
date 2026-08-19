# metabarcoding_tools_0-2c.R
# Utility functions for phyloseq-based metabarcoding workflows.
# Revision 0-2c: consolidated 0-2a robustness fixes with layered metadata handling.
# Includes improved argument validation, OTU-table orientation handling,
# reproducibility, portability, and removal of unnecessary package dependencies.


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

.assert_phyloseq <- function(ps, allow_empty = FALSE) {
  if (!inherits(ps, "phyloseq")) {
    stop("Input must be a phyloseq object.", call. = FALSE)
  }

  if (!allow_empty) {
    if (phyloseq::ntaxa(ps) == 0L) {
      stop("Phyloseq object has zero taxa.", call. = FALSE)
    }
    if (phyloseq::nsamples(ps) == 0L) {
      stop("Phyloseq object has zero samples.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.otu_taxa_rows <- function(ps) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  otu
}

.restore_otu_orientation <- function(mat_taxa_rows, taxa_rows) {
  if (taxa_rows) mat_taxa_rows else t(mat_taxa_rows)
}

.validate_positive_scalar <- function(x, name, allow_zero = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!ok || (!allow_zero && x <= 0) || (allow_zero && x < 0)) {
    comparator <- if (allow_zero) ">= 0" else "> 0"
    stop(name, " must be a single finite numeric value ", comparator, ".", call. = FALSE)
  }
  invisible(TRUE)
}
# -----------------------------------------------------------------------------
# Combine ASV and Species names into a single label for plotting
# -----------------------------------------------------------------------------
add_taxon_to_asv_id <- function(physeq, rank = NULL, sep = ".") {
  .assert_phyloseq(physeq)

  tax <- phyloseq::tax_table(physeq)

  # Default to the lowest/final taxonomic rank
  if (is.null(rank)) {
    rank <- colnames(tax)[ncol(tax)]
  }

  if (!rank %in% colnames(tax)) {
    stop(
      "Taxonomic rank '", rank, "' not found in tax_table.",
      call. = FALSE
    )
  }

  taxon <- as.character(tax[, rank])
  asv <- phyloseq::taxa_names(physeq)

  # Avoid creating unusable IDs from missing taxonomy
  missing_taxon <- is.na(taxon) | taxon == ""
  taxon[missing_taxon] <- "unresolved"

  phyloseq::taxa_names(physeq) <- paste0(taxon, sep, asv)

  physeq
}

# -----------------------------------------------------------------------------
# Empty taxa / samples and abundance filtering
# -----------------------------------------------------------------------------

drop_empty_taxa_samples <- function(ps, verbose = TRUE, context = "Empty filtering") {
  .assert_phyloseq(ps)

  n_taxa_before <- phyloseq::ntaxa(ps)
  n_samples_before <- phyloseq::nsamples(ps)

  empty_taxa <- phyloseq::taxa_sums(ps) == 0
  n_empty_taxa <- sum(empty_taxa)

  if (n_empty_taxa > 0L) {
    ps <- phyloseq::prune_taxa(!empty_taxa, ps)
  }

  if (phyloseq::ntaxa(ps) == 0L) {
    stop("All taxa have zero abundance after pruning.", call. = FALSE)
  }

  empty_samples <- phyloseq::sample_sums(ps) == 0
  n_empty_samples <- sum(empty_samples)

  if (n_empty_samples > 0L) {
    ps <- phyloseq::prune_samples(!empty_samples, ps)
  }

  if (phyloseq::nsamples(ps) == 0L) {
    stop("All samples have zero abundance after pruning.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(
      context, ": removed ", n_empty_taxa, " empty taxa and ",
      n_empty_samples, " empty samples."
    )
  }

  list(
    phyloseq = ps,
    n_taxa_before = n_taxa_before,
    n_taxa_after = phyloseq::ntaxa(ps),
    n_samples_before = n_samples_before,
    n_samples_after = phyloseq::nsamples(ps),
    n_empty_taxa_removed = n_empty_taxa,
    n_empty_samples_removed = n_empty_samples
  )
}


filter_low_abundance <- function(ps_counts, threshold = 0.01, min_samples = 1L,
                                 verbose = TRUE) {
  .assert_phyloseq(ps_counts)
  .validate_positive_scalar(threshold, "threshold", allow_zero = TRUE)

  if (threshold > 1) {
    stop("threshold is a relative-abundance cutoff and must be <= 1.", call. = FALSE)
  }
  if (!is.numeric(min_samples) || length(min_samples) != 1L ||
      !is.finite(min_samples) || min_samples < 1 || min_samples %% 1 != 0) {
    stop("min_samples must be a positive integer.", call. = FALSE)
  }
  min_samples <- as.integer(min_samples)

  empty_pre <- drop_empty_taxa_samples(
    ps_counts,
    verbose = verbose,
    context = "Before low-abundance filtering"
  )
  ps_counts <- empty_pre$phyloseq

  if (min_samples > phyloseq::nsamples(ps_counts)) {
    stop(
      "min_samples (", min_samples, ") exceeds the number of retained samples (",
      phyloseq::nsamples(ps_counts), ").",
      call. = FALSE
    )
  }

  taxa_rows <- phyloseq::taxa_are_rows(ps_counts)
  otu_counts <- .otu_taxa_rows(ps_counts)

  sample_totals <- colSums(otu_counts, na.rm = TRUE)
  if (any(sample_totals <= 0)) {
    stop("Internal error: empty samples remain after pruning.", call. = FALSE)
  }

  otu_rel <- sweep(otu_counts, 2L, sample_totals, "/")

  # A taxon is retained if it reaches the threshold in at least min_samples.
  keep_cell <- otu_rel >= threshold
  keep_taxa <- rowSums(keep_cell, na.rm = TRUE) >= min_samples
  n_taxa_removed_low_abundance <- sum(!keep_taxa)

  if (!any(keep_taxa)) {
    stop(
      "Filtering would remove all taxa. No taxon reaches threshold = ",
      threshold, " in at least ", min_samples, " sample(s).",
      call. = FALSE
    )
  }

  otu_counts_filtered <- otu_counts[keep_taxa, , drop = FALSE]
  otu_rel_filtered <- otu_rel[keep_taxa, , drop = FALSE]

  # Values below the threshold are removed even for taxa retained elsewhere.
  low_cell <- otu_rel_filtered < threshold
  otu_counts_filtered[low_cell] <- 0
  otu_rel_filtered[low_cell] <- 0

  keep_samples <- colSums(otu_counts_filtered, na.rm = TRUE) > 0
  n_samples_removed_low_abundance <- sum(!keep_samples)

  if (!any(keep_samples)) {
    stop("Filtering would remove all samples.", call. = FALSE)
  }

  otu_counts_filtered <- otu_counts_filtered[, keep_samples, drop = FALSE]
  otu_rel_filtered <- otu_rel_filtered[, keep_samples, drop = FALSE]

  keep_taxa_after <- rowSums(otu_counts_filtered, na.rm = TRUE) > 0
  n_taxa_removed_empty_after_low <- sum(!keep_taxa_after)

  if (!any(keep_taxa_after)) {
    stop("Filtering would remove all taxa after empty-sample pruning.", call. = FALSE)
  }

  otu_counts_filtered <- otu_counts_filtered[keep_taxa_after, , drop = FALSE]
  otu_rel_filtered <- otu_rel_filtered[keep_taxa_after, , drop = FALSE]

  kept_taxa_names <- rownames(otu_counts_filtered)
  kept_sample_names <- colnames(otu_counts_filtered)

  ps_counts_filtered <- phyloseq::prune_taxa(kept_taxa_names, ps_counts)
  ps_counts_filtered <- phyloseq::prune_samples(kept_sample_names, ps_counts_filtered)

  ps_rel_filtered <- phyloseq::prune_taxa(kept_taxa_names, ps_counts)
  ps_rel_filtered <- phyloseq::prune_samples(kept_sample_names, ps_rel_filtered)

  phyloseq::otu_table(ps_counts_filtered) <- phyloseq::otu_table(
    .restore_otu_orientation(otu_counts_filtered, taxa_rows),
    taxa_are_rows = taxa_rows
  )
  phyloseq::otu_table(ps_rel_filtered) <- phyloseq::otu_table(
    .restore_otu_orientation(otu_rel_filtered, taxa_rows),
    taxa_are_rows = taxa_rows
  )

  if (isTRUE(verbose)) {
    message(
      "Low-abundance filtering: removed ", n_taxa_removed_low_abundance,
      " taxa by threshold, ", n_samples_removed_low_abundance,
      " samples that became empty, and ", n_taxa_removed_empty_after_low,
      " additional empty taxa."
    )
    message(
      "Final object: ", phyloseq::ntaxa(ps_counts_filtered), " taxa and ",
      phyloseq::nsamples(ps_counts_filtered), " samples retained."
    )
  }

  list(
    counts = ps_counts_filtered,
    relative = ps_rel_filtered,
    n_taxa_before = empty_pre$n_taxa_before,
    n_taxa_after = phyloseq::ntaxa(ps_counts_filtered),
    n_samples_before = empty_pre$n_samples_before,
    n_samples_after = phyloseq::nsamples(ps_counts_filtered),
    n_empty_taxa_removed_before = empty_pre$n_empty_taxa_removed,
    n_empty_samples_removed_before = empty_pre$n_empty_samples_removed,
    n_taxa_removed_low_abundance = n_taxa_removed_low_abundance,
    n_samples_removed_low_abundance = n_samples_removed_low_abundance,
    n_empty_taxa_removed_after_low_abundance = n_taxa_removed_empty_after_low,
    threshold = threshold,
    min_samples = min_samples
  )
}


# -----------------------------------------------------------------------------
# Metadata helpers
# -----------------------------------------------------------------------------

.read_metadata_file <- function(file) {
  if (is.null(file)) return(NULL)
  if (!is.character(file) || length(file) != 1L || !nzchar(file)) {
    stop("Metadata file must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(file)) {
    stop("Metadata file not found: ", file, call. = FALSE)
  }

  # check.names = FALSE is important for metadata templates containing names
  # such as '*sample_name'. Empty strings are retained as NA where possible.
  utils::read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

.resolve_metadata_id_column <- function(metadata, sample_id_col = "sample_name") {
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame.", call. = FALSE)
  }

  if (!is.null(sample_id_col) && sample_id_col %in% names(metadata)) {
    return(sample_id_col)
  }

  candidates <- c("sample_name", "SampleID", "sample_id", "sample", "Sample")
  hits <- candidates[candidates %in% names(metadata)]
  if (length(hits) == 1L) return(hits)

  stop(
    "Could not identify the sample-ID column in metadata. Specify sample_id_col explicitly.",
    call. = FALSE
  )
}

append_sample_metadata <- function(phyloseq, metadata,
                                   sample_id_col = "sample_name",
                                   conflict = c("fill_missing", "replace", "keep"),
                                   drop_unmatched_metadata = TRUE,
                                   verbose = TRUE) {
  .assert_phyloseq(phyloseq)
  conflict <- match.arg(conflict)

  if (is.character(metadata) && length(metadata) == 1L) {
    metadata <- .read_metadata_file(metadata)
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame or path to a CSV file.", call. = FALSE)
  }

  sample_id_col <- .resolve_metadata_id_column(metadata, sample_id_col)
  ids <- as.character(metadata[[sample_id_col]])

  if (anyNA(ids) || any(!nzchar(ids))) {
    stop("Metadata contains missing or empty sample IDs in column '", sample_id_col, "'.",
         call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    duplicates <- unique(ids[duplicated(ids)])
    stop(
      "Metadata contains duplicated sample IDs: ",
      paste(utils::head(duplicates, 10L), collapse = ", "),
      if (length(duplicates) > 10L) " ..." else "",
      call. = FALSE
    )
  }

  ps_ids <- phyloseq::sample_names(phyloseq)
  matched <- match(ps_ids, ids)
  n_matched <- sum(!is.na(matched))

  if (n_matched == 0L) {
    stop("None of the metadata sample IDs match the phyloseq sample names.", call. = FALSE)
  }

  unmatched_ps <- ps_ids[is.na(matched)]
  unmatched_md <- setdiff(ids, ps_ids)

  # Preserve all phyloseq samples and their current metadata. Metadata are
  # aligned to sample_names rather than relying on input row order.
  if (!is.null(phyloseq::sample_data(phyloseq, errorIfNULL = FALSE))) {
    current <- as.data.frame(phyloseq::sample_data(phyloseq), stringsAsFactors = FALSE)
  } else {
    current <- data.frame(row.names = ps_ids)
  }
  current <- current[ps_ids, , drop = FALSE]

  incoming_cols <- setdiff(names(metadata), sample_id_col)
  incoming <- metadata[matched, incoming_cols, drop = FALSE]
  rownames(incoming) <- ps_ids

  # For unmatched phyloseq samples, matched indexing creates NA rows. This is
  # intentional: metadata attachment must never silently remove samples.
  overlap <- intersect(names(current), names(incoming))
  new_cols <- setdiff(names(incoming), names(current))

  for (nm in new_cols) {
    current[[nm]] <- incoming[[nm]]
  }

  for (nm in overlap) {
    if (conflict == "replace") {
      replace_idx <- !is.na(incoming[[nm]])
      current[[nm]][replace_idx] <- incoming[[nm]][replace_idx]
    } else if (conflict == "fill_missing") {
      empty_current <- is.na(current[[nm]]) |
        (is.character(current[[nm]]) & !nzchar(trimws(current[[nm]])))
      fill_idx <- empty_current & !is.na(incoming[[nm]])
      current[[nm]][fill_idx] <- incoming[[nm]][fill_idx]
    }
    # conflict == "keep": retain existing column unchanged.
  }

  phyloseq::sample_data(phyloseq) <- phyloseq::sample_data(current)

  if (isTRUE(verbose)) {
    message(
      "Metadata appended: ", n_matched, "/", length(ps_ids),
      " phyloseq samples matched; ", length(incoming_cols), " metadata columns considered."
    )
    if (length(unmatched_ps) > 0L) {
      message("  Phyloseq samples without metadata: ", length(unmatched_ps), ".")
    }
    if (!isTRUE(drop_unmatched_metadata) && length(unmatched_md) > 0L) {
      message("  Metadata rows not represented in phyloseq: ", length(unmatched_md), ".")
    }
  }

  phyloseq
}

read_sample_metadata_table <- function(file,
                                       sample_id_col = NULL,
                                       sep = ",",
                                       normalize_hyphens = FALSE,
                                       id_substitutions = NULL) {

  # Read one metadata table independently, before any merging.
  # Examples:
  #   CSV:       sep = ","
  #   TSV:       sep = "\t"
  #   semicolon: sep = ";"
  #
  # Sample names can subsequently be modified manually, e.g.:
  #   sample_names(qc.map) <- gsub(
  #     ".ITS_S1", "", sample_names(qc.map), fixed = TRUE
  #   )

  if (is.null(file) || !file.exists(file)) {
    stop("Metadata file not found: ", file, call. = FALSE)
  }

  x <- utils::read.table(
    file,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    quote = "\"",
    comment.char = ""
  )

  if (nrow(x) == 0L) {
    stop("Metadata file contains no rows: ", file, call. = FALSE)
  }

  # Use first column as sample ID unless another column is specified.
  if (is.null(sample_id_col)) {
    ids <- as.character(x[[1L]])
    x <- x[-1L]
  } else {
    if (!sample_id_col %in% names(x)) {
      stop(
        "Sample ID column '", sample_id_col,
        "' not found in ", file, ".",
        call. = FALSE
      )
    }

    ids <- as.character(x[[sample_id_col]])
    x[[sample_id_col]] <- NULL
  }

  if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
    stop(
      "Missing/empty sample IDs in ", file, ".",
      call. = FALSE
    )
  }

  # Optional convenience conversion.
  if (isTRUE(normalize_hyphens)) {
    ids <- gsub("-", ".", ids, fixed = TRUE)
  }

  # Optional ordered gsub rules.
  #
  # Example:
  # id_substitutions = list(
  #   list(
  #     pattern = ".ITS_S1",
  #     replacement = "",
  #     fixed = TRUE
  #   ),
  #   list(
  #     pattern = "-",
  #     replacement = ".",
  #     fixed = TRUE
  #   )
  # )

  if (!is.null(id_substitutions)) {
    if (!is.list(id_substitutions)) {
      stop(
        "id_substitutions must be NULL or a list of gsub rule lists.",
        call. = FALSE
      )
    }

    for (rule in id_substitutions) {
      if (!is.list(rule) || is.null(rule$pattern)) {
        stop(
          "Each id_substitutions rule must contain at least 'pattern'.",
          call. = FALSE
        )
      }

      replacement <- if (is.null(rule$replacement)) {
        ""
      } else {
        rule$replacement
      }

      fixed <- if (is.null(rule$fixed)) {
        FALSE
      } else {
        isTRUE(rule$fixed)
      }

      ignore.case <- if (is.null(rule$ignore.case)) {
        FALSE
      } else {
        isTRUE(rule$ignore.case)
      }

      ids <- gsub(
        pattern = rule$pattern,
        replacement = replacement,
        x = ids,
        fixed = fixed,
        ignore.case = ignore.case
      )
    }
  }

  # Check uniqueness after all transformations.
  if (anyDuplicated(ids)) {

    dup <- unique(ids[duplicated(ids)])

    stop(
      "Duplicated sample IDs in ", file,
      " after sample-name transformations: ",
      paste(head(dup, 10L), collapse = ", "),
      if (length(dup) > 10L) " ..." else "",
      call. = FALSE
    )
  }

  rownames(x) <- ids
  phyloseq::sample_data(x)
}


merge_sample_metadata <- function(qc_metadata,
                                  biological_metadata = NULL,
                                  generate_pseudo = FALSE,
                                  pseudo = list(
                                    location = c("DummyLoc1", "DummyLoc2", "DummyLoc3"),
                                    treatment = c("DummyTreat1", "DummyTreat2", "DummyTreat3")
                                  ),
                                  pseudo_prefix = "pseudo_",
                                  seed = NULL,
                                  verbose = TRUE,
                                  preview_n = 10L) {
  # Merge already-loaded metadata objects after any desired sample-name cleanup.
  # QC metadata define the sample universe.

  if (!inherits(qc_metadata, "sample_data")) {
    stop("qc_metadata must be a phyloseq::sample_data object.", call. = FALSE)
  }
  if (!is.null(biological_metadata) && !inherits(biological_metadata, "sample_data")) {
    stop("biological_metadata must be NULL or a phyloseq::sample_data object.", call. = FALSE)
  }

  if (!is.numeric(preview_n) || length(preview_n) != 1L ||
      is.na(preview_n) || preview_n < 0) {
    stop("preview_n must be a single non-negative number.", call. = FALSE)
  }
  preview_n <- as.integer(preview_n)

  md <- as.data.frame(qc_metadata, stringsAsFactors = FALSE)
  qc_ids <- rownames(md)

  bio <- NULL
  bio_ids <- character(0)
  overlap_ids <- character(0)
  qc_only_ids <- qc_ids
  bio_only_ids <- character(0)

  if (!is.null(biological_metadata)) {
    bio <- as.data.frame(biological_metadata, stringsAsFactors = FALSE)
    bio_ids <- rownames(bio)

    overlap_ids <- intersect(qc_ids, bio_ids)
    qc_only_ids <- setdiff(qc_ids, bio_ids)
    bio_only_ids <- setdiff(bio_ids, qc_ids)

    duplicate_cols <- intersect(names(md), names(bio))

    if (length(duplicate_cols) > 0L) {

      new_qc_names <- paste0("QC.", duplicate_cols)

      if (any(new_qc_names %in% names(md)) ||
          any(new_qc_names %in% names(bio))) {
        stop(
          "Cannot rename duplicated QC metadata columns because one or more ",
          "'QC.<column>' names already exist: ",
          paste(
            new_qc_names[
              new_qc_names %in% c(names(md), names(bio))
            ],
            collapse = ", "
          ),
          ".",
          call. = FALSE
        )
      }

      names(md)[match(duplicate_cols, names(md))] <- new_qc_names

      if (isTRUE(verbose)) {
        message(
          "Metadata columns present in both tables: ",
          paste(duplicate_cols, collapse = ", "),
          ". QC versions renamed to: ",
          paste(new_qc_names, collapse = ", "),
          "."
        )
      }
    }

    for (nm in names(bio)) {
      md[[nm]] <- NA
      if (length(overlap_ids) > 0L) {
        md[overlap_ids, nm] <- bio[overlap_ids, nm]
      }
    }
  }

  all_ids <- union(qc_ids, bio_ids)
  metadata_check <- data.frame(
    sample_name = all_ids,
    qc_metadata = all_ids %in% qc_ids,
    biological_metadata = all_ids %in% bio_ids,
    stringsAsFactors = FALSE
  )
  metadata_check$status <- ifelse(
    metadata_check$qc_metadata & metadata_check$biological_metadata,
    "overlap",
    ifelse(metadata_check$qc_metadata, "QC only", "biological only")
  )
  metadata_check <- metadata_check[
    order(
      match(metadata_check$status, c("overlap", "QC only", "biological only")),
      metadata_check$sample_name
    ),
    , drop = FALSE
  ]
  rownames(metadata_check) <- NULL

  if (isTRUE(verbose)) {
    message("Metadata check:")
    message("  QC metadata samples:           ", length(qc_ids))
    message("  Biological metadata samples:   ", length(bio_ids))
    message("  Present in both:               ", length(overlap_ids))
    message("  QC only:                       ", length(qc_only_ids))
    message("  Biological only:               ", length(bio_only_ids))

    if (is.null(biological_metadata)) {
      message("  Biological metadata: not supplied; skipped.")
    } else if (preview_n > 0L) {
      if (length(qc_only_ids) > 0L) {
        cat("\nQC-only samples (first ",
            min(preview_n, length(qc_only_ids)), " of ",
            length(qc_only_ids), "):\n", sep = "")
        print(utils::head(qc_only_ids, preview_n))
      }

      if (length(bio_only_ids) > 0L) {
        cat("\nBiological-metadata-only samples (first ",
            min(preview_n, length(bio_only_ids)), " of ",
            length(bio_only_ids), "):\n", sep = "")
        print(utils::head(bio_only_ids, preview_n))
      }

      if (length(qc_only_ids) == 0L && length(bio_only_ids) == 0L) {
        message("  Sample IDs match completely between QC and biological metadata.")
      }
    }
  }

  if (isTRUE(generate_pseudo)) {
    if (!is.list(pseudo) || is.null(names(pseudo)) || any(!nzchar(names(pseudo)))) {
      stop(
        "pseudo must be a named list, e.g. list(location = c('DummyLoc1', 'DummyLoc2')).",
        call. = FALSE
      )
    }
    if (!is.character(pseudo_prefix) || length(pseudo_prefix) != 1L || is.na(pseudo_prefix)) {
      stop("pseudo_prefix must be a single character string.", call. = FALSE)
    }
    if (!is.null(seed)) set.seed(seed)

    generated <- character(0)
    for (nm in names(pseudo)) {
      values <- pseudo[[nm]]
      if (length(values) == 0L || all(is.na(values))) {
        stop("Pseudo-metadata values for '", nm, "' are empty.", call. = FALSE)
      }

      target <- paste0(pseudo_prefix, nm)
      if (target %in% names(md)) {
        stop("Pseudo metadata column already exists: ", target, call. = FALSE)
      }

      md[[target]] <- sample(values, nrow(md), replace = TRUE)
      generated <- c(generated, target)
    }

    if (isTRUE(verbose)) {
      message("Pseudo metadata generated: ", paste(generated, collapse = ", "), ".")
    }
  }

  data_map <- phyloseq::sample_data(md)

  structure(
    list(
      metadata = data_map,
      overlap = metadata_check,
      only_in_QC = qc_only_ids,
      only_in_biological = bio_only_ids
    ),
    class = c("sample_metadata_layers", "list")
  )
}


load_sample_metadata_layers <- function(metadata_file = "samples_metadata.csv",
                                        biological_metadata_file = NULL,
                                        sample_id_col = NULL,
                                        biological_sample_id_col = NULL,
                                        normalize_hyphens = FALSE,
                                        qc_id_substitutions = NULL,
                                        biological_id_substitutions = NULL,
                                        generate_pseudo = FALSE,
                                        pseudo = list(
                                          location = c("DummyLoc1", "DummyLoc2", "DummyLoc3"),
                                          treatment = c("DummyTreat1", "DummyTreat2", "DummyTreat3")
                                        ),
                                        pseudo_prefix = "pseudo_",
                                        seed = NULL,
                                        verbose = TRUE,
                                        preview_n = 10L) {
  # Convenience wrapper around read_sample_metadata_table() + merge_sample_metadata().
  # For maximum transparency/flexibility, the two-step workflow is recommended.

  qc_map <- read_sample_metadata_table(
    file = metadata_file,
    sample_id_col = sample_id_col,
    normalize_hyphens = normalize_hyphens,
    id_substitutions = qc_id_substitutions
  )

  if (isTRUE(verbose)) {
    message(
      "Technical/QC metadata loaded: ",
      phyloseq::nsamples(qc_map), " samples, ",
      ncol(as.data.frame(qc_map)), " columns."
    )
  }

  bio_map <- NULL
  if (!is.null(biological_metadata_file)) {
    bio_map <- read_sample_metadata_table(
      file = biological_metadata_file,
      sample_id_col = biological_sample_id_col,
      normalize_hyphens = normalize_hyphens,
      id_substitutions = biological_id_substitutions
    )
  }

  merge_sample_metadata(
    qc_metadata = qc_map,
    biological_metadata = bio_map,
    generate_pseudo = generate_pseudo,
    pseudo = pseudo,
    pseudo_prefix = pseudo_prefix,
    seed = seed,
    verbose = verbose,
    preview_n = preview_n
  )
}

# Backwards-compatible alias for code that used the old helper name.
# The function returns a list; use $metadata to obtain the sample_data object.
load_sample_metadata <- load_sample_metadata_layers

# -----------------------------------------------------------------------------
# Taxonomy helpers
# -----------------------------------------------------------------------------
select_taxonomic_rank <- function(
    physeq,
    target_taxa = 20,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  if (!is.numeric(target_taxa) ||
      length(target_taxa) != 1L ||
      !is.finite(target_taxa) ||
      target_taxa < 1) {
    stop("target_taxa must be a positive number.", call. = FALSE)
  }

  tax <- as.data.frame(
    phyloseq::tax_table(physeq),
    stringsAsFactors = FALSE
  )

  if (ncol(tax) == 0L) {
    stop("No taxonomic ranks found in tax_table.", call. = FALSE)
  }

  ## Count resolved categories at each available taxonomic rank
  rank_summary <- data.frame(
    rank = names(tax),
    n_taxa = vapply(
      tax,
      function(x) {
        x <- as.character(x)
        x <- x[!is.na(x) & x != ""]
        length(unique(x))
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  ## Ignore ranks without any assignments
  rank_summary <- rank_summary[
    rank_summary$n_taxa > 0,
    ,
    drop = FALSE
  ]

  if (nrow(rank_summary) == 0L) {
    stop(
      "No populated taxonomic ranks available.",
      call. = FALSE
    )
  }

  ## Select rank closest to requested number of categories
  rank_summary$difference <- abs(
    rank_summary$n_taxa - target_taxa
  )

  best <- rank_summary[
    which.min(rank_summary$difference),
    ,
    drop = FALSE
  ]

  selected_rank <- best$rank[[1]]

  if (isTRUE(verbose)) {
    message(
      "Taxonomic rank selected for plotting: ",
      selected_rank,
      " (", best$n_taxa[[1]], " categories; target = ",
      target_taxa, ")."
    )
  }

  selected_rank
}

tax_glom_postcluster_species <- function(
    physeq,
    sp_pattern = "_sp| spc| sp\\.| spc\\.",
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  tax_df <- as.data.frame(
    phyloseq::tax_table(physeq),
    stringsAsFactors = FALSE
  )

  if (!"species" %in% names(tax_df)) {
    stop(
      "Taxonomy table must contain a 'species' column.",
      call. = FALSE
    )
  }

  species_values <- as.character(tax_df$species)

  # Identify unresolved species-level assignments.
  # Examples after taxonomy propagation:
  # Trifolium_spc, Asteraceae_spc, Viridiplantae_spc
  is_placeholder <- !is.na(species_values) &
    grepl(
      sp_pattern,
      species_values,
      ignore.case = TRUE,
      perl = TRUE
    )

  # Keep each unresolved postcluster separately by numbering identical
  # placeholder labels. This preserves an approximation of hidden species
  # diversity instead of collapsing all "Trifolium_spc" into one taxon.
  if (any(is_placeholder)) {

    placeholder_idx <- which(is_placeholder)
    placeholder_values <- species_values[placeholder_idx]

    within_placeholder_index <- ave(
      seq_along(placeholder_idx),
      placeholder_values,
      FUN = seq_along
    )

    tax_df$species[placeholder_idx] <- paste0(
      placeholder_values,
      " ",
      within_placeholder_index
    )
  }

  phyloseq::tax_table(physeq) <-
    phyloseq::tax_table(as.matrix(tax_df))

  # Resolved taxa with the same species assignment are merged.
  # Numbered unresolved postclusters remain separate.
  result <- phyloseq::tax_glom(
    physeq,
    taxrank = "species"
  )

  if (isTRUE(verbose)) {
    message(
      "Postcluster-aware species aggregation: ",
      phyloseq::ntaxa(physeq), " taxa -> ",
      phyloseq::ntaxa(result), " taxa; ",
      sum(is_placeholder),
      " unresolved postclustered units retained separately."
    )
  }

  result
}

select_analysis_unit <- function(
    physeq,
    method = c("asv", "taxonomy", "postcluster"),
    rank = NULL,
    postcluster_threshold = 0,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  method <- match.arg(
    tolower(method),
    choices = c("asv", "taxonomy", "postcluster")
  )

  postcluster_available <- !is.null(postcluster_threshold) &&
    !is.na(postcluster_threshold) &&
    postcluster_threshold > 0

  n_before <- phyloseq::ntaxa(physeq)

  if (method == "asv") {
    message(
      "Analysis unit: ",
      if (postcluster_available) {
        "pASVs (ASVs were postclustered during bioinformatics; no further aggregation)."
      } else {
        "ASVs (no taxonomic aggregation)."
      }
    )

    result <- physeq

    if (isTRUE(verbose)) {
      message(
        "Analysis unit: ",
        if (postcluster_available) "postclustered ASVs (pASVs)." else "ASVs."
      )
    }

  } else if (method == "taxonomy") {

    if (is.null(rank) || !nzchar(rank)) {
      stop(
        "A taxonomic rank must be supplied for method = 'taxonomy'.",
        call. = FALSE
      )
    }

    available_ranks <- colnames(phyloseq::tax_table(physeq))

    if (!rank %in% available_ranks) {
      stop(
        "Taxonomic rank '", rank, "' not found. Available ranks: ",
        paste(available_ranks, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    result <- phyloseq::tax_glom(
      physeq,
      taxrank = rank
    )

    phyloseq::taxa_names(result) <-
      as.character(phyloseq::tax_table(result)[, rank])

    if (isTRUE(verbose)) {
      message(
        "Analysis unit: taxonomic aggregation at ",
        rank,
        " level."
      )
    }

  } else if (method == "postcluster") {

    if (!postcluster_available) {
      stop(
        "Postcluster analysis was selected, but no postclustering was ",
        "performed during the bioinformatics pipeline. ",
        "Choose analysis_unit = 'asv' or 'taxonomy'.",
        call. = FALSE
      )
    }

    result <- tax_glom_postcluster_species(
      physeq,
      verbose = verbose
    )

    if (isTRUE(verbose)) {
      message(
        "Bioinformatics postclustering was performed at ",
        postcluster_threshold,
        "% identity."
      )
      message(
        "Analysis unit: postcluster-aware approximate species units."
      )
    }
  }

  if (isTRUE(verbose)) {
    message(
      "Taxa before: ", n_before,
      " | taxa after: ", phyloseq::ntaxa(result)
    )
  }

  result
}

remove_unresolved_taxa <- function(phyloseq, marker = NULL, verbose = TRUE) {
  .assert_phyloseq(phyloseq)

  # Backward compatibility with older scripts that set marker globally.
  if (is.null(marker)) {
    if (exists("marker", envir = parent.frame(), inherits = TRUE)) {
      marker <- get("marker", envir = parent.frame(), inherits = TRUE)
      warning(
        "Using global 'marker' is deprecated. Pass marker explicitly, e.g. ",
        "remove_unresolved_taxa(ps, marker = 'ITS2').",
        call. = FALSE
      )
    } else {
      stop("marker must be supplied: 'ITS2', 'COI', or '16S'.", call. = FALSE)
    }
  }

  marker <- toupper(as.character(marker)[1L])

  if (!marker %in% c("ITS2", "COI", "16S")) {
    stop(
      "Unsupported marker: ", marker,
      ". Use 'ITS2', 'COI', or '16S'.",
      call. = FALSE
    )
  }

  tax <- as.data.frame(
    phyloseq::tax_table(phyloseq),
    stringsAsFactors = FALSE
  )

  required <- switch(
    marker,
    ITS2 = c("kingdom", "phylum", "species"),
    COI = c("phylum", "class"),
    `16S` = c("kingdom", "phylum", "order", "genus")
  )

  missing_cols <- setdiff(required, names(tax))

  if (length(missing_cols) > 0L) {
    stop(
      "Taxonomy table lacks required rank(s) for ", marker, ": ",
      paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }

  # Store each removal criterion separately.
  reasons <- list()

  if (marker == "ITS2") {

    reasons$wrong_phylum <-
      is.na(tax$phylum) |
      tax$phylum != "p:Streptophyta"

    reasons$unresolved_taxon <-
      !is.na(tax$species) &
      tax$species %in% c(
        "p:Streptophyta_spc_spc_spc_spc",
        "k:Viridiplantae_spc_spc_spc_spc_spc"
      )

    reasons$missing_kingdom <-
      is.na(tax$kingdom) |
      tax$kingdom == ""

  } else if (marker == "COI") {

    reasons$unresolved_phylum <-
      !is.na(tax$phylum) &
      tax$phylum %in% "k:Animalia_spc"

    reasons$unresolved_class <-
      !is.na(tax$class) &
      tax$class %in% c(
        "p:Insecta_spc",
        "p:Arachnida_spc",
        "p:Gastropoda_spc"
      )

  } else if (marker == "16S") {

    reasons$non_bacterial <-
      is.na(tax$kingdom) |
      tax$kingdom != "d:Bacteria"

    reasons$chloroplast <-
      !is.na(tax$phylum) &
      tax$phylum %in% "p:Cyanobacteria/Chloroplast"

    reasons$unresolved_genus <-
      !is.na(tax$genus) &
      tax$genus %in% "d:Bacteria_spc_spc_spc_spc"

    reasons$mitochondria <-
      !is.na(tax$order) &
      tax$order %in% "f:Mitochondria"

    reasons$missing_kingdom <-
      is.na(tax$kingdom) |
      tax$kingdom == ""
  }

  # Convert possible NAs to FALSE.
  reasons <- lapply(
    reasons,
    function(x) {
      x[is.na(x)] <- FALSE
      x
    }
  )

  reason_matrix <- do.call(cbind, reasons)

  remove <- rowSums(reason_matrix) > 0
  keep <- !remove

  n_total <- nrow(tax)
  n_removed <- sum(remove)
  n_kept <- sum(keep)

  if (isTRUE(verbose)) {

    cat("\nUnresolved-taxon filtering:", marker, "\n")
    cat("  Taxa before filtering:", n_total, "\n")
    cat("  Taxa removed:         ", n_removed, "\n")
    cat("  Taxa retained:        ", n_kept, "\n")

    cat("\nRemoval reasons:\n")

    reason_counts <- sort(
      colSums(reason_matrix),
      decreasing = TRUE
    )

    for (nm in names(reason_counts)) {
      cat(
        "  ",
        format(nm, width = max(nchar(names(reason_counts)))),
        ": ",
        reason_counts[[nm]],
        "\n",
        sep = ""
      )
    }

    if (n_removed > 0L) {

      n_multiple <- sum(rowSums(reason_matrix) > 1)

      if (n_multiple > 0L) {
        cat(
          "\n  Note: ",
          n_multiple,
          " taxa matched more than one removal criterion.\n",
          sep = ""
        )
      }
    }

    cat("\n")
  }

  if (!any(keep)) {
    warning(
      "All taxa were removed by unresolved-taxon filtering.",
      call. = FALSE
    )
  }

  phyloseq::prune_taxa(keep, phyloseq)
}

propagate_incomplete_taxonomy <- function(phyloseq) {
  .assert_phyloseq(phyloseq)

  tax <- as.data.frame(phyloseq::tax_table(phyloseq), stringsAsFactors = FALSE)
  tax[is.na(tax)] <- ""
  taxranks <- names(tax)

  if (length(taxranks) < 2L) return(phyloseq)

  for (i in 2:length(taxranks)) {
    empty_current <- tax[[i]] == ""
    previous_known <- tax[[i - 1L]] != ""
    fill_rows <- empty_current & previous_known

    if (any(fill_rows)) {
      tax[[i]][fill_rows] <- paste0(tax[[i - 1L]][fill_rows], "_spc")
    }
  }

  phyloseq::tax_table(phyloseq) <- phyloseq::tax_table(as.matrix(tax))
  phyloseq
}


replace_tax_prefixes <- function(phyloseq) {
  .assert_phyloseq(phyloseq)

  tax <- as(phyloseq::tax_table(phyloseq), "matrix")
  tax <- gsub("^\\w:", "", tax)
  tax <- gsub("_spc_.*", "_spc", tax)
  tax <- gsub("_", " ", tax)
  tax <- gsub(";$", "", tax)

  phyloseq::tax_table(phyloseq) <- phyloseq::tax_table(tax)
  phyloseq
}

# -----------------------------------------------------------------------------
# Rarefaction
# -----------------------------------------------------------------------------

rarefaction_subsample <- function(x, sample.size, replace = FALSE) {
  if (!is.numeric(x) || anyNA(x) || any(x < 0) || any(x %% 1 != 0)) {
    stop("x must contain non-negative integer counts without NA values.", call. = FALSE)
  }
  if (!is.numeric(sample.size) || length(sample.size) != 1L ||
      !is.finite(sample.size) || sample.size < 0 || sample.size %% 1 != 0) {
    stop("sample.size must be a non-negative integer.", call. = FALSE)
  }

  sample.size <- as.integer(sample.size)
  rarvec <- numeric(length(x))
  total <- sum(x)

  if (total <= 0 || sample.size == 0L) return(rarvec)
  if (!replace && sample.size > total) {
    stop(
      "sample.size (", sample.size, ") exceeds total counts (", total,
      ") when replace = FALSE.",
      call. = FALSE
    )
  }

  if (replace) {
    # Multinomial sampling avoids constructing a vector containing every read.
    return(as.numeric(stats::rmultinom(1L, size = sample.size, prob = x)))
  }

  # Sequential multivariate-hypergeometric sampling. This has memory use O(taxa)
  # rather than O(number of reads), unlike expanding all observations first.
  remaining_total <- total
  remaining_draws <- sample.size

  if (length(x) == 1L) {
    rarvec[1L] <- sample.size
    return(rarvec)
  }

  for (i in seq_len(length(x) - 1L)) {
    if (remaining_draws == 0L) break
    if (x[i] == 0L) {
      remaining_total <- remaining_total - x[i]
      next
    }

    draw_i <- stats::rhyper(
      1L,
      m = x[i],
      n = remaining_total - x[i],
      k = remaining_draws
    )
    rarvec[i] <- draw_i
    remaining_draws <- remaining_draws - draw_i
    remaining_total <- remaining_total - x[i]
  }
  rarvec[length(x)] <- remaining_draws
  rarvec
}


# -----------------------------------------------------------------------------
# GBIF occurrence helpers
# -----------------------------------------------------------------------------

prepare_species_occurences <- function(
    gbif_file = "/Users/ra39huv/TMP/Data_processing/_resources/GBIF_20230724.csv") {
  if (!file.exists(gbif_file)) {
    stop("GBIF occurrence file not found: ", gbif_file, call. = FALSE)
  }

  gbif_records <- utils::read.table(
    gbif_file, sep = "\t", row.names = 1, fill = TRUE,
    stringsAsFactors = FALSE
  )

  if (ncol(gbif_records) < 3L) {
    stop("GBIF occurrence file must contain at least three data columns.", call. = FALSE)
  }
  gbif_records <- gbif_records[, seq_len(3L), drop = FALSE]
  names(gbif_records) <- c("value", "SciName", "Country")

  gbif_records <- gbif_records[
    gbif_records$Country != "countryCode" &
      !is.na(gbif_records$SciName) & gbif_records$SciName != "",
    , drop = FALSE
  ]

  gbif_wide <- stats::reshape(
    gbif_records,
    idvar = "SciName",
    timevar = "Country",
    direction = "wide"
  )
  gbif_wide <- gbif_wide[!is.na(gbif_wide$SciName), , drop = FALSE]
  rownames(gbif_wide) <- gbif_wide$SciName
  gbif_wide
}


get_species_occurences_wide <- function(phyloseq, gbif_wide = NULL, gbif_file = NULL) {
  .assert_phyloseq(phyloseq)

  if (is.null(gbif_wide)) {
    gbif_wide <- if (is.null(gbif_file)) {
      prepare_species_occurences()
    } else {
      prepare_species_occurences(gbif_file)
    }
  }

  otu <- .otu_taxa_rows(phyloseq)
  abund <- sort(rowSums(otu, na.rm = TRUE), decreasing = TRUE)
  species <- names(abund)

  matches <- gbif_wide[match(species, rownames(gbif_wide)), , drop = FALSE]
  out <- data.frame(
    Species = species,
    Abundance = as.numeric(abund),
    GBIF.match = rownames(matches),
    matches,
    row.names = NULL,
    check.names = FALSE
  )
  out
}


get_country_codes <- function(
    countrycode_file = "/Users/ra39huv/TMP/Data_processing/_resources/GBIF_countrycodes.csv") {
  if (!file.exists(countrycode_file)) {
    stop("GBIF country-code file not found: ", countrycode_file, call. = FALSE)
  }

  cc <- utils::read.table(
    countrycode_file, sep = ",", fill = TRUE, header = TRUE,
    row.names = 2, stringsAsFactors = FALSE
  )
  rownames(cc)[rownames(cc) == "Na"] <- "NA"
  cc
}


get_species_occurences_long <- function(phyloseq, notInGBIF = TRUE,
                                        gbif_wide = NULL, gbif_file = NULL,
                                        country_codes = NULL,
                                        countrycode_file = NULL) {
  .assert_phyloseq(phyloseq)

  taxa_matches <- get_species_occurences_wide(
    phyloseq, gbif_wide = gbif_wide, gbif_file = gbif_file
  )

  occurrence_cols <- grep("^value\\.", names(taxa_matches), value = TRUE)
  if (length(occurrence_cols) == 0L) {
    stop("No country occurrence columns ('value.*') found in GBIF data.", call. = FALSE)
  }

  long_list <- lapply(occurrence_cols, function(col) {
    data.frame(
      Species = taxa_matches$Species,
      Abundance = taxa_matches$Abundance,
      GBIF.match = taxa_matches$GBIF.match,
      Country = sub("^value\\.", "", col),
      value = taxa_matches[[col]],
      stringsAsFactors = FALSE
    )
  })
  taxa_matches_long <- do.call(rbind, long_list)
  taxa_matches_long <- taxa_matches_long[!is.na(taxa_matches_long$value), , drop = FALSE]

  if (is.null(country_codes)) {
    country_codes <- if (is.null(countrycode_file)) {
      get_country_codes()
    } else {
      get_country_codes(countrycode_file)
    }
  }

  wanted_cc_cols <- intersect(c("name", "region", "sub.region"), names(country_codes))
  cc_match <- country_codes[match(taxa_matches_long$Country, rownames(country_codes)),
                            wanted_cc_cols, drop = FALSE]
  taxa_matches_long2 <- cbind(taxa_matches_long, cc_match)

  if ("region" %in% names(taxa_matches_long2)) {
    taxa_matches_long2 <- taxa_matches_long2[order(taxa_matches_long2$region), , drop = FALSE]
  }

  if (isTRUE(notInGBIF)) {
    matched_species <- unique(taxa_matches_long2$Species)
    all_taxa <- phyloseq::taxa_names(phyloseq)
    missing_taxa <- setdiff(all_taxa, matched_species)

    if (length(missing_taxa) > 0L) {
      tax <- as.data.frame(phyloseq::tax_table(phyloseq), stringsAsFactors = FALSE)
      otu <- .otu_taxa_rows(phyloseq)
      species_labels <- if ("species" %in% names(tax)) {
        as.character(tax[missing_taxa, "species"])
      } else {
        missing_taxa
      }

      not_in_gbif <- data.frame(
        Species = species_labels,
        Abundance = as.numeric(rowSums(otu[missing_taxa, , drop = FALSE], na.rm = TRUE)),
        GBIF.match = NA_character_,
        Country = NA_character_,
        value = 10,
        stringsAsFactors = FALSE
      )
      for (nm in wanted_cc_cols) not_in_gbif[[nm]] <- NA_character_

      taxa_matches_long2 <- rbind(taxa_matches_long2, not_in_gbif)
    }
  }

  rownames(taxa_matches_long2) <- NULL
  taxa_matches_long2
}


# Correctly spelled aliases while retaining the historical function names.
prepare_species_occurrences <- prepare_species_occurences
get_species_occurrences_wide <- get_species_occurences_wide
get_species_occurrences_long <- get_species_occurences_long


# -----------------------------------------------------------------------------
# Sample labelling
# -----------------------------------------------------------------------------

label_low_throughput <- function(phyloseq,
                                 threshold,
                                 verbose = TRUE,
                                 preview_n = 10L) {
  .assert_phyloseq(phyloseq)
  .validate_positive_scalar(threshold, "threshold", allow_zero = TRUE)

  reads <- phyloseq::sample_sums(phyloseq)
  low <- reads < threshold

  n_total <- length(reads)
  n_low <- sum(low)

  flagged <- data.frame(
    sample = phyloseq::sample_names(phyloseq)[low],
    reads = as.numeric(reads[low]),
    threshold = rep(threshold, sum(low)),
    stringsAsFactors = FALSE
  )

  if (nrow(flagged) > 0L) {
    flagged <- flagged[order(flagged$reads), , drop = FALSE]
    rownames(flagged) <- NULL
  }

  if (isTRUE(verbose)) {
    cat("\nLow-throughput sample check\n")
    cat("  Threshold:       <", threshold, " reads\n")
    cat("  Samples checked: ", n_total, "\n")
    cat(
      "  Samples flagged: ",
      n_low,
      " (",
      round(100 * n_low / n_total, 1),
      "%)\n",
      sep = ""
    )

    if (n_low > 0L) {
      cat(
        "  Flagged range:   ",
        min(flagged$reads),
        " - ",
        max(flagged$reads),
        " reads\n"
      )

      if (preview_n > 0L) {
        cat(
          "\nLowest-throughput samples (first ",
          min(preview_n, n_low),
          " of ",
          n_low,
          "):\n",
          sep = ""
        )

        print(
          utils::head(flagged, preview_n),
          row.names = FALSE
        )
      }
    } else {
      cat("  No samples below threshold.\n")
    }

    cat("\n")
  }

  if (n_low > 0L) {
    names_new <- phyloseq::sample_names(phyloseq)
    names_new[low] <- paste0(names_new[low], "|LT", threshold)
    phyloseq::sample_names(phyloseq) <- names_new
  }

  structure(
    list(
      data = phyloseq,
      low_throughput = flagged
    ),
    class = c("low_throughput_result", "list")
  )
}

label_sample_by_host <- function(phyloseq, hostcol, projcol = NULL, idcol = NULL,
                                 seed = NULL) {
  .assert_phyloseq(phyloseq)
  md <- as.data.frame(phyloseq::sample_data(phyloseq), stringsAsFactors = FALSE)

  if (!hostcol %in% names(md)) {
    stop("hostcol '", hostcol, "' not found in sample_data.", call. = FALSE)
  }

  hostlist <- as.character(md[[hostcol]])
  projlist <- if (is.null(projcol) || identical(projcol, "")) {
    rep("", nrow(md))
  } else {
    if (!projcol %in% names(md)) {
      stop("projcol '", projcol, "' not found in sample_data.", call. = FALSE)
    }
    as.character(md[[projcol]])
  }

  if (is.null(idcol) || identical(idcol, "")) {
    if (!is.null(seed)) set.seed(seed)
    if (nrow(md) > 99999L) {
      stop("Cannot generate unique five-digit IDs for more than 99,999 samples.", call. = FALSE)
    }
    ids <- sample.int(99999L, nrow(md), replace = FALSE)
  } else if (length(idcol) == 1L && is.character(idcol) && idcol %in% names(md)) {
    ids <- md[[idcol]]
  } else if (length(idcol) == nrow(md)) {
    ids <- idcol
  } else {
    stop(
      "idcol must be NULL/'', a sample_data column name, or a vector with one ID per sample.",
      call. = FALSE
    )
  }

  id_text <- if (is.numeric(ids) && all(!is.na(ids)) && all(ids %% 1 == 0)) {
    sprintf("%05d", as.integer(ids))
  } else {
    as.character(ids)
  }

  fields <- if (all(projlist == "" | is.na(projlist))) {
    paste(hostlist, id_text, sep = "|")
  } else {
    paste(hostlist, projlist, id_text, sep = "|")
  }

  if (anyDuplicated(fields)) {
    stop("Generated sample labels are not unique.", call. = FALSE)
  }

  phyloseq::sample_names(phyloseq) <- fields
  phyloseq
}


# -----------------------------------------------------------------------------
# Contamination diagnostics
# -----------------------------------------------------------------------------

id_cont_lh <- function(phyloseq, distri = TRUE, neg = character(), pos = character(),
                       throughput_threshold = 2500, plot = TRUE) {
  .assert_phyloseq(phyloseq)
  .validate_positive_scalar(throughput_threshold, "throughput_threshold", allow_zero = TRUE)

  all_samples <- phyloseq::sample_names(phyloseq)
  neg <- intersect(as.character(neg), all_samples)
  pos <- intersect(as.character(pos), all_samples)

  sample_names_main <- setdiff(all_samples, neg)
  if (length(sample_names_main) == 0L) {
    stop("No non-negative-control samples remain.", call. = FALSE)
  }

  sam_samples <- phyloseq::prune_samples(sample_names_main, phyloseq)
  sam_otu <- .otu_taxa_rows(sam_samples)
  sam_rel <- sweep(sam_otu, 2L, colSums(sam_otu), "/")

  d1 <- rowMeans(sam_otu > 0, na.rm = TRUE)
  d2 <- rowSums(sam_rel, na.rm = TRUE)
  d3 <- rowMeans(sam_rel, na.rm = TRUE)

  result <- list(
    sample_prevalence = d1,
    sample_relative_sum = d2,
    sample_relative_mean = d3,
    negative_controls = neg,
    positive_controls = pos
  )

  if (length(neg) > 0L) {
    neg_samples <- phyloseq::prune_samples(neg, phyloseq)
    neg_otu <- .otu_taxa_rows(neg_samples)
    neg_totals <- colSums(neg_otu)
    neg_rel <- sweep(neg_otu, 2L, neg_totals, "/")

    d1n <- rowMeans(neg_otu > 0, na.rm = TRUE)
    d2n <- rowSums(neg_rel, na.rm = TRUE)
    d3n <- rowMeans(neg_rel, na.rm = TRUE)

    result$negative_prevalence <- d1n
    result$negative_relative_sum <- d2n
    result$negative_relative_mean <- d3n
    result$sample_minus_negative_relative_mean <- d3 - d3n
  }

  sample_totals <- colSums(sam_otu)
  high <- sample_totals > throughput_threshold
  low <- sample_totals < throughput_threshold

  if (any(high)) result$high_throughput_taxon_mean <- rowMeans(sam_rel[, high, drop = FALSE])
  if (any(low)) result$low_throughput_taxon_mean <- rowMeans(sam_rel[, low, drop = FALSE])

  if (isTRUE(plot)) {
    graphics::plot(
      d2, d1,
      xlab = "Summed relative abundance across samples",
      ylab = "Sample prevalence",
      main = "Potential widespread low-abundance contaminants"
    )

    if (length(neg) > 0L) {
      graphics::plot(
        d3, d3n,
        xlab = "Mean relative abundance in samples",
        ylab = "Mean relative abundance in negative controls",
        main = "Samples vs negative controls"
      )
    }

    if (any(high) && any(low)) {
      graphics::plot(
        result$low_throughput_taxon_mean,
        result$high_throughput_taxon_mean,
        xlab = "Mean relative abundance - low throughput",
        ylab = "Mean relative abundance - high throughput",
        main = "Low vs high throughput samples"
      )
    }
  }

  invisible(result)
}


# -----------------------------------------------------------------------------
# Trait mapping
# -----------------------------------------------------------------------------

map_interactions_trait <- function(phyloseq, traittable, traitcols = NULL, speccol) {
  .assert_phyloseq(phyloseq)

  if (missing(speccol) || length(speccol) != 1L) {
    stop("speccol must identify the species/taxon column in traittable.", call. = FALSE)
  }
  if (is.character(speccol)) {
    if (!speccol %in% names(traittable)) {
      stop("speccol '", speccol, "' not found in traittable.", call. = FALSE)
    }
    species_values <- traittable[[speccol]]
  } else {
    species_values <- traittable[[speccol]]
  }

  if (is.null(traitcols) || length(traitcols) == 0L || identical(traitcols, "")) {
    traitcols <- setdiff(names(traittable), if (is.character(speccol)) speccol else names(traittable)[speccol])
  }
  traits <- traittable[, traitcols, drop = FALSE]

  non_numeric <- !vapply(traits, is.numeric, logical(1))
  if (any(non_numeric)) {
    stop(
      "All selected trait columns must be numeric. Non-numeric: ",
      paste(names(traits)[non_numeric], collapse = ", "), ".",
      call. = FALSE
    )
  }

  otu <- .otu_taxa_rows(phyloseq)
  taxa <- rownames(otu)
  samples <- colnames(otu)

  trait_matrix <- matrix(
    NA_real_, nrow = length(taxa), ncol = ncol(traits),
    dimnames = list(taxa, names(traits))
  )
  results_nomatch <- character()
  results_multiple <- character()

  for (i in seq_along(taxa)) {
    matches <- which(!is.na(species_values) & species_values == taxa[i])
    if (length(matches) == 1L) {
      trait_matrix[i, ] <- as.numeric(traits[matches, , drop = TRUE])
    } else if (length(matches) > 1L) {
      results_multiple <- c(results_multiple, taxa[i])
    } else {
      results_nomatch <- c(results_nomatch, taxa[i])
    }
  }

  mapping <- matrix(
    NA_real_, nrow = length(samples), ncol = ncol(trait_matrix),
    dimnames = list(samples, colnames(trait_matrix))
  )
  uncertainty <- mapping
  uncertainty_weighted <- mapping

  for (j in seq_along(samples)) {
    abundance <- otu[, j]
    weighted <- trait_matrix * abundance

    mapping[j, ] <- colSums(weighted, na.rm = TRUE)
    uncertainty[j, ] <- colMeans(is.na(weighted))

    known <- !is.na(trait_matrix)
    uncertainty_weighted[j, ] <- colSums(known * abundance, na.rm = TRUE)
  }

  list(
    mapping = as.data.frame(mapping),
    uncertainity = as.data.frame(uncertainty),
    uncertainity_weighted = as.data.frame(uncertainty_weighted),
    nomatch = results_nomatch,
    multiple_matches = results_multiple
  )
}

# -----------------------------------------------------------------------------
# Handle control samples
# -----------------------------------------------------------------------------

inspect_controls <- function(
    physeq,
    control_column = "Type",
    negative = "negative",
    positive = "positive",
    rank = NULL,
    output_dir = "plots",
    top_n_negative = 20L,
    positive_display_min = 0.0001,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  md <- as.data.frame(
    phyloseq::sample_data(physeq),
    stringsAsFactors = FALSE
  )

  if (!control_column %in% names(md)) {
    if (isTRUE(verbose)) {
      message(
        "Control inspection skipped: metadata column '",
        control_column,
        "' not found."
      )
    }

    return(invisible(list(
      negative_samples = character(0),
      positive_samples = character(0),
      biological_samples = phyloseq::sample_names(physeq),
      negative_enrichment = NULL
    )))
  }

  control_type <- as.character(md[[control_column]])

  negative_samples <- rownames(md)[
    !is.na(control_type) & control_type == negative
  ]

  positive_samples <- rownames(md)[
    !is.na(control_type) & control_type == positive
  ]

  biological_samples <- rownames(md)[
    is.na(control_type) |
      !(control_type %in% c(negative, positive))
  ]

  has_negative <- length(negative_samples) > 0L
  has_positive <- length(positive_samples) > 0L

  if (!has_negative && !has_positive) {
    if (isTRUE(verbose)) {
      message(
        "No positive or negative controls detected in metadata column '",
        control_column,
        "'. Control inspection skipped."
      )
    }

    return(invisible(list(
      negative_samples = character(0),
      positive_samples = character(0),
      biological_samples = biological_samples,
      negative_enrichment = NULL
    )))
  }

  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  tax <- as.data.frame(
    phyloseq::tax_table(physeq),
    stringsAsFactors = FALSE
  )

  if (is.null(rank) || !rank %in% names(tax)) {
    rank <- tail(names(tax), 1L)
  }


  # Common relative-abundance scale used in control plots.
  abundance_breaks <- c(
    0,
    0.00001,  # 0.001%
    0.0001,   # 0.01%
    0.001,    # 0.1%
    0.01,     # 1%
    0.1,      # 10%
    1         # 100%
  )

  abundance_labels <- c(
    "0",
    "0.001%",
    "0.01%",
    "0.1%",
    "1%",
    "10%",
    "100%"
  )


  # ========================================================
  # NEGATIVE CONTROLS
  # ========================================================

  negative_enrichment <- NULL

  if (has_negative) {

    # ------------------------------------------------------
    # Sequencing depth: negatives versus biological samples

    compare_samples <- c(
      negative_samples,
      biological_samples
    )

    depth <- data.frame(
      sample = compare_samples,
      reads = as.numeric(
        phyloseq::sample_sums(physeq)[compare_samples]
      ),
      sample_group = c(
        rep("negative", length(negative_samples)),
        rep("biological", length(biological_samples))
      ),
      stringsAsFactors = FALSE
    )

    depth_plot <- depth[
      depth$reads > 0,
      ,
      drop = FALSE
    ]

    if (nrow(depth_plot) > 0L) {

      p.depth <- ggplot2::ggplot(
        depth_plot,
        ggplot2::aes(
          x = sample_group,
          y = reads
        )
      ) +
        ggplot2::geom_violin(
          trim = FALSE,
          scale = "width"
        ) +
        ggplot2::geom_boxplot(
          width = 0.12,
          outlier.shape = NA
        ) +
        ggplot2::geom_jitter(
          width = 0.10,
          height = 0,
          size = 2
        ) +
        ggplot2::scale_y_log10() +
        ggplot2::labs(
          x = NULL,
          y = "Total reads per sample",
          title = "Sequencing depth: negative vs biological samples"
        ) +
        ggplot2::theme_classic()

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          "controls_negative_sequencing_depth.pdf"
        ),
        plot = p.depth,
        width = 5,
        height = 5
      )
    }


    # ------------------------------------------------------
    # Taxa enriched in negative controls

    ps.compare <- phyloseq::prune_samples(
      compare_samples,
      physeq
    )

    ps.compare.rel <- phyloseq::transform_sample_counts(
      ps.compare,
      function(x) {
        if (sum(x) == 0) x else x / sum(x)
      }
    )

    otu <- as(
      phyloseq::otu_table(ps.compare.rel),
      "matrix"
    )

    if (!phyloseq::taxa_are_rows(ps.compare.rel)) {
      otu <- t(otu)
    }

    neg_idx <- colnames(otu) %in% negative_samples
    bio_idx <- colnames(otu) %in% biological_samples

    if (any(neg_idx) && any(bio_idx)) {

      mean_negative <- rowMeans(
        otu[, neg_idx, drop = FALSE],
        na.rm = TRUE
      )

      mean_biological <- rowMeans(
        otu[, bio_idx, drop = FALSE],
        na.rm = TRUE
      )

      prevalence_negative <- rowMeans(
        otu[, neg_idx, drop = FALSE] > 0,
        na.rm = TRUE
      )

      prevalence_biological <- rowMeans(
        otu[, bio_idx, drop = FALSE] > 0,
        na.rm = TRUE
      )

      nonzero <- otu[otu > 0]

      pseudocount <- if (length(nonzero) > 0L) {
        min(nonzero, na.rm = TRUE) / 2
      } else {
        1e-12
      }

      log2_enrichment <- log2(
        (mean_negative + pseudocount) /
          (mean_biological + pseudocount)
      )

      negative_enrichment <- data.frame(
        taxon = rownames(otu),
        mean_negative = mean_negative,
        mean_biological = mean_biological,
        prevalence_negative = prevalence_negative,
        prevalence_biological = prevalence_biological,
        log2_enrichment = log2_enrichment,
        stringsAsFactors = FALSE
      )

      # Only taxa actually detected in negative controls.
      negative_enrichment <- negative_enrichment[
        negative_enrichment$prevalence_negative > 0,
        ,
        drop = FALSE
      ]

      negative_enrichment <- negative_enrichment[
        order(
          negative_enrichment$log2_enrichment,
          decreasing = TRUE
        ),
        ,
        drop = FALSE
      ]

      # Show the taxa most enriched in negative controls.
      plot_enrichment <- negative_enrichment[
        negative_enrichment$log2_enrichment > 0,
        ,
        drop = FALSE
      ]

      plot_enrichment <- utils::head(
        plot_enrichment,
        top_n_negative
      )

      if (nrow(plot_enrichment) > 0L) {

        # Same taxon order in both panels.
        plot_enrichment$taxon <- factor(
          plot_enrichment$taxon,
          levels = rev(plot_enrichment$taxon)
        )


        ## Panel 1: enrichment relative to biological samples

        p.neg.fold <- ggplot2::ggplot(
          plot_enrichment,
          ggplot2::aes(
            x = taxon,
            y = log2_enrichment
          )
        ) +
          ggplot2::geom_col() +
          ggplot2::coord_flip() +
          ggplot2::geom_hline(
            yintercept = 0,
            linetype = "dashed"
          ) +
          ggplot2::labs(
            x = NULL,
            y = "log2 enrichment",
            title = "Overrepresentation in negatives"
          ) +
          ggplot2::theme_classic()


        ## Panel 2: mean relative abundance in negatives and biological samples

        abundance_long <- rbind(
          data.frame(
            taxon = as.character(plot_enrichment$taxon),
            sample_group = "negative",
            mean_abundance = plot_enrichment$mean_negative
          ),
          data.frame(
            taxon = as.character(plot_enrichment$taxon),
            sample_group = "biological",
            mean_abundance = plot_enrichment$mean_biological
          )
        )

        abundance_long$taxon <- factor(
          abundance_long$taxon,
          levels = levels(plot_enrichment$taxon)
        )

        p.neg.abundance <- ggplot2::ggplot(
          abundance_long,
          ggplot2::aes(
            x = taxon,
            y = mean_abundance,
            fill = sample_group
          )
        ) +
          ggplot2::geom_col(
            position = "dodge"
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_y_continuous(
            trans = scales::pseudo_log_trans(
              base = 10,
              sigma = 0.000005
            ),
            limits = c(0, 1),
            breaks = abundance_breaks,
            labels = abundance_labels
          ) +
          ggplot2::geom_hline(
            yintercept = c(
              0.00001,  # 0.001%
              0.0001,   # 0.01%
              0.001,    # 0.1%
              0.01      # 1%
            ),
            linetype = "dashed"
          ) +
          ggplot2::labs(
            x = NULL,
            y = "Mean relative abundance",
            fill = NULL,
            title = "Mean abundance"
          ) +
          ggplot2::theme_classic() +
          ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
          )


        # Combine panels when patchwork is available.
        if (requireNamespace("patchwork", quietly = TRUE)) {

          p.neg <- p.neg.fold +
            p.neg.abundance +
            patchwork::plot_layout(
              widths = c(1.15, 1.35)
            )

          ggplot2::ggsave(
            filename = file.path(
              output_dir,
              "controls_negative_enriched_taxa.pdf"
            ),
            plot = p.neg,
            width = 11,
            height = max(
              5,
              0.25 * nrow(plot_enrichment) + 2
            )
          )

        } else {

          ggplot2::ggsave(
            filename = file.path(
              output_dir,
              "controls_negative_enriched_taxa.pdf"
            ),
            plot = p.neg.fold,
            width = 6,
            height = max(
              5,
              0.25 * nrow(plot_enrichment) + 2
            )
          )

          ggplot2::ggsave(
            filename = file.path(
              output_dir,
              "controls_negative_enriched_taxa_abundance.pdf"
            ),
            plot = p.neg.abundance,
            width = 6,
            height = max(
              5,
              0.25 * nrow(plot_enrichment) + 2
            )
          )

          if (isTRUE(verbose)) {
            message(
              "Package 'patchwork' not installed: negative-control panels ",
              "were written as separate PDFs."
            )
          }
        }
      }
    }
  }


  # ========================================================
  # POSITIVE CONTROLS
  # ========================================================

  if (has_positive) {

    # Positive controls are inspected on their own. Relative abundance makes
    # low-level unexpected taxa easier to compare between non-empty samples.

    ps.pos <- phyloseq::prune_samples(
      positive_samples,
      physeq
    )

    ps.pos.rel <- phyloseq::transform_sample_counts(
      ps.pos,
      function(x) {
        if (sum(x) == 0) x else x / sum(x)
      }
    )

    # Remove taxa absent from all positive controls.
    ps.pos.rel <- phyloseq::prune_taxa(
      phyloseq::taxa_sums(ps.pos.rel) > 0,
      ps.pos.rel
    )

    pos.melt <- phyloseq::psmelt(
      ps.pos.rel
    )

    # Keep taxa reaching the display threshold in at least one positive sample.
    taxa_to_plot <- unique(
      pos.melt[[rank]][
        pos.melt$Abundance >= positive_display_min
      ]
    )

    pos.melt <- pos.melt[
      pos.melt[[rank]] %in% taxa_to_plot &
        pos.melt$Abundance >= positive_display_min,
      ,
      drop = FALSE
    ]

    if (nrow(pos.melt) > 0L) {

      p.pos <- ggplot2::ggplot(
        pos.melt,
        ggplot2::aes(
          x = .data[[rank]],
          y = Abundance
        )
      ) +
        ggplot2::geom_boxplot(
          outlier.shape = NA
        ) +
        ggplot2::geom_jitter(
          width = 0.15,
          size = 1.5
        ) +
        ggplot2::scale_y_continuous(
          trans = scales::pseudo_log_trans(
            base = 10,
            sigma = 0.000005
          ),
          limits = c(0, 1),
          breaks = abundance_breaks,
          labels = abundance_labels
        ) +
        ggplot2::geom_hline(
          yintercept = c(
            0.00001,  # 0.001%
            0.0001,   # 0.01%
            0.001,    # 0.1%
            0.01      # 1%
          ),
          linetype = "dashed"
        ) +
        ggplot2::labs(
          x = NULL,
          y = "Relative abundance",
          title = "Taxonomic composition of positive controls"
        ) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 7
          )
        )

      n_taxa_plot <- length(
        unique(pos.melt[[rank]])
      )

      ggplot2::ggsave(
        filename = file.path(
          output_dir,
          "controls_positive_relative.pdf"
        ),
        plot = p.pos,
        width = min(
          11,
          max(
            7,
            0.18 * n_taxa_plot + 4
          )
        ),
        height = 5
      )
    }
  }


  # ========================================================
  # SUMMARY
  # ========================================================

  if (isTRUE(verbose)) {

    message("Control inspection:")

    message(
      "  Negative controls:  ",
      length(negative_samples)
    )

    message(
      "  Positive controls:  ",
      length(positive_samples)
    )

    message(
      "  Biological samples: ",
      length(biological_samples)
    )

    if (has_negative) {
      message(
        "  Negative-control plots: sequencing depth and enriched taxa."
      )
    }

    if (has_positive) {
      message(
        "  Positive-control plot: relative taxonomic composition."
      )
    }

    message(
      "  Plots written to: ",
      output_dir,
      "/"
    )
  }

  invisible(
    list(
      negative_samples = negative_samples,
      positive_samples = positive_samples,
      biological_samples = biological_samples,
      negative_enrichment = negative_enrichment
    )
  )
}

# ========================================================
# Hill number helper functions  
# ========================================================

inspect_sample_coverage <- function(
    physeq,
    output_file = "plots/sample_coverage_curves.pdf",
    endpoint_factor = 2,
    nboot = 20,
    conf = 0.95,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    stop(
      "Package 'iNEXT' is required for coverage-based diversity analysis.",
      call. = FALSE
    )
  }

  if (!is.numeric(endpoint_factor) ||
      length(endpoint_factor) != 1L ||
      !is.finite(endpoint_factor) ||
      endpoint_factor < 1) {
    stop(
      "endpoint_factor must be a numeric value >= 1.",
      call. = FALSE
    )
  }

  otu <- .otu_taxa_rows(physeq)
  otu <- otu[rowSums(otu, na.rm = TRUE) > 0, , drop = FALSE]
  otu <- otu[, colSums(otu, na.rm = TRUE) > 0, drop = FALSE]

  if (nrow(otu) == 0L || ncol(otu) == 0L) {
    stop(
      "No non-empty samples and taxa available for coverage analysis.",
      call. = FALSE
    )
  }

  # One abundance vector per sample.
  abundance_list <- lapply(seq_len(ncol(otu)), function(i) otu[, i])
  names(abundance_list) <- colnames(otu)

  sample_reads <- vapply(abundance_list, sum, numeric(1))


  ## Observed coverage
  info <- iNEXT::DataInfo(
    abundance_list,
    datatype = "abundance"
  )

  coverage_table <- data.frame(
    Sample = as.character(info$Assemblage),
    reads = as.numeric(info$n),
    observed_taxa = as.numeric(info$S.obs),
    observed_coverage = as.numeric(info$SC),
    stringsAsFactors = FALSE
  )


  ## Completeness curves
  # iNEXT evaluates coverage from rarefaction through extrapolation.
  # By default we allow extrapolation up to 2x the observed sample size.

  endpoint <- ceiling(
    endpoint_factor * max(sample_reads)
  )

  coverage_inext <- iNEXT::iNEXT(
    abundance_list,
    q = 0,
    datatype = "abundance",
    endpoint = endpoint,
    se = nboot > 0,
    conf = conf,
    nboot = nboot
  )


  ## Suggested common coverage
  # iNEXT's automatic coverage standardization gives the common completeness
  # level attainable across assemblages under its interpolation/extrapolation
  # framework.

  common_est <- iNEXT::estimateD(
    abundance_list,
    q = 0,
    datatype = "abundance",
    base = "coverage",
    level = NULL,
    nboot = 0
  )

  suggested_coverage <- min(
    as.numeric(common_est$SC),
    na.rm = TRUE
  )

  coverage_table$requires_extrapolation <-
    coverage_table$observed_coverage < suggested_coverage

  coverage_table <- coverage_table[
    order(coverage_table$observed_coverage),
    ,
    drop = FALSE
  ]

  rownames(coverage_table) <- NULL


  ## Diagnostic plot
  # Type 2 = sample completeness curve: coverage versus sample size.
  # Solid portions represent interpolation and dashed portions extrapolation.

  coverage_plot <- iNEXT::ggiNEXT(
    coverage_inext,
    type = 2,
    se = nboot > 0,
    color.var = "Assemblage"
  ) +
    ggplot2::geom_hline(
      yintercept = suggested_coverage,
      linetype = "dotted"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1)
    ) +
    ggplot2::labs(
      x = "Sequencing depth",
      y = "Sample coverage",
      title = "Sample completeness curves",
      subtitle = paste0(
        "Dotted line: suggested common coverage = ",
        round(100 * suggested_coverage, 1),
        "%"
      ),
      colour = "Sample"
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(
    output_file,
    plot = coverage_plot,
    width = 9,
    height = 6
  )


  ## Console summary
  if (isTRUE(verbose)) {

    n_interpolation <- sum(
      !coverage_table$requires_extrapolation,
      na.rm = TRUE
    )

    message("Sample coverage overview:")
    message(
      "  Samples: ", nrow(coverage_table),
      " | median observed coverage: ",
      round(
        100 * stats::median(
          coverage_table$observed_coverage,
          na.rm = TRUE
        ),
        1
      ),
      "%"
    )
    message(
      "  Minimum observed coverage: ",
      round(
        100 * min(
          coverage_table$observed_coverage,
          na.rm = TRUE
        ),
        1
      ),
      "%"
    )
    message(
      "  Suggested common coverage: ",
      round(100 * suggested_coverage, 1),
      "%"
    )
    message(
      "  Reached without extrapolation by ",
      n_interpolation, "/",
      nrow(coverage_table),
      " samples."
    )
    message(
      "  Inspect the completeness curves before accepting this value, ",
      "especially if a few poorly covered samples lower the common target."
    )
    message(
      "  Diagnostic plot written to: ",
      output_file
    )
  }

  invisible(
    list(
      table = coverage_table,
      suggested_coverage = suggested_coverage,
      inext = coverage_inext,
      plot = coverage_plot
    )
  )
}

estimate_hill_diversity <- function(
    physeq,
    coverage,
    q = c(0, 1, 2),
    nboot = 50,
    conf = 0.95,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    stop(
      "Package 'iNEXT' is required for coverage-based diversity analysis.",
      call. = FALSE
    )
  }

  if (!is.numeric(coverage) ||
      length(coverage) != 1L ||
      !is.finite(coverage) ||
      coverage <= 0 ||
      coverage > 1) {
    stop(
      "coverage must be a single value > 0 and <= 1.",
      call. = FALSE
    )
  }

  if (!is.numeric(q) || length(q) == 0L || any(!is.finite(q)) || any(q < 0)) {
    stop(
      "q must contain one or more non-negative Hill-number orders.",
      call. = FALSE
    )
  }

  otu <- .otu_taxa_rows(physeq)
  otu <- otu[rowSums(otu, na.rm = TRUE) > 0, , drop = FALSE]
  otu <- otu[, colSums(otu, na.rm = TRUE) > 0, drop = FALSE]

  if (nrow(otu) == 0L || ncol(otu) == 0L) {
    stop(
      "No non-empty samples and taxa available for diversity estimation.",
      call. = FALSE
    )
  }

  abundance_list <- lapply(
    seq_len(ncol(otu)),
    function(i) otu[, i]
  )
  names(abundance_list) <- colnames(otu)

  ## Coverage-standardized Hill numbers
  hill <- iNEXT::estimateD(
    abundance_list,
    q = q,
    datatype = "abundance",
    base = "coverage",
    level = coverage,
    nboot = nboot,
    conf = conf
  )

  # Standardize column names for easier downstream use.
  rename_if_present <- function(x, old, new) {
    if (old %in% names(x)) {
      names(x)[names(x) == old] <- new
    }
    x
  }

  hill <- rename_if_present(hill, "Assemblage", "Sample")
  hill <- rename_if_present(hill, "Order.q", "q")
  hill <- rename_if_present(hill, "qD", "diversity")
  hill <- rename_if_present(hill, "SC", "coverage")
  hill <- rename_if_present(hill, "Method", "method")
  hill <- rename_if_present(hill, "qD.LCL", "diversity_LCL")
  hill <- rename_if_present(hill, "qD.UCL", "diversity_UCL")
  hill <- rename_if_present(hill, "m", "sample_size")

  hill$Sample <- as.character(hill$Sample)

  ## Add biological/QC metadata
  md <- as.data.frame(
    phyloseq::sample_data(physeq),
    stringsAsFactors = FALSE
  )
  md$Sample <- phyloseq::sample_names(physeq)

  md$Sample <- rownames(md)

  match_idx <- match(hill$Sample, md$Sample)

  if (anyNA(match_idx)) {
    warning(
      "Metadata could not be matched for ",
      sum(is.na(match_idx)),
      " Hill-diversity rows. Example sample IDs: ",
      paste(utils::head(hill$Sample[is.na(match_idx)], 5), collapse = ", "),
      call. = FALSE
    )
  }

  metadata_cols <- setdiff(names(md), "Sample")

  for (nm in metadata_cols) {
    hill[[nm]] <- md[[nm]][match_idx]
  }
 
  if (isTRUE(verbose)) {
    message(
      "  Metadata matched for ",
      sum(!is.na(match_idx)), "/", length(match_idx),
      " diversity estimates."
    )
  }

  # Restore original sample order.
  sample_order <- phyloseq::sample_names(physeq)
  hill$.sample_order <- match(hill$Sample, sample_order)

  hill <- hill[
    order(hill$.sample_order, hill$q),
    ,
    drop = FALSE
  ]

  hill$.sample_order <- NULL
  rownames(hill) <- NULL

  ## Human-readable diversity labels
  hill$Hill_number <- ifelse(
    hill$q == 0,
    "q = 0 (richness)",
    ifelse(
      hill$q == 1,
      "q = 1 (exp. Shannon)",
      ifelse(
        hill$q == 2,
        "q = 2 (inv. Simpson)",
        paste0("q = ", hill$q)
      )
    )
  )

  if (isTRUE(verbose)) {
    message(
      "Coverage-standardized Hill diversity calculated for ",
      length(unique(hill$Sample)),
      " samples at ",
      round(100 * coverage, 1),
      "% coverage."
    )
    message(
      "  q = 0: richness | q = 1: exp. Shannon | ",
      "q = 2: inverse Simpson"
    )

    if ("method" %in% names(hill)) {
      method_counts <- table(hill$method)
      message(
        "  Estimation: ",
        paste(
          paste(names(method_counts), method_counts, sep = "="),
          collapse = ", "
        )
      )
    }
  }

  hill
}

plot_hill_diversity <- function(
    hill_table,
    group,
    color = NULL,
    output_file = "plots/sample_hill_diversity.pdf",
    verbose = TRUE) {

  if (!is.data.frame(hill_table)) {
    stop("hill_table must be a data.frame.", call. = FALSE)
  }

  required <- c("Sample", "q", "diversity", "Hill_number")
  missing_cols <- setdiff(required, names(hill_table))

  if (length(missing_cols) > 0L) {
    stop(
      "hill_table lacks required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!group %in% names(hill_table)) {
    stop(
      "Grouping variable '", group, "' not found.",
      call. = FALSE
    )
  }

  if (!is.null(color) && !color %in% names(hill_table)) {
    stop(
      "Colour variable '", color, "' not found.",
      call. = FALSE
    )
  }

  ## Remove rows lacking variables required for plotting
  keep <- !is.na(hill_table[[group]]) &
          !is.na(hill_table$diversity)

  if (!is.null(color)) {
    keep <- keep & !is.na(hill_table[[color]])
  }

  plot_df <- hill_table[keep, , drop = FALSE]

  if (nrow(plot_df) == 0L) {
    stop(
      "No complete observations available for Hill-diversity plotting.",
      call. = FALSE
    )
  }

  ## Keep q panels in numerical order
  hill_levels <- unique(
    plot_df$Hill_number[order(plot_df$q)]
  )

  plot_df$Hill_number <- factor(
    plot_df$Hill_number,
    levels = hill_levels
  )

  ## Plot dimensions
  n_groups <- dplyr::n_distinct(plot_df[[group]], na.rm = TRUE)
  plot_width <- min(16, max(6, 1.1 * n_groups + 3))
  plot_height <- 10

  ## Plot
  if (is.null(color)) {

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[group]],
        y = diversity
      )
    ) +
      ggplot2::geom_boxplot(
        outlier.shape = NA,
        width = 0.65
      ) +
      ggplot2::geom_jitter(
        width = 0.15,
        alpha = 0.5,
        size = 2
      )

  } else {

    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(
        x = .data[[group]],
        y = diversity,
        color = .data[[color]]
      )
    ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(group = interaction(
          .data[[group]],
          .data[[color]]
        )),
        outlier.shape = NA,
        position = ggplot2::position_dodge(width = 0.7),
        width = 0.6
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.12,
          dodge.width = 0.7
        ),
        alpha = 0.6,
        size = 2
      )
  }

  p <- p +
    ggplot2::facet_wrap(
      ~ Hill_number,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = group,
      y = "Effective number of taxa",
      color = color,
      title = "Coverage-standardized Hill diversity"
    )

  ggplot2::ggsave(
    output_file,
    plot = p,
    width = plot_width,
    height = plot_height
  )

  if (isTRUE(verbose)) {
    message(
      "Hill-diversity plot grouped by '", group, "'",
      if (!is.null(color)) paste0(" and coloured by '", color, "'") else "",
      " written to: ", output_file
    )
  }

  invisible(p)
}



#############################
## Network Plotting
make_bipartite_network <- function(
    physeq,
    group,
    rank = NULL,
    min_abundance = 0,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  md <- as.data.frame(
    phyloseq::sample_data(physeq),
    stringsAsFactors = FALSE
  )

  if (!group %in% names(md)) {
    stop(
      "Grouping variable '", group, "' not found in sample metadata.",
      call. = FALSE
    )
  }

  if (!is.numeric(min_abundance) ||
      length(min_abundance) != 1L ||
      !is.finite(min_abundance) ||
      min_abundance < 0 ||
      min_abundance > 1) {
    stop(
      "min_abundance must be between 0 and 1.",
      call. = FALSE
    )
  }

  ## Optionally aggregate taxonomy
  ps <- physeq

  if (!is.null(rank)) {

    ranks <- colnames(phyloseq::tax_table(ps))

    if (!rank %in% ranks) {
      stop(
        "Taxonomic rank '", rank, "' not found. Available ranks: ",
        paste(ranks, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    ps <- phyloseq::tax_glom(
      ps,
      taxrank = rank
    )

    phyloseq::taxa_names(ps) <- as.character(
      phyloseq::tax_table(ps)[, rank]
    )
  }

  ## Relative abundance
  ps <- phyloseq::transform_sample_counts(
    ps,
    function(x) {
      if (sum(x) == 0) x else x / sum(x)
    }
  )

  otu <- .otu_taxa_rows(ps)

  groups <- as.character(
    phyloseq::sample_data(ps)[[group]]
  )
  names(groups) <- phyloseq::sample_names(ps)

  valid <- !is.na(groups) & groups != ""
  groups <- groups[valid]
  otu <- otu[, names(groups), drop = FALSE]

  if (length(groups) == 0L) {
    stop(
      "No samples with non-missing values for '",
      group, "'.",
      call. = FALSE
    )
  }

  ## Mean relative abundance within each biological group
  group_levels <- unique(groups)

  network_matrix <- vapply(
    group_levels,
    function(g) {
      rowMeans(
        otu[, groups == g, drop = FALSE],
        na.rm = TRUE
      )
    },
    numeric(nrow(otu))
  )

  rownames(network_matrix) <- rownames(otu)
  colnames(network_matrix) <- group_levels

  ## Remove weak links
  network_matrix[
    network_matrix < min_abundance
  ] <- 0

  ## Remove disconnected taxa and groups
  network_matrix <- network_matrix[
    rowSums(network_matrix) > 0,
    colSums(network_matrix) > 0,
    drop = FALSE
  ]

  if (nrow(network_matrix) == 0L ||
      ncol(network_matrix) == 0L) {
    stop(
      "No interactions remain after network filtering.",
      call. = FALSE
    )
  }

  # bipartite expects lower level in rows and higher level in columns.
  network_matrix <- t(network_matrix)

  if (isTRUE(verbose)) {
    message(
      "Bipartite network: ",
      nrow(network_matrix), " ", group, " categories x ",
      ncol(network_matrix), " taxa | ",
      sum(network_matrix > 0), " links retained."
    )
  }

  network_matrix
}

make_bipartite_network <- function(
    physeq,
    group,
    rank = NULL,
    min_abundance = 0,
    verbose = TRUE) {

  .assert_phyloseq(physeq)

  md <- as.data.frame(
    phyloseq::sample_data(physeq),
    stringsAsFactors = FALSE
  )

  if (!group %in% names(md)) {
    stop(
      "Grouping variable '", group,
      "' not found in sample metadata.",
      call. = FALSE
    )
  }

  if (!is.numeric(min_abundance) ||
      length(min_abundance) != 1L ||
      !is.finite(min_abundance) ||
      min_abundance < 0 ||
      min_abundance > 1) {
    stop(
      "min_abundance must be between 0 and 1.",
      call. = FALSE
    )
  }

  ps <- physeq


  ## Aggregate taxonomy if requested
  if (!is.null(rank)) {

    tax_ranks <- colnames(
      phyloseq::tax_table(ps)
    )

    if (!rank %in% tax_ranks) {
      stop(
        "Taxonomic rank '", rank,
        "' not found. Available ranks: ",
        paste(tax_ranks, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    ps <- phyloseq::tax_glom(
      ps,
      taxrank = rank
    )

    rank_names <- as.character(
      phyloseq::tax_table(ps)[, rank]
    )

    missing_rank <- is.na(rank_names) | rank_names == ""

    if (any(missing_rank)) {
      rank_names[missing_rank] <- phyloseq::taxa_names(ps)[missing_rank]
    }

    phyloseq::taxa_names(ps) <- make.unique(rank_names)
  }


  ## Store taxonomy belonging to network taxa
  taxon_metadata <- as.data.frame(
    phyloseq::tax_table(ps),
    stringsAsFactors = FALSE
  )

  taxon_metadata$network_taxon <- phyloseq::taxa_names(ps)


  ## Convert to relative abundance
  ps <- phyloseq::transform_sample_counts(
    ps,
    function(x) {
      if (sum(x) == 0) x else x / sum(x)
    }
  )

  otu <- .otu_taxa_rows(ps)

  groups <- as.character(
    phyloseq::sample_data(ps)[[group]]
  )

  names(groups) <- phyloseq::sample_names(ps)

  valid <- !is.na(groups) & groups != ""

  groups <- groups[valid]
  otu <- otu[, names(groups), drop = FALSE]

  if (length(groups) == 0L) {
    stop(
      "No samples with non-missing values for '",
      group, "'.",
      call. = FALSE
    )
  }


  ## Mean relative abundance within groups
  group_levels <- unique(groups)

  network_matrix <- vapply(
    group_levels,
    function(g) {
      rowMeans(
        otu[, groups == g, drop = FALSE],
        na.rm = TRUE
      )
    },
    numeric(nrow(otu))
  )

  rownames(network_matrix) <- rownames(otu)
  colnames(network_matrix) <- group_levels


  ## Remove weak links
  network_matrix[
    network_matrix < min_abundance
  ] <- 0


  ## Remove disconnected taxa and groups
  network_matrix <- network_matrix[
    rowSums(network_matrix) > 0,
    colSums(network_matrix) > 0,
    drop = FALSE
  ]

  if (nrow(network_matrix) == 0L ||
      ncol(network_matrix) == 0L) {
    stop(
      "No interactions remain after network filtering.",
      call. = FALSE
    )
  }


  ## Rows = sample groups, columns = taxa
  network_matrix <- t(network_matrix)

  taxon_metadata <- taxon_metadata[
    match(
      colnames(network_matrix),
      taxon_metadata$network_taxon
    ),
    ,
    drop = FALSE
  ]

  attr(network_matrix, "taxon_metadata") <- taxon_metadata
  attr(network_matrix, "group_variable") <- group
  attr(network_matrix, "taxonomic_rank") <- rank


  if (isTRUE(verbose)) {
    message(
      "Bipartite network: ",
      nrow(network_matrix), " ", group, " categories x ",
      ncol(network_matrix), " taxa | ",
      sum(network_matrix > 0), " links retained."
    )
  }

  network_matrix
}



plot_bipartite_network <- function(
    network,
    type = c("matrix", "fr", "bipartite", "web"),
    labels = TRUE,
    tax_color_rank = NULL,
    output_dir = "plots",
    filename_prefix = "network",
    web_horizontal = TRUE,
    verbose = TRUE) {

  type <- match.arg(type)

  if (!is.matrix(network)) {
    network <- as.matrix(network)
  }

  if (nrow(network) == 0L ||
      ncol(network) == 0L) {
    stop(
      "Network matrix contains no nodes.",
      call. = FALSE
    )
  }

  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  output_file <- file.path(
    output_dir,
    paste0(filename_prefix, "_", type, ".pdf")
  )


  # ========================================================
  # NODE COLOURS

  taxon_metadata <- attr(
    network,
    "taxon_metadata"
  )

  taxon_group <- rep(
    "Taxon",
    ncol(network)
  )

  names(taxon_group) <- colnames(network)


  if (!is.null(tax_color_rank)) {

    if (is.null(taxon_metadata)) {
      stop(
        "Taxonomic metadata are not stored in the network object. ",
        "Recreate the network using make_bipartite_network().",
        call. = FALSE
      )
    }

    if (!tax_color_rank %in% names(taxon_metadata)) {
      stop(
        "Taxonomic colour rank '",
        tax_color_rank,
        "' not found. Available ranks: ",
        paste(
          setdiff(
            names(taxon_metadata),
            "network_taxon"
          ),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    taxon_group <- as.character(
      taxon_metadata[[tax_color_rank]]
    )

    taxon_group[
      is.na(taxon_group) |
        taxon_group == ""
    ] <- "Unassigned"

    names(taxon_group) <-
      taxon_metadata$network_taxon
  }


  ## Colours shared between graph and classical web
  tax_levels <- unique(
    taxon_group[colnames(network)]
  )

  tax_cols <- if (exists(
    "palette_discrete",
    mode = "function"
  )) {
    palette_discrete(length(tax_levels))
  } else {
    grDevices::hcl.colors(
      length(tax_levels),
      palette = "Dark 3"
    )
  }

  names(tax_cols) <- tax_levels


  # ========================================================
  # 1. INTERACTION MATRIX

  if (type == "matrix") {

    matrix_df <- as.data.frame(
      as.table(network),
      stringsAsFactors = FALSE
    )

    names(matrix_df) <- c(
      "group",
      "taxon",
      "weight"
    )

    group_order <- names(
      sort(
        rowSums(network),
        decreasing = TRUE
      )
    )

    taxon_order <- names(
      sort(
        colSums(network),
        decreasing = TRUE
      )
    )

    matrix_df$group <- factor(
      matrix_df$group,
      levels = rev(group_order)
    )

    matrix_df$taxon <- factor(
      matrix_df$taxon,
      levels = taxon_order
    )

    p <- ggplot2::ggplot(
      matrix_df,
      ggplot2::aes(
        x = taxon,
        y = group,
        fill = weight
      )
    ) +
      ggplot2::geom_tile(
        colour = "white",
        linewidth = 0.35
      ) +
      ggplot2::scale_fill_viridis_c(
        trans = scales::pseudo_log_trans(
          base = 10,
          sigma = 0.001
        )
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 55,
          hjust = 1
        ),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = "Taxon",
        y = NULL,
        fill = "Interaction\nstrength"
      )

    ggplot2::ggsave(
      output_file,
      p,
      width = min(
        18,
        max(8, ncol(network) * 0.55)
      ),
      height = min(
        14,
        max(4.5, nrow(network) * 0.65 + 2)
      )
    )
  }


  # ========================================================
  # 2. FRUCHTERMAN-REINGOLD GRAPH

  if (type == "fr") {

    if (!requireNamespace("igraph", quietly = TRUE) ||
        !requireNamespace("ggraph", quietly = TRUE)) {
      stop(
        "Packages 'igraph' and 'ggraph' are required.",
        call. = FALSE
      )
    }

    graph <- igraph::graph_from_biadjacency_matrix(
      network,
      weighted = TRUE
    )

    igraph::V(graph)$level <- ifelse(
      igraph::V(graph)$type,
      "Taxon",
      "Group"
    )

    igraph::V(graph)$node_group <- "Group"

    tax_idx <- igraph::V(graph)$level == "Taxon"

    igraph::V(graph)$node_group[tax_idx] <-
      taxon_group[
        igraph::V(graph)$name[tax_idx]
      ]

    igraph::V(graph)$strength <- igraph::strength(
      graph,
      weights = igraph::E(graph)$weight
    )

    set.seed(1)

    coords <- igraph::layout_with_fr(
      graph,
      weights = igraph::E(graph)$weight,
      niter = 1000
    )

    p <- ggraph::ggraph(
      graph,
      layout = "manual",
      x = coords[, 1],
      y = coords[, 2]
    ) +
      ggraph::geom_edge_link(
        ggplot2::aes(width = weight),
        alpha = 0.25,
        lineend = "round"
      ) +
      ggraph::scale_edge_width(
        range = c(0.2, 3),
        guide = "none"
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(
          colour = node_group,
          shape = level,
          size = strength
        )
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          Group = 15,
          Taxon = 16
        ),
        guide = "none"
      ) +
      ggplot2::scale_size_continuous(
        range = c(3, 7),
        guide = "none"
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          Group = "grey25",
          tax_cols
        )
      ) +
      ggplot2::labs(
        colour = if (is.null(tax_color_rank))
          NULL
        else tax_color_rank
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin = ggplot2::margin(
          15, 20, 15, 20
        )
      )

    if (labels) {
      p <- p +
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          repel = TRUE,
          size = 3,
          show.legend = FALSE
        )
    }

    ggplot2::ggsave(
      output_file,
      p,
      width = min(
        14,
        max(
          8,
          sqrt(igraph::vcount(graph)) * 2.4
        )
      ),
      height = min(
        12,
        max(
          7,
          sqrt(igraph::vcount(graph)) * 2
        )
      )
    )
  }


  # ========================================================
  # 3. IGRAPH BIPARTITE LAYOUT

  if (type == "bipartite") {

    if (!requireNamespace("igraph", quietly = TRUE) ||
        !requireNamespace("ggraph", quietly = TRUE)) {
      stop(
        "Packages 'igraph' and 'ggraph' are required.",
        call. = FALSE
      )
    }

    graph <- igraph::graph_from_biadjacency_matrix(
      network,
      weighted = TRUE
    )

    igraph::V(graph)$level <- ifelse(
      igraph::V(graph)$type,
      "Taxon",
      "Group"
    )

    igraph::V(graph)$node_group <- "Group"

    tax_idx <- igraph::V(graph)$level == "Taxon"

    igraph::V(graph)$node_group[tax_idx] <-
      taxon_group[
        igraph::V(graph)$name[tax_idx]
      ]

    igraph::V(graph)$strength <- igraph::strength(
      graph,
      weights = igraph::E(graph)$weight
    )

    coords <- igraph::layout_as_bipartite(
      graph,
      hgap = 1,
      vgap = 1,
      maxiter = 500
    )

    # Rotate to two vertical columns, which gives labels more room.
    coords <- coords[, c(2, 1)]

    p <- ggraph::ggraph(
      graph,
      layout = "manual",
      x = coords[, 1],
      y = coords[, 2]
    ) +
      ggraph::geom_edge_link(
        ggplot2::aes(width = weight),
        alpha = 0.22,
        lineend = "round"
      ) +
      ggraph::scale_edge_width(
        range = c(0.2, 3),
        guide = "none"
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(
          colour = node_group,
          shape = level,
          size = strength
        )
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          Group = 15,
          Taxon = 16
        ),
        guide = "none"
      ) +
      ggplot2::scale_size_continuous(
        range = c(3, 7),
        guide = "none"
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          Group = "grey25",
          tax_cols
        )
      ) +
      ggplot2::labs(
        colour = if (is.null(tax_color_rank))
          NULL
        else tax_color_rank
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin = ggplot2::margin(
          15, 50, 15, 50
        )
      )

    if (labels) {
      p <- p +
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          repel = TRUE,
          size = 3,
          show.legend = FALSE
        )
    }

    ggplot2::ggsave(
      output_file,
      p,
      width = 10,
      height = min(
        16,
        max(
          7,
          max(nrow(network), ncol(network)) * 0.45
        )
      )
    )
  }


  # ========================================================
  # 4. CLASSICAL BIPARTITE PACKAGE WEB

  if (type == "web") {

    if (!requireNamespace("bipartite", quietly = TRUE)) {
      stop(
        "Package 'bipartite' is required for type = 'web'.",
        call. = FALSE
      )
    }

    ## Taxa are the higher level (columns); groups are the lower level (rows).
    higher_cols <- unname(
      tax_cols[taxon_group[colnames(network)]]
    )
    lower_cols <- rep("grey25", nrow(network))

    ## plotweb requires character vectors for labels.
    higher_labels <- if (labels) {
      colnames(network)
    } else {
      rep("", ncol(network))
    }

    lower_labels <- if (labels) {
      rownames(network)
    } else {
      rep("", nrow(network))
    }

    plot_width <- if (web_horizontal) {
      min(16, max(8, nrow(network) * 0.75 + 6))
    } else {
      min(20, max(10, ncol(network) * 0.75))
    }

    plot_height <- if (web_horizontal) {
      min(18, max(8, ncol(network) * 0.55))
    } else {
      min(12, max(7, nrow(network) * 0.8 + 4))
    }

    grDevices::pdf(
      output_file,
      width = plot_width,
      height = plot_height
    )

    bipartite::plotweb(
      network,
      sorting = "normal",
      empty = FALSE,

      higher_labels = higher_labels,
      lower_labels = lower_labels,

      # Automatically reduces font size to avoid label overlap.
      text_size = "auto",
      spacing = 0.3,

      box_size = 0.08,
      lab_distance = 0.04,

      higher_color = higher_cols,
      lower_color = lower_cols,

      link_color = "higher",
      link_alpha = 0.3,

      curved_links = FALSE,
      horizontal = web_horizontal,

      mar = c(1, 1, 1, 1)
    )

    grDevices::dev.off()
  }


  if (isTRUE(verbose)) {
    message(
      "Network plot (",
      type,
      ") written to: ",
      output_file
    )
  }

  invisible(output_file)
}