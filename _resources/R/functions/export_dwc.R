
# -----------------------------------------------------------------------------
# Darwin Core export helpers
# -----------------------------------------------------------------------------
# These functions are called by 08_export_dwc.R.
# They implement the DwC field derivation and assembly logic so
# the template script stays readable for non-bioinformaticians.
#
# Main entry point: export_dwc()
# Internal helpers are prefixed .dwc_ and not intended for direct use.
# -----------------------------------------------------------------------------


# Derive the resolved taxonRank for each taxon in a taxonomy data.frame.
#
# After replace_tax_prefixes(), unresolved ranks end in " spc" (plain
# propagation) or " spc N" (numbered postcluster unit from
# tax_glom_postcluster_species()). The resolved rank is the rightmost
# column whose value is neither empty nor matches either pattern.
#
# Returns a character vector of rank names, one per taxon row, or NA
# where no rank is resolved.
.dwc_taxon_rank <- function(tax_df) {
  rank_cols <- names(tax_df)
  vapply(seq_len(nrow(tax_df)), function(i) {
    row <- as.character(tax_df[i, ])
    is_unresolved <- is.na(row) | row == "" |
      grepl(" spc( \\d+)?$", row, perl = TRUE)
    resolved <- which(!is_unresolved)
    if (length(resolved) == 0L) return(NA_character_)
    rank_cols[max(resolved)]
  }, character(1))
}


# Derive scientificName for each taxon from the taxonomy data.frame.
#
# Returns the value at the resolved rank. For postclustered placeholder
# units (e.g. "Trifolium spc 1"), strips the " spc N" suffix so the
# exported name is the parent taxon ("Trifolium"), not the internal
# pipeline label. The postcluster context is captured separately in
# identificationRemarks.
#
# tax_df     : taxonomy data.frame (rows = taxa, cols = ranks)
# taxon_rank : character vector of resolved rank names (from .dwc_taxon_rank)
.dwc_scientific_name <- function(tax_df, taxon_rank) {
  vapply(seq_len(nrow(tax_df)), function(i) {
    rk <- taxon_rank[i]
    if (is.na(rk)) return(NA_character_)
    val <- as.character(tax_df[i, rk])
    if (is.na(val) || val == "") return(NA_character_)
    # Strip " spc N" postcluster suffix to recover the parent name.
    val <- sub(" spc( \\d+)?$", "", val, perl = TRUE)
    trimws(val)
  }, character(1))
}


# Extract one metadata column from a sample data.frame by DwC term name.
#
# Uses the col_map list (term -> column name) supplied by the user in the
# template. Returns NA_character_ for every sample when the column is absent,
# and emits a diagnostic message so the user knows what is missing.
#
# md      : data.frame of sample metadata (rows = samples)
# term    : character; the DwC term name used as key in col_map
# col_map : named list mapping DwC terms to metadata column names
.dwc_meta_col <- function(md, term, col_map) {
  col <- col_map[[term]]
  if (is.null(col) || (length(col) == 1L && is.na(col))) {
    return(rep(NA_character_, nrow(md)))
  }
  if (!col %in% names(md)) {
    message("  DwC note: metadata column '", col, "' (mapped to '", term,
            "') not found in sample_data — left empty.")
    return(rep(NA_character_, nrow(md)))
  }
  as.character(md[[col]])
}


# Parse combined lat/lon strings into signed decimal degrees.
#
# Handles the three formats that commonly appear in MIxS templates and
# lab notebooks:
#   "48.1351 N 11.5820 E"  (INSDC/MIxS format, with compass letters)
#   "48.1351N 11.5820E"    (compact variant)
#   "48.1351, 11.5820"     (already signed decimal, no compass)
#
# Returns a list(lat, lon), each a numeric vector of length(x).
# Entries that cannot be parsed are NA.
.dwc_parse_latlon <- function(x) {
  lat <- rep(NA_real_, length(x))
  lon <- rep(NA_real_, length(x))
  compass_re <- "^([0-9.]+)\\s*([NS])\\s+([0-9.]+)\\s*([EW])$"

  for (i in seq_along(x)) {
    s <- trimws(x[i])
    if (is.na(s) || s == "") next

    # Try "DD.DDD N EE.EEE E" style first.
    m <- regmatches(s, regexec(compass_re, s,
                               perl = TRUE, ignore.case = TRUE))[[1]]
    if (length(m) == 5L) {
      lat_val <- as.numeric(m[2])
      lon_val <- as.numeric(m[4])
      lat[i]  <- if (toupper(m[3]) == "S") -lat_val else lat_val
      lon[i]  <- if (toupper(m[5]) == "W") -lon_val else lon_val
      next
    }

    # Try "DD.DDD, EE.EEE" (signed decimal, comma or space separated).
    parts <- suppressWarnings(as.numeric(strsplit(s, "[,\\s]+")[[1]]))
    if (length(parts) == 2L && !anyNA(parts)) {
      lat[i] <- parts[1]
      lon[i] <- parts[2]
    }
  }

  list(lat = lat, lon = lon)
}


# Build and write Darwin Core tables from a processed phyloseq object.
#
# This is the main function called from 08_export_dwc.R. It produces three
# files in the tables/ directory:
#   dwc_occurrence.csv     — Occurrence core (one row per taxon x sample)
#   dwc_dna_extension.csv  — DNA-derived data extension (one row per taxon)
#   dwc_fields_guide.txt   — plain-text guide to empty columns
#
# Arguments
# ---------
# physeq            : phyloseq object with count data and sample metadata.
#                     Use data.ps.filter (section 04): controls removed,
#                     taxonomy cleaned, metadata merged.
# institution_code  : character; your institution code for occurrenceID,
#                     e.g. "LMU". Register at https://www.gbif.org/grscicoll
# project           : character; project name (set in section 01 header).
# marker            : character; amplicon marker, e.g. "ITS2", "16S", "COI".
# col_map           : named list mapping DwC terms to sample_data() column
#                     names. Set to NA for terms not present in your metadata.
# target_gene       : character; DwC DNA extension target_gene value.
# target_subfragment: character; DwC DNA extension target_subfragment value.
# pipeline_postcluster : integer; postcluster identity threshold from pipeline
#                     (0 = disabled). Used in identificationRemarks.
# tax_threshold     : numeric; classification identity threshold from config.
# sintax_cutoff     : numeric; SINTAX bootstrap cutoff from config.
# use_blast_sintax_combination : integer (0/1); from config.
# output_dir        : character; directory for output files (default "tables").
# verbose           : logical; print progress messages.
#
# Returns invisibly: list(occurrence = dwc_occ, dna = dwc_dna)
export_dwc <- function(
    physeq,
    institution_code,
    project,
    marker,
    col_map,
    target_gene,
    target_subfragment       = NA_character_,
    pipeline_postcluster     = 0,
    tax_threshold            = NA,
    sintax_cutoff            = NA,
    use_blast_sintax_combination = 0,
    output_dir               = "tables",
    verbose                  = TRUE) {

  .assert_phyloseq(physeq)

  if (!is.character(institution_code) || !nzchar(institution_code)) {
    stop("institution_code must be a non-empty character string.", call. = FALSE)
  }
  if (!is.list(col_map) || is.null(names(col_map))) {
    stop("col_map must be a named list.", call. = FALSE)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # ------------------------------------------------------------------
  # Long occurrence table: one row per taxon x sample, non-zero only.
  # psmelt() includes zero-abundance rows; we drop those — an absence
  # is not a DwC occurrence.
  if (verbose) message("DwC export: building Occurrence core ...")

  occ_long <- phyloseq::psmelt(physeq)
  occ_long <- occ_long[occ_long$Abundance > 0L, , drop = FALSE]

  if (nrow(occ_long) == 0L) {
    stop("DwC export: no non-zero occurrences found in the phyloseq object.",
         call. = FALSE)
  }

  # ------------------------------------------------------------------
  # Taxonomy: derive taxonRank and scientificName per taxon.
  # tax_table columns are fixed ranks; values after replace_tax_prefixes()
  # are plain names or " spc"-suffixed placeholders.
  tax_df <- as.data.frame(phyloseq::tax_table(physeq),
                           stringsAsFactors = FALSE)

  taxon_rank_vec            <- .dwc_taxon_rank(tax_df)
  names(taxon_rank_vec)     <- rownames(tax_df)

  sci_name_vec              <- .dwc_scientific_name(tax_df, taxon_rank_vec)
  names(sci_name_vec)       <- rownames(tax_df)

  # ------------------------------------------------------------------
  # Sample metadata.
  md          <- as.data.frame(phyloseq::sample_data(physeq),
                                stringsAsFactors = FALSE)
  md$.Sample  <- rownames(md)

  # Lat/Lon: prefer explicit separate columns; fall back to parsing a
  # combined string (e.g. "*lat_lon" from the pipeline QC metadata).
  lat_col <- col_map[["decimalLatitude"]]
  lon_col <- col_map[["decimalLongitude"]]

  if (!is.null(lat_col) && !is.na(lat_col) && lat_col %in% names(md) &&
      !is.null(lon_col) && !is.na(lon_col) && lon_col %in% names(md)) {
    lat_vec <- as.numeric(md[[lat_col]])
    lon_vec <- as.numeric(md[[lon_col]])
  } else {
    combined_raw <- .dwc_meta_col(md, "latlon_combined", col_map)
    parsed       <- .dwc_parse_latlon(combined_raw)
    lat_vec      <- parsed$lat
    lon_vec      <- parsed$lon
    if (verbose && any(!is.na(lat_vec))) {
      message("  DwC: parsed lat/lon from combined '",
              col_map[["latlon_combined"]], "' column.")
    }
  }
  names(lat_vec) <- rownames(md)
  names(lon_vec) <- rownames(md)

  # Match sample index once, used for every metadata column below.
  sample_idx <- match(occ_long$Sample, md$.Sample)

  # ------------------------------------------------------------------
  # Pipeline citation — used in both Occurrence core (references field)
  # and DNA extension (otu_class_appr field).
  # Leonhardt, Peters & Keller (2022) Phil. Trans. R. Soc. B 377: 20210171.
  pipeline_citation <- paste0(
    "Leonhardt SD, Peters B, Keller A (2022) ",
    "Do amino and fatty acid profiles of pollen provisions correlate with ",
    "bacterial microbiomes in the mason bee Osmia bicornis? ",
    "Philosophical Transactions of the Royal Society B 377 (1853): 20210171. ",
    "https://doi.org/10.1098/rstb.2021.0171"
  )

  # ------------------------------------------------------------------
  # Flag postclustered units: species column ends in " spc N".
  sp_col         <- if ("species" %in% names(occ_long)) occ_long$species else ""
  is_postcluster <- grepl(" spc \\d+$", as.character(sp_col), perl = TRUE)

  # Identification remarks string: records pipeline settings per row.
  id_remarks <- paste0(
    "marker=", marker,
    if (!is.na(tax_threshold))  paste0("; tax_threshold=",  tax_threshold)  else "",
    if (!is.na(sintax_cutoff))  paste0("; sintax_cutoff=",  sintax_cutoff)  else "",
    "; pipeline_postcluster=", pipeline_postcluster,
    ifelse(taxon_rank_vec[occ_long$OTU] != "species" &
             !is.na(taxon_rank_vec[occ_long$OTU]),
           paste0("; resolved_to=", taxon_rank_vec[occ_long$OTU]), ""),
    ifelse(is_postcluster, "; postclustered_unit=TRUE", ""),
    if (isTRUE(use_blast_sintax_combination == 1))
      "; classification=BLAST_LCA+SINTAX" else "; classification=SINTAX"
  )

  # ------------------------------------------------------------------
  # Assemble Occurrence core data.frame.
  dwc_occ <- data.frame(

    occurrenceID         = paste(institution_code, project,
                                 occ_long$Sample, occ_long$OTU, sep = ":"),
    institutionCode      = institution_code,
    datasetName          = project,
    basisOfRecord        = "MaterialSample",

    eventID              = occ_long$Sample,
    materialSampleID     = .dwc_meta_col(md, "materialSampleID",
                                          col_map)[sample_idx],
    eventDate            = .dwc_meta_col(md, "eventDate",
                                          col_map)[sample_idx],
    samplingProtocol     = .dwc_meta_col(md, "samplingProtocol",
                                          col_map)[sample_idx],

    country              = .dwc_meta_col(md, "country",
                                          col_map)[sample_idx],
    locality             = .dwc_meta_col(md, "locality",
                                          col_map)[sample_idx],
    decimalLatitude      = lat_vec[occ_long$Sample],
    decimalLongitude     = lon_vec[occ_long$Sample],
    geodeticDatum        = ifelse(!is.na(lat_vec[occ_long$Sample]),
                                  "WGS84", NA_character_),

    scientificName       = sci_name_vec[occ_long$OTU],
    kingdom              = occ_long$kingdom,
    phylum               = occ_long$phylum,
    class                = occ_long$class,
    order                = occ_long$order,
    family               = occ_long$family,
    genus                = occ_long$genus,
    taxonRank            = taxon_rank_vec[occ_long$OTU],

    occurrenceStatus     = "present",
    organismQuantity     = as.integer(occ_long$Abundance),
    organismQuantityType = "DNA sequence reads",
    associatedTaxa       = .dwc_meta_col(md, "associatedTaxa",
                                          col_map)[sample_idx],

    identificationRemarks = id_remarks,

    # DwC references field: cites the pipeline used to generate the data.
    references = pipeline_citation,

    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  # Remove " spc"-placeholder values from taxonomy rank columns —
  # these are pipeline artefacts for phyloseq, not valid taxon names.
  spc_re <- " spc( \\d+)?$"
  for (col in c("kingdom", "phylum", "class", "order", "family", "genus")) {
    if (col %in% names(dwc_occ)) {
      dwc_occ[[col]] <- ifelse(
        grepl(spc_re, dwc_occ[[col]], perl = TRUE),
        NA_character_,
        dwc_occ[[col]]
      )
    }
  }

  # ------------------------------------------------------------------
  # DNA-derived data extension: one row per unique taxon.
  if (verbose) message("DwC export: building DNA-derived data extension ...")

  dna_taxa <- unique(occ_long$OTU)

  # Sequences are only available when analysis_unit = "asv" and refseq()
  # is stored in the phyloseq object. For postclustered or taxonomy-
  # aggregated objects the field is left empty; sequences are in asvs.merge.fa.
  seq_vec        <- rep(NA_character_, length(dna_taxa))
  names(seq_vec) <- dna_taxa

  rs_slot <- tryCatch(
    phyloseq::refseq(physeq, errorIfNULL = FALSE),
    error = function(e) NULL
  )
  if (inherits(rs_slot, "XStringSet")) {
    rs   <- as.character(rs_slot)
    hits <- names(rs) %in% dna_taxa
    seq_vec[names(rs)[hits]] <- rs[hits]
  }

  # Bioinformatics method strings for the DNA extension.
  # otu_class_appr is the standard DwC field for bioinformatics provenance.
  otu_class_appr <- paste0(
    "Keller metabarcoding pipeline (https://github.com/chiras/metabarcoding_pipeline); ",
    "VSEARCH UNOISE3 denoising",
    if (pipeline_postcluster > 0)
      paste0("; postclustered at ", pipeline_postcluster, "% identity")
    else "",
    "; chimera removal uchime3_denovo; ",
    "cite: ", pipeline_citation
  )

  tax_class_appr <- paste0(
    "VSEARCH SINTAX; cutoff=", sintax_cutoff,
    if (isTRUE(use_blast_sintax_combination == 1))
      "; combined BLAST LCA + SINTAX" else ""
  )

  dwc_dna <- data.frame(
    DNA_sequence_ID         = dna_taxa,
    DNA_sequence            = seq_vec,
    target_gene             = target_gene,
    target_subfragment      = target_subfragment,
    # Primer fields: filled by user from config.txt if needed.
    pcr_primer_name_forward = NA_character_,
    pcr_primer_name_reverse = NA_character_,
    pcr_primer_forward      = NA_character_,
    pcr_primer_reverse      = NA_character_,
    otu_class_appr          = otu_class_appr,
    otu_seq_comp_appr       = paste0("VSEARCH usearch_global; id=",
                                      if (!is.na(tax_threshold))
                                        tax_threshold / 100 else "NA"),
    otu_db                  = NA_character_,
    tax_class_appr          = tax_class_appr,
    # taxonID deferred: run verify_taxonomy.py after deposition.
    tax_class_id            = NA_character_,
    stringsAsFactors        = FALSE,
    check.names             = FALSE
  )

  # ------------------------------------------------------------------
  # Write output files.
  occ_file   <- file.path(output_dir, "dwc_occurrence.csv")
  dna_file   <- file.path(output_dir, "dwc_dna_extension.csv")
  guide_file <- file.path(output_dir, "dwc_fields_guide.txt")

  write.csv(dwc_occ, occ_file, row.names = FALSE, na = "")
  write.csv(dwc_dna, dna_file, row.names = FALSE, na = "")

  # ------------------------------------------------------------------
  # Fields guide: tell the user what is empty and what fills it.
  empty_occ <- names(dwc_occ)[
    vapply(dwc_occ, function(x) all(is.na(x)), logical(1))
  ]
  empty_dna <- names(dwc_dna)[
    vapply(dwc_dna, function(x) all(is.na(x)), logical(1))
  ]

  combined_col <- col_map[["latlon_combined"]]

  guide_lines <- c(
    paste("Darwin Core fields guide —", project),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "=== REQUIRED for GBIF submission ===",
    "  eventDate         : ISO 8601 date (YYYY-MM-DD or YYYY-MM).",
    "                      Supply in biological metadata.",
    "  decimalLatitude   : Signed decimal degrees (e.g. 48.1351).",
    "  decimalLongitude  : Signed decimal degrees (e.g. 11.5820).",
    if (!is.null(combined_col) && !is.na(combined_col))
      paste0("                      OR supply '", combined_col,
             "' as a combined string and re-run.")
    else "",
    "",
    "=== RECOMMENDED for GBIF submission ===",
    "  samplingProtocol  : e.g. 'Malaise trap', 'Floral swab'.",
    "  country           : ISO 3166-1 alpha-2 or country name.",
    "  locality          : More specific location description.",
    "  associatedTaxa    : Host, e.g. 'host: Apis mellifera'.",
    "",
    "=== DEFERRED (fill after deposition) ===",
    "  tax_class_id      : taxonID — GBIF backbone or UNITE SH accession.",
    "                      Run verify_taxonomy.py after sequences are deposited.",
    "  pcr_primer_*      : Primer names and sequences from your config.txt.",
    "  otu_db            : Reference database name and version.",
    "",
    "=== EMPTY OCCURRENCE COLUMNS IN THIS EXPORT ===",
    if (length(empty_occ) > 0) paste(" ", empty_occ) else "  (none)",
    "",
    "=== EMPTY DNA EXTENSION COLUMNS IN THIS EXPORT ===",
    if (length(empty_dna) > 0) paste(" ", empty_dna) else "  (none)"
  )

  writeLines(guide_lines, guide_file)

  if (verbose) {
    message("DwC export complete.")
    message("  Occurrence rows   : ", nrow(dwc_occ))
    message("  Unique taxa       : ", nrow(dwc_dna))
    message("  Empty occ cols    : ", length(empty_occ))
    message("  Empty DNA cols    : ", length(empty_dna))
    message("  Files written to  : ", output_dir, "/")
    if (length(empty_occ) > 0 || length(empty_dna) > 0) {
      message("  See dwc_fields_guide.txt for what biological metadata fills them.")
    }
  }

  invisible(list(occurrence = dwc_occ, dna = dwc_dna))
}
