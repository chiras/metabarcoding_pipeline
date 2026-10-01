############################################################
# 08 - DARWIN CORE EXPORT 
# good for interoperability and metadata public deposition

# THIS HAS NOT BEEN TESTED YET. Please check the output files carefully before using. 
# ERRORS in the output files are likely if your metadata columns are not mapped correctly to DwC terms.
# SCRIPT MIGHT BREAK.

#
# Exports results as Darwin Core tables for public deposition
# to GBIF, OBIS, or similar biodiversity repositories.
#
# DwC standard       : https://dwc.tdwg.org/terms/
# DNA-derived data   : https://rs.gbif.org/extension/gbif/1.0/
#                      dna_derived_data_2021-07-05.xml
#
# Output files (tables/):
#   dwc_occurrence.csv     — one row per taxon x sample detected
#   dwc_dna_extension.csv  — one row per unique taxon (sequences,
#                            primers, bioinformatics provenance)
#   dwc_fields_guide.txt   — plain-text list of empty columns and
#                            what biological metadata fills them
#
# This script uses data.ps.filter from section 04 (controls
# already removed, taxonomy cleaned, metadata merged).
# Add biological metadata (dates, coordinates, host) in
# section 02 via merge_sample_metadata() before running this.
############################################################


# ==========================================================
# 1. INSTITUTION CODE
# ==========================================================

# Your institution's code, used to build globally unique
# occurrence IDs. Register at https://www.gbif.org/grscicoll
# Examples: "LMU", "ZFMK", "NATURALIS", "SMNS"
dwc_institution_code <- "MSB"


# ==========================================================
# 2. METADATA COLUMN MAPPING
# ==========================================================

# Map your sample_data() column names to DwC terms.
# Set a term to NA if that information is not in your metadata.
#
# The columns on the right come from your biological metadata
# file loaded in section 02. Column names with * are from the
# pipeline QC metadata (samples_metadata.csv).

dwc_col_map <- list(
  eventDate         = "*collection_date",  # ISO 8601: YYYY-MM-DD
  decimalLatitude   = NA,      # if you have separate lat column, put name here
  decimalLongitude  = NA,      # if you have separate lon column, put name here
  latlon_combined   = "*lat_lon",  # pipeline combined string, parsed automatically
  country           = "*geo_loc_name",
  locality          = NA,      # e.g. "locality" or "site_name" in your metadata
  materialSampleID  = "source_material_id",
  associatedTaxa    = "*host", # host plant/animal, e.g. "Apis mellifera"
  samplingProtocol  = NA       # e.g. "sampling_method" in your metadata
)


# ==========================================================
# 3. MARKER / TARGET GENE
# ==========================================================

# These are derived automatically from the marker variable
# set in section 01. Only change if your marker is not listed.
dwc_target_gene <- switch(
  toupper(marker),
  ITS2  = "ITS2",
  FITS  = "ITS1-ITS2",
  `16S` = "16S rRNA",
  COI   = "COI",
  marker   # fallback: use the marker string as-is
)

dwc_target_subfragment <- switch(
  toupper(marker),
  ITS2  = "ITS2",
  FITS  = "ITS",
  `16S` = "V4",
  COI   = "5P",
  NA_character_
)


# ==========================================================
# 4. EXPORT
# ==========================================================

export_dwc(
  physeq                       = data.ps.filter,
  institution_code             = dwc_institution_code,
  project                      = project,
  marker                       = marker,
  col_map                      = dwc_col_map,
  target_gene                  = dwc_target_gene,
  target_subfragment           = dwc_target_subfragment,
  pipeline_postcluster         = pipeline_postcluster,
  tax_threshold                = tax_threshold,
  sintax_cutoff                = sintax_cutoff,
  use_blast_sintax_combination = if (exists("use_blast_sintax_combination"))
                                   use_blast_sintax_combination else 0,
  output_dir                   = "tables"
)
