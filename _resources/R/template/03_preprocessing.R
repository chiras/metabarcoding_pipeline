############################################################
# 03 - PREPROCESSING AND QC
############################################################


# ==========================================================
# TAXONOMY PREPROCESSING

## Propagate incomplete taxonomy
# Propagate lineage for such taxa classified only to higher taxonomic levels.
data.ps <- propagate_incomplete_taxonomy(data.ps)

## Remove unresolved and non-target taxa
# Filtering rules depend on the marker and are reported to the console.
data.ps <- remove_unresolved_taxa(data.ps, marker = marker)

## Clean taxonomic labels
data.ps <- replace_tax_prefixes(data.ps)

## Add taxonomic assignment to ASV IDs
# Example: ASV123 -> Trifolium_pratense.ASV123
data.ps <- add_taxon_to_asv_id(data.ps)


# ==========================================================
# SAMPLE QC

## Label low-throughput samples
# Samples below the threshold are retained but labelled "|LT<threshold>".
# Adjust the threshold if required.

lt.result <- label_low_throughput(
  data.ps,
  threshold = 1000
)

data.ps <- lt.result$data

# Full list of flagged samples:
# lt.result$low_throughput


# ----------------------------------------------------------
# Optional sample labelling

## Add metadata information to sample names
# Can be useful for plots or manual inspection.
# data.ps <- label_sample_by_host(data.ps, "host", "project")


# ==========================================================
# PREPROCESSED DATA

## Inspect resulting object
data.ps


