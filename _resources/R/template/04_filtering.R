############################################################
# 04 - FILTERING AND AGGREGATION
############################################################


# ==========================================================
# ANALYSIS UNIT AND TAXONOMIC LEVEL

# Choose the biological unit used for downstream analyses:
#
# "asv"         Keep ASVs/pASVs as separate units.
# "taxonomy"    Aggregate sequences at the selected taxonomic level.
# "postcluster" Aggregate resolved species while retaining unresolved
#               postclustered units as approximate species.
#
# Good starting points:
# Bacteria 16S  -> "asv", or genus/family for broader patterns
# Plant ITS2    -> species; "postcluster" if many taxa are unresolved
# Fungal ITS    -> species/genus; ASVs/pASVs if poorly resolved
# Arthropod COI -> species; ASVs/pASVs if poorly resolved
#
# "postcluster" is only available when postclustering was already performed
# during bioinformatics. If uncertain, compare alternative approaches.

analysis_unit <- "postcluster"               # "asv", "taxonomy", or "postcluster"

# If you want to change the level of taxonomic aggregation, specify the rank here. Only relevant if analysis_unit = "taxonomy".
# Valid ranks: kingdom, phylum, class, order, family, genus, species
# analysis_rank <- "genus"    # e.g. "species", "genus", "family"

data.ps.taxa <- select_analysis_unit(
  data.ps,
  method = analysis_unit,
  rank = analysis_rank,
  postcluster_threshold = pipeline_postcluster
)

# ==========================================================
# CONTROL INSPECTION

# Controls can help to judge low-abundance signals before filtering.
#
# Negative controls:
# Low-level reads are common because metabarcoding is highly sensitive.
# Taxa dominating true samples may also dominate negatives at much lower
# absolute abundance. Raw read counts are therefore most informative here.
#
# Positive controls:
# Known positive-control taxa can reveal cross-sample leakage into non-empty
# samples. Their relative abundance in biological samples provides a useful
# empirical indication of the background level.
#
# Adjust these labels only if your metadata use different names.

control_column <- "Type"
negative_label <- "negative"
positive_label <- "positive"

control.results <- inspect_controls(
  data.ps.taxa,
  control_column = control_column,
  negative = negative_label,
  positive = positive_label,
  rank = analysis_rank,
  output_dir = "plots"
)


## Remove controls before downstream analyses
data.ps.taxa <- prune_samples(!sample_names(data.ps.taxa) %in% c(control.results$negative_samples, control.results$positive_samples), data.ps.taxa)
data.ps.taxa <- prune_taxa(taxa_sums(data.ps.taxa) > 0, data.ps.taxa)

# ==========================================================
# OPTIONAL MANUAL TAXON FILTERING

## Remove known irrelevant taxa
# Add project-specific exclusions only if required.
#
# data.ps.taxa <- subset_taxa(
#   data.ps.taxa,
#   !species %in% c(
#     "Microcycas_calocoma",
#     "Ephedra_major"
#   )
# )

# ==========================================================
# ABUNDANCE TRANSFORMATION AND FILTERING

## Relative abundance
# Standardizes sequencing depth across samples and makes abundances comparable between samples.
data.ps.rel <- transform_sample_counts(data.ps.taxa, function(x) x / sum(x))

## Filter low-abundance taxa
# Low-abundance detections may represent genuine rare taxa, but are also more
# susceptible to PCR/sequencing errors, cross-sample leakage and contamination.
# Controls inspected above can help choose a suitable filtering level.
#
# threshold: minimum relative abundance within a sample (0.01 = 1%).
# Detections below this threshold are removed from that sample, while the taxon
# can still be retained if it passes the threshold in other samples.
#
# min_samples: number of samples in which a taxon must reach the threshold.
# Higher values reduce isolated detections, which may represent sequencing
# errors, but can also remove genuinely rare or sample-specific taxa.
#
# Practical starting points:
# 0.001 (0.1%) = permissive; retains more rare taxa, but also more potential noise
# 0.01  (1%)   = moderate; useful general starting point
# 0.05  (5%)   = strict; focuses mainly on dominant taxa
#
# There is no universal cutoff. Consider controls, expected biology and the
# study question and sample replication design. If uncertain, compare settings 
# and check whether conclusions are robust.

abundance_threshold <- 0.01
abundance_min_samples <- 1

filtered <- filter_low_abundance(
  data.ps.taxa,
  threshold = abundance_threshold,
  min_samples = abundance_min_samples
)

data.ps.filter <- filtered$counts
data.ps.rel.filter <- filtered$relative

# ==========================================================
# FILTERED DATA OBJECTS

## Inspect resulting objects
data.ps.taxa
data.ps.filter
data.ps.rel.filter

