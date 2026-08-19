# ==========================================================
# SAMPLE DIVERSITY

## Coverage-standardized  Hill diversity
# Uses iNEXT to compare diversity at equal sample completeness rather than equal
# read depth. Hill numbers: q=0 richness, q=1 common taxa, q=2 dominant taxa.

## 1. Inspect sample coverage
# First inspect sample completeness before choosing a common coverage level.
# A suitable target should be reached by most samples without requiring strong
# extrapolation. Coverage is preferable to equal read depth because it compares
# samples at similar completeness.

coverage.results <- inspect_sample_coverage(data.ps.filter)
coverage.results$table

# Suggested starting point based on the observed samples:
coverage.results$suggested_coverage

## 2. Estimate coverage-standardized Hill diversity
# Hill numbers:
# q = 0  richness, sensitive to rare taxa
# q = 1  exponential Shannon diversity, weights taxa by abundance
# q = 2  inverse Simpson diversity, emphasizes dominant taxa
#
# Set the target coverage after inspecting the coverage plot above.

diversity_coverage <- coverage.results$suggested_coverage
# diversity_coverage <- 0.95

hill.results <- estimate_hill_diversity(
  data.ps.filter,
  coverage = diversity_coverage,
  q = c(0, 1, 2)
)

write.csv(hill.results, "tables/sample_hill_diversity.csv", row.names = FALSE)


## 3. Plot Hill diversity by biological group
# Compare coverage-standardized diversity among a metadata grouping variable.

diversity_group <- "host"
diversity_color <- "host"   # set NULL for no second grouping variable

plot_hill_diversity(
  hill.results,
  group = diversity_group,
  color = diversity_color,
  output_file = "plots/sample_hill_diversity.pdf"
)

# ==========================================================
# ORDINATION

## NMDS overview
# Bray-Curtis ordination provides a first view of compositional differences
# among samples without requiring predefined biological groups.

# ==========================================================
# SAMPLE ORDINATION

## Ordination of community composition
# NMDS is a good general choice. For large datasets, PCoA is faster and
# deterministic. Both use Bray-Curtis dissimilarity here.

ordination_method <- "NMDS"         # "NMDS" or "PCoA"
ordination_color <- "host"          # set NULL for no colour grouping
ordination_shape <- "host"          # set NULL for no shape grouping
ordination_labels <- FALSE          # TRUE can help identify unusual samples


## Calculate ordination
ordination <- ordinate(
  data.ps.rel.filter,
  method = ordination_method,
  distance = "bray",
  k = 2,
  trymax = if (ordination_method == "NMDS") 50 else NULL
)


## Plot ordination
p <- plot_ordination(
  data.ps.rel.filter,
  ordination,
  color = ordination_color,
  shape = ordination_shape
) +
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "NMDS - Bray-Curtis")

if (ordination_labels) {
  p <- p +
    geom_text(
      aes(label = rownames(p$data)),
      size = 2.5,
      vjust = -0.7,
      show.legend = FALSE
    )
}


ggsave(
  paste0("plots/sample_ordination_", ordination_method, ".pdf"),
  p,
  width = 7,
  height = 6
)


