# ==========================================================
# INTERACTION NETWORKS

## Build interaction network
# Samples are grouped by one metadata variable and linked to detected taxa.
# Relative abundances are averaged within groups before constructing the network.

network_group <- "host"
network_rank <- analysis_rank
network_min_abundance <- 0.01

network <- make_bipartite_network(
  data.ps.rel.filter,
  group = network_group,
  rank = network_rank,
  min_abundance = network_min_abundance
)


## Plot interaction networks
# matrix: compact overview, particularly useful for larger networks
# fr: force-directed structure based on interactions
# bipartite: explicit two-level igraph network
# web: classical ecological bipartite network

network_labels <- TRUE
network_tax_color <- "family"     # taxonomic rank for colours; NULL for none

for (network_plot in c("matrix", "fr", "bipartite", "web")) {
  plot_bipartite_network(
    network,
    type = network_plot,
    labels = network_labels,
    tax_color_rank = network_tax_color
  )
}

