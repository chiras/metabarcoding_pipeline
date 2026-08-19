############################################################
# 05 - BASIC PLOTS
############################################################


# ==========================================================
# TAXON ABUNDANCE

## Most abundant taxa by metadata group
# Shows the cumulative relative abundance of the most abundant taxa, partitioned
# by one biological metadata variable. Change group_column as required.

top_n_taxa <- 20
group_column <- "host"

plot_df <- psmelt(data.ps.rel.filter)

if (!group_column %in% names(plot_df)) {
  stop(
    "Metadata column '", group_column,
    "' not found. Available columns: ",
    paste(names(sample_data(data.ps.rel.filter)), collapse = ", ")
  )
}

# Sum relative abundance across samples for each taxon and group.
plot_df <- plot_df %>%
  group_by(.data[[analysis_rank]], .data[[group_column]]) %>%
  summarise(
    cumulative_abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Select the most abundant taxa across all groups.
top_taxa <- plot_df %>%
  group_by(.data[[analysis_rank]]) %>%
  summarise(
    total_abundance = sum(cumulative_abundance),
    .groups = "drop"
  ) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = top_n_taxa)

plot_df <- plot_df %>%
  filter(.data[[analysis_rank]] %in% top_taxa[[analysis_rank]])

# Order taxa by their total abundance.
plot_df[[analysis_rank]] <- factor(
  plot_df[[analysis_rank]],
  levels = rev(top_taxa[[analysis_rank]])
)

pdf("plots/cumulative_relative_abundance.pdf", width = 8, height = 7)

print(
  ggplot(plot_df, aes(x = .data[[analysis_rank]], y = cumulative_abundance, fill = .data[[group_column]])) +
    geom_col() +
    scale_fill_manual(
        values = palette_discrete(
            dplyr::n_distinct(plot_df[[group_column]], na.rm = TRUE)
        )
    ) +
    coord_flip() +
    theme_bw() +
    labs(x = NULL,
      y = "Cumulative relative abundance",
      fill = group_column,
      title = "Most abundant taxa"
    )
)
dev.off()


# ==========================================================
# SAMPLE COMPOSITION

## Select taxonomic level for visualization
# Automatically chooses the taxonomic rank giving approximately the requested
# number of categories. Alternatively, set plot_rank manually (e.g. "genus").

plot_target_taxa <- 20
plot_rank <- select_taxonomic_rank(data.ps.rel.filter, target_taxa = plot_target_taxa)
# plot_rank <- "genus"


## Relative abundance across samples
data.ps.plot <- tax_glom(data.ps.rel.filter, taxrank = plot_rank)
taxa_names(data.ps.plot) <- as.character(tax_table(data.ps.plot)[, plot_rank])

data.melt <- psmelt(data.ps.plot)
data.melt$Sample <- factor(data.melt$Sample, levels = sample_names(data.ps.plot))

pdf("plots/sample_relative_abundance.pdf",
    width = min(14, max(7, nsamples(data.ps.plot) / 8)),
    height = 6)

print(
  ggplot(
    data.melt,
    aes(x = Sample, y = Abundance, fill = OTU)
  ) +
    geom_col(width = 1) +
    scale_fill_manual(
    values = palette_discrete(
        dplyr::n_distinct(data.melt$OTU, na.rm = TRUE)
    )
    )+
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      axis.text.x = if (nsamples(data.ps.plot) <= 60)
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
      else element_blank(),
      axis.ticks.x = if (nsamples(data.ps.plot) <= 60)
        element_line()
      else element_blank(),
      legend.key.height = grid::unit(0.35, "cm")
    ) +
    labs(
      x = NULL,
      y = "Relative abundance",
      fill = plot_rank,
      title = paste("Sample composition at", plot_rank, "level")
    )
)
dev.off()

## Distribution of taxa across samples
# Choose taxonomic and metadata variables used to structure the plot.
# Multiple variables can be supplied to either dimension.

## Distribution of taxa across samples
# Taxonomic ranks structure the rows; metadata variables structure the columns.

distribution_tax_ranks <- c("order", "family")
distribution_facet_x <- c("host")

## Optional pseudo-log scaling
# Useful for displaying both dominant and low-abundance taxa. Unlike a true
# log transformation, pseudo-log scaling can display zero abundance.
distribution_log_scale <- TRUE

data.melt <- psmelt(data.ps.rel.filter)

p <- ggplot(
  data.melt,
  aes(
    x = OTU,
    y = Abundance,
    fill = .data[[distribution_tax_ranks[length(distribution_tax_ranks)]]]
  )
) +
  geom_col() +
  facet_grid(
    rows = vars(!!!rlang::syms(distribution_tax_ranks)),
    cols = vars(!!!rlang::syms(distribution_facet_x)),
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    x = "Taxa",
    y = "Relative abundance",
    fill = distribution_tax_ranks[length(distribution_tax_ranks)]
  )


if (distribution_log_scale) {
  p <- p +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10, sigma = 1),
      breaks = c(0.01, 0.1, 0.25, 0.5,0.75, 1),
      labels = c("1%", "10%", "25%", "50%", "75%", "100%"),
      limits = c(0, 1),
      expand = c(0, 0)

    )
} else {
  p <- p +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::percent,
      expand = c(0, 0)
    )
}

ggsave(
  "plots/taxon_distribution_across_samples.pdf",
  p,
  width = 14,
  height = max(6, ntaxa(data.ps.rel.filter) / 5)
)

