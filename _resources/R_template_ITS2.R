# Loading in data
## Taxonomy
data.tax <- tax_table(as.matrix(read.table("taxonomy.vsearch", header=T,row.names=1,fill=T,sep=",")))

## Community table
data.otu <- otu_table(read.table("asv_table.merge.txt"), taxa_are_rows=T)

## if phylogeny included
# data.tre <- read.tree("asvs.tre")

## Sample metadata (second line optional if sample names include "-"):
data.map <- 	sample_data(read.table("samples_metadata.csv", header=T, row.names=1,  sep=";", fill=T))
sample_names(data.map) <- gsub("-",".",sample_names(data.map))

## check metadata vs. samples in sequencing data consistency
sample_names(data.map)[!(sample_names(data.map) %in% sample_names(data.otu))]
sample_names(data.otu)[!(sample_names(data.otu) %in% sample_names(data.map))]

## merge the three tables to a single phyloseq object
(data.comp <- merge_phyloseq(data.otu,data.tax,data.map))
## if phylogeny included: 
# (data.comp <- merge_phyloseq(data.otu,data.tax,data.map,data.tre))


# create fake metadata if needed
data.comp <- fill_pseudo_metadata(data.comp)

# preprocessing data pt.1 :
## given hierarchical classification options at the end, we have to propagate the taxonomy over taxonomic levels to not throw out stuff only classified to higher tax levels
data.comp <- propagate_incomplete_taxonomy(data.comp)

## filtering irrelevant taxa, zb. unresolved, algae, fungi etc
data.comp.filter <- remove_unresolved_taxa(data.comp)

#(data.comp.filter = subset_taxa(data.comp.filter, kingdom!=""))

## Make taxa labels nice for plots
data.comp.filter <- replace_tax_prefixes(data.comp.filter)
taxa_names(data.comp.filter) <- interaction(taxa_names(data.comp.filter),tax_table(data.comp.filter)[,length(colnames(tax_table(data.comp.filter)))])

### Check the names
tail(tax_table(data.comp.filter))

## Multiple ASVs might represent the same species, here they are collated
(data.species <- tax_glom(data.comp.filter,taxrank="species"))
taxa_names(data.species) <- tax_table(data.species)[,"species"]

# alternatively try if you used postclustering: 
# data.species2 <- tax_glom_species_filtered(data.comp.filter,rank="species")


## (optional) Rename samples by type different metadata factors
#(data.species <- label_sample_by_host(data.species,"host","project"))

## (optional) Label samples with low throughput with LT
(data.species <- label_low_throughput(data.species , 500))
sample_names(data.species)

# Transform to relative data
data.species.rel = transform_sample_counts(data.species, function(x) x/sum(x))

## (optional) consider looking at controls and
control_samples = c("negative","positive")

pdf("plots/controls_absolute.pdf", width=5, height=5)

controls <- subset_samples(data.species, Type %in% control_samples)
controls <- prune_taxa(taxa_sums(controls)>100, controls)
controls.melt <- psmelt(controls)
ggplot(controls.melt, aes(x=species, y=Abundance, col=Type))+geom_boxplot()+  
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
dev.off()

pdf("plots/controls_relative.pdf", width=15, height=5)

controls.rel <- subset_samples(data.species.rel, Type %in% control_samples)
controls.rel <- prune_taxa(taxa_sums(controls.rel)>0.1, controls.rel)
controls.melt.rel <- psmelt(controls.rel)
controls.melt.rel  <- controls.melt.rel %>% filter(Abundance > 0.001)

ggplot(controls.melt.rel, aes(x=species, y=Abundance, col=Type))+geom_boxplot()+scale_y_log10()+ 
  geom_hline(yintercept = c(0.01, 0.1, 0.001, 0.0001), 
             linetype = "dashed", color = "red") +
  annotate("text", x = 10, y = 0.01, 
           label = "0.01", hjust = 1, vjust = -0.5, color = "red") +
  annotate("text", x = 10, y = 0.1, 
           label = "0.1", hjust = 1, vjust = -0.5, color = "red") +
  annotate("text", x = 10, y = 0.001, 
           label = "0.001", hjust = 1, vjust = -0.5, color = "red") +
  annotate("text", x = 10, y = 0.0001, 
           label = "0.0001", hjust = 1, vjust = -0.5, color = "red")+  
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
dev.off()

data.species <- subset_samples(data.species.rel, !(Type %in% control_samples))
data.species <- subset_samples(data.species.rel, !(Type %in% control_samples))

# low abundance filtering
otu_table(data.species.rel)[otu_table(data.species.rel)<0.01 ]<-0
otu_table(data.species)[otu_table(data.species.rel)<0.01 ]<-0
data.species.filter		= prune_taxa(taxa_sums(data.species)>0, data.species)
(data.species.rel.filter = prune_taxa(rowSums(otu_table(data.species.rel))>0, data.species.rel))

# define paramters for plot definitions
ntaxa <- length(taxa_names(data.species.rel.filter))



# # check GBIF location records
# species_records_long <- get_species_occurences_long(data.species.rel.filter, notInGBIF=T)

# species_records_long$regCountry <- (interaction(species_records_long$region,species_records_long$sub.region))
# species_records_long <- species_records_long[order(species_records_long$region),]
# species_records_long$regCountry <- factor(species_records_long$regCountry, levels = unique(sort(as.character(species_records_long$regCountry))))

# species_records_long <- cbind(species_records_long,tax_table(data.species.rel.filter)[species_records_long$Species,])

# species_records_long$SpecAbund <- interaction(species_records_long$Species,sprintf("%.3f", round(species_records_long$Abundance, digits=4)), sep=" | ")

# pdf("plots/species_occurence.pdf", width=15, height=ntaxa/5)

# ggplot(species_records_long, aes(fill=regCountry, y=value, x=SpecAbund, alpha=log(Abundance*1000,10))) + 
#   facet_grid(order+family~region, space = "free", scales = "free",switch = "y")+
#   geom_bar(position="stack", stat="identity")+
#   theme(strip.text.y.left = element_text(angle = 0))+
#   theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) + 
#   xlab("Species | cumulative relative Abundance")+
#   ylab("GBIF record abundance (log)")+
#   theme(legend.position="bottom")+
#   coord_flip()+ 
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
#                      breaks = scales::trans_breaks("log10", function(x) 10^x), 
#                      labels = scales::trans_format("log10", scales::math_format(10^.x))
#   )#+ annotation_logticks(sides = "b", )  
# dev.off()


# Filtering taxa that are not relevant 
data.species.rel.filter <- subset_taxa(data.species.rel.filter, !species %in% c("Microcycas calocoma","Ephedra major","Impatiens acuminata"))


# First diversity and community metrics and graphs
# Distribution of major taxa, accumulated over all samples
par(mar=c(4,15,1,1), mfrow=c(1,1))
barplot(t(as.data.frame(sort(taxa_sums(data.species.rel.filter), decreasing=T)[1:20])), las=2, horiz=T)

## with ggplot
plot_df <- data.frame(cbind((taxa_sums(data.species.rel.filter))))
colnames(plot_df) <- c("cummulative")

plot_df$species <- rownames(plot_df)
plot_df <- plot_df[order(plot_df$cummulative, decreasing = T),]
plot_df_long <-  pivot_longer(plot_df[1:20,], cols=c(1:1), names_to = "Source", values_to = "relabundance")

pdf("plots/cum_rel_abundance.pdf", width=5, height=5)
ggplot(plot_df_long, aes(fill=Source, y=relabundance, x=reorder(species, -relabundance))) +
  geom_bar(position="stack", stat="identity")+coord_flip()+theme_bw()+ scale_fill_grey(start = 0.3, end = .9) + theme(legend.position = c(0.8, 0.2)) + theme(legend.background = element_rect(colour="black", linewidth=0.5, linetype="solid"))+
  labs(x="Taxa",y="Relative abundance",fill="Source")
dev.off()

## Distribution over samples
data.melt <- psmelt(data.species.rel.filter)
data.melt$Sample <- as.factor(data.melt$Sample)

pdf("plots/sample_rel_abundance.pdf", width=20, height=ntaxa/5)
ggplot(data.melt, aes(OTU, Abundance, fill= family)) +
  facet_grid(order+family~location, space = "free", scales = "free",switch = "y")+
  theme_bw()+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_bar(stat="identity")+
  theme(legend.position="bottom")+
  theme(strip.text.y.left = element_text(angle = 0))+
  theme(strip.text.x = element_text(angle = 90))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),axis.title.x = element_text(family = "sans", size = 15)) + xlab("Species")+
  coord_flip()+ scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))

dev.off()


## Richness (Observed) / Shannon Diversity (replace x/col by group in metadata)
pdf("plots/sample_diversity.pdf", width=25, height=10)
p <- plot_richness(data.species.filter,x= "treatment" , col= "location", measures=c("Observed","Shannon"))+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1) ) + 
  geom_boxplot()+
  scale_color_npg()+
  geom_point(size=4, alpha=0.4,position = position_dodge(width = 0.75)) 

p$layers <- p$layers[-1]
p
dev.off()

## Ordination
data.species.filter.nmds <-  ordinate(data.species.rel.filter, method="NMDS", "bray",k=3, trymax=500)

pdf("plots/sample_ordination.pdf", width=12, height=10)
plot_ordination(data.species.rel.filter, data.species.filter.nmds , color="treatment", shape="location")+
  geom_point(size=6)+
  theme_bw()+
  scale_color_npg()#+
  #geom_label(label=sample_names(data.species.rel.filter))
dev.off()

## if phylogeny included:
# pdf("plots/sample_ordination_PCOA_unifrac_uw.pdf", width=12, height=10)

# unifrac.dist <- UniFrac(data.species.rel.filter)
# ordi = ordinate(data.species.rel.filter, "PCoA", "unifrac", weighted=F)
# plot_ordination(data.species.rel.filter, ordi, color="Plot")
# dev.off()

## Networks (replace id by group if you want to have them merged by metadata)
netmat <- t(otu_table(data.species.rel.filter))
#netmat <- otu_table(merge_samples(data.species.rel.filter,group="id"))

### (optional) keep only major links, otherwise it often becomes overwhelming
otu_table(netmat)[otu_table(netmat)<0.05 ]<-0
netmat		= prune_taxa(taxa_sums(netmat)>0, netmat)

### plotting
pdf("plots/sample_network.pdf", width=25, height=10)
plotweb(data.frame(t(otu_table(netmat))))
dev.off()
