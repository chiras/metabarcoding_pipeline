### pseudo fill metadata
fill_pseudo_metadata <- function(phyloseq){
  if (sum(!is.na(sample_data(phyloseq)$location)) == 0){
    sample_data(phyloseq)$location <- sample(c("DummyLoc1","DummyLoc2","DummyLoc3"),length(sample_data(phyloseq)$location), replace=T)
  }

  if (sum(!is.na(sample_data(phyloseq)$treatment)) == 0){
    sample_data(phyloseq)$treatment <- sample(c("DummyTreat1","DummyTreat2","DummyTreat3"),length(sample_data(phyloseq)$treatment), replace=T)
  }
  return(phyloseq)
  
}
  
tax_glom_species_filtered <- function(physeq, rank = "species", sp_pattern = "_sp| spc| sp\\.| spc\\.") {
  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(physeq))
  
  # Identify placeholder species
  is_placeholder <- str_detect(tax_df[[rank]], regex(sp_pattern, ignore_case = TRUE))
  
  # Enumerate placeholders by genus
  tax_df_enumerated <- tax_df
  if ("genus" %in% colnames(tax_df)) {
    tax_df_enumerated[is_placeholder, rank] <- tax_df[is_placeholder, ] %>%
      group_by(genus) %>%
      mutate(
        !!rank := paste0(!!sym(rank), " ", row_number())
      ) %>%
      pull(!!rank)
  } else {
    stop("Taxonomy table must contain a 'genus' column to enumerate.")
  }
  
  # Replace the taxonomy in the phyloseq object
  tax_table(physeq) <- tax_table(as.matrix(tax_df_enumerated))
  
  # Redo identification now that names are unique
  tax_df <- tax_df_enumerated
  is_placeholder <- str_detect(tax_df[[rank]], regex(sp_pattern, ignore_case = TRUE))
  
  # Split physeq
  physeq_agglom <- prune_taxa(!is_placeholder, physeq)
  physeq_leave <- prune_taxa(is_placeholder, physeq)
  
  # Agglomerate only non-placeholder species
  physeq_agglommed <- tax_glom(physeq_agglom, taxrank = rank)
  
  # Merge both
  merged <- merge_phyloseq(physeq_agglommed, physeq_leave)
  
  return(merged)
}

### define taxa filtering
remove_unresolved_taxa <- function(phyloseq){
  if (marker == "ITS2"){
    (phyloseq = subset_taxa(phyloseq, phylum=="p:Streptophyta"))
    (phyloseq = subset_taxa(phyloseq, species!="p:Streptophyta_spc_spc_spc_spc"))
    (phyloseq = subset_taxa(phyloseq, species!="k:Viridiplantae_spc_spc_spc_spc_spc"))
    (phyloseq = subset_taxa(phyloseq, kingdom!=""))
  }
  
  if (marker == "COI"){
    (phyloseq = subset_taxa(phyloseq, phylum!="k:Animalia_spc"))
    (phyloseq = subset_taxa(phyloseq, class!="p:Insecta_spc"))
    (phyloseq = subset_taxa(phyloseq, class!="p:Arachnida_spc"))
    (phyloseq = subset_taxa(phyloseq, class!="p:Gastropoda_spc"))
  }
  
  if (marker == "16S"){
    (phyloseq = subset_taxa(phyloseq, kingdom=="d:Bacteria"))
    (phyloseq = subset_taxa(phyloseq, phylum!="p:Cyanobacteria/Chloroplast"))
    (phyloseq = subset_taxa(phyloseq, genus!="d:Bacteria_spc_spc_spc_spc"))
    (phyloseq = subset_taxa(phyloseq, order!="f:Mitochondria"))
    (phyloseq = subset_taxa(phyloseq, kingdom!=""))
  }
  
  return(phyloseq)
}

### uneven rarefaction_subsample
rarefaction_subsample <- function(x, sample.size, replace=FALSE){
  rarvec <- numeric(length(x))
  if(sum(x) <= 0){
    # Protect against, and quickly return an empty vector,
    # if x is already an empty count vector
    return(rarvec)
  }
  if(replace){
    # resample with replacement
    suppressWarnings(subsample <- sample(1:length(x), sample.size, replace=TRUE, prob=x))
  } else {
    # resample without replacement
    obsvec <- apply(data.frame(OTUi=1:length(x), times=x), 1, function(x){
      rep_len(x["OTUi"], x["times"])
    })
    obsvec <- unlist(obsvec, use.names=FALSE)
    # use `sample` for subsampling. Hope that obsvec doesn't overflow.
    suppressWarnings(subsample <- sample(obsvec, sample.size, replace=FALSE))
  }
  # Tabulate the results (these are already named by the order in `x`)
  sstab <- table(subsample)
  # Assign the tabulated random subsample values to the species vector
  rarvec[as(names(sstab), "integer")] <- sstab
  # Return abundance vector. Let replacement happen elsewhere.
  return(rarvec)
}

#### get species country occurences
# curl -i -H "Content-Type: application/json" -H "Accept: application/json" -X POST -d @filter.json  https://api.gbif.org/v1/occurrence/search?year=1800,1899
# wget https://api.gbif.org/v1/occurrence/download/request/0056302-230530130749713.zip # All
# wget https://api.gbif.org/v1/occurrence/download/request/0096050-230530130749713.zip # Germany
# wget https://api.gbif.org/v1/occurrence/download/request/0096276-230530130749713.zip # Austria
# cat 0*.csv -> combined.csv
# cut -f 10,16 combined.csv > species_country.csv
# cat species_country.csv | sort | uniq -c | sed -e "s/^[[:space:]]*//" -e "s/[[:space:]]/  /" > species_country.sorted.csv
# cat species_country.sorted.csv | sed -e "/            /d" -e 's/"//g' -e "/^1 $/d" -e "s/'//g" > species_country.sorted2.csv
# cat -n species_country.sorted2.csv | sed -e "s/ × / x /g" -e "s/^[[:space:]]*/T/" -e "/       $/d"> species_country.sorted3.csv
# cp species_country.sorted2.csv ../../_resources/GBIF_latest.csv

prepare_species_occurences <- function(){
  gbif_records <- read.table("/Users/ra39huv/TMP/Data_processing/_resources/GBIF_20230724.csv", sep="\t", row.names=1, fill=T)
  colnames(gbif_records) <- c("value","SciName","Country")
  gbif_records <- gbif_records[gbif_records$Country !="countryCode",]
  gbif_records <- gbif_records[gbif_records$SciName !="",]
  gbif_wide <- reshape(gbif_records, idvar = "SciName", timevar = "Country", direction = "wide")
  gbif_wide <- gbif_wide[is.na(gbif_wide$SciName)==F,]
  row.names(gbif_wide) <- gbif_wide$SciName
  return(gbif_wide)
}
  
get_species_occurences_wide <- function(phyloseq){
    gbif_wide <- prepare_species_occurences()
    taxa_abundance <- as.data.frame(sort(rowSums(otu_table(phyloseq)), decreasing=T))
    taxa_matches <- cbind(row.names(taxa_abundance),taxa_abundance,gbif_wide[rownames(taxa_abundance),])
    colnames(taxa_matches)[1:3] <- c("Species","Abundance","GBIF.match")        
    return(taxa_matches)
}

get_country_codes<- function(){
  cc <- read.table("/Users/ra39huv/TMP/Data_processing/_resources/GBIF_countrycodes.csv", sep=",", fill=T, header=T, row.names=2)
  rownames(cc)[rownames(cc)=="Na"]<-"NA"
  return(cc)
}

get_species_occurences_long <- function(phyloseq, notInGBIF = T){
  taxa_matches <- get_species_occurences_wide(phyloseq)
  taxa_matches_long <- melt(taxa_matches, id.vars=c("Species","Abundance","GBIF.match"), variable.name = "Country")
  taxa_matches_long$Country <- gsub("value.","",taxa_matches_long$Country)
  taxa_matches_long=taxa_matches_long[is.na(taxa_matches_long$value)==F,]

  cc <- get_country_codes()
  
  taxa_matches_long2 <- cbind(taxa_matches_long,cc[taxa_matches_long$Country, c("name","region","sub.region")])
  taxa_matches_long2 <- taxa_matches_long2[order(taxa_matches_long2$region),]
  
  if (notInGBIF){
    not.in.gbif <- data.frame(Species=as.character(data.frame(tax_table(phyloseq))[!(taxa_names(phyloseq) %in% taxa_matches_long2$Species),"species"]), Abundance=as.numeric(rowSums(otu_table(phyloseq)[!(taxa_names(phyloseq) %in% taxa_matches_long2$Species),])), GBIF.match=NA, Country=NA, value=10, name=NA , region = NA,sub.region=NA)
    taxa_matches_long2 <- rbind(taxa_matches_long2,not.in.gbif)
  }
  
  return(taxa_matches_long2)
}




#### Custom functions
#Propagate taxonomie for tax_glom to work
propagate_incomplete_taxonomy <- function(phyloseq){
	taxranks <- colnames(tax_table(phyloseq))
	for (i in 2:length(taxranks)){
	  if (sum(sum(tax_table(phyloseq)[,taxranks[i]]==""))>1){
		tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i]]<-paste(tax_table(phyloseq)[tax_table(phyloseq)[,taxranks[i]]=="",taxranks[i-1]],"_spc",sep="")
	  }
	}
	  return(phyloseq)
}



#Make taxa labels nice for plots
replace_tax_prefixes <- function(phyloseq){
  tmp_tax_table <- apply(tax_table(phyloseq), c(1, 2),function(y) gsub("^\\w:","",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_spc_.*","_spc",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub("_"," ",y))
  tmp_tax_table <- apply(tmp_tax_table, c(1, 2),function(y) gsub(";$","",y))
  tax_table(phyloseq)<- tmp_tax_table
  return(phyloseq)
}

#label low throughput samples
label_low_throughput <- function(phyloseq, threshold){
  sample_names(phyloseq)[sample_sums(phyloseq)<threshold]<-paste(sample_names(phyloseq)[sample_sums(phyloseq)<threshold],"|LT",threshold,sep="")
  return(phyloseq)
}

label_sample_by_host <- function(phyloseq, hostcol, projcol="", idcol=""){
  hostlist <- sample_data(phyloseq)[[hostcol]]
  projlist <- sample_data(phyloseq)[[projcol]]
  print(length(idcol))
  if(length(idcol)!=length(hostlist)){
   	uniqIDs <- sample(1:99999,length(hostlist) ,replace=F)
   	uniqIDs <- sprintf(paste(hostlist, projlist,"%05d", sep="|"), uniqIDs)
  }else{
  	uniqIDs <- sprintf(paste(hostlist, projlist,"%05d", sep="|"), idcol)
  }

  sample_names(phyloseq)<-uniqIDs
  return(phyloseq)
}

id_cont_lh <-function(phyloseq, distri=T, neg, pos){
  neg_samples=c()
  if (length(neg)>0){
    neg_samples= subset_samples(test, sample_names(test) %in% neg)
    neg_samples.rel = transform_sample_counts(neg_samples, function(x) x/sum(x))
    sam_samples= subset_samples(test, !(sample_names(test) %in% neg))
    sam_samples.rel = transform_sample_counts(sam_samples, function(x) x/sum(x))

        tail(sort(rowSums(otu_table(neg_samples))))
  }
  pos_samples=c()
  if (length(pos)>0){
    pos_samples=pos
  }


  ## checks
  ### distribution of continuous distribution over majority of samples
  rowSums(otu_table(test)>0)
  rowSums(otu_table(test))

  d1 <- ((rowSums(otu_table(sam_samples)>0))/length(otu_table(sam_samples)[1,]))
  d2 <- rowSums(otu_table(sam_samples.rel))
  d3 <- rowSums(otu_table(sam_samples.rel))/length(otu_table(sam_samples)[1,])

  d1n <- ((rowSums(otu_table(neg_samples)>0))/length(otu_table(neg_samples)[1,]))
  d2n <- rowSums(otu_table(neg_samples.rel))
  d3n <- rowSums(otu_table(neg_samples.rel))/length(otu_table(neg_samples.rel)[1,])

  #
  # names(d1)[d1-d1n<0]
  #
  # sort(d1-d1n)
  #
  # plot(d1n~d1)
  # plot(d2n~d2)

    plot(d1~d2 ) # i would expect contamination to be in many samples, but low abundances, i.e. high y and low x (upper left corner)
    plot(d3n~d3) # i would expect contamination to be low in d2 and high in d2n, NOT cross-contamination though between samples

    # points(d2n~d1n ,add=T, col="red")



  #### relation to low abundance samples

   abu_h <- subset_samples(sam_samples.rel, colSums(otu_table(sam_samples)) > 2500)
   abu_l <- subset_samples(sam_samples.rel, colSums(otu_table(sam_samples)) < 2500)

    plot(rowMeans(otu_table(abu_h))~rowMeans(otu_table(abu_l)), pch=19, alpha=0.3)
      # I would expect low throughput samples to have higher values for contamination taxa

    tail(sort(round(rowMeans(otu_table(abu_h))-rowMeans(otu_table(abu_l)), digits=4)))
    par(mar=c(4,15,1,1), mfrow=c(1,1))

    barplot(head(sort(round(rowMeans(otu_table(sam_samples.rel))-rowMeans(otu_table(neg_samples.rel)), digits=4)), n=20), horiz=T, las=2)
      # i would expect contamination to occur more in negative controls


        # plot(c(otu_table(its2.species.rel)[10,])~ c(colSums(otu_table(its2.species))), type="n", xlim=c(0,150000), ylim=c(0,1))

        # require(RColorBrewer)
        # col_vector = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],length(otu_table(its2.species.rel)[,1]), replace=T)
        # #
    # for (i in 1:length(otu_table(its2.species.rel)[,1])){
    #   points(c(otu_table(its2.species.rel)[i,])~ c(colSums(otu_table(its2.species))),col=col_vector[i], add=T)
    # }
    #

    # check lt flag in name

  ### presence in negative controls
  ### presence in positive controls


}


map_interactions_trait <- function(phyloseq, traittable, traitcols="", speccol){
	otu <- otu_table(phyloseq)
	taxa <- taxa_names(phyloseq)
	samples <- sample_names(phyloseq)

	#traittable <- chemie
	#speccol <-"host"
	#traitcols = c(5:31)

	traits <- traittable[, traitcols]
	#traits$spc <- traittable[, speccol]

	tmp_matrix = matrix(ncol=length(traits[1,]), nrow=length(taxa), data=NA, byrow=F)
	rownames(tmp_matrix)<- taxa;		colnames(tmp_matrix)<- names(traits)

	results = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results)<- names(traits)
	results_uncertainity = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results_uncertainity)<- names(traits)
	results_uncertainity_weighted = data.frame(matrix(ncol=length(traits[1,]), nrow=0))
	colnames(results_uncertainity_weighted)<- names(traits)

	results_nomatch = c()

	for (i in 1:length(taxa)){
		matches <- traittable[, speccol] == taxa[i]
		if (sum(matches) == 1){
			tmp_matrix[taxa[i],] <- as.matrix(traits[matches,])#as.matrix(traits[matches,1:length(traits[1,])-1])
		}else if(sum(matches) > 1){
			# not sure what to do is multiple hits
		}else {
			results_nomatch =c(results_nomatch, taxa[i])
		}
	}
	template_matrix <- tmp_matrix

	for (j in 1:length(samples)){
		conc_matrix = template_matrix

#		if (sum(row.names(otu[,j])!=row.names(conc_matrix))>1){
#			print "Taxa inconsistency"
#		}else{
			prob_matrix <- matrix(ncol=length(tmp_matrix[1,]) , nrow=length(tmp_matrix[,1])  , data=otu[,j], byrow=F)
			mult_matrix <- conc_matrix * prob_matrix
			tmp_results <- as.data.frame(t(colSums	(mult_matrix, na.rm=T)))
			tmp_uncertainity <- as.data.frame(t(colSums	(is.na(mult_matrix)))/length(taxa))
			tmp_uncertainity_weighted <- colSums(as.numeric(!is.na(mult_matrix)) * prob_matrix)
			rownames(tmp_results) <- samples[j]
			rownames(tmp_uncertainity) <- samples[j]
			results = rbind(results, tmp_results)
			results_uncertainity = rbind(results_uncertainity, tmp_uncertainity)
			results_uncertainity_weighted = rbind(results_uncertainity_weighted, tmp_uncertainity_weighted)
#		} # end if taxa inconsistency


	} # end for j
	rownames(results_uncertainity_weighted)<- samples

	result_list=list()
	result_list$mapping <- results
	result_list$uncertainity <- results_uncertainity
	result_list$uncertainity_weighted <- results_uncertainity_weighted
	result_list$nomatch <- results_nomatch

  return(result_list)
}
