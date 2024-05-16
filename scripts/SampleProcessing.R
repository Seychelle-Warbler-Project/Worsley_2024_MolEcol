setwd("/users/sarah/OneDrive - University of East Anglia/Seychelles_warbler/Analysis4_toS2022/MergedQIIME/Age_paper/Age/paper/FinalDraft/MolEcol/revisions/DRYAD/SampleProcessing")

#load packages
  library(phyloseq)
  library(ape)
  library(ggplot2)
  library(plyr)
  library(decontam)
  library(vegan)
  library(FSA)
  library(microbiome)
  library(tidyverse)

#ASV table
  asv_table <- read.csv ("ASV_table.csv", row.names=1) #read in asv table with feature names as rownames
  str (asv_table) #should be 62598 features, 1015 samples
  head (asv_table)

  asv_table <- as.matrix (asv_table) #make into a matrix

#taxonomy
  taxonomy <- read.csv ("Taxonomy.csv", row.names=1)
  str (taxonomy) #should be 62598 observations, 7 taxonomic groupings, feature names as row names.
  taxonomy <- as.matrix (taxonomy)

#load filtered metadata with sample names as rownames (1015 samples)
  metadata<-read.csv("Metadata.csv", row.names = 1)
  str(metadata)  
  head(metadata)
  
#read in tree as a phyloseq object
  phy_tree <- read_tree ("tree.nwk")


#### MERGE INTO PHYLOSEQ OBJECT ####
#import all as phyloseq objects
  ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxonomy)
  META <- sample_data(metadata)
  head(META)

#check that the ASV and sample names are consistent across objects (e.g have dashes been replaced with dots?)
  str(taxa_names(TAX))
  str(taxa_names(ASV))
  str(taxa_names(phy_tree))

  str(sample_names(ASV))
  str(sample_names(META))

#merge into phyloseq object
  physeq <- phyloseq(ASV, TAX, META, phy_tree)
  physeq #check there are the right numbers of samples/features/variables etc... 

  table(sample_data(physeq)$SampleType4Way)

### 969 faecal samples
### 21 Collection controls
### 15 extraction controls
### 10 Positive controls - zymo biomics

#check the number of reads per sample
  sample_reads<-data.frame(reads=sample_sums(physeq))
  sample_reads$ID<-rownames(sample_reads)

  head(sample_reads)
  sum(sample_reads$reads) # 81911398 reads in total across all samples
  mean(sample_reads$reads) #80700.88
  sd(sample_reads$reads) # 62600.86

  
  
##################
### FILTERING ####
##################
  
### filter to remove instances where features are not assigned as bacteria/ unassigned ###

  table(tax_table(physeq)[,"Kingdom"]) 
  
  # removes all labelled as unassigned at kingdom level (73 ASVs)
  physeq2<-subset_taxa(physeq, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned")) 
  table(tax_table(physeq2)[,"Kingdom"])
  
  #remove those features labelled as Archaea (439 ASVs)
  physeq3<-subset_taxa(physeq2, !Kingdom %in% "Archaea") 
  table(tax_table(physeq3)[,"Kingdom"]) 
  
 
### check the positive controls ###
  
  str(metadata)
  PosControls<- subset_samples(physeq3, SampleType4Way == "PositiveControl")
  PosControlsASVs<- data.frame(otu_table(PosControls))
  head(PosControlsASVs)
  PosControlsReads<- subset(PosControlsASVs, rowSums(PosControlsASVs, na.rm = TRUE) > 0)
  PosControlsReads$TotalDetected<- rowSums(PosControlsReads)
  PosControlsTaxonomy<- merge(PosControlsReads, taxonomy, by="row.names")
  head(PosControlsTaxonomy)
  
  PosControlsTaxonomy %>%
    group_by(Family,Genus) %>%
    summarise(TotalReads = sum(TotalDetected))
  
  # There are 8 Genera that have high number of reads in pos controls
  # (remaining 6 have <70 reads in total each so probably errors)
  # The 8 genera are bacillus, Escherichia-Shigella, Enterococcus,
  # Lactobacillus, Listeria, Pseudomonas, Staphylococcus and one unknown Enterobacteriaceae
  # If you blast the unclassified Enterobacteriaceae- this is Salmonella
  # These match the 8 species in the positive zymo control.
 

### remove the positive controls ###
  
  physeq4<- subset_samples(physeq3, SampleType4Way != "PositiveControl") 
  table(sample_data(physeq4)$SampleType4Way)
  

### remove ASVs that are unassigned at phylum level
  
  #print features per phylum- can see that 797 are phylum unassigned (have a blank name)
  table(tax_table(physeq4)[,"Phylum"])
  # removes all labelled as unassigned at phylum level
  physeq5<-subset_taxa(physeq4, !is.na(Phylum) & !Phylum %in% c("", "Unassigned"))
  table(tax_table(physeq5)[,"Phylum"]) #check removed- gives the number of features per phylum
  
  #check adjusted numbers - now 61289 ASVs, 1005 samples
  physeq5
  table(sample_data(physeq5)$SampleType4Way)
  

# Check if any samples have 0 reads following removal of ASVs unassigned to phylum
  sample_data(physeq5)$SampleReads<- sample_sums(physeq5)
  DataPhyseq<- data.frame(sample_data(physeq5))
  summary(DataPhyseq$SampleReads) #min is 3 reads
  
  
################
### DECONTAM ###
################
  
## plot read counts per sample
  df <- as.data.frame(sample_data(physeq5)) # Put sample_data into a data.frame
  df$LibrarySize <- sample_sums(physeq5)
  df <- df[order(df$LibrarySize),]
  df$Index <- seq(nrow(df))
  ggplot(data=df, aes(x=Index, y=LibrarySize, color=SampleType4Way)) + geom_point(size=0.4)
  #extraction controls fall at lower end of library size scale
  #collection controls are distributed throughout
  
## Filter based on prevalence ##
  
# 1) Filter based on ASVs in extraction controls (there are 15 of these- approx 2 per kit)
  
  physeq5 # 61289 ASVs
  sample_data(physeq5)$is.negLab <- ifelse(sample_data(physeq5)$SampleType4Way == "ExtractionControl", TRUE, FALSE)
  
  decontamdf.prev.Lab <- isContaminant(physeq5, method="prevalence", neg="is.negLab")
  table(decontamdf.prev.Lab$contaminant) # 32 contaminants identified from extraction controls
  ContamsLabControls<- decontamdf.prev.Lab[decontamdf.prev.Lab$contaminant=="TRUE",]
  
  ContamsLabControls<- merge(ContamsLabControls, taxonomy, by="row.names")
  str(ContamsLabControls)
  # write.csv(ContamsLabControls, "ContaminantsLabControls.csv") 
  # Lachnospiraceae particularly prevalent- known to be part of the kitome (https://bmcmicrobiol.biomedcentral.com/articles/10.1186/s12866-020-01839-y)
  
  # Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa <- transform_sample_counts(physeq5, function(abund) 1*(abund>0))
  ps.pa.Lab <- prune_samples(sample_data(ps.pa)$is.negLab== "TRUE", ps.pa)
  ps.pa.pos <- prune_samples(sample_data(ps.pa)$is.negLab== "FALSE", ps.pa)
  
  # Make data.frame of prevalence in positive and negative samples
  df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.Lab),
                      contaminant=decontamdf.prev.Lab$contaminant)
  ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
  #Remove contaminating taxa
  taxaNames<- taxa_names(physeq5)
  str(taxaNames)
  IDsLabContam<- as.vector(ContamsLabControls$Row.names)
  taxaNames2<- taxaNames[!taxaNames %in% IDsLabContam]
  str(taxaNames2)
  physeq6<-prune_taxa(taxaNames2, physeq5) 
  physeq6 # 61257 asvs
  
  summary(sample_sums(physeq6))
  
  
  
# 2) Remove ASVs that are collection contaminants (21 collection controls)
  
  sample_data(physeq6)$is.negCollection <- sample_data(physeq6)$SampleType4Way == "CollectionControl"
  
  decontamdf.prev.Collection <- isContaminant(physeq6, method="prevalence", neg="is.negCollection")
  table(decontamdf.prev.Collection$contaminant) #5983 contaminants
  ContamsCollectionControls<- decontamdf.prev.Collection[decontamdf.prev.Collection$contaminant=="TRUE",]
  
  ContamsCollectionControls<- merge(ContamsCollectionControls, taxonomy, by="row.names")
  str(ContamsCollectionControls)
  ggplot(ContamsCollectionControls,aes(Phylum, prev)) + geom_point() +
    theme(axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1))
  #write.csv(ContamsCollectionControls, "ContaminantsCollectionControls.csv")
  
  # Make phyloseq object of presence-absence in negative controls and true samples
  ps.pa.Collection1 <- transform_sample_counts(physeq6, function(abund) 1*(abund>0))
  ps.pa.Collection2 <- prune_samples(sample_data(ps.pa.Collection1)$is.negCollection== "TRUE", ps.pa.Collection1)
  ps.pa.posCollection <- prune_samples(sample_data(ps.pa.Collection1)$is.negLab== "FALSE", ps.pa.Collection1)
  
  # Make data.frame of prevalence (number of samples taxa appear in) for positive and negative samples
  df.paCollection <- data.frame(pa.pos=taxa_sums(ps.pa.posCollection), pa.neg=taxa_sums(ps.pa.Collection2),
                                contaminant=decontamdf.prev.Collection$contaminant)
  ggplot(data=df.paCollection, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
    xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
  
  #Remove contaminating taxa
  
  ContamsCollectionToFilter<- ContamsCollectionControls
  str(ContamsCollectionToFilter) # 5983
  taxaNamesAll<- taxa_names(physeq6)
  str(taxaNamesAll)#61257
  IDsCollectionContam<- as.vector(ContamsCollectionToFilter$Row.names)
  length(IDsCollectionContam)#5983
  taxaNamesToKeep<- taxaNamesAll[!taxaNamesAll %in% IDsCollectionContam]
  str(taxaNamesToKeep)
  physeq7<-prune_taxa(taxaNamesToKeep, physeq6) 
  physeq7 # now 55274 ASVs remaining
  
  
  #check if any samples with 0 reads following decontam
  sample_data(physeq7)$SampleReads<- sample_sums(physeq7)
  DataPhyseq<- data.frame(sample_data(physeq7))
  ZeroCounts<- DataPhyseq[DataPhyseq$SampleReads==0,]
  ZeroCounts #0 samples now has zero reads 

  table(sample_data(physeq7)$SampleType4Way)
  
  
#### remove control samples ###
  
  physeq8<- subset_samples(physeq7, SampleType4Way == "Faecal") 
  physeq8 <- prune_taxa(taxa_sums(physeq8) > 0,physeq8) # get rid of ASVs that have 0 reads after removing controls
  physeq8 # 51360 ASVs, 969 faecal samples
  table(sample_data(physeq8)$SampleType4Way)
  
  
  
########################## 
##### Filter by reads ####
##########################
  
#### SAMPLE COMPLETENESS using iNEXT ####

#make the ASV abundance table into a matrix and check by printing first 2 rows/columns
  abund <- as(otu_table(physeq8), "matrix") 
  abund[1:2,1:2]
  
#convert to a dataframe
  abund2 <- as.data.frame(abund)
  str(abund2)
  
#iNEXT only takes numeric values, so change all values in the dataframe to numeric values instead of integers.
  df2 <- mutate_all(abund2, function(x) as.numeric(x)) 
  str(df2)
  
# Install and load inext
  library(iNEXT)
  #iNEXT= https://cran.r-project.org/web/packages/iNEXT/iNEXT.pdf
  
#Run inext
  #q=0 specifies that the function should use species richness (rather than shannon (1) or simpson (2) indices) for the rarefaction. 
  #Datatype= abundance because you have raw abundances. 
  #endpoint=20000 specifies the number of reads (sample size) that should be used as an endpoint for the rarefaction/extrapolation. 
  
  inext<-iNEXT(df2, q=0, datatype="abundance", endpoint=20000)
  dput(inext, "inextResults.txt")
  
  inextresults<- dget("inextResults.txt")
  
  #plot rarefaction curve
  rarefaction<- ggiNEXT(inextresults, type=1, se=TRUE, grey= TRUE) + theme(legend.position = "none")+ xlab("Sequencing depth") + ylab("Observed ASVs")
  rarefaction
  
  #Plot sample completeness curve
  completeness<-ggiNEXT(inextresults, type=2, se=TRUE, grey=TRUE)+scale_shape_manual(values=rep(20,164))+ theme(legend.position = "none") +xlab("Read count") +ylab("Sample completeness")
  completeness + geom_vline(xintercept=10000, alpha=0.5, linetype=2)
  # completeness begins to plateau around 7000-10000 reads
  
  
  
## Identify samples with <8'000 reads
  SampleReadsphyseq8<- data.frame(sample_sums(physeq8))
  SampleReadsphyseq8$Sample.ID<- row.names(SampleReadsphyseq8) 
  colnames(SampleReadsphyseq8)<- c("reads","sample")
  SampleReadsphyseq8 <-  SampleReadsphyseq8[order(SampleReadsphyseq8$reads),]
  SampleReadsphyseq8[c(1:50),]
  samplesLowReads<- SampleReadsphyseq8[SampleReadsphyseq8$reads <8000,]
  samplesLowReadsID<- samplesLowReads$sample 
  length(samplesLowReadsID)# 27 faecal samples have <8000
  
  physeq8Data<- data.frame(sample_data(physeq8))
  physeq8DataLowReads<- physeq8Data[physeq8Data$sample.id %in% samplesLowReadsID,] # almost all samples had very little material/low qubit
  str(physeq8DataLowReads)

  
### remove samples with less than 8000 reads 
  physeq9<-prune_samples(sample_sums(physeq8)>=8000, physeq8)
  physeq9 <- prune_taxa(taxa_sums(physeq9) > 0,physeq9) # get rid of ASVs that have 0 reads after removing samples
  physeq9 # removes 27 samples - 942 samples, 51267 taxa
  
  
### Filter out ASVs with low number of reads ###

# Plot distribution of total reads for each ASV

  # find max abundance for each ASV across samples
  
  maxAb<- data.frame(otu_table(physeq9))
  maxAb$max <- do.call(pmax, maxAb[1:936])
  maxAb2<- maxAb[maxAb$max >0 & maxAb$max <100,]
  hist(maxAb2$max, breaks=50)
  ticks<-seq(from=0, to=100, by=10)
  axis(1, at=ticks)
  # a lot of low frequency ASVs
  # There were also several spurious ASVs in the positive mock controls with <70 reads
  # remove ASVs with <50 reads in total across all samples

### Filter ASVs with <50 reads in total across all samples
  
  # Compute prevalence of each feature (total number of samples in which a taxon appears at least once)
  prevdf = apply(X = otu_table(physeq9),
                 MARGIN = ifelse(taxa_are_rows(physeq9), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeq9),
                      tax_table(physeq9))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
  # Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq9, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / phyloseq::nsamples(physeq9),color=Phylum)) +
    # Include threshold: total abundance across all samples set to 50
    geom_vline(xintercept = 50, alpha = 0.5, linetype = 2)+  geom_point(size = 0.3, alpha = 0.7) +
    geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2)+
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
# ASVs with <50 reads in total tend to be found in <1% of samples (~9 samples)
  physeq9

# Define abundance threshold as 50 total reads across samples
  
  abundanceThreshold<-50
  
# Execute the filter, using `prune_taxa()` function
  
  head(prevdf1)
  
  KeepTaxa<-prevdf1[prevdf1$TotalAbundance>=abundanceThreshold,]
  str(KeepTaxa)
  head(KeepTaxa)
  KeepTaxaNames<- rownames(KeepTaxa)
  
  physeq10<- prune_taxa(KeepTaxaNames, physeq9)
  physeq10 # 23151 taxa when 50 used as cutoff (51267 ASVs before filtering)
  
  sum(sample_sums(physeq9)) # 68841738
  sum(sample_sums(physeq10)) # 68343354
  (68343354/68841738)*100 # 99.27% of reads retained
  

  
#Find mean number of ASVs across samples prior to rarefaction
  
  range(sample_sums(physeq10)) #range of read numbers across samples: 8166-962984
  observationThreshold <- 1
  sumASVs<- (apply(otu_table(physeq10) > observationThreshold, 2, sum))
  str(sumASVs)
  sumASVs <- data.frame(sumASVs)
  str(sumASVs)
  mean(sumASVs[,1]) #mean ASVs per sample = 225.621
  sd(sumASVs[,1]) #Standard deviation of ASVs per sample = 162.484
  range(sumASVs[,1]) # range 5-933
  print(sd(sumASVs[,1])/sqrt(length((sumASVs[,1])))) # SE = 5.294
  
  ReadsAfterFiltering<- sample_sums(physeq10)
  head(ReadsAfterFiltering)
  length(ReadsAfterFiltering)
  range(ReadsAfterFiltering) 
  mean(ReadsAfterFiltering) # 72551.33 (mean reads per sample)
  sd(ReadsAfterFiltering) #55495.42 (sd reads)
  
  
### extract the filtered taxonomy,asv and metadata files for other downstream analysis
  
  filteredMeta<- data.frame(sample_data(physeq10))
  write.csv(filteredMeta, "FilteredMetadata.csv")
  
  filteredASV<- data.frame(otu_table(physeq10))
  write.csv(filteredASV, "FilteredASVTable.csv")
  
  filteredTax<- data.frame(tax_table(physeq10))
  write.csv(filteredTax, "FilteredTaxonomyTable.csv")
  
  
#######################
#### Repeatability ####
#######################

### 1) alpha diversity measures ###

#ASV table
  asvFiltered <- read.csv ("FilteredASVTable.csv", row.names=1) #read in asv table with feature names as rownames
  str (asvFiltered) #should be 942 samples, 23151
  head (asvFiltered)
  
  asvFiltered_table <- as.matrix (asvFiltered) #make into a matrix
  
#taxonomy
  taxonomyFiltered <- read.csv ("FilteredTaxonomyTable.csv", row.names=1)
  str (taxonomyFiltered) #should be 23151 observations, 7 taxonomic groupings, feature names as row names.
  taxonomyFiltered_matrix <- as.matrix (taxonomyFiltered)
  
#load filtered metadata with sample names as rownames
  metadataFiltered<-read.csv("FilteredMetadata.csv", row.names = 1)
  str(metadataFiltered)  
  head(metadataFiltered)
  
#read in tree as a phyloseq object
  phy_tree <- read_tree ("tree.nwk")
  
  
#### MERGE INTO PHYLOSEQ OBJECT ####
  #import all as phyloseq objects
  ASV_Filtered <- otu_table(asvFiltered_table, taxa_are_rows = TRUE)
  TAX_Filtered <- tax_table(taxonomyFiltered_matrix)
  META_Filtered <- sample_data(metadataFiltered)

#merge into phyloseq object
  physeqRepeatability <- phyloseq(ASV_Filtered, TAX_Filtered, META_Filtered, phy_tree)
  physeqRepeatability  #23151 ASVs, 942 samples
  
  
  
# RAREFY READS TO MIN SAMPLING DEPTH: 8,000 reads #
  
  #rarefy to 8000 and set seed before rarefying (58392), so that results are reproducible  
  physeqRare<-rarefy_even_depth(physeqRepeatability, 8000, rngseed = 58392)
  
  # 73 ASVs removed after subsampling- leaves 23078 ASVs across 942 samples
  physeqRare
  sample_sums(physeqRare)
  
  
##### CALCULATE ALPHA DIVERSITY METRICS ####
  
  # calculate diversity metrics using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRare, split=TRUE, measures= c("Shannon", "Observed"))
  head(richnessEstRare)
  
  #add alpha diversity metrics to metadata
  physeqRareMeta <- data.frame(sample_data(physeqRare))
  head(physeqRareMeta)
  physeqRareMeta$Shannon <- richnessEstRare$Shannon
  physeqRareMeta$Observed <- richnessEstRare$Observed

  str(data.frame(physeqRareMeta))

  write.csv(physeqRareMeta, "MetaForMatrix.csv")
  
  
### IDENTIFY SAMPLES EXTRACTED TWICE ######
  
  MetaForMatrix<- read.csv("MetaForMatrix.csv", row.names = 1)
  head(MetaForMatrix)
  str(MetaForMatrix)
  
  #extract those sequenced twice
  row_namesSeqDup <- as.vector(grep("b", rownames(physeqRareMeta), value=TRUE)) 
  str(row_namesSeqDup) #69 seq reps
  TubesSeqDup<- MetaForMatrix[MetaForMatrix$sample.id %in% row_namesSeqDup,]
  TubesIDsSeqDup<- TubesSeqDup$TubeNumber #tube numbers of duplicates (should match in both duplicates)
  
  #extract those where faecal sample was extracted and sequenced twice
  row_namesExDup <- as.vector(grep("R", rownames(physeqRareMeta), value=TRUE)) 
  str(row_namesExDup) #10 seq reps
  TubesExDup<- MetaForMatrix[MetaForMatrix$sample.id %in% row_namesExDup,]
  TubesIDsExDup<- TubesExDup$TubeNumber #tube numbers of duplicates (should match in both duplicates)
  
  TubeDupAll<- c(TubesIDsSeqDup, TubesIDsExDup)
  length(TubeDupAll)
  
  AllDups<- MetaForMatrix[MetaForMatrix$TubeNumber %in% TubeDupAll,] #Filter metadata to only include duplicated tube numbers  
  sort(AllDups$TubeNumber) #all duplicates there for 79 tubes (158 samples total)
  length(AllDups$TubeNumber)
  #write.csv(AllDups, "SequenceDuplicates.csv")
  
  
### CONSISTENCY OF ALPHA METRICS ###
  
#Make a distance matrix with euclidean distances of Shannon diversity
  str(AllDups)
  Shannon <- AllDups$Shannon
  samples <- AllDups$TubeNumberUnique
  names(Shannon) <- samples
  distMatrix <- vegdist(Shannon, method="euclidean")
  distMatrix <- as.matrix(distMatrix)
  distMatrix[1:3,1:3]
  
  
#extract as a dataframe (just upper triangle- one way comparisons)
  str(AllDups)
  pairwiseDist <- t(combn(colnames(distMatrix), 2))
  pairwiseDist<- data.frame(pairwiseDist, dist=distMatrix[pairwiseDist])
  head(pairwiseDist)
  colnames(pairwiseDist)<- c("ID1","ID2", "dist")
  str(pairwiseDist)
  
# Add in sequence run number
  
  SeqRun1<- AllDups[,c("TubeNumberUnique","SequencingRun")]
  head(SeqRun1)
  colnames(SeqRun1)<-c("ID1", "Run1")
  pairwiseDist<- merge(pairwiseDist, SeqRun1, by="ID1", all.x=T)
  head(pairwiseDist)
  
  SeqRun2<- AllDups[,c("TubeNumberUnique","SequencingRun")]
  head(SeqRun2)
  colnames(SeqRun2)<-c("ID2", "Run2")
  pairwiseDist<- merge(pairwiseDist, SeqRun2, by="ID2", all.x=T)
  head(pairwiseDist)

  
#Copy the two Id columns so that we can check for tube no. matches in seq dup: ID1seq, ID2seq
  pairwiseDist$ID1seq<- pairwiseDist$ID1
  pairwiseDist$ID2seq<- pairwiseDist$ID2
  str(pairwiseDist)
  
  #remove "b" from ID1seq and ID2seq to check for sequencing replicates (note, the second time sequenced has b after the tube number)
  pairwiseDist$ID1seq<-gsub("b","",as.character(pairwiseDist$ID1seq))
  pairwiseDist$ID2seq<-gsub("b","",as.character(pairwiseDist$ID2seq))
  str(pairwiseDist)
  
  
#Copy the two Id columns so that we can check for tube no. matches in extraction dup: ID1seq, ID2seq
  pairwiseDist$ID1extr<- pairwiseDist$ID1
  pairwiseDist$ID2extr<- pairwiseDist$ID2
  str(pairwiseDist)
  
  #remove "R" from ID1seq and ID2seq to check for sequencing replicates (note, the second time extracted has R after the tube number)
  pairwiseDist$ID1extr<-gsub("R","",as.character(pairwiseDist$ID1extr))
  pairwiseDist$ID2extr<-gsub("R","",as.character(pairwiseDist$ID2extr))
  str(pairwiseDist)
  
#name type of comparison
  pairwiseDist$Comparison<-ifelse(pairwiseDist$ID1seq==pairwiseDist$ID2seq, "SeqDuplicate",
                                            ifelse(pairwiseDist$ID1extr==pairwiseDist$ID2extr, "ExtractionDuplicate","DifferentSample"))
  str(pairwiseDist)
  
  table(pairwiseDist$Comparison) # 69 sequence duplicate comparisons, 10 extraction duplicate comparisons, 12324 between samples
  
  
#new column to say if duplicates of the same sample sequenced in the same or different run

  pairwiseDist$Comparison2<-ifelse(pairwiseDist$ID1seq==pairwiseDist$ID2seq,
                                             ifelse(pairwiseDist$Run1==pairwiseDist$Run2,
                                                    "DuplicateSameRun", "DuplicateDifferentRun"), pairwiseDist$Comparison)
  str(pairwiseDist)
  table(pairwiseDist$Comparison2) # different run=46, same run = 23 (for 69 sequenced twice)
  
#plot to check between versus within samples

  samplePairs<- c("Different samples", "Extraction duplicates", "Sequencing duplicates")
  AlphaPlot<- ggplot(pairwiseDist, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())
  AlphaPlot
  
#plot to check sequencing duplicates across same and different run
  pairwiseDist$Comparison2<- as.factor(pairwiseDist$Comparison2)
  levels(pairwiseDist$Comparison2)
  sampleCategory<- c("Different samples", "Duplicate different run", "Duplicate same run", "Extraction duplicate")
  AlphaPlot2<- ggplot(pairwiseDist, aes(Comparison2, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=14)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=sampleCategory) + theme(element_blank())
  AlphaPlot2
  
# Kruskal-Wallis test
  
  kruskal.test(dist ~ Comparison2, data = pairwiseDist)
  
  
#perform Dunn's Test with Bonferroni correction for p-values: 
# no sig difference for duplicates sequenced in same or different run in terms of shannon
# no diff between any of duplicates

  dunnAlpha<- dunnTest(dist ~ Comparison2,
           data=pairwiseDist,
           method="bonferroni")
  
  dput(dunnAlpha$res, "DunnAlphaRepeatability.txt")

  
  
### CONSISTENCY OF BETA METRICS ###
  
  physeqRepeatability #use the unrarefied reads
  
# filter rare taxa from phyloseq object (those with low prevalence)
  
  # Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(physeqRepeatability),
                 MARGIN = ifelse(taxa_are_rows(physeqRepeatability), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(physeqRepeatability),
                      tax_table(physeqRepeatability))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
  # Plot the unique phyla: each dot will be a feature- total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeqRepeatability, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeqRepeatability),color=Phylum)) +
    # Include a threshold- here 5% = approx 18 samples
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2)+  geom_point(size = 0.3, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
  # Define prevalence threshold as 5% samples (abundance threshold already set as <50 reads)
  prevalenceThreshold<-5*(942/100)
  
  # Execute the prevalence filter, using `prune_taxa()` function
  head(prevdf1)
  
  KeepTaxa1prev<-prevdf1[prevdf1$Prevalence>=prevalenceThreshold,]
  str(KeepTaxa1prev)
  ASVNamesToKeep<- rownames(KeepTaxa1prev)
  
  physeqRepeatability
  physeqRepfiltered<- prune_taxa(ASVNamesToKeep, physeqRepeatability)
  physeqRepfiltered # 714 taxa, 942 samples

  
#Filter the phyloseq object to contain just the sequencing and repeat extractions
  
  TubeNosDups<- AllDups$TubeNumberUnique #tube nos of seq dups and repeat extractions
  
  physeqFilteredDups <- prune_samples(sample_data(physeqRepfiltered)$TubeNumberUnique %in% TubeNosDups, physeqRepfiltered)
  physeqFilteredDups #158 samples (69 sequencing duplicates + 10 extraction dup)
  
#CLR transformation of ASV abundances and sistance matrix
  physeq_clrDups <- microbiome::transform(physeqFilteredDups, "clr")
  
  #function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  #Extract OTU Matrix and Sample Data  
  clr_v<-vegan_otu(physeq_clrDups)
  clr_s<-as(sample_data(physeq_clrDups),"data.frame")  
  
#distance matrix
  distMatrixBeta <- vegdist(clr_v, method="euclidean")
  str(distMatrixBeta)
  distMatrixBeta <- as.matrix(distMatrixBeta)
  distMatrixBeta[1:3,1:3]
  
#extract as a dataframe (just upper triangle- one way comparisons)
  pairwiseDistBeta <- t(combn(colnames(distMatrixBeta), 2))
  head(pairwiseDistBeta)
  pairwiseDistBeta<- data.frame(pairwiseDistBeta, dist=distMatrixBeta[pairwiseDistBeta])
  head(pairwiseDistBeta)
  str(pairwiseDistBeta)
  colnames(pairwiseDistBeta)<- c("sampleID1", "sampleID2", "dist")
  
# get unique tube IDs
  TubesBeta<- AllDups[,c("sample.id","TubeNumberUnique")]
  head(TubesBeta)
  colnames(TubesBeta)<- c("sampleID1","ID1")
  pairwiseDistBeta<- merge(pairwiseDistBeta, TubesBeta, by="sampleID1", all.x=T)
  head(pairwiseDistBeta)
  
  colnames(TubesBeta)<- c("sampleID2","ID2")
  pairwiseDistBeta<- merge(pairwiseDistBeta, TubesBeta, by="sampleID2", all.x=T)
  head(pairwiseDistBeta)
  
  
# Add in sequence run number
  
  SeqRun1<- AllDups[,c("TubeNumberUnique","SequencingRun")]
  head(SeqRun1)
  colnames(SeqRun1)<-c("ID1", "Run1")
  pairwiseDistBeta<- merge(pairwiseDistBeta, SeqRun1, by="ID1", all.x=T)
  head(pairwiseDistBeta)
  
  SeqRun2<- AllDups[,c("TubeNumberUnique","SequencingRun")]
  head(SeqRun2)
  colnames(SeqRun2)<-c("ID2", "Run2")
  pairwiseDistBeta<- merge(pairwiseDistBeta, SeqRun2, by="ID2", all.x=T)
  head(pairwiseDistBeta)

  
#Copy the two Id columns so that we can check for tube no. matches in seq dup: ID1seq, ID2seq
  pairwiseDistBeta$ID1seq<- pairwiseDistBeta$ID1
  pairwiseDistBeta$ID2seq<- pairwiseDistBeta$ID2
  str(pairwiseDistBeta)
  
  #remove "b" from ID1seq and ID2seq to check for sequencing replicates (note, the second time sequenced has b after the tube number)
  pairwiseDistBeta$ID1seq<-gsub("b","",as.character(pairwiseDistBeta$ID1seq))
  pairwiseDistBeta$ID2seq<-gsub("b","",as.character(pairwiseDistBeta$ID2seq))
  str(pairwiseDistBeta)
  
  
#Copy the two Id columns so that we can check for tube no. matches in extraction dup: ID1seq, ID2seq
  pairwiseDistBeta$ID1extr<- pairwiseDistBeta$ID1
  pairwiseDistBeta$ID2extr<- pairwiseDistBeta$ID2
  str(pairwiseDistBeta)
  
  #remove "R" from ID1seq and ID2seq to check for sequencing replicates (note, the second time extracted has R after the tube number)
  pairwiseDistBeta$ID1extr<-gsub("R","",as.character(pairwiseDistBeta$ID1extr))
  pairwiseDistBeta$ID2extr<-gsub("R","",as.character(pairwiseDistBeta$ID2extr))
  str(pairwiseDistBeta)
  
#name type of comparison
  pairwiseDistBeta$Comparison<-ifelse(pairwiseDistBeta$ID1seq==pairwiseDistBeta$ID2seq, "SeqDuplicate",
                                  ifelse(pairwiseDistBeta$ID1extr==pairwiseDistBeta$ID2extr, "ExtractionDuplicate","DifferentSample"))
  str(pairwiseDistBeta)
  
  table(pairwiseDistBeta$Comparison) # 69 sequence duplicate comparisons, 10 extraction duplicate comparisons, 12324 between samples
  
  
#new column to say if duplicates of the same sample sequenced in the same or different run
  
  pairwiseDistBeta$Comparison2<-ifelse(pairwiseDistBeta$ID1seq==pairwiseDistBeta$ID2seq,
                                   ifelse(pairwiseDistBeta$Run1==pairwiseDistBeta$Run2,
                                          "DuplicateSameRun", "DuplicateDifferentRun"), pairwiseDistBeta$Comparison)
  str(pairwiseDistBeta)
  table(pairwiseDistBeta$Comparison2) # different run=46, same run = 23 (for 69 sequenced twice)
  
  
#plot to check between versus within samples
  
  samplePairs<- c("Different samples", "Extraction duplicates", "Sequencing duplicates")
  BetaPlot<- ggplot(pairwiseDistBeta, aes(Comparison, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=16)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=samplePairs) + theme(element_blank())
  BetaPlot
  
#plot to check sequencing duplicates across same and different run
  pairwiseDistBeta$Comparison2<- as.factor(pairwiseDistBeta$Comparison2)
  levels(pairwiseDistBeta$Comparison2)
  sampleCategory<- c("Different samples", "Duplicate different run", "Duplicate same run", "Extraction duplicate")
  BetaPlot2<- ggplot(pairwiseDistBeta, aes(Comparison2, dist)) +geom_jitter(col="darkgrey", width=0.2, shape=21, size=2, stroke=1)+ 
    geom_boxplot(alpha=0.1, outlier.shape=NA) + stat_summary(fun=mean, geom="point", shape=23, size=4, stroke=2) +
    xlab("\nPairwise comparison") +ylab("Euclidean distance\n") +theme_minimal() +
    theme(axis.title = element_text(size=20), axis.text = element_text(size=14)) +theme (axis.line = element_line(colour="black"))+
    scale_x_discrete(labels=sampleCategory) + theme(element_blank())
  BetaPlot2
  
# Kruskal-Wallis test
  
  kruskal.test(dist ~ Comparison2, data = pairwiseDistBeta)
  
  
#perform Dunn's Test with Bonferroni correction for p-values: 
  # no sig difference for duplicates sequenced in same or different run
  
  betaDunn<-dunnTest(dist ~ Comparison2,
           data=pairwiseDistBeta,
           method="bonferroni")
  
  dput(betaDunn$res, "DunnBetaRepeatability.txt") 
  