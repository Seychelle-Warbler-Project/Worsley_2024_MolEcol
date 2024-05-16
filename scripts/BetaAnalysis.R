setwd("/users/sarah/OneDrive - University of East Anglia/Seychelles_warbler/Analysis4_toS2022/MergedQIIME/Age_paper/Age/paper/FinalDraft/MolEcol/revisions/DRYAD/Beta")

#load packages
  
  library(phyloseq)
  library(ape)
  library(ggplot2)
  library(vegan)
  library(microbiome)
  library(ggordiplots)
  library(BiodiversityR)
  library(data.table)
  library(dplyr)
  library(plyr)
  #library(philr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
  library(car)
  library(broom.mixed)
  library(ggeffects)
  library(gllvm)
  library(gplots)
  library(readr)
  library(tidyr)



### load sample data ###

# decontam used to remove contaminants identified in extraction blanks and collection controls
# samples with <8000 reads removed
# ASVs with <50 reads in total across all samples removed


#### load filtered metadata with sample names as rownames (942 samples) ####
  metadataFiltered<-read.csv("FilteredMetadata.csv", row.names = 1)
  str(metadataFiltered)  
  head(metadataFiltered)
  

# standardise territory quality relative to mean in that season
  
  MeansTQ_FP<- aggregate(TQFinal_SameType ~ FieldPeriodID, data=metadataFiltered, FUN=mean) # mean of TQ
  colnames(MeansTQ_FP)<- c("FieldPeriodID", "MeanTQ_FP")
  metadataFiltered2<- merge(metadataFiltered, MeansTQ_FP, by="FieldPeriodID", all.x=T)
  str(metadataFiltered2)
  metadataFiltered2$MeanCentredTQ<- metadataFiltered2$TQFinal_SameType- metadataFiltered2$MeanTQ_FP
  
  
  ggplot(metadataFiltered2, aes(as.factor(FieldPeriodID), MeanTQ_FP)) + geom_point()
  ggplot(metadataFiltered2, aes(as.factor(FieldPeriodID), MeanCentredTQ)) + geom_point()


# make status categories
  
  Unknown<- c("SEEN1", "SEEN2")
  Sub<- c("AB", "ABX", "H")
  Dom <- c("BrF", "BrM")
  Juvenile <- c("CH", "FL", "OFL")
  metadataFiltered2$Status2<- ifelse(metadataFiltered2$Status %in% Unknown, "Unknown",
                                ifelse(metadataFiltered2$Status %in% Sub, "Subordinate", 
                                       ifelse(metadataFiltered2$Status %in% Dom, "Dominant", 
                                              ifelse(metadataFiltered2$Status %in% Juvenile, "Juvenile", metadataFiltered2$Status))))
  
  metadataFiltered2$Status2<- as.factor(metadataFiltered2$Status2)

  
# make ageclass category
  
  metadataFiltered2$Ageclass2<- factor(ifelse(metadataFiltered2$Ageclass=="CH", "N", 
                                         ifelse(metadataFiltered2$AgeYears <=0.25, "FL",
                                                ifelse(metadataFiltered2$AgeYears <= 0.5, "OFL", 
                                                       ifelse(metadataFiltered2$AgeYears <1, "SA",
                                                              ifelse(metadataFiltered2$AgeYears <3, "A 1-3 yrs",
                                                                     ifelse(metadataFiltered2$AgeYears <6, "A 3-6 yrs", "Old A")))))))
  
  
  
  table(metadataFiltered2$Ageclass2)
  ggplot(metadataFiltered2, aes(Ageclass2, AgeYears)) + geom_point()
  
  metadataFiltered2$Ageclass2<- factor(metadataFiltered2$Ageclass2, levels=c("N", "FL", "OFL", "SA", "A 1-3 yrs", "A 3-6 yrs", "Old A"))

 
# make category for if sample was taken in the terminal year or not
  
  metadataFiltered2$YrsToLastSeen<- metadataFiltered2$DaysToLastSeen/365.25
  str(metadataFiltered2)
  metadataFiltered2$TerminalYear<- ifelse(metadataFiltered2$YrsToLastSeen <=1 & metadataFiltered2$DiedEver =="yes", "yes", "no")
  table(metadataFiltered2$TerminalYear) # 218 samples taken in the terminal year, 724 not
  
  

  
  
##### make phyloseq object: pre-processed to remove contaminants/spurious ASVs  #####
  
#ASV table
  asvFiltered <- read.csv ("FilteredASVTable.csv", row.names=1) #read in asv table with feature names as rownames
  str (asvFiltered) #should be 942 samples, 23151
  head (asvFiltered)
  
  asvFiltered_table <- as.matrix (asvFiltered) #make into a matrix
  
#taxonomy
  taxonomyFiltered <- read.csv ("FilteredTaxonomyTable.csv", row.names=1)
  str (taxonomyFiltered) #should be 23151 observations, 7 taxonomic groupings, feature names as row names.
  taxonomyFiltered_matrix <- as.matrix (taxonomyFiltered)
  
#metadata
  str(metadataFiltered2)
  head(metadataFiltered2)
  rownames(metadataFiltered2)<- metadataFiltered2$sample.id

#read in tree as a phyloseq object
  phy_tree <- read_tree ("tree.nwk")

#import all as phyloseq objects
  ASV <- otu_table(asvFiltered_table, taxa_are_rows = TRUE)
  TAX <- tax_table(taxonomyFiltered_matrix)
  META <- sample_data(metadataFiltered2)

#merge into phyloseq object
  physeq <- phyloseq(ASV, TAX, META, phy_tree)
  physeq 

  sample_data(physeq)$SampleReads<- sample_sums(physeq)
  

  
###################
#### FILTERING ####
###################

#### Filter the phyloseq object to remove the sequencing and repeat extractions ####
  
  #extract those sequenced twice
  physeqMeta<- data.frame(sample_data(physeq))
  row_namesSeqDup <- as.vector(grep("b", rownames(physeqMeta), value=TRUE)) 
  row_namesExDup <- as.vector(grep("R", rownames(physeqMeta), value=TRUE)) 
  IDsDups<- c(row_namesExDup,row_namesSeqDup)
  Tubes<- physeqMeta[physeqMeta$sample.id %in% IDsDups,]
  TubesIDs<- Tubes$TubeNumber #tube numbers of duplicates (should match in both duplicates)
  
  AllDups<- physeqMeta[physeqMeta$TubeNumber %in% TubesIDs,] #Filter metadata to only include duplicated tube numbers 
  str(AllDups) #158
  AllDupsReads<- AllDups[,c("sample.id", "TubeNumber", "SampleReads")]
  str(AllDupsReads)
  
  #sort by tube number and reads (so the lowest read number is first)
  AllDupsReads<- AllDupsReads[
    with(AllDupsReads, order(TubeNumber,SampleReads)),
  ]
  
  #Keep only the lowest read no (first entry)- these will be removed from phyloseq object
  LowestReads<- AllDupsReads[!duplicated(AllDupsReads[,2]),]
  head(AllDupsReads)
  head(LowestReads)
  IDsLowest<- LowestReads$sample.id
  
  #Filter repeat with lowest reads from phyloseq obj
  physeq
  physeq2<- subset_samples(physeq, !sample.id %in% IDsLowest) 
  physeq2 #863 samples (repeats with the lowest read number removed)
  
  
#### remove floaters that have no territory quality ####
  
  physeq3<- subset_samples(physeq2, !Status == "FLOAT") #12 floaters removed
  
  physeq3 #851 samples remain
  
  
#### deal with catch repeats ####
  
  catchDupSorted<- read.csv("CatchDuplicatesSorted.csv")
  # duplicates from same catch sorted depending on source
  # (prioritise samples from tray, then bag, then other source, or highest read no if matched source)
  
  catchDupSorted<- catchDupSorted[,c("sample.id", "CatchDupKeep")]
  CatchToRemove<- subset(catchDupSorted, CatchDupKeep=="N")
  CatchToRemoveID<- CatchToRemove$sample.id
  
  physeq4<- subset_samples(physeq3, !sample.id %in% CatchToRemoveID)#793 samples
  physeq4
  length(unique(sample_data(physeq4)$CatchID)) # 789- 5 have NA (as not caught but sampled)- 789 includes the 1 NA + 4 = 793 unique catch ids
  length(unique(sample_data(physeq4)$BirdID)) # samples are from 439 birds, 50% of birds have >1 sample
  
  
  
################################################
### 1) Adults- age and terminal effects  #######
################################################
  
#### filter to adults only ####
  
  PhyseqA<- subset_samples(physeq4, AgeYears >= 1)
  PhyseqA
  str(sample_data(PhyseqA))
  
  #write_rds(PhyseqA, "PhyseqA.rds")
  PhyseqA<-readRDS("PhyseqA.rds")

  table(sample_data(PhyseqA)$Ageclass2) # 192 = Y adult, 172 = M adults, 98 = O adults

  
######## filter rare taxa from phyloseq object (those with low prevalence) #########
  
  # Compute prevalence of each feature (total number of samples in which a taxon appears at least once), store as data.frame
  prevdf = apply(X = otu_table(PhyseqA),
                 MARGIN = ifelse(taxa_are_rows(PhyseqA), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts for each phylum to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(PhyseqA),
                      tax_table(PhyseqA))
  head(prevdf)
  str(prevdf)
  
  plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))}) 
  #average prevelance of features within each phylum and the sum of feature prevalence within each phylum
  
  # Plot the unique phyla: each dot will be a feature
  #total abundance of that feature across samples on x axis and the prevalance (the fraction of all samples it occurs in on the y axis).
  prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(PhyseqA, "Phylum"))
  ggplot(prevdf1, aes(TotalAbundance, Prevalence / phyloseq::nsamples(PhyseqA),color=Phylum)) +
    # Include a threshold- here 5% 
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2)+  geom_point(size = 0.3, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) + theme(legend.position="none")
  
  # Define prevalence threshold as 5% samples (abundance threshold already set as <50 reads)
  PhyseqA
  prevalenceThresholdA<-5*(462/100) # ~23 samples
  
  # Execute the prevalence filter, using `prune_taxa()` function
  KeepTaxa1prev<-prevdf1[prevdf1$Prevalence>=prevalenceThresholdA,]
  str(KeepTaxa1prev)
  ASVNamesToKeep<- rownames(KeepTaxa1prev)
  
  PhyseqA
  sum(sample_sums(PhyseqA)) #33901412
  
  physeqBetaA<- prune_taxa(ASVNamesToKeep, PhyseqA)
  physeqBetaA  #462 samples, 674 taxa
  sum(sample_sums(physeqBetaA)) #26194392
  
  26194392/33901412 #77% of ASVs retained
  
  Prevalence<- (KeepTaxa1prev$Prevalence/462)*100
  range(Prevalence) # 5-84%
  mean(Prevalence) #16.14
  print(sd(Prevalence)/sqrt(length((Prevalence)))) #0.50
  
  Reads1<- sample_sums(PhyseqA)
  Reads2<- sample_sums(physeqBetaA)
  ReadsPerSample<- (Reads2/Reads1)*100
  mean(ReadsPerSample)
  print(sd(ReadsPerSample)/sqrt(length((ReadsPerSample)))) #1.12
  
  
  
############ CLR transform abundances ##########
  
  physeq_clr_A <- microbiome::transform(physeqBetaA, "clr")
  physeq_clr_A #462 samples, 674 taxa
  #saveRDS(physeq_clr_A, "physeq_clr_A.rds")
  
  ## Extract ASV Matrix and Sample Data  
  
  #function to extract an ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  clr_Matrix_A<-vegan_otu(physeq_clr_A)
  clr_SampleData_A<-as(sample_data(physeq_clr_A),"data.frame")
  str(clr_SampleData_A) #462 samples
  
  ## correct format of sample data
  clr_SampleData_A$SexEstimate <- as.factor(clr_SampleData_A$SexEstimate)
  clr_SampleData_A$Season <- as.factor(clr_SampleData_A$Season)
  clr_SampleData_A$FieldPeriodID <- as.factor(clr_SampleData_A$FieldPeriodID)
  clr_SampleData_A$BirdID <- as.factor(clr_SampleData_A$BirdID)
  clr_SampleData_A$Ageclass2 <- as.factor(clr_SampleData_A$Ageclass2)
  clr_SampleData_A$TerminalYear<- as.factor(clr_SampleData_A$TerminalYear)
  
  str(clr_SampleData_A)
  
  table(clr_SampleData_A$TerminalYear) #terminal samples: 356 no, 106 samples yes (within a year of death date)

  
################  permanova analysis ##############
  
## a) permanova no interaction terms 
  
  
  perm <- how(nperm = 9999)
  set.seed(3782)
  setBlocks(perm) <- with(clr_SampleData_A, BirdID)
  permanovaAdult<- adonis2(clr_Matrix_A ~ AgeYears + SexEstimate + MeanCentredTQ +
                             MeanTQ_FP + Season +
                             MinutesSinceSunrise +TimeAt4Degrees + TerminalYear,
                           data=clr_SampleData_A, permutations = perm, method = "euclidean", by= "margin", parallel = 2)
  permanovaAdult
  
  # terminal year isn't significant, age is marginally significant
  
  ResultsPermanovaAdult<- data.frame(permanovaAdult)
  write.csv(ResultsPermanovaAdult, "PermanovaAdult_terminal.csv")
  
  
  
## b) check with terminal year*Age and sex*Age interaction
  
  perm <- how(nperm = 9999)
  set.seed(3783)
  setBlocks(perm) <- with(clr_SampleData_A, BirdID)
  permanovaAdult_interaction<- adonis2(clr_Matrix_A ~ AgeYears + SexEstimate + MeanCentredTQ +
                                         MeanTQ_FP + Season +
                                         MinutesSinceSunrise +TimeAt4Degrees + TerminalYear +
                                         AgeYears*TerminalYear + SexEstimate*TerminalYear,
                                       data=clr_SampleData_A, permutations = perm, method = "euclidean", by= "margin", parallel = 2)
  permanovaAdult_interaction ### no significant interactions
  
  ResultsPermanovaAdult_interaction<- data.frame(permanovaAdult_interaction)
  write.csv(ResultsPermanovaAdult_interaction, "PermanovaAdult_interaction.csv")
  

## c) check results if use phyloseq object with no filtering 
  
  PhyseqA # 23151 ASVs, 462 samples
  
  physeqACLR_Rare<- microbiome::transform(PhyseqA, "clr")
  clr_Matrix_ARare<-vegan_otu(physeqACLR_Rare)
  clr_SampleData_ARare<-as(sample_data(physeqACLR_Rare),"data.frame")
  str(clr_SampleData_ARare) #462 samples
  
  ## correct format of sample data
  clr_SampleData_ARare$SexEstimate <- as.factor(clr_SampleData_ARare$SexEstimate)
  clr_SampleData_ARare$Season <- as.factor(clr_SampleData_ARare$Season)
  clr_SampleData_ARare$FieldPeriodID <- as.factor(clr_SampleData_ARare$FieldPeriodID)
  clr_SampleData_ARare$BirdID <- as.factor(clr_SampleData_ARare$BirdID)
  clr_SampleData_ARare$Ageclass2 <- as.factor(clr_SampleData_ARare$Ageclass2)
  clr_SampleData_ARare$TerminalYear<- as.factor(clr_SampleData_ARare$TerminalYear)
  
  str(clr_SampleData_ARare)
  
  
  perm <- how(nperm = 9999)
  set.seed(3785)
  setBlocks(perm) <- with(clr_SampleData_ARare, BirdID)
  permanovaAdult_rare<- adonis2(clr_Matrix_ARare ~ AgeYears + SexEstimate + MeanCentredTQ +
                             MeanTQ_FP + Season +
                             MinutesSinceSunrise +TimeAt4Degrees + TerminalYear,
                           data=clr_SampleData_ARare, permutations = perm, method = "euclidean", by= "margin", parallel = 2)
  permanovaAdult_rare 
  
  ResultsPermanovaAdult_NoFilter<- data.frame(permanovaAdult_rare)
  #write.csv(ResultsPermanovaAdult_NoFilter, "PermanovaAdult_NO_FILTER.csv") ### results the same
  
  
  
############ Principal Components Analysis (PCA) ###############
  
  
  clr_pca_A<-rda(clr_Matrix_A, scale=FALSE)
  
  screeplot(clr_pca_A) # PC1 gives most information, followed by 2.
  
  sig <- PCAsignificance(clr_pca_A, axes = 8) 
  sig # PC1 = 10.89, PC2 = 4.35, PC3 = 2.26, PC4 = 1.92
  
  
  #PC scores
  pca_scores<-scores(clr_pca_A, choices=c(1,2,3,4), scaling=1)
  sample_scores<- data.frame(pca_scores$sites)
  str(sample_scores)
  summary(sample_scores)
  sample_scores<-data.frame(setDT(sample_scores, keep.rownames = TRUE)[])
  str(sample_scores)
  colnames(sample_scores)<-c("sample.id", "PC1", "PC2","PC3", "PC4")
  
  
  
### Ageclass  ###
  
  str(clr_SampleData_A)
  clr_SampleData_A$Ageclass2<- factor(clr_SampleData_A$Ageclass2, levels=c("A 1-3 yrs", "A 3-6 yrs", "Old A"))
  
  plotPCA_age<-gg_ordiplot(clr_pca_A,groups=clr_SampleData_A$Ageclass2,ellipse = TRUE,plot=FALSE, pt.size=1, scaling=1, choices = c(3,4)) 
  plotPCA_age
  
  
  ord.data <- plotPCA_age$df_ord
  str(ord.data)
  head(ord.data)
  colnames(ord.data)<- c("x","y","Ageclass2")
  
  
  
  #centroids
  
  sample_scores$Ageclass2<- clr_SampleData_A$Ageclass2
  
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(Ageclass2), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(Ageclass2), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(Ageclass2), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(Ageclass2), summarize, mean=mean(PC4)))
  centroidsAge<- merge(centroidsPC1, centroidsPC2, by="Ageclass2")
  colnames(centroidsAge)<- c("Ageclass2", "x", "y")
  ord.data$AgeYears<- clr_SampleData_A$AgeYears
  
  ellipseAge<-data.frame(plotPCA_age$df_ellipse)
  colnames(ellipseAge)<-c("Ageclass2", "x", "y")
  ellipseAge$Ageclass2<-as.factor(ellipseAge$Ageclass2)
  levels(ellipseAge$Ageclass2)
  ellipseAge$Ageclass2<- factor(ellipseAge$Ageclass2,
                                levels=c("A 1-3 yrs","A 3-6 yrs","Old A"))
  
  levels(ord.data$Ageclass2)
  levels(centroidsAge$Ageclass2)
  levels(ellipseAge$Ageclass2)
  
  ord.data$BirdID<- as.factor(clr_SampleData_A$BirdID)
  
  PlotAge<- ggplot(data = ord.data, aes(x = x, y = y, fill = Ageclass2)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=21, alpha=0.7, colour="black")+
    geom_point(data=centroidsAge, aes(x,y), shape=23, stroke=2, size=5, colour="black")+
    scale_fill_manual(values=c("peachpuff","orangered3","royalblue2"),
                      labels=c("1-3 years", "3-6 years", "> 6 years"), name="Age (years)") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 10.89%") + ylab("PC2 4.35%") + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE)) 
  
  PlotAge 

  
  #+ geom_path(data = ellipseAge, 
  #aes(x = x, y = y, colour = Ageclass2))+
  #scale_colour_manual(values=c("red", "goldenrod1", "lightpink1", "darkorchid3","skyblue1", "royalblue2", "gray")) 
  

  # boxplot of scores
  
  ggplot(ord.data, aes(Ageclass2, x)) + geom_boxplot() + ylab("PC1 scores")
  ggplot(ord.data, aes(Ageclass2, y)) + geom_boxplot() + ylab("PC2 scores")
  
  
  
  #betadisper
  distMatrixBeta <- vegdist(clr_Matrix_A, method="euclidean")
  BDAge<-betadisper(distMatrixBeta, clr_SampleData_A$Ageclass2)
  set.seed(2699)
  permutest(BDAge, permutations = 9999) 
  boxplot(BDAge,ylab="Distance to centroid\n", xlab="\nAge class", cex.axis=1.5, cex.lab=1.5)
  
  
  
####### PCA of Season ######
  
  
  plotPCA_Season<-gg_ordiplot(clr_pca_A,groups=clr_SampleData_A$Season,ellipse = TRUE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCA_Season
  
  ord.dataSeason <- plotPCA_Season$df_ord
  str(ord.dataSeason)
  head(ord.dataSeason)
  colnames(ord.dataSeason)<- c("x","y","Season")
  
  #centroids
  
  sample_scores$Season<- clr_SampleData_A$Season
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(Season), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(Season), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(Season), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(Season), summarize, mean=mean(PC4)))
  centroidsSeason<- merge(centroidsPC1, centroidsPC2, by="Season")
  colnames(centroidsSeason)<- c("Season", "x", "y")
  
  
  PlotSeason<- ggplot(data = ord.dataSeason, aes(x = x, y = y, fill = Season)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=21, alpha=0.8, colour="grey55")+
    geom_point(data=centroidsSeason, aes(x,y), shape=23, stroke=2, size=5, colour="grey5")+
    scale_fill_manual(values=c("goldenrod", "darkgreen")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC1 10.89%") + ylab("PC2 4.35%") + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotSeason
  
  
  ### plot of scores
  
  ggplot(ord.dataSeason, aes(Season, x)) + geom_boxplot() + ylab("PC1 scores")
  ggplot(ord.dataSeason, aes(Season, y)) + geom_boxplot() + ylab("PC2 scores")
  
  
  #betadisper
  distMatrixBeta <- vegdist(clr_Matrix_A, method="euclidean")
  BDSeason<-betadisper(distMatrixBeta, clr_SampleData_A$Season)
  set.seed(2699)
  permutest(BDSeason, permutations = 9999) 
  boxplot(BDSeason,ylab="Distance to centroid\n", xlab="\nSeason", cex.axis=1.5, cex.lab=1.5)
  
  
  
##### PCA of minutes since sunrise #####
  
  clr_SampleData_A$Time_AmPm<- ifelse(clr_SampleData_A$MinutesSinceSunrise <360, "AM","PM")
  
  plotPCA_time<-gg_ordiplot(clr_pca_A,groups=clr_SampleData_A$Time_AmPm,ellipse = TRUE,plot=FALSE, pt.size=2, scaling=1, choices = c(3,4)) 
  plotPCA_time
  
  ord.data_time<- plotPCA_time$df_ord
  colnames(ord.data_time)<- c("x","y","Time")
  
  #centroids
  
  sample_scores$Time<- clr_SampleData_A$Time_AmPm
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC4)))
  centroidsTime<- merge(centroidsPC3, centroidsPC4, by="Time")
  colnames(centroidsTime)<- c("Time", "x", "y")
  centroidsTime$MinutesSinceSunrise<- c(220, 560)
  
  
  ord.data_time$Minutes<- clr_SampleData_A$MinutesSinceSunrise
  
  
  PlotTime<- ggplot(data = ord.data_time, aes(x = x, y = y, colour = Minutes)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=19, alpha=1) +
    geom_point(data=centroidsTime, aes(x,y, fill=Time), shape=23, stroke=2, size=5, colour="black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("\nPC3 2.26%") + ylab("PC4 1.93%") + 
    scale_colour_gradient(low = "darkorchid4", high = "gold") +
    scale_fill_manual(values=c("darkorchid4", "goldenrod")) +
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTime
  
  
  ### plot scores
  
  ggplot(ord.data_time, aes(Time, x)) + geom_boxplot() + ylab("PC3")
  ggplot(ord.data_time, aes(Time, y)) + geom_boxplot() + ylab("PC4")
  
  
  
  #betadisper
  BDtime<-betadisper(distMatrixBeta, clr_SampleData_A$Time_AmPm)
  set.seed(2699)
  permutest(BDtime, permutations = 9999) 
  boxplot(BDtime,ylab="Distance to centroid\n", xlab="\nTime", cex.axis=1.5, cex.lab=1.5)
  
  
  
##### PCA of TQ ######
  
  summary(clr_SampleData_A$MeanTQ_FP)
  
  clr_SampleData_A$MeanTQ_category<- ifelse(clr_SampleData_A$MeanTQ_FP < 17087, "low",
                                          ifelse(clr_SampleData_A$MeanTQ_FP <36602, "medium", "high"))
  
  plotPCA_TQ<-gg_ordiplot(clr_pca_A,groups=clr_SampleData_A$MeanTQ_category,ellipse = TRUE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCA_TQ
  
  ord.dataTQ <- plotPCA_TQ$df_ord
  colnames(ord.dataTQ)<- c("x","y","TQ")
  
  #centroids
  
  sample_scores$TQ<- clr_SampleData_A$MeanTQ_category
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(TQ), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(TQ), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(TQ), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(TQ), summarize, mean=mean(PC4)))
  centroidsTQ<- merge(centroidsPC1, centroidsPC2, by="TQ")
  colnames(centroidsTQ)<- c("TQ", "x", "y")
  
  ord.dataTQ$MeanTQ<- clr_SampleData_A$MeanTQ_FP
  
  centroidsTQ$MeanTQ<- c(17136,27929,36602)
  centroidsTQ$TerritoryQuality<- factor(centroidsTQ$TQ, levels=c("low", "medium", "high"))
  ord.dataTQ$TerritoryQuality<- factor(ord.dataTQ$TQ, levels=c("low", "medium", "high"))
  
  PlotTQ<- ggplot(data = ord.dataTQ, aes(x = x, y = y, fill = TerritoryQuality)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=21, alpha=0.8, col="grey55") +
    geom_point(data=centroidsTQ, aes(x,y, fill=TQ), shape=23, stroke=2, size=5, colour="black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_colour_gradient2(low="lightcyan1", mid="palegreen3", high="#043927", limits=c(9000, 55000), midpoint = 25000, n.breaks=20)+
    scale_fill_manual(values=c("purple4", "grey55", "goldenrod")) +
    xlab("\nPC1 10.89%") + ylab("PC2 4.35%") + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTQ
  
  #scores 
  
  ggplot(ord.dataTQ, aes(TerritoryQuality, x)) + geom_boxplot() + ylab("PC1 scores")
  ggplot(ord.dataTQ, aes(TerritoryQuality, y)) + geom_boxplot() + ylab("PC2 scores")
  
  
  
##### Time to freeze ######
  
  summary(clr_SampleData_A$TimeAt4Degrees)
  
  clr_SampleData_A$TimeAt4Degrees_cat<- ifelse(clr_SampleData_A$TimeAt4Degrees < 30, "< 30 days",
                                            ifelse(clr_SampleData_A$TimeAt4Degrees <= 60, "30-60 days", ">60 days"))
  
  clr_SampleData_A$TimeAt4Degrees_cat<- factor(clr_SampleData_A$TimeAt4Degrees_cat, levels=c("< 30 days", "30-60 days", ">60 days"))
  
  plotPCA_Freeze<-gg_ordiplot(clr_pca_A,groups=clr_SampleData_A$TimeAt4Degrees_cat,ellipse = TRUE,plot=FALSE, pt.size=2, scaling=1, choices = c(1,2)) 
  plotPCA_Freeze
  
  ord.dataFreeze <- plotPCA_Freeze$df_ord
  colnames(ord.dataFreeze)<- c("x","y","Time")
  
  #centroids
  
  sample_scores$Time<- clr_SampleData_A$TimeAt4Degrees_cat
  
  centroidsPC1<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC1)))
  centroidsPC2<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC2)))
  centroidsPC3<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC3)))
  centroidsPC4<-data.frame(ddply(sample_scores, .(Time), summarize, mean=mean(PC4)))
  centroidsTime<- merge(centroidsPC1, centroidsPC2, by="Time")
  colnames(centroidsTime)<- c("Time", "x", "y")
  
  ord.dataFreeze$Time<- clr_SampleData_A$TimeAt4Degrees_cat

  
  PlotTimeFreeze<- ggplot(data = ord.dataFreeze, aes(x = x, y = y, fill = Time)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=21, alpha=0.8, col="grey55") +
    geom_point(data=centroidsTime, aes(x,y, fill=Time), shape=23, stroke=2, size=5, colour="black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_fill_manual(values=c("navyblue", "darkgreen", "lightsteelblue3"),
                      labels= c("< 30 days", "30-60 days", "> 60 days"),
                      name= "Time at 4\u00B0C") +
    xlab("\nPC1 10.89%") + ylab("PC2 4.35%") + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTimeFreeze
  
  #scores 
  
  ggplot(ord.dataTQ, aes(TerritoryQuality, x)) + geom_boxplot() + ylab("PC1 scores")
  ggplot(ord.dataTQ, aes(TerritoryQuality, y)) + geom_boxplot() + ylab("PC2 scores")
  
  
#### continuous gradient for time at 4 degrees #####
  

  ord.dataFreeze$TimeCont<- as.numeric(clr_SampleData_A$TimeAt4Degrees)
  str(ord.dataFreeze)
  
  PlotTimeFreeze2<- ggplot(data = ord.dataFreeze, aes(x = x, y = y, colour = TimeCont)) + 
    theme_bw() + geom_point(size = 3, stroke=0.5, shape=19, alpha=1) +
    geom_point(data=centroidsTime, aes(x,y, fill=Time), shape=23, stroke=2, size=5, colour="black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_colour_gradientn(colours = c("#E8FF70","#CAF270","#417E4E","#042721"),
                           values = c(0,0.2,0.3,0.5,0.7, 1))+
    scale_fill_manual(values=c("#E8FF70","#94CD69", "#417E4E")) +
    xlab("\nPC1 10.89%") + ylab("PC2 4.35%") + 
    theme(axis.title = element_text(size=24), axis.text = element_text(size=22),
          legend.title = element_text(size=18), legend.text = element_text(size=18))+
    theme(legend.spacing.y = unit(0.3, 'cm'))  +
    guides(fill = guide_legend(byrow = TRUE))
  
  PlotTimeFreeze2
  
  


  
  
  

######################################################
### 2) GAMM/GLLVM ANALYSIS - individual taxa   #######
######################################################

  
################ a) gllvm #####################
  
  
  PhyseqA<- readRDS("PhyseqA.rds")
  PhyseqA # 462 samples
  
  CoreASV<- core(PhyseqA, detection = 0.01, prevalence = 30/100) # 96 ASVs are shared at 30% 
  
# changes in genus-level abundances with age
  genus_betaA<- tax_glom(PhyseqA, taxrank="Genus", NArm=FALSE) 
  genus_betaA #1829 genera
  #write_rds(genus_betaA, "genus_betaA.rds")
  genus_betaA<- readRDS("genus_betaA.rds")
  
  
  
#### filter to core taxa (i.e. >50% prevalence and 0.01 abundance minimum) ####

  physeqBeta50_genusA<- core(genus_betaA, detection = 0.01, prevalence = 50/100)
  physeqBeta50_genusA # 54 genera, 462 samples
  taxa(physeqBeta50_genusA)
  
  sum(sample_sums(genus_betaA)) 
  sum(sample_sums(physeqBeta50_genusA)) 
  
  (sum(sample_sums(physeqBeta50_genusA)))/(sum(sample_sums(genus_betaA))) #63% of reads retained
  
  
### average relative abundance of core taxa ###
  
  Core_RA<- microbiome::transform(physeqBeta50_genusA, "compositional")
  Core_RAmelt<- psmelt(Core_RA)
  Core_RAmelt$Abundance<- Core_RAmelt$Abundance*100 # change to percent
  head(Core_RAmelt)
  
  SError <- function(x) sd(x)/sqrt(length(x))
  
  AverageRA_Core<- Core_RAmelt %>%
    group_by(OTU) %>%
    summarise(Mean=mean(Abundance), SE=SError(Abundance))
  
  View(AverageRA_Core)
  
  mean(AverageRA_Core$Mean) #1.85%
  SError(AverageRA_Core$SE) #0.03%
  median(AverageRA_Core$Mean) #0.88%
  range(AverageRA_Core$Mean) # 0.2- 16.5% RA
  
  CoreTaxonomy<- unique(Core_RAmelt[,c("OTU","Order", "Family", "Genus")])
  AverageRA_Core2<- merge(AverageRA_Core, CoreTaxonomy, by="OTU")
  View(AverageRA_Core2)
  
#### CLR transform ####
  
  genus_clr_scaled <- microbiome::transform(physeqBeta50_genusA, "clr")
  
#############
  
  taxa_names(genus_clr_scaled)
  
##Data preparation
  
  
  ys <- data.frame(t(otu_table(genus_clr_scaled)))
  
  rowSums(ys)
  
  names(ys) <-taxa_names(genus_clr_scaled)
  
  head(ys)[1:5,1:5]
  
  Xs<-data.frame(sample_data(genus_clr_scaled))
  head(Xs)
  
  
  
#### correct format of sample data ####
  
  Xs$SexEstimate <- as.factor(Xs$SexEstimate)
  Xs$Season <- as.factor(Xs$Season)
  Xs$FieldPeriodID <- as.factor(Xs$FieldPeriodID)
  Xs$BirdID <- as.factor(Xs$BirdID)
  Xs$Ageclass2 <- as.factor(Xs$Ageclass2)
  Xs$TerminalYear<- as.factor(Xs$TerminalYear)
  
  
  
### centre and scale continuous variables  ###
  
  Xs$AgeYears.sc<- arm::rescale(Xs$AgeYears, binary.inputs = "full")
  Xs$MeanTQ_FP.sc<- arm::rescale(Xs$MeanTQ_FP, binary.inputs = "full")
  Xs$CentredTQ.sc<- arm::rescale(Xs$MeanCentredTQ, binary.inputs = "full")
  Xs$TimeOfDay.sc<- arm::rescale(Xs$MinutesSinceSunrise, binary.inputs = "full")
  Xs$TimeToFreeze.sc <- arm::rescale(Xs$TimeAt4Degrees, binary.inputs = "full")
  
  
  names(Xs)
  
  Xs<-Xs[,c("Season", "SexEstimate", "TerminalYear","AgeYears.sc",
            "TimeOfDay.sc","TimeToFreeze.sc",
            "CentredTQ.sc", "MeanTQ_FP.sc", "BirdID")]
  names(Xs)
  
  
  
  
  
  
### Model fitting with gllvm ###
  
  
  fit_model <- gllvm(ys, Xs, 
                     num.lv = 2,
                     starting.val='zero',
                     formula = ~  AgeYears.sc +
                       SexEstimate +
                       Season + TerminalYear + MeanTQ_FP.sc+
                       CentredTQ.sc +
                       TimeToFreeze.sc+
                       TimeOfDay.sc,
                     family = "gaussian",
                     row.eff = ~(1|BirdID))
  
  
  
  coefplot(fit_model, cex.ylab = 0.5 ,mar = c(5,12,2,1), which.Xcoef = c(1))
  
  #par(mfrow = c(3, 2), mar = c(4, 4, 2, 1))
  #plot(fit_model, var.colors = 1)
  
  summary(fit_model)$AICc # for 1 lv= 118185, 2 lv = 117325.3, 
  
  
######## Extract estimates #########
  
  df<-coef(fit_model)
  est_df<-data.frame(df$Intercept)
  est_df2<-data.frame(df$Xcoef) 
  est_df3<-merge(est_df, est_df2, by = 0)
  
  est_df3
  
  # order genera
  
  row.names(est_df3)<-est_df3$Row.names
  est_df3<-est_df3[colnames(ys),]
  
  #put est_df3 into long format
  
  names(est_df3)[1]<- "ASV"
  names(est_df3)[2]<- "Intercept"
  
  estimates_df <- gather(est_df3, Treatment, Estimate, names(est_df3)[2]:names(est_df3)[ncol(est_df3)], factor_key=TRUE)
  head(estimates_df)
  
  
###### extract confindence intervals  #######
  
  
  confint_df<-data.frame(confint(fit_model))
  
  
  confint_df<-rbind(confint_df[rownames(confint_df) %like% "Xcoef", ],
                    confint_df[rownames(confint_df) %like% "Intercept", ])
  
  head(confint_df)
  
# add a column with correct variable level
  variables<- colnames(est_df3)[3:ncol(est_df3)]
  variables<-c(variables, "Intercept")
  variables1<-rep(variables, nrow(est_df))
  variables2<-variables1[order(match(variables1, variables))]
  
  confint_df$Treatment<-variables2
  
  
# column with taxa names. Should be automatically in the correct order but double check
  
  confint_df$ASV<-rep(colnames(ys), length(unique(confint_df$Treatment)))
  
  
  merged<-merge(estimates_df, confint_df, by = c("Treatment", "ASV"))
  
  
  final_estimates<- merged
  names(final_estimates)[4]<-"CI_lower"
  names(final_estimates)[5]<-"CI_upper"
  
  
  head(final_estimates)
  unique(final_estimates$Treatment)
  
  
## add significance
  
  final_estimates$Significant<- !data.table::between(0, final_estimates$CI_lower, final_estimates$CI_upper)
  
  
  final_estimates$ASV<-factor(final_estimates$ASV)
  head(final_estimates)
  
  
####  rename ambiguous genera ####
  GenusNames1<-data.frame(tax_table(genus_clr_scaled))
  head(GenusNames1)
  GenusNames1$ASV<- rownames(GenusNames1)
  GenusNames2<-GenusNames1[,c("ASV","Order", "Family", "Genus", "Phylum")]
  head(GenusNames2)
  
  length(unique(GenusNames2$ASV)) #54 unique genera
  length(unique(GenusNames2$Genus)) #only 38 unique names - the unknown genera are not given a name, also need to rename uncultured
  
  
  GenusNames2$Family<- ifelse(GenusNames2$Family== "uncultured bacterium", "", GenusNames2$Family)
  
  GenusNames2$Family<- ifelse( GenusNames2$Family=="", NA, GenusNames2$Family)
  
  GenusNames3<- GenusNames2 %>% 
    mutate(Family = coalesce(Family,Order)) # give name of order if no family
  
  GenusNames3$Family
  
  
  GenusNames3$Number<- 1:54
  GenusNames3$NewName<- paste(GenusNames3$Family,"_",GenusNames3$Number)
  head(GenusNames3)
  
  GenusNames3$Genus<- ifelse(GenusNames3$Genus=="", NA, GenusNames3$Genus)
  
  GenusNames4<- GenusNames3 %>% 
    mutate(Genus = coalesce(Genus,NewName)) # give family name if no genus
  
  head(GenusNames4)
  #View(GenusNames4)
  
  GenusNames4$Genus<- ifelse(GenusNames4$Genus== "uncultured bacterium", GenusNames4$NewName, GenusNames4$Genus)
  GenusNames4$Genus<- ifelse(GenusNames4$Genus== "uncultured", GenusNames4$NewName, GenusNames4$Genus)
  GenusNames4$Genus<- ifelse(GenusNames4$Genus== "uncultured soil bacterium", GenusNames4$NewName, GenusNames4$Genus)
  
  length(unique(GenusNames4$Genus))#54 unique names
  
  
  final_estimates2<- merge(final_estimates, GenusNames4, by="ASV") # add in taxa names
  head(final_estimates2)
  levels(final_estimates2$Treatment)

  
  final_estimatesAge<- subset(final_estimates2, Treatment=="AgeYears.sc")
  head(final_estimatesAge)
  NewAge<-final_estimatesAge[order(final_estimatesAge$Estimate),] # order taxa according to estimate for age
  NewAge$OrderGenus<- 1:54
  head(NewAge)
  NewAgeMatch<- NewAge[,c("ASV", "OrderGenus")]
  final_estimates3<- merge(final_estimates2,NewAgeMatch)
  
  str(final_estimates3$OrderGenus)
  subset(final_estimates3,OrderGenus==1) #check
  
  final_estimates3$Genus<- as.factor(final_estimates3$Genus)
  final_estimates3$Genus<- factor(final_estimates3$Genus, levels=unique(final_estimates3$Genus[order(final_estimates3$OrderGenus,final_estimates3$Genus)]), ordered=TRUE)
  
  #rename for graph
  namesGen<- unique(final_estimates3[,c("ASV", "Genus")])
  originalNames<-levels(final_estimates3$Genus)
  originalNames
  NewNames<-  c("Ruminococcaceae genus 1","Gordonia","Lactococcus", "Parabacteroides", 
                "Desulfovibrio","Lachnoclostridium","Candidatus Soleaferrea", "Ruminococcaceae genus 2",
                "Saccharimonadaceae genus 1", "Lachnospiraceae genus 1", "Christensenellaceae genus 1", 
                "Akkermansia", "Rhodococcus", "Mycobacterium", "Enterobacteriaceae genus 1",
                "Corynebacterium genus 1" ,"Nocardioides", "Rhodospirillaceae genus 1", 
                "Isosphaeraceae genus 1", "Rhodobacteraceae genus 1", "Saccharimonadales genus 1",
                "Thermomicrobiales JG30-KF-CM45 genus 1", "Friedmanniella" ,"Thermomicrobiales JG30-KF-CM45 genus 2",
                "Bacteroides","Amnibacterium" , "Actinomycetospora","Paracoccus", "Enterococcus" ,
                "Marmoricola", "Rhizobiaceae genus 1", "Stenotrophomonas" ,"Saccharimonadales genus 2",
                "Curtobacterium", "Microbacteriaceae genus 1", "Acetobacteraceae genus 1","Beijerinckiaceae genus 1", 
                 "Devosia" , "Aureimonas","Methylobacterium" ,  "Rhizobium","Gemmata"  ,"Singulisphaera"  , 
                "Thermomicrobiales JG30-KF-CM45 genus 3" ,"Solirubrobacterales 67-14 genus 1", "Kocuria", 
                "Geodermatophilus", "Microbacterium",  "Sphingomonas","Quadrisphaera",  "Pseudonocardia",
                "Kineococcus","Pantoea" ,"Micromonosporaceae genus 1" )
  Namedf<- data.frame(originalNames, NewNames)
  colnames(Namedf)<- c("Genus", "NewNames")
  head(Namedf)
  namesGen2<- merge(namesGen, Namedf, by="Genus")
  head(namesGen2)
  namesGen3<- namesGen2[,c(2,3)]
  
  final_estimates4<- merge(final_estimates3, namesGen3, by="ASV")
  head(final_estimates4)
  levels(final_estimates4$Treatment)
  
  
  final_estimates4$NewNames = factor(final_estimates4$NewNames, levels=unique(final_estimates4$NewNames[order(final_estimates4$OrderGenus)]), ordered=TRUE)
  str(final_estimates4)
  
  levels(final_estimates4$Treatment)
  final_estimates4$Treatment <- factor(final_estimates4$Treatment,
                                       levels = c("Intercept", "AgeYears.sc", "SexEstimate1", "SeasonMinor","TerminalYearyes",
                                                  "MeanTQ_FP.sc","CentredTQ.sc", "TimeToFreeze.sc", "TimeOfDay.sc"), 
                                       labels = c("Intercept","Age", "Sex (male)", "Season (minor)","Terminal year (yes)",
                                                  "Mean TQ","Mean-centered TQ","Time at 4\u00B0C", "Time of day"))
  
  
#### make plots #####
  
# plot for age and terminal yr
  plot1<- c("Age", "Terminal year (yes)")

  final_estimatesplot<- subset(final_estimates4, Treatment %in% plot1)
  head(final_estimatesplot)
  
  
  p1<-ggplot(final_estimatesplot, aes(x=Estimate, y=NewNames, fill=Phylum, col=Significant)) +
    geom_linerange(aes(xmin=CI_lower, xmax=CI_upper), linewidth=1)+
    geom_point(shape=21)+
    scale_color_manual(values=c("grey70", "grey3"))+
    geom_vline(xintercept = 0, linetype="dashed")+
    facet_grid(. ~ Treatment)+
    ylab("")
  p1+ theme_bw()+theme(panel.grid.minor = element_blank())+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=20),
          legend.text = element_text(size=16), legend.title = element_text(size=20))
  
  
# plot for other predictors
  final_estimatesplot2<- subset(final_estimates4, Treatment!="Age")
  final_estimatesplot2<- subset(final_estimatesplot2, Treatment!="Terminal year (yes)")
  final_estimatesplot2<- subset(final_estimatesplot2, Treatment!="Intercept")
  
  head(final_estimatesplot2)
  
  p2<-ggplot(final_estimatesplot2, aes(x=Estimate, y=NewNames, fill=Phylum, col=Significant)) +
    geom_linerange(aes(xmin=CI_lower, xmax=CI_upper), linewidth=1)+
    geom_point(shape=21)+
    scale_color_manual(values=c("grey70", "grey3"))+
    geom_vline(xintercept = 0, linetype="dashed") + 
    facet_grid(. ~ Treatment)+
    ylab("")
  p2+ theme_bw()+theme(panel.grid.minor = element_blank())+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=20),
          legend.text = element_text(size=16), legend.title = element_text(size=20),
          strip.text = element_text(size=12)) +
    theme(legend.position="bottom", legend.box = "vertical")
  
  
  
#Find mean abundance of significant taxa
  
  head(final_estimatesplot)
  
  SigAgeCoreGen<- subset(final_estimatesplot, Treatment=="Age")
  SigAgeCoreGen2<- subset(SigAgeCoreGen, Significant=="TRUE")
  
  SigIDs<- SigAgeCoreGen2$ASV
  
  
  genus_betaA
  
  genus_betaAcomp<-transform(genus_betaA, "compositional") # converts to relative abundance between 0 and 1
  sample_sums(genus_betaAcomp)
  
  genus_betaAlong<-psmelt(genus_betaAcomp)
  head(genus_betaAlong)
  
  genus_betaAlong$Abundance2<- genus_betaAlong$Abundance*100
  BetaMean<- genus_betaAlong %>%
    group_by(OTU, Genus, Family) %>%
    dplyr::summarize(Max = max(Abundance2), Mean = mean(Abundance2, na.rm=TRUE), SD= sd(Abundance),
                     SE = print(sd(Abundance)/sqrt(length((Abundance)))))
  
  head(BetaMean)
  AbundancesOrdered<- BetaMean[order(BetaMean$Mean, decreasing = TRUE),]
  print(AbundancesOrdered, n=20)
  
  
  SigTaxa<- subset(genus_betaAlong, OTU %in% SigIDs)
  head(SigTaxa) 
  SigTaxa$Abundance<- SigTaxa$Abundance*100
  
  SigTaxa %>%
    group_by(OTU, Genus, Family) %>%
    dplyr::summarize(Mean = mean(Abundance, na.rm=TRUE), SD= sd(Abundance),
                     SE = print(sd(Abundance)/sqrt(length((Abundance)))))
  
  
  str(SigTaxa)
  SigTaxa %>%
    group_by(OTU, Ageclass2, Genus, Family) %>%
    dplyr::summarize(Mean = mean(Abundance, na.rm=TRUE), SD= sd(Abundance))
  
  
  
  
  
################### b) GAMM ANALYSIS - individual taxa in adults #################
  
  
  genus_betaA<- readRDS("genus_betaA.RDS")
  
  
#### filter to core taxa (i.e. >50% prevalence and 0.01 abundance minimum) ####
  
  physeqBeta50_genusA<- core(genus_betaA, detection = 0.01, prevalence = 50/100)
  physeqBeta50_genusA # 54 genera, 462 samples
  taxa(physeqBeta50_genusA)
  
  sum(sample_sums(genus_betaA)) 
  sum(sample_sums(physeqBeta50_genusA)) 
  
  (sum(sample_sums(physeqBeta50_genusA)))/(sum(sample_sums(genus_betaA))) #63% of reads retained
  

#### CLR transform ####
  
  Genus_clr50A <- microbiome::transform(physeqBeta50_genusA, "clr")
  
  
#### convert to long and rename ambiguous genera ####
  
  Genus_clr50longA<-psmelt(Genus_clr50A)
  head(Genus_clr50longA)
  
  length(unique(Genus_clr50longA$OTU)) #54 unique genera
  length(unique(Genus_clr50longA$Genus)) #only 38 unique names - the unknown genera are not given a name, also need to rename uncultured
  
  taxaNamesA<-Genus_clr50longA[,c("OTU", "Order", "Family", "Genus")]
  unique(taxaNamesA)
  length(unique(taxaNamesA$OTU))
  
  
  Genus_clr50longA2<- Genus_clr50longA
  
  Genus_clr50longA2$Family<- ifelse(Genus_clr50longA2$Family== "uncultured bacterium", "", Genus_clr50longA2$Family)
  
  Genus_clr50longA2$Family<- ifelse( Genus_clr50longA2$Family=="", NA, Genus_clr50longA2$Family)
  
  Genus_clr50longA2<- Genus_clr50longA2 %>% 
    mutate(Family = coalesce(Family,Order))
  
  taxaNamesA2<-Genus_clr50longA2[,c("OTU", "Order", "Family", "Genus")]
  taxaNamesUniqueA<-unique(taxaNamesA2)
  
  taxaNamesUniqueA$Number<- 1:54
  taxaNamesUniqueA$NewName<- paste(taxaNamesUniqueA$Family,"_",taxaNamesUniqueA$Number)
  head(taxaNamesUniqueA)
  taxaNamesReplaceA<- taxaNamesUniqueA[,c("OTU", "NewName")]
  Genus_clr50longA3<- merge(Genus_clr50longA2, taxaNamesReplaceA, by="OTU")
  head(Genus_clr50longA3)
  
  Genus_clr50longA3$Genus<- ifelse( Genus_clr50longA3$Genus=="", NA, Genus_clr50longA3$Genus)
  
  Genus_clr50longA3<- Genus_clr50longA3 %>% 
    mutate(Genus = coalesce(Genus,NewName))
  
  head(Genus_clr50longA3)
  
  Genus_clr50longA3$Genus<- ifelse(Genus_clr50longA3$Genus== "uncultured bacterium", Genus_clr50longA3$NewName, Genus_clr50longA3$Genus)
  Genus_clr50longA3$Genus<- ifelse(Genus_clr50longA3$Genus== "uncultured", Genus_clr50longA3$NewName, Genus_clr50longA3$Genus)
  Genus_clr50longA3$Genus<- ifelse(Genus_clr50longA3$Genus== "uncultured soil bacterium", Genus_clr50longA3$NewName, Genus_clr50longA3$Genus)
  
  taxaNamesA<-Genus_clr50longA3[,c("OTU", "Order", "Family", "Genus")]
  unique(taxaNamesA)
  length(unique(taxaNamesA$Genus))
  
  
  
#### correct format of sample data ####
  
  Genus_clr50longA3$SexEstimate <- as.factor(Genus_clr50longA3$SexEstimate)
  Genus_clr50longA3$Season <- as.factor(Genus_clr50longA3$Season)
  Genus_clr50longA3$FieldPeriodID <- as.factor(Genus_clr50longA3$FieldPeriodID)
  Genus_clr50longA3$BirdID <- as.factor(Genus_clr50longA3$BirdID)
  Genus_clr50longA3$Ageclass2 <- as.factor(Genus_clr50longA3$Ageclass2)
  Genus_clr50longA3$TerminalYear<- as.factor(Genus_clr50longA3$TerminalYear)
  
  
### centre and scale continuous variables  ###
  
  Genus_clr50longA3$AgeYears.sc<- arm::rescale(Genus_clr50longA3$AgeYears, binary.inputs = "full")
  Genus_clr50longA3$MeanTQ_FP.sc<- arm::rescale(Genus_clr50longA3$MeanTQ_FP, binary.inputs = "full")
  Genus_clr50longA3$CentredTQ.sc<- arm::rescale(Genus_clr50longA3$MeanCentredTQ, binary.inputs = "full")
  Genus_clr50longA3$TimeOfDay.sc<- arm::rescale(Genus_clr50longA3$MinutesSinceSunrise, binary.inputs = "full")
  Genus_clr50longA3$TimeToFreeze.sc <- arm::rescale(Genus_clr50longA3$TimeAt4Degrees, binary.inputs = "full")
  
  str(Genus_clr50longA3)
  

  
##### loop to fit GAMMs to each genus - prints results to PDF #####
  range(Genus_clr50longA3$AgeYears.sc)
  mean(Genus_clr50longA3$AgeYears) #4.234792
  2*sd(Genus_clr50longA3$AgeYears) #6.173092
  
  mean(Genus_clr50longA3$MinutesSinceSunrise) #367.4524
  2*sd(Genus_clr50longA3$MinutesSinceSunrise) #340.2907
  
  mean(Genus_clr50longA3$MeanTQ_FP) #26994.68
  2*sd(Genus_clr50longA3$MeanTQ_FP) # 25712.83
  
  resultsA<- list()
  age_estimatesA<-list()
  Time_estimatesA<- list()
  TerminalYearA<- list()
  TQA<- list()
  
  taxanamesA<-unique(Genus_clr50longA3$Genus)
  
  
  pdf("taxaPlotsA.pdf")
  
  for (i in 1:length(taxanamesA)){
    tryCatch({ #catch errors
      print(i)
      print(taxanamesA[i])
      
      taxa1A<-subset(Genus_clr50longA3, Genus == taxanamesA[i])
      metadataA<-taxa1A

      # fit gamm    
      m_taxaA <- mgcv::bam(Abundance~
                             s(AgeYears.sc, bs="cr", k=6) + 
                             s(TimeOfDay.sc, bs="cr") + 
                             s(TimeToFreeze.sc, bs="cr") +
                             SexEstimate +
                             Season +
                             TerminalYear +
                             CentredTQ.sc+
                             MeanTQ_FP.sc +
                             s(BirdID, bs="re"),
                           data=metadataA,
                           family = gaussian)
      
      print(summary(m_taxaA))
      
      predA <- data.frame(ggpredict(m_taxaA, terms = c("AgeYears.sc[all]")))  # this gives overall predictions for the model
      predA$Age<- (predA$x*6.173092)+4.234792
      
      PredPlotA<- ggplot(predA) +
        geom_point(data=metadataA, aes(x=AgeYears, y=Abundance)) +
        theme_bw() +
        geom_ribbon(aes(x = Age, ymin = conf.low, ymax = conf.high), 
                    fill = "lightgrey", alpha = 0.5) +
        geom_line(aes(x = Age, y = predicted), linewidth=2) +          
        labs(x = "\nAge", y = "Abundance\n") + ggtitle(taxanamesA[i]) +
        scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18))+
        theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))+
        theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(PredPlotA)
      
      
      
      predA2 <- data.frame(ggpredict(m_taxaA, terms = c("TimeOfDay.sc[all]")))  # this gives overall predictions for the model
      predA2$Time<- (predA2$x*340.2907)+367.4524
      
      PredPlotA2<- ggplot(predA2) +
        geom_point(data=metadataA, aes(x=MinutesSinceSunrise, y=Abundance)) +
        geom_ribbon(aes(x = Time, ymin = conf.low, ymax = conf.high), 
                    fill = "lightgrey", alpha = 0.5) +
        geom_line(aes(x = Time, y = predicted), linewidth=2) +          
        labs(x = "\nMinutes since sunrise", y = "Abundance\n") + ggtitle(taxanamesA[i]) +
        theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) +
        theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(PredPlotA2)
      
      
      predA3 <- data.frame(ggpredict(m_taxaA, terms = c("TerminalYear")))  # this gives overall predictions for the model
      
      PredPlotA3<- ggplot(predA3) +
        geom_boxplot(data=metadataA, aes(x=TerminalYear,y=Abundance))+
        geom_point(data=metadataA, aes(x=TerminalYear, y=Abundance),
                   position= position_jitter(width=0.1), alpha=0.3) +
        geom_linerange(aes(x=x, ymin=conf.low, ymax=conf.high), linewidth=2)+
        geom_point(aes(x = x, y = predicted), col="red", size=4) + 
        theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))+
        labs(x = "\nTerminal year", y = "Abundance\n") + ggtitle(taxanamesA[i]) +
        theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(PredPlotA3)
      
      predA4 <- data.frame(ggpredict(m_taxaA, terms = c("MeanTQ_FP.sc[all]")))  # this gives overall predictions for the model
      predA4$TQ<- (predA4$x*25712.83)+26994.68
      
      PredPlotA4<- ggplot(predA4) +
        geom_point(data=metadataA, aes(x=MeanTQ_FP, y=Abundance)) +
        geom_ribbon(aes(x = TQ, ymin = conf.low, ymax = conf.high), 
                    fill = "lightgrey", alpha = 0.5) +
        geom_line(aes(x = TQ, y = predicted), linewidth=2) +          
        labs(x = "\nMean territory quality", y = "Abundance\n") + ggtitle(taxanamesA[i]) +
        theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) +
        theme(plot.margin = unit(c(0.1,0.4,0.1,0.1), "cm"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      print(PredPlotA4)
      
      resultsA[[i]]<- summary(m_taxaA)
      #gam.check(m_taxa)
      
      ## add effect size and p value

      summary<-summary(m_taxaA)
      textplot(capture.output(summary(m_taxaA)))

      resultsdf1<- data.frame(degFreedom= summary$s.table[1,1], effectSize= summary$s.table[1,3], P_val=summary$s.table[1,4])
      resultsdf2<- data.frame(degFreedom= summary$s.table[2,1], effectSize= summary$s.table[2,3], P_val=summary$s.table[2,4])
      resultsdf3<- data.frame(estimate= summary$p.table[4,1], SE= summary$p.table[4,2], P_val=summary$p.table[4,4])
      resultsdf4<- data.frame(estimate= summary$p.table[6,1], SE= summary$p.table[6,2], P_val=summary$p.table[6,4])
      
      age_estimatesA[[i]]<-resultsdf1
      Time_estimatesA[[i]]<- resultsdf2
      TerminalYearA[[i]]<- resultsdf3
      TQA[[i]]<- resultsdf4
      
    }, error=function(e){})
  }
  
  dev.off()
  
  names(age_estimatesA)<-taxanamesA
  names(resultsA)<-taxanamesA
  names(Time_estimatesA)<-taxanamesA
  names(TerminalYearA)<-taxanamesA
  names(TQA)<-taxanamesA
  
  age_estimates_dfA<-data.frame(do.call(rbind, age_estimatesA))
  age_estimates_dfA
  sigAge_A<- subset(age_estimates_dfA, P_val<0.05)
  sigAge_A
  str(sigAge_A) #7/167 are significant with age= 8%



  
  
  
  
  