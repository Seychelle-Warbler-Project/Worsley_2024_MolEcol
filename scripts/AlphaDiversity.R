setwd("/users/sarah/OneDrive - University of East Anglia/Seychelles_warbler/Analysis4_toS2022/MergedQIIME/Age_paper/Age/paper/FinalDraft/MolEcol/revisions/DRYAD/Alpha")

#load packages
  library(phyloseq)
  library(ape)
  library(ggplot2)
  library(vegan)
  library(picante)
  library(car)
  library(DHARMa)
  library(dplyr)
  library(ggeffects)
  library(arm)
  library(mgcv)
  library(mgcv.helper)
  library(multcomp)
  library(visreg)
  library(ggpubr)
  library(patchwork)

### make phyloseq object: tables pre-processed to remove contaminants/spurious ASVs ###

# decontam used to remove contaminants identified in extraction blanks and collection controls
# samples with <8000 reads removed
# ASVs with <50 reads in total across all samples removed


#### load filtered metadata with sample names as rownames (942 samples) ####
  metadataFiltered<-read.csv("FilteredMetadata.csv", row.names = 1)
  str(metadataFiltered)  
  head(metadataFiltered)


# standardise territory quality relative to mean in that season
  
  MeansTQ_FP<- aggregate(TQFinal_SameType ~ FieldPeriodID, data=metadataFiltered, FUN=mean)
  colnames(MeansTQ_FP)<- c("FieldPeriodID", "MeanTQ_FP")

  metadataFiltered2<- merge(metadataFiltered, MeansTQ_FP, by="FieldPeriodID", all.x=T)
  str(metadataFiltered2)
  
  metadataFiltered2$MeanCentredTQ<- metadataFiltered2$TQFinal_SameType - metadataFiltered2$MeanTQ_FP
  
  ggplot(metadataFiltered2, aes(as.factor(FieldPeriodID), MeanTQ_FP)) + geom_point()
  ggplot(metadataFiltered2, aes(as.factor(FieldPeriodID), MeanCentredTQ)) + geom_point()
  ggplot(metadataFiltered2, aes(as.factor(FieldPeriodID), TQFinal_SameType)) + geom_point()
  


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

#load filtered metadata with sample names as rownames (942 samples)
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

  
#### filter to adults only ####

  PhyseqAdults<- subset_samples(physeq4, AgeYears >= 1)
  PhyseqAdults # 462 samples, 23151
  table(sample_data(PhyseqAdults)$Ageclass2) # 192 = Y adult, 172 = M adults, 98 = O adults

  PhyseqA <- prune_taxa(taxa_sums(PhyseqAdults) > 0, PhyseqAdults) # get rid of ASVs that have 0 reads after all the filtering steps above
  PhyseqA #17159 taxa , 462 samples
  
  

####################  
### RAREFY READS ###
####################

#rarefy to 8000 and set seed before rarefying (58392), so that results are reproducible  
  physeqRare<-rarefy_even_depth(PhyseqA, 8000, rngseed = 58392)

# 581 ASVs removed after subsampling- leaves 16578 ASVs across 462 samples
  physeqRare
  sample_sums(physeqRare)
  
  observationThreshold <- 1
  sumASVs<- (apply(otu_table(physeqRare) > observationThreshold, 2, sum))
  str(sumASVs)
  sumASVs <- data.frame(sumASVs)
  str(sumASVs)
  mean(sumASVs[,1]) #mean ASVs per sample = 161.58 in rarefied dataset
  sd(sumASVs[,1]) #Standard deviation of ASVs per sample = 110.46
  print(sd(sumASVs[,1])/sqrt(length((sumASVs[,1])))) #SE = 5.139
  range(sumASVs[,1]) # range 3-731



############################################
##### CALCULATE ALPHA DIVERSITY METRICS ####
############################################ 

# Calculate diversity metrics using estimate_richness()
  richnessEstRare<-estimate_richness(physeqRare, split=TRUE, measures= c("Shannon", "Observed"))
  head(richnessEstRare)

#add alpha diversity metrics to metadata
  physeqRareMeta <- data.frame(sample_data(physeqRare))
  head(physeqRareMeta)
  physeqRareMeta$Shannon <- richnessEstRare$Shannon
  physeqRareMeta$Observed <- richnessEstRare$Observed

  str(data.frame(physeqRareMeta))

  #write.csv(physeqRareMeta, "AlphaDiversityMetadata.csv")


############################  
### ALPHA DIVERSITY DATA ###
############################  

  alphaDat<- read.csv("AlphaDiversityMetadata.csv", row.names = 1)
  str(alphaDat)

  alphaDat$BirdID<- as.factor(alphaDat$BirdID)
  
  table(alphaDat$Status2)




# frequency per bird #

  FreqBird<- data.frame(table(alphaDat$BirdID))
  colnames(FreqBird)<- c("BirdID", "Frequency")
  str(FreqBird)
  mean(FreqBird$Frequency) #1.69
  print(sd(FreqBird$Frequency)/sqrt(length(FreqBird$Frequency))) #0.05 SE
  
  
  FreqBirdMulti<- subset(FreqBird, Frequency >1)
  length(FreqBirdMulti$BirdID) #129 have > 1 sample
  129/273 # 47% of adult birds have >1 sample
  FreqBirdMulti2<- subset(FreqBird, Frequency >2)
  length(FreqBirdMulti2$BirdID) # 43 have more than 2
  43/273
  max(FreqBirdMulti2$Frequency) # max is 6

# number of samples from birds >6 #

  str(alphaDat)
  AlphaDatModel6<-subset(alphaDat, (AgeYears) >= 6)
  str(AlphaDatModel6) #98 samples
  BirdsWithSamplesOver6<-unique(AlphaDatModel6$BirdID) #66
  BirdsOver6Died<- subset(AlphaDatModel6, TerminalYear=="yes")
  length(unique(BirdsOver6Died$BirdID)) #24
  FreqBirdsSenescent<- data.frame(table(AlphaDatModel6$BirdID))
  BirdsSenecentMulti<- subset(FreqBirdsSenescent, Freq>1)
  str(BirdsSenecentMulti) #24 individuals have samples taken longitudinally post senescence only
  
  AlphaDat_longitudinal<- subset(alphaDat, BirdID %in% BirdsWithSamplesOver6)
  FreqBirdsOver6<-data.frame(table(AlphaDat_longitudinal$BirdID)) # including longitudinal samples pre/post senescence
  head(FreqBirdsOver6)
  BirdsOver6Multi<-subset(FreqBirdsOver6, Freq>1)
  str(BirdsOver6Multi)
  37/66 #56% of birds have >1 sample
  max(BirdsOver6Multi$Freq)
  max(alphaDat$AgeYears)
  
  BirdsMultiIDs<- BirdsOver6Multi$Var1
  length(BirdsMultiIDs)
  BirdsMultiDat<-subset(AlphaDat_longitudinal, BirdID %in% BirdsMultiIDs)
  BirdsMultiLess6<- subset(BirdsMultiDat, AgeYears <6)
  length(unique(BirdsMultiLess6$BirdID)) # 19 have samples pre and post senescence
  19/66 #29%
  BirdMultiMore6<- subset(BirdsMultiDat, AgeYears >6)
  length(unique(BirdMultiMore6$BirdID))
  
  
  str(AlphaDat_longitudinal)
  AlphaDat_longitudinal$SampleDate<- as.Date(AlphaDat_longitudinal$SampleDate, format="%d/%m/%Y")
  AlphaDat_longitudinalDates<- AlphaDat_longitudinal[,c("BirdID", "AgeYears", "SampleDate")]
  library(dplyr)
  Earliest<- AlphaDat_longitudinalDates %>%
    group_by(BirdID) %>%
    filter(SampleDate ==min(SampleDate))
  Earliest<- data.frame(Earliest)
  colnames(Earliest)<- c("BirdID", "AgeYearsMin", "DateMin")
  
  Latest<- AlphaDat_longitudinalDates %>%
    group_by(BirdID) %>%
    filter(SampleDate ==max(SampleDate))
  Latest<- data.frame(Latest)
  colnames(Latest)<- c("BirdID", "AgeYearsMax", "DateMax")
  
  DatesSpanned<- merge(Earliest, Latest, by="BirdID")
  head(DatesSpanned)
  str(DatesSpanned) #66 birds
  DatesSpanned$TimePassed<- difftime(DatesSpanned$DateMax,DatesSpanned$DateMin, units="days")
  DatesSpanned$YearPassed<- as.numeric(DatesSpanned$TimePassed/365.25)
  mean(DatesSpanned$YearPassed) #1.4 years
  print(sd(DatesSpanned$YearPassed)/sqrt(length((DatesSpanned$YearPassed))))
  min(DatesSpanned$YearPassed) #0
  max(DatesSpanned$YearPassed) #4.8
  hist(DatesSpanned$YearPassed)
  
  table(alphaDat$TerminalYear) #236 no 106 yes
  TermSamples<- subset(alphaDat, TerminalYear=="yes")
  range(TermSamples$AgeYears)
  
# Age stats for paper
  
  min(alphaDat$AgeYears)
  max(alphaDat$AgeYears)
  mean(alphaDat$AgeYears)
  sd(alphaDat$AgeYears)
  
#females
  alphaDat$SexEstimate
  alphaDatF<- subset(alphaDat, SexEstimate==0)
  str(alphaDatF) #211 female samples
  length(unique(alphaDatF$BirdID))
  
  min(alphaDatF$AgeYears)
  max(alphaDatF$AgeYears)
  mean(alphaDatF$AgeYears)
  sd(alphaDatF$AgeYears)
  print(sd(alphaDatF$AgeYears)/sqrt(length(alphaDatF$AgeYears))) #0.21 SE

  # terminal samples F
  
  alphaDatFterm<- subset(alphaDatF, TerminalYear=="yes")
  str(alphaDatFterm)
  range(alphaDatFterm$AgeYears) # 1-15.6
  mean(alphaDatFterm$AgeYears) # 5.09
  print(sd(alphaDatFterm$AgeYears)/sqrt(length(alphaDatFterm$AgeYears))) # SE 0.6
  min(alphaDatFterm$AgeYears) # 1
  max(alphaDatFterm$AgeYears) #15.6
  
  # number of samples from F birds >6
  
  str(alphaDatF)
  alphaDatFOld<-subset(alphaDatF, (AgeYears) >= 6)
  str(alphaDatFOld) #43 samples
  length(unique(alphaDatFOld$BirdID)) #30
  BirdsWithSamplesOver6F<-unique(alphaDatFOld$BirdID)
  alphaDatFOld_long<- subset(alphaDatF, BirdID %in% BirdsWithSamplesOver6F)
  FreqBirdsOver6F<-data.frame(table(alphaDatFOld_long$BirdID))
  head(FreqBirdsOver6F)
  BirdsOver6MultiF<-subset(FreqBirdsOver6F, Freq>1)
  str(BirdsOver6MultiF)
  14/30 #47% of F birds over 6 have >1 sample
  max(BirdsOver6MultiF$Freq)
  
  
  
#males
  alphaDat$SexEstimate
  alphaDatM<- subset(alphaDat, SexEstimate==1)
  str(alphaDatM) #251 males
  length(unique(alphaDatM$BirdID))
  
  min(alphaDatM$AgeYears)
  max(alphaDatM$AgeYears)
  mean(alphaDatM$AgeYears)
  sd(alphaDatM$AgeYears)
  print(sd(alphaDatM$AgeYears)/sqrt(length(alphaDatM$AgeYears))) #0.20 SE
  
  
  alphaDatMterm<- subset(alphaDatM, TerminalYear=="yes")
  str(alphaDatMterm)
  range(alphaDatMterm$AgeYears) #1-17.2
  mean(alphaDatMterm$AgeYears)
  print(sd(alphaDatMterm$AgeYears)/sqrt(length(alphaDatMterm$AgeYears))) 
  
  # number of samples from M birds >6
  
  str(alphaDatM)
  alphaDatMOld<-subset(alphaDatM, (AgeYears) >= 6)
  str(alphaDatMOld) #55 samples
  length(unique(alphaDatMOld$BirdID)) #36
  BirdsWithSamplesOver6M<-unique(alphaDatMOld$BirdID)
  alphaDatMOld_long<- subset(alphaDatM, BirdID %in% BirdsWithSamplesOver6M)
  FreqBirdsOver6M<-data.frame(table(alphaDatMOld_long$BirdID))
  head(FreqBirdsOver6M)
  BirdsOver6MultiM<-subset(FreqBirdsOver6M, Freq>1)
  str(BirdsOver6MultiM)
  23/36 #63% of F birds over 6 have >1 sample

  
### numbers with estimated birth date ###
  
  
  str(alphaDat)
  BirdsBD<- alphaDat[,c("BirdID", "BirthDate", "BirthDateEstimated")]
  BirdsBD<- unique(BirdsBD)
  str(BirdsBD) # 273 birds
  
  table(BirdsBD$BirthDateEstimated) # 100 false, 173 true (63%)
  
  EstSamplesPerFP<-alphaDat %>% 
    group_by(SampleYear) %>%
    count(BirthDateEstimated) %>%
    mutate(prop = prop.table(n))
  
  BD_Yr<- ggplot(alphaDat, aes(x = as.factor(SampleYear),fill = as.factor(BirthDateEstimated))) + 
    geom_bar(position = "fill")+
    xlab("Sample year") + ylab("Proportion with estimated birth date") +
    scale_fill_manual(values=c("grey25", "grey60")) +
    guides(fill=guide_legend(title="Birth date estimated"))
  
  Sexlabel<-c("Female", "Male")
  BD_Sex<- ggplot(alphaDat, aes(x = as.factor(SexEstimate),fill = as.factor(BirthDateEstimated))) + 
    geom_bar(position = "fill") +
    xlab("Sex")+
    theme(axis.title.y =element_blank(),
          axis.text.y =element_blank(),
          axis.ticks.y =element_blank())+
    scale_fill_manual(values=c("grey25", "grey60")) +
    guides(fill=guide_legend(title="Birth date estimated"))+
    scale_x_discrete(labels= Sexlabel)
  
  BD_Season<- ggplot(alphaDat, aes(x = as.factor(Season),fill = as.factor(BirthDateEstimated))) + 
    geom_bar(position = "fill")+
    xlab("Season")+ ylab("Proportion with estimated birth date")+
    scale_fill_manual(values=c("grey25", "grey60")) +
    guides(fill=guide_legend(title="Birth date estimated"))
  
  Agelabel<-c("1-3", "3-6", ">6")
  BD_Age<- ggplot(alphaDat, aes(x = as.factor(Ageclass2),fill = as.factor(BirthDateEstimated))) + 
    geom_bar(position = "fill")+
    theme(axis.title.y =element_blank(),
          axis.text.y =element_blank(),
          axis.ticks.y =element_blank())+
    xlab("Age (years)")+ 
    scale_fill_manual(values=c("grey25", "grey60")) +
    guides(fill=guide_legend(title="Birth date estimated"))+
    scale_x_discrete(labels= Agelabel)
  
  
  plot(BD_Yr + BD_Sex + BD_Season + BD_Age) + plot_layout(guides = "collect")

  
#### Time at 4 degrees stats ####
  
# Info about the variable for paper
  
  hist(alphaDat$TimeAt4Degrees)
  
  range(alphaDat$TimeAt4Degrees) #0-115 days (~4 months max)
  
  alphaDat %>%
    group_by(Season) %>%
    reframe(RangeDays= range(TimeAt4Degrees)) # major= 0-115 days, minor= 7-74 days
  
  
  mean(alphaDat$TimeAt4Degrees) #43.76 days mean across all samples
  
  alphaDat %>%
    group_by(Season) %>%
    summarise(MeanDays= mean(TimeAt4Degrees)) # major= 47.4 days, minor= 36.1 days
  
  SE <- function(x) sd(x)/sqrt(length(x))
  SE(alphaDat$TimeAt4Degrees) # SE= 1.083 all samples
  
  alphaDat %>%
    group_by(Season) %>%
    summarise(SEDays= SE(TimeAt4Degrees)) # 1.42 major, 1.33 minor SE
  
  
  median(alphaDat$TimeAt4Degrees) # 41 days all samples
  
  alphaDat %>%
    group_by(Season) %>%
    summarise(MedianDays= median(TimeAt4Degrees)) # major= 46 days, minor= 35 days
  

### check for correlations with other variables
  
  CorAge<- ggscatter(alphaDat, x = "TimeAt4Degrees", y = "AgeYears", 
            add = "reg.line", conf.int = TRUE, color="grey40",
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Time at 4 degrees", ylab = "AgeYears") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  CorTimeDay<- ggscatter(alphaDat, x = "TimeAt4Degrees", y = "MinutesSinceSunrise", 
                         add = "reg.line", conf.int = TRUE, color = "grey40",
                         cor.coef = TRUE, cor.method = "pearson",
                         xlab = "Time at 4 degrees", ylab = "Time of day (minutes since sunrise)") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  
  CorMeanTQ<- ggscatter(alphaDat, x = "TimeAt4Degrees", y = "MeanTQ_FP", 
            add = "reg.line", conf.int = TRUE, color="grey40",
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Time at 4 degrees", ylab = "Mean territory quality") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  CorCentredTQ<- ggscatter(alphaDat, x = "TimeAt4Degrees", y = "MeanCentredTQ", 
            add = "reg.line", conf.int = TRUE, color="grey40",
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "Time at 4 degrees", ylab = "Centred territory quality") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  
  

  alphaDat$SeasonBinary<- ifelse(alphaDat$Season=="Major", 1, 0)
  cor.test(alphaDat$SeasonBinary, alphaDat$TimeAt4Degrees) # 0.23 point biseral correlation
  CorSeason<- ggscatter(alphaDat, x = "SeasonBinary", y = "TimeAt4Degrees",
                           cor.coef = TRUE, cor.method = "pearson",
                           xlab = "Season", ylab = "Time at 4 degrees", point=FALSE) +
    scale_x_continuous(labels=c("Minor", "Major"), breaks=0:1)+
    geom_jitter(width=0.1, colour="grey40") +
    geom_boxplot(aes(group=Season), alpha=0.4)+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))

 
  
  cor.test(alphaDat$SexEstimate, alphaDat$TimeAt4Degrees) #0.13 sex point biseral correlation
  CorSex<- ggscatter(alphaDat, x = "SexEstimate", y = "TimeAt4Degrees",
                        cor.coef = TRUE, cor.method = "pearson",
                        xlab = "Sex", ylab = "Time at 4 degrees", point=FALSE) +
    scale_x_continuous(labels=c("Female", "Male"), breaks=0:1)+
    geom_jitter(width=0.1, colour="grey40") +
    geom_boxplot(aes(group=SexEstimate), alpha=0.4)+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))

  
  
  
  plot((CorAge + CorTimeDay +CorMeanTQ) / (CorCentredTQ +CorSeason + CorSex)) 

  
  
  
  
  
  


  
##################### 
### SHANNON MODEL ###
#####################

  alphaDatModel<- alphaDat
  str(alphaDatModel) #462


### Make categorical variables factors ###
  
  alphaDatModel$SexEstimate<- ifelse(alphaDatModel$SexEstimate == 0, "F", "M")
  alphaDatModel$SexEstimate<- as.factor(alphaDatModel$SexEstimate)
  alphaDatModel$Season<- as.factor(alphaDatModel$Season)
  alphaDatModel$BirdID<- as.factor(alphaDatModel$BirdID)
  alphaDatModel$Status2<- as.factor(alphaDatModel$Status2)
  alphaDatModel$SurvivedToNextSeason<- as.factor(alphaDatModel$SurvivedToNextSeason)
  alphaDatModel$SampleDate<- as.Date(alphaDatModel$SampleDate,format="%d/%m/%Y")
  alphaDatModel$TerminalYear<- as.factor(alphaDatModel$TerminalYear)
  table(alphaDatModel$TerminalYear) # 106 samples taken in the terminal year, 356 not

  TerminalIDs<- subset(alphaDatModel, TerminalYear=="yes")
  TerminalBirds<- length(unique(TerminalIDs$BirdID)) #91 birds
  

#plot of sampling regime

  SamplePlot<- ggplot(alphaDatModel, aes(AgeYears, BirdID)) +
    geom_vline(xintercept = 1, linetype="dashed", color = "darkblue")+
    geom_vline(xintercept = 3, linetype="dashed", color = "darkblue")+
    geom_vline(xintercept = 6, linetype="dashed", color = "darkblue")+
    geom_line(aes(group=BirdID), col="grey50")+
    geom_point(aes(col=TerminalYear)) +
    scale_colour_manual(values=c("black", "goldenrod"), name="Terminal year")+
    scale_x_continuous(limits=c(0,18), n.breaks=20)+
    xlab("\n Age (years)") +
    theme_bw() +
    theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())+
    theme(axis.text = element_text(size=16), axis.title = element_text(size=20))+
    theme(legend.text = element_text(size=16), legend.title = element_text(size=18)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  SamplePlot


### rescale the variables
  alphaDatModel$AgeYears.sc<- rescale(alphaDatModel$AgeYears, binary.inputs = "full")
  alphaDatModel$MeanTQ_FP.sc<- rescale(alphaDatModel$MeanTQ_FP, binary.inputs = "full")
  alphaDatModel$CentredTQ.sc<- rescale(alphaDatModel$MeanCentredTQ, binary.inputs = "full")
  alphaDatModel$TimeOfDay.sc<- rescale(alphaDatModel$MinutesSinceSunrise, binary.inputs = "full")
  alphaDatModel$TimeToFreeze.sc <- rescale(alphaDatModel$TimeAt4Degrees, binary.inputs = "full")
 
  
### need to make sex and terminal year ordered factors to be able to test interactions in mgcv ###

  alphaDatModel$SexEstimate<- factor(alphaDatModel$SexEstimate, levels = c(
    "F", "M"), ordered = TRUE) 
  alphaDatModel$TerminalYear<- factor(alphaDatModel$TerminalYear, levels = c(
    "yes", "no"), ordered = TRUE) 

  str(alphaDatModel)
  

#### a) check the vifs with no interaction terms present ####
  
  car::vif(lmer(Shannon~ AgeYears.sc+ TimeOfDay.sc + TimeToFreeze.sc + SexEstimate+Season+ 
             CentredTQ.sc + MeanTQ_FP.sc + TerminalYear + (1|BirdID), data=alphaDatModel))

  gammModelS_A<- mgcv::gam(Shannon ~ s(AgeYears.sc, bs="cr", k=10) + 
                           s(TimeOfDay.sc, bs="cr") + 
                           s(TimeToFreeze.sc, bs="cr") +
                           SexEstimate + 
                           Season +
                           CentredTQ.sc +
                           MeanTQ_FP.sc + 
                           TerminalYear +
                           s(BirdID, bs="re"),
                         data=alphaDatModel, 
                         family=gaussian,
                         method="REML")

  summary(gammModelS_A) #summary of model
  mgcv::gam.vcomp(gammModelS_A) # variance components incl. for random effects
  visreg(gammModelS_A, "AgeYears.sc")
  visreg(gammModelS_A, "MeanTQ_FP.sc")

  
  

  

## function to calculate vifs from gam from sam clifford github (couldn't get mgcv.helper::vif.gam() to work) ##
  
  gamvif <- function(mod){
    # mod is an mgcv object
    mod.sum <- summary(mod)
    s2 <- mod$sig2 # estimate of standard deviation of residuals
    X <- model.matrix(mod) # data used to fit the model
    n <- nrow(X) # how many observations were used in fitting?
    v <- -1 # omit the intercept term, it can't inflate variance
    varbeta <- mod.sum$p.table[v,2]^2 # variance in estimates
    varXj <- apply(X=X[,row.names(mod.sum$p.table)[v]],MARGIN=2, var) # variance of all the explanatory variables
    VIF <- varbeta/(s2/(n-1)*1/varXj) # the variance inflation factor, obtained by rearranging
    # var(beta_j) = s^2/(n-1) * 1/var(X_j) * VIF_j
    VIF.df <- data.frame(variable=names(VIF),
                         vif=VIF, 
                         row.names=NULL)
    
    return(VIF.df)
  }

  gamvif(gammModelS_A) # all VIFS <2 for parametric terms
  concurvity(gammModelS_A) # when random effects excluded all fine


## model diagnostics ##

  gam.check(gammModelS_A)
# the k-index:  1 for Age, 1.10 for Time of day, 1 for time to freeze
# The further below 1 this is, the more likely it is that there is missed pattern left in the residuals.

#Low p-values may indicate that residuals are not randomly distrubuted and
# that the basis dimension, k, has been set too low, especially if the reported edf is close to k', the maximum possible EDF for the term.
# p-values are 0.47 Age, 0.98 Time of day, 0.52 time to freeze=  all looks fine

# plots look fine




### b) interaction terms for sex and terminal year ###

  gammModelS_A2<- mgcv::gam(Shannon ~ s(AgeYears.sc, bs="cr", k=10) +
                            s(AgeYears.sc, bs="cr", by=TerminalYear, m=1, k=10) + 
                            s(AgeYears.sc, bs="cr", by=SexEstimate, m=1, k=10)+
                            s(TimeOfDay.sc, bs="cr") + 
                            s(TimeToFreeze.sc, bs="cr") +
                            SexEstimate + 
                            Season +
                            CentredTQ.sc + 
                            MeanTQ_FP.sc + 
                            TerminalYear +
                            s(BirdID, bs="re"),
                          data=alphaDatModel, 
                          family=gaussian,
                          method="REML")

  summary(gammModelS_A2) #summary of model- neither sex*age or terminal sample*age sig
  AIC(gammModelS_A2, gammModelS_A) # interactions not sig and higher AIC



### c) remove sex interaction term and retest ###

  gammModelS_A3<- mgcv::gam(Shannon ~ s(AgeYears.sc, bs="cr", k=10) +
                            s(AgeYears.sc, bs="cr", by=TerminalYear, m=1, k=10) + 
                            s(TimeOfDay.sc, bs="cr") + 
                            s(TimeToFreeze.sc, bs="cr") +
                            SexEstimate + 
                            Season +
                            CentredTQ.sc + 
                            MeanTQ_FP.sc + 
                            TerminalYear +
                            s(BirdID, bs="re"),
                          data=alphaDatModel, 
                          family=gaussian,
                          method="REML")

  summary(gammModelS_A3) #summary of model- terminal sample*age not sig

  AIC(gammModelS_A3, gammModelS_A) # interaction not sig and higher AIC

  
# export results with interactions 
  resultsS_Ainteractions<- summary(gammModelS_A2)
  sink("results_terminalShannonA_interactions.txt")
  resultsS_Ainteractions$p.table
  resultsS_Ainteractions$s.table
  resultsS_Ainteractions$r.sq
  mgcv::gam.vcomp(gammModelS_A2)
  sink()
  

# export results of model with no interaction
  resultsS_Aterminal<- summary(gammModelS_A)
  sink("results_terminalShannonA.txt")
  resultsS_Aterminal$p.table
  resultsS_Aterminal$s.table
  resultsS_Aterminal$r.sq
  mgcv::gam.vcomp(gammModelS_A) # variance components incl. for random effects
  sink()



#### c) plot the results ####

# i) for time to freeze:

# Extract the prediction data frame
  predFreeze <- data.frame(ggpredict(gammModelS_A, terms = c("TimeToFreeze.sc[all]")))  # this gives overall predictions for the model
  str(predFreeze)
  head(predFreeze)

#rescale
  mean(alphaDatModel$TimeAt4Degrees) # 43.75541
  2*sd(alphaDatModel$TimeAt4Degrees) # 46.56452

  predFreeze$TimeToFreeze_rescale<- (predFreeze$x*46.56452)+43.75541
  head(predFreeze)

# Plot the predictions 

  ShannonFreeze<- ggplot(predFreeze) + 
    geom_point(data = alphaDatModel, aes(x = TimeAt4Degrees, y = Shannon),col="grey35") + 
    theme_bw()+
    geom_ribbon(aes(x = TimeToFreeze_rescale, ymin = conf.low, ymax = conf.high), 
                fill = "lightgrey", alpha = 0.5) +
    geom_line(aes(x = TimeToFreeze_rescale, y = predicted), linewidth=2) +          
    labs(x = "\nTime at 4 degrees (days)", y = "Shannon diversity\n") +
    theme(axis.text = element_text(size=14), axis.title = element_text(size=20)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  ShannonFreeze


  
  


#####################    
### Observed ASVs ###
#####################

  mean(alphaDatModel$Observed)
  var(alphaDatModel$Observed) # variance much greater than mean: overdispersed count data- use neg.bin


####### Negative binomial GAMM model- Observed ASVs ########

### terminal effect in adults #####

  str(alphaDatModel) #462
  length(unique(alphaDatModel$BirdID)) #273

  gammModelO_adults<- mgcv::gam(Observed ~ s(AgeYears.sc, bs="cr", k=10) + 
                                s(TimeOfDay.sc, bs="cr") + 
                                s(TimeToFreeze.sc, bs="cr") +
                                SexEstimate + 
                                Season +
                                CentredTQ.sc + 
                                MeanTQ_FP.sc +
                                TerminalYear +
                                s(BirdID, bs="re"),
                              data=alphaDatModel, 
                              family="nb", 
                              method="REML")

  summary(gammModelO_adults) #summary of model
  mgcv::gam.vcomp(gammModelO_adults) # variance components incl. for random effects

  gamvif(gammModelO_adults) # all VIFS <2 for parametric terms
  concurvity(gammModelO_adults) # when random effects excluded all less than 0.2 so fine


## model diagnostics ##

  gam.check(gammModelO_adults)
# the k-index:  0.97 for Age, 1.03 for Time of day, 0.94 for time to freeze
# The further below 1 this is, the more likely it is that there is missed pattern left in the residuals.

#Low p-values may indicate that residuals are not randomly distributed and
# that the basis dimension, k, has been set too low, especially if the reported edf is close to k', the maximum possible EDF for the term.
# p-values are 0.60 Age, 0.91 Time of day, 0.34 time to freeze - all fine

# diagnostic plots are fine


### b) with interaction terms ###

  str(alphaDatModel) # check factors for interactions are ordered 

  gammModelO2_adults<- mgcv::gam(Observed ~ s(AgeYears.sc, bs="cr", k=10)+
                                 s(AgeYears.sc, bs="cr", by=TerminalYear, m=1, k=10) + 
                                 s(AgeYears.sc, bs="cr", by=SexEstimate, m=1, k=10) + 
                                 s(TimeOfDay.sc, bs="cr") + 
                                 s(TimeToFreeze.sc, bs="cr") +
                                 SexEstimate + 
                                 CentredTQ.sc + 
                                 MeanTQ_FP.sc +
                                 Season +
                                 TerminalYear +
                                 s(BirdID, bs="re"),
                               data=alphaDatModel, 
                               family=nb, 
                               method="REML")

  summary(gammModelO2_adults) #summary of model
  gam.check(gammModelO2_adults)

  AIC(gammModelO_adults, gammModelO2_adults) #  AICs are worse with interactions and neither sig



### c) drop sex interaction term and retest ###

  gammModelO3_adults<- mgcv::gam(Observed ~ s(AgeYears.sc, bs="cr", k=10)+
                                 s(AgeYears.sc, bs="cr", by=TerminalYear, m=1, k=10) + 
                                 s(TimeOfDay.sc, bs="cr") + 
                                 s(TimeToFreeze.sc, bs="cr") +
                                 SexEstimate + 
                                 CentredTQ.sc + 
                                 MeanTQ_FP.sc +
                                 Season +
                                 TerminalYear +
                                 s(BirdID, bs="re"),
                               data=alphaDatModel, 
                               family=nb, 
                               method="REML")

  summary(gammModelO3_adults) #summary of model
  gam.check(gammModelO3_adults)
  AIC(gammModelO_adults, gammModelO3_adults) # AIC worse and term not sig


## output results of model without interaction terms for supp##
  
  resultsObserved_Ainteractions<- summary(gammModelO2_adults)
  sink("ObservedResults_Ainteractions.txt")
  resultsObserved_Ainteractions$p.table
  resultsObserved_Ainteractions$s.table
  resultsObserved_Ainteractions$r.sq
  mgcv::gam.vcomp(gammModelO2_adults)
  sink()
  
## output results of model without interaction terms ##

  resultsObserved_A<- summary(gammModelO_adults)
  sink("ObservedResults_adultsOnly.txt")
  resultsObserved_A$p.table
  resultsObserved_A$s.table
  resultsObserved_A$r.sq
  mgcv::gam.vcomp(gammModelO_adults)
  sink()


#### c) plot the results ####

## for sex

# Extract the prediction data frame
  predSexObs <- data.frame(ggpredict(gammModelO_adults, terms = c("SexEstimate")))  # this gives overall predictions for the model
  str(predSexObs)
  head(predSexObs)
  colnames(predSexObs)<- c("SexEstimate", "predicted", "std.error", "conf.low" , "conf.high", "group")


# Plot the predictions 

  SexPredictionPlotObserved<- ggplot(predSexObs) + 
  geom_violin(data= alphaDatModel,aes(x = SexEstimate, y = Observed))+
  geom_point(data = alphaDatModel,aes(x = SexEstimate, y = Observed),
             position= position_jitter(width=0.1), col="grey55", alpha=0.3 ) + 
  theme_bw()+
  theme(axis.title.y = element_blank())+
  geom_linerange(aes(x=SexEstimate,ymin=conf.low, ymax=conf.high), col="black", linewidth=2)+
  geom_point(aes(x=SexEstimate, y=predicted), col="red", shape=19, size=5) +
  labs(x = "\nSex", y = "Observed ASVs\n")+
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  SexPredictionPlotObserved

