setwd("/users/sarah/OneDrive - University of East Anglia/Seychelles_warbler/Analysis4_toS2022/MergedQIIME/Age_paper/Age/paper/FinalDraft/MolEcol/DRYAD/Similarity")

#load packages
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(brms)
  library(rstan)
  library(ggdist)
  library(ggplot2)
  library(bayesplot)
  library(patchwork)


#1. Making microbiome pairwise similarity matrices:

# load phyloseq object- clr transformed abundances for just adults 
# (note this is generated in previous beta analysis r script)- ASVs in <5% of sample removed
  
  physeq_clr_A<- readRDS("physeq_clr_A.rds") #462 samples, 624 taxa
  physeq_clr_A
  
  length(unique((sample_data(physeq_clr_A)$BirdID))) #273 birds
  
  FreqSamples<- data.frame(table(sample_data(physeq_clr_A)$BirdID))
  head(FreqSamples)
  FreqSamplesMulti<- subset(FreqSamples, Freq>1)
  length(FreqSamplesMulti$Var1) #129 birds have more than one sample
  
  sample_data<-data.frame(sample_data(physeq_clr_A))
  
  str(sample_data)
  table(sample_data$TerminalYear) # 106 yes, 356 no
  table(sample_data$Ageclass2)
  
  sample_data<- sample_data[,c("sample.id","BirdID","SexEstimate","Ageclass2", "SampleDate", "TerminalYear", "AgeYears", "DiedEver", "TimeToDeath", "Season")]
                    
  saveRDS(sample_data,"sample_data_from_phyloseq.rds")
  
  sample_data<-readRDS("sample_data_from_phyloseq.rds")
  
                    
#make a key for the order of sample names and their associated individual IDs.
  
  key<-data.frame(ID=sample_data(physeq_clr_A)$BirdID, Sample_name=sample_data(physeq_clr_A)$sample.id)
                    
  

###### Make Aitchisons similarity matrix ########

 
#function to extract ASV matrix
  vegan_otu <- function(physeq) {
    OTU <- otu_table(physeq)
    if (taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    return(as(OTU, "matrix"))
  }
  
  
#Extract ASV Matrix and calculated euclidean distances
  clr_Matrix<-vegan_otu(physeq_clr_A)
  
  AitM<- as.matrix(vegdist(clr_Matrix, method="euclidean"))

                    
#Note that distance matrix has rownames and colnames in the same order as key
  all(rownames(AitM)==key$Sample_name)

#Eyeball distances across individuals
  AitM[1:3]
  hist(AitM)
  
#convert distances to similarity values (ranging from 0 to 1 where 0 is the same, 1 is the most different pair)
  max(AitM)
  AitM<- 1-(AitM/max(AitM))
  
  AitM[1:3] #check
  hist(AitM)
  all(rownames(AitM)==key$Sample_name)
                    
#Save matrix to ready matrices folder
  saveRDS(AitM,"AitM.rds")
  

  
###### 2. Making continuous temporal distance matrices ######
                    
### Temporal distance matrix 
### describes the distance in days between microbiome samples based on the sample collection dates recorded in sample_data

#Transform dates into a numeric variable
  sample_data$date<-as.Date(as.character(sample_data$SampleDate), format="%d/%m/%Y")
  sample_data$date_numeric<-as.numeric(sample_data$date)
  str(sample_data)
  
#Create data frame with each sample name (character) and sampling time (numeric) 
  SampleTime_frame<-sample_data[,c("sample.id","date_numeric")]
  str(SampleTime_frame)
  rownames(SampleTime_frame)

#Create an empty matrix to fill with distances
  TEMPM<-array(0,c(nrow(SampleTime_frame),nrow(SampleTime_frame)))

#Derive matrix with time distances between each sample using abs()-function
  for (i in 1:nrow(SampleTime_frame)){
  for (j in 1:nrow(SampleTime_frame))
    {TEMPM[i,j]=abs(SampleTime_frame$date_numeric[i] -SampleTime_frame$date_numeric[j])
  }
}

#Note that Temporal distance matrix has rownames and colnames in the same order as key

#Name rown amd colnames with sample names 
  rownames(TEMPM)<-key$Sample_name
  colnames(TEMPM)<-key$Sample_name

#Save matrix
  saveRDS(TEMPM,"TEMPM.rds")


  
###### 3. Making binary similarity matrices (here, age similarity) ######
  
## The resulting matrix will have a value of "1"= same age or "0"=different age for each pair


#Create data frame with each sample name (character) and their Age (Character)
  Age_frame<-sample_data[,c("sample.id","Ageclass2")]
  str(Age_frame)
  Age_frame$Ageclass2<-as.character(Age_frame$Ageclass2)

#Create an empty numeric matrix to fill with distances
  AGEM<-array(0,c(nrow(Age_frame),nrow(Age_frame)))

#Derive matrix with binary Age similarity between each sample
  for(i in 1:nrow(Age_frame)){
  for(j in 1:nrow(Age_frame)){
    if(Age_frame$Ageclass2[i]==Age_frame$Ageclass2[j]){
      AGEM[i,j]= 1
      } else{
        AGEM[i,j]= 0
      }
  }
  } 

#Note that AGE similarity matrix has rownames and colnames in the same order as key
  all(rownames(AGEM)==key$Sample_name)
  rownames(AGEM)<-key$Sample_name
  colnames(AGEM)<-key$Sample_name

#Save matrix
  saveRDS(AGEM,"AGEM.rds")
  

  
###### 4. Making binary similarity matrices (here, season similarity) ######
  
## The resulting matrix will have a value of "1"= same season (Major-Major, Minor-Minor) or "0"=different season for each pair
  
  
#Create data frame with each sample name (character) and their Age (Character)
  Season_frame<-sample_data[,c("sample.id","Season")]
  str(Season_frame)
  Season_frame$Season <-as.character(Season_frame$Season)
  
  #Create an empty numeric matrix to fill with distances
  SEASONM<-array(0,c(nrow(Season_frame),nrow(Season_frame)))
  
  #Derive matrix with binary Age similarity between each sample
  for(i in 1:nrow(Season_frame)){
    for(j in 1:nrow(Season_frame)){
      if(Season_frame$Season[i]==Season_frame$Season[j]){
        SEASONM[i,j]= 1
      } else{
        SEASONM[i,j]= 0
      }
    }
  } 
  
  #Note that AGE similarity matrix has rownames and colnames in the same order as key
  all(rownames(SEASONM)==key$Sample_name)
  rownames(SEASONM)<-key$Sample_name
  colnames(SEASONM)<-key$Sample_name
  
  #Save matrix to ready matrices folder
  saveRDS(SEASONM,"SEASONM.rds")

  
###### 5. Making binary similarity matrices (here, sex similarity) ########
  
##### The resulting matrix will have a value of "1"= same sex or "0"=different sex for each pair
  
  
#Create data frame with each sample name (character) and their Age (Character)
  Sex_frame<-sample_data[,c("sample.id","SexEstimate")]
  str(Sex_frame)
  Sex_frame$SexEstimate<-as.character(Sex_frame$SexEstimate)
  
#Create an empty numeric matrix to fill with distances
  SEXM<-array(0,c(nrow(Sex_frame),nrow(Sex_frame)))
  
  #Derive matrix with binary Age similarity between each sample
  for(i in 1:nrow(Sex_frame)){
    for(j in 1:nrow(Sex_frame)){
      if(Sex_frame$SexEstimate[i]==Sex_frame$SexEstimate[j]){
        SEXM[i,j]= 1
      } else{
        SEXM[i,j]= 0
      }
    }
  } 
  
#Note that AGE similarity matrix has rownames and colnames in the same order as key
  all(rownames(SEXM)==key$Sample_name)
  rownames(SEXM)<-key$Sample_name
  colnames(SEXM)<-key$Sample_name
  
#Save matrix to ready matrices folder
  saveRDS(SEXM,"SEXM.rds")
  
  

  
###### 6. Making combination-factor matrices (here, ageclass combination of a pair) ######
  
# The resulting matrix will have for each individual pair a value of:
# "YY"= both young adult, "MM"= both middle adult, or "OO"=both old adult.
# Also YM, MO, YO comparisons

  
#Create data frame with each sample name (character) and their Ageclass (Character)
  Ageclass_frame<-sample_data[,c("sample.id","Ageclass2")]
  #rename factor levels
  Ageclass_frame$Ageclass2<- recode_factor(Ageclass_frame$Ageclass2, "A 1-3 yrs" = "Y", "A 3-6 yrs" = "M", "Old A" ="O")
  levels(Ageclass_frame$Ageclass2)
  Ageclass_frame$Ageclass2<-as.character(Ageclass_frame$Ageclass2)
  str(Ageclass_frame)
  
#Create an empty character matrix to fill with characters
  AGECLASSM<-array(as.character(NA),c(nrow(Ageclass_frame),nrow(Ageclass_frame)))
  
  for(i in 1:nrow(Ageclass_frame)){
    for(j in 1:nrow(Ageclass_frame)){ 
      if(Ageclass_frame$Ageclass2[i]=="Y" & Ageclass_frame$Ageclass2[i]==Ageclass_frame$Ageclass2[j]){
        AGECLASSM[i,j]= "YY"}
      if(Ageclass_frame$Ageclass2[i]=="M" & Ageclass_frame$Ageclass2[i]==Ageclass_frame$Ageclass2[j]){
        AGECLASSM[i,j]= "MM"}
      if(Ageclass_frame$Ageclass2[i]=="O" & Ageclass_frame$Ageclass2[i]==Ageclass_frame$Ageclass2[j]){
        AGECLASSM[i,j]= "OO"}
      if(Ageclass_frame$Ageclass2[i]=="Y" & Ageclass_frame$Ageclass2[j]=="M"){
        AGECLASSM[i,j]= "YM"}
      if(Ageclass_frame$Ageclass2[i]=="M" & Ageclass_frame$Ageclass2[j]=="Y"){
        AGECLASSM[i,j]= "YM"}
      if(Ageclass_frame$Ageclass2[i]=="Y" & Ageclass_frame$Ageclass2[j]=="O"){
        AGECLASSM[i,j]= "YO"}
      if(Ageclass_frame$Ageclass2[i]=="O" & Ageclass_frame$Ageclass2[j]=="Y"){
        AGECLASSM[i,j]= "YO"}
      if(Ageclass_frame$Ageclass2[i]=="M" & Ageclass_frame$Ageclass2[j]=="O"){
        AGECLASSM[i,j]= "MO"}
      if(Ageclass_frame$Ageclass2[i]=="O" & Ageclass_frame$Ageclass2[j]=="M"){
        AGECLASSM[i,j]= "MO"}
    }
  }
  
  rownames(AGECLASSM)<-key$Sample_name
  colnames(AGECLASSM)<-key$Sample_name  
  
  saveRDS(AGECLASSM,"AGECLASSM.rds")
  
  
  
###### 6. Making combination-factor matrices (here, terminal sample comparison) #####
  
#The resulting matrix will have for each individual pair a value of:
#"TT"= both terminal; "TN"=terminal v not terminal; "NN"= not terminal not terminal

  #Create data frame with each sample name (character) and their sex
  str(sample_data)
  Terminal_frame<-sample_data[,c("sample.id","TerminalYear")]
  str(Terminal_frame)
  
  #Create an empty character matrix to fill with characters
  TERMM<-array(as.character(NA),c(nrow(Terminal_frame),nrow(Terminal_frame)))
  
  for(i in 1:nrow(Terminal_frame)){
    for(j in 1:nrow(Terminal_frame)){ 
      if(Terminal_frame$TerminalYear[i]=="yes" & Terminal_frame$TerminalYear[i]==Terminal_frame$TerminalYear[j]){
        TERMM[i,j]= "TT"}
      if(Terminal_frame$TerminalYear[i]=="no" & Terminal_frame$TerminalYear[i]==Terminal_frame$TerminalYear[j]){
        TERMM[i,j]= "NN"}
      if( Terminal_frame$TerminalYear[i]!=Terminal_frame$TerminalYear[j]){
        TERMM[i,j]= "NT"}
    }
  }
  
  rownames(TERMM)<-key$Sample_name
  colnames(TERMM)<-key$Sample_name
  
  
  #Save matrix to ready matrices folder
  saveRDS(TERMM,"TERMM.rds")  

  
  
##### These matrices can be "unraveled" into a long dyadic data frame ######
  
  ## where each row depicts one dyad
  ## columns mark the various calculated distances between members of this pair
  ## Note that for the sake of constructing a multimembership random effect in dyadic glms later
  ## also need to add columns giving the identity of sample/individual A and B in each dyad

##### first unravel the existing matrices into one dyadic data frame:
  
#Read in the matrices if not in already:
  AitM<-readRDS("AitM.rds")
  TEMPM<-readRDS("TEMPM.rds")
  AGEM<-readRDS("AGEM.rds")
  SEASONM<- readRDS("SEASONM.rds")
  SEXM<-readRDS("SEXM.rds")
  AGECM<-readRDS("AGECLASSM.rds")
  TERMM<- readRDS("TERMM.rds")

#unravel the matrices into vectors matching the lower quantile of each matrix. 

#From numeric matrices, this can be done by making a list (c()) of the distance object (dist()) derived from the matrix. 
  # as.dist() by default includes only the lower quantile of the matrix and excludes the diagonal.
  
#From categorical matrices, this can be done by making a list (c()) of the lower quantile of the matrix with lower.tri() -function.

  ait<-c(as.dist(AitM)) # similarity matrix
  temp<-c(as.dist(TEMPM)) # number of days between sampling points
  age<-c(as.dist(AGEM)) # same age or different (0 or 1)
  season<- c(as.dist(SEASONM)) # same season or different (0 or 1)
  sex<-c(as.dist(SEXM)) # same sex or different (0 or 1)
  ageClass<-c(AGECM[lower.tri(AGECM)]) # ageclass combination
  terminal<-c(TERMM[lower.tri(TERMM)]) # type of terminal/non-terminal comparison

#Combine these vectors into a data frame
  
  data.dyad<-data.frame(Microbiome_similarity=ait,Temporal_distance=temp,
                        age_similarity=age,sex_similarity=sex,
                        season_similarity= season,
                        ageclass_combination=ageClass,
                        terminal_combination=terminal) 
  str(data.dyad)


#Add sample IDs in each dyad as separate columns into the data frame and exclude self-comparisons

# extracting sample-combinations present in the matrices
  list<-expand.grid(key$Sample_name,key$Sample_name) 
  str(list)

# This created sample-to-same-sample pairs as well. Get rid of these:
  list<-list[which(list$Var1!=list$Var2),] 

# this still has both quantiles in--> add 'unique' key 
  list$key <- apply(list, 1, function(x)paste(sort(x), collapse='')) 
  list<-subset(list, !duplicated(list$key)) 
# sanity check that the sample name combinations are in the same exact order as the lower quantile value vector of the matrices
  i=73
  AitM[which(rownames(AitM)==list$Var1[i]),which(colnames(AitM)==list$Var2[i])]==ait[i]

# add the names of both samples participating in each dyad into the data frame
  data.dyad$IDA<-list$Var2
  data.dyad$IDB<-list$Var1
  head(data.dyad)
  

#Individual IDs in each dyad

  str(sample_data)
  IndivID<- sample_data[,c("sample.id", "BirdID")]
  colnames(IndivID)<- c("IDA", "BirdIDA")
  data.dyad2<- merge(data.dyad,IndivID, by="IDA" )
  colnames(IndivID)<- c("IDB", "BirdIDB")
  data.dyad3<- merge(data.dyad2,IndivID, by="IDB" )
  head(data.dyad3)
  str(data.dyad3)
  

# Make sure you have got rid of all same-sample comparisons
  data.dyadFinal<-data.dyad3[which(data.dyad3$IDA!=data.dyad3$IDB),] 
  str(data.dyadFinal) #106491
  
  hist(data.dyadFinal$Microbiome_similarity)


#label whether comparison is same bird or different bird comparison

  
  data.dyadFinal$Individual_combination<- ifelse(data.dyadFinal$BirdIDA==data.dyadFinal$BirdIDB, "Same", "Different")
  str(data.dyadFinal)
  table(data.dyadFinal$Individual_combination)

 
  saveRDS(data.dyadFinal, "data.dyadFinal.rds")


  
  
#################################
######### MODELLING #############
#################################


  data.dyad<-readRDS("data.dyadFinal.rds")
  head(data.dyad)

  str(data.dyad) #106491 pairs

  
  
##### 1) within individuals- comparisons of samples from the same age class ######
  
# eliminate all between individual comparisons
  
  data.dyad_within<-data.dyad[which(data.dyad$BirdIDA==data.dyad$BirdIDB),] 
  str(data.dyad_within)# 271 same individual comparisons
  length(unique(data.dyad_within$BirdIDA)) # 129 unique birds
  
  
# eliminate comparisons between different age groups 
  # (limit to comparisons of individuals within the same age group for age model- YY, MM, OO)
  
  data.dyad_withinAge<- subset(data.dyad_within, age_similarity==1)
  str(data.dyad_withinAge) #155 same individual comparisons
  head(data.dyad_withinAge)
  
  table(data.dyad_withinAge$ageclass_combination) # 56 YY, 58 MM, 41 OO comparisons
  data.dyad_withinAge %>%
    group_by(ageclass_combination, season_similarity) %>%
    dplyr::summarise(count=n())
  table(data.dyad_withinAge$terminal_combination)
  
  
  Old<- subset(data.dyad_withinAge, ageclass_combination=="OO")
  head(Old)
  oldIDs<- c(Old$IDA, Old$IDB)
  length(unique(oldIDs)) # 56 samples in OO comparisons
  oldBirdIDs<- c(Old$BirdIDA, Old$BirdIDB)
  
  
  
  length(unique(oldBirdIDs)) # 24 birds in OO comparisons (with self)
  mean(Old$Temporal_distance) # mean = 554.6585 days - 1.52 yrs
  range(Old$Temporal_distance) # 15 days to 1455 days (3.98 yrs)
  
  Mid<- subset(data.dyad_withinAge, ageclass_combination=="MM")
  MidIDs<- c(Mid$IDA, Mid$IDB)
  length(unique(MidIDs)) # 82 samples in MM comparisons
  MidBirdIDs<- c(Mid$BirdIDA, Mid$BirdIDB)
  length(unique(MidBirdIDs)) # 36 birds in MM comparisons (with self)
  mean(Mid$Temporal_distance)
  range(Mid$Temporal_distance)
  
  Young<- subset(data.dyad_withinAge, ageclass_combination=="YY")
  YIDs<- c(Young$IDA, Young$IDB)
  length(unique(YIDs)) # 89 samples in YY comparisons
  YBirdIDs<- c(Young$BirdIDA, Young$BirdIDB)
  length(unique(YBirdIDs)) # 41 birds in YY comparisons (with self)
  mean(Young$Temporal_distance)
  range(Young$Temporal_distance)
  
  
#scale all predictors to range between 0-1 if they are not already naturally on that scale
  
  #define scaling function:
  range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
  
  data.dyad_withinAge$Temporal_distance<- range.use(data.dyad_withinAge$Temporal_distance, min.use = 0, max.use = 1)
  
  str(data.dyad_withinAge)
  
  
#Make age combination into a factor for age model
  data.dyad_withinAge$ageclass_combination<-factor(data.dyad_withinAge$ageclass_combination, levels=c("OO","MM","YY"))
  levels(data.dyad_withinAge$ageclass_combination)
  
#Make BirdID into a factor for age model
  data.dyad_withinAge$BirdIDA<-factor(data.dyad_withinAge$BirdIDA)
  data.dyad_withinAge$BirdIDB<-factor(data.dyad_withinAge$BirdIDB)
  
  
### Preliminary plotting to get an idea of the data.
  
  hist(data.dyad_withinAge$Microbiome_similarity, breaks=20) 
  min(data.dyad_withinAge$Microbiome_similarity) #0.08
  max(data.dyad_withinAge$Microbiome_similarity) #0.665
  
  # Beta distribution most appropriate as can range between 0 and 1. 
  

  ggplot(data.dyad_withinAge,
         aes(ageclass_combination, Microbiome_similarity))+
    geom_boxplot()+
    geom_point()
  
  ggplot(data.dyad_withinAge,
         aes(ageclass_combination, Temporal_distance))+
    geom_boxplot() +
    geom_point()
  
  ggplot(data.dyad_withinAge,
         aes(Temporal_distance, Microbiome_similarity))+
    geom_point() +
    geom_smooth()
  

  
  
###  Construct the models with brms package ###
  
  # Let the intercept vary across a multimembership random effect of sample IDs in each pair

  
### MODEL 1 - stability test WITHIN individuals and age groups ###
  
  
# distance of samples taken from same individual within the same age group
# (hypothesis- does similarity decrease in old age group due to lower stability within individuals?)
  

### a) Use default priors in first run
  
  set.seed(921)
  model1_within<-brm(Microbiome_similarity ~ 1 + ageclass_combination
                     + Temporal_distance
                     + (1|mm(IDA,IDB)),  
                     data = data.dyad_withinAge, 
                     family= "Beta",
                     warmup = 1000, iter = 3000, 
                     cores = 4, chains = 4,
                     init=0)
  
  prior_summary(model1_within) # 1 divergent transition, Bulk ESS too low, missing info low
  summary(model1_within)
  

  
## model diagnostics
  
  #traceplots
  plot(model1_within) # phi and random intercept not very good
  
  # Posterior predictive checks can be performed using pp_check function
  
  pp_check(model1_within, ndraws=100) #not bad
  
  pp_check(model1_within, type = "error_hist", ndraws = 11)
  pp_check(model1_within, type = "scatter_avg", ndraws = 100)
  pp_check(model1_within, x = 'Temporal_distance', type='error_scatter_avg_vs_x')
  pp_check(model1_within, type = "stat_2d")
  pp_check(model1_within, type = "loo_pit_qq", ndraws = 100)
  
  
  
### b) set priors to regularise model and penalise extreme estimates ###

  #plot out gamma distribution
  c(prior(gamma(0.01, 0.01)),  # brms default
    prior(gamma(1, 0.1))) %>%  # alternative phi (moves away from low values- constrained to positive real values, and is concentrated primarily in the double-digit range indicating greater precision)
    parse_dist() %>% 
    
    ggplot(aes(xdist = .dist_obj, y = prior)) + 
    stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
    scale_x_continuous(expression(italic(p)(phi)), breaks = 0:4 * 25) +
    scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
    labs(title = "Phi Prior") +
    coord_cartesian(xlim = c(0, 110))
  
# set new priors #
  priors <- c(set_prior("normal(0,1)", class="b"),
              set_prior("gamma(1,0.1)", class="phi"),
              set_prior("student_t(3,0,2.5)", class="Intercept"),
              set_prior("student_t(3,0,2.5)",class="sd", group="mmIDAIDB"))
  
  
  set.seed(902)
  model1_within2<-brm(Microbiome_similarity ~ 1 + ageclass_combination
                     + Temporal_distance
                     + (1|mm(IDA,IDB)),  
                     data = data.dyad_withinAge, 
                     family= "Beta", prior=priors,
                     warmup = 2000, iter = 4000,
                     cores = 4, chains = 4,
                     init=0)
  
  prior_summary(model1_within2)
  summary(model1_within2) #much better
  

  saveRDS(model1_within2, "model1_withinAge.rds")
  
  model1_within2<-readRDS("model1_withinAge.rds")
  
  
## model diagnostics
  
  #traceplots
  plot(model1_within2)
  
  # Posterior predictive checks can be performed using pp_check function
  
  pp_check(model1_within2, ndraws=100) 
  # no major systematic discrepancies of data from what can be predicted from model

  
  pp_check(model1_within2, type = "error_hist", ndraws = 10)
  pp_check(model1_within2, type = "scatter_avg", ndraws = 100)
  pp_check(model1_within2, x = 'Temporal_distance', type='error_scatter_avg_vs_x')
  pp_check(model1_within2, type = "stat_2d")
  pp_check(model1_within2, type = "loo_pit_qq", ndraws = 100)

  
  
##### plot the results of the model ####
  

  resrand<-summary(model1_within2)$random
  resdf<-summary(model1_within2)$fixed
  resdf<-as.data.frame(resdf)
  resdf<-resdf[c("Estimate","l-95% CI","u-95% CI")]
  resdf<-resdf[2:nrow(resdf),]
  resdf$Predictor<-rownames(resdf)
  colnames(resdf)<-c("Estimate","lCI","uCI","Predictor")
  str(resdf)
  
  resdf$Predictor<-factor(resdf$Predictor, levels=c("Temporal_distance", "ageclass_combinationMM","ageclass_combinationYY"))
  
  ticks<-rev(resdf$Predictor)
  
  plot_mod1<-ggplot(resdf,aes(x=Estimate,y=Predictor))+
    geom_linerange(aes(xmin = lCI, xmax = uCI),linewidth=2.5, colour="grey45")+
    geom_point(size=4,colour="black", shape=21, fill="white")+
    theme_bw()+
    theme(legend.position='none',text = element_text(size=22))+
    scale_y_discrete(labels= c("Temporal distance", "YY", "MM"))+
    labs(x="\nEffect on microbiome similarity",y="")+
    scale_x_continuous(limits=c(-1,1), breaks=c( -1,-0.5, 0, 0.5, 1))+
    geom_vline(xintercept=0, linetype="dashed")
  
  plot_mod1
  
  

  
  
  
  
  
######## MODEL 2 - BETWEEN individual comparisons #########
  
  data.dyad<- readRDS("data.dyadFinal.rds")
  
# eliminate all self comparisons
  str(data.dyad) #106491
  
  data.dyad_between<-data.dyad[which(data.dyad$BirdIDA!=data.dyad$BirdIDB),] 
  str(data.dyad_between) #106220
  
  
##### One sample at random per individual ######
  
  
  sample_data<-readRDS("sample_data_from_phyloseq.rds")
  str(sample_data)
  
#### take one sample at random-remember to set seed so reproducible
  
  set.seed(53223)
  RandomSample<- sample_data %>%
    group_by(BirdID) %>%
    sample_n(1)
  
  RandomSample<- data.frame(RandomSample)
  
  RandomSample %>%
    group_by(Ageclass2, TerminalYear) %>%
    dplyr::summarise(count=n()) 
  
  p1<-ggplot(sample_data, aes(AgeYears, as.factor(BirdID),col=TerminalYear)) + geom_point()
  p2<-ggplot(RandomSample, aes(AgeYears, as.factor(BirdID), col=TerminalYear)) + geom_point()
  
  p1+p2
  
  
  
  # YA = 83 n, 32 y; MA 74 n, 27 y ; OA 40 n, 17 y
  
  RandomSampleIDs<- RandomSample$sample.id
  
  str(data.dyad_between) #106220
  data.dyad_between$IDB<- as.character(data.dyad_between$IDB)
  data.dyad_between$IDA<- as.character(data.dyad_between$IDA)
  
  data.dyad_between2<- subset(data.dyad_between, IDA %in% RandomSampleIDs)
  data.dyad_between3<- subset(data.dyad_between2, IDB %in% RandomSampleIDs)
  
  str(data.dyad_between3) #37128
  
  data.dyadfinal_RandomSample<- data.dyad_between3[complete.cases(data.dyad_between3), ]
  str(data.dyadfinal_RandomSample) #37128
  
  data.dyadfinal_RandomSample$IDB<- as.factor(data.dyadfinal_RandomSample$IDB)
  data.dyadfinal_RandomSample$IDA<- as.factor(data.dyadfinal_RandomSample$IDA)
  
  saveRDS(data.dyadfinal_RandomSample, "data.dyad_RandomSample.rds")
  

  
##### for model #####
  
  str(data.dyadfinal_RandomSample)
  head(data.dyadfinal_RandomSample)
  
# eliminate comparisons between age groups
  # (i.e. limit to comparisons of individuals within the same age group- YY, MM, OO)
  data.dyad_betweenAge<- subset(data.dyadfinal_RandomSample, age_similarity==1)
  str(data.dyad_betweenAge)  #13201
  head(data.dyad_betweenAge)
  
#eliminate comparisons that are more than a year apart
  data.dyad_betweenAge2<- subset(data.dyad_betweenAge, Temporal_distance<=365.25)
  str(data.dyad_betweenAge2) #4867
  table(data.dyad_betweenAge2$ageclass_combination) # 2349 YY, 1883 MM, 635 OO
  table(data.dyad_betweenAge2$season_similarity) # 2113, 2754 - more or less equal 
  data.dyad_betweenAge2 %>%
    group_by(ageclass_combination, season_similarity) %>%
    dplyr::summarise(count=n()) # similar proportions across age classes (i.e. young not more likely to be in same season)
  
  range(data.dyad_betweenAge2$Microbiome_similarity) #0.04- 0.75 similarity scores, mean is 0.36
  
  
#Make BirdIDs into factors
  data.dyad_betweenAge2$BirdIDA<-factor(data.dyad_betweenAge2$BirdIDA)
  data.dyad_betweenAge2$BirdIDB<-factor(data.dyad_betweenAge2$BirdIDB)
  
#Make ageclass combination into a factor 
  data.dyad_betweenAge2$ageclass_combination<-factor(data.dyad_betweenAge2$ageclass_combination,
                                                     levels=c("OO","MM","YY"))
  
#Make terminal combination into a factor 
  data.dyad_betweenAge2$terminal_combination<-factor(data.dyad_betweenAge2$terminal_combination,
                                                     levels=c("NN","NT","TT"))
  
  str(data.dyad_betweenAge2)

  
  data.dyadfinal<- data.dyad_betweenAge2[complete.cases(data.dyad_betweenAge2), ]
  str(data.dyadfinal) #4867
  
  data.dyadfinal %>%
    group_by(ageclass_combination, terminal_combination) %>%
    dplyr::summarise(count=n())
  
  table(data.dyadfinal$terminal_combination) # 2570 NN, 1912 NT, 385 TT (75 TT are in OO, 135 are in MM, 175 in YY)
  ggplot(data.dyadfinal, aes(terminal_combination, Microbiome_similarity)) + geom_boxplot()
  
#save dataset
  
  saveRDS(data.dyadfinal, "data.dyad_betweenAgeclass_oneSample.rds")
  
  data.dyadfinal<-readRDS("data.dyad_betweenAgeclass_oneSample.rds")

 
###### run the model #####
  
## NOTE this was done on a High-perfomance Computing Cluster at UEA using R 3.6.2 because of memory constraints
  
## R script for running dyadic model (BUT RUN ON CLUSTER): ##
  
  data.dyad_between<-readRDS("data.dyad_betweenAgeclass_oneSample.rds")
  
  
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
  
  range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
  
#scale temporal distance
  data.dyad_between$Temporal_distance<- range.use(data.dyad_between$Temporal_distance, min.use = 0, max.use = 1)
  
  
###  Construct the models with brms package ###
# Let the intercept vary across a multimembership random effect of sample IDs/individual IDs in each pair
# use normal(0,1) for b and gamma(1,0.1) priors to regularise model
  
  priors<- c(set_prior("normal(0,1)", class="b"),
             set_prior("student_t(3,0,2.5)", class="sd", group="mmIDAIDB"),
             set_prior("student_t(3,0,2.5)",class="Intercept"),
             set_prior("gamma(1,0.1)", class="phi"))
  
  
  set.seed(123)
  
  model1<-brm(bf(Microbiome_similarity~ 1 + ageclass_combination +
                   sex_similarity + terminal_combination +
                   Temporal_distance + (1|mm(IDA,IDB)), decomp= "QR"),
              data = data.dyad_between,
              prior=priors,
              family= "Beta",
              warmup = 2000, iter = 8000, thin=2,
              cores = 4, chains = 4,
              init=0)
  
  saveRDS(model1, "modeldist_Between.rds")
  
  pdf("modeldist_Between.pdf")
  pairs(model1)
  plot(model1)
  as_draws_df(model1) %>%
    mcmc_acf(pars = vars(b_Intercept:sd_mmIDAIDB__Intercept))
  pp_check(model1, ndraws=100)
  pp_check(model1, type = "loo_pit_qq", ndraws=100)
  dev.off()
  
  
  output<- readRDS("modeldist_Between.rds")
  summary(output)
  prior_summary(output)
  pp_check(output, type="stat_2d")
  pairs(output)
  
  
  cov2cor(vcov(output))
  
  
##### plot the results of the model ######
  
  resdf2<-summary(output)$fixed
  resrand2<- summary(output)$rand
  resdf2<-as.data.frame(resdf2)
  resdf2<-resdf2[c("Estimate","l-95% CI","u-95% CI")]
  resdf2<-resdf2[2:nrow(resdf2),]
  resdf2$Predictor<-rownames(resdf2)
  colnames(resdf2)<-c("Estimate","lCI","uCI","Predictor")
  str(resdf2)
  
  
  resdf2$Predictor<-factor(resdf2$Predictor, levels=c("Temporal_distance","sex_similarity","terminal_combinationNT", "terminal_combinationTT", "ageclass_combinationYY", "ageclass_combinationMM"))
  
  ticks<-rev(resdf2$Predictor)
  
  plot_mod2<-ggplot(resdf2,aes(x=Estimate,y=Predictor))+
    geom_linerange(aes(xmin = lCI, xmax = uCI),linewidth=2.5, colour="grey45")+
    geom_point(size=2.5,colour="black", shape=21, fill="white")+
    theme_bw()+
    theme(legend.position='none',text = element_text(size=22))+
    labs(x="\nEffect on microbiome similarity",y="")+
    scale_y_discrete(labels= c("Temporal distance","Sex similarity", "NT", "TT", "YY", "MM"))+
    scale_x_continuous(limits=c(-0.5,0.5), breaks=c(-0.5, 0, 0.5))+
    geom_vline(xintercept=0, linetype="dashed")
  
  plot_mod2
  
  
  
  
  
  
  
  
  
  


  

    