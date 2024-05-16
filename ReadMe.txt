This folder contains data and R script needed to reproduce the results in the paper:
Worsley S.F., Davies C.S., Lee C.Z., Mannarelli, ME. et al. Longitudinal gut microbiome dynamics in relation to age and senescence in a wild animal population. Molecular Ecology

All variable names in these data should be self-explanatory, but please don't hesitate to send an email to the corresponding author if anything is unclear.


Files included:

### Scripts ###

#SampleProcessing.R
This R script contains the code necessary to reproduce all initial pre-processing steps using outputs from QIIMEII (e.g. decontamination, sample filtering, extraction repeatability tests). The outputs from this code (Filtered ASV, taxonomy and metadata files) are used in all downstream analysis.

#AlphaDiversity.R
This R script contains the code necessary to reproduce alpha diversity analysis. It uses the filtered ASV, taxonomy and metadata files generated in SampleProcessing.R

#BetaAnalysis.R
This R script contains the code necessary to reproduce the results of the analysis investigating the association between senescence and gut microbiome beta diversity (i.e. PERMANOVA, PCA, and differential abundance tests). It uses the filtered ASV, taxonomy and metadata files generated in SampleProcessing.R

#SimilarityAnalysis.R
This R script contains the code necessary to reproduce the results of bayesian dyadic models investigating how GM personalisation and stability change with age. It uses a phyloseq object generated in BetaAnalysis.R


### Data ###

# ASV_table.csv
This file contains the raw abundances of Amplicon Sequence Variants (ASV) inferred for each faecal sample using QIIMEII. This ASV table is used to generate phyloseq objects. The paired-end sequences for each sample used to generate the table have been deposited with the ENA under the accession numbers under the study accession numbers PRJEB45408 (samples taken in 2017 and 2018) and PRJEB47095 (samples taken in 2019 and 2020) and PRJEB67634 (samples taken in 2021 and 2022)

# Taxonomy.csv
This file contains the taxonomy of ASVs in ASV_table.csv. It was generated as part of the QIIMEII pipeline. This file is used to generate phyloseq objects.

#tree.nwk
ASV rooted tree file generated as part of the QIIMEII pipeline.

# Metadata.csv
This file contains the sample metadata. Columns correspond to:

sample.id = Unique identifier for each faecal sample
BirdID = Unique identifier for each bird
TubeNumber = Tube number for each faecal sample given at the time of sampling
OriginalTubeNumberUnique = Unique tube number for sequenced faecal samples- (suffix "R" means a repeat extraction of the same sample was sequenced, suffix "b" refers to the second time a sample was sequenced to asses repeatability within and across sequencing runs)
FieldPeriodID = Unique identifier for the field period sample was taken in (164 = major17, 166= minor18, 167= major18, 171= minor19, 173= major19, 174= minor20, 175 and 176= major21, 177=minor22, 178=major22)
Season= Major or minor breeding season
SequencingRun= sequencing run for that sample
CatchID = Unique identifier for catch
OccasionID = Unique identifier for data associated with catch
BTO = Unique BTO ring number for each Bird
FieldRing = Unique colour ring combination for each Bird
CatchTerritoryID = Territory ID where bird was caught/sampled (note not necessarily breed group territory).
SampleDate = date faecal sample taken
SampleYear = year faecal sample taken (2017-2022)
FreezerDate = date samples were frozen at -80 degrees C.
TimeAt4Degrees = number of days faecal sample stored at 4 degrees C prior to freezing.
SampleType = Sample type (Control or Faecal)
SampleType4way = Sample type in further detail. Collection control; Extraction control (blank swab) or PCR negative (water control in PCR); Positive control (mock community); Faecal sample.
ExtractionDate = date of DNA extraction from faecal sample
Qubit = Concentration (ng/μl) of DNA extracted from sample
BirthDate = estimated or true birth date for each bird
SexEstimate = sex of bird (1= male, 0= female)
SampleTime = time of day bird was caught/sampled
MinutesSinceSunrise = time bird was caught as minutes since sunrise (06:00 am)
Ageclass = ageclass of bird (FL= fledgling, OFL= old fledgling, SA= Sub-adult, A= adults)
AgeDays = Age in days
AgeYears = Age in years
Eyecolour= eyecolour of bird (G= grey, LB= light brown, RB= red brown).
BreedTerritoryNumber = Territory number of that bird's breed group
BreedTerritoryID = ID for that Territory number in database
BreedGroupID = Unique ID for that bird's breed group in a particular season
Status = status of bird in population census
TQcorrected = quality of breed territory calculated using the "corrected" method in the database (note there are missing values for some seasons)
averageTQ_SameType = quality of territory (as above). For territories with missing qualities an average of the preceding and following seasons of the same type (i.e. major or minor) was taken.
TQnotes_SameType = notes display is an average of other seasons (NA if not)
FieldPeriodID_LastSeen= last field period a bird was observed (note this will be 2023 for birds that survived until the end of the study)
LastPeriodSeenDate = date of the end of the field period when bird last seen
LastBreedYear= year the bird was last observed (2023 for birds that survived)
DaysToLastSeen = difference between SampleDate and LastPeriodSeenDate columns
DiedEver= binary column indicating whether the bird died during the study (yes or no)
BirthDateEstimated = if the birth date was estimated to refine the age of new birds hatched within a year for which the exact hatch date is unknown. Estimates are based on a combination of behavioural observations and eye colour (TRUE= estimated, FALSE= exact dates known)



# FilteredASVTable.csv
A filtered ASV table generated using the SampleProcessing.R script (following decontam etc...) which is used in all downstream analysis

# FilteredTaxonomyTableBacteria.csv
A filtered taxonomy table generated using the SampleProcessing.R script (following decontam etc...) which is used in all downstream analysis

# FilteredMetadataBacteria.csv
A filtered metadata table generated using the SampleProcessing.R script (following decontam etc...) which is used in all downstream analysis

#CatchDuplicatesSorted.csv
A csv file defining which sample to keep in analyses when there are duplicate samples from the same catch (decision made based on sample source). File is used in alpha and beta scripts to filter catch duplicates.