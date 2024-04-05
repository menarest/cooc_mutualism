#' ----
#' title: Lepidoptera data processing - co-occurrences
#' authors: Esteban Menares & Hugo Saiz
#' date: 16.11.2021
#' ----

#' Project: TrophCost - co-occurrence analysis
#' Processing step No: 2/4


#' **Overall workflow data processing for Lepidoptera data:** 
#' 1. Data processing 1 - trophic-links
#' 2. Data processing 2 - co-occurrences 
#' 3. Data processing 3 - species correlation matrices
#' 4. Data processing 4 - flower visitors (data validation)
  

#' **General Aim**: add co-occurrences into the species-level dataset 
#' of imago trophic-links data for each region. We don't have 
#' co-occurrence data for larva. 


#' **cols to add...**
#' co_occurrence:       Co-occurrence value (obs - exp)/(obs + exp) = association value
#' exp_rii:             Mean of expected RII index across all randomisations                           using a pairwise null model implemented through the                            "curveball algorithm"
#' sites_plant:         Number of sites where the plant species occur 
#' sites_lepi:          Number of sites where the butterfly species occur
#' obs_cooccur:         Observed number of sites where both species co-occur
#' exp_cooccur:         Expected number of sites having both species
#' prob_cooccur:        Probability that both species occur at a site 
#' p_lt:                Probability of co-occurrence at a frequency lower 
#'                      than the observed frequency
#' p_gt:                Probability of co-occurrence at a frequency greater 
#'                      than the observed frequency
#' std_effects:         Standardized effect sizes = (obs - exp)/No. of sites


# Set up --------------------------------------------------------------

## ---- Requirements
library(tidyverse)  # ggplots and data wrangling 
library(data.table) # loading and wrangling data
library(cooccur)    # species co-occurrence
# devtools::install_github("GotelliLab/EcoSimR")
library(EcoSimR) # for co-occurrence null models
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)



# 1. Reading data ---------------------------------------------------------

# long format dataset with all possible pairs of plants-lepidoptera per region
imago_SCH <- fread(file = "data/processed/troph_int_imago_SCH.csv") 
imago_ALB <- fread(file = "data/processed/troph_int_imago_ALB.csv") 

# site_spp matrix
cooccur_all <- fread(file = "data/raw/occur.csv") 

# NOTES: 
# - the datasets imago_SCH and imago_ALB were outputted in ~/troph-cost/scripts/wrangling/pairs_1_trophic_links.R. 
# - The dataset cooccur_all was cleaned up and prepared in section 1) in ~/troph-cost/scripts/wrangling/data_tidying.R  




# 2. Data tidying ---------------------------------------------------------

## ---- remove NAs
cooccur_all <- na.omit(cooccur_all)



## ---- remove sites of region HAI to match trophic-link matrix
dim(cooccur_all)
cooccur_all <- cooccur_all[!grepl("HEG", EPID), ] 
dim(cooccur_all)




## ---- save unique plants and lepi names from imago per region 

# imago_SCH
names_lepi_imago_SCH <- unique(imago_SCH$lepidoptera_sp) # names imago spp 
names_plant_imago_SCH <- unique(imago_SCH$plant_sp)      # names plants spp

# imago_ALB
names_lepi_imago_ALB <- unique(imago_ALB$lepidoptera_sp) # names imago spp 
names_plant_imago_ALB <- unique(imago_ALB$plant_sp)      # names plants spp




## ---- subset per region plants and butterflies from trophic-link matrices

# imago_SCH
cooccur_imago_SCH <- 
  cooccur_all %>%
  select(EPID, all_of(c(names_lepi_imago_SCH, names_plant_imago_SCH))) 
str(cooccur_imago_SCH)

# imago_ALB
cooccur_imago_ALB <- 
  cooccur_all %>%
  select(EPID, all_of(c(names_lepi_imago_ALB, names_plant_imago_ALB))) 
str(cooccur_imago_ALB)
# remove unnecessary objects from environment 
rm(cooccur_all) 




## ---- make EPID as row.names so that they are not treated as data
# imago_SCH
class(cooccur_imago_SCH)
cooccur_imago_SCH <- as.data.frame(cooccur_imago_SCH) # set up as data frame, cooccur needs data frame
class(cooccur_imago_SCH)
row.names(cooccur_imago_SCH) <- cooccur_imago_SCH$EPID # make EPIDs as row names
cooccur_imago_SCH$EPID <- NULL # remove EPID column

# imago_ALB
class(cooccur_imago_ALB)
cooccur_imago_ALB <- as.data.frame(cooccur_imago_ALB) # set up as data frame, cooccur needs data frame
class(cooccur_imago_ALB)
row.names(cooccur_imago_ALB) <- cooccur_imago_ALB$EPID # make EPIDs as row names
cooccur_imago_ALB$EPID <- NULL # remove EPID column


#' Check that the dimensions of cooccur_imago datasets are equal to the 
#' sum of unique names of plants + lepidoptera in imago datasets
length(c(names_plant_imago_SCH, 
         names_lepi_imago_SCH)) == dim(cooccur_imago_SCH)[2]
length(c(names_plant_imago_ALB, 
         names_lepi_imago_ALB)) == dim(cooccur_imago_ALB)[2]

# remove unnecessary objects from environment 
rm(names_lepi_imago_SCH, names_lepi_imago_ALB,
   names_plant_imago_SCH, names_plant_imago_ALB)




## ---- make P/A matrices 
# imago_SCH
cooccur_imago_SCH_PA <- cooccur_imago_SCH
cooccur_imago_SCH_PA[cooccur_imago_SCH_PA >0] <- 1

# imago_ALB
cooccur_imago_ALB_PA <- cooccur_imago_ALB
cooccur_imago_ALB_PA[cooccur_imago_ALB_PA >0] <- 1

# remove unnecessary objects from environment
rm(cooccur_imago_SCH, cooccur_imago_ALB) 




## ---- make numeric as integers
cooccur_imago_SCH_PA <- cooccur_imago_SCH_PA %>% 
  mutate(across(where(is.numeric), as.integer))

cooccur_imago_ALB_PA <- cooccur_imago_ALB_PA %>% 
  mutate(across(where(is.numeric), as.integer))

str(cooccur_imago_SCH_PA)
str(cooccur_imago_ALB_PA)




## ---- subset co-occurrences per region
# imago ALB
dim(cooccur_imago_ALB_PA)
cooccur_imago_ALB_PA <- cooccur_imago_ALB_PA[
  grepl("AEG", row.names(cooccur_imago_ALB_PA)), ]
dim(cooccur_imago_ALB_PA)

# imago SCH
dim(cooccur_imago_SCH_PA)
cooccur_imago_SCH_PA <- cooccur_imago_SCH_PA[
  grepl("SEG", row.names(cooccur_imago_SCH_PA)), ]
dim(cooccur_imago_SCH_PA)




## ---- remove cols of species without occurrences (colSum = 0)
# Region ALB
cooccur_imago_ALB_PA_filtered <- 
  cooccur_imago_ALB_PA %>% 
  select(
    where(
      ~ sum(.x) !=0 ))

# Region SCH
cooccur_imago_SCH_PA_filtered <- 
  cooccur_imago_SCH_PA %>% 
  select(
    where(
      ~ sum(.x) !=0))



## ---- remove species without occurrences from trophic-links matrices

#' first save the names of butterflies and plants species in cooccur_imago_ALB_PA and .SCH
cols_cooccur_imago_ALB_PA_filtered <-
  cooccur_imago_ALB_PA_filtered %>%
  select(everything()) %>%
  names()

cols_cooccur_imago_SCH_PA_filtered <-
  cooccur_imago_SCH_PA_filtered %>%
  select(everything()) %>%
  names()

#' then use names to filter species from imago that don't occur
imago_ALB_filtered <- imago_ALB[
  plant_sp %in% cols_cooccur_imago_ALB_PA_filtered,] # common plants
imago_ALB_filtered <- imago_ALB_filtered[
  lepidoptera_sp %in% cols_cooccur_imago_ALB_PA_filtered,] # common lepidoptera

imago_SCH_filtered <- imago_SCH[
  plant_sp %in% cols_cooccur_imago_SCH_PA_filtered,] # common plants
imago_SCH_filtered <- imago_SCH_filtered[
  lepidoptera_sp %in% cols_cooccur_imago_SCH_PA_filtered,] # common lepidoptera




## ---- Check that the dimensions of cooccurrence and imago datasets are equal
#' Use the sum of unique names of plants + lepidoptera of imago datasets

#' make a vector with unique names
names_lepi_imago_ALB <- unique(imago_ALB_filtered$lepidoptera_sp) 
names_plant_imago_ALB <- unique(imago_ALB_filtered$plant_sp)

names_lepi_imago_SCH <- unique(imago_SCH_filtered$lepidoptera_sp) 
names_plant_imago_SCH <- unique(imago_SCH_filtered$plant_sp)

length(c(names_lepi_imago_ALB, 
         names_plant_imago_ALB)) == dim(cooccur_imago_ALB_PA_filtered)[2]

length(c(names_lepi_imago_SCH, 
         names_plant_imago_SCH)) == dim(cooccur_imago_SCH_PA_filtered)[2]



## ---- remove unnecessary objects from environment 
rm(cooccur_imago_SCH_PA, cooccur_imago_ALB_PA, 
   imago_SCH, imago_ALB)





# 3. Run co-occurrence analysis ------------------------------------------

#' **NOTE**: the analysis was done in two ways. 
#' First, manually calculating the observed and expected co-occurrence and the association metric (RII) between pairs, 

#' Second, using the cooccur package which calculates co-occurrences using the hypergeometric distribution to derive a p-values and obtained the standardized size effects for each pair (Veech 2013).


# 3.1 RII calculation ------------------------------------------------------

#' **Manual calculation workflow**
#' define a function to: 
#' Step 1: calculate observed co-occurrences
#' Step 2: calculate expected co-occurrences
#' Step 3: calculate association matrix (RII index = obs-exp/obs+exp)
#' Step 4: define a function to randomize the observe matrix and calculate p-values



#' **STEP 1: calculate observed co-occurrences**
#' obs = observed number of sites having both species
#' Do a particular pair of species co-occur at a specific site?



#' **STEP 2: calculate expected co-occurrences**
#' exp = expected number of sites having both species
#' Obtained through basic probability theory (Bowers & Brown, 1982, Veech, 2006, 2013; Araujo et al., 2011) 
#' exp = (N1/N) x (N2/N) x N  = (N1 * N2) / N
#' N1 = number of sites that have species 1
#' N2 = number of sites that have species 2
#' N = total number of sites. 
#' **NOTE** that N1/N and N2/N are the marginal probabilities of occurrence of species 1 and 2, respectively, also known as their incidence rates.



#' **STEP 3: calculate association matrix (co-occurrences)**
#' "The most straightforward way to measure co-occurrence between two species is by the observed number of times that the two species co-occur relative to the expected number of times that the two species do not co-occur (Sanderson, 2000; Sfenthourakis et al., 2004, 2006; Veech, 2006, 2013; Pitta et al., 2012). Also known as the "natural metric" (Sfenthourakis et al., 2004, 2006; Pitta et al., 2012)." (Veech 2013) 

#' To calculate co-occurrences, we adapted the RII index from Armas et al. (2004) "RII is a comparative index to measure the Relative Interaction Intensity of plant interactions. RII has values ranging from -1 to 1, is symmetrical around zero and is negative for competition and positive for facilitation. RII can be calculated for any type of net interaction (from competitive exclusion to symbiosis) and in the range of 0.6 ≤ alpha ≤ 1.6 it has the lowest variability. RII has a “logarithm” behavior, being asymptotic to (y) lines 1 and -1 as a sigmoid function. More-over, the index can be scaled up and used to measure multispecific interactions at the community level". (Armas et al., 2004)

#' We use the formula: 
#' assoc = (obs - exp) / (obs + exp) 

#' 0 < assoc ≤ +1: Positive co-occurrence or aggregated, two species co-occur together at more locations than would be expected if each were randomly distributed relative to the other species.

#' -1 ≤ assoc < 0: Negative co-occurrence or segregated, two species co-occur at fewer locations than expected if each were randomly distributed relative to the other species. 

#' assoc = 0: Random co-occurrence, the two species are distributed randomly (or independently) of one another

#' Thus, a pair of species could potentially have a random, positive, or negative association – and the type of association can indirectly indicate (or provide supporting evidence for) the ecological process(es) or factor(s) leading to the particular pattern of co-occurrence.



## 3.1.1 ALB -------------------------------------------------------------

#' **NOTE**: we work on a site (rows) by species (cols) matrix.
#' Therefore, after analysis: N = (No. of sites) = (No. of rows)


## ---- create a vector with occurrences per species
occur_imago_ALB <- colSums(cooccur_imago_ALB_PA_filtered) 
#View(occur_imago_ALB)


## ---- observed co-occurrences

# create an empty array

obs_ALB <- array(0,c(dim(cooccur_imago_ALB_PA_filtered)[2],
                     dim(cooccur_imago_ALB_PA_filtered)[2]),
                 dimnames = list(cols_cooccur_imago_ALB_PA_filtered,
                                 cols_cooccur_imago_ALB_PA_filtered)) # create empty array

# create a for loop 

for (i in 1:length(occur_imago_ALB)){ 
  for (j in 1:length(occur_imago_ALB)){
    obs_ALB[i,j] <- 
      sum(cooccur_imago_ALB_PA_filtered[,i] + 
            cooccur_imago_ALB_PA_filtered[,j] == 2) # calculate obs co-occurrences 
  }
}
#View(obs_ALB)

# set diagonal of a matrix to zero

diag(obs_ALB) <- 0 # we don't need the obs occurrences of the same species
#View(obs_ALB)




## ---- expected co-occurrences 

# create empty array 

exp_ALB <- array(0,c(dim(cooccur_imago_ALB_PA_filtered)[2],
                      dim(cooccur_imago_ALB_PA_filtered)[2]),
                 dimnames = list(cols_cooccur_imago_ALB_PA_filtered,
                                 cols_cooccur_imago_ALB_PA_filtered)) # create empty array

# create for loop

for (i in 1:length(occur_imago_ALB)){
  for (j in 1:length(occur_imago_ALB)){
    exp_ALB[i,j] <- #' calculate exp = (N1 * N2) / N
      occur_imago_ALB[i] * occur_imago_ALB[j] / 
      dim(cooccur_imago_ALB_PA_filtered)[1] # 1 = dim of rows
  }
}
#View(exp_ALB)

# set diagonal of a matrix to zero

diag(exp_ALB) <- 0
#View(exp_ALB)




## ---- associations matrix

# create empty array

assoc_ALB <- array(0,c(dim(cooccur_imago_ALB_PA_filtered)[2],
                       dim(cooccur_imago_ALB_PA_filtered)[2]),
                   dimnames = list(cols_cooccur_imago_ALB_PA_filtered,
                                   cols_cooccur_imago_ALB_PA_filtered)) # create empty array

# create for loop

for (i in 1:length(occur_imago_ALB)){
  for (j in 1:length(occur_imago_ALB)){
    assoc_ALB[i,j] <- 
      (obs_ALB[i,j] - exp_ALB[i,j]) / # (Obs - Exp)/(Obs + Exp)
      (obs_ALB[i,j] + exp_ALB[i,j])   # values -1 to +1
  }
}
#View(assoc_ALB)

# change diagonal to zero to avoid the error "missing value where TRUE/FALSE needed" in the next loop
diag(assoc_ALB) <- NA #TODO change to 0 if using threshold




## 3.1.2 SCH -----------------------------------------------------------

## ---- create a vector with occurrences per species
occur_imago_SCH <- colSums(cooccur_imago_SCH_PA_filtered) 
#View(occur_imago_SCH)



## ---- observed co-occurrences
obs_SCH <- array(0,c(dim(cooccur_imago_SCH_PA_filtered)[2],
                     dim(cooccur_imago_SCH_PA_filtered)[2]),
                 dimnames = list(cols_cooccur_imago_SCH_PA_filtered,
                                 cols_cooccur_imago_SCH_PA_filtered)) # create empty array

for (i in 1:length(occur_imago_SCH)){ 
  for (j in 1:length(occur_imago_SCH)){
    obs_SCH[i,j] <- 
      sum(cooccur_imago_SCH_PA_filtered[,i] + cooccur_imago_SCH_PA_filtered[,j] == 2) # calculate obs co-occurrences 
  }
}
#View(obs_SCH)

# set diagonal of a matrix to zero

diag(obs_SCH) <- 0 # we don't need the co-occurrence of the same species
#View(obs_SCH)





## ---- expected co-occurrences 
exp_SCH <- array(0,c(dim(cooccur_imago_SCH_PA_filtered)[2],
                     dim(cooccur_imago_SCH_PA_filtered)[2]),
                 dimnames = list(cols_cooccur_imago_SCH_PA_filtered,
                                 cols_cooccur_imago_SCH_PA_filtered)) # create empty array

for (i in 1:length(occur_imago_SCH)){
  for (j in 1:length(occur_imago_SCH)){
    exp_SCH[i,j] <- # calculate exp = (N1 * N2) / N
      occur_imago_SCH[i] * occur_imago_SCH[j] / 
      dim(cooccur_imago_SCH_PA_filtered)[1] # 1 = dim of rows
  }
}
#View(exp_SCH)

# set diagonal of a matrix to zero
diag(exp_SCH) <- 0
#View(exp_SCH)




## ---- associations matrix
assoc_SCH <- array(0,c(dim(cooccur_imago_SCH_PA_filtered)[2],
                      dim(cooccur_imago_SCH_PA_filtered)[2]),
                  dimnames = list(cols_cooccur_imago_SCH_PA_filtered,
                                  cols_cooccur_imago_SCH_PA_filtered)) # create empty array

for (i in 1:length(occur_imago_SCH)){
  for (j in 1:length(occur_imago_SCH)){
    assoc_SCH[i,j] <- 
      (obs_SCH[i,j] - exp_SCH[i,j]) / # (Obs - Exp)/(Obs + Exp)
      (obs_SCH[i,j] + exp_SCH[i,j])   # values -1 to +1
  }
}
#View(assoc_SCH)

# change diagonal to zero to avoid the error "missing value where TRUE/FALSE needed" in the next loop
diag(assoc_SCH) <- NA #TODO change to 0 if using threshold

# clean global environment

rm(exp_ALB, exp_SCH, obs_ALB, obs_SCH, occur_imago_ALB, occur_imago_SCH)


## Save association matrices ---

# write_rds(assoc_ALB, 'data/processed/assoc_RII_ALB.rds')
# write_rds(assoc_SCH, 'data/processed/assoc_RII_SCH.rds')






# 3.2 Cooccur package --------------------

#' **NOTE**: 
#' we work on a species by site (rows, cols) matrix (transposed 
#' matrix). Therefore, after analysis: N = (No. of sites) = (No. of rows)

#' After running the analysis manually, we also run the analysis with the 
#' cooccur package and check the values of obs and exp obtained manually.
#' We also get the values for: 
#' sp1_inc      = Number of sites (or samples) that have species 1 (sites_plant)
#' sp2_inc      = Number of sites that have species 2 (sites_lepi)
#' prob_cooccur = Probability that both species occur at a site. 
#' p_lt         = Probability that the two species would co-occur at a 
#' frequency less than the observed number of co-occurrence sites if the 
#' two species were distributed randomly (independently) of one another.
#' p_gt         = Probability of co-occurrence at a frequency greater than
#' the observed frequency

#' "p_lt and p_gt can be obtained analytically under the condition where a 
#' species probability of occurrence at each site is equal to its observed 
#' frequency among all the sites." (Veech 2014). 

#' "if p_lt < 0.05, then the species pair co-occurs at a frequency lower 
#' than we would expect to find by chance. If p_gt < 0.05, the pair 
#' co-occurs at a rate higher than we would expect to find by chance. 
#' Since we’ve stored only significant interactions, either p_lt or p_gt
#' will be less than 0.05 for each row."
#' **NOTE END**

 


## 3.2.1 ALB --------------------------------------------------------------

class(cooccur_imago_ALB_PA_filtered)

# turn into matrix
matrix_cooccur_ALB <- as.matrix(cooccur_imago_ALB_PA_filtered)
class(matrix_cooccur_ALB)

# transpose matrix to have a spp by site matrix
matrix_cooccur_ALB_transposed <- t(matrix_cooccur_ALB)

# calculate co-occurrences with cooccur package
cooccur_ALB <- 
  cooccur(
    mat = matrix_cooccur_ALB_transposed,
    type = "spp_site", # rows = spp, cols = sites
    thresh = FALSE, # remove species pairs that are expected to have less than 1 co-occurrence
    spp_names = TRUE, # if species names are in the rows or colums
    true_rand_classifier = 0.1, # default, truly random associations that deviate from their exp by < 10%
    only_effects = FALSE, # if TRUE, shows only effect sizes
    eff_standard = TRUE, # if FALSE, effects = obs - exp. If TRUE, standardized effects
    eff_matrix = FALSE) # if TRUE, effect sizes returned in a distance matrix
class(cooccur_ALB)

#' return an analysis-wide count of the number of species combinations classified as 
#' positive, negative, or random. 
summary(cooccur_ALB)

# access and save the results from cooccur object 
cooccur_imago_ALB_package <- cooccur_ALB$results

# Calculate standardized effect sizes = (obs - exp)/No. of sites
cooccur_ALB_effects <- effect.sizes(mod = cooccur_ALB, 
                                    standardized = TRUE,
                                    matrix = FALSE)






## 3.2.2 SCH ----------------------------------------------------------

class(cooccur_imago_SCH_PA_filtered)

# turn into matrix
matrix_cooccur_SCH <- as.matrix(cooccur_imago_SCH_PA_filtered)
class(matrix_cooccur_SCH)

# transpose matrix
matrix_cooccur_SCH_transposed <- t(matrix_cooccur_SCH)

# calculate co-occurrences with cooccur package
cooccur_SCH <- 
  cooccur(
    mat = matrix_cooccur_SCH_transposed,
    type = "spp_site", # rows = spp, cols = sites
    thresh = FALSE, # remove species pairs that are expected to have less than 1 co-occurrence
    spp_names = TRUE, # if species names are in the rows or colums
    true_rand_classifier = 0.1, # default, truly random associations that deviate from their exp by < 10%
    only_effects = FALSE, # if TRUE, shows only effect sizes
    eff_standard = TRUE, # if FALSE, effects = obs - exp. If TRUE, standardized effects
    eff_matrix = FALSE) # if TRUE, effect sizes returned in a distance matrix
class(cooccur_SCH)

#' return an analysis-wide count of the number of species combinations classified as 
#' positive, negative, or random. 
summary(cooccur_SCH)


# access and save the results from cooccur object 
cooccur_imago_SCH_package <- cooccur_SCH$results

# Calculate standardized effect sizes = (obs - exp)/No. of sites
cooccur_SCH_effects <- effect.sizes(mod = cooccur_SCH, 
                                    standardized = TRUE,
                                    matrix = FALSE)







# 4. Wrangling data  ----------------------------------------------------

##  4.1 ALB -------------------------------------------------------------

## ---- Wrangling manually calculated co-occurrences

# Plant-Lepi co-occurrences: filter out all Lepi-Lepi and Plant-Plant pairs 
# to select only Plant-Lepi pairs

assoc_ALB_filtered <-
  assoc_ALB[rownames(assoc_ALB) %in% names_lepi_imago_ALB,
            colnames(assoc_ALB) %in% names_plant_imago_ALB]
#View(assoc_ALB_filtered) 

# Melt matrix into long format 
assoc_ALB_filtered <- reshape2::melt(assoc_ALB_filtered, 
                                     na.rm = TRUE, 
                                     value.name = "co_occurrence")

# Rename variables as in imago_ALB_filtered
names(assoc_ALB_filtered) <- c("lepidoptera_sp",
                               "plant_sp",
                               "co_occurrence")

# Check that the dimensions are the same
dim(assoc_ALB_filtered)[1] == dim(imago_ALB_filtered)[1]

#' **NOTE**: this we replaced with the part of 5. Update, we removed after
# Join datasets
imago_ALB_filtered <- dplyr::left_join(imago_ALB_filtered,
                                       assoc_ALB_filtered,
                                       by = c("plant_sp",
                                              "lepidoptera_sp"))




## ---- Wrangling cooccur package calculated co-occurrences

# make a copy of results to work on the self-co-occurrences [update]

cooccur_imago_ALB_package_new <- cooccur_imago_ALB_package

# Harmonize names of both datasets
colnames(cooccur_imago_ALB_package)[which(
  names(cooccur_imago_ALB_package) == "sp1_name")] <- "lepidoptera_sp"

colnames(cooccur_imago_ALB_package)[which(
  names(cooccur_imago_ALB_package) == "sp2_name")] <- "plant_sp"

colnames(cooccur_imago_ALB_package)[which(
  names(cooccur_imago_ALB_package) == "sp1_inc")] <- "sites_lepi"

colnames(cooccur_imago_ALB_package)[which(
  names(cooccur_imago_ALB_package) == "sp2_inc")] <- "sites_plant"

# Matching of both co-occurrence datasets by plant and butterfly species
imago_ALB <- merge(imago_ALB_filtered[, c("plant_sp", 
                                          "lepidoptera_sp", 
                                          "troph_pair",
                                          "troph_interaction",
                                          "co_occurrence")], 
                   cooccur_imago_ALB_package[, c("plant_sp", 
                                                 "lepidoptera_sp",
                                                 "sites_plant",
                                                 "sites_lepi",
                                                 "obs_cooccur",
                                                 "exp_cooccur",
                                                 "prob_cooccur",
                                                 "p_lt",
                                                 "p_gt")], 
                   by = c("plant_sp",
                          "lepidoptera_sp"), 
                   all.x = TRUE)



## ---- Adding standardized size effects 

# make a copy of results to work on the self-co-occurrences [update]

cooccur_ALB_effects_new <- cooccur_ALB_effects

# Harmonize names of both datasets

colnames(cooccur_ALB_effects)[which(
  names(cooccur_ALB_effects) == "sp1")] <- "lepidoptera_sp"

colnames(cooccur_ALB_effects)[which(
  names(cooccur_ALB_effects) == "sp2")] <- "plant_sp"

colnames(cooccur_ALB_effects)[which(
  names(cooccur_ALB_effects) == "effects")] <- "std_effects"

# Merge both datasets

imago_ALB <- merge(imago_ALB[, c("plant_sp", 
                                 "lepidoptera_sp", 
                                 "troph_pair",
                                 "troph_interaction",
                                 "co_occurrence", 
                                 "sites_plant",
                                 "sites_lepi",
                                 "obs_cooccur",
                                 "exp_cooccur",
                                 "prob_cooccur",
                                 "p_lt",
                                 "p_gt")], 
                   cooccur_ALB_effects[, c("plant_sp", 
                                           "lepidoptera_sp",
                                           "std_effects")], 
                   by = c("plant_sp",
                          "lepidoptera_sp"), 
                   all.x = TRUE)



##  4.2 SCH -------------------------------------------------------------

## ---- Manually calculated RII

# Filter out all Lepi-Lepi and Plant-Plant pairs i.e. select only Plant-Lepi pairs 

assoc_SCH_filtered <- 
  assoc_SCH[rownames(assoc_SCH) %in% names_lepi_imago_SCH,
            colnames(assoc_SCH) %in% names_plant_imago_SCH] 
#View(assoc_SCH_filtered) 

# Melt matrix into long format 
assoc_SCH_filtered <- reshape2::melt(assoc_SCH_filtered, 
                                     na.rm = TRUE, 
                                     value.name = "co_occurrence")

# Rename variables as in imago_ALB_filtered
names(assoc_SCH_filtered) <- c("lepidoptera_sp", 
                               "plant_sp",
                               "co_occurrence")

# Check that the dimensions are the same
dim(assoc_SCH_filtered)[1] == dim(imago_SCH_filtered)[1]

#' **NOTE**: this we replaced with the part of 5. Update, we removed after
# Join datasets
imago_SCH_filtered <- dplyr::left_join(imago_SCH_filtered,
                                       assoc_SCH_filtered,
                                       by = c("plant_sp",
                                              "lepidoptera_sp"))



## ---- Wrangling cooccur package calculated co-occurrences

# make a copy of results to work on the self-co-occurrences [update]

cooccur_imago_SCH_package_new <- cooccur_imago_SCH_package

# Harmonize names of both datasets
colnames(cooccur_imago_SCH_package)[which(
  names(cooccur_imago_SCH_package) == "sp1_name")] <- "lepidoptera_sp"

colnames(cooccur_imago_SCH_package)[which(
  names(cooccur_imago_SCH_package) == "sp2_name")] <- "plant_sp"

colnames(cooccur_imago_SCH_package)[which(
  names(cooccur_imago_SCH_package) == "sp1_inc")] <- "sites_lepi"

colnames(cooccur_imago_SCH_package)[which(
  names(cooccur_imago_SCH_package) == "sp2_inc")] <- "sites_plant"

# Matching of both co-occurrence datasets by plant and butterfly species
imago_SCH <- merge(imago_SCH_filtered[, c("plant_sp", 
                                          "lepidoptera_sp", 
                                          "troph_pair",
                                          "troph_interaction",
                                          "co_occurrence")], 
                   cooccur_imago_SCH_package[, c("plant_sp", 
                                                 "lepidoptera_sp",
                                                 "sites_plant",
                                                 "sites_lepi",
                                                 "obs_cooccur",
                                                 "exp_cooccur",
                                                 "prob_cooccur",
                                                 "p_lt",
                                                 "p_gt")], 
                   by = c("plant_sp",
                          "lepidoptera_sp"), 
                   all.x = TRUE)



## ---- Adding standardized size effects 

# make a copy of results to work on the self-co-occurrences [update]

cooccur_SCH_effects_new <- cooccur_SCH_effects

# Harmonize names of both datasets
colnames(cooccur_SCH_effects)[which(
  names(cooccur_SCH_effects) == "sp1")] <- "lepidoptera_sp"

colnames(cooccur_SCH_effects)[which(
  names(cooccur_SCH_effects) == "sp2")] <- "plant_sp"

colnames(cooccur_SCH_effects)[which(
  names(cooccur_SCH_effects) == "effects")] <- "std_effects"

# Merge both datasets
imago_SCH <- merge(imago_SCH[, c("plant_sp",  
                                 "lepidoptera_sp", 
                                 "troph_pair",
                                 "troph_interaction",
                                 "co_occurrence", 
                                 "sites_plant",
                                 "sites_lepi",
                                 "obs_cooccur",
                                 "exp_cooccur",
                                 "prob_cooccur",
                                 "p_lt",
                                 "p_gt")], 
                   cooccur_SCH_effects[, c("plant_sp", 
                                           "lepidoptera_sp",
                                           "std_effects")], 
                   by = c("plant_sp",
                          "lepidoptera_sp"), 
                   all.x = TRUE)







# 5. [UPDATE]: Add pairwise null models for RII significance -----------------

#' **NOTE**: 
#' get_rii: takes a site (rows) by species (cols) matrix and calculates RII. 
#' EcoSimR::sim9 and sim9_single: generate randomised matrices with species  in rows and sites in cols, we need to transpose the matrix!

## 5.1 Define function to calculate Armas's RII  ------------------------------

# this function works with a sites by species matrix. It also takes a species by site matrix by it transpose it first and then run the analysis. 

get_rii <- function(m,
                    matrix_type = "spp_site",
                    filter_names = FALSE,
                    names_taxa1,
                    names_taxa2,
                    supress_sp_names = FALSE) {
  
  # defining matrix type
  if (matrix_type == "site_spp") {
    mat <- m
  }
  if (matrix_type == "spp_site") {
    mat <- t(m)
  }
  
  # get total occurrences per species
  occur <- colSums(mat)
  
  # get total number of species
  n_species = ncol(mat)
  
  # calculate observed co-occurrences
  obs <- matrix(
    0,
    nrow = ncol(mat),
    ncol = ncol(mat),
    dimnames = list(colnames(mat), colnames(mat))
  )
  
  for (i in 1:length(occur)) {
    for (j in 1:length(occur)) {
      obs[i, j] <- sum(mat[, i] + mat[, j] == 2)
    }
  }
  
  # set diagonal to zero
  diag(obs) <- 0
  
  # calculate expected co-occurrences
  exp <- matrix(
    0,
    nrow = ncol(mat),
    ncol = ncol(mat),
    dimnames = list(colnames(mat), colnames(mat))
  )
  
  for (i in 1:length(occur)) {
    for (j in 1:length(occur)) {
      exp[i, j] <- occur[i] * occur[j] / nrow(mat)
    }
  }
  
  # set diagonal to zero
  diag(exp) <- 0
  
  # calculate associations i.e. RII index per species
  assoc <- matrix(
    0,
    nrow = ncol(mat),
    ncol = ncol(mat),
    dimnames = list(colnames(mat), colnames(mat))
  )
  
  for (i in 1:length(occur)) {
    for (j in 1:length(occur)) {
      assoc[i, j] <- (obs[i, j] - exp[i, j]) / (obs[i, j] + exp[i, j])
    }
  }
  
  # set diagonal to NA
  diag(assoc) <- NA
  
  if (filter_names) {
    if (missing(names_taxa1) || 
        missing(names_taxa2) || 
        length(names_taxa1) == 0 || 
        length(names_taxa2) == 0){
      stop("When filter_names is TRUE, names_taxa1 and names_taxa2 must be provided.")
    }
    
    # Check if vectors of taxa names are valid subsets
    if (!all(names_taxa1 %in% colnames(mat))) {
      stop(
        "names_taxa1 contain invalid species names or check if the the provided matrix matches the specified matrix_type."
      )
    }
    
    if (!all(names_taxa2 %in% colnames(mat))) {
      stop(
        "names_taxa2 contain invalid species names or check if the the provided matrix matches the specified matrix_type."
      )
    }
    
    # Filter plants and butterflies pairs only
    filtered_assoc <- assoc[rownames(assoc) %in% names_taxa1,
                            colnames(assoc) %in% names_taxa2]
    
    # Melt matrix into long format
    melted_assoc <- reshape2::melt(
      filtered_assoc,
      na.rm = TRUE,
      varnames = c("taxa1", "taxa2"),
      value.name = "obs_rii"
    )
    
    if (supress_sp_names) {
      melted_assoc <- dplyr::arrange(melted_assoc, taxa1, taxa2)
      melted_assoc <- dplyr::select(melted_assoc, -c(taxa1, taxa2))
      return(melted_assoc$obs_rii)
      
    } else {
      return(dplyr::arrange(melted_assoc, taxa1, taxa2))
    }
    
  } else {
    filtered_assoc <- assoc
    
    # Melt matrix into long format
    melted_assoc <- reshape2::melt(
      filtered_assoc,
      na.rm = TRUE,
      varnames = c("taxa1", "taxa2"),
      value.name = "obs_rii"
    )
    
    # Convert taxa1 and taxa2 to character
    melted_assoc$taxa1 <- as.character(melted_assoc$taxa1)
    melted_assoc$taxa2 <- as.character(melted_assoc$taxa2)
    
    # Extract values only from upper or lower triangle
    upper_triangle <-
      melted_assoc[melted_assoc$taxa1 <= melted_assoc$taxa2,]
    unique_pairs <-
      unique(upper_triangle[, c("taxa1", "taxa2", "obs_rii")])
    
    if (supress_sp_names) {
      unique_pairs <- dplyr::arrange(unique_pairs, taxa1, taxa2)
      unique_pairs <- dplyr::select(unique_pairs, -c(taxa1, taxa2))
      return(unique_pairs$obs_rii)
      
    } else {
      return(dplyr::arrange(unique_pairs, taxa1, taxa2))
    }
  }
}

# TODO: this function could be faster and more efficient if we remove the loops and do more vectorized operations... 


## 5.2 Run some test for the function ---------------------------------------

# save vectors with lepi and plant names

# lepi <- cooccur_imago_ALB_PA_filtered %>% 
#   as_tibble() %>% 
#   select(Aphantopus_hyperantus:Vanessa_cardui) %>% 
#   names()
# 
# plant <- cooccur_imago_ALB_PA_filtered %>% 
#   as_tibble() %>% 
#   select(Achillea_millefolium:last_col()) %>% 
#   names()
# 
# # test without removing names of species, with filtering 
# 
# test_rii_WITH_names <- get_rii(
#   m = cooccur_imago_ALB_PA_filtered,
#   matrix_type = "site_spp",
#   filter_names = TRUE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = FALSE
# )
# 
# # test with removing names i.e. returning only rii values as vector
# 
# test_rii_NO_names <- get_rii(
#   m = cooccur_imago_ALB_PA_filtered,
#   matrix_type = "site_spp",
#   filter_names = TRUE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = TRUE
# )
# 
# # test removing names of species, without filtering 
# 
# test_rii_no_names_no_filtering <- get_rii(
#   m = cooccur_imago_ALB_PA_filtered,
#   matrix_type = "site_spp",
#   filter_names = FALSE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = TRUE
# )
# 
# # test with removing name i.e. returning only rii values for compatibility with null models of EcoSim9
# 
# test_rii_no_filtering <- get_rii(
#   m = cooccur_imago_ALB_PA_filtered,
#   matrix_type = "site_spp",
#   filter_names = FALSE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = FALSE
# )
# 
# 
# ## ---- test with sample matrix from Lavender et al. 2019
# 
# test_mat <- read.table("resources_troph-cost/co-occurrence/Lavender_et_al_2019/DataS1_sampleMatrix.txt", sep = ",")
# 
# test_mat <- as.matrix(test_mat)
# rownames(test_mat) <- paste0("sp", 1:nrow(test_mat))
# colnames(test_mat) <- paste0("site", 1:ncol(test_mat))
# test_mat
# 
# # save a vector with species names 
# 
# plant <- rownames(test_mat[1:15,])
# lepi <- rownames(test_mat[16:20,])
# 
# # run analysis 
# 
# rii_test_mat <- get_rii(
#   m = test_mat,
#   matrix_type = "spp_site",
#   filter_names = FALSE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = FALSE
# ) 
# 
# class(rii_test_mat)
# 
# # supressing species names 
# 
# rii_test_mat_supress_sp_names <- get_rii(
#   m = test_mat,
#   matrix_type = "spp_site",
#   filter_names = FALSE,
#   names_taxa1 = plant,
#   names_taxa2 = lepi,
#   supress_sp_names = TRUE
# ) 
# 
# class(rii_test_mat_supress_sp_names)
# 
# ## ---- compare veech's pairwise approach results with Lavender's sample matrix
# 
# veech_test <- cooccur::cooccur(test_mat, 
#                                type = "spp_site", 
#                                thresh = FALSE, 
#                                spp_names = TRUE,
#                                prob = "hyper")
# 
# veech_test$results %>% view()

# OK

# END OF TESTS FOR RII FUNCTION





## 5.3 Calculate pairwise null models for RII --------------------------------

# We need to be able to compare this metric (RII) to something, to see if we obtain significance values. We will randomized the matrix then calculate the mean RII for all the generated null models and: 

# 1. calculate obs co-occurrences for each species pair in the observed matrix
# 2. calculate expected co-occurrences for each species pair in the observed matrix
# 3. calculate RII index using the obs and exp co-occurrences in the observed matrix. 
# 4. create random matrices using the sim9 “curveball algorithm” using the rule of thumb of >4N (more than four times the total number of 1s in the network) for the number of swaps (Neal et al. 2023)
# 4. repeat steps 1-3 for each pair on each randomised matrix. 
# 5. calculate the two-tailed probabilities to obtain the p-values for negative and positive co-occurrences for each pair.

# Based on: 
# https://rdrr.io/cran/EcoSimR/man/sim9.html
# https://rdrr.io/cran/EcoSimR/man/sim9_single.html 
# https://github.com/GotelliLab/EcoSimR/issues/54 
# https://github.com/emhart/Misc_Func/blob/master/NullModelsWks.R#L23


null_model_rii <- function(m, n_reps = 1, matrix_type = "spp_site", filter_names = FALSE, names_taxa1, names_taxa2) {
  
  n_species <- nrow(m)
  
  if (filter_names) {
    obs <- get_rii(m, matrix_type = "spp_site", filter_names = TRUE, names_taxa1, names_taxa2)
    null <- matrix(NA, nrow = length(names_taxa1)*length(names_taxa2), ncol = n_reps)
    
    for (i in 1:n_reps) {
      cat("Progress:", i, "/", n_reps, "\n")
      flush.console()
      
      mat <- EcoSimR::cooc_null_model(
        speciesData = m,
        algo = "sim9",
        nReps = 1,
        algoOpts = list(burn_in = 5 * sum(m == 1)),
        saveSeed = FALSE,
        suppressProg = TRUE
      )$Randomized.Data
      
      rii_values <- get_rii(mat,
                            matrix_type = "spp_site",
                            filter_names = TRUE,
                            names_taxa1, 
                            names_taxa2,
                            supress_sp_names = TRUE)
      
      if (length(rii_values) != length(names_taxa1)*length(names_taxa2)) {
        stop("The length of returned values in get_rii does not match the expected number of pairs.")
      }
      null[, i] <- rii_values
    }
    
    results <- cbind(obs, data.frame(
      exp_rii = rowMeans(null),
      p_lower = sapply(obs$obs_rii, function(x) mean(x >= null)),
      p_upper = sapply(obs$obs_rii, function(x) mean(x <= null)))) |>
      as.data.frame()
    
  } else {
    obs <- get_rii(m, filter_names = FALSE)
    null <- matrix(NA, nrow = choose(n_species, 2), ncol = n_reps)
    
    for (i in 1:n_reps) {
      cat("Progress:", i, "/", n_reps, "\n")
      flush.console()
      
      mat <- EcoSimR::cooc_null_model(
        speciesData = m,
        algo = "sim9",
        nReps = 1,
        algoOpts = list(burn_in = 5 * sum(m == 1)),
        saveSeed = FALSE,
        suppressProg = TRUE
      )$Randomized.Data
      
      rii_values <- get_rii(mat, 
                            matrix_type = "spp_site",
                            supress_sp_names = TRUE)
      
      if (length(rii_values) != choose(n_species, 2)) {
        stop("The length of returned values in get_rii does not match the expected number of pairs.")
      }
      null[, i] <- rii_values
    }
    
    results <- cbind(obs, data.frame(
      exp_rii = rowMeans(null),
      p_lower = sapply(obs$obs_rii, function(x) mean(x >= null)),
      p_upper = sapply(obs$obs_rii, function(x) mean(x <= null)))) |>
      as.data.frame()
  }
  
  rm(null)
  gc()
  
  return(results)
}

# check results of null models with the RII calculated with the get_rii function and manually with loops (section 3.1)

# Wrangle data 

# lepi <- cooccur_imago_ALB_PA_filtered %>%
#   as_tibble() %>%
#   select(Aphantopus_hyperantus:Vanessa_cardui) %>%
#   names()
# 
# plant <- cooccur_imago_ALB_PA_filtered %>%
#   as_tibble() %>%
#   select(Achillea_millefolium:last_col()) %>%
#   names()
# 
# rii <- get_rii(t(cooccur_imago_ALB_PA_filtered),
#                matrix_type = "spp_site",
#                filter_names = TRUE,
#                names_taxa1 = plant,
#                names_taxa2 = lepi,
#                supress_sp_names = FALSE)
 
# res_rii_alb <- null_model_rii(t(cooccur_imago_ALB_PA_filtered), 
#                               n_reps = 10,
#                               filter_names = FALSE,
#                               names_taxa1 = plant, 
#                               names_taxa2 = lepi)
# 
# res_rii_alb_filt <- null_model_rii(t(cooccur_imago_ALB_PA_filtered), 
#                                    n_reps = 10,
#                                    filter_names = TRUE,
#                                    names_taxa1 = plant, 
#                                    names_taxa2 = lepi)

# stopifnot(imago_ALB$co_occurrence == res_rii_alb_filt$obs_rii) 
# stopifnot(imago_ALB$co_occurrence == rii$obs_rii)
# stopifnot(rii$obs_rii == res_rii_alb_filt$obs_rii) 

# test_taxa1_plant <- null_model_rii(m, n_reps = 1, filter_names = TRUE, names_taxa1 = plant, names_taxa2 = lepi)
# 
# imago_ALB_filtered %>% 
#   left_join(
#     test_taxa1_plant,
#     by = c('plant_sp' = 'taxa1', 'lepidoptera_sp' = 'taxa2')
#   ) %>% 
#   left_join(
#     rii,
#     by = c('plant_sp' = 'taxa1', 'lepidoptera_sp' = 'taxa2')
#   ) %>% 
#   left_join(assoc_ALB_filtered, 
#             by = c("plant_sp", 
#                    "lepidoptera_sp")
#   ) %>% 
#   as_tibble() %>% view()



## 5.4 Run analysis  -----------------------------------------------------

## ---- get names of lepi and plants in ALB

lepi <- cooccur_imago_ALB_PA_filtered %>%
  as_tibble() %>%
  select(Aphantopus_hyperantus:Vanessa_cardui) %>%
  names()

plant <- cooccur_imago_ALB_PA_filtered %>%
  as_tibble() %>%
  select(Achillea_millefolium:last_col()) %>%
  names()

# Run RII and pairwise null models calculation 

res_rii_alb_filt <- null_model_rii(t(cooccur_imago_ALB_PA_filtered), 
                              n_reps = 1000,
                              filter_names = TRUE,
                              names_taxa1 = plant, 
                              names_taxa2 = lepi)



## ---- get names of lepi and plants in SCH

lepi <- cooccur_imago_SCH_PA_filtered %>%
  as_tibble() %>%
  select(Anthocharis_cardamines:Thymelicus_sylvestris) %>%
  names()

plant <- cooccur_imago_SCH_PA_filtered %>%
  as_tibble() %>%
  select(Achillea_millefolium:last_col()) %>%
  names()

# Run RII and pairwise null models calculation 

res_rii_sch_filt <- null_model_rii(t(cooccur_imago_SCH_PA_filtered), 
                                   n_reps = 1000,
                                   filter_names = TRUE,
                                   names_taxa1 = plant, 
                                   names_taxa2 = lepi)

## 5.5 Join results  -----------------------------------------------------

# ALB

imago_ALB <- imago_ALB %>% 
  as_tibble() %>% 
  left_join(
    res_rii_alb_filt %>% 
      as_tibble(),
    by = c('plant_sp' = 'taxa1', 'lepidoptera_sp' = 'taxa2')
  ) %>%
  select(-co_occurrence,
         co_occurrence = obs_rii) %>% 
  relocate(co_occurrence, exp_rii, p_lower, p_upper,
           .after = sites_lepi) 

# SCH

imago_SCH <- imago_SCH %>% 
  as_tibble() %>% 
  left_join(
    res_rii_sch_filt %>% 
      as_tibble(),
    by = c('plant_sp' = 'taxa1', 'lepidoptera_sp' = 'taxa2')
  ) %>%
  select(-co_occurrence,
         co_occurrence = obs_rii) %>% 
  relocate(co_occurrence, exp_rii, p_lower, p_upper,
           .after = sites_lepi)



# clean global environment

rm(
  matrix_cooccur_ALB,
  matrix_cooccur_ALB_transposed,
  matrix_cooccur_SCH,
  matrix_cooccur_SCH_transposed,
  cols_cooccur_imago_ALB_PA_filtered,
  cols_cooccur_imago_SCH_PA_filtered,
  imago_ALB_filtered,
  imago_SCH_filtered,
  assoc_ALB_filtered,
  assoc_SCH_filtered,
  cooccur_imago_ALB_package,
  cooccur_imago_SCH_package,
  cooccur_ALB_effects,
  cooccur_SCH_effects,
  rii
)






# 6. [Update] Plant-Plant and Lepi-Lepi co-occurrences -------------------

##  6.1 ALB ---------------------------------------------------------------

## ---- Wrangling manually calculated co-occurrences

# Plant-PLant co-occurrences: filter out all Lepi-Lepi and Plant-Plant pairs 
# to select only Plant-Lepi pairs 

assoc_ALB_lepi <- 
  assoc_ALB[rownames(assoc_ALB) %in% names_lepi_imago_ALB,
            colnames(assoc_ALB) %in% names_lepi_imago_ALB]

assoc_ALB_plant <- 
  assoc_ALB[rownames(assoc_ALB) %in% names_plant_imago_ALB,
            colnames(assoc_ALB) %in% names_plant_imago_ALB]

# Melt matrices into long format 

assoc_ALB_lepi_new <- reshape2::melt(assoc_ALB_lepi, na.rm = TRUE, 
                                     value.name = "co_occurrence")
assoc_ALB_plant_new <- reshape2::melt(assoc_ALB_plant, na.rm = TRUE, 
                                     value.name = "co_occurrence")

# Rename variables 

names(assoc_ALB_lepi_new) <- c("sp1_name", 
                               "sp2_name",
                               "co_occurrence")

names(assoc_ALB_plant_new) <- c("sp1_name", 
                                "sp2_name",
                                "co_occurrence")

# Matching of lepi-lepi co-occurrence datasets

cooccur_lepi_ALB <-
  assoc_ALB_lepi_new %>% 
  as_tibble() %>% 
  inner_join(
    cooccur_imago_ALB_package_new %>% 
      as_tibble() %>% 
      select(-c("sp1", "sp2")),
    by =  c("sp1_name", "sp2_name")) %>% 
  # Add standardized size effects 
  inner_join(
    cooccur_ALB_effects_new %>% 
      as_tibble(),
    by =  c("sp1_name" = "sp1", 
            "sp2_name" = "sp2")) %>% 
  # add a new region col and a lepi-lepi to specify the type of co-occurrence
  mutate(region = "ALB",
         type = "lepi_lepi") %>% 
  # reorder cols
  relocate(region,
           type, 
           sp1_name, 
           sp2_name, 
           co_occurrence:last_col()) %>% 
  arrange(sp1_name) 

# Matching of plant-plant co-occurrence datasets

cooccur_plant_ALB <-
  assoc_ALB_plant_new %>% 
  as_tibble() %>% 
  inner_join(
    cooccur_imago_ALB_package_new %>% 
      as_tibble() %>% 
      select(-c("sp1", "sp2")),
    by =  c("sp1_name", "sp2_name")) %>% 
  # Add standardized size effects 
  inner_join(
    cooccur_ALB_effects_new %>% 
      as_tibble(),
    by =  c("sp1_name" = "sp1", 
            "sp2_name" = "sp2")) %>% 
  # add a new region col and a plant-plant to specify the type of co-occurrence
  mutate(region = "ALB",
         type = "plant_plant") %>% 
  # reorder cols
  relocate(region,
           type, 
           sp1_name, 
           sp2_name, 
           co_occurrence:last_col()) %>% 
  arrange(sp1_name) 


# clean global environment 

rm(assoc_ALB, 
   assoc_ALB_lepi, 
   assoc_ALB_lepi_new, 
   assoc_ALB_plant,
   assoc_ALB_plant_new, 
   cooccur_ALB, 
   cooccur_ALB_effects_new,
   cooccur_imago_ALB_package_new)





##  6.1 SCH -------------------------------------------------------------

## ---- Wrangling manually calculated co-occurrences

# Plant-PLant co-occurrences: filter out all Lepi-Lepi and Plant-Plant pairs 
# to select only Plant-Lepi pairs 

assoc_SCH_lepi <- 
  assoc_SCH[rownames(assoc_SCH) %in% names_lepi_imago_SCH,
            colnames(assoc_SCH) %in% names_lepi_imago_SCH]

assoc_SCH_plant <- 
  assoc_SCH[rownames(assoc_SCH) %in% names_plant_imago_SCH,
            colnames(assoc_SCH) %in% names_plant_imago_SCH]

# Melt matrices into long format 

assoc_SCH_lepi_new <- reshape2::melt(assoc_SCH_lepi, na.rm = TRUE, 
                                     value.name = "co_occurrence")
assoc_SCH_plant_new <- reshape2::melt(assoc_SCH_plant, na.rm = TRUE, 
                                      value.name = "co_occurrence")

# Rename variables 

names(assoc_SCH_lepi_new) <- c("sp1_name", 
                               "sp2_name",
                               "co_occurrence")

names(assoc_SCH_plant_new) <- c("sp1_name", 
                                "sp2_name",
                                "co_occurrence")

# Matching of lepi-lepi co-occurrence datasets

cooccur_lepi_SCH <-
  assoc_SCH_lepi_new %>% 
  as_tibble() %>% 
  inner_join(
    cooccur_imago_SCH_package_new %>% 
      as_tibble() %>% 
      select(-c("sp1", "sp2")),
    by =  c("sp1_name", "sp2_name")) %>% 
  # Add standardized size effects 
  inner_join(
    cooccur_SCH_effects_new %>% 
      as_tibble(),
    by =  c("sp1_name" = "sp1", 
            "sp2_name" = "sp2")) %>% 
  # add a new region col and a lepi-lepi to specify the type of co-occurrence
  mutate(region = "SCH",
         type = "lepi_lepi") %>% 
  # reorder cols
  relocate(region,
           type, 
           sp1_name, 
           sp2_name, 
           co_occurrence:last_col()) %>% 
  arrange(sp1_name) 

# Matching of plant-plant co-occurrence datasets

cooccur_plant_SCH <-
  assoc_SCH_plant_new %>% 
  as_tibble() %>% 
  inner_join(
    cooccur_imago_SCH_package_new %>% 
      as_tibble() %>% 
      select(-c("sp1", "sp2")),
    by =  c("sp1_name", "sp2_name")) %>% 
  # Add standardized size effects 
  inner_join(
    cooccur_SCH_effects_new %>% 
      as_tibble(),
    by =  c("sp1_name" = "sp1", 
            "sp2_name" = "sp2")) %>% 
  # add a new region col and a plant-plant to specify the type of co-occurrence
  mutate(region = "SCH",
         type = "plant_plant") %>% 
  # reorder cols
  relocate(region,
           type, 
           sp1_name, 
           sp2_name, 
           co_occurrence:last_col()) %>% 
  arrange(sp1_name) 

# clean global environment 

rm(assoc_SCH, 
   assoc_SCH_lepi, 
   assoc_SCH_lepi_new, 
   assoc_SCH_plant,
   assoc_SCH_plant_new, 
   cooccur_SCH,
   cooccur_SCH_effects_new,
   cooccur_imago_SCH_package_new)

# put all cooccur tables together into one

self_cooccur <- 
  bind_rows(cooccur_plant_ALB,
            cooccur_lepi_ALB,
            cooccur_plant_SCH,
            cooccur_lepi_SCH)





# 7. Save data ----------------------------------------------------------

# write_csv(cooccur_imago_ALB_PA_filtered %>%
#             rownames_to_column(var = "EPID"),
#           file = "data/processed/occur_imago_ALB_PA.csv",
#           na = "NA")
# 
# write_csv(cooccur_imago_SCH_PA_filtered %>%
#             rownames_to_column(var = "EPID"),
#           file = "data/processed/occur_imago_SCH_PA.csv",
#           na = "NA")
# 
# write_csv(imago_ALB, file = "data/raw/imago_ALB.csv", na = "NA")
# write_csv(imago_SCH, file = "data/raw/imago_SCH.csv", na = "NA")
# write_csv(self_cooccur, file = "data/raw/self_cooccur.csv")
