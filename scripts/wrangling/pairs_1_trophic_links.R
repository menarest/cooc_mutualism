#' ----
#' title: Lepidoptera data processing - trophic-links
#' author: Esteban Menares
#' date: 11.10.2021
#' ---

#' Project: TrophCost - co-occurrence analysis
#' Processing step No: 1/4
  

#' **Overall workflow data processing for Lepidoptera data**: 
#' 1. Data processing 1 - trophic-links
#' 2. Data processing 2 - co-occurrences 
#' 3. Data processing 3 - species correlation matrices
#' 4. Data processing 4 - flower visitors (data validation)
  
  
  
#' **General Aim**: to make a species-level dataset for larva and 
#' imago data for both regions combined (overlap of both regions)
#' and decide if to analyse either combined or each region per
#' separate

#' **Columns** (for each dataset, one dataset per region):
#' plant_sp:            Plant species name
#' lepidoptera_sp:      Lepidoptera species name
#' troph_pair:          Plant-Lepidoptera Pair (abbreviated unique name)
#' troph_interaction:   Trophic-link strength from literature



## ---- Requirements
library(tidyverse) # ggplots and data wrangling 
library(data.table) # loading and wrangling data
library(cowplot) # formatting ggplots 
library(Hmisc) # descriptive statistics 




#### 1. Reading data ####

# NOTE: these datasets were cleaned up and prepared in section 3) in ~/troph-cost/scripts/wrangling/data_tidying.R and can be found in BExIS with the numbers: 

# Imago ALB: 31734
# Larva ALB: 31735
# Imago SCH: 31736
# Larva SCH: 31737

# These datasets were saved as long versions. To run this script, you would need to convert from long to wide format. Run this line when loading each dataset to turn it into wide format again: 
# %>% pivot_wider(names_from = "lepidoptera_species", values_from = "int_strength")

larva_SCH <- fread(file = "resources/trophic-interactions/Larva_BB.csv")
dim(larva_SCH)
str(larva_SCH)

larva_ALB <- fread(file = "resources/trophic-interactions/Larva_BW.csv")
dim(larva_ALB)
str(larva_ALB)

imago_SCH <- fread(file = "resources/trophic-interactions/Imago_BB.csv")
dim(imago_SCH)
str(imago_SCH)

imago_ALB <- fread(file = "resources/trophic-interactions/Imago_BW.csv")
dim(imago_ALB)
str(imago_ALB)




#### 2. Data integrity ####


## ---- Classify integers columns to class numeric 

larva_SCH <-
  larva_SCH %>%
  mutate(across(Anthocharis_cardamines:last_col(),
                ~ as.numeric(.)))
str(larva_SCH)

larva_ALB <-
  larva_ALB %>%
  mutate(across(Anthocharis_cardamines:last_col(),
                ~ as.numeric(.)))
str(larva_ALB)

imago_SCH <-
  imago_SCH %>%
  mutate(across(Anthocharis_cardamines:last_col(),
                ~ as.numeric(.)))
str(imago_SCH)

imago_ALB <-
  imago_ALB %>%
  mutate(across(Anthocharis_cardamines:last_col(),
                ~ as.numeric(.)))
str(imago_ALB)


## ---- Check if butterflies and plants names are same in all 4 datasets
names(larva_SCH) == names(larva_ALB) &
  names(larva_ALB) == names(imago_SCH) &
  names(imago_SCH) == names(imago_ALB) 

larva_SCH$plant_genus == larva_ALB$plant_genus &
  larva_ALB$plant_genus == imago_SCH$plant_genus &
  imago_SCH$plant_genus == imago_ALB$plant_genus

larva_SCH$plant_species == larva_ALB$plant_species &
  larva_ALB$plant_species == imago_SCH$plant_species &
  imago_SCH$plant_species == imago_ALB$plant_species




#### 3. Data wrangling ####



## ---- Remove unwanted cols
larva_SCH <- select(larva_SCH,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
larva_ALB <- select(larva_ALB,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
imago_SCH <- select(imago_SCH,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
imago_ALB <- select(imago_ALB,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 




## ---- remove rows with NA
larva_SCH_filter <- larva_SCH[!is.na(Anthocharis_cardamines),]
larva_ALB_filter <- larva_ALB[!is.na(Anthocharis_cardamines),]
imago_SCH_filter <- imago_SCH[!is.na(Anthocharis_cardamines),]
imago_ALB_filter <- imago_ALB[!is.na(Anthocharis_cardamines),]




## ---- save the names of removed cols (butterflies) from each datasets
rm_cols_larva_SCH <- larva_SCH_filter %>% 
  select_if(~any(is.na(.))) %>% 
  names()

rm_cols_larva_ALB <- larva_ALB_filter %>% 
  select_if(~any(is.na(.))) %>% 
  names()

rm_cols_imago_SCH <- imago_SCH_filter %>% 
  select_if(~any(is.na(.))) %>% 
  names()

rm_cols_imago_ALB <- imago_ALB_filter %>% 
  select_if(~any(is.na(.))) %>% 
  names()




## ---- remove cols with NAs
larva_SCH_filter <- larva_SCH_filter %>% 
  select_if(~any(!is.na(.)))

larva_ALB_filter <- larva_ALB_filter %>% 
  select_if(~any(!is.na(.)))

imago_SCH_filter <- imago_SCH_filter %>% 
  select_if(~any(!is.na(.)))

imago_ALB_filter <- imago_ALB_filter %>% 
  select_if(~any(!is.na(.)))




## ---- Overlap cols from both regions for larva and imago datasets
# Filter both plants and butterflies that are in both regions to 
# get a matrix with same amount of plant_sp (rows) and lepidoptera_sp 
# (cols) between both regions SCH - ALB for larva and imago independently

# check if cols removed from rm_cols_larva_ALB are in rm_cols_larva_SCH
rm_cols_larva_ALB %in% rm_cols_larva_SCH
# rm_cols_imago_ALB = 0, hence no cols with NA in imago_ALB, 
# therefore, we use rm_cols_imago_SCH to remove from imago_ALB all 9 
# butterflies removed from imago_SCH to get tow datasets of equal cols dim


# combine removed cols (larva butterfly) into a vector with unique names
rm_cols_larva <- c(rm_cols_larva_SCH, rm_cols_larva_ALB) %>% 
  sort() %>% 
  unique()

# since all 15 butterflies from combined vector rm_cols_larva were 
# already removed from larva_SCH, we only remove them from larva_ALB



## ---- intersect cols from larva_SCH and larva_ALB (cols in both)
larva_ALB_filter <- dplyr::select(larva_ALB_filter, -any_of(rm_cols_larva))
imago_ALB_filter <- dplyr::select(imago_ALB_filter, -any_of(rm_cols_imago_SCH))

# check col names are the same
stopifnot(names(larva_ALB_filter) == names(larva_SCH_filter))
stopifnot(names(imago_ALB_filter) == names(imago_SCH_filter))



## ---- intersect rows from both regions for larva and imago
# Combine plant genus and species into "plant_sp" to make unique names
larva_SCH_filter <- unite(larva_SCH_filter, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = FALSE)

larva_ALB_filter <- unite(larva_ALB_filter, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = FALSE)

imago_SCH_filter <- unite(imago_SCH_filter, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = FALSE)

imago_ALB_filter <- unite(imago_ALB_filter, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = FALSE)

# get the intersect of plant_sp from both datasets
common_plants_larva <- intersect(larva_SCH_filter$plant_sp, 
                                 larva_ALB_filter$plant_sp)  
length(common_plants_larva)

common_plants_imago <- intersect(imago_SCH_filter$plant_sp, 
                                 imago_ALB_filter$plant_sp)  
length(common_plants_imago)

# remove non-common rows from both datasets

larva_SCH_filter <- larva_SCH_filter[larva_SCH_filter$plant_sp %in% common_plants_larva,] # give you common rows
larva_ALB_filter <- larva_ALB_filter[larva_ALB_filter$plant_sp %in% common_plants_larva,] # give you common rows
stopifnot(larva_SCH_filter$plant_sp == larva_ALB_filter$plant_sp) # to check if rows are same

# STOP! check manually a few values with original datasets from excel 

imago_SCH_filter <- imago_SCH_filter[imago_SCH_filter$plant_sp %in% common_plants_imago,] # give you common rows 
imago_ALB_filter <- imago_ALB_filter[imago_ALB_filter$plant_sp %in% common_plants_imago,] # give you common rows 
stopifnot(imago_SCH_filter$plant_sp == imago_ALB_filter$plant_sp) # to check if rows are same

# STOP! check manually a few values with original datasets from excel 




## ---- Extract first character of each plant genus and 3 first characters of plant species
p_genus <- str_sub(larva_SCH_filter$plant_genus, 1, 1) # extract first character of genus
p_species <- str_sub(larva_SCH_filter$plant_species, 1, 3) # extract first 3 characters of species
No <- seq(1:length(larva_SCH_filter$plant_species)) # add a unique number to the plant names
plant <- paste("P", No, sep = "") # paste P for plant with each number
plant <- paste(plant, p_genus, p_species, sep = "") # paste all together
is.vector(plant) # to check if is vector 
plant # to see vector
rm(p_genus, p_species, No) # remove unnecessary vectors from environment 




## ---- Create short unique col names: larva
# save the old col names first
old_names_larva <- larva_SCH_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_larva, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_larva, # extract first 3 characters of species
                          "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_larva)) # add a unique number to the Lepidoptera names
lepi_larva <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_larva <-  paste(lepi_larva, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_larva)
lepi_larva # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(larva_SCH_filter[,4:46]) == old_names_larva




## ---- Create short unique col names: imago
# save the old col names first
old_names_imago <- imago_SCH_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_imago, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_imago, # extract first 3 characters of species
                        "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_imago)) # add a unique number to the Lepidoptera names
lepi_imago <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_imago <-  paste(lepi_imago, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_imago)
lepi_imago # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(imago_SCH_filter[,4:52]) == old_names_imago




## ---- Change names of butterflies for extracted ones
larva_SCH_filter %>% data.table::setnames(old = old_names_larva, new = lepi_larva)
larva_ALB_filter %>% data.table::setnames(old = old_names_larva, new = lepi_larva)
imago_SCH_filter %>% data.table::setnames(old = old_names_imago, new = lepi_imago)
imago_ALB_filter %>% data.table::setnames(old = old_names_imago, new = lepi_imago)




## ---- Add a new column plants to the data set
larva_SCH_filter$plants <- plant
larva_ALB_filter$plants <- plant
imago_SCH_filter$plants <- plant
imago_ALB_filter$plants <- plant




## ---- Reorder columns by name using select() and keep only one cell plant names 
larva_SCH_filter <- larva_SCH_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L43Vcar)

larva_ALB_filter <- larva_ALB_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L43Vcar)

imago_SCH_filter <- imago_SCH_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L49Vcar)

imago_ALB_filter <- imago_ALB_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L49Vcar)

# Remove unnecesary data from environment
rm(larva_SCH, larva_ALB, imago_SCH, imago_ALB)




# Creating unique plant-butterfly pairs: Larva

## ---- Create pairs of plant-butterfly larva_SCH_filter
# Pivot data to reorganize values into plant-butterfly pairs
LBB <- pivot_longer(larva_SCH_filter, cols = L1Acar:L43Vcar, 
             names_to ="pairs", 
             values_to = "interaction_strength")

# Combine plant-butterfly cells into one single "troph_pair" cell
LBB <- unite(LBB, plants, pairs, col = "troph_pair", sep = "-")




## ---- Create pairs of plant-butterfly larva_ALB_filter
# Pivot data to reorganize values into plant-butterfly pairs
LBW <- pivot_longer(larva_ALB_filter, cols = L1Acar:L43Vcar,
                  names_to ="pairs", 
                  values_to = "interaction_strength")

# Combine plant-butterfly cells into one single "pairs" cell
LBW <- unite(LBW, plants, pairs, col = "troph_pair", sep = "-")




## ---- check if pairs of plant-butterfly larva_SCH_filter == larva_ALB_filter
stopifnot(LBB$troph_pair == LBW$troph_pair)




## ---- Combine both datasets LBB and LBW into one final "larva"
larva <- bind_cols(LBB[,c("plant_sp",
                          "troph_pair", 
                          "interaction_strength")], 
                   LBW[,"interaction_strength"])

# rename variables uniquely
larva <- larva %>% 
  rename(troph_link_SCH = interaction_strength...3, 
         troph_link_ALB = interaction_strength...4)




## ---- Add a column "lepidoptera_sp" with the full names 
# of butterflies to match new row dimensions: 

larva <-
  larva %>%
  separate(
    col = troph_pair,
    sep = "-",
    into = c("troph_pair_P", "troph_pair_L")
  ) %>%
  left_join(bind_cols(old_names_larva, lepi_larva),
            by = c("troph_pair_L" = "...2")) %>%
  unite(troph_pair_P, troph_pair_L,
        col = "troph_pair",
        sep = "-") %>%
  rename(lepidoptera_sp = `...1`) %>%
  # reorder cols
  relocate(lepidoptera_sp, .after = plant_sp)



## ---- remove unnecessary objects from environment
rm(LBB, LBW, larva_SCH_filter, larva_ALB_filter)




#' Preparing the species-level trophic-interactions matrix:
#' In order to identify target pairs of species for conservation 
#' we need to remove all genus-level observations from the matrices



## ---- Separate col "plant_sp" into plant_genus and plant_species
larva <- separate(larva, plant_sp, sep = "_", 
         into = c("plant_genus", "plant_species"))




## ---- remove rows with genus-level plants (i.e. genus_sp)
larva <- as.data.table(larva)
larva <- larva[plant_species != "sp", ]




# replace 0.1 values for 0 to overwrite all genus-level interactions
# previously coded as 0.1 which were not removed by removing genus_sp
with(larva, table(troph_link_ALB, troph_link_SCH))
larva$troph_link_SCH <- replace(larva$troph_link_SCH, 
                              larva$troph_link_SCH == 0.1, 0)
larva$troph_link_ALB <- replace(larva$troph_link_ALB, 
                              larva$troph_link_ALB == 0.1, 0)

# check unique values from vectors to confirm that there are no 0.1
unique(larva$troph_link_SCH)
unique(larva$troph_link_ALB)



## ---- harmonize values of trophLinks in both regions
# since the level of information for trophic-links is different depending
# on the region: 
# troph_link_SCH: 0 and 1
# troph_link_ALB: 0, 0.2, 0.4, 0.6, 0.8, 1
# To be able to compare both regions, we transform values of ALB >0 to 1

larva$troph_link_ALB <- 
  replace(larva$troph_link_ALB, larva$troph_link_ALB > 0, 1)
unique(larva$troph_link_ALB) # check values



## ---- unite plant_genus and plant_species back into plant_sp 
larva <- unite(larva, plant_genus, plant_species, 
               col = "plant_sp", sep = "_", remove = TRUE)





# Creating unique plant-butterfly pairs: Imago

## ---- Create pairs of plant-butterfly imago_SCH
# Pivot data to reorganize values into plant-butterfly pairs
IBB <- pivot_longer(imago_SCH_filter, cols = L1Acar:L49Vcar, 
                    names_to ="pairs", 
                    values_to = "interaction_strength")

# Combine plant-butterfly cells into one single "troph_pairs" cell
IBB <- unite(IBB, plants, pairs, col = "troph_pair", sep = "-")




## ---- Create pairs of plant-butterfly imago_ALB
# Pivot data to reorganize values into plant-butterfly pairs
IBW <- pivot_longer(imago_ALB_filter, cols = L1Acar:L49Vcar, 
                    names_to ="pairs", 
                    values_to = "interaction_strength")

# Combine plant-butterfly cells into one single "pairs" cell
IBW <- unite(IBW, plants, pairs, col = "troph_pair", sep = "-")



## ---- check if pairs of plant-butterfly imago_SCH == imago_ALB
stopifnot(IBB$troph_pair == IBW$troph_pair)




## ---- Combine both datasets LBB and LBW into one "larva"
imago <- bind_cols(IBB[,c("plant_sp",
                          "troph_pair", 
                          "interaction_strength")], 
                   IBW[,"interaction_strength"])

# rename variables uniquely
imago <- imago %>% 
  rename(troph_link_SCH = interaction_strength...3, 
         troph_link_ALB = interaction_strength...4)


## ---- Add a column "lepidoptera_sp" with the full names of butterflies 
# to match new row dimensions: 

imago <-
  imago %>%
  separate(
    col = troph_pair,
    sep = "-",
    into = c("troph_pair_P", "troph_pair_L")
  ) %>%
  left_join(bind_cols(old_names_imago, lepi_imago),
            by = c("troph_pair_L" = "...2")) %>%
  unite(troph_pair_P, troph_pair_L,
        col = "troph_pair",
        sep = "-") %>%
  rename(lepidoptera_sp = `...1`) %>%
  # reorder cols
  relocate(lepidoptera_sp, .after = plant_sp)

  

## ---- remove unnecessary objects from environment
rm(IBB, IBW, imago_SCH_filter, imago_ALB_filter)




#' Preparing the species-level trophic-interactions matrix:
#' In order to identify target pairs of species for conservation 
#' we need to remove all genus-level observations from the matrices



## ---- Separate col "plant_sp" into plant_genus and plant_species
imago <- separate(imago, plant_sp, sep = "_", 
                  into = c("plant_genus", "plant_species"))




## ---- remove rows with genus-level plants (i.e. genus_sp)
imago <- as.data.table(imago)
imago <- imago[plant_species != "sp", ]




# replace 0.1 values for 0 to overwrite all genus-level interactions
# previously coded as 0.1 which were not removed by removing genus_sp
with(imago, table(troph_link_ALB, troph_link_SCH))

imago$troph_link_SCH <- replace(imago$troph_link_SCH, 
                              imago$troph_link_SCH == 0.1, 0)
imago$troph_link_ALB <- replace(imago$troph_link_ALB, 
                              imago$troph_link_ALB == 0.1, 0)

# check unique values from vectors to confirm that there are no 0.1
unique(imago$troph_link_SCH)
unique(imago$troph_link_ALB)



## ---- unite plant_genus and plant_species back into plant_sp 
imago <- unite(imago, plant_genus, plant_species, 
               col = "plant_sp", sep = "_", remove = TRUE)

# change troph_link_cols to numeric

larva <-
  larva %>%
  mutate(
    troph_link_SCH = as.character(troph_link_SCH),
    troph_link_ALB = as.character(troph_link_ALB)
  )

imago <-
  imago %>%
  mutate(
    troph_link_SCH = as.character(troph_link_SCH),
    troph_link_ALB = as.character(troph_link_ALB)
  )




#### 4. Data exploration and visualization ####

##### 4.2 Larva ####

# Larva interaction data was collected for possible other purposes, but not used since we did not have occurrence larva data.

## ---- Visualise data using scatterplots of pairs

p <- 
  ggplot(larva %>% 
           as_tibble() %>% 
           mutate(troph_link_SCH = factor(troph_link_SCH,
                                          levels = c("0", "1")),
                  troph_link_ALB = factor(troph_link_ALB,
                                          levels = c("0", "1"))), 
         aes(x = troph_link_SCH, 
             y = troph_link_ALB)) + 
  geom_count() + 
  scale_size_area(max_size = 15) +
  theme_cowplot(12) + # for nice formatting
  labs(x = "Trophic-links region SCH (North-East)", 
       y = "Trophic-links region ALB (South-West)",
       title ="Overlap of trophic-link pairs (plant-butterfly LARVA) per region")
p

# cowplot::ggsave2("overlap-larva-SCH-ALB.pdf",
#                  plot = ggplot2::last_plot(),
#                  path = "output/plots/",
#                  width = 7,
#                  height = 4,
#                  scale = 1)


## ---- evaluate quality of data

# make a contingency table to evaluate quality of data
df <- with(larva, table(troph_link_ALB, troph_link_SCH))

# compare 0,0; 0,1; 1,0; and 1,1 between ALB and SCH

# get the total 

sum(df)
# 6116 + 18 + 35 + 23 # 6192 OLD

# 0,0 + 1,1 ---

df[1,1] + df[2,2] # 6139
(df[1,1] + df[2,2])/sum(df)*100 # 99.14 % 


# 0,1 + 1,0 --- 

df[2,1] + df[1,2] # 53
(df[2,1] + df[1,2])/sum(df)*100 # 0.86 % 


# 0,0 ---

df[1,1] # 6116
df[1,1]/sum(df)*100 # 98.78 %


# 1,1 ---

df[2,2] # 23
df[2,2]/sum(df)*100 # 0.37 %




##### 4.2 Imago ####


## ---- Visualise data using scatterplots of pairs

imago %>% as_tibble() %>% count(troph_link_ALB)
imago %>% as_tibble() %>% count(troph_link_SCH)

p <- 
  ggplot(imago %>% 
           as_tibble() %>% 
           mutate(troph_link_SCH = 
                    factor(troph_link_SCH,
                           levels = c("0", "0.2", "0.6", "1")),
                  troph_link_ALB = 
                    factor(troph_link_ALB,
                           levels = c("0", "0.2", "0.4", "0.6", "0.8", "1"))), 
         aes(x = troph_link_SCH, 
             y = troph_link_ALB)) + 
  geom_count() + 
  scale_size_area(max_size = 15) +
  theme_cowplot(12) + # for nice formatting
  labs(x = "Trophic-links region SCH (North-East)", 
       y = "Trophic-links region ALB (South-West)",
       title ="Overlap of trophic-link pairs (plant-butterfly IMAGO) per region")
p

# cowplot::ggsave2("overlap-imago-SCH-ALB.pdf",
#                  plot = ggplot2::last_plot(),
#                  path = "output/plots/",
#                  width = 7,
#                  height = 4,
#                  scale = 1)


## ---- evaluate quality of data
with(imago, table(troph_link_ALB, troph_link_SCH))

# data has different levels, so we need to transform to P/A 

# make a copy of imago
imago_PA <- imago

# replace values to 0 and 1
imago_PA$troph_link_SCH <- replace(imago_PA$troph_link_SCH, 
                                imago_PA$troph_link_SCH > 0, 1)
imago_PA$troph_link_ALB <- replace(imago_PA$troph_link_ALB, 
                                imago_PA$troph_link_ALB > 0, 1)

# check unique values to confirm that there are only 0 and 1
unique(imago_PA$troph_link_SCH)
unique(imago_PA$troph_link_ALB)

# make contingency table
df <- with(imago_PA, table(troph_link_ALB, troph_link_SCH))
df

# compare 0,0; 0,1; 1,0; and 1,1 between ALB and SCH

# get the total 

sum(df) # 7056 

# 0,0 + 1,1 ---

df[1,1] + df[2,2] # 6087
(df[1,1] + df[2,2])/sum(df)*100 # 86.27 % 

# 0,1 + 1,0 --- 

df[2,1] + df[1,2] # 969
(df[2,1] + df[1,2])/sum(df)*100 # 13.73 % 

# 0,0 ---

df[1,1] # 5527
df[1,1]/sum(df)*100 # 78.33 %

# 1,1 ---

df[2,2] # 560
df[2,2]/sum(df)*100 # 7.94 %





#### 5. Data re-wrangling ####

#' Based on the data exploration results, in specific the low level of 
# overlap between trophic-interaction pairs at both regions, we decided to 
# analyse the data for each region per separate. This mean we will produce 
# four datasets:
#' larva_SCH
#' larva_ALB
#' imago_SCH
#' imago_ALB

#' **Columns**: 
#' - Plant species (full name)
#' - Lepidoptera species (full name)
#' - Plant-Lepidoptera Pair (abbreviated unique name)
#' - Trophic-link strength at region




## ---- First remove all objects from environment to start over again
rm(list = ls())  #will clear all objects includes hidden objects




##### 5.1 Read and tidy data ####

## ---- Reading data

larva_SCH <- fread(file = "resources_troph-cost/trophic-interactions/Larva_BB.csv")
larva_ALB <- fread(file = "resources_troph-cost/trophic-interactions/Larva_BW.csv")
imago_SCH <- fread(file = "resources_troph-cost/trophic-interactions/Imago_BB.csv")
imago_ALB <- fread(file = "resources_troph-cost/trophic-interactions/Imago_BW.csv")



## ---- Classify integers columns to class numeric 
larva_SCH <-
  larva_SCH %>%
  mutate(across(Anthocharis_cardamines:last_col(), ~ as.numeric(.)))
str(larva_SCH)

larva_ALB <-
  larva_ALB %>%
  mutate(across(Anthocharis_cardamines:last_col(), ~ as.numeric(.)))
str(larva_ALB)

imago_SCH <-
  imago_SCH %>%
  mutate(across(Anthocharis_cardamines:last_col(), ~ as.numeric(.)))
str(imago_SCH)

imago_ALB <-
  imago_ALB %>%
  mutate(across(Anthocharis_cardamines:last_col(), ~ as.numeric(.)))
str(imago_ALB)



## ---- Remove unwanted cols
larva_SCH <- select(larva_SCH,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
larva_ALB <- select(larva_ALB,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
imago_SCH <- select(imago_SCH,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 
imago_ALB <- select(imago_ALB,-c("plant_order",
                                 "plant_family",
                                 "plant_qualifiers",
                                 "plant_synonym",
                                 "plant_de_name")) 



## ---- remove rows with NA
larva_SCH_filter <- larva_SCH[!is.na(Anthocharis_cardamines),]
larva_ALB_filter <- larva_ALB[!is.na(Anthocharis_cardamines),]
imago_SCH_filter <- imago_SCH[!is.na(Anthocharis_cardamines),]
imago_ALB_filter <- imago_ALB[!is.na(Anthocharis_cardamines),]



## ---- remove cols with NAs
larva_SCH_filter <- larva_SCH_filter %>% select_if(~any(!is.na(.)))
larva_ALB_filter <- larva_ALB_filter %>% select_if(~any(!is.na(.)))
imago_SCH_filter <- imago_SCH_filter %>% select_if(~any(!is.na(.)))
imago_ALB_filter <- imago_ALB_filter %>% select_if(~any(!is.na(.)))



##### 5.2 Change names ####


## ---- Create short unique plant names (row) for each region: SCH
p_genus <- str_sub(larva_SCH_filter$plant_genus, 1, 1) # extract first character of genus
p_species <- str_sub(larva_SCH_filter$plant_species, 1, 3) # extract first 3 characters of species
No <- seq(1:length(larva_SCH_filter$plant_species)) # add a unique number to the plant names
plant_SCH <- paste("P", No, sep = "") # paste P for plant with each number
plant_SCH <- paste(plant_SCH, p_genus, p_species, sep = "") # paste all together
is.vector(plant_SCH) # to check if is vector 
plant_SCH # to see vector
rm(p_genus, p_species, No) # remove unnecessary vectors from environment 




## ---- Create short unique plant names (row) for each region: ALB
p_genus <- str_sub(larva_ALB_filter$plant_genus, 1, 1) # extract first character of genus
p_species <- str_sub(larva_ALB_filter$plant_species, 1, 3) # extract first 3 characters of species
No <- seq(1:length(larva_ALB_filter$plant_species)) # add a unique number to the plant names
plant_ALB <- paste("P", No, sep = "") # paste P for plant with each number
plant_ALB <- paste(plant_ALB, p_genus, p_species, sep = "") # paste all together
is.vector(plant_ALB) # to check if is vector 
plant_ALB # to see vector
rm(p_genus, p_species, No) # remove unnecessary vectors from environment 




## ---- Create short unique butterfly names (cols) for each dataset: larva_SCH
# save the old col names first
old_names_larva_SCH <- larva_SCH_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_larva_SCH, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_larva_SCH, # extract first 3 characters of species
                            "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_larva_SCH)) # add a unique number to the Lepidoptera names
lepi_larva_SCH <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_larva_SCH <-  paste(lepi_larva_SCH, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_larva_SCH)
lepi_larva_SCH # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(larva_SCH_filter[,3:45]) == old_names_larva_SCH




## ---- Create short unique butterfly names (cols) for each dataset: larva_ALB
# save the old col names first
old_names_larva_ALB <- 
  larva_ALB_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_larva_ALB, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_larva_ALB, # extract first 3 characters of species
                            "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_larva_ALB)) # add a unique number to the Lepidoptera names
lepi_larva_ALB <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_larva_ALB <-  paste(lepi_larva_ALB, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_larva_ALB)
lepi_larva_ALB # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(larva_ALB_filter[,3:58]) == old_names_larva_ALB




## ---- Create short unique butterfly names (cols) for each dataset: imago_SCH
# save the old col names first
old_names_imago_SCH <- imago_SCH_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_imago_SCH, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_imago_SCH, # extract first 3 characters of species
                            "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_imago_SCH)) # add a unique number to the Lepidoptera names
lepi_imago_SCH <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_imago_SCH <-  paste(lepi_imago_SCH, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_imago_SCH)
lepi_imago_SCH # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(imago_SCH_filter[,3:51]) == old_names_imago_SCH




## ---- Create short unique butterfly names (cols) for each dataset: imago_ALB
# save the old col names first
old_names_imago_ALB <- 
  imago_ALB_filter %>% 
  select("Anthocharis_cardamines":"Vanessa_cardui") %>% 
  colnames()  # save the old butterfly names first

lepi_genus <- str_sub(old_names_imago_ALB, 1, 1) # extract first character of genus
lepi_species <- str_extract(old_names_imago_ALB, # extract first 3 characters of species
                            "(?<=[:punct:])\\w{3}") 
No <- seq(1:length(old_names_imago_ALB)) # add a unique number to the Lepidoptera names
lepi_imago_ALB <- paste("L", No, sep = "") # paste L for Lepidoptera wiht each number
lepi_imago_ALB <-  paste(lepi_imago_ALB, lepi_genus, lepi_species, sep = "") # paste all together
is.vector(lepi_imago_ALB)
lepi_imago_ALB # to see vector
rm(lepi_genus, lepi_species, No) # remove unnecessary vectors from environment 

# check if the names are exactly equal
colnames(imago_ALB_filter[,3:60]) == old_names_imago_ALB




## ---- Change names of butterflies for extracted ones
larva_SCH_filter %>% data.table::setnames(old = old_names_larva_SCH, 
                                          new = lepi_larva_SCH)
larva_ALB_filter %>% data.table::setnames(old = old_names_larva_ALB, 
                                          new = lepi_larva_ALB)
imago_SCH_filter %>% data.table::setnames(old = old_names_imago_SCH, 
                                          new = lepi_imago_SCH)
imago_ALB_filter %>% data.table::setnames(old = old_names_imago_ALB, 
                                          new = lepi_imago_ALB)




## ---- Add a new column plants to the data set
larva_SCH_filter$plants <- plant_SCH
larva_ALB_filter$plants <- plant_ALB
imago_SCH_filter$plants <- plant_SCH
imago_ALB_filter$plants <- plant_ALB




## ---- Combine plant genus and species into "plant_sp" to make unique names
larva_SCH_filter <- unite(larva_SCH_filter, plant_genus, plant_species, 
                          col = "plant_sp", sep = "_", remove = FALSE)

larva_ALB_filter <- unite(larva_ALB_filter, plant_genus, plant_species, 
                          col = "plant_sp", sep = "_", remove = FALSE)

imago_SCH_filter <- unite(imago_SCH_filter, plant_genus, plant_species, 
                          col = "plant_sp", sep = "_", remove = FALSE)

imago_ALB_filter <- unite(imago_ALB_filter, plant_genus, plant_species, 
                          col = "plant_sp", sep = "_", remove = FALSE)




## ---- Reorder columns by name using select() and keep only one cell plant names 
larva_SCH_filter <- larva_SCH_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L43Vcar)

larva_ALB_filter <- larva_ALB_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L56Vcar)

imago_SCH_filter <- imago_SCH_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L49Vcar)

imago_ALB_filter <- imago_ALB_filter %>% 
  select(plant_sp, -c(plant_genus, plant_species), plants, L1Acar:L58Vcar)

# Remove unnecessary data from environment
rm(larva_SCH, larva_ALB, imago_SCH, imago_ALB)





##### 5.3 Create pairs ####

#' Create unique plant-butterfly pairs for all regions and life stages combinations


## ---- Create pairs of plant-butterfly 

# Pivot data to reorganize values into plant-butterfly pairs
larva_SCH_pivot <- pivot_longer(larva_SCH_filter, 
                                cols = L1Acar:L43Vcar, 
                                names_to ="pairs", 
                                values_to = "troph_interaction")

larva_ALB_pivot <- pivot_longer(larva_ALB_filter, 
                                cols = L1Acar:L56Vcar,
                                names_to ="pairs", 
                                values_to = "troph_interaction")

imago_SCH_pivot <- pivot_longer(imago_SCH_filter, 
                                cols = L1Acar:L49Vcar, 
                                names_to ="pairs", 
                                values_to = "troph_interaction")

imago_ALB_pivot <- pivot_longer(imago_ALB_filter, 
                                cols = L1Acar:L58Vcar, 
                                names_to ="pairs", 
                                values_to = "troph_interaction")

# Combine plant-butterfly cells into one single "troph_pair" cell
larva_SCH_pivot <- unite(larva_SCH_pivot, plants, pairs, 
                         col = "troph_pair", sep = "-")

larva_ALB_pivot <- unite(larva_ALB_pivot, plants, pairs, 
                         col = "troph_pair", sep = "-")

imago_SCH_pivot <- unite(imago_SCH_pivot, plants, pairs, 
                         col = "troph_pair", sep = "-")

imago_ALB_pivot <- unite(imago_ALB_pivot, plants, pairs, 
                         col = "troph_pair", sep = "-")




## ---- Add a column "lepidoptera_sp" 
# Use the full names of butterflies to match new row dimensions: 
# (rows of larva_SCH_pivot)/(length of old_names_larva) = repetitions
rep <- dim(larva_SCH_pivot)[1]/length(old_names_larva_SCH)
larva_SCH_pivot$lepidoptera_sp <- rep(old_names_larva_SCH, rep)

rep <- dim(larva_ALB_pivot)[1]/length(old_names_larva_ALB)
larva_ALB_pivot$lepidoptera_sp <- rep(old_names_larva_ALB, rep)

rep <- dim(imago_SCH_pivot)[1]/length(old_names_imago_SCH)
imago_SCH_pivot$lepidoptera_sp <- rep(old_names_imago_SCH, rep)

rep <- dim(imago_ALB_pivot)[1]/length(old_names_imago_ALB)
imago_ALB_pivot$lepidoptera_sp <- rep(old_names_imago_ALB, rep)




## ---- reorder by column name
larva_SCH_pivot <- larva_SCH_pivot[,c("plant_sp",
                                      "lepidoptera_sp",
                                      "troph_pair",
                                      "troph_interaction")]

larva_ALB_pivot <- larva_ALB_pivot[,c("plant_sp",
                                      "lepidoptera_sp",
                                      "troph_pair",
                                      "troph_interaction")]

imago_SCH_pivot <- imago_SCH_pivot[,c("plant_sp",
                                      "lepidoptera_sp",
                                      "troph_pair",
                                      "troph_interaction")]

imago_ALB_pivot <- imago_ALB_pivot[,c("plant_sp",
                                      "lepidoptera_sp",
                                      "troph_pair",
                                      "troph_interaction")]


## ---- remove unnecessary objects from environment
rm(larva_SCH_filter, larva_ALB_filter, imago_SCH_filter, imago_ALB_filter)





##### 5.4 Species-level matrices ####

#' In order to identify target pairs of species for conservation 
#' we need to remove all genus-level observations from the matrices



## ---- Separate col "plant_sp" into plant_genus and plant_species 
# to be able to filter plant_genus = sp
larva_SCH_pivot <- separate(larva_SCH_pivot, plant_sp, sep = "_", 
                            into = c("plant_genus", "plant_species"))

larva_ALB_pivot <- separate(larva_ALB_pivot, plant_sp, sep = "_", 
                            into = c("plant_genus", "plant_species"))

imago_SCH_pivot <- separate(imago_SCH_pivot, plant_sp, sep = "_", 
                            into = c("plant_genus", "plant_species"))

imago_ALB_pivot <- separate(imago_ALB_pivot, plant_sp, sep = "_", 
                            into = c("plant_genus", "plant_species"))



## ---- remove rows with genus-level plants (i.e. genus_sp)
larva_SCH_pivot <- as.data.table(larva_SCH_pivot)
larva_ALB_pivot <- as.data.table(larva_ALB_pivot)
imago_SCH_pivot <- as.data.table(imago_SCH_pivot)
imago_ALB_pivot <- as.data.table(imago_ALB_pivot)

larva_SCH_pivot <- larva_SCH_pivot[plant_species != "sp", ]
larva_ALB_pivot <- larva_ALB_pivot[plant_species != "sp", ]
imago_SCH_pivot <- imago_SCH_pivot[plant_species != "sp", ]
imago_ALB_pivot <- imago_ALB_pivot[plant_species != "sp", ]




## ---- Overwrite all genus-level interactions (0.1 => 0)
# to remove remaining genus-level interactions 
# which were not removed by removing genus_sp

# check if there are troph_interaction values = 0.1
table(larva_SCH_pivot$troph_interaction)
table(larva_ALB_pivot$troph_interaction)
table(imago_SCH_pivot$troph_interaction)
table(imago_ALB_pivot$troph_interaction)

# replace these values 0.1 with 0
larva_SCH_pivot$troph_interaction <- 
  replace(larva_SCH_pivot$troph_interaction, 
          larva_SCH_pivot$troph_interaction == 0.1, 0)

larva_ALB_pivot$troph_interaction <- 
  replace(larva_ALB_pivot$troph_interaction, 
          larva_ALB_pivot$troph_interaction == 0.1, 0)

imago_SCH_pivot$troph_interaction <- 
  replace(imago_SCH_pivot$troph_interaction, 
          imago_SCH_pivot$troph_interaction == 0.1, 0)

imago_ALB_pivot$troph_interaction <- 
  replace(imago_ALB_pivot$troph_interaction, 
          imago_ALB_pivot$troph_interaction == 0.1, 0)

# confirm that there are no 0.1 values
table(larva_SCH_pivot$troph_interaction)
table(larva_ALB_pivot$troph_interaction)
table(imago_SCH_pivot$troph_interaction)
table(imago_ALB_pivot$troph_interaction)




## ---- unite plant_genus and plant_species back into plant_sp 
larva_SCH_pivot <- unite(larva_SCH_pivot, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = TRUE)

larva_ALB_pivot <- unite(larva_ALB_pivot, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = TRUE)

imago_SCH_pivot <- unite(imago_SCH_pivot, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = TRUE)

imago_ALB_pivot <- unite(imago_ALB_pivot, plant_genus, plant_species, 
                         col = "plant_sp", sep = "_", remove = TRUE)



##### 5.5 [UPDATE] Tidy final table, update ids and check up names --------

# Aim: 
# put together the unfiltered imago trophic interaction tables for both regions
# (unfiltered for occurrences). Then: 

# 1) Add a region col, 
# 2) Update the troph_pair id to the newer one with leading zeros and region
# 3) Bind both tables 


# harmonize the names of plant species as in '1. Ocurrences'

# ALB

imago_ALB <-
  imago_ALB_pivot %>% 
  as_tibble() %>% 
  # add a col "region" 
  mutate(region = 'ALB') %>% 
  # reorder cols
  relocate(region, plant_sp:troph_interaction) %>% 
  # separate col troph_pair to be able to add leading zero to number in id
  extract(., troph_pair, 
          into = paste0('troph_pair', 1:5),
          # getting all characters except a dash till the first dash
          regex = '(^P)([[:digit:]]+)(.{4}-L)([[:digit:]]+)(.{4})') %>% 
  # add leading zeros to numerical groups
  mutate(
    troph_pair2 = str_pad(troph_pair2, 3, side = "left", pad = "0"),
    troph_pair4 = str_pad(troph_pair4, 3, side = "left", pad = "0")) %>% 
  # collapse cells from troph_pair cols back into one single col
  unite(.,
        troph_pair1, troph_pair2, troph_pair3, troph_pair4, troph_pair5,
        col = "troph_pair",
        sep = "") %>% 
  # separate troph_pair id into plant and butterfly id
  separate(troph_pair, 
           into = c('plant_sp_id', 'lepi_sp_id'), 
           remove = FALSE) %>%
  # correct mistaken cell entries and col names
  mutate(
    # correct wrongly written id 
    lepi_sp_id = str_replace(lepi_sp_id, 
                             '(.*)(INA$)', # assign groups to modify 
                             '\\1Iio'), # modify only group 2
    # add region into lepi id 
    lepi_sp_id = paste0("ALB",as.character(lepi_sp_id)),
    # add region into plant id 
    plant_sp_id = paste0("ALB",as.character(plant_sp_id))) %>% 
  # paste id cols back to one single pair id
  unite(., col = "troph_pair", 
        plant_sp_id, lepi_sp_id, 
        sep = "-",
        remove = FALSE) %>% 
  select(-plant_sp_id, -lepi_sp_id) %>% 
  # tidy genus-level plant species "aggr."
  filter(!str_detect(plant_sp, "aggr.")) %>% 
  mutate(
    plant_sp = str_replace(plant_sp, 
                           "Capsella_bursa-pastoris",
                           "Capsella_bursa_pastoris"),
    plant_sp = str_replace(plant_sp,
                           "Silene_flos-cuculi",
                           "Silene_flos_cuculi"))
  
# SCH

imago_SCH <-
  imago_SCH_pivot %>% 
  as_tibble() %>% 
  # add a col "region" 
  mutate(region = 'SCH') %>% 
  # reorder cols
  relocate(region, plant_sp:troph_interaction) %>% 
  # separate col troph_pair to be able to add leading zero to number in id
  extract(., troph_pair, 
          into = paste0('troph_pair', 1:5),
          # getting all characters except a dash till the first dash
          regex = '(^P)([[:digit:]]+)(.{4}-L)([[:digit:]]+)(.{4})') %>% 
  # add leading zeros to numerical groups
  mutate(
    troph_pair2 = str_pad(troph_pair2, 3, side = "left", pad = "0"),
    troph_pair4 = str_pad(troph_pair4, 3, side = "left", pad = "0")) %>% 
  # collapse cells from troph_pair cols back into one single col
  unite(.,
        troph_pair1, troph_pair2, troph_pair3, troph_pair4, troph_pair5,
        col = "troph_pair",
        sep = "") %>% 
  # separate troph_pair id into plant and butterfly id
  separate(troph_pair, 
           into = c('plant_sp_id', 'lepi_sp_id'), 
           remove = FALSE) %>%
  # correct mistaken cell entries and col names
  mutate(
    # correct wrongly written id 
    lepi_sp_id = str_replace(lepi_sp_id, 
                             '(.*)(INA$)', # assign groups to modify 
                             '\\1Iio'), # modify only group 2
    # add region into lepi id 
    lepi_sp_id = paste0("SCH",as.character(lepi_sp_id)),
    # add region into plant id 
    plant_sp_id = paste0("SCH",as.character(plant_sp_id))) %>% 
  # paste id cols back to one single pair id
  unite(., col = "troph_pair", 
        plant_sp_id, lepi_sp_id, 
        sep = "-",
        remove = FALSE) %>% 
  select(-plant_sp_id, -lepi_sp_id) %>% 
  # tidy genus-level plant species "aggr."
  filter(!str_detect(plant_sp, "aggr.")) %>% 
  mutate(
    plant_sp = str_replace(plant_sp, 
                           "Capsella_bursa-pastoris",
                           "Capsella_bursa_pastoris"),
    plant_sp = str_replace(plant_sp,
                           "Silene_flos-cuculi",
                           "Silene_flos_cuculi"))

  
# bind ALB and SCH together

troph_ints_raw <-
  bind_rows(imago_ALB, imago_SCH)




#### 6. Saving data ####

## ---- write .csv files for each region and life-stage per separate

# fwrite(larva_SCH_pivot, file = "data/processed/archive/troph_int_larva_SCH.csv")
# fwrite(larva_ALB_pivot, file = "data/processed/archive/troph_int_larva_ALB.csv")
# fwrite(imago_SCH, file = "data/processed/troph_int_imago_SCH.csv")
# fwrite(imago_ALB, file = "data/processed/troph_int_imago_ALB.csv")
# write_csv(troph_ints_raw, 'data/raw/troph_links.csv')

# clean global environment 

rm(imago_ALB, imago_SCH, troph_ints_raw)



### END OF SCRIPT ###