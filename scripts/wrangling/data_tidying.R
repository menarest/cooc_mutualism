
# title: Lepidoptera data tidying for raw data
# author: Esteban Menares
# date: Aug 2021


# NOTE: many of the datasets were obtained from BExIS (BExIS (https://www.bexis.uni-jena.de/). Use the indicated ID (and version) to search the original datasets on their website. 


# Setup ---------------------------------------------------------------

library(tidyverse) # for data wrangling
library(data.table) # for data wrangling 

source('scripts/source_script.R')


# 1. OCCURRENCES -------------------------------------------------------

# read in data

# lepi: table with butterfly abundances. They counted all species and their individual numbers within 2.5 m either side and 5 m in front of the scientists on transects of 300 m length within 30 min. A transect was divided in 50m of 5 min intervals. 3 surveys on all plots.

# BExIS ID # 12526 (1.8.19)

# Börschig C, Klein A-M, von Wehrden H, Krauss J (2013) Traits of butterfly communities change from specialist to generalist characteristics with increasing land-use intensity. Basic and Applied Ecology 14:547–554. https://doi.org/10.1016/j.baae.2013.09.002

lepi <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/12526_Lepidoptera-2008-OpenData/12526.xlsx",
                    sheet = "12526", 
                    na = "NA")

# Plant cover table: In spring of the years 2008-2020, we sampled all species in an area of 4m x 4m and estimated the percentage cover of each species relative to the whole 4 m x 4 m plot. The vegetation records are done in early summer, normally in May and beginning of June. NA = Missing data, because already mown or not accessible because of grazing animals

# BExIS ID # 23586 (1.3.1)

plants <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/23586_Vegetation_BExIS-OpenData/23586.xlsx",
                    sheet = "data_2008-2017", 
                    na = "NA")

# tidy data

# Lepi

lepi_new <-
  lepi %>%
  # remove moths
  filter(but_moth == "B") %>% 
  arrange(species, EPID, date) %>% 
  # sum up abundances over all intervals i.e. total abundance per survey,  EPID, species
  group_by(species, EPID, survey) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>% 
  # average the summed up abundances over all surveys i.e. total abundance for the whole year per species and EPID
  group_by(species, EPID) %>%
  summarise(abundance = mean(abundance), .groups = "drop") %>% 
  mutate(# add underscore to names
    species = str_replace_all(species, ' ', "_")) %>%
  # put species names as cols
  pivot_wider(
    .,
    names_from = species,
    values_from = abundance,
    values_fill = 0
  ) %>% 
  # sum values of Pyrgus_alveus + Pyrgus_alveus_Komplex because is hard to separate
  mutate(Pyrgus_alveus_Komplex = Pyrgus_alveus + Pyrgus_alveus_Komplex) %>%
  # drop unneeded cols
  select(-Pyrgus_alveus) %>% 
  rename(Leptidea_sinapis_complex = `Leptidea_sinapis/reali`,
         Melitaea_athalia = `Melitaea_athalia-Komplex`,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex) %>%
  # round values to have no decimals
  mutate(across(where(is.numeric), round, 0)) %>%
  # arrange by plot
  arrange(EPID) 

# plants

plant_new <-
  plants %>%
  filter(Year == "2008") %>% 
  # per plant species, average all cover values for all sites
  group_by(Species, Useful_EP_PlotID) %>%
  summarise(abundance = mean(Cover), .groups = "drop") %>%
  rename(species = Species,
         EPID = Useful_EP_PlotID) %>%
  pivot_wider(
    .,
    names_from = species,
    values_from = abundance,
    values_fill = 0
  ) %>%
  # remove unnecessary observations
  select(
    -c(
      Brassicaceae_sp,
      Asteraceae_sp,
      Caryophyllaceae_sp,
      Orchidacea_sp,
      Baumkeimling_sp,
      Unknown,
      Rhinanthus_aggr.,
      # Tree species
      Acer_sp,
      Betula_pendula,
      Carpinus_betulus,
      Fraxinus_excelsior,
      Pinus_sylvestris,
      Populus_tremula,
      Prunus_avium,
      Prunus_sp,
      Prunus_spinosa,
      Quercus_robur,
      Tilia_sp,
      # Shrubs species
      Juniperus_communis,
      Crataegus_sp,
      # fern species
      Ophioglossum_vulgatum
    )
  ) %>%
  # rename species names to have only genus_species
  rename(
    Achillea_millefolium = Achillea_millefolium_aggr.,
    Alchemilla_vulgaris = Alchemilla_vulgaris_aggr.,
    Allium_oleraceum = Allium_cf_oleraceum,
    Arabis_hirsuta = Arabis_hirsuta_aggr.,
    Bromus_hordeaceus = Bromus_hordeaceus_aggr.incl_B_commutatus,
    Carex_muricata = Carex_muricata_aggr.,
    Cerastium_pumilum = Cerastium_pumilum_aggr.,
    Euphrasia_rostkoviana = Euphrasia_rostkoviana_aggr.,
    Euphrasia_sp = Euphrasia_sp_cf,
    Festuca_ovina = Festuca_ovina_aggr.,
    Festuca_rubra = Festuca_rubra_aggr.,
    Galium_mollugo = Galium_mollugo_aggr.,
    Hordeum_vulgare = Hordeum_vulgare_aggr.,
    Leucanthemum_vulgare = Leucanthemum_vulgare_aggr.,
    Medicago_sativa = Medicago_sativa_aggr.,
    Muscari_neglectum = Muscari_neglectum_aggr.,
    Odontites_vernus = Odontites_vernus_aggr.,
    Ononis_spinosa = Ononis_repens_spinosa_aggr.,
    Phleum_pratense = Phleum_pratense_aggr.,
    Pimpinella_saxifraga = Pimpinella_saxifraga_aggr.,
    Poa_pratensis = Poa_pratensis_aggr.,
    Polygala_comosa = Polygala_comosa_aggr.,
    Potentilla_verna = Potentilla_verna_aggr.,
    Rosa_canina = Rosa_canina_aggr.,
    Vicia_cracca = Vicia_cracca_aggr.,
    Vicia_sativa = Vicia_sativa_aggr.,
    Vicia_tetrasperma = Vicia_tetrasperma_aggr.
    )

# join both df

occur_new <- 
  lepi_new %>% 
  inner_join(
    plant_new, 
    by = "EPID"
  )

# write_csv(occur_new, 'data/raw/occur.csv', quote = "none") # use quote = none to avoid "" on NAs)





# 2. LANDUSE -------------------------------------------------------------

#' **Aim**: to create a dataset containing all information for experimental 
#' sites including lad use data

# Read in data

vogt <- read_delim('resources_troph-cost/338774_Vogt_et_at_2019_Landuse-OpenData/338774_Vogt_et_at_2019.csv')

obs <- 
  read_csv('data/raw/cooccurAll.csv')

# Data tidying

vogt <-
  vogt %>% 
  # rename variable as in other datasets
  rename(EPID = PlotID) %>% 
  # filter out region HAI and EPIDs that don't have occurrences
  filter(StudyRegion != 'HAI') %>% 
  # add a leading zero 
  mutate(
    EPID = EPID %>%
      str_replace('(AEG|SEG)([1-9]{1})$', # create groups for regex
                  '\\10\\2')) %>%  # add a leading zero between groups
  # create a new ID col from unite Year and EPID
  unite(.,
        EPID, 
        Year,
        col = 'ID_new',
        sep = '_',
        remove = FALSE) %>% 
  # cleaning up...
  select(-ID) %>% 
  relocate(ID_new, StudyRegion:DescSeeds) %>% 
  rename(ID = ID_new) %>% 
  # set dates as date (year, month, day)
  mutate(across(c(Date, DateMowing1, DateMowing2, DateMowing3,DateMowing4),
    ~ dmy(.x)))

# write_csv(vogt, file = 'data/raw/sites_landuse.csv')

# Create vectos with EPIDs for site SCH and ALB where species have occurrences

EP_to_filter <- 
  obs %>% 
  filter(
    str_detect(EPID, 'SEG|AEG')) %>% 
  column_to_rownames('EPID') %>% 
  select(
    where(
      ~ sum(.x) !=0)) %>% 
  row.names()

# filter sites from year 2007 and 2008 and that have occurrences

sites <- 
  vogt %>% 
  filter(Year == 2008 | Year == 2007)  %>% 
  filter(EPID %in% EP_to_filter)

# write_csv(sites, file = 'data/raw/sites_landuse_2007-2008.csv')







# 3. TROPHIC INTERACTIONS -------------------------------------------------

# NOTE: this data were extracted from Ebert (2005) for ALB, and Richert and Brauner (2018) for SCH, then collected in a single excel file which then was here cleaned up. Cleaned versions with a single file per region for larva and adults were stored in BExIS as long format and can be found using the ID: 

# Imago ALB: 31734
# Larva ALB: 31735
# Imago SCH: 31736
# Larva SCH: 31737

# we then use them in ~troph-cost/scripts/wrangling/pairs_1_trophic_links.R in wide format to perform the analyses of co-occurrence and other posterior analyses. 

## 3.1 Larva_BB (Berlin Brandenburg = BB = SCH) 

Larva_BB <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Larva_BB", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.2 Larva_BW (Baden-Württemberg = BW = ALB) 
Larva_BW <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Larva_BW", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.3 Imago_BB  (Berlin Brandenburg = BB = SCH) 

Imago_BB <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Imago_BB", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.4 Imago_BW (Baden-Württemberg = BW = ALB) 

Imago_BW <- 
  readxl::read_xlsx("/Users/estebanmenaresbarraza/Documents/BTU-PhD/troph-cost/resources_troph-cost/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Imago_BW", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

# check the number of plants and butterflies (dimensions without NAs)

dim_no_NA <- function(x){
  x %>%
    as_tibble() %>% 
    select(Anthocharis_cardamines:last_col()) %>%
    select_if(~any(!is.na(.))) %>% 
    filter(!is.na(Anthocharis_cardamines)) %>% 
    dim()
} 

dim_no_NA(Imago_BW)
dim_no_NA(Imago_BB)
dim_no_NA(Larva_BW)
dim_no_NA(Larva_BB)

# check the number of plants and butterflies (dimensions without NAs, no Zero)

dim_no_NA_no_zero <- function(x, ref_NA_row = ref_NA_row){
  x %>%
    as_tibble() %>% 
    select(all_of(ref_NA_row):last_col()) %>%
    select(where( ~ any(!is.na(.)))) %>% 
    filter(!is.na(ref_NA_row)) %>% 
    select(where(~ sum(., na.rm = TRUE) > 0)) %>% 
    filter(rowSums(across(everything())) > 0) %>%
    dim()
} 

dim_no_NA_no_zero(Imago_BW, ref_NA_row = "Anthocharis_cardamines")
dim_no_NA_no_zero(Imago_BB, ref_NA_row = "Anthocharis_cardamines")
dim_no_NA_no_zero(Larva_BW, ref_NA_row = "Anthocharis_cardamines")
dim_no_NA_no_zero(Larva_BB, ref_NA_row = "Anthocharis_cardamines")


# replace any "," or "/" for ";" to avoid readability problems when exporting as csv

Larva_BB <- Larva_BB %>% 
  mutate(across(where(is.character), ~str_replace_all(.x, "/|,", ";")))

Larva_BW <- Larva_BW %>% 
  mutate(across(where(is.character), ~str_replace_all(.x, "/|,", ";")))

Imago_BB <- Imago_BB %>% 
  mutate(across(where(is.character), ~str_replace_all(.x, "/|,", ";")))

Imago_BW <- Imago_BW %>% 
  mutate(across(where(is.character), ~str_replace_all(.x, "/|,", ";")))

## ---- Save data

# Save data with wide format

write_excel_csv(Larva_BB, file = "resources_troph-cost/trophic-interactions/Larva_BB.csv")
write_excel_csv(Larva_BW, file = "resources_troph-cost/trophic-interactions/Larva_BW.csv")
write_excel_csv(Imago_BB, file = "resources_troph-cost/trophic-interactions/Imago_BB.csv")
write_excel_csv(Imago_BW, file = "resources_troph-cost/trophic-interactions/Imago_BW.csv")

# Save data with long format

Larva_BB %>%
  pivot_longer(cols = Anthocharis_cardamines:last_col(),
               names_to = "lepidoptera_species",
               values_to = "int_strength") %>% 
  write_excel_csv(., file = "resources_troph-cost/trophic-interactions/Larva_BB_long.csv")

Larva_BW %>%
  pivot_longer(cols = Anthocharis_cardamines:last_col(),
               names_to = "lepidoptera_species",
               values_to = "int_strength") %>% 
  write_excel_csv(., file = "resources_troph-cost/trophic-interactions/Larva_BW_long.csv")

Imago_BB %>%
  pivot_longer(cols = Anthocharis_cardamines:last_col(),
               names_to = "lepidoptera_species",
               values_to = "int_strength") %>% 
  write_excel_csv(., file = "resources_troph-cost/trophic-interactions/Imago_BB_long.csv")

Imago_BW %>%
  pivot_longer(cols = Anthocharis_cardamines:last_col(),
               names_to = "lepidoptera_species",
               values_to = "int_strength") %>% 
  write_excel_csv(., file = "resources_troph-cost/trophic-interactions/Imago_BW_long.csv")







# 4. MOST VISITED PLANTS --------------------------------------------------

popular_plants_SCH <- read_xlsx("data/raw/popular_plants_SCH.xlsx")
dim(popular_plants_SCH)
str(popular_plants_SCH)

# write_csv(popular_plants_SCH, file = "data/popular_plants_SCH.csv")




# 5. FLOWER AVAILABILITY ---------------------------------------------------

#' **Aim:**
#' per species and EP, calculate average of flowering units across all 
#' dates and turn dataset into long format dataset: obs at each site for
#' each plant and butterfly species.


## ---- Load data

# BExIS dataset 4981 (v2) "Flower availability 2008 Alb-korrigiert"
flowers_ALB <- readxl::read_xlsx("resources_troph-cost/4981_2_FlowerAvailability_ALB-2008-OpenData/4981_2_data.xlsx")
str(flowers_ALB)

# BExIS dataset 4964 (v2) "Flower availability 2008 Schorfheide"
flowers_SCH <- readxl::read_xlsx("resources_troph-cost/4964_2_FlowerAvailability_SCH-2008-OpenData/4964_2_data.xlsx")
str(flowers_SCH)


## ---- Tidy data 

# Group data per EP and plant species, and calculate mean across all dates

#ALB
flowers_ALB_filtered <- flowers_ALB %>% 
  group_by(EP_ID, Species, .drop = FALSE) %>%
  summarise(flower_units = 
              mean(Flowering_unit))

#SCH
flowers_SCH_filtered <- flowers_SCH %>% 
  group_by(EP_ID, Species, .drop = FALSE) %>%
  summarise(flower_units = 
              mean(Flowering_unit))

# Change variable names

#ALB
colnames(flowers_ALB_filtered)[which(
  names(flowers_ALB_filtered) == "Species")] <- "plant_sp"

colnames(flowers_ALB_filtered)[which(
  names(flowers_ALB_filtered) == "EP_ID")] <- "EPID"
# SCH
colnames(flowers_SCH_filtered)[which(
  names(flowers_SCH_filtered) == "Species")] <- "plant_sp"

colnames(flowers_SCH_filtered)[which(
  names(flowers_SCH_filtered) == "EP_ID")] <- "EPID"

# Make plant_sp names
#ALB
flowers_ALB_filtered$plant_sp <- 
  str_replace_all(flowers_ALB_filtered$plant_sp, " ", "_")
#SCH
flowers_SCH_filtered$plant_sp <- 
  str_replace_all(flowers_SCH_filtered$plant_sp, " ", "_")

# Trim empty spaces at the right side of plant_sp names
#ALB
flowers_ALB_filtered$plant_sp <- 
  str_trim(flowers_ALB_filtered$plant_sp, side = "right")
#SCH
flowers_SCH_filtered$plant_sp <- 
  str_trim(flowers_SCH_filtered$plant_sp, side = "right")

# Check unique names EPID
unique(flowers_ALB_filtered$EPID)
unique(flowers_SCH_filtered$EPID)

# Replace values of EPID

#ALB
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG1")] <- c("AEG01")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG2")] <- c("AEG02")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG3")] <- c("AEG03")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG4")] <- c("AEG04")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG6")] <- c("AEG06")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG7")] <- c("AEG07")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG8")] <- c("AEG08")
flowers_ALB_filtered$EPID[flowers_ALB_filtered$EPID == c("AEG9")] <- c("AEG09")

#SCH
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG1")] <- c("SEG01")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG2")] <- c("SEG02")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG3")] <- c("SEG03")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG4")] <- c("SEG04")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG6")] <- c("SEG06")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG7")] <- c("SEG07")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG8")] <- c("SEG08")
flowers_SCH_filtered$EPID[flowers_SCH_filtered$EPID == c("SEG9")] <- c("SEG09")

# Check unique names EPID again
unique(flowers_ALB_filtered$EPID)
unique(flowers_SCH_filtered$EPID)

# Sort dataset by EPID and plant_sp
flowers_ALB_filtered <- arrange(flowers_ALB_filtered, EPID, plant_sp)
flowers_SCH_filtered <- arrange(flowers_SCH_filtered, EPID, plant_sp)


## ---- Save data
# fwrite(flowers_ALB_filtered, file = "data/raw/flower_units_ALB.csv")
# fwrite(flowers_SCH_filtered, file = "data/raw/flower_units_SCH.csv")




# 6. FLOWER VISITATION --------------------------------------------------


## ---- Read data 

# BExIS ID 10160 (v3)
flower_visitor <- 
  read_delim(file = "resources_troph-cost/10160_FlowerVisitorInteractions-2008-OpenData/10160_3_data.csv",
           delim = ';')
str(flower_visitor)


## ---- filter only Lepidoptera species
flower_visitor <- 
  flower_visitor %>% 
  filter(Order == "Lepidoptera")



## ---- harmonize plant species names
# check names
unique(flower_visitor$Plant)

flower_visitor$Plant <- 
  str_replace_all(flower_visitor$Plant, " ", "_")

flower_visitor$Plant <- 
  str_replace_all(flower_visitor$Plant, "_ssp._orientalis", "")

flower_visitor$Plant <- 
  str_replace_all(flower_visitor$Plant, "_agg.", "")

flower_visitor$Plant <- 
  str_replace_all(flower_visitor$Plant, "_cf._", "_")

flower_visitor$Plant <- 
  str_trim(flower_visitor$Plant, side = "right")

# check names again 
unique(flower_visitor$Plant)



## ---- harmonize lepidoptera species names
# check names 
unique(flower_visitor$Species)

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, " ", "_")

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, "/minos_agg.", "")

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, "/argyrognomon/ideas_agg.", "")

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, "/argyrognomom/ideas", "")

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, "_agg.", "")

flower_visitor$Species <- 
  str_replace_all(flower_visitor$Species, "sp.", "sp")

# check names again 
unique(flower_visitor$Species)

# remove family-level observations
flower_visitor <- 
  flower_visitor[!grepl("_undet.", flower_visitor$Species),]

# check EP_ID names

unique(flower_visitor$EP_ID)

## ---- tidy dataset

flower_visitor <- 
  flower_visitor %>% 
  # filter out EP_IDs from region HAI
  filter(
    !str_detect(EP_ID, 'HEG')) %>% 
  # Add a leading 0 to numeric IDs to have all names with 2 digits
  mutate(
    EP_ID = EP_ID %>%
      str_replace('(AEG|SEG)([1-9]{1})$', # create groups for regex
                  '\\10\\2')) %>%         # add zero between groups
  # arrange data 
  arrange(EP_ID, Plant, Species) %>% 
  # create a region col
  mutate(
    Region = 
      if_else(
        str_detect(EP_ID, 'AEG.*') == TRUE,
        'ALB',
        'SCH')) %>% 
  distinct() %>% 
  # put region at the front
  relocate(Region, EP_ID:Total) %>% 
  # rename cols
  rename(plant_sp = Plant,
         lepidoptera_sp = Species)

# Check if each row is a unique observation

flower_visitor %>% 
  distinct()

## ---- save data
# write_csv(flower_visitor, "data/raw/flower_visitation.csv", quote = "none")






# 7. RED LIST  -----------------------------------------------------------

# In order to type information of red lists for each plant and lepidoptera
# species we will: 
# Create a dataframe with all plant and lepidoptera species (unique 
# distinct) values which have a trophic interaction AND occurr in a 
# sampling site. The dataframe should consider both study regions.

# requierements
library(tidyverse)

# load data
imago_ALB <- read_csv(file = "data/imago_ALB.csv")
imago_SCH <- read_csv(file = "data/imago_SCH.csv")

# get the unique distinct names of plants and butterflies to create vectors
plants_ALB <- unique(imago_ALB$plant_sp)
lepi_ALB <- unique(imago_ALB$lepidoptera_sp)

plants_SCH <- unique(imago_SCH$plant_sp)
lepi_SCH <- unique(imago_SCH$lepidoptera_sp)

region_ALB <- rep("ALB", length(c(plants_ALB, lepi_ALB)))
region_SCH <- rep("SCH", length(c(plants_SCH, lepi_SCH)))

org_type <- c(rep("plant_sp", length(plants_ALB)), 
              rep("lepidoptera_sp", length(lepi_ALB)), 
              rep("plant_sp", length(plants_SCH)), 
              rep("lepidoptera_sp", length(lepi_SCH)))

# create a data frame containing all vectors
red_list <- data.frame(region = c(region_ALB, 
                                  region_SCH),
                       taxa = org_type, 
                       species = c(plants_ALB, 
                                   lepi_ALB, 
                                   plants_SCH, 
                                   lepi_SCH), 
                       red_list_detail = NA,
                       red_list = NA, 
                       distribution = NA)

# save data frame

# write_csv(red_list, file = "data/raw/red_list_info.csv", quote = "none")

# Work continued in ~/scripts/wrangling/scores1_red_list.R




# 8. CROPS -----------------------------------------------------------

# variable:	description

# region:	study region (ALB/SCH) representing the respective state of germany (Baden-Württemberg/Brandenburg)
# crop_name:	crop scientific name
# crop_group:	crop grouping common name according to data provided by Rader et al. (2020) whenever a single species presented many varieties or subsp. If no group found, the same name as in crop_en
# crop_en:	crop common name in english 
# crop_de:	crop common name in german 
# UAA_2008: utilized agricultural area in year 2008
# UAA_2019: utilized agricultural area in year 2019
# mean_UAA_2015_2020:	mean utilized agricultural area of a crop in hectars for years 2015-2020
# prop_area_state:	proportion of a crop's mean_UAA_2015_2020 from the total state's area
# prop_area_UAA:	proportion of a crop's mean_UAA_2015_2020 from the total UAA of the state
# comments:	details and descriptions for each raw. 

crops <- 
  readxl::read_excel('resources_troph-cost/crops_pollinators/crops.xlsx',
                     sheet = 'data',
                     na = 'NA') %>% 
  as_tibble() 

crops <-
  crops %>% 
  # round values 
  mutate(
    UAA_2008 = round(UAA_2008, 1),
    UAA_2019 = round(UAA_2019, 1),
    mean_UAA_2015_2020 = round(mean_UAA_2015_2020, 1),
    prop_area_state = round(prop_area_state, 5),
    prop_area_UAA = round(prop_area_UAA, 5)) 

# write_excel_csv(crops, file = 'data/raw/crops.csv', na = 'NA')

# Work continued in script UAA_crops_regions.R...




# 9. PAIRS ---------------------------------------------------------------

# Read in data

library(tidyverse)

species <- 
  read_csv('data/raw/species.csv')

pairs_ALB <- 
  read_csv('data/raw/imago_ALB.csv')

pairs_SCH <- 
  read_csv('data/raw/imago_SCH.csv')

# ALB

pairs_ALB <-
  pairs_ALB %>% 
  # separate troph_pair id into plant and butterfly id
  separate(troph_pair, 
           into = c('plant_sp_id', 'lepi_sp_id'), 
           remove = TRUE) %>%
  # paste id cols back to one single pair id
  unite(., col = "troph_pair_id", 
        plant_sp_id, lepi_sp_id, 
        sep = "-",
        remove = FALSE) %>% 
  # reorder cols to put troph_pair at the beginning
  relocate(troph_pair_id, plant_sp:last_col())

# SCH

pairs_SCH <-
  pairs_SCH %>% 
  # separate troph_pair id into plant and butterfly id
  separate(troph_pair, 
           into = c('plant_sp_id', 'lepi_sp_id'), 
           remove = TRUE) %>%
  # paste id cols back to one single pair id
  unite(., col = "troph_pair_id", 
        plant_sp_id, lepi_sp_id, 
        sep = "-",
        remove = FALSE) %>% 
  # reorder cols to put troph_pair at the beginning
  relocate(troph_pair_id, plant_sp:last_col())


## ---- create a single table with all information

pairs_ALB_new <- 
  pairs_ALB %>% 
  # add a var called region
  mutate(region = "ALB") %>% 
  # put at front
  relocate(region, troph_pair_id:last_col())

pairs_SCH_new <- 
  pairs_SCH %>% 
  # add a var called region
  mutate(region = "SCH") %>% 
  # put at front
  relocate(region, troph_pair_id:last_col())

# combine tables into one called "pairs"

pairs <- 
  bind_rows(pairs_ALB_new, pairs_SCH_new)

# remove unnecessary tables

rm(pairs_SCH, pairs_ALB, pairs_ALB_new, pairs_SCH_new)


## ---- Check and harmonize names of pairs and species tables

# 1. Plants ALB

# get plant species of region ALB from species table

species_plant_sp <-
  species %>% 
  filter(region == "ALB" & taxa == "plant") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

# get plant species of region ALB from pairs table

pairs_plant_sp <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(plant_sp) %>% 
  distinct() %>% 
  pull()

setdiff(species_plant_sp, pairs_plant_sp)
setdiff(pairs_plant_sp, species_plant_sp) 

# ALl GOOD!

# 2. Lepi ALB

# get lepi species of region ALB from species table

species_lepi_sp <-
  species %>% 
  filter(region == "ALB" & taxa == "lepidoptera") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

# get lepi species of region ALB from pairs table

pairs_lepi_sp <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(lepidoptera_sp) %>% 
  distinct() %>% 
  pull()

setdiff(species_lepi_sp, pairs_lepi_sp) 
# ALL GOOD

setdiff(pairs_lepi_sp, species_lepi_sp) 
# ALL GOOD


# 3. Plants SCH

# get plant species of region SCH from species table

species_plant_sp <-
  species %>% 
  filter(region == "SCH" & taxa == "plant") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

# get plant species of region SCH from pairs table

pairs_plant_sp <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(plant_sp) %>% 
  distinct() %>% 
  pull()

setdiff(species_plant_sp, pairs_plant_sp)
# no differences

setdiff(pairs_plant_sp, species_plant_sp) 
# no differences

# 4. Lepi SCH

# get lepi species of region SCH from species table

species_lepi_sp <-
  species %>% 
  filter(region == "SCH" & taxa == "lepidoptera") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

# get lepi species of region SCH from pairs table

pairs_lepi_sp <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(lepidoptera_sp) %>% 
  distinct() %>% 
  pull()

setdiff(species_lepi_sp, pairs_lepi_sp) 
# no differences

setdiff(pairs_lepi_sp, species_lepi_sp) 
# no differences


# round non-integer values to 5 decimal places  

pairs_new <- pairs %>% 
  mutate(
    across(c(sites_plant, sites_lepi, obs_cooccur, flower_visit), 
    as.integer)) %>% 
  mutate(
    across(
      where(is.double), ~sprintf("%.5f", .x))) 

## ---- save data 

write_excel_csv(pairs_new, 'data/raw/pairs.csv', na = "NA")





# 10. SPECIES AND SCORES_SP -----------------------------------------------

# Species table: a table extracted manually with all the species names, their family, genus and name, german names, synonyms, and comments for each of the species found in the study. Species names and synonyms were checked using the Lepiforum (Rennwald and Rodeland n.d., https://lepiforum.org/) and the Global Lepidoptera Index (Beccaloni et al. 2022). Leptidea sinapis, L. reali and L. juvernica were considered as L. sinapis complex according to Dincă et al. (2011). Pyrgus alveus, P. trebevicensis, P. accrete were considered as P. alveus complex according to Tshikolovets (2011). Species names and synonyms were checked using the International Plant Names Index (IPNI 2023) and Plants of the World Online database (POWO 2023). 


# Setup

library(tidyverse)

species <- 
  readr::read_csv('data/raw/species.csv') #%>% 
  # # rename var
  # rename(species_name = species)

pairs <- 
  read_csv('data/raw/pairs.csv')


# Data tidying 

## ---- Add the species id to the species table 

# ALB

species_new_ALB <- 
  species %>% 
  filter(region == "ALB") %>% 
  select(-species_id) %>% 
  # add plant species id
  left_join(
    pairs %>% 
      filter(region == "ALB") %>% 
      select(plant_sp, plant_sp_id) %>% 
      distinct(),
    by = c("species_name" = "plant_sp")) %>%
  # add lepidoptera species id
  left_join(
    pairs %>% 
      filter(region == "ALB") %>% 
      select(lepidoptera_sp, lepi_sp_id) %>% 
      distinct(),
    by = c("species_name" = "lepidoptera_sp")) %>% 
  # unite both cols into one
  unite(., col = "species_id", 
        plant_sp_id, lepi_sp_id, 
        sep = "",
        remove = TRUE,
        na.rm = TRUE)

# check data ALB

plant_check <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(plant_sp) %>% 
  distinct() %>% 
  pull()

plant_check_id <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(plant_sp_id) %>% 
  distinct() %>% 
  pull()

lepi_check <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(lepidoptera_sp) %>% 
  distinct() %>% 
  pull()

lepi_check_id <- 
  pairs %>% 
  filter(region == "ALB") %>% 
  select(lepi_sp_id) %>% 
  distinct() %>% 
  pull()

pairs_check <- c(plant_check, lepi_check)
pairs_check_id <- c(plant_check_id, lepi_check_id)

species_check <- 
  species_new_ALB %>% 
  filter(region == "ALB") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

species_check_id <- 
  species_new_ALB %>% 
  filter(region == "ALB") %>% 
  select(species_id) %>% 
  distinct() %>% 
  pull()

identical(pairs_check, species_check)
identical(pairs_check_id, species_check_id)

# SCH

species_new_SCH <-
  species %>% 
  filter(region == "SCH") %>% 
  select(-species_id) %>% 
  # add plant species id
  left_join(
    pairs %>% 
      filter(region == "SCH") %>% 
      select(plant_sp, plant_sp_id) %>% 
      distinct(),
    by = c("species_name" = "plant_sp")) %>%
  # add lepidoptera species id
  left_join(
    pairs %>% 
      filter(region == "SCH") %>% 
      select(lepidoptera_sp, lepi_sp_id) %>% 
      distinct(),
    by = c("species_name" = "lepidoptera_sp")) %>% 
  # unite both cols into one
  unite(., col = "species_id", 
        plant_sp_id, lepi_sp_id, 
        sep = "",
        remove = TRUE,
        na.rm = TRUE)

# check data SCH

plant_check <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(plant_sp) %>% 
  distinct() %>% 
  pull()

plant_check_id <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(plant_sp_id) %>% 
  distinct() %>% 
  pull()

lepi_check <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(lepidoptera_sp) %>% 
  distinct() %>% 
  pull()

lepi_check_id <- 
  pairs %>% 
  filter(region == "SCH") %>% 
  select(lepi_sp_id) %>% 
  distinct() %>% 
  pull()

pairs_check <- c(plant_check, lepi_check)
pairs_check_id <- c(plant_check_id, lepi_check_id)

species_check <- 
  species_new_SCH %>% 
  filter(region == "SCH") %>% 
  select(species_name) %>% 
  distinct() %>% 
  pull()

species_check_id <- 
  species_new_SCH %>% 
  filter(region == "SCH") %>% 
  select(species_id) %>% 
  distinct() %>% 
  pull()

identical(pairs_check, species_check)
identical(pairs_check_id, species_check_id)

# Combine cases from both tables back into one single species table 

species_new <- 
  # combine both tables
  bind_rows(species_new_ALB, species_new_SCH) %>% 
  # move species_id to the front
  relocate(region, species_id:last_col())

# remove unnecessary objects

rm(species_new_ALB, species_new_SCH, lepi_check, pairs_check, plant_check, 
   species_check)

# final clean up species table

species_new <-
  species_new %>% 
  # # tidy de_name var
  # separate(col = de_name, 
  #          sep = ',', 
  #          into = c("de_name1", "de_name2", "de_name3", "de_name4")) %>%
  relocate(species_id, region, taxa:last_col()) %>% 
  # select vars for species df
  select(species_id:de_name4, comment) %>%
  # clean up comment col
  mutate(comment = comment %>% 
           str_replace("^red_list.+",
                       NA_character_) %>% 
           str_replace("; red_list.+",
                       "")) 
# because 
species_new <-
  species_new %>% 
    mutate(comment = str_replace_all(comment, ",", ";"))

# write data into csv

write_excel_csv(species_new, "data/raw/species.csv", na = NA_character_)


# PART 2: SCORES.

# This part cannot be run again because the old species table had the info about red list species on it which I added manually and then removed with this script. I separated then species, a table containing only the species ids, sci names, synonyms and german names, from scores, a table containing species ids, sci names and different "scores" or information about certain  biological or conservation aspects. DO NOT run again this. 

# Separate species table into species and scores table i.e. remove data that
# is of a different level of observation
# 
scores <-
  species %>%
  # reorder vars
  relocate(species_id, region:last_col()) %>%
  # select vars for scores df
  select(species_id:species_name, red_list:last_col())

# Score crop_pest was done manually 

# write data into csv

# write_csv(scores, "data/raw/scores.csv", na = NA_character_)






# 11. SELF CO-OCCURRENCES -------------------------------------------------

# tidy table of lepi-lepi and plant-plant species co-occurrences

# 1) Rename butterfly names with "komplex" and "sinapis reali" with:
#   - Leptidea_sinapis_reali >>> Leptidea_sinapis_complex
#   - Melitaea_athalia_Komplex >>> Melitaea_athalia
#   - Pyrgus_alveus_Komplex >>> Pyrgus_alveus_complex
# 2) remove genus-level plant species "aggr."
# 2) Update the troph_pair id to the newer one with leading zeros and region

# read in data 

self_cooccur <- 
  read_csv('data/raw/self_cooccur.csv')

# data tidying 

self_cooccur <- 
  self_cooccur %>% 
  # rename butterfly names 
  mutate(sp1_name = sp1_name %>% 
           str_replace("Leptidea_sinapis_reali", 
                       "Leptidea_sinapis_complex") %>% 
           str_replace("Melitaea_athalia_Komplex",
                       "Melitaea_athalia") %>% 
           str_replace("Pyrgus_alveus_Komplex", 
                       "Pyrgus_alveus_complex"),
         sp2_name = sp2_name %>% 
           str_replace("Leptidea_sinapis_reali", 
                       "Leptidea_sinapis_complex") %>% 
           str_replace("Melitaea_athalia_Komplex",
                       "Melitaea_athalia") %>% 
           str_replace("Pyrgus_alveus_Komplex", 
                       "Pyrgus_alveus_complex")) %>% 
  # remove genus-level plant species "aggr."
  filter(!str_detect(sp1_name, "aggr.")) %>% 
  filter(!str_detect(sp2_name, "aggr.")) %>%
  # add ids with leading zeros and region for each col of species name
  left_join(
    species %>% 
      select(region, species_id, species_name),
    by = c("sp1_name" = "species_name",
           "region" = "region")) %>%
  left_join(
    species %>% 
      select(region, species_id, species_name),
    by = c("sp2_name" = "species_name",
           "region")) %>% 
  # rename joined vars
  rename(sp1_id = species_id.x,
         sp2_id = species_id.y) %>% 
  # reorder cols
  relocate(region:sp2_name, sp1_id:sp2_id, co_occurrence:last_col()) 

# Save data

# write_csv(self_cooccur, 'data/raw/self_cooccur.csv')




# 12. LUI -------------------------------------------------------------

## 12.1  LUI ------------------------------------------------------------

lui <- 
  read_delim('resources_troph-cost/landuse/LuiData/LUI_default components set_regional_overall_2022-09-07.txt')

lui <- 
  lui %>% 
  BEplotZeros(dat = ., column = "PLOTID", colname = "PLOTID") %>% 
  as_tibble() %>% 
  rename(EPID = PLOTID) %>% 
  arrange(EPID) %>% 
  select(-c(YEAR, EXPLO))

# write_csv(lui, 'data/raw/lui.csv')


## 12.2 LUI regional --------------------------------------------------

lui_reg <- 
  read_delim('resources_troph-cost/landuse/LUI_regional/LUI_default components set_regional_overall_2023-02-09.txt')

lui_reg <-
  lui_reg %>% 
  BEplotZeros(dat = ., column = "PLOTID", colname = "PLOTID") %>% 
  as_tibble() %>% 
  rename(EPID = PLOTID) %>% 
  arrange(EPID) %>% 
  select(-c(YEAR, EXPLO))
  
# write_csv(lui_reg, 'data/raw/lui_reg.csv')
  


# 13. SITES ------------------------------------------------------

sites <-
  read_delim('resources_troph-cost/EPs/20826_6_Dataset_EPs_info/20826_6_data.csv') %>%
  # add a leading zero
  mutate(EP_PlotID = EP_PlotID %>%
           str_replace_all('(AEG|SEG|HEG)([1-9]{1})$', # create groups for regex
                           '\\10\\2'))  # add a leading zero between groups
  
env_vars <- read_delim('resources_troph-cost/EPs/EP_env_vars/plots.csv')

env_vars2 <- 
  read_delim('resources_troph-cost/EPs/20826_6_Dataset_EPs_info/20826_6_data.csv') %>% 
  # add a leading zero 
  mutate(
    EP_PlotID = EP_PlotID %>%
      str_replace_all('(AEG|SEG|HEG)([1-9]{1})$', # create groups for regex
                      '\\10\\2')) %>%  # add a leading zero between groups
  rename_with(tolower, .cols = everything())
  
landuse <- read_delim('data/raw/sites_landuse.csv') %>% 
  filter(Year == 2008)

LUI_2008 <- read_csv('data/raw/lui_reg.csv')

env_vars3 <- 
  read_delim('resources_troph-cost/31018_5_env_and_land-use_covariates_grassland_EPs/31018_5_data.csv') %>% 
  # add a leading zero 
  mutate(
    EP_PlotID = EP_PlotID %>%
      str_replace_all('(AEG|SEG|HEG)([1-9]{1})$', # create groups for regex
                      '\\10\\2'))   # add a leading zero between groups

env_vars4 <- 
  read_delim('resources_troph-cost/18148_2_Dataset/18148_2_data.csv') %>% 
  distinct() %>% 
  # add a leading zero 
  mutate(
    EPID = EP_PLOTID %>%
      str_replace_all('(AEG|SEG|HEG)([1-9]{1})$', # create groups for regex
                      '\\10\\2'))   # add a leading zero between groups
  

# tidy data

sites <-
  sites %>% 
  select(EPID = EP_PlotID, 
         region = Exploratory, 
         Landuse, 
         ActivePlot,
         RW, 
         HW, 
         SoilTypeWRB, 
         Elevation, 
         Slope) %>% 
  filter(Landuse == "G", region != "HAI", ActivePlot == "yes") %>% 
  select(-c(Landuse, ActivePlot)) %>% 
  as_tibble() %>% 
    rename_with(tolower, .cols = c(region:Slope)) %>% 
  arrange(region) %>% 
    left_join(
      env_vars %>% 
        select(EPID = plotID, 
               year = datetime,
               precip = precipitation_radolan,
               rain_days = precipitation_radolan_rain_days),
      by = "EPID"
    ) %>% 
    inner_join(
      landuse %>% 
        filter(Year == 2008) %>% 
        select(EPID, size_unit = SizeManagementUnit_ha, 
               DateMowing1:DateMowing3),
      by = "EPID"
    ) %>% 
  inner_join(
    LUI_2008,
    by = "EPID"
  ) %>% 
    left_join(
      env_vars2 %>% 
        filter(activeplot == "yes") %>% 
        select(EPID = ep_plotid, vip, mip, 
               biotope_id = biotopeid, 
               biotop_name = biotopname, 
               conservation_type = conservationtype,
               org_farming = ecologicalfarming),
      by = "EPID"
    ) %>% 
    left_join(
      env_vars3 %>% 
        select(EPID = EP_PlotID,
               twi = TWI, # for soil moisture
               ph = Soil.pH, 
               soil_depth = Soil.depth,
               clay = Soil.clay.content,
               sand = Soil.sand.content),
      by = "EPID"
    ) %>% 
  left_join(
    env_vars4 %>% 
      select(EPID,
             Arable.500 = A500,
             Forest.500 = F500,
             Grassland.500	= G500,
             Semi.Natural.500	= N500,
             Roads.500 = R500, 
             Woodland.500	= T500,
             Urban.500 = U500, 
             Water.500 = W500,
             Arable.1000 = A1000,
             Forest.1000 = F1000,
             Grassland.1000	= G1000,
             Semi.Natural.1000	= N1000,
             Roads.1000 = R1000, 
             Woodland.1000	= T1000,
             Urban.1000 = U1000, 
             Water.1000 = W1000),
    by = "EPID"
  )

# write_csv(sites, 'data/raw/sites.csv')




# 14.  QGIS UAA per site ----------------------------------------------

# read in data

alb_3km <-
  readxl::read_excel('resources_troph-cost/UAA/UAA_alb_buffer_3km.xlsx') %>% 
  as_tibble()

alb_5km <-
  readxl::read_excel('resources_troph-cost/UAA/UAA_alb_buffer_5km.xlsx') %>% 
  as_tibble() 

sch_3km <-
  readxl::read_excel('resources_troph-cost/UAA/UAA_sch_buffer_3km.xlsx') %>% 
  as_tibble()

sch_5km <-
  readxl::read_excel('resources_troph-cost/UAA/UAA_sch_buffer_5km.xlsx') %>% 
  as_tibble()

# tidy data

alb_3km <-
  alb_3km %>% 
  select(-geometry) %>% 
  mutate(buffer = "3km") %>% 
  relocate(buffer, .after = Y) %>% 
  mutate(Elevation = as.numeric(Elevation),
         Slope = as.numeric(Slope))

alb_5km <-
  alb_5km %>% 
  select(-geometry) %>% 
  mutate(buffer = "5km") %>% 
  relocate(buffer, .after = Y) %>% 
  mutate(Elevation = as.numeric(Elevation),
         Slope = as.numeric(Slope))

sch_3km <- 
  sch_3km %>% 
  select(-geometry) %>% 
  mutate(buffer = "3km") %>% 
  relocate(buffer, .after = Y) 

sch_5km <-
  sch_5km %>% 
  select(-geometry) %>% 
  mutate(buffer = "5km") %>% 
  relocate(buffer, .after = Y) 

UAA <- 
  bind_rows(
    alb_3km,
    alb_5km,
    sch_3km,
    sch_5km
  )

UAA <- 
  UAA %>% 
  select(-c(HISTO_NODA, HISTO_NODATA)) 

# write_csv(UAA, 'resources_troph-cost/UAA/UAA_QGIS.csv')




# 15. REG DIST PLANTS -----------------------------------------------------

# Original data of plant occurrences at regional level was obtained via mail to Rudolf May from Bundesamt für Naturschutz (BfN), and Flora and Vegetation Deutschlands (FloraWeb). This was used to calculate the regional distribution of plants by calculating the total area of TK25 which the species occupy as a proportion of the total regional area. 

#' *STEPS*
#'  - read files
#'  - tidy names
#'  - set minimum share of federal state in the total area (XX_PART_AREA):
#'  for this we take just any species with area > 0
#'  - set time periods: from 1980 on 
#'  - set which statuses: any
#'  - per species, calculate the total area of TK25 which the species occupy

# Read files 

regdist_sch <- 
  readxl::read_excel("resources_troph-cost/regional_distribution/plants/export_menares_bb.xlsx") %>% 
  as_tibble() %>% 
  # rename wrongly named var
  rename(BB_PART_AREA = BW_PART_AREA) %>% 
  # change decimal separator from , to .
  mutate(GRID_AREA = GRID_AREA %>% 
           str_replace(., ",", ".") %>% 
           as.numeric(.),
         BB_PART_AREA = BB_PART_AREA %>% 
           str_replace(., ",", ".") %>% 
           as.numeric(.))

regdist_alb <- 
  readxl::read_excel("resources_troph-cost/regional_distribution/plants/export_menares_bw.xlsx") %>% 
  as_tibble() %>% 
  # change decimal separator from , to .
  mutate(GRID_AREA = GRID_AREA %>% 
           str_replace(., ",", ".") %>% 
           as.numeric(.),
         BW_PART_AREA = BW_PART_AREA %>% 
           str_replace(., ",", ".") %>% 
           as.numeric(.))

# Explore and tidy files 

names(regdist_sch)
names(regdist_alb)

str(regdist_sch)
str(regdist_alb)

regdist_sch %>% distinct(ZEITRAUM, ZEITRAUM_TEXT)
regdist_alb %>% distinct(ZEITRAUM, ZEITRAUM_TEXT)

# check duplicated values 

regdist_sch %>% 
  distinct(SUCHNR, SCI_NAME) %>%
  add_count(SCI_NAME) %>% 
  filter(n > 1) 

regdist_alb %>% 
  distinct(SUCHNR, SCI_NAME) %>% 
  add_count(SCI_NAME) %>% 
  filter(n > 1) # Muscari neglectum has 2 numbers, why? -> because in FloraWeb the SUCHNR of agg. lead to M. neglectum, so both are the same. 

summary(regdist_sch)
summary(regdist_alb)


# Explore the total area of grid and total share for each unique grid codes

# sum of GRID_AREA all unique GRIDCODE - SCH
regdist_sch %>% 
  distinct(GRIDCODE, .keep_all = TRUE) %>% 
  summarise(sum = sum(GRID_AREA))

# sum of BB_PART_AREA all unique GRIDCODE - SCH
total_area_BB <- 
  regdist_sch %>% 
  distinct(GRIDCODE, .keep_all = TRUE) %>% 
  summarise(sum = sum(BB_PART_AREA)) %>% 
  pull()

# sum of GRID_AREA all unique GRIDCODE - ALB
regdist_alb %>% 
  distinct(GRIDCODE, .keep_all = TRUE) %>% 
  summarise(sum = sum(GRID_AREA))

# sum of BW_PART_AREA all unique GRIDCODE - ALB

total_area_BW <- 
  regdist_alb %>% 
  distinct(GRIDCODE, .keep_all = TRUE) %>% 
  summarise(sum = sum(BW_PART_AREA)) %>% 
  pull()


# SCH

regdist_sch %>% 
  # set which statuses 
  filter(ZEITRAUM == 3) %>% 
  mutate(BB_PART_AREA = round(BB_PART_AREA, digits = 3)) %>%
  group_by(SCI_NAME) %>% 
  summarise(TOTAL_GRID_AREA = sum(GRID_AREA),
            TOTAL_SHARE_REG = sum(BB_PART_AREA),
            PROP_SHARE_REG = TOTAL_SHARE_REG/total_area_BB) %>%
  mutate(PROP_SHARE_REG = round(PROP_SHARE_REG, digits = 3)) %>% View()

# ALB

regdist_alb %>% 
  # set which statuses 
  filter(ZEITRAUM == 3) %>% 
  mutate(BW_PART_AREA = round(BW_PART_AREA, digits = 3)) %>%
  group_by(SCI_NAME) %>% 
  summarise(TOTAL_GRID_AREA = sum(GRID_AREA),
            TOTAL_SHARE_REG = sum(BW_PART_AREA),
            PROP_SHARE_REG = TOTAL_SHARE_REG/total_area_BW) %>% 
  mutate(PROP_SHARE_REG = round(PROP_SHARE_REG, digits = 3)) %>% View()

# check what's wrong with Rosa canina

regdist_alb %>% 
  filter(ZEITRAUM == 3 & SCI_NAME == "Rosa canina") %>% 
  select(GRID_AREA) %>% range()

regdist_alb %>% 
  filter(ZEITRAUM == 3 & SCI_NAME == "Rosa canina") %>% 
  select(BW_PART_AREA) %>% range()

regdist_alb %>% 
  filter(ZEITRAUM == 3 & SCI_NAME == "Rosa canina") %>% 
  count(SUCHNR)

regdist_alb %>% 
  filter(ZEITRAUM == 3 & SCI_NAME == "Rosa canina") %>% 
  count(GRIDCODE)

# It seems that this species has duplicated entries of GRIDCODE -> because in FloraWeb the SUCHNR of agg. doesn't exist, so I got duplicated values for R. canina. -> remove one, they are the same

regdist_alb %>% 
  filter(ZEITRAUM == 3) %>% 
  add_count(SCI_NAME, GRIDCODE) %>% 
  filter(n > 1) %>% 
  distinct(SUCHNR, SCI_NAME)

# species Muscari neglectum is under two different SUCHNR, Polygala comosa and Rosa canina have duplicated entries for each GRIDCODE. We should then sum only the area of each unique GRIDCODE per species 

sch <- 
regdist_sch %>% 
  # set which statuses 
  filter(ZEITRAUM == 3) %>% 
  distinct(SCI_NAME, GRIDCODE, .keep_all = TRUE) %>% 
  group_by(SCI_NAME) %>% 
  summarise(TOTAL_GRID_AREA = sum(GRID_AREA),
            TOTAL_SHARE_REG = sum(BB_PART_AREA),
            PROP_SHARE_REG = TOTAL_SHARE_REG/total_area_BB,
            REGION = "SCH") %>%
  mutate(PROP_SHARE_REG = round(PROP_SHARE_REG, digits = 3)) %>% 
  select(REGION, SCI_NAME, PROP_SHARE_REG)

# ALB

alb <- 
regdist_alb %>% 
  # set which statuses 
  filter(ZEITRAUM == 3) %>% 
  distinct(SCI_NAME, GRIDCODE, .keep_all = TRUE) %>% 
  group_by(SCI_NAME) %>% 
  summarise(TOTAL_GRID_AREA = sum(GRID_AREA),
            TOTAL_SHARE_REG = sum(BW_PART_AREA),
            PROP_SHARE_REG = TOTAL_SHARE_REG/total_area_BW,
            REGION = "ALB") %>%
  mutate(PROP_SHARE_REG = round(PROP_SHARE_REG, digits = 3)) %>% 
  select(REGION, SCI_NAME, PROP_SHARE_REG)

# join both datasets

reg_dist_plants <- bind_rows(alb, sch)

# fix names

reg_dist_plants <-
  reg_dist_plants %>% 
  mutate(SCI_NAME = SCI_NAME %>% 
           str_replace_all(., " ", "_") %>% 
           str_replace(., "_agg\\.$", "_aggr.")) 
    

# save data

# write_csv(reg_dist_plants, file = "data/raw/reg_dist_plants.csv")





# 15. ECONOMIC MODELLING -----------------------------------------------------

#' *Steps*: 
# - read each optimisation file
# - add column name with the name of the file (i.e. optimisation)
# - wrangle all files and merge
# - save each individual optimisation file merged with the all scores table 
# - save all files as one big table as well 
# - do the same for scenarios 

library(tidyverse)

## ---- read in and tidy data

# read big table with all scores and measures 

all_measures_alb <-
  read_delim("resources_troph-cost/economic/optimisations/EP_all_measures_scores_ALB_cost_indicators.csv") %>% 
  rename_with(tolower, .cols = everything()) %>% 
  relocate(region, .before = measure_id) %>% 
  relocate(epid, id, .after = region)



## 15.1 Unconstrained optimisations  ------------------------------------

# list of paths for unconstrained optimisations

paths <- file.path("resources_troph-cost/economic/optimisations",
                   list.files("resources_troph-cost/economic/optimisations",
                              pattern = "^Result_"))

# create a list with all optimisations

optis_list <-
  paths %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>%
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(optimisation =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>%
      separate_wider_regex(optimisation, 
                           c(budget = "^.{0,4}", "_", optimisation = ".*")) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(optis_list, names(optis_list), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations",
                      "/merged_unconstrained-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis <-
  optis_list %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  mutate(budget = if_else(budget == "n", NA_character_, budget) %>% 
           as.integer(budget)) %>% 
  filter(ecovalue != 0) %>% 
  mutate(budget = sum(cost),
         .by = optimisation) 

# save the combined long file 

write_csv(optis, "data/raw/optimisations.csv")


## 15.2 [OLD] Constrained optimisations for a 1458€ budget  -----------------

# OLD...archived!!

## ---- read each optimization file 

# list of paths for constrained optimisations for a 1458€ budget 
paths_1458 <- file.path("resources_troph-cost/economic/optimisations/1458",
                        list.files("resources_troph-cost/economic/optimisations/1458",
                                   pattern = "^Result_"))


optis_list_1458 <-
  paths_1458 %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>%
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(optimisation =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>%
      separate_wider_regex(optimisation, 
                           c(budget = "^.{0,4}", "_", optimisation = ".*")) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(optis_list_1458, names(optis_list_1458), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/1458",
                      "/merged-1458-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_1458 <- 
  optis_list_1458 %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  filter(ecovalue != 0)

# save the combined long file 

write_csv(optis_1458, "data/raw/optimisations_1458.csv")



## 15.3 Constrained optimisations for a 1344€ budget  ----------------------

## ---- read each optimization file 

# list of paths for constrained optimisations for a 1344€ budget 
paths_1344 <- file.path("resources_troph-cost/economic/optimisations/1344",
                        list.files("resources_troph-cost/economic/optimisations/1344",
                                   pattern = "^Result_"))


optis_list_1344 <-
  paths_1344 %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>%
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(optimisation =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>%
      separate_wider_regex(optimisation, 
                           c(budget = "^.{0,4}", "_", optimisation = ".*")) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(optis_list_1344, names(optis_list_1344), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/1344",
                      "/merged-1344-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_1344 <- 
  optis_list_1344 %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  filter(ecovalue != 0)

# save the combined long file 

write_csv(optis_1344, "data/raw/optimisations_1344.csv")





## 15.4 Constrained optimisations for a 991€ budget  ----------------------

## ---- read each optimization file 

# list of paths for constrained optimisations for a 991€ budget 
paths_991 <- file.path("resources_troph-cost/economic/optimisations/991",
                        list.files("resources_troph-cost/economic/optimisations/991",
                                   pattern = "^Result_"))


optis_list_991 <-
  paths_991 %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>%
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(optimisation =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>%
      separate_wider_regex(optimisation, 
                           c(budget = "^.{0,4}", "_", optimisation = ".*")) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(optis_list_991, names(optis_list_991), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/991",
                      "/merged-991-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_991 <- 
  optis_list_991 %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  filter(ecovalue != 0)

# save the combined long file 

write_csv(optis_991, "data/raw/optimisations_991.csv")




## 15.4 Unconstrained optimisations for scenarios ---------------------------

## ---- read each optimization file 

# list of paths for unconstrained optimisations

paths_scenarios <-
  file.path(
    "resources_troph-cost/economic/optimisations/scenarios",
    list.files(
      "resources_troph-cost/economic/optimisations/scenarios",
      pattern = "^Result_"
    )
  )

# create a list with all optimisations

scenarios_list <-
  paths_scenarios %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>% 
  
# loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(opti =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>% 
      mutate(optimisation = str_extract(opti, ".*[^_[:digit:]]"),
             budget = str_extract(opti, "[:digit:]{0,4}$"),
             .before = epid) %>%  
      select(-opti) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(scenarios_list, names(scenarios_list), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/scenarios",
                      "/merged_scenario_unconstrained-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_scenarios <-
  scenarios_list %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  mutate(budget = if_else(budget == "n", NA_character_, budget) %>% 
           as.integer(budget)) %>% 
  filter(ecovalue != 0) %>% 
  mutate(budget = sum(cost),
         .by = optimisation) 

# save the combined long file 

write_csv(optis_scenarios, "data/raw/optis_scenarios.csv")



## 15.5 Optimisations scenarios for a 991€ budget ----------------------------

## ---- read each optimization file 

# list of paths for unconstrained optimisations

paths_scenarios_991 <-
  file.path(
    "resources_troph-cost/economic/optimisations/scenarios_budget991",
    list.files(
      "resources_troph-cost/economic/optimisations/scenarios_budget991",
      pattern = "^Result_"
    )
  )

# create a list with all optimisations

scenarios_list_991 <-
  paths_scenarios_991 %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>% 
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(opti =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>% 
      mutate(optimisation = str_extract(opti, ".*[^_[:digit:]]"),
             budget = str_extract(opti, "[:digit:]{0,4}$"),
             .before = epid) %>%  
      select(-opti) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(scenarios_list_991, names(scenarios_list_991), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/scenarios_budget991",
                      "/merged_scenario_991-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_scenarios_991 <-
  scenarios_list_991 %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  mutate(budget = if_else(budget == "n", NA_character_, budget) %>% 
           as.integer(budget)) %>% 
  filter(ecovalue != 0) %>% 
  mutate(budget = sum(cost),
         .by = optimisation) 

# save the combined long file 

write_csv(optis_scenarios_991, "data/raw/optis_scenarios_991.csv")




## 15.6 Optimisations scenarios for a 1344 € budget --------------------------

## ---- read each optimization file 

# list of paths for unconstrained optimisations

paths_scenarios_1344 <-
  file.path(
    "resources_troph-cost/economic/optimisations/scenarios_budget1344",
    list.files(
      "resources_troph-cost/economic/optimisations/scenarios_budget1344",
      pattern = "^Result_"
    )
  )

# create a list with all optimisations

scenarios_list_1344 <-
  paths_scenarios_1344 %>%
  set_names(basename(.) %>%
              str_remove_all(., "Result_") %>%
              str_remove_all(., ".csv")) %>% 
  
  # loop over each file
  map(
    ~ read.csv(.x, sep = ";", dec = ",") %>%
      as_tibble() %>%
      select(-X) %>%
      
      # clean names
      rename_with(.fn = tolower, .cols = everything()) %>%
      rename_with(~ gsub("x.", "", .x, fixed = TRUE)) %>%
      rename_with(~ gsub(".", "", .x, fixed = TRUE)) %>%
      rename(measure_id = matches("MID"),
             epid = matches("FID")) %>%
      
      # add the file name as a col name and extract budget and indicator name
      mutate(opti =
               basename(.x) %>%
               str_remove(., "Result_") %>%
               str_remove(., ".csv"), 
             .before = epid) %>% 
      mutate(optimisation = str_extract(opti, ".*[^_[:digit:]]"),
             budget = str_extract(opti, "[:digit:]{0,4}$"),
             .before = epid) %>%  
      select(-opti) %>% 
      relocate(budget, .after = totalcost) %>% 
      
      # merge the files with all measures with each optimisation in the list
      inner_join(all_measures_alb, .,
                 by = c("epid", "measure_id")) %>% 
      rename(cost = totalcost)) 

## ---- save all individual merged files

walk2(scenarios_list_1344, names(scenarios_list_1344), ~ {
  file_name <- paste0("resources_troph-cost/economic/optimisations/scenarios_budget1344",
                      "/merged_scenario_1344-", .y, ".csv")
  write_csv(.x, file_name)
})


## ---- combine all files into one

optis_scenarios_1344 <-
  scenarios_list_1344 %>%
  list_rbind() %>%
  # add a unique identifier to each row
  unite(
    optimisation,
    epid,
    col = "id",
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) %>%
  relocate(id, .before = region) %>% 
  mutate(budget = if_else(budget == "n", NA_character_, budget) %>% 
           as.integer(budget)) %>% 
  filter(ecovalue != 0) %>% 
  mutate(budget = sum(cost),
         .by = optimisation) 

# save the combined long file 

write_csv(optis_scenarios_1344, "data/raw/optis_scenarios_1344.csv")
