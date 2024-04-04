
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
  readxl::read_xlsx("resources/12526_Lepidoptera-2008-OpenData/12526.xlsx",
                    sheet = "12526", 
                    na = "NA")

# Plant cover table: In spring of the years 2008-2020, we sampled all species in an area of 4m x 4m and estimated the percentage cover of each species relative to the whole 4 m x 4 m plot. The vegetation records are done in early summer, normally in May and beginning of June. NA = Missing data, because already mown or not accessible because of grazing animals

# BExIS ID # 23586 (1.3.1)

plants <- 
  readxl::read_xlsx("resources/23586_Vegetation_BExIS-OpenData/23586.xlsx",
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





# 3. TROPHIC INTERACTIONS -------------------------------------------------

# NOTE: this data were extracted from Ebert (2005) for ALB, and Richert and Brauner (2018) for SCH, then collected in a single excel file which then was here cleaned up. Cleaned versions with a single file per region for larva and adults were stored in BExIS as long format and can be found using the ID: 

# Imago ALB: 31734
# Larva ALB: 31735
# Imago SCH: 31736
# Larva SCH: 31737

# we then use them in ~troph-cost/scripts/wrangling/pairs_1_trophic_links.R in wide format to perform the analyses of co-occurrence and other posterior analyses. 

## 3.1 Larva_BB (Berlin Brandenburg = BB = SCH) 

Larva_BB <- 
  readxl::read_xlsx("resources/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Larva_BB", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.2 Larva_BW (Baden-Württemberg = BW = ALB) 
Larva_BW <- 
  readxl::read_xlsx("resources/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Larva_BW", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.3 Imago_BB  (Berlin Brandenburg = BB = SCH) 

Imago_BB <- 
  readxl::read_xlsx("resources/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Imago_BB", na = "NA") %>% 
  # rename butterfly names 
  rename(Leptidea_sinapis_complex = Leptidea_sinapis_reali,
         Melitaea_athalia = Melitaea_athalia_Komplex,
         Pyrgus_alveus_complex = Pyrgus_alveus_Komplex)

## 3.4 Imago_BW (Baden-Württemberg = BW = ALB) 

Imago_BW <- 
  readxl::read_xlsx("resources/trophic-interactions/trophic-links_matrix_lepidoptera_plants.xlsx", sheet = "Imago_BW", na = "NA") %>% 
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







# 5. FLOWER AVAILABILITY ---------------------------------------------------

#' **Aim:**
#' per species and EP, calculate average of flowering units across all 
#' dates and turn dataset into long format dataset: obs at each site for
#' each plant and butterfly species.


## ---- Load data

# BExIS dataset 4981 (v2) "Flower availability 2008 Alb-korrigiert"
flowers_ALB <- readxl::read_xlsx("resources/4981_2_FlowerAvailability_ALB-2008-OpenData/4981_2_data.xlsx")
str(flowers_ALB)

# BExIS dataset 4964 (v2) "Flower availability 2008 Schorfheide"
flowers_SCH <- readxl::read_xlsx("resources/4964_2_FlowerAvailability_SCH-2008-OpenData/4964_2_data.xlsx")
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
  read_delim(file = "resources/10160_FlowerVisitorInteractions-2008-OpenData/10160_3_data.csv",
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

# round non-integer numbers to 3 decimal places

pairs <- pairs %>% 
  mutate(
    across(c(sites_plant, sites_lepi, obs_cooccur, flower_visit), 
    as.integer)) %>% 
  mutate(
    across(
      where(is.double), ~round(.x, digits = 3))) 

## ---- save data 

write_excel_csv(pairs, 'data/raw/pairs.csv', na = "NA")