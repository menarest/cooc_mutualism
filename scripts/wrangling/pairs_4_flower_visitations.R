#' ----
#' title: Lepidoptera data processing - species correlations
#' author: Esteban Menares
#' date: 14.02.2022
#' ----

#' Project: TrophCost - co-occurrence analysis
#' Processing step No: 4/4

  
#' **Overall workflow data processing for Lepidoptera data**
#' 1. Data processing 1 - trophic-links
#' 2. Data processing 2 - co-occurrences 
#' 3. Data processing 3 - species correlation matrices
#' 4. Data processing 4 - flower visitations (data validation)
 

#' **Description**
#' We use the flower-visitor interaction dataset from BExIS (ID 10160)
#' for data validation. This data give us information of interactions 
#' between plants and butterflies. The level of information provided 
#' by this data is higher than co-occurrences but lower than trophic-
#' interactions in terms of sampling effort, geographical extension and considers only one year instead of multiple years.

#' Data collected between May - August 2008, 50 grassland EP per exploratory, 
#' sampling different land use types simultaneously. Some plots sampled 
#' repeatedly over the year. All flowering plant species (an independent 
#' estimate of the total number of flower units per species) and all 
#' flower visitor-plant interactions were be recorded on a 600 m2 transect
#'  around each plot sampled for six hours. A ‘flower unit’ is defined 
#'  as a unit of one or more flowers demanding insects to fly from one 
#'  unit to another. Insects were caught with a dip net or with help of 
#'  an exhauster. Only insects sitting directly in the centre of the 
#'  flowers seeming to feed on pollen or nectar were caught. Insects 
#'  sitting on the petals were not taken into consideration.

# Weiner, Christiane; Werner, Michael; Linsenmair, Karl Eduard; Blüthgen, Nico (2016): flower visitor interactions 2008. Version 3. Biodiversity Exploratories Information System. Dataset. https://www.bexis.uni-jena.deddm/data/Showdata/10160?version=3

#'   Plant species and number of interactions between butterflies and 
#'   plant species: each cell shows the number of individuals of a 
#'   butterfly species on a certain flower species.


#' **General Aim** 
#' Add the flower visitation data into the main dataset. Each value 
#' represent the sum of individuals of a butterfly species on a given 
#' flower species across all intercepts. 
 

#' **cols to add...**
#' flower_visit:        Flower visitation record (data validation)



# 1. Setup -----------------------------------------------------------------

# Requirements

library(tidyverse)  # ggplots and data wrangling 


# Reading in data 

# long format dataset with all possible pairs of plants-lepidoptera ALB

imago_ALB <- 
  read_csv(file = "data/raw/imago_ALB.csv")
dim(imago_ALB)
str(imago_ALB)

# long format dataset with all possible pairs of plants-lepidoptera SCH

imago_SCH <- 
  read_csv(file = "data/raw/imago_SCH.csv") 
str(imago_SCH)

# long format dataset with butterfly's flower visitations

visit_all <- 
  read_csv(file = "data/raw/flower_visitation.csv")
str(visit_all)


#' NOTES: 
#' - the datasets imago_SCH and imago_ALB were outputted in ~/troph-cost/scripts/wrangling/pairs_3_species_correlations.R
#' - the dataset visit_all was cleaned up and prepared in section 6) in ~/troph-cost/scripts/wrangling/data_tidying.R. 




# 2. Data wrangling -----------------------------------------------------

# get the total sum of butterflies visits per region and plant species across all EPs, dates and intercepts 

# region ALB

visit_ALB <-
  visit_all %>% 
  # Group data 
  group_by(Region, plant_sp, lepidoptera_sp) %>%
  # calculate mean flower visitation
  summarise(flower_visit = sum(Total),
            .groups = 'drop') %>% 
  arrange(plant_sp, lepidoptera_sp) %>%
  mutate(flower_visit = 
           round(flower_visit, 2)) %>% 
  filter(Region == 'ALB')

# region SCH

visit_SCH <-
  visit_all %>%
  # Group data
  group_by(Region, plant_sp, lepidoptera_sp) %>%
  # calculate mean flower visitation
  summarise(flower_visit = sum(Total),
            .groups = 'drop') %>%
  arrange(plant_sp, lepidoptera_sp) %>%
  mutate(flower_visit =
           round(flower_visit, 2)) %>%
  filter(Region == 'SCH')

# Check if each row is a unique observation

visit_ALB %>% 
  distinct()

visit_ALB %>% 
  select(plant_sp, lepidoptera_sp) %>% 
  distinct()

visit_SCH %>% 
  distinct()

visit_SCH %>% 
  select(plant_sp, lepidoptera_sp) %>% 
  distinct()

# yes!

# Join datasets with main dataset

# ALB

imago_ALB_new <- 
  imago_ALB %>% 
  left_join(
    visit_ALB %>% 
      select(plant_sp,
             lepidoptera_sp,
             flower_visit),
    by = c('plant_sp', 
           'lepidoptera_sp')) %>% 
  mutate(
    across(
      c(co_occurrence:p_upper, prob_cooccur:flower_cor_eff_size), 
      ~ round(., 5))) 

# SCH

imago_SCH_new <-
  imago_SCH %>% 
  left_join(
    visit_SCH %>% 
      select(plant_sp,
             lepidoptera_sp,
             flower_visit),
    by = c('plant_sp', 
           'lepidoptera_sp')) %>% 
  mutate(
    across(
      c(co_occurrence:p_upper, prob_cooccur:flower_cor_eff_size), 
      ~ round(., 5))) 



# 5. Save data ------------------------------------------------------------
 
# write_csv(imago_ALB_new, file = "data/raw/imago_ALB.csv")
# write_csv(imago_SCH_new, file = "data/raw/imago_SCH.csv")

#' **NOTE**: continue in ~/scripts/wrangling/data_tidying.R in point
#' 9. PAIRS to tidy both datasets into one 



## END OF SCRIPT ##