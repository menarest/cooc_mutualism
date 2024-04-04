#' ----
#' title: Lepidoptera data processing - species correlations
#' author: Esteban Menares
#' date: 26.11.2021
#' ----

#' Project: TrophCost - co-occurrence analysis
#' Processing step No: 3/4

  
#' **Overall workflow data processing for Lepidoptera data:**
#' 1. Data processing 1 - trophic-links
#' 2. Data processing 2 - co-occurrences 
#' 3. Data processing 3 - species correlation matrices
#' 4. Data processing 4 - flower visitations (data validation)
 
#' **Description**
#' We want to have an estimate of the relationship between plants and 
#' butterflies abundances per site. Are certain butterfly species more
#' correlated with certain plants? To answer this we run, in addition to 
#' the co-occurrence analysis (P/A data), a species correlation analysis. 
#' We can explore this by correlating the abundance of butterflies to the 
#' abundance of plants, but we can also estimate it more precisely using 
#' the abundance of flowers or flower availability. 

#' Flower availability is estimated by determining the number of 'flower 
#' units' across the grassland types of the Exploratories as a basis for 
#' analyzing flower-visitor interactions networks.


#' **General Aim** 
#' Add correlations between plants and imago butterfly species into the 
#' species-level dataset per site using data for plants abundances and 
#' flower abundances (availability).
 

#' **Workflow** 
#' For each butterfly species:

#' **Calculate species correlations with plant abundances => sp_corr**
#'      - load data and check data 
#'      - check NAs
#'      - filter plants that are not present in trophic-links dataset
#'      - check NAs
#'      - make a pairwise correlation with the BE list abundances
#'      - filter plants and butterflies not present in trophic-links dataset
#'      - check and remove cols of species without occurrences 
#'      - do analysis for each individual region
#'      - per pairs, correlation between abundances across all plots
#'      - use Spearman correlation (non-parametric) Rho (r) 
#'      - calculate effect size

#' **Calculate species correlations with flower availability => flower_corr** 
#' - Step 1: we take the butterfly abundances from the corr_imago_region 
#'        datasets, which are already checked and evaluated for: 
#'      - check NAs
#'      - subset plants that do not occur in both
#'      - separate per EPID
#'      - remove cols of species without occurrences 
#' - Step 2: we take the flower abundances from the flower_region dataset
#'      - load and check data
#'      - log + 1 transform to avoid high range fo values
#'      - transform from long to wide format (sites by species)
#'      - check NAs
#'      - filter out plants not present in trophic-links dataset
#'      - separate per EPID
#'      - remove cols of species without occurrences 
#' - Step 3: we combine both of them into one single dataset
#'      - combine both previous datasets
#'      - do analysis for each individual region


#' **cols to add...**
#' plant_cor:           Spearman's rank correlation coefficient for the abundance of butterflies and plants species
#' plant_cor_p:         p-value for the Spearman's rank correlation of butterfly and plant abundances
#' plant_cor_eff_size:  effect size (r^2) of the Spearman's rank correlation of butterfly and plant abundances
#' flower_cor:          Spearman's rank correlation coefficients for the abundance of butterfly and flowers (log(x + 1) transformed)
#' flower_cor_p:        p-value for the spearman correlation of butterflies and flowers abundances
#' flower_cor_eff_size: Effect sizes of the correlation of butterflies and flower abundances (r^2)



# Set up --------------------------------------------------------------

#' ---- Requirements
library(tidyverse)  # ggplots and data wrangling 
library(data.table) # loading and wrangling data
library(cowplot)    # formatting ggplots 
library(Hmisc)      # correlations coefficients and p-values



# 1. Reading data ---------------------------------------------------------

# long format dataset with all possible pairs of plants-lepidoptera ALB
imago_ALB <- fread(file = "data/raw/imago_ALB.csv") 
dim(imago_ALB)
str(imago_ALB)

# long format dataset with all possible pairs of plants-lepidoptera SCH
imago_SCH <- fread(file = "data/raw/imago_SCH.csv") 
dim(imago_SCH)
str(imago_SCH)

# wide format dataset with all occurrences per plot of plants and butterflies
corr_all <- fread(file = "data/raw/occur.csv") 
str(corr_all)

# long format dataset with flower units per plant species and plot
flowers_ALB <- fread(file = "data/raw/flower_units_ALB.csv")
flowers_SCH <- fread(file = "data/raw/flower_units_SCH.csv")

# NOTES: 
# - the datasets imago_SCH and imago_ALB were outputted in ~/troph-cost/scripts/wrangling/pairs_2_co-occurrence.R. 
# - the dataset corr_all was cleaned up and prepared in section 1) in ~/troph-cost/scripts/wrangling/data_tidying.R. 
# - the datasets flowers_ALB and flowers_SCH were cleaned up and prepared in section 5) in ~/troph-cost/scripts/wrangling/data_tidying.R. 




# 2. Plant abundance data -------------------------------------------------

## 2.1 Data tidying and wrangling -----------------------------------------

## ---- remove NAs
corr_all <- na.omit(corr_all)



## ---- remove sites of region HAI to match trophic-link matrix
dim(corr_all)
corr_all <- corr_all[!grepl("HEG", EPID), ] 
dim(corr_all)



## ---- per region, save unique plants and lepi names from imago datasets  
# to use them for sub-setting occurrence/abundance data

# imago_SCH
names_lepi_imago_SCH <- unique(imago_SCH$lepidoptera_sp) # names imago spp 
names_plant_imago_SCH <- unique(imago_SCH$plant_sp)      # names plants spp

# imago_ALB
names_lepi_imago_ALB <- unique(imago_ALB$lepidoptera_sp) # names imago spp 
names_plant_imago_ALB <- unique(imago_ALB$plant_sp)      # names plants spp




## ---- per region, subset only plants and butterflies (cols) that are 
# present in trophic-links

# imago_SCH
corr_imago_SCH <- 
  corr_all %>%
  select(EPID, 
         all_of(c(names_lepi_imago_SCH, names_plant_imago_SCH))) 
str(corr_imago_SCH)

# imago_ALB
corr_imago_ALB <- 
  corr_all %>%
  select(EPID, 
         all_of(c(names_lepi_imago_ALB, names_plant_imago_ALB))) 
str(corr_imago_ALB)





## ---- make EPID as row.names so that they are not treated as data
# imago_SCH
class(corr_imago_SCH)
corr_imago_SCH <- as.data.frame(corr_imago_SCH) # set up as data frame, cooccur needs data frame
class(corr_imago_SCH)
row.names(corr_imago_SCH) <- corr_imago_SCH$EPID # make EPIDs as row names
corr_imago_SCH$EPID <- NULL # remove EPID column

# imago_ALB
class(corr_imago_ALB)
corr_imago_ALB <- as.data.frame(corr_imago_ALB) # set up as data frame, cooccur needs data frame
class(corr_imago_ALB)
row.names(corr_imago_ALB) <- corr_imago_ALB$EPID # make EPIDs as row names
corr_imago_ALB$EPID <- NULL # remove EPID column


# Check that the dimensions of cooccur_imago datasets are equal to the 
# sum of unique names of plants + lepidoptera in main (imago) datasets
length(c(names_plant_imago_SCH, 
         names_lepi_imago_SCH)) == dim(corr_imago_SCH)[2]
length(c(names_plant_imago_ALB, 
         names_lepi_imago_ALB)) == dim(corr_imago_ALB)[2]

# remove unnecessary objects from environment 
rm(names_lepi_imago_SCH, names_lepi_imago_ALB,
   names_plant_imago_SCH, names_plant_imago_ALB)




## ---- make integers as numeric
corr_imago_SCH <- 
  corr_imago_SCH %>% 
  mutate(
    across(
      where(is.integer), as.numeric))
corr_imago_ALB <- 
  corr_imago_ALB %>% 
  mutate(
    across(
      where(is.integer), as.numeric))
str(corr_imago_SCH)
str(corr_imago_ALB)




## ---- subset obs (rows) per region
# imago_ALB

dim(corr_imago_ALB)
corr_imago_ALB <- 
  corr_imago_ALB[grepl("AEG", row.names(corr_imago_ALB)), ]
dim(corr_imago_ALB)

# imago_SCH

dim(corr_imago_SCH)
corr_imago_SCH <- 
  corr_imago_SCH[grepl("SEG", row.names(corr_imago_SCH)), ]
dim(corr_imago_SCH)



## ---- remove cols of species without occurrences (colSum = 0)

# Region ALB 
corr_imago_ALB_filtered <- 
  corr_imago_ALB %>%
  select(
    where(~ sum(.x) !=0))

# Region SCH 
corr_imago_SCH_filtered <- 
  corr_imago_SCH %>%
  select(
    where(~ sum(.x) !=0))



## ---- remove those species without occurrences from trophic-links (imago)

# first save the names of butterflies and plants species
cols_corr_imago_ALB_filtered <- 
  corr_imago_ALB_filtered %>%
  select(
    everything()) %>%
  names()

cols_corr_imago_SCH_filtered <- 
  corr_imago_SCH_filtered %>%
  select(
    everything()) %>%
  names()

# then use those names to remove species from imago that don't occur 
imago_ALB_filtered <- imago_ALB[
  plant_sp %in% cols_corr_imago_ALB_filtered, ] # common plants

imago_ALB_filtered <- imago_ALB_filtered[
  lepidoptera_sp %in% cols_corr_imago_ALB_filtered, ] # common lepidoptera

imago_SCH_filtered <- imago_SCH[
  plant_sp %in% cols_corr_imago_SCH_filtered, ] # common plants

imago_SCH_filtered <- imago_SCH_filtered[
  lepidoptera_sp %in% cols_corr_imago_SCH_filtered, ] # common lepidoptera

# then make a vector with unique names to check and compare dims
names_lepi_imago_ALB <- unique(imago_ALB_filtered$lepidoptera_sp)
names_plant_imago_ALB <- unique(imago_ALB_filtered$plant_sp)

names_lepi_imago_SCH <- unique(imago_SCH_filtered$lepidoptera_sp)
names_plant_imago_SCH <- unique(imago_SCH_filtered$plant_sp)

# Check: the dimensions of corr_imago_ALB_filtered are equal to the sum
# of unique names of plants + Lepidoptera of imago_ALB_filtered
length(
  c(names_lepi_imago_ALB, names_plant_imago_ALB)) == 
  dim(corr_imago_ALB_filtered)[2] # dim(...)[2] to select cols

# Check: the same for corr_imago_SCH_filtered and imago_SCH_filtered
length(
  c(names_lepi_imago_SCH, names_plant_imago_SCH)) == 
  dim(corr_imago_SCH_filtered)[2]



## ---- remove unnecessary objects from environment
rm(corr_imago_SCH, corr_imago_ALB, imago_ALB, imago_SCH, 
   names_lepi_imago_ALB, names_lepi_imago_SCH, 
   names_plant_imago_ALB, names_plant_imago_SCH)





## 2.2 Run analysis --------------------------------------------------------

## ---- Create function for formatting a correlation matrix into a table 
# with 4 columns containing:

# Column 1 : row names (variable 1 for the correlation test)
# Column 2 : column names (variable 2 for the correlation test)
# Column 3 : the correlation coefficients
# Column 4 : the p-values of the correlations

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values

flattenCorrMatrix <- 
  function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    lepidoptera_sp = rownames(cormat)[row(cormat)[ut]],
    plant_sp = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# more info: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software

# The output of the function rcorr() is a list containing the following 
# elements : - r : the correlation matrix - n : the matrix of the number of 
# observations used in analyzing each pair of variables - P : the p-values 
# corresponding to the significance levels of correlations.

### ALB --------------------------------------------------------------------

# create and extract correlation coefficients and p-values
cor_mat_plant_ALB <- 
  rcorr(
    as.matrix(corr_imago_ALB_filtered), 
    type = "spearman")

cor_df_plant_ALB <- 
  flattenCorrMatrix(
    cor_mat_plant_ALB$r, 
    cor_mat_plant_ALB$P)

# rename cols 
colnames(cor_df_plant_ALB)[which(
  names(cor_df_plant_ALB) == "cor")] <- "plant_cor"

colnames(cor_df_plant_ALB)[which(
  names(cor_df_plant_ALB) == "p")] <- "plant_cor_p"

# load the universal matrix
imago_ALB <- 
  read_csv(file = "data/raw/imago_ALB.csv")

# combine datasets: add correlations into main dataset
imago_ALB <- 
  left_join(imago_ALB, 
            cor_df_plant_ALB, 
            by = c("plant_sp", 
                   "lepidoptera_sp"))

# add effect sizes: r2 or coefficient of determination
imago_ALB$plant_cor_eff_size <- (imago_ALB$plant_cor)^2



### SCH -------------------------------------------------------------------

# create and extract correlation coefficients and p-values
cor_mat_plant_SCH <- rcorr(as.matrix(corr_imago_SCH_filtered), 
                           type = "spearman")
cor_df_plant_SCH <- flattenCorrMatrix(cor_mat_plant_SCH$r, 
                                      cor_mat_plant_SCH$P)

# rename cols 
colnames(cor_df_plant_SCH)[which(
  names(cor_df_plant_SCH) == "cor")] <- "plant_cor"

colnames(cor_df_plant_SCH)[which(
  names(cor_df_plant_SCH) == "p")] <- "plant_cor_p"

# load the universal matrix
imago_SCH <- fread(file = "data/raw/imago_SCH.csv")

# combine datasets: add correlations into main dataset
imago_SCH <- dplyr::left_join(imago_SCH, 
                              cor_df_plant_SCH, 
                              by = c("plant_sp", 
                                     "lepidoptera_sp"))

# add effect sizes: r2 or coefficient of determination
imago_SCH$plant_cor_eff_size <- (imago_SCH$plant_cor)^2




# 3. Flower abundance (availability) data -------------------------------

## 3.1 Data wrangling: butterflies --------------------------------------

## ---- subset, per region, only cols with butterflies 

# imago_ALB
butterflies_imago_ALB <- 
  corr_imago_ALB_filtered %>%
  select(Aphantopus_hyperantus:Vanessa_cardui)
str(butterflies_imago_ALB)

# imago_SCH
butterflies_imago_SCH <- 
  corr_imago_SCH_filtered %>%
  select(Anthocharis_cardamines:Thymelicus_sylvestris)
str(butterflies_imago_SCH)



## 3.1 Data wrangling: flowers ------------------------------------------

## ---- Check data
str(flowers_ALB)
str(flowers_SCH)

range(flowers_ALB$flower_units)
range(flowers_SCH$flower_units)

# Log + 1 transform values to account for high range of values
flowers_ALB$log_flower_units <- log(flowers_ALB$flower_units + 1)
flowers_SCH$log_flower_units <- log(flowers_SCH$flower_units + 1)

# Check data again 
range(flowers_ALB$log_flower_units)
range(flowers_SCH$log_flower_units)



## ---- Long to wide format: sites per plant species (rows, cols)
flowers_ALB <- 
  flowers_ALB %>% 
  select(-flower_units) %>% 
  pivot_wider(names_from = plant_sp, 
              values_from = log_flower_units)

flowers_SCH <- 
  flowers_SCH %>% 
  select(-flower_units) %>% 
  pivot_wider(names_from = plant_sp, 
              values_from = log_flower_units)



## ---- Remove NAs
flowers_ALB <- na.omit(flowers_ALB)
flowers_SCH <- na.omit(flowers_SCH)



## ---- Subset, per region, plants not present in trophic-links datasets

# first make a subset of unique distinct names from imago_region_filtered
plants_imago_ALB <- unique(imago_ALB_filtered$plant_sp)
plants_imago_SCH <- unique(imago_SCH_filtered$plant_sp)


# second use names to filter species that are not present in imago_region
# imago_ALB
flowers_ALB_filtered <- 
  flowers_ALB %>%
  select(EPID, 
         any_of(plants_imago_ALB))
str(flowers_ALB_filtered)

# imago_SCH
flowers_SCH_filtered <- 
  flowers_SCH %>%
  select(EPID, 
         any_of(plants_imago_SCH))
str(flowers_SCH_filtered)



## ---- Make EP as row.names so that they are not treated as data

# ALB
class(flowers_ALB_filtered)
flowers_ALB_filtered <- as.data.frame(flowers_ALB_filtered)
class(flowers_ALB_filtered)
row.names(flowers_ALB_filtered) <- flowers_ALB$EPID # make EPIDs as row names
flowers_ALB_filtered$EPID <- NULL # remove EPID column

# SCH
class(flowers_SCH_filtered)
flowers_SCH_filtered <- as.data.frame(flowers_SCH_filtered)
class(flowers_SCH_filtered)
row.names(flowers_SCH_filtered) <- flowers_SCH$EPID # make EPIDs as row names
flowers_SCH_filtered$EPID <- NULL # remove EPID column


## ---- Remove cols of species without occurrences (colSum = 0)

# Check if any col sums zero 
all(colSums(flowers_ALB_filtered) != 0)
all(colSums(flowers_SCH_filtered) != 0)


#' **NOTE** we DO NOT remove species without occurrences from trophic-links
#' dataset, but we will code them later as NA in the main dataset 
#' when we match both datasets per species names



## ---- Add butterfly abundances and filter EPIDs that exist in both datasets

# first save vectors with the EPID from flowers and imago datasets
EPID_flowers_ALB <- row.names(flowers_ALB_filtered)
EPID_butterflies_ALB <- row.names(butterflies_imago_ALB)

EPID_flowers_SCH <- row.names(flowers_SCH_filtered)
EPID_butterflies_SCH <- row.names(butterflies_imago_SCH)


# check which EPID that appear in x but not in y: setdiff(x, y)
setdiff(EPID_flowers_ALB, EPID_butterflies_ALB)
setdiff(EPID_butterflies_ALB, EPID_flowers_ALB)

setdiff(EPID_flowers_SCH, EPID_butterflies_SCH)
setdiff(EPID_butterflies_SCH, EPID_flowers_SCH)


# join the two datasets 
# we keep only EPID with a match in flower_region and butterflies_region

#ALB
lepi_flowers_ALB <- 
  merge(butterflies_imago_ALB, 
        flowers_ALB_filtered, 
        by = "row.names", 
        copy = FALSE, 
        keep = FALSE, 
        na_matches = "na")

# SCH
lepi_flowers_SCH <- 
  merge(butterflies_imago_SCH, 
        flowers_SCH_filtered, 
        by = "row.names", 
        copy = FALSE, 
        keep = FALSE, 
        na_matches = "na")


## ---- Prepare data frames to run matrix operations

row.names(lepi_flowers_ALB) <- lepi_flowers_ALB$Row.names # make EPIDs as row names
lepi_flowers_ALB$Row.names <- NULL # remove EPID column

row.names(lepi_flowers_SCH) <- lepi_flowers_SCH$Row.names # make EPIDs as row names
lepi_flowers_SCH$Row.names <- NULL # remove EPID column



## 3.2 Run analysis -------------------------------------------------------

### ALB -------------------------------------------------------------------

# create and extract correlation coefficients and p-values
cor_mat_flowers_ALB <- rcorr(as.matrix(lepi_flowers_ALB), 
                             type = "spearman")
cor_df_flowers_ALB <- flattenCorrMatrix(cor_mat_flowers_ALB$r, 
                                        cor_mat_flowers_ALB$P)

# rename cols 
colnames(cor_df_flowers_ALB)[which(
  names(cor_df_flowers_ALB) == "cor")] <- "flower_cor"

colnames(cor_df_flowers_ALB)[which(
  names(cor_df_flowers_ALB) == "p")] <- "flower_cor_p"

# combine datasets: add correlations into main dataset
imago_ALB <- dplyr::left_join(imago_ALB, 
                              cor_df_flowers_ALB, 
                              by = c("plant_sp", 
                                     "lepidoptera_sp"))

# add effect sizes: r2 or coefficient of determination
imago_ALB$flower_cor_eff_size <- (imago_ALB$flower_cor)^2





### SCH -------------------------------------------------------------------

# create and extract correlation coefficients and p-values
cor_mat_flowers_SCH <- rcorr(as.matrix(lepi_flowers_SCH), 
                             type = "spearman")
cor_df_flowers_SCH <- flattenCorrMatrix(cor_mat_flowers_SCH$r, 
                                        cor_mat_flowers_SCH$P)

# rename cols 
colnames(cor_df_flowers_SCH)[which(
  names(cor_df_flowers_SCH) == "cor")] <- "flower_cor"

colnames(cor_df_flowers_SCH)[which(
  names(cor_df_flowers_SCH) == "p")] <- "flower_cor_p"

# combine datasets: add correlations into main dataset
imago_SCH <- dplyr::left_join(imago_SCH, 
                              cor_df_flowers_SCH, 
                              by = c("plant_sp", 
                                     "lepidoptera_sp"))

# add effect sizes: r2 or coefficient of determination
imago_SCH$flower_cor_eff_size <- (imago_SCH$flower_cor)^2



# Correct for multiple testing --------------------------------------------

imago_ALB <- imago_ALB %>% 
  mutate(
    # apply multiple testing correction (Benjamini & Hochberg's FDR)
    plant_cor_p = p.adjust(plant_cor_p, method = "fdr"),
    flower_cor_p = p.adjust(flower_cor_p, method = "fdr"))

imago_SCH <- imago_SCH %>% 
  mutate(
    # apply multiple testing correction (Benjamini & Hochberg's FDR)
    plant_cor_p = p.adjust(plant_cor_p, method = "fdr"),
    flower_cor_p = p.adjust(flower_cor_p, method = "fdr"))



# 6. Save data ------------------------------------------------------------

# write_csv(imago_ALB, file = "data/raw/imago_ALB.csv", na = 'NA')
# write_csv(imago_SCH, file = "data/raw/imago_SCH.csv", na = 'NA')

## END OF SCRIPT ##