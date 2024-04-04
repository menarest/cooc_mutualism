#' ----
#' title: Lepidoptera data analysis - matrix correlations
#' author: Esteban Menares
#' date: 30.06.2022
#' ----



# Set up  -------------------------------------------------------------

library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter, dplyr::select)
library(cooccur)
library(GGally)

source('scripts/source_script.R')

# read in data

pairs <- 
  read_csv('data/raw/pairs.csv')


# Filtered tables with species occurrences per site 

occur_alb <-
  read_csv('data/processed/occur_imago_ALB_PA.csv') %>%
  column_to_rownames("EPID") %>%
  as.matrix()

occur_sch <- 
  read_csv('data/processed/occur_imago_SCH_PA.csv') %>% 
  column_to_rownames("EPID") %>% 
  as.matrix()



# Trophic Interactions Matrices --------------------------------------------

### ALB heatmap ----

pairs %>%
  filter(region == "ALB") %>% 
  select(plant_sp, lepidoptera_sp, troph_interaction) %>%
  mutate(
    plant_sp = factor(plant_sp,
                      levels = sort(unique(plant_sp))),
    lepidoptera_sp = factor(lepidoptera_sp, 
                            levels = rev(sort(unique(lepidoptera_sp)))),
    troph_interaction = factor(as.character(troph_interaction),
                               levels = rev(sort(unique(troph_interaction))))
  ) %>% 
  ggplot(    
    aes(
      x = plant_sp,
      y = lepidoptera_sp,
      # use this if you want to reorder by interaction
      # x = reorder(plant_sp, as.numeric(troph_interaction)),
      # y = reorder(lepidoptera_sp, as.numeric(troph_interaction)),
      fill = troph_interaction)) +
  geom_tile(color = "black", size = 0.05) +
  # to keep tiles squares 
  coord_fixed() +
  scale_fill_brewer(palette = "YlGnBu", direction = -1) +
  # remove extra space
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  #set a base size for all fonts
  theme_grey(base_size = 5) +
  # rotate axis labels 90° counterclockwise and align them vertically at 
  # their end (hjust = 1) and their centers horizontally 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 6),
        plot.margin = unit(c(0, 1, 0, 1), "mm")) +
  labs(
    title = 'Trophic Interaction Matrix - ALB',
    x = "Plant",
    y = "Lepidoptera",
    fill = 'Interaction \nstrength')

# save plot 

ggsave(
  filename = "heatmap_troph_int_ALB.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 5.5,
  width = 7.5
) 

### SCH heatmap ----

pairs %>%
  filter(region == "SCH") %>%
  select(plant_sp, lepidoptera_sp, troph_interaction) %>%
  mutate(
    plant_sp = factor(plant_sp,
                      levels = sort(unique(plant_sp))),
    lepidoptera_sp = factor(lepidoptera_sp, 
                            levels = rev(sort(unique(lepidoptera_sp)))),
    troph_interaction = factor(as.character(troph_interaction),
                               levels = rev(sort(unique(troph_interaction))))
  ) %>%
  ggplot(    
    aes(
      x = plant_sp,
      y = lepidoptera_sp,
      # use this if you want to reorder by interaction
      # x = reorder(plant_sp, as.numeric(troph_interaction)),
      # y = reorder(lepidoptera_sp, as.numeric(troph_interaction)),
      fill = troph_interaction)) +
  geom_tile(color = "black", size = 0.05) +
  # to keep tiles squares 
  coord_fixed() +
  # choose a color-blind colour palette 
  scale_fill_brewer(palette = "YlGnBu", direction = -1) +
  #remove extra space
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  #set a base size for all fonts
  theme_grey(base_size = 5) +
  # rotate axis labels 90° counterclockwise and align them vertically at 
  # their end (hjust = 1) and their centers horizontally 
  theme(axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 6)) +
  labs(
    title = 'Trophic Interaction Matrix - SCH',
    x = "Plant",
    y = "Lepidoptera",
    fill = 'Interaction \nstrength')

# save plot 

ggsave(
  filename = "heatmap_troph_int_SCH.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 5,
  width = 7.5
) 


# Prediction Troph Int vs Co-occurrence --------------------------------

### ALB, troph int > 0, RII with p-value ≤ 0.2 -----

pairs %>% 
  filter(region == "ALB") %>%
  mutate(
    RII_0.2 = 
      case_when(
        troph_interaction > 0 & p_upper <= 0.2 ~ "TP",
        troph_interaction == 0 & p_upper > 0.2 & p_upper <= 1 ~ "TN",
        troph_interaction > 0 & p_upper > 0.2 & p_upper <= 1 ~ "FN",
        troph_interaction == 0 & p_upper <= 0.2 ~ "FP"
      )
  ) %>%
  mutate(
    RII_0.2 = 
      factor(as.character(RII_0.2)),
    lepidoptera_sp = factor(lepidoptera_sp, 
                            levels = rev(sort(unique(lepidoptera_sp))))
  ) %>% 
  ggplot(
    aes(
      x = plant_sp,
      y = lepidoptera_sp,
      fill = RII_0.2)) +
  geom_tile(color="white", size=0.1) +
  coord_fixed() +
  # choose a color-blind colour palette 
  scale_fill_brewer(palette = "Accent", direction = -1) +
  #remove extra space
  scale_y_discrete(expand=c(0, 0)) +
  scale_x_discrete(expand = c(0,0)) +
  #set a base size for all fonts
  theme_grey(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)) +
  labs(
    title = 'Predicted Associations vs True Interactions Matrix, ALB, RII with p ≤ 0.2',
    x = "Plants",
    y = "Lepidoptera",
    fill = 'Detection')

# save plot 

ggsave(
  filename = "heatmap_accuracy_species_ALB.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 5,
  width = 7.5
) 
  



### SCH, troph int > 0, RII with p-value ≤ 0.2 -----

pairs %>% 
  filter(region == "SCH") %>%
  mutate(
    p_upper = if_else(co_occurrence < 0, 1, p_upper),
    RII_0.2 = 
      case_when(
        troph_interaction > 0 & p_upper <= 0.2 ~ "TP",
        troph_interaction == 0 & p_upper > 0.2 & p_upper <= 1 ~ "TN",
        troph_interaction > 0 & p_upper > 0.2 & p_upper <= 1 ~ "FN",
        troph_interaction == 0 & p_upper <= 0.2 ~ "FP"
      )
  ) %>%
  mutate(
    RII_0.2 = 
      factor(as.character(RII_0.2)),
    lepidoptera_sp = factor(lepidoptera_sp, 
                            levels = rev(sort(unique(lepidoptera_sp))))
  ) %>% 
  ggplot(
    aes(
      x = plant_sp,
      y = lepidoptera_sp,
      fill = RII_0.2)) +
  geom_tile(color = "white", size = 0.1) +
  coord_fixed() +
  # choose a color-blind colour palette 
  scale_fill_brewer(palette = "Accent", direction = -1) +
  #remove extra space
  scale_y_discrete(expand=c(0, 0)) +
  scale_x_discrete(expand = c(0,0)) +
  #set a base size for all fonts
  theme_grey(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 12),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)) +
  labs(
    title = 'Predicted Associations vs True Interactions Matrix, SCH, RII with p ≤ 0.2',
    x = "Plants",
    y = "Lepidoptera",
    fill = 'Detection')

# save plot 

ggsave(
  filename = "heatmap_accuracy_species_SCH.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 5,
  width = 7.5
)



# Number of predictions ----------------------------------------------------

pairs %>% 
  # filter(region == "ALB") %>%
  select(region, plant_sp, lepidoptera_sp, troph_interaction, co_occurrence,
  ) %>% 
  mutate(
    RII_0.5 = 
      case_when(
        troph_interaction > 0 & co_occurrence >= 0.5 ~ 1,    # True positive 
        troph_interaction == 0 & co_occurrence < 0.5 ~ 0,    # True negative
        troph_interaction > 0 & co_occurrence < 0.5  ~ 0.75, # False negative
        troph_interaction == 0 & co_occurrence >= 0.5 ~ 0.25 # False positive
      )
  ) %>%
  # mutate(
  #   RII_0.5 =
  #     factor(as.character(RII_0.5))) %>%
  group_by(region) %>%
  summarise(
    TP = sum(RII_0.5 == 1),
    TN = sum(RII_0.5 == 0),
    FN = sum(RII_0.5 == 0.75),
    FP = sum(RII_0.5 == 0.25),
    TPpercent = sum(RII_0.5 == 1)/n()*100,
    Tnpercent = sum(RII_0.5 == 0)/n()*100,
    FNpercent = sum(RII_0.5 == 0.75)/n()*100,
    FPpercent = sum(RII_0.5 == 0.25)/n()*100
  )






# Package cooccur -----------------------------------------------------


### ALB ----------------------------------------------------------------

class(occur_alb)

# transpose matrix to have a spp by site matrix
occur_alb_t <- t(occur_alb)

# calculate co-occurrences with cooccur package
cooccur_ALB <- cooccur(mat = occur_alb_t,
                       type = "spp_site", # rows = spp, cols = sites
                       thresh = FALSE, # remove species pairs that are expected to have less than 1 co-occurrence
                       spp_names = TRUE, # if species names are in the rows or colums
                       true_rand_classifier = 0.1, # default, truly random associations that deviate from their exp by < 10%
                       only_effects = FALSE, # if TRUE, shows only effect sizes
                       eff_standard = TRUE, # if FALSE, effects = obs - exp. If TRUE, standardized effects
                       eff_matrix = FALSE) # if TRUE, effect sizes returned in a distance matrix
class(cooccur_ALB)

#' return an analysis-wide count of the number of species combinations 
#' classified as positive, negative, or random. 
summary(cooccur_ALB) # BUT this gives the overall including plant-plant, and
# lepi-lepi interactions, we need to calculate it manually. 


# RETURN A TABLE OF RESULTS
#' Return the entire table with all species pairs and their co-occurrence 
#' statistics. 
#' The table has values for “p_lt” and “p_gt” which are the probabilities 
#' that the two species co-occur with a frequency less than (lt) or 
#' greater than (gt) expected by chance. 
prob.table(cooccur_ALB)

# access and save the results from cooccur object 
all_cooccur <- cooccur_ALB$results

# output a pairwise probability table containing significant species combinations only
print(cooccur_ALB) 


# PLOT A VISUAL REPRESENTATIONS OF SPECIES CO-OCCURRENCES
#' Visually interpret these co-occurrence results.
#' Create a lower triangle heat map visually indicating the significant 
# positive, negative, and random co-occurrence patterns among all species 
plot(cooccur_ALB) # add "plotrand = TRUE" to include completely random species


# EXTRACT INFORMATION FOR A FOCAL SPECIES
#' extracts the significant positive and negative association data for a 
# single species
#pair(mod = cooccur_ALB, spp = "Geospiza scandens")


# EXTRACT INFORMATION FOR ALL SPECIES
#' summarize for each species the percent of its associations that 
#' are positive, negative, or random.
pair.attributes(cooccur_ALB)

# create a ranked (by percent significant associations) bar plot showing 
# the percentage of positive, negative, and random associations for each 
# species
# pair.profile(cooccur_ALB)
pair.profile.EM(cooccur_ALB)

# plot the observed versus expected number of co-occurrence sites for each 
# species pair
obs.v.exp(cooccur_ALB) 


### SCH -----------------------------------------------------------------


class(occur_sch)

# transpose matrix
occur_sch_transposed <- t(occur_sch)

# calculate co-occurrences with cooccur package
cooccur_SCH <- cooccur(mat = occur_sch_transposed,
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


# output a pairwise probability table containing significant species combinations only
print(cooccur_SCH) 


# PLOT A VISUAL REPRESENTATIONS OF SPECIES CO-OCCURRENCES
#' Visually interpret these co-occurrence results.
#' Create a lower triangle heat map visually indicating the significant positive, 
#' negative, and random co-occurrence patterns among all species 
plot(cooccur_SCH) # add "plotrand = TRUE" to include completely random species


# EXTRACT INFORMATION FOR A FOCAL SPECIES
#' extracts the significant positive and negative association data for a single species
#pair(mod = cooccur_SCH, spp = "Geospiza scandens")


# EXTRACT INFORMATION FOR ALL SPECIES
#' summarize for each species the percent of its associations that 
#' are positive, negative, or random.
pair.attributes(cooccur_SCH)

# create a ranked (by percent significant associations) bar plot showing the percentage 
# of positive, negative, and random associations for each species
# pair.profile(cooccur_SCH)
pair.profile.EM(cooccur_SCH)

# plot the observed versus expected number of co-occurrence sites for each species pair
obs.v.exp(cooccur_SCH) 
