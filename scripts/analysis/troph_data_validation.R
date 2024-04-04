#' ----
#' title: Lepidoptera data analysis - data validation
#' author: Esteban Menares
#' date: 13.02.2023
#' ----


#' **AIM**: 
#' Compare the pairs that do have visits and check whether those pairs have significant negative co-occurrences. For that:
#'   - check if each pair has negative RII along with the p_lt that I obtained with the cooccur package.
#'   - check the standardised effect sizes (= (obs - exp)/No. of sites) and p_lt that I also obtained from the cooccur package.
#'   - check plant and flower neg correlations using pairs with p-value ≤ 0.2
#'   - correlation between visits and RII, probabilistic, and Spearman. 
#'   - correlation between visits and interactions


# Set up --------------------------------------------------------------

library(tidymodels)
library(tidyverse)
conflicted::conflicts_prefer(dplyr::select(), dplyr::filter())
library(sjPlot)
library(GGally) # add-on to ggplot 
library(ggpubr) # stacking ggplots 
library(MASS) # function glm.nb 
library(AER) # for a formal test of overdispersion
library(DHARMa) # model checks
library(AICcmodavg) # model selection and overdispersion (c hat) test

theme_set(theme_classic(6))

source("scripts/source_script.R")

# load data

pairs <- read_csv("data/raw/pairs.csv")


# Data exploration and analysis -------------------------------------------

# save pairs df 

pairs_orig <- pairs
str(pairs)

# create a subset with only obs where flower_visit is not NA 

pairs <-
  pairs %>% 
  filter(!is.na(flower_visit)) %>%
  select(region,
         troph_pair_id,
         troph_interaction,
         rii = co_occurrence,
         p_lower,
         p_upper,
         cooccur_stdeff = std_effects,
         p_lt,
         p_gt,
         plant_cor,
         plant_cor_p,
         flower_cor, 
         flower_cor_p,
         flower_visit)


# 1. Check negative co-occurrences -------------------------------------------

pairs %>% 
  
  # reclassify variables with (significant) negative co-occurrences
  mutate(
    neg_rii = if_else(rii < 0 & p_lower <= 0.2, 1, 0),
    neg_prob = if_else(cooccur_stdeff < 0 & p_lt <= 0.2, 1, 0),
    neg_plant_cor = if_else(plant_cor < 0 & plant_cor_p <= 0.2, 1, 0),
    neg_flower_cor = if_else(flower_cor < 0 & flower_cor_p <= 0.2, 1, 0)) %>% 
  
  # count the number of negative co-occurrences and the percentage 
  summarise(
    n = n(),
    n_rii = sum(neg_rii, na.rm = TRUE),
    n_prob = sum(neg_prob, na.rm = TRUE),
    n_plant_cor = sum(neg_plant_cor, na.rm = TRUE),
    n_flower_cor = sum(neg_flower_cor, na.rm = TRUE),
    perc_rii = n_rii/n()*100,
    perc_prob = n_prob/n()*100,
    perc_plant_cor = n_plant_cor/n()*100,
    perc_flower_cor = n_flower_cor/n()*100,
    .by = region
  ) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 1))) %>% 
  write.table(., "output/cooccur_validation/neg_cooccur_methods.txt",
              sep = ",",
              quote = FALSE,
              row.names = FALSE)



# 2. Check correlations visits ~ co-occurrence associations -----------------

# create models for each association method - ALB

modList <- list()

## ---- ALB

# plant cor
modList[["Plant cor - ALB"]] <- lm_plant_cor_alb <- 
  lm(plant_cor ~ flower_visit, data = pairs,
     subset = 
       pairs$plant_cor_p <= 0.2 & pairs$plant_cor > 0 & pairs$region == "ALB")

# flower cor
modList[["Flower cor - ALB"]] <- lm_flower_cor_alb <- 
  lm(flower_cor ~ flower_visit, data = pairs,
     subset = 
       pairs$flower_cor_p <= 0.2 & pairs$flower_cor > 0 & pairs$region == "ALB")

# probabilistic
modList[["Probabilistic - ALB"]] <- lm_prob_alb <- 
  lm(cooccur_stdeff ~ flower_visit, data = pairs,
     subset = 
       pairs$p_gt <= 0.2 & pairs$cooccur_stdeff > 0 & pairs$region == "ALB")

# RII
modList[["RII - ALB"]] <- lm_rii_alb <- 
  lm(rii ~ flower_visit, data = pairs,
     subset = pairs$p_upper <= 0.2 & pairs$rii > 0 & pairs$region == "ALB")

## ---- SCH

# plant cor
modList[["Plant cor - SCH"]] <- lm_plant_cor_sch <- 
  lm(plant_cor ~ flower_visit, data = pairs,
     subset = 
       pairs$plant_cor_p <= 0.2 & pairs$plant_cor > 0 & pairs$region == "SCH")
# Concl: not enough degrees of freedom, undetermined model. Not assigned! remove

# flower cor
modList[["Flower cor - SCH"]] <- lm_flower_cor_sch <- 
  lm(flower_cor ~ flower_visit, data = pairs,
     subset = 
       pairs$flower_cor_p <= 0.2 & pairs$flower_cor > 0 & pairs$region == "SCH")
# Concl: not enough degrees of freedom, undetermined model. Assigned but remove!

modList$`Flower cor - SCH` <- NULL

# probabilistic
modList[["Probabilistic - SCH"]] <- lm_prob_sch <- 
  lm(cooccur_stdeff ~ flower_visit, data = pairs,
     subset = 
       pairs$p_gt <= 0.2 & pairs$cooccur_stdeff > 0 & pairs$region == "SCH")

# RII
modList[["RII - SCH"]] <- lm_rii_sch <- 
  lm(rii ~ flower_visit, data = pairs,
     subset = pairs$p_upper <= 0.2 & pairs$rii > 0 & pairs$region == "SCH")

# export results into word tables

modList %>% names()

tab_model(
  modList[1],
  show.se = TRUE,
  dv.labels = c("Plant correlation - ALB"), 
  pred.labels = c("Intercept",
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_plant_cor.doc"
)

tab_model(
  modList[2],
  show.se = TRUE,
  dv.labels = c("Flower correlation - ALB"), 
  pred.labels = c("Intercept", 
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_flower_cor.doc"
)

tab_model(
  modList[c(3,5)],
  show.se = TRUE,
  dv.labels = c("Probabilistic - ALB",
                "Probabilistic - SCH"), 
  pred.labels = c("Intercept", 
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_probabilistic.doc"
)

tab_model(
  modList[c(4,6)],
  show.se = TRUE,
  dv.labels = c("RII - ALB",
                "RII - SCH"),
  pred.labels = c("Intercept", 
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_rii.doc"
)

# or per region summary
modList %>% names()

tab_model(
  modList[c(1:4)],
  show.se = FALSE,
  collapse.ci = TRUE,
  dv.labels = c("Plant cor - ALB",
                "Flower cor - ALB",    
                "Probabilistic - ALB",
                "RII - ALB"),
  pred.labels = c("Intercept", 
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_cooccur_methods_alb.doc"
)

tab_model(
  modList[c(5:6)],
  show.se = FALSE,
  collapse.ci = TRUE,
  dv.labels = c("Probabilistic - SCH",
                "RII - SCH"),
  pred.labels = c("Intercept", 
                  "Flower visits"),
  file = "output/cooccur_validation/flower_visits_vs_cooccur_methods_sch.doc"
)

# generate a ggplot for flower visits with (significant) positive associations 

pairs %>%
  ggplot() +
  facet_wrap("region", 
             scales = "fixed") +
  
  # add plant_cor 
  geom_point(data = pairs %>% 
               filter(plant_cor_p <= 0.2, plant_cor > 0 & region != "SCH"),
             aes(x = flower_visit, y = plant_cor, col = "Plant cor")) +
  geom_smooth(data = pairs %>% 
                filter(plant_cor_p <= 0.2, plant_cor > 0  & region != "SCH"),
              method = "lm",
              aes(x = flower_visit, y = plant_cor, col = "Plant cor"), 
              se = FALSE) + 
  
  # add flower_cor
  geom_point(data = pairs %>% 
               filter(flower_cor_p <= 0.2, flower_cor > 0 & region != "SCH"),
             aes(x = flower_visit, y = flower_cor, 
                 col = "Flower cor")) +
  geom_smooth(data = pairs %>% 
                filter(flower_cor_p <= 0.2, flower_cor > 0 & region != "SCH"),
              method = "lm",
              aes(x = flower_visit, y = flower_cor, 
                  col = "Flower cor"),
              se = FALSE) + 
  
  # add probabilistic 
  geom_point(data = pairs %>%
               filter(p_gt <= 0.2, cooccur_stdeff > 0),
             aes(x = flower_visit, y = cooccur_stdeff, col = "Probabilistic")) +
  geom_smooth(data = pairs %>%
                filter(p_gt <= 0.2, cooccur_stdeff > 0),
              method = "lm",
              aes(x = flower_visit, y = cooccur_stdeff, col = "Probabilistic"),
              se = FALSE) + 
  
  # add RII
  geom_point(data = pairs %>%
               filter(p_upper <= 0.2, rii > 0),
             aes(x = flower_visit, y = rii, col = "RII")) +
  geom_smooth(data = pairs %>%
                filter(p_upper <= 0.2, rii > 0),
              method = "lm",
              aes(x = flower_visit, y = rii, col = "RII"),
              se = FALSE) + 
  
  # add scales 
  scale_color_manual(name = "Association method: ",
                     values = c("Plant cor" = "seagreen",
                                "Flower cor" = "slateblue",
                                "Probabilistic" = "gold",
                                "RII" = "skyblue")) + 
  labs(x = "Flower visit count", y = "Association strenght") +
  theme(legend.position = "bottom") 

# save ggplot

ggsave(
  filename = "correlation-flower-visits-vs-assoc-methods.pdf",
  plot = last_plot(),
  path = "output/cooccur_validation/",
  dpi = "retina",
  height = 10,
  width = 10,
  units = "cm"
)




# 3. Check correlation: visits ~ ecological interactions --------------------

# ALB
cor_alb <- cor.test(
  data = pairs,
  ~ flower_visit + troph_interaction,
  alternative = "greater", # corresponds to positive association
  method = "spearman", # for non-normal values i.e. rank-based association
  exact = FALSE, # p-values are computed using asymptotic t approximation
  subset = pairs$region == "ALB"
)

# SCH
cor_sch <- cor.test(
  data = pairs,
  ~ flower_visit + troph_interaction,
  alternative = "greater", # corresponds to positive association
  method = "spearman", # for non-normal values i.e. rank-based association
  exact = FALSE, # p-values are computed using asymptotic t approximation
  subset = pairs$region == "SCH"
)

# save results

cor_alb %>% 
  tidy() %>% 
  bind_rows(cor_sch %>% tidy()) %>% 
  mutate(region = c("ALB", "SCH"),
         .before = estimate) %>%   
  mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>% 
  write.table(., 'output/cooccur_validation/visits_vs_troph_int.txt',
              sep = ",",
              quote = FALSE,
              row.names = FALSE)
