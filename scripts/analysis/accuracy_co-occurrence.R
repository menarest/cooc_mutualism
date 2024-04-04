#' ----
#' title: accuracy analysis for co-occurrence results 
#' author: Esteban Menares
#' date: 30.06.2022
#' ----

#' *AIM*: 
#' evaluate accuracy of co-occurrence methods and compare their results

# Based on: 

# Freilich et al. 2018 https://doi.org/10.1002/ecy.2142
# Goberna & Verdú 2022 https://doi.org/10.1016/j.soilbio.2021.108534
# Brisson et al. 2019 https://doi.org/10.3389/fmicb.2019.00585


# We will analyse the sensitivity and specificity ratio 
# of the co-occurrence algorithm results obtained with P/A and abundances 
# using different cutoffs:

# P/A: RII ≥ 0.5, 0.8, 1
# P/A: p_gt ≤ 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1
# Abundances: p-value ≤ 0.01, 0.05, 0.1, 0,2, 0,5, 0,8, 1 

# For this we first calculate the following:

# TP: number of links that are correctly identified as ecological interactions
# TN: number of links that are correctly identified non-interactions
# FP: number of links that are incorrectly identified as interactions
# FN: number of links that are incorrectly identified as non-interactions

# Sensitivity = TP / TP + FN i.e Probability of detection of a true link

# Specificity = TN / TN + FP i.e. True negative rate

# since plant-lepidoptera interaction is a mutualistic interaction the sign 
# of the relationship should only be positive. Therefore, we need to first 
# filter the negative interactions found by the co-occurrence method or just
# focus on the positive ones. So...

#' *STEPS*

#' 1. Keep only positive associations values from each method
#'   - for cooccur pckg results, use p_gt and values
#'  - for rii method, filter out negative co-occurrence values 
#'  - for Spearman method, filter out negative r values and use p-values
#' 2. per method, summarise data of TP, TN, FP, FN by cutoff 
#' 3. evaluate cutoff results and choose one cuttof to work with
#' 4. with selected cutoff, compare all methods with sensibility and specificity
#' 5. run a chi2 test on the sensibility values per method


# Set up  -------------------------------------------------------------

library(tidyverse)
library(ggpubr) # for extra ggplot functionalities
library(scales) # for extra ggplot functionalities
theme_set(theme_classic(7))
source('scripts/source_script.R')

# read in data

pairs <- read_csv('data/raw/pairs.csv')


# Rename variables  -------------------------------------------------------

# No. ("real") ecological trophic links -----------------------------------

pairs %>% 
  group_by(region) %>% 
  summarise(
    potential_links = 
      length(unique(plant_sp)) * length(unique(lepidoptera_sp)),
    realized_links = sum(troph_interaction > 0),
    connectance = potential_links/realized_links)


# Summary cooccur package plant-lepi interactions -----------------------

pairs %>% 
  group_by(region) %>% 
  summarise(
    n_plant_sp = n_distinct(plant_sp),
    n_lepi_sp = n_distinct(lepidoptera_sp),
    n_spp = n_plant_sp + n_lepi_sp,
    n_positive = sum(p_gt <= 0.05),
    n_negative = sum(p_lt <= 0.05), # negative associations are not of interest
    n_random = n() - (n_positive + n_negative),
    non_random = ((n_positive + n_negative) / n() * 100),
    total_pairs = n()
  )


# 1) P/A based results --------------------------------------------------

## 1.1) p_gt (cooccur) --------------------------------------------------

# explore number of negative and positive associations

pairs %>% 
  group_by(region) %>% 
  summarise(
    gt_0.01 = sum(p_gt <= 0.01),
    lt_0.01 = sum(p_lt <= 0.01),
    gt_0.05 = sum(p_gt <= 0.05),
    lt_0.05 = sum(p_lt <= 0.05),
    gt_0.1 = sum(p_gt <= 0.1),
    lt_0.1 = sum(p_lt <= 0.1),
    gt_0.2 = sum(p_gt <= 0.2),
    lt_0.2 = sum(p_lt <= 0.2),
    p_0.01 = sum(p_gt <= 0.01 & p_lt <= 0.01),
    p_0.05 = sum(p_gt <= 0.05 & p_lt <= 0.05),
    p_0.1 = sum(p_gt <= 0.1 & p_lt <= 0.1),
    p_0.2 = sum(p_gt <= 0.2 & p_lt <= 0.2),
    p_0.5 = sum(p_gt <= 0.5 & p_lt <= 0.5),
    p_1 = sum(p_gt <= 1 & p_lt <= 1)
  ) 

# no cases where a single pair has at the same time low p_gt and low p_lt, 
# i.e. a negative and positive associations are mutually exclusive

# now lets classify detections as True positives, True negative, False negative or False positive

pairs_p_gt <-
  pairs %>% 
  select(region, 
         troph_pair_id, 
         troph_interaction,
         p_gt,
         p_lt,
         std_effects
  ) %>% 
  mutate(
    p_gt = if_else(std_effects < 0, 2, p_gt),
    p_gt_0.05 =
      case_when(
        troph_interaction > 0 & p_gt <= 0.05 ~ "TP",
        troph_interaction == 0 & p_gt > 0.05 & p_gt <= 1 ~ "TN",
        # troph_interaction > 0 & p_gt > 0.05 & p_gt <= 1 ~ "FN",
        # troph_interaction == 0 & p_gt <= 0.05 ~ "FP",
        TRUE ~ "other"
      ),
    p_gt_0.1 = 
      case_when(
        troph_interaction > 0 & p_gt <= 0.1 ~ "TP",
        troph_interaction == 0 & p_gt > 0.1 & p_gt <= 1 ~ "TN",
        # troph_interaction > 0 & p_gt > 0.1 & p_gt <= 1 ~ "FN",
        # troph_interaction == 0 & p_gt <= 0.1 ~ "FP",
        TRUE ~ "other"
      ),
    p_gt_0.2 = 
      case_when(
        troph_interaction > 0 & p_gt <= 0.2 ~ "TP",
        troph_interaction == 0 & p_gt > 0.2 & p_gt <= 1 ~ "TN",
      #   troph_interaction > 0 & p_gt > 0.2 & p_gt <= 1 ~ "FN",
      #   troph_interaction == 0 & p_gt <= 0.2 ~ "FP",
        TRUE ~ "other"
      ),
    p_gt_0.5 = 
      case_when(
        troph_interaction > 0 & p_gt <= 0.5 ~ "TP",
        troph_interaction == 0 & p_gt > 0.5 & p_gt <= 1 ~ "TN",
        # troph_interaction > 0 & p_gt > 0.5 & p_gt <= 1 ~ "FN",
        # troph_interaction == 0 & p_gt <= 0.5 ~ "FP",
        TRUE ~ "other"
      ),
    p_gt_0.05 = factor(p_gt_0.05),
    p_gt_0.1 = factor(p_gt_0.1),
    p_gt_0.2 = factor(p_gt_0.2),
    p_gt_0.5 = factor(p_gt_0.5)
  ) %>% 
  pivot_longer(
    cols = p_gt_0.05:p_gt_0.5,
    names_to = "cutoff",
    values_to = "detected"
  )

### summary table ----

tab_pairs_p_gt <-
  pairs_p_gt %>%
  group_by(region, cutoff) %>%
  summarise(
    # Calculate n of known trophic links
    known_link = sum(troph_interaction > 0),
    # Calculate n of known non-trophic links
    known_nonlink = sum(troph_interaction == 0),
    # Calculate n of detected links per cutoff
    TP = sum(detected == "TP"),
    TN = sum(detected == "TN"),
    FP = known_nonlink - TN,
    FN = known_link - TP,
    # FP = sum(detected == "FP"),
    # FN = sum(detected == "FN"),
    # other = sum(detected == "other"),
    .groups = "drop"
  ) %>%
  mutate(sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) 

tab_pairs_p_gt

### counts per detection ----

# create stacked bar plot with counts per detection type (FN, FP, TN, TP)

p_p_gt <- 
  tab_pairs_p_gt %>% 
  pivot_longer(
    cols = TP:FN,
    names_to = "detected",
    values_to = "n"
  ) %>% 
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x = cutoff, 
             y = n,
             fill = detected)) +
  geom_col() + 
  facet_wrap("region") + 
  # increase the number of breaks
  scale_y_continuous(n.breaks = 5) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1) +
  labs(
    x = "Probabilistic method p-value",
    y = "Number of pairs") + 
  theme(legend.title = element_blank(),
        legend.position = "top")
p_p_gt

# save ggplot

ggsave(
  filename = "cutoff_p_gt.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)

### network sensitivity ----

# remember... 
# Sensitivity = TP / TP + FN i.e Probability of detection of a true link

# because the co-occurrence method included also negative associations, we will recalculate FN and sensitivity using only positive associations.

p_net_sensitivity_p_gt <-
  tab_pairs_p_gt %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FN,TP),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FN" = "Known interaction",
                               "TP" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level probabilistic method",
    y = "Count of links") + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_sensitivity_p_gt

# save ggplot

ggsave(
  filename = "sensitivity_p_gt.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5, 
  units = "cm"
)



### network specificity ----

# remember... Specificity = TN / TN + FP i.e. True negative rate

# because the co-occurrence method included also negative associations, 
# we will recalculate FP and specificity using only positive associations.

p_net_specificity_p_gt <-
  tab_pairs_p_gt %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FP,TN),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FP" = "Known non-interactions",
                               "TN" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level probabilistic method",
    y = "Count of non-links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_specificity_p_gt

# save ggplot
 
ggsave(
  filename = "specificity_p_gt.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)


## 1.2) [UPDATE] RII + Null Models  ------------------------------------------

pairs_rii <- pairs %>% 
  select(region, 
         troph_pair_id, 
         troph_interaction,
         co_occurrence,
         p_lower,
         p_upper
  ) %>% 
  mutate(
    p_upper = if_else(co_occurrence < 0, 2, p_upper),
    RII_0.05 =
      case_when(
        troph_interaction > 0 & p_upper <= 0.05 ~ "TP",
        troph_interaction == 0 & p_upper > 0.05 & p_upper <= 1 ~ "TN",
        # troph_interaction > 0 & p_upper > 0.05 & p_upper <= 1 ~ "FN",
        # troph_interaction == 0 & p_upper <= 0.05 ~ "FP",
        TRUE ~ "other"
      ),
    RII_0.1 = 
      case_when(
        troph_interaction > 0 & p_upper <= 0.1 ~ "TP",
        troph_interaction == 0 & p_upper > 0.1 & p_upper <= 1 ~ "TN",
        # troph_interaction > 0 & p_upper > 0.1 & p_upper <= 1 ~ "FN",
        # troph_interaction == 0 & p_upper <= 0.1 ~ "FP",
        TRUE ~ "other"
      ),
    RII_0.2 = 
      case_when(
        troph_interaction > 0 & p_upper <= 0.2 ~ "TP",
        troph_interaction == 0 & p_upper > 0.2 & p_upper <= 1 ~ "TN",
        # troph_interaction > 0 & p_upper > 0.2 & p_upper <= 1 ~ "FN",
        # troph_interaction == 0 & p_upper <= 0.2 ~ "FP",
        TRUE ~ "other"
      ),
    RII_0.5 = 
      case_when(
        troph_interaction > 0 & p_upper <= 0.5 ~ "TP",
        troph_interaction == 0 & p_upper > 0.5 & p_upper <= 1 ~ "TN",
        # troph_interaction > 0 & p_upper > 0.5 & p_upper <= 1 ~ "FN",
        # troph_interaction == 0 & p_upper <= 0.5 ~ "FP",
        TRUE ~ "other"
      ),
    RII_0.05 = factor(RII_0.05),
    RII_0.1 = factor(RII_0.1),
    RII_0.2 = factor(RII_0.2),
    RII_0.5 = factor(RII_0.5)
  ) %>% 
  pivot_longer(
    cols = RII_0.05:RII_0.5,
    names_to = "cutoff",
    values_to = "detected"
  ) %>% 
  mutate(
    cutoff = fct_relevel(cutoff, c("RII_0.05", "RII_0.1", "RII_0.2", "RII_0.5")))


### summary table ----

tab_pairs_rii <-
  pairs_rii %>%
  group_by(region, cutoff) %>%
  summarise(
    # Calculate n of known trophic links
    known_link = sum(troph_interaction > 0),
    # Calculate n of known non-trophic links
    known_nonlink = sum(troph_interaction == 0),
    # Calculate n of detected links per cutoff
    TP = sum(detected == "TP"),
    TN = sum(detected == "TN"),
    FP = known_nonlink - TN,
    FN = known_link - TP,
    # FP = sum(detected == "FP"),
    # FN = sum(detected == "FN"),
    # other = sum(detected == "other"),
    .groups = "drop"
  ) %>%
  mutate(sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) 

tab_pairs_rii

### counts per detection ----

# create stacked bar plot with counts per detection type (FN, FP, TN, TP)

p_rii <-
  tab_pairs_rii %>% 
  pivot_longer(
    cols = TP:FN,
    names_to = "detected",
    values_to = "n"
  ) %>% 
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x = cutoff, 
             y = n,
             fill = detected)) +
  geom_col() + 
  facet_wrap("region") + 
  # increase the number of breaks
  scale_y_continuous(n.breaks = 5) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1) +
  labs(
    x = "RII pairwise null models alpha",
    y = "Number of pairs",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top")
p_rii

# save ggplot

ggsave(
  filename = "cutoff_rii_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)

### network sensitivity ----

# because the co-occurrence method included also negative associations, 
# we will recalculate FN and sensitivity using only positive associations.

p_net_sensitivity_rii <-
  tab_pairs_rii %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FN,TP),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FN" = "Known interaction",
                               "TP" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level RII pairwise null models",
    y = "Count of links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_sensitivity_rii

# save ggplot

ggsave(
  filename = "sensitivity_rii_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)


### network specificity ----

# remember... Specificity = TN / TN + FP i.e. True negative rate

# because the co-occurrence method included also negative associations, 
# we will recalculate FN and specificity using only positive associations.

p_net_specificity_rii <- 
  tab_pairs_rii %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FP,TN),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FP" = "Known non-interactions",
                               "TN" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level RII pairwise null models",
    y = "Count of non-links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_specificity_rii

# save ggplot

ggsave(
  filename = "specificity_rii_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)





# 2) Abundance-based results ------------------------------------

## 2.1) (corrected) plant_p_cor (Spearman) ----

pairs_p_cor <-
  pairs %>% 
  select(region, 
         troph_pair_id, 
         troph_interaction,
         plant_cor,
         plant_cor_p
  ) %>% 
  mutate(
    # make negative associations non-significant 
    plant_cor_p = if_else(plant_cor < 0, 2, plant_cor_p),
    p_0.05 =
      case_when(
        troph_interaction > 0 & plant_cor_p <= 0.05 ~ "TP",
        troph_interaction == 0 & plant_cor_p > 0.05 & plant_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & plant_cor_p > 0.05 & plant_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & plant_cor_p <= 0.05 ~ "FP",
        TRUE ~ "other"
      ),  
    p_0.1 =
      case_when(
        troph_interaction > 0 & plant_cor_p <= 0.1 ~ "TP",
        troph_interaction == 0 & plant_cor_p > 0.1 & plant_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & plant_cor_p > 0.1 & plant_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & plant_cor_p <= 0.1 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.2 =
      case_when(
        troph_interaction > 0 & plant_cor_p <= 0.2 ~ "TP",
        troph_interaction == 0 & plant_cor_p > 0.2 & plant_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & plant_cor_p > 0.2 & plant_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & plant_cor_p <= 0.2 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.5 =
      case_when(
        troph_interaction > 0 & plant_cor_p <= 0.5 ~ "TP",
        troph_interaction == 0 & plant_cor_p > 0.5 & plant_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & plant_cor_p > 0.5 & plant_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & plant_cor_p <= 0.5 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.05 = factor(p_0.05),
    p_0.1 = factor(p_0.1),
    p_0.2 = factor(p_0.2),
    p_0.5 = factor(p_0.5)
  ) %>% 
  pivot_longer(
    cols = p_0.05:p_0.5,
    names_to = "cutoff",
    values_to = "detected"
  )

### summary table ----

tab_pairs_p_cor <-
  pairs_p_cor %>%
  group_by(region, cutoff) %>%
  summarise(
    # Calculate n of known trophic links
    known_link = sum(troph_interaction > 0),
    # Calculate n of known non-trophic links
    known_nonlink = sum(troph_interaction == 0),
    # Calculate n of detected links per cutoff
    TP = sum(detected == "TP"),
    TN = sum(detected == "TN"),
    FP = known_nonlink - TN,
    FN = known_link - TP,
    # FP = sum(detected == "FP"),
    # FN = sum(detected == "FN"),
    # other = sum(detected == "other"),
    .groups = "drop"
  ) %>%
  mutate(sensitivity = TP / (TP + FN),
         specificity = TN / (TN + FP)) 

tab_pairs_p_cor

### counts per detection ---- 

# create stacked bar plot with counts per detection type (FN, FP, TN, TP)

p_p_cor <-
  tab_pairs_p_cor %>% 
  pivot_longer(
    cols = TP:FN,
    names_to = "detected",
    values_to = "n"
  ) %>% 
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x = cutoff, 
             y = n,
             fill = detected)) +
  geom_col() + 
  facet_wrap("region") + 
  # increase the number of breaks
  scale_y_continuous(n.breaks = 5) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1) +
  labs(
    x = "Plant Spearman correlation alpha (corrected)",
    y = "Number of pairs",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top")
p_p_cor

# save ggplot

ggsave(
  filename = "cutoff_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)

### network sensitivity ----

# because the co-occurrence method included also negative associations, 
# we will recalculate FN and sensitivity using only positive associations.

p_net_sensitivity_p_cor <-
  tab_pairs_p_cor %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FN,TP),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FN" = "Known interaction",
                               "TP" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level Plant correlation (corrected)",
    y = "Count of links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_sensitivity_p_cor

# save ggplot

ggsave(
  filename = "sensitivity_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)


### network specificity ----

# remember... Specificity = TN / TN + FP i.e. True negative rate

p_net_specificity_p_cor <-
  tab_pairs_p_cor %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FP,TN),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FP" = "Known non-interactions",
                               "TN" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level Plant correlation (corrected)",
    y = "Count of non-links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_specificity_p_cor

# save ggplot

ggsave(
  filename = "specificity_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)





## 2.2) (corrected) flower_p_cor (Spearman) ------------------------------

# we did not have data of flowering availability (estimated by flower units) 
# for all species pairs, hence we should first filter out those pairs and
# recalculate known and not-known interactions

pairs %>% 
  select(region, 
         troph_pair_id, 
         troph_interaction,
         flower_cor,
         flower_cor_p
  ) %>% 
  filter(!is.na(flower_cor)) %>% 
  group_by(region) %>% 
  summarise(known_link = sum(troph_interaction > 0))

# calculate TP, TN... for different cutoff levels using flower abundance data

flower_p_cor <-
  pairs %>% 
  select(region, 
         troph_pair_id, 
         troph_interaction,
         flower_cor,
         flower_cor_p
  ) %>% 
  filter(!is.na(flower_cor)) %>% 
  mutate(
    # make negative associations non-significant 
    flower_cor_p = if_else(flower_cor < 0, 2, flower_cor_p),
    p_0.05 =
      case_when(
        troph_interaction > 0 & flower_cor_p <= 0.05 ~ "TP",
        troph_interaction == 0 & flower_cor_p > 0.05 & flower_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & flower_cor_p > 0.05 & flower_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & flower_cor_p <= 0.05 ~ "FP",
        TRUE ~ "other"
      ),  
    p_0.1 =
      case_when(
        troph_interaction > 0 & flower_cor_p <= 0.1 ~ "TP",
        troph_interaction == 0 & flower_cor_p > 0.1 & flower_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & flower_cor_p > 0.1 & flower_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & flower_cor_p <= 0.1 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.2 =
      case_when(
        troph_interaction > 0 & flower_cor_p <= 0.2 ~ "TP",
        troph_interaction == 0 & flower_cor_p > 0.2 & flower_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & flower_cor_p > 0.2 & flower_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & flower_cor_p <= 0.2 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.5 =
      case_when(
        troph_interaction > 0 & flower_cor_p <= 0.5 ~ "TP",
        troph_interaction == 0 & flower_cor_p > 0.5 & flower_cor_p <= 1 ~ "TN",
        # troph_interaction > 0 & flower_cor_p > 0.5 & flower_cor_p <= 1 ~ "FN",
        # troph_interaction == 0 & flower_cor_p <= 0.5 ~ "FP",
        TRUE ~ "other"
      ),
    p_0.05 = factor(p_0.05),
    p_0.1 = factor(p_0.1),
    p_0.2 = factor(p_0.2),
    p_0.5 = factor(p_0.5)
  ) %>% 
  pivot_longer(
    cols = p_0.05:p_0.5,
    names_to = "cutoff",
    values_to = "detected"
  )

### summary table ----

tab_flower_p_cor <-
  flower_p_cor %>%
  group_by(region, cutoff) %>%
  summarise(
    # Calculate n of known trophic links
    known_link = sum(troph_interaction > 0),
    # Calculate n of known non-trophic links
    known_nonlink = sum(troph_interaction == 0),
    # Calculate n of detected links per cutoff
    TP = sum(detected == "TP"),
    TN = sum(detected == "TN"),
    FP = known_nonlink - TN,
    FN = known_link - TP,
    # FP = sum(detected == "FP"),
    # FN = sum(detected == "FN"),
    # other = sum(detected == "other"),
    .groups = "drop"
  ) %>%
  mutate(sensitivity = TP / known_link,
         specificity = TN / known_nonlink) 

tab_flower_p_cor

### counts per detection ---- 

# create stacked bar plot with counts per detection type (FN, FP, TN, TP)

p_flower_p_cor <-
  tab_flower_p_cor %>% 
  pivot_longer(
    cols = TP:FN,
    names_to = "detected",
    values_to = "n"
  ) %>% 
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>% 
  ggplot(aes(x = cutoff, 
             y = n,
             fill = detected)) +
  geom_col() + 
  facet_wrap("region") + 
  # increase the number of breaks
  scale_y_continuous(n.breaks = 5) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1) +
  labs(
    x = "Flower Spearman correlation alpha (corrected)",
    y = "Number of pairs",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top")


# save ggplot

ggsave(
  filename = "cutoff_flower_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)

### network sensitivity ----

# because the co-occurrence method included also negative associations, 
# we will recalculate FN and sensitivity using only positive associations.

p_net_sensitivity_flower_p_cor <-
  tab_flower_p_cor %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FN,TP),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FN" = "Known interaction",
                               "TP" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level Flower correlation (corrected)",
    y = "Count of links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_sensitivity_flower_p_cor


# save ggplot

ggsave(
  filename = "sensitivity_flower_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)


### network specificity ----

# remember... Specificity = TN / TN + FP i.e. True negative rate

p_net_specificity_flower_p_cor <- 
  tab_flower_p_cor %>%
  mutate(cutoff = str_remove(cutoff, "^[^0-9]*")) %>% 
  pivot_longer(
    cols = c(FP,TN),
    names_to = "detected",
    values_to = "n"
  ) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = cutoff,
      y = n,
      fill = detected)) +
  geom_col() + 
  # facet per region
  facet_wrap("region") + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FP" = "Known non-interactions",
                               "TN" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), ")")),
            position = position_stack(vjust = 0.5), 
            size = 1.8) +
  labs(
    x = "Significance level Flower correlation (corrected)",
    y = "Count of non-links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )
p_net_specificity_flower_p_cor


# save ggplot

ggsave(
  filename = "specificity_flower_p_cor.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)





# 3) Comparing all methods -----------------------------------------------


###  Counts per detection type -------------------------------------------


ggarrange(
  p_p_cor +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )),
  p_p_gt  +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )),
  p_rii  +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )),
  labels = c("(a)", "(b)", "(c)"),
  hjust = -0.1,
  ncol = 3,
  nrow = 1,
  common.legend = TRUE,
  legend = "top",
  align = "hv"
)

ggsave(
  filename = "cutoff_detection_pairs_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 12,
  width = 18,
  units = "cm"
)


# TODO: continue with corrections of Noelle:
#' - change axis labels of individual plots: replace p-value for significance level
#' show only numbers in each axis (remove p_gt_0.05...)
#' increase font sizes of axis labels, and bar labels?




### Sensitivity all ----

ggarrange(p_net_sensitivity_p_cor,
          p_net_sensitivity_p_gt, 
          p_net_sensitivity_rii, 
          labels = c("(a)", "(b)", "(c)"),
          font.label = list(size = 10),
          hjust = -0.1,
          vjust = 1,
          ncol = 1, 
          nrow = 3, 
          common.legend = TRUE,
          legend = "top")

# save ggplot

ggsave(
  filename = "sensitivity_network_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 18,
  width = 8.5,
  units = "cm"
)



### Specificity all ----

ggarrange(p_net_specificity_p_cor,
          p_net_specificity_p_gt, 
          p_net_specificity_rii, 
          labels = c("(d)", "(e)", "(f)"),
          font.label = list(size = 10),
          hjust = -0.1,
          vjust = 1,
          ncol = 1, 
          nrow = 3,
          common.legend = TRUE,
          legend = "top")

# save ggplot

ggsave(
  filename = "specificity_network_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 18,
  width = 8.5,
  units = "cm"
)




## Cor: Plant vs Flower - Sensitivity ----

ggarrange(p_net_sensitivity_p_cor, 
          p_net_sensitivity_flower_p_cor, 
          labels = c("(a)", "(b)"),
          font.label = list(size = 10),
          hjust = -0.1,
          vjust = 1,
          ncol = 1, 
          nrow = 2,
          common.legend = TRUE,
          legend = "top")

# save ggplot

ggsave(
  filename = "sensitivity_plant_vs_flower.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 18,
  width = 8.5,
  units = "cm"
)

## Cor: Plant vs Flower - Specificity ----

ggarrange(p_net_specificity_p_cor, 
          p_net_specificity_flower_p_cor, 
          labels = c("(c)", "(d)"),
          font.label = list(size = 10),
          hjust = -0.1,
          vjust = 1,
          ncol = 1, nrow = 2,
          common.legend = TRUE,
          legend = "top")

# save ggplot

ggsave(
  filename = "specificity_plant_vs_flower.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 18,
  width = 8.5,
  units = "cm"
)



## Summary table methods ---------------------------------------------------

tab_methods <-
  bind_rows(
    tab_pairs_rii,
    tab_pairs_p_gt,
    tab_pairs_p_cor %>% mutate(method = "Plant cor"),
    tab_flower_p_cor %>% mutate(method = "Flower cor")) %>% 
  mutate(
    method = case_when(
      str_detect(cutoff, "RII") ~ "RII",
      str_detect(cutoff, "p_gt") ~ "Probabilistic",
      .default = method
    )
  ) %>% 
  mutate(cutoff = str_remove(cutoff, "RII_"),
         cutoff = str_remove(cutoff, "p_gt_"),
         cutoff = str_remove(cutoff, "p_")) %>% 
  select(-c(known_link, known_nonlink)) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>% 
  relocate(method, .before = region) %>% 
  write.table(., "output/plots/accuracy/tab_cooccur_methods_TP_TN.txt",
              sep = ",",
              quote = FALSE,
              row.names = FALSE)



# 4) Comparing/selecting co-occurrence method ------------------------------

# based on the last results we filter contingency tables of co-occurrences using cutoff of p_gt ≤ 0.2, p_cor ≤ 0.2 and RII ≥ 0.2 and then compare these results using a Chi-squared to see which method allows to predict more interactions of butterflies and plants

# real interactions

pairs %>% 
  filter(troph_interaction > 0) %>% 
  count(region)

# get values of TP per method (contingency table)

tab_summary <-
  tab_pairs_p_cor %>%
  select(region:TN) %>%
  filter(cutoff == "p_0.2") %>%
  bind_rows(tab_pairs_p_gt %>%
              select(region:TN) %>%
              filter(cutoff == "p_gt_0.2")) %>%
  bind_rows(tab_pairs_rii %>%
              select(region:TN) %>%
              filter(cutoff == "RII_0.2")) %>%
  arrange(region) 



## 4.1) Plot cut-off by method ------------------------------------------------

## ---- sensitivity 

tab_summary %>%
  # add a new FP col including only positive interactions
  mutate(FN = known_link - TP) %>%
  pivot_longer(cols = c(FN, TP),
               names_to = "detected",
               values_to = "n") %>%
  mutate(
    cutoff = case_match(
      cutoff,
      "p_0.2" ~ "Correlation",
      "p_gt_0.2" ~ "Probabilistic",
      "RII_0.2" ~ "RII null models")) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = fct_reorder(cutoff, n), 
      y = n,
      fill = detected)) +
  # facet per region
  facet_wrap(~region) + 
  geom_col() + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FN" = "Known interactions",
                               "TP" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1.5) +
  labs(
    x = "Co-occurrence method",
    y = "Count of links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )

# save plot 

ggsave(
  filename = "sensitivity_0.2_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)
  
## --- specificity

tab_summary %>%
  # add a new FP col including only positive interactions
  mutate(FP = known_nonlink - TN) %>%
  pivot_longer(cols = c(FP, TN),
               names_to = "detected",
               values_to = "n") %>%
  mutate(
    cutoff = case_match(
      cutoff,
      "p_0.2" ~ "Correlation",
      "p_gt_0.2" ~ "Probabilistic",
      "RII_0.2" ~ "RII null models")) %>% 
  # group data and calculate percentage of FN and TP per cutoff and region
  group_by(region, cutoff) %>% 
  mutate(perc = n/sum(n)*100) %>% 
  ungroup() %>%
  # Stacked bar plot
  ggplot(
    aes(
      x = fct_reorder(cutoff, n), 
      y = n,
      fill = detected)) +
  # facet per region
  facet_wrap(~region) + 
  geom_col() + 
  # increase the number of breaks 
  scale_y_continuous(n.breaks = 5) +
  # customize label keys 
  scale_fill_manual(values=c("cadetblue", "sandybrown"),
                    labels = c("FP" = "Known non-interactions",
                               "TN" = "Detected")) +
  geom_text(aes(x = cutoff, y = n,
                label = paste0(n, " (", round(perc, 1), "%)")),
            position = position_stack(vjust = 0.5), 
            size = 1.5) +
  labs(
    x = "Co-occurrence method",
    y = "Count of non-links",
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top"
  )

ggsave(
  filename = "specificity_0.2_null.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy",
  dpi = "retina",
  height = 10,
  width = 8.5,
  units = "cm"
)


## 4.2) Chi-squared tests ----------------------------------------------------

### Sensitivity ----

# separate data per region

tab_alb <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(cutoff:known_link, TP) %>% 
  mutate(FN = known_link - TP) %>% 
  select(-known_link) %>% 
  column_to_rownames("cutoff") 

tab_sch <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(cutoff:known_link, TP)  %>% 
  mutate(FN = known_link - TP) %>% 
  select(-known_link) %>% 
  column_to_rownames("cutoff") 

# run a Chi-squared test of independence 

# assumptions:
# 1) cases randomly sampled from the population
# 2) non-occurrences are included

# chi-squared test

(chisq_alb <- chisq.test(tab_alb))
(chisq_sch <- chisq.test(tab_sch))

# proportion tables for interpretation of results

prop.table(as.matrix(tab_alb), 2)
prop.table(as.matrix(tab_sch), 2)

# obtain observed, expected and residuals 

chisq_alb$observed
chisq_alb$expected
chisq_alb$residuals # Pearson residuals

chisq_sch$observed
chisq_sch$expected
chisq_sch$residuals # Pearson residuals


### Sensitivity no Spearman: RII vs. Probabilistic ---------------------------

# repeat analysis without Spearman to see if RII and Cooccur differ significantly

tab_alb_noSpearman <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(cutoff:known_link, TP) %>% 
  mutate(FN = known_link - TP) %>% 
  select(-known_link) %>%
  filter(cutoff != "p_0.2") %>% 
  column_to_rownames("cutoff") 

tab_sch_noSpearman <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(cutoff:known_link, TP)  %>% 
  mutate(FN = known_link - TP) %>% 
  select(-known_link) %>% 
  filter(cutoff != "p_0.2") %>% 
  column_to_rownames("cutoff") 

# run chi-test 

(chisq_alb_noSpearman <- chisq.test(tab_alb_noSpearman))
(chisq_sch_noSpearman <- chisq.test(tab_sch_noSpearman))

prop.table(as.matrix(tab_alb_noSpearman), 2)
prop.table(as.matrix(tab_sch_noSpearman), 2)

# obtain observed, expected and residuals 

chisq_alb_noSpearman$observed
chisq_alb_noSpearman$expected
chisq_alb_noSpearman$residuals # Pearson residuals

chisq_sch_noSpearman$observed
chisq_sch_noSpearman$expected
chisq_sch_noSpearman$residuals # Pearson residuals

# In both regions, significant difference for sensitivity between RII with null models with plant abundance and analytical p-values obtained from probabilistic method at a cutoff of p-value = 0.2. 



### Specificity -------------------------------------------------------------

tab_alb <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(cutoff, known_nonlink, TN) %>% 
  mutate(FP = known_nonlink - TN) %>% 
  select(-known_nonlink) %>% 
  column_to_rownames("cutoff")

tab_sch <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(cutoff, known_nonlink, TN) %>%
  mutate(FP = known_nonlink - TN) %>%
  select(-known_nonlink) %>%
  column_to_rownames("cutoff")

# chi-squared test

(chisq_alb <- chisq.test(tab_alb))
(chisq_sch <- chisq.test(tab_sch))

# proportion tables for interpretation of results

prop.table(as.matrix(tab_alb), 2)
prop.table(as.matrix(tab_sch), 2)



### Specificity RII vs Cooccur ------------------------------------------------

tab_alb <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(cutoff, known_nonlink, TN) %>% 
  mutate(FP = known_nonlink - TN) %>% 
  select(-known_nonlink) %>% 
  filter(cutoff != "p_0.2") %>% 
  column_to_rownames("cutoff")

tab_sch <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(cutoff, known_nonlink, TN) %>%
  mutate(FP = known_nonlink - TN) %>%
  select(-known_nonlink) %>%
  filter(cutoff != "p_0.2") %>% 
  column_to_rownames("cutoff")

# chi-squared test

(chisq_alb <- chisq.test(tab_alb))
(chisq_sch <- chisq.test(tab_sch))

# proportion tables for interpretation of results

prop.table(as.matrix(tab_alb), 2)
prop.table(as.matrix(tab_sch), 2)

# obtain observed, expected and residuals 

chisq_alb$observed
chisq_alb$expected
chisq_alb$residuals # Pearson residuals

chisq_sch$observed
chisq_sch$expected
chisq_sch$residuals # Pearson residuals

# Concl: the specificity obtained with the RII index was significantly higher than with the probabilistic method only for the region ALB but not for the region SCH. 


### Plant vs flower sensitivity --------------------------------------------

tab_summary <-
  tab_pairs_p_cor %>%
  select(region:TN) %>%
  filter(cutoff == "p_0.2") %>%
  bind_rows(tab_flower_p_cor %>%
              select(region:TN) %>%
              filter(cutoff == "p_0.2") %>% 
              mutate(cutoff = "p_flower_0.2")) %>% 
  mutate(FN = known_link - TP) %>% 
  select(-c(known_link, known_nonlink, TN))

# ALB

tab_alb <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(-region) %>%
  column_to_rownames("cutoff")

# SCH

tab_sch <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(-region) %>%
  column_to_rownames("cutoff") 

# run chi-squared test

(chisq_alb <- chisq.test(tab_alb))
(chisq_sch <- chisq.test(tab_sch))

# proportion tables for interpretation of results

prop.table(as.matrix(tab_alb), 1)
prop.table(as.matrix(tab_sch), 1)

# obtain observed, expected and residuals 

chisq_alb$observed
chisq_alb$expected
chisq_alb$residuals # Pearson residuals

chisq_sch$observed
chisq_sch$expected
chisq_sch$residuals # Pearson residuals


# No significant difference between Spearman correlation with the plant abundance or the flower abundance at a cutoff of p-value = 0.2 in both regions. 



### Plant vs flower Specificity ------------------------------------------------
  
tab_summary <-
  tab_pairs_p_cor %>%
  select(region, cutoff, TN, FP) %>%
  filter(cutoff == "p_0.2") %>%
  bind_rows(tab_flower_p_cor %>%
              select(region, cutoff, TN, FP) %>%
              filter(cutoff == "p_0.2") %>% 
              mutate(cutoff = "p_flower_0.2"))

# ALB

tab_alb <-
  tab_summary %>%
  filter(region == "ALB") %>%
  select(-region) %>%
  column_to_rownames("cutoff")

# SCH

tab_sch <-
  tab_summary %>%
  filter(region == "SCH") %>%
  select(-region) %>%
  column_to_rownames("cutoff") 

# chi-squared test

(chisq_alb <- chisq.test(tab_alb))
(chisq_sch <- chisq.test(tab_sch))

# proportion tables for interpretation of results

prop.table(as.matrix(tab_alb), 2)
prop.table(as.matrix(tab_sch), 2)

# obtain observed, expected and residuals 

chisq_alb$observed
chisq_alb$expected
chisq_alb$residuals # Pearson residuals

chisq_sch$observed
chisq_sch$expected
chisq_sch$residuals # Pearson residuals

# Concl: the specificity obtained with the RII index was significantly higher than with the probabilistic method only for the region ALB but not for the region SCH. 





# 5) Summarising final networks ---------------------------------------------

# add method col to each summary table 

tab_pairs_p_cor <- tab_pairs_p_cor %>% mutate(method = "Plant cor")
tab_pairs_p_gt <- tab_pairs_p_gt %>% mutate(method = "Probabilistic")
tab_pairs_rii <- tab_pairs_rii %>% mutate(method = "RII null models")
tab_flower_p_cor <- tab_flower_p_cor %>% mutate(method = "Flower cor")


pairs %>% 
  group_by(region) %>% 
  summarise(
    method = "Empirical",
    potential_links = 
      length(unique(plant_sp)) * length(unique(lepidoptera_sp)),
    realized_links = sum(troph_interaction > 0),
    connectance = realized_links/potential_links) %>%
  bind_rows(
    tab_pairs_rii %>% 
      filter(str_detect(cutoff, "0.2")) %>% 
      mutate(
        region,
        method, 
        potential_links = known_link + known_nonlink,
        realized_links = TP,
        connectance = realized_links/potential_links,
        .keep = "none")) %>% 
  bind_rows(
    map(
      list(tab_pairs_p_gt,
           tab_pairs_p_cor,
           tab_flower_p_cor),
      ~ .x %>% filter(str_detect(cutoff, "0.2")) %>% 
        mutate(
          region,
          method, 
          potential_links = known_link + known_nonlink,
          realized_links = TP,
          connectance = realized_links/potential_links,
          .keep = "none")) %>% 
      list_rbind()) %>% 
  mutate(across(where(is.numeric), ~round(.x, digits = 3))) %>% 
  relocate(method, .before = region) %>% 
  arrange(region) %>% 
  write.table(., "output/plots/accuracy/tab_summary_network.txt",
              sep = ",",
              quote = FALSE,
              row.names = FALSE)


# Conclusion ----

# The RII pairwise null models method performed better in terms of sensitivity but worst for specificity. The probabilistic method performed better for specificity. Spearman method with multiple testing correction using Benjamini & Hochberg's FDR.The RII method, performed better in specificity compared to the probabilistic method, but only for the region ALB and not in region SCH. In general, and since we are more interested to know the "real interactions" than the "real non interactions", we selected the probabilistic method for our futher analysis. 