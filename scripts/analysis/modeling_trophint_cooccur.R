#' ----
#' title: Lepidoptera data analysis - modelling trophic interactions as a 
#' function of co-occurrences using P/A data 
#' author: Esteban Menares
#' date: 01.11.2022
#' ----


#' **AIM**: to test the predictive capacity of co-occurrences to predict 
#' trophic interactions of plant and Lepidoptera in grasslands in Germany. 


# Set up --------------------------------------------------------------

library(foreign)
library(MASS) # to fit the proportional odds logistic regression
library(Hmisc)
library(ggpubr)
library(reshape2) 
library(tidyverse) # data wrangling and visualization 
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
library(performance) # for model validation
library(ggeffects)   # for model predictions visualization 
library(multcompView) # for post-hoc plot 



# read in data

pairs <- 
  read_csv('data/raw/pairs.csv')

source('scripts/source_script.R')

# set theme for plots

theme_set(theme_sci(8)) 



# Data wrangling and exploration --------------------------------------

# using abundances at p-value ≤ 0.2

pairs_new <-
  pairs %>% 
  select(region,
         troph_pair_id, 
         troph_interaction, 
         co_occurrence,
         p_upper) %>% 
  mutate(
    p_upper = if_else(co_occurrence < 0, 1, p_upper),
    troph_interaction = factor(troph_interaction)) %>% 
  filter(p_upper <= 0.2) %>% 
  rename(troph = troph_interaction, 
         cooc = co_occurrence)


# create subsets per region for independent analysis using inverted ANOVA

alb <-
  pairs_new %>% 
  filter(region == 'ALB') 

sch <- 
  pairs_new %>% 
  filter(region == 'SCH')

## one at a time, table region and troph_int

lapply(pairs_new[, c("region", "troph")], table)

## three way cross tabs (xtabs) and flatten the table

ftable(xtabs(~ region + troph, data = pairs_new))

# drop unused levels of sch 

sch <- 
  sch %>% 
  mutate(troph = factor(troph)) %>% 
  droplevels(.$troph)

lapply(sch[, "troph"], table)

# examine the continuous variable co_occurrence

summary(alb$cooc)
summary(sch$cooc)
sd(alb$cooc)
sd(sch$cooc)

# examine the distribution of co_occurrence at every level of trophic interaction broken down by region

ggplot(pairs_new, 
       aes(x = troph, 
           y = cooc)) +
  geom_boxplot(size = 0.1, outlier.alpha = 0.1) +
  geom_point(alpha = 0.5, shape = 21) +
  facet_grid(~region, margins = FALSE) +
  scale_y_continuous(limits = c(0, 1)) + 
  labs(
    x = "Trophic interaction strength",
    y = "Co-occurrence association strength"
  ) + 
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed") + 
  geom_vline(xintercept = 3.5, color = "black", linetype = "dashed") +
  coord_flip()

# save plot 

ggsave(
  filename = "cooccur_by_troph_int_per_region.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 5,
  width = 7.5
) 




# Ordered Logistic Regression -----------------------------------------

# Below we use the polr command from the MASS package to estimate an 
# ordered logistic regression model. The command name comes from 
# proportional odds logistic regression, highlighting the proportional 
# odds assumption in our model. We specify Hess=TRUE to have the model 
# return the observed information matrix from optimization (called the 
# Hessian) which is used to get standard errors.
# https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/

# In R’s polr the ordinal logistic regression model is parameterized as
# logit (P(Y \le j)) = \beta_{j0} – \eta_{1}x_1 – \cdots – \eta_{p} x_p.


## fit ordered logit model and store results
m <- polr(troph ~ cooc + region, data = pairs_new, Hess=TRUE)
m_alb <- polr(troph ~ cooc, data = alb, Hess=TRUE)
m_sch <- polr(troph ~ cooc, data = sch, Hess=TRUE)

## view a summary of the model
summary(m)

# logit(P(Y ≤ 0.2)) = 1.1722 - (-1.409)*cooc - 0.8244 *SCH
# logit(P(Y ≤ 0.4)) = 2.0348 - (-1.409)*cooc - 0.8244 *SCH
# logit(P(Y ≤ 0.6)) = 2.1537 - (-1.409)*cooc - 0.8244 *SCH
# logit(P(Y ≤ 0.8)) = 3.0606 - (-1.409)*cooc - 0.8244 *SCH
# logit(P(Y ≤ 1))   = 4.6890 - (-1.409)*cooc - 0.8244 *SCH

summary(m_alb)

# logit(P(Y ≤ 0.2)) = 1.1722 - (-1.409)*cooc 
# logit(P(Y ≤ 0.4)) = 2.0348 - (-1.409)*cooc 
# logit(P(Y ≤ 0.6)) = 2.1537 - (-1.409)*cooc 
# logit(P(Y ≤ 0.8)) = 3.0606 - (-1.409)*cooc 
# logit(P(Y ≤ 1))   = 4.6890 - (-1.409)*cooc 
summary(m_sch)

# logit(P(Y ≤ 0.2)) = 0.7997 - (-0.4274)*cooc 
# logit(P(Y ≤ 0.6)) = 1.6274 - (-0.4274)*cooc 
# logit(P(Y ≤ 1))   = 3.4264 - (-0.4274)*cooc 

# Calculate p-values 

# One way to calculate a p-value in this case is by comparing the t-value 
# against the standard normal distribution, like a z test. Of course this 
# is only true with infinite degrees of freedom, but is reasonably 
# approximated by large samples, becoming increasingly biased as sample 
# size decreases. 

## Calculate p-values: store table
(ctable <- coef(summary(m)))
(ctable_alb <- coef(summary(m_alb)))
(ctable_sch <- coef(summary(m_sch)))

p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p_alb <- pnorm(abs(ctable_alb[, "t value"]), lower.tail = FALSE) * 2
p_sch <- pnorm(abs(ctable_sch[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p-value" = p))
(ctable_alb <- cbind(ctable_alb, "p-value" = p_alb))
(ctable_sch <- cbind(ctable_sch, "p-value" = p_sch))

## export results as table 
round(ctable_alb, digits = 3) %>% 
write.table(., "output/plots/accuracy/tab_p-values_ORL_ALB.txt",
              sep = ",",
              quote = FALSE,
              row.names = TRUE)

round(ctable_sch, digits = 3) %>% 
  write.table(., "output/plots/accuracy/tab_p-values_ORL_SCH.txt",
              sep = ",",
              quote = FALSE,
              row.names = TRUE)


# We can also get confidence intervals for the parameter estimates. These 
# can be obtained either by profiling the likelihood function or by using 
# the standard errors and assuming a normal distribution. 
# If the 95% CI does not cross 0, the parameter estimate is statistically 
# significant.

# default method gives profiled CIs
(ci <- confint(m)) 
(ci_alb <- confint(m_alb))
(ci_sch <- confint(m_sch))

# CIs assuming normality
confint.default(m) 
confint.default(m_alb) 
# 2.5 %     97.5 %
# cooc -2.280718 -0.5364507

confint.default(m_sch)
# 2.5 %   97.5 %
# cooc -1.882379 1.027607


# only the parameter estimate for the ALB region is significant, since it 
# doesn't cross zero. The estimates in the output are given in units of 
# ordered logits, or ordered log odds, so we can interpret this as:

# ALB: for a unit increase in co-occurrence, we would not expect a 1.22 
# increase in the expected value of trophic interaction in the log odds 
# scale... so what?
# https://stats.oarc.ucla.edu/r/faq/ologit-coefficients/ for interpretation..

# Instead of interpreting the odds of being in the th category or less, 
# we can interpret the odds of being greater than the th category by 
# exponentiating  itself. 

# The coefficients scaled in terms of logs are difficult to interpret. 
# Another way to interpret logistic regression models is to convert the 
# coefficients into odds ratios (OR). To get the OR and confidence intervals, 
# we just exponentiate the estimates and confidence intervals.

## OR and CI
exp(cbind(OR = coef(m), ci))
exp(cbind(OR = coef(m_alb), t(ci_alb)))
exp(cbind(OR = coef(m_sch), t(ci_sch)))

# These coefficients are called proportional odds ratios 

# OR = 1, no association 
# OR > 1, positively associated
# OR < 1, negatively associated

### Results interpretation ----

(1-0.2444891) * 100 # 75.55109
# For every unit -increase- in the co-occurrence value of a pair, the odds of being -more- likely to have a strong trophic interaction is 75.55 % lower

1/0.2444891 # 4.090162
# For every unit -decrease- in the co-occurrence value of a pair, the odds of being -less- likely to have a strong trophic interaction is multiplied by 4.1 times 


### Assumptions ---------------------------------------------------------

# check the section proportional odds assumption for more details
# https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/

# One of the assumptions underlying ordinal logistic (and ordinal probit) 
# regression is that the relationship between each pair of outcome groups 
# is the same. This is called the proportional odds assumption or the 
# parallel regression assumption. Because the relationship between all pairs
# of groups is the same, there is only one set of coefficients.

sf <- function(y) {
  c('Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)),
    'Y>=4' = qlogis(mean(y >= 4)),
    'Y>=5' = qlogis(mean(y >= 5)),
    'Y>=6' = qlogis(mean(y >= 6)))
}

(s <- with(alb, summary(as.numeric(troph) ~ cooc, fun=sf)))

glm(I(as.numeric(troph) >= 2) ~ cooc, family="binomial", data = alb)
glm(I(as.numeric(troph) >= 3) ~ cooc, family="binomial", data = alb)
glm(I(as.numeric(troph) >= 4) ~ cooc, family="binomial", data = alb)
glm(I(as.numeric(troph) >= 5) ~ cooc, family="binomial", data = alb)

# use the values in this table (df= s) to help us assess whether the 
# proportional odds assumption is reasonable for our model.
# for each range of cooc, subtract the predicted value for troph greater 
# than or equal to 2 (second factor level value i.e. 0.2) and troph greater 
# than or equal to 3 (third factor level value i.e. 0.4), and so on.
# Compare the values 

s[, 7] <- s[, 7] - s[, 6]
s[, 6] <- s[, 6] - s[, 5]
s[, 5] <- s[, 5] - s[, 4]
s[, 4] <- s[, 4] - s[, 3]
s[, 3] <- s[, 3] - s[, 3]
s # print

# plot the results of table. 
# The command which=1:n is a list of values indicating levels of y should 
# be included in the plot. The command pch=1:3 selects the markers to use, 
# and is optional, as are xlab='logit' which labels the x-axis, and main=' '
# which sets the main label for the graph to blank.

# If the proportional odds assumption holds, for each predictor variable,
# distance between the symbols for each set of categories of the dependent 
# variable, should remain similar. To help demonstrate this, we normalized 
# all the first set of coefficients to be zero so there is a common 
# reference point.

plot(s, which=1:6, pch=1:6, xlab='logit', main=' ', xlim=range(s[,3:7]))


### Model visualization -------------------------------------------------

newdat <- data.frame(
  cooc = rep(seq(from = 0, to = 1, length.out = 100), 4))

newdat <- cbind(newdat, predict(m_alb, newdat, type = "probs"))

##show first few rows
head(newdat)

lnewdat <-
  newdat %>% 
    as_tibble()  %>% 
    pivot_longer(
      cols = `0`:`1`,
      names_to = "Level",
      values_to = "Probability"
    )

## view first few rows
head(lnewdat)

p_orl_alb <-
  lnewdat %>% 
  ggplot(
    aes(x = cooc, 
        y = Probability, 
        colour = Level)) +
  geom_line() +
  scale_color_viridis_d() +
  labs(
    x = "Co-occurrence association strength",
    colour = "Interaction \nstrength"
  )
p_orl_alb

ggsave(
   filename = "OLR_cooccur_ALB.pdf",
   plot = last_plot(),
   path = "output/plots/accuracy/",
   dpi = "retina",
   height = 3,
   width = 3.5
 ) 


# Inverted ANOVA ------------------------------------------------------

# If you use only one continuous predictor, we could “flip” the model 
# around so that, say, cooc was the outcome variable and troph was the 
# predictor variable. Then you could run a one-way ANOVA. This isn’t a 
# bad thing to do if you only have one predictor variable (from the 
# logistic model), and it is continuous.
# https://stats.oarc.ucla.edu/r/dae/ordinal-logistic-regression/

## ALB unfiltered data ---------------------------------------------------- 

## ---- ALB - Modeling
  
lm1 <- aov(formula = cooc ~ troph, data = alb)
summary(lm1)

## export results as table 
capture.output(summary(lm1),file="output/plots/accuracy/tab_ANOVA_ALB.doc")

## ---- Post-hoc analysis 

# What is the effect of the treatment on the value ?
ANOVA = aov(lm1)

# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x = ANOVA, "troph", conf.level = 0.95)
TUKEY

pdf(file = "output/plots/accuracy/tukey_LM_coocur_troph_int_ALB.pdf",
     width = 5.5,
     height = 3.5,
     pointsize = 8)

# Tuckey test representation :
plot(TUKEY , las = 1 , col = "brown") 

dev.off()


## ---- ALB - Model validation
  
check_model(lm1)

# check autocorrelation for linear models
check_autocorrelation(aov(formula = cooc ~ troph, 
                         data = alb))
# Warning: Auto-correlated residuals detected (p < .001)!!!


# Check Predicted Distribution of Residuals and Predicted Distribution of Response
result <- check_distribution(lm1)
result

plot(result)

# **ALB - Model visualisation**
  
dat <- 
  ggpredict(lm1, terms = "troph") # gives the predicted values for Y based on terms Xn

p_lm_alb <- 
  dat %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                alpha = 0.5
                ) +
  labs(x = "Trophic interaction strength",
       y = "Co-occurrence association strength"
      ) + 
  coord_cartesian(ylim = c(0.3, 0.7)) 

p_lm_alb

# save plot 
 
ggsave(
  filename = "LM_coocur_troph_int_ALB.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 3,
  width = 3.5
) 



## SCH unfiltered data ----

# **ALB - Modeling** 
  
lm2 <- aov(formula = cooc ~ troph, 
          data = sch)
summary(lm2)

## export results as table 
capture.output(summary(lm2),file="output/plots/accuracy/tab_ANOVA_SCH.doc")

# **ALB - Model validation**
  
check_model(lm2)

# check autocorrelation for linear models
check_autocorrelation(aov(formula = cooc ~ troph, 
                         data = sch))
# OK: Residuals appear to be independent and not autocorrelated (p = 0.300).


# Check Predicted Distribution of Residuals and Predicted Distribution of Response
result <- check_distribution(lm2)
result

plot(result)

# **ALB - Model visualisation** 

dat <- ggpredict(lm2, terms = c("troph")) # gives the predicted values for Y based on terms Xn


p_lm_sch <- dat %>%
  ggplot(aes(x = x,
             y = predicted)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                alpha = .5) +
  labs(x = "Trophic interaction strength",
       y = "Co-occurrence association strength") + 
  coord_cartesian(ylim = c(0.3, 0.7)) 

p_lm_sch

# save plot 

ggsave(
  filename = "LM_coocur_troph_int_SCH.pdf",
  plot = last_plot(),
  path = "output/plots/accuracy/",
  dpi = "retina",
  height = 3,
  width = 3.5,
) 

