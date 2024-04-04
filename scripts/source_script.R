#' ---------------------------------
#' title: My compilation of functions
#' author: Esteban Menares
#' version: 2023.11.21
#' --------------------------------


#' **NOTE**: this is my private compilations of functions written by my self or by others, either taken from the internet, colleagues, or even adapted from other packages::functions




# 1. Descriptive Statistics -------------------------------------------------

summary_stats <- function(data, var) {
  data |> summarize(
    min = min({{ var }}, na.rm = TRUE),
    mean = mean({{ var }}, na.rm = TRUE),
    median = median({{ var }}, na.rm = TRUE),
    max = max({{ var }}, na.rm = TRUE),
    n = n(),
    n_miss = sum(is.na({{ var }})),
    sem = sd({{ var }})/sqrt(length({{ var }})),
    .groups = "drop"
  )
}



# 2. RMSE ----------------------------------------------------------------
# The root mean square error (RMSE) allows us to measure how far predicted values are from observed values in a regression analysis. In other words, how concentrated the data around the line of best fit. Lesser the RMSE, the better fit is our model.

# https://www.r-bloggers.com/2021/07/how-to-calculate-root-mean-square-error-rmse-in-r/

# RMSE = √[ Σ(Pi – Oi)2 / n ]

# where:

# Σ symbol indicates “sum”
# Pi is the predicted value for the ith observation in the dataset
# Oi is the observed value for the ith observation in the dataset
# n is the sample size

# e.g. RMSE <- sqrt(mean((data$actual - data$predicted)^2))

# a custom function

RMSE <- function(mod, obs) {
  sqrt(mean((mod - obs) ^ 2))
}

# mod = for model (fitted) values, 
# obs = for observed (true) values.






# 3. Keygen --------------------------------------------------------------

# Generate unique identifier for a single row

#' *NOTE*: to generate unique identifiers in bulk during data tidying 
#' process use the digest package with mutate and apply functions: e.g. 

# df %>% mutate(id = apply(., 1, digest::digest))

# or 

# df %>% mutate(id = apply(., 1, uuid::UUIDgenerate))

# or

generate_id <-
  function(n_char = 5, 
           my_date = Sys.Date()) {
    require(lubridate)
    c(letters,
      LETTERS,
      rep(0:9, n_char)) %>% 
      sample(size = n_char) %>% 
      str_c(collapse = '') %>% 
      str_c(lubridate::as_datetime(my_date), sep = "-")
  }

# 3. Correlation matrix to table -----------------------------------------


# Formatting a correlation matrix into a table with 4 columns

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

# Usage: 

# create and extract correlation coefficients and p-values
# cor_mat <- 
#   rcorr(
#     as.matrix(df), 
#     type = "spearman")
# 
# cor_df <- 
#   flattenCorrMatrix(
#     cor_mat$r, 
#     cor_mat$P)
# 
# # rename cols 
# colnames(cor_df)[
#   which(
#     names(cor_df) == "cor")] <- "plant_cor"
# 
# colnames(cor_df)[
#   which(
#     names(cor_df) == "p")] <- "plant_cor_p"



# 4. Remove cols with some or all NA ----------------------------------------

# cols with some NA rows

some_are_na <- 
  function(x) {
    !all(
      is.na(x))
  } 

# cols without any NA row

none_is_na <- 
  function(x) {
    !any(
      is.na(x))
  } 

# examples: 

# df %>% 
#   select(where(~none_is_na(.)))

# df %>% 
#   select(where(~some_are_na(.)))




# 5. Add leading zeros to EPID -------------------------------------------

#' Change Exploratories plot names without zeros to plot names with zeros
# https://github.com/biodiversity-exploratories-synthesis/Synthesis_useful_functions/blob/9cd963b4bd45f3f6f700df33ae9e6d4221c99275/BEplotZeros%20.R 

#' This function adds the sortable plot names to an existing dataset.
#' Based on a column with Bexis plot names (e.g. AEG1), the function creates a column
#' with sortable plot names (e.g. AEG01). 
#' Note: This function does the opposite of the "BEplotNonZeros.R" function
#'
#' @dat a dataset (data.frame) with at least one column with Bexis plot names (e.g. "AEG1)
#' @column the column name of the column containing the Bexis plot names
#' @plotnam the desired name of the new column with sortable plot names. Default is "PlotSTD"
#' @return the same dataset with an extra column containing the sortable plot names (e.g. "AEG01)
#' @author Caterina Penone
#' 
#' @examples
#' #create a dataset with a plot name column
#' dat <- data.frame(Plot_name = c("AEG1", "AEG2", "HEW4", "SEG8", "SEW10"), Values=1:5)
#' dat <- BEplotZeros(dat, column = "Plot_name", plotnam = "Sorted_plot_name")
#' 
#' @export
BEplotZeros <- function (dat, column, colname = "PlotSTD") {
  dat <- as.data.frame(dat)
  funz <- function(x)
    ifelse((nchar(as.character(x)) == 4),
           gsub("(.)$", "0\\1", x),
           as.character(x))
  dat[, colname] <- sapply(dat[, column], funz)
  
  return(dat)
}





# 6. Plot themes -------------------------------------------------------------

# My own theme for 7.5 in (190 mm) plots adapted from ggthemes::theme_few(). 
# This theme is based on the rules and examples from Stephen Few's Show Me the Numbers and "Practical Rules for Using Color in Charts".


theme_sci <- function (base_size = 5,
          base_family = "")
{
  gray <- "#4D4D4D"
  black <- "#000000"
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      line = element_line(colour = gray),
      rect = element_rect(fill = "white",
                          colour = NA),
      text = element_text(colour = black),
      axis.ticks = element_line(colour = gray),
      legend.key = element_rect(colour = NA),
      panel.border = element_rect(colour = gray),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted",
                                        colour = "gray50"),
      strip.background = element_rect(fill = "white",
                                      colour = NA)
    )
}



# 7.  Pairs profile  --------------------------------------------------

# adapted from cooccur::pair.profile()

pair.profile.EM <-
  function (mod) {
    ptab <- mod$results
    spp <- unique(c(ptab$sp1, ptab$sp2))
    profiles <- data.frame(
      spp = c(),
      rand = c(),
      pos = c(),
      neg = c()
    )
    for (i in as.character(spp)) {
      spptab <- ptab[ptab$sp1 == i | ptab$sp2 == i,]
      pos <- nrow(spptab[spptab$p_gt < 0.05,])
      neg <- nrow(spptab[spptab$p_lt < 0.05,])
      rand <- nrow(spptab) - (pos + neg)
      total <- rand + pos + neg
      if (total == 0) {
        profiles[i, "pos"] <- 0
        profiles[i, "neg"] <- 0
        profiles[i, "rand"] <- 100
        profiles[i, "spp"] <- i
      }
      else {
        profiles[i, "pos"] <- pos / total * 100
        profiles[i, "neg"] <- neg / total * 100
        profiles[i, "rand"] <- rand / total * 100
        profiles[i, "spp"] <- i
      }
    }
    profiles$cooccur <- profiles$pos + profiles$neg
    profiles <- profiles[with(profiles, order(cooccur)),]
    profiles$spp <-
      factor(x = profiles$spp,
             ordered = T,
             levels = profiles$spp)
    profiles$cooccur <- NULL
    pos <- mod$positive
    neg <- mod$negative
    rand <- mod$random
    total <- pos + neg + rand
    pos <- mod$positive / total * 100
    neg <- mod$negative / total * 100
    rand <- mod$random / total * 100
    if ("sp1_name" %in% colnames(ptab)) {
      spp_key <- unique(data.frame(
        sppnum = c(ptab$sp1, ptab$sp2),
        sppname = c(
          as.character(ptab$sp1_name),
          as.character(ptab$sp2_name)
        )
      ))
      names <-
        merge(
          x = data.frame(order = 1:length(row.names(profiles)),
                         profiles),
          y = spp_key,
          by.x = "spp",
          by.y = "sppnum",
          all.x = T
        )
      names <- names[with(names, order(order)),]
      names$order <- NULL
      names$spp <- NULL
      names$sppname <- factor(
        x = names$sppname,
        ordered = T,
        levels = names$sppname
      )
      sidebar <- data.frame(
        pos = pos,
        neg = neg,
        rand = rand,
        sppname = "All Species"
      )
      names <- rbind(names, sidebar)
      md <- reshape2::melt(names, id = (c("sppname")))
      sppname <- md$sppname
      value <- md$value
      variable <- md$variable
      p <-
        ggplot(data = md, aes(x = sppname, y = value, fill = variable)) +
        geom_bar(stat = "identity",
                 position = "stack",
                 width = 1)
      p <- p + scale_fill_manual(
        values = c("light blue", "#FFCC66",
                   "dark gray"),
        name = "",
        labels = c("positive", "negative",
                   "random")
      ) + theme(
        plot.title = element_text(
          vjust = 2,
          size = 20,
          face = "bold"
        ),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(
          size = 6,
          hjust = 0,
          vjust = 1,
          angle = -45
        ),
        panel.background = element_rect(fill = "white",
                                        colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + xlab("") +
        ylab("Percent of pairings")
      p <-
        p + ggtitle("Species Association Profile") + scale_y_continuous(expand = c(0.005,
                                                                                   0)) + scale_x_discrete(expand = c(0.005, 0))
      p <- p + geom_bar(
        stat = "identity",
        data = md[(md$sppname ==
                     "All Species"),],
        aes(sppname),
        alpha = 0,
        size = 0.1,
        color = "white",
        width = 1
      )
      p
    }
    else {
      sidebar <- data.frame(
        pos = pos,
        neg = neg,
        rand = rand,
        spp = "All Species"
      )
      profiles <- rbind(profiles, sidebar)
      md <- reshape2::melt(profiles, id = (c("spp")))
      p <-
        ggplot(data = md, aes(x = spp, y = value, fill = variable)) +
        geom_bar(stat = "identity",
                 position = "stack",
                 width = 1)
      p <- p + scale_fill_manual(
        values = c("light blue", "#FFCC66",
                   "dark gray"),
        name = "",
        labels = c("positive", "negative",
                   "random")
      ) + theme(
        plot.title = element_text(
          vjust = 2,
          size = 20,
          face = "bold"
        ),
        legend.text = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(
          size = 6,
          hjust = 0,
          vjust = 1
        ),
        panel.background = element_rect(fill = "white",
                                        colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + xlab("") +
        ylab("Percent of pairings")
      p <-
        p + ggtitle("Species Association Profile") + scale_y_continuous(expand = c(0.005,
                                                                                   0)) + scale_x_discrete(expand = c(0.005, 0))
      p <- p + geom_bar(
        stat = "identity",
        data = md[(md$spp ==
                     "All Species"),],
        aes(spp),
        alpha = 0,
        size = 0.1,
        color = "white",
        width = 1
      ) + theme(axis.text.x = element_text(
        hjust = 0,
        vjust = 1,
        angle = -45
      ))
      p
    }
  }



# 8. Diagnostic plots for binomial data ----------------------------------

# function for diagnostic plots for binomial distributed data shared by Joe Kolowski's GLMM course at SMSC 2023.

# model: your regression model (the value returned from a binomail glm)
# bin.size: the number of points to put in a bin
# xvar: the (numeric!) vector you would like to plot your residuals against. Defaults to fitted(model)
# the other settings are likely to be useful at their default values

# examples: (Parasites in cod data set) 
# plotDiagBinom(g, col=rgb(0,0,0,.3), bin.size=15, xlab="fitted values")
# plotDiagBinom(g, xvar= parasite$Length , col=rgb(0,0,0,.3), bin.size=15, xlab="Length")
# plotDiagBinom(g, xvar= parasite$Depth , col=rgb(0,0,0,.3), bin.size=15, xlab="Depth")
#

plotDiagBinom <-
  function(model,
           bin.size = 10,
           xvar = predict(model, type = "link"),
           fun = mean,
           resid.type = "pearson",
           xlab = "",
           h = 0,
           ylim = c(-1.1, 1.1),
           ...)
  {
    theDf <-
      as.data.frame(cbind(residuals(model, type = resid.type), xvar))
    names(theDf) <- c("residuals", "xvar")
    n <- nrow(theDf)
    idx <- order(theDf$xvar)
    theDf <- theDf[idx, ]
    theDf$group.id <-
      unlist(lapply(1:round(n / bin.size), function(x)
        rep(x, bin.size)))[1:n]
    
    avg.resid <- tapply(theDf$residuals, theDf$group.id, fun)
    avg.fit <- tapply(theDf$xvar, theDf$group.id, mean)
    
    plot(
      avg.resid ~ avg.fit ,
      xlab = xlab,
      ylab = "Pearson residuals",
      pch = 16,
      ylim = ylim,
      ...
    )
    abline(h = h)
  }

getBinnedBinomialData <- function(x, y, bin.size)
{
  theDf <- data.frame(x = x, y = y)
  n <- nrow(theDf)
  idx <- order(theDf$x)
  theDf <- theDf[idx, ]
  theDf$group.id <-
    unlist(lapply(1:round(n / bin.size), function(x)
      rep(x, bin.size)))[1:n]
  
  avg.y <- tapply(theDf$y, theDf$group.id, mean)
  lower.y <-
    tapply(
      theDf$y,
      theDf$group.id,
      FUN = function(x)
        quantile(x, probs = .25)
    )
  upper.y <-
    tapply(
      theDf$y,
      theDf$group.id,
      FUN = function(x)
        quantile(x, probs = .75)
    )
  avg.x <- tapply(theDf$x, theDf$group.id, mean)
  
  res <-
    data.frame(
      avg.y = avg.y,
      avg.x = avg.x,
      lower.y = lower.y,
      upper.y = upper.y
    )
  return(res)
}





# 9. Pairs panel ------------------------------------------------------

# function for diagnostic plots for binomial distributed data shared by Joe Kolowski's GLMM course at SMSC 2023.

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor.fac=.9, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r.sign <- sign(cor(x, y, use = "pairwise.complete.obs"))
  r <- abs(cor(x, y, use = "pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  cex.cor <- cex.cor.fac/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r, col=ifelse(r.sign < 0, "blue", "red"))
}

## put histograms on the diagonal

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}




# 10. Ordered coefplot -------------------------------------------------

# Created by Ben Bolker. Used for ordering estimates from a model from low to high coming as output from dotwhisker::dwplot(). 
# https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/bolker_mixedmodels/utils.R

# Example: 

### coefficient plots

## basic coefficient plot
# dotwhisker::dwplot(m1, effects="fixed") + geom_vline(xintercept = 0, lty = 2)

## ordered coefficient plot (from dwplot_ordered())
# dwplot_ordered(m1, effects = "fixed")

## compare plots
# dwplot_ordered(list(intercept_only = m1, random_slopes = m2), effects = "fixed")

dwplot_ordered <- function(model, ..., dodge = 0.5) {
  if (!require("broom.mixed")) stop("please install broom.mixed")
  multi_model <- is.list(model)
  tt <- if (!multi_model) tidy(model, ...) else purrr::map_dfr(model, tidy, .id = "model", ...)
  tt <- (tt
         %>% filter(term != "(Intercept)")
         %>% mutate(across(term, fct_reorder, estimate))
  )
  gg <- (ggplot(tt, aes(y=term, x = estimate,
                        xmin = estimate - 2*std.error,
                        xmax = estimate + 2*std.error)))
  if (multi_model) {
    gg <- gg + aes(colour = model)
    pd <- position_dodge(width = dodge)
  } else pd <- position_identity()
  gg1 <- (gg
          + geom_pointrange(position = pd)
          + geom_vline(xintercept = 0, lty = 2)
  )
  return(gg1)
}



# 11. Multidiversity --------------------------------------------------

# This function was created by Eric Allan 
# https://github.com/eric-allan/multidiversity/blob/master/.multidiv

#### a function to calculate multidiversity
#### threshold can be a proportion of the maximum e.g. 0.5 or can be the median or mean
#### it can also be a vector of the same length as the number of diversities to
#### allow different thresholds and therefore different weightings for different diversities, setting the threshold to NA drops the group from the calculation of multidiversity
#### threshold = FALSE calculates the mean of the scaled diversities or functions
#### scaling can be by any function specified in "sc", i.e. max, mean, sd etc., max is the default
#### centering by the mean is possible with cent = TRUE, to take z-scores of the diversities/processes, use sc="sd" and cent = TRUE 
#### "by" specifies a vector of the same length as nrow(x) to be used to split the data and use different
#### thresholds for the groups in "by"
#### "weights" allows different weightings for the different groups and should be a vector of the same length as the number of groups

multidiv <-
  function(x,
           threshold = FALSE,
           sc = "max",
           cent = FALSE,
           by = FALSE,
           weights = FALSE) {
    result <- matrix(nrow = nrow(x), ncol = 2)
    
    if (any(by != FALSE)) {
      xs <- split(x, by)
      xst <- list()
      
      for (i in 1:length(unique(by))) {
        xst[[i]] <-
          scale(xs[[i]],
                scale = apply(xs[[i]], 2, match.fun(sc), na.rm = T),
                center = cent)
      }
      x.stand <- do.call("rbind", xst)
    }
    
    else{
      x.stand <-
        scale(x,
              scale = apply(x, 2, match.fun(sc), na.rm = T),
              center = cent)
    }
    
    if (length(threshold) == 1) {
      gm <- apply(x, 1, function(x)
        (sum(complete.cases(x))))
    }
    
    else{
      x2 <-
        sweep(x, 2, threshold, "*")  ## remove diversities with NA threshold from calc. of how many measured
      gm <- apply(x2, 1, function(x)
        (sum(complete.cases(x))))
    }
    
    if (FALSE %in% threshold) {
      ### prevent error message if threshold is vector
      m <- apply(x.stand, 1, mean, na.rm = T)
    }
    
    else{
      if (any(c("median", "mean") %in% threshold)) {
        tf <- match.fun(threshold)
        
        x.thresh <-
          apply(x.stand, 2, function(x)
            (1 * (x > tf(x, na.rm = T))))
        gg <- apply(x.thresh, 1, sum, na.rm = T)
        m <- gg / gm
        
      }
      
      else{
        x.thresh <-
          1 * sweep(x.stand, 2, threshold, ">")  ### does each variable pass threshold?
        
        if (any(weights != FALSE)) {
          weights2 <-
            matrix(rep(weights, nrow(x.thresh)),
                   nrow = nrow(x.thresh),
                   byrow = TRUE)
          x.thresh <- x.thresh * weights2   ### multiply by weights
        }
        else{
        }
        gg <- apply(x.thresh, 1, sum, na.rm = T)
        m <- gg / gm
      }
    }
    
    result[, 1] <- m
    result[, 2] <- gm
    
    colnames(result) <- c("m", "groups measured")
    
    return(result)
  }



# 12. Normalising/Scaling and Inverse-Scaling ------------------------------

# based on https://medium.com/swlh/data-normalisation-with-r-6ef1d1947970

## 12.1 Min-Max ----

# Min-Max Normalization transforms x to x’ by converting each value of features to a range between 0 and 1, and this is also known as (0–1) Normalization. If the data has negative values the range would have been between -1 and 1.

norm_minmax <- function(x, na.rm = FALSE) {
  (x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
}

## 12.2 Robust Scalar ----

# Robust Scalar transforms x to x’ by subtracting each value of features by the median and dividing it by the interquartile range between the 1st quartile (25th quantile) and the 3rd quartile (75th quantile). 
# Can be used when we have outliers in the data

robust_scalar <- function(x) {
  (x - median(x)) / (quantile(x, probs = .75) - quantile(x, probs = .25))
}

## 12.3 Mean Normalization ----

# Mean Normalization transforms x to x’ the same way as Min-Max Normalization; one different thing is each value of the feature is first subtracted by the sample mean.

mean_norm_minmax <- function(x){
  (x- mean(x)) /(max(x)-min(x))
}

## 12.4 Inverse min-max scaling ----

inverse_minmax <- function(x, na.rm = FALSE) {
  x * (max(x, na.rm = na.rm) - min(x, na.rm = na.rm)) + min(x, na.rm = na.rm)
}





# 13.  Calculate species interactions per site -------------------------------

# Authors: Esteban Menares & Nöelle Schenk 

## ---- Define a function to calculate interactions per site

interactions_by_site <- function(occurrence_matrix,
                                 interactions_df,
                                 butterfly_species_col,
                                 plant_species_col,
                                 interaction_code_col,
                                 interaction_strength_col) {
  
    
    # Make sure mentioned columns are in interactions table and NOT in occurrence matrix
    if (!all(c(butterfly_species_col, 
               plant_species_col, 
               interaction_code_col, 
               interaction_strength_col) %in% names(interactions_df))) {
      stop("mentioned names are not column names in interactions_df. Did you confound occurrence_matrix and interactions_df?")
    }
    
    # Check if occurrence_matrix is a data frame or matrix
    if (is.data.frame(occurrence_matrix)) {
      site_names <- rownames(occurrence_matrix)
      occurrence_matrix <- as.matrix(occurrence_matrix)
    } else if (is.matrix(occurrence_matrix)) {
      site_names <- rownames(occurrence_matrix)
    } else {
      stop("occurrence_matrix must be a data frame or matrix.")
    }
    
    # Check for NA values in occurrence_matrix
    if (any(is.na(occurrence_matrix))) {
      stop("NA values found in occurrence_matrix.")
    }
    
    # Check for NA values in interactions_df
    if (any(is.na(interactions_df))) {
      stop("NA values found in interactions_df.")
    }
    
    # record initial number of sites and interactions
    n_sites <- nrow(occurrence_matrix)
    n_interactions <- nrow(interactions_df)
    
    # create result matrix: 
    # columns are interactions (names are interaction codes), rows are sites
    # Default value is 0, if both species occur, value is overwritten by interaction strength
    interactions_per_site <- matrix(0, nrow = n_sites, ncol = n_interactions)
    colnames(interactions_per_site) <- interactions_df[[interaction_code_col]]
    rownames(interactions_per_site) <- rownames(occurrence_matrix)
    
    for (i in 1:n_sites) {
      for (j in 1:n_interactions) {
        
        # getting current row of the interactions table
        butterfly_species <- interactions_df[[butterfly_species_col]][j]
        plant_species <- interactions_df[[plant_species_col]][j]
        interaction_code <- interactions_df[[interaction_code_col]][j] 
        interaction_strength <- interactions_df[[interaction_strength_col]][j]
        
        if (butterfly_species %in% colnames(occurrence_matrix) &&
            plant_species %in% colnames(occurrence_matrix)) {
          butterfly_occurrence <- occurrence_matrix[i, butterfly_species]
          plant_occurrence <- occurrence_matrix[i, plant_species]
          
          # Check if both occurrences are greater than zero (see Noelle comment below)
          if (butterfly_occurrence > 0 && plant_occurrence > 0) {
            interaction_strength <- as.numeric(interactions_df[j, interaction_strength_col])
            
            interactions_per_site[i, j] <- interaction_strength
          }
        }
      }
    }
    return(interactions_per_site)
  }


## ---- Toy data and test function

# Sample butterfly and plant occurrence matrix
# butterfly_matrix <-
#   matrix(
#     c(1, 2, 0, 3, 0, 0, 1, 0, 1, 2), # added a zero here, to have a case where Lep and plants both do not occur
#     nrow = 5,
#     ncol = 2,
#     dimnames = list(
#       c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5"),
#       c("Lepi1", "Lepi2")
#     )
#   )
# plant_matrix <-
#   matrix(
#     c(2, 1, 3, 0, 2, 0, 1, 2, 0, 1),
#     nrow = 5,
#     ncol = 2,
#     dimnames = list(
#       c("Site 1", "Site 2", "Site 3", "Site 4", "Site 5"),
#       c("Plant1", "Plant2")
#     )
#   )

# Combine butterfly and plant occurrence matrices, make sure rownames are matched correctly using match() for matrices

# occurrence_matrix <-
#   cbind(butterfly_matrix,
#         plant_matrix[match(rownames(butterfly_matrix), 
#                            rownames(plant_matrix)), ])

# interactions_df2 <- data.frame(
#   but_species = c("Lepi1", "Lepi2", "Lepi1", "Lepi2", "Lepi3"),
#   pla_species = c("Plant1", "Plant2", "Plant2", "Plant3", "Plant1"),
#   int_code = c(
#     "Lepi1-Plant1",
#     "Lepi2-Plant2",
#     "Lepi1-Plant2",
#     "Lepi2-Plant3",
#     "Lepi3-Plant1"
#   ),
#   int_strength = c(0.5, 0.3, 0.2, 0.4, 0.1)
# )


# EXAMPLE: Call the function with the provided data and column names
# interactions_per_site <-
#   interactions_by_site(
#     occurrence_matrix,
#     interactions_df2,
#     butterfly_species_col = "but_species",
#     plant_species_col = "pla_species",
#     interaction_code_col = "int_code",
#     interaction_strength_col = "int_strength"
#   )




# 14. Calculate network-level indices ---------------------------------

# use this function to calculate network-level (or other) indices from a list of networks ("webs_alb" and "webs_sch" in the examples) using bipartite, generate null models, compare obs vs. expected values, and calculate different final values including (and not restricted to this score): 
# obs_nestedness (of each network)
# null_nestedness (mean of null models)
# final_nestedness (obs - mean_null)
# p_value (proportion of non-random values from 1000 null models)
# LCI (lower CI)
# UCI (upper CI)
# z_scores ((obs - mean_null) / SD_null)

## ---- Turn data frames into incidence matrices in case you have a dataframe and not matrices. Check ?bipartite::frame2webs for details. 

# webs_alb <- bipartite::frame2webs(
#   int_ALB,
#   varnames = c("plant_sp",
#                "lepidoptera_sp",
#                "site",
#                "interaction_strength"),
#   type.out = "list",
#   emptylist = TRUE
# )
# 
# webs_sch <- bipartite::frame2webs(
#   int_SCH,
#   varnames = c("plant_sp",
#                "lepidoptera_sp",
#                "site",
#                "interaction_strength"),
#   type.out = "list",
#   emptylist = TRUE
# )

## ---- Create a function to calculate nestedness (or any other indices) and related statistics for a single network

calculate_nestedness <- function(network) {
  # Calculate observed nestedness using NODF index
  obs <- unlist(networklevel(network, index = "NODF"))
  # Generate null models based on Vazquez's method
  nulls_vazquez <- bipartite::nullmodel(network, N = 999, method = "vaznull")
  # Calculate NODF for null models
  null_vazquez_NODF <- unlist(purrr::map(nulls_vazquez, networklevel, index = "NODF"))
  # calculate mean nestedness of null models
  mean_null <- mean(null_vazquez_NODF)
  # calculate final nestedness = obs - mean(nulls). Higher value is better
  final_nestedness <- obs - mean_null
  # Calculate p-value
  praw <- sum(null_vazquez_NODF > obs) / length(null_vazquez_NODF)
  p_value <- if_else(praw > 0.5, 1 - praw, praw)
  # Compute 95% CIs
  LCI <- quantile(null_vazquez_NODF, c(0.025), na.rm = TRUE)
  UCI <- quantile(null_vazquez_NODF, c(0.975), na.rm = TRUE)
  # Compute standard deviation of z-scores
  sd_null <- sd(null_vazquez_NODF)
  # Compute z-scores for the network
  z_scores <- (obs - mean_null) / sd_null
  
  data.frame(
    obs_nestedness = obs,
    null_nestedness = mean_null,
    final_nestedness = final_nestedness,
    p_value = p_value,
    LCI = LCI,
    UCI = UCI,
    z_scores = z_scores
  )
}

## ---- Use purrr::map to iterate through each element in the webs list and combine the results row-wise

# The following step could be run also with one list. This step can be modified to output a list instead of a dataframe. 

# nestedness_list_ALB <- 
#   purrr::map(webs_alb, ~calculate_nestedness(.x), .progress = TRUE)
# nestedness_list_SCH <- 
#   purrr::map(webs_sch, ~calculate_nestedness(.x), .progress = TRUE)
# 
# # Combine the results into a single data frame using list_rbind
# nestedness_ALB <- purrr::list_rbind(nestedness_list_ALB)
# nestedness_SCH <- purrr::list_rbind(nestedness_list_SCH)

## ---- Reset row names
# row.names(nestedness_ALB) <- NULL
# row.names(nestedness_SCH) <- NULL

## ---- Add a "EPID" (sites) column to store the names of the networks
# Here it takes the names of each site (i.e. each network) to name the data frames.

# nestedness_ALB$EPID <- names(webs_alb)
# nestedness_SCH$EPID <- names(webs_sch)

## ---- Combine both dfs
# nestedness <- bind_rows(nestedness_ALB, nestedness_SCH)
# nestedness <- nestedness %>% relocate(EPID, .before = obs_nestedness)
# 
# network_scores <- nestedness %>% 
#   left_join(
#     scores_site_final %>% 
#       select(EPID, n_lepi, n_plant, sum_troph_int_pa, connectance_pa),
#     by = "EPID"
#   ) %>% 
#   relocate(EPID, n_lepi, n_plant, sum_troph_int_pa, connectance_pa, 
#            .after = EPID)




# 15. Calculate R2  ---------------------------------------------------

get_R2 <- function(model) {
  # Check if the input model has 'deviance' and 'null.deviance' components
  if (!all(c("deviance", "null.deviance") %in% names(model))) {
    stop("Input model must have 'deviance' and 'null.deviance' components.")
  }
  
  # Calculate 1 - deviance/null.deviance
  `R^2` <- 1 - model$deviance / model$null.deviance
  
  return(`R^2`)
}




# 16. Clean objects names ------------------------------------------------------
# function for cleaning objects names. Shared by Brian Evans's R-GIS course at SMSC 2023.


# Clean object names:

clean_names <- 
  function(df) {
    df %>% 
      set_names(
        names(.) %>% 
          str_replace_all("[:blank:]", "_") %>% 
          str_replace_all("([a-z])([A-Z]){1}", "\\1_\\2") %>% 
          tolower())
  }





# 17. LonLat to utm -----------------------------------------------------------
# function for cleaning objects names. Shared by Brian Evans's R-GIS course at SMSC 2023.


# Function for getting an epsg code from lonlat data (from Lovelace, chapter 7,
# https://geocompr.robinlovelace.net/reproj-geo-data.html):

lonlat_to_utm <-
  function(lonlat) {
    utm <-
      (floor((lonlat[1] + 180) / 6) %% 60) + 1
    if(lonlat[2] > 0) {
      utm + 32600
    } else{
      utm + 32700
    }
  }




# 18. Find your UTM zone --------------------------------------------------
# function for cleaning objects names. Shared by Brian Evans's R-GIS course at SMSC 2023.

# You can calculate your UTM zone number using the following steps:

# Add 180 to your your target longitude;
# Divide the resultant value by 6;
# Round up to the nearest whole number (here, I’m using the ceiling() function to round up)

my_utm_zone <-
  function(lon) {
    ceiling(
      (180 + lon)/6)
  }

# For example, the longitudinal center of the District of Columbia is approximately, -77° – we can calculate the UTM zone number of the city with:

# my_utm_zone(-77)

# Note: The full UTM zone would be described as 18N because it is north of the Equator. If the city were south of the Equator, this would be 18S. Knowing the UTM zone is only part of the battle though – we need to know the CRS associated with that zone. For this we use  the function, lonlat_to_utm(), to determine the EPSG code associated with our target UTM zone! 




# 19. Dimensions of dataframe without Zero or NAs -----------------------------

## ---- Dimensions without NAs

# ref_NA_row indicates the reference column to use to remove rows with NAs 

dim_no_NA <- function(x, ref_NA_row = ref_NA_row){
  x %>%
    as_tibble() %>% 
    select(all_of(ref_NA_row):last_col()) %>%
    select(where( ~ any(!is.na(.)))) %>% 
    filter(!is.na(ref_NA_row)) %>% 
    dim()
} 


## ---- Dimensions no NAs no Zero 

# ref_NA_row indicates the reference column to use to remove rows with NAs 

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
