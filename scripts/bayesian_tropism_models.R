################################################################################
## Script Name:        Quantify changes in tropism and pathogenicity over time
## Purpose:           
## Author:             James Baxter
## Date Created:      2026-09-18
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)

############################## DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)
library(brms)
library(ape)
library(tidyverse)
library(tidybayes)
library(bayesplot)

################################### DATA #######################################
# Read and inspect data
data <- read_csv('./data/00_master_dataset_wide.csv')
tree <- read.tree('./trees/HA.fas.contree')

# Create the covariance matrix from the phylogenetic tree
A <- vcv.phylo(tree) 

# Standardize predictors for better prior specification
data_std <- data %>%
  mutate(day5_specificity_tissue = (day5_lung_mean - day5_nasal_mean) / (day5_lung_mean + day5_nasal_mean)) %>%
  mutate(across(c(starts_with('day'),
                  max_body_weight_loss_percent,
                  isolation_year),
                .fns = ~ scale(.x)[,1],
                .names = '{.col}_std')) %>%
  #select(strain, ends_with('std')) %>%
  filter(strain %in% rownames(A))

# Reorder A to match data order
A <- A[data_std$strain, data_std$strain]


################################### MAIN #######################################
# Q1: Does the trophism of H9N2 change through time?

# Q1a: lung mean day 5 ~ time

# Q1b: nasal mean day 5 ~ time

# Q1c: lung mean day 5 ~ time (no phylogenetic covariance)

# Q1d: nasal mean day 5 ~ time (no phylogenetic covariance)





# Q2: What is the consequence of the change in trophism on pathogenicity?
model_tropism_pathogenicity <- brm(
  max_body_weight_loss_percent_std ~ s(day5_specificity_tissue, k = 10) + 
    (1|gr(strain, cov = A)),
  data = data_std,
  data2 = list(A = A),
  family = gaussian(),
 # prior = prior_spec,
  chains = 4,
  iter = 5000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "cmdstanr" 
)

posterior_prediction <- predictions(model_tropism_pathogenicity,
                                        #newdata = 'balanced',
                                    by ='day5_specificity_tissue',
                                    conf_level = 0.9,
                                    re_formula = NA,
                                    type = "response", # posterior draws of the expected value 
                                    #newdata = NULL # Unit-level predictions for each observed value in the dataset (empirical distribution)
) %>% 
  as_tibble() %>%
  mutate(across(c(estimate, conf.low, conf.high), .fns = ~.x * sd(data$max_body_weight_loss_percent) + mean(data$max_body_weight_loss_percent))) %>% 
  ggplot(aes(x = day5_specificity_tissue)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.7, fill = 'grey') + 
  geom_line(aes(y = estimate)) +
  geom_point(data = data_std, aes(x = day5_specificity_tissue, y = max_body_weight_loss_percent)) +
  theme_bw() + 
  ylab('BW % loss') + 
  xlab('Specificity')


avg_predictions(model_tropism_pathogenicity,
                newdata = 'balanced',
                by ='isolation_year_std',
                conf_level = 0.9,
                type = "response", # posterior draws of the expected value 
                #newdata = NULL # Unit-level predictions for each observed value in the dataset (empirical distribution)
) %>% 
  as_tibble() %>%
  mutate(isolation_year = isolation_year_std * sd(data$isolation_year) + mean(data$isolation_year)) %>%
  mutate(across(c(estimate, conf.low, conf.high), .fns = ~.x * sd(data$max_body_weight_loss_percent) + mean(data$max_body_weight_loss_percent))) %>% 
  ggplot(aes(x = isolation_year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.7, fill = 'grey') + 
  geom_line(aes(y = estimate)) +
  geom_point(data = data_std, aes(x = isolation_year, y = max_body_weight_loss_percent)) +
  theme_bw() + 
  ylab('BW % loss') + 
  xlab('year')
