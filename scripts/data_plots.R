################################################################################
## Script Name:        <INSERT_SCRIPT_NAME_HERE>
## Purpose:            <BRIEFLY_DESCRIBE_SCRIPT_PURPOSE>
## Author:             James Baxter
## Date Created:       2026-09-21
################################################################################

############################### SYSTEM OPTIONS #################################
options(
  scipen = 6,     # Avoid scientific notation
  digits = 7      # Set precision for numerical display
)
memory.limit(30000000)

############################### DEPENDENCIES ###################################
# Load required libraries
library(tidyverse)
library(magrittr)


################################### DATA #######################################
# Read and inspect data
data <- read_csv('./data/00_master_dataset_wide.csv')

# Create the covariance matrix from the phylogenetic tree
A <- vcv.phylo(tree) 

# Standardize predictors for better prior specification
data_std <- data %>%
  mutate(across(c(starts_with('day'),
                  isolation_year),
                .fns = ~ scale(.x)[,1],
                .names = '{.col}_std')) %>%
  #select(strain, ends_with('std')) %>%
  filter(strain %in% rownames(A))

################################### MAIN #######################################
# Main analysis or transformation steps
data %>%
  select(starts_with('day5'),
         isolation_year) %>%
  mutate(day5_specificity_tissue = (day5_lung_mean - day5_nasal_mean) / (day5_lung_mean + day5_nasal_mean)) %>%
  pivot_longer(cols = starts_with('day5'),
               values_to = 'titre',
               names_to = 'tissue') %>%
  mutate(tissue = str_split_i(tissue, '_', 2)) %>%
  ggplot(aes(x = isolation_year, 
             y = titre)) +
  geom_point() +
  facet_wrap(~tissue, scales = 'free_y') + 
  theme_bw()


  
  
################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################