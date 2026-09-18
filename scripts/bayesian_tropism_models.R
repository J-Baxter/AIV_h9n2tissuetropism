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
  mutate(across(c(starts_with('day'),
                  isolation_year),
                .fns = ~ scale(.x)[,1],
                .names = '{.col}_std')) %>%
  #select(strain, ends_with('std')) %>%
  filter(strain %in% rownames(A))

# Reorder A to match data order
A <- A[data_std$strain, data_std$strain]


################################### MAIN #######################################

# Q1: Does the trophism of H9N2 change through time?



# Q2: What is the consequence of the change in trophism on pathogenicity?
# 

