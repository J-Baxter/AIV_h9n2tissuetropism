
Define priors for the model
prior_spec <- c(
  # Intercept - centered on observed mean
  prior(normal(0, 2), class = "Intercept", lb = 0.00000000000001),
  
  # Slope - effect of year on trait (standardized units)
  # Small effect expected unless strong selection for change
  prior(normal(0, 1), class = "b", coef = "isolation_year_std"),
  
  # Phylogenetic correlation strength (0-1, where 1 = full phylogenetic correlation)
  # You can specify different priors here depending on expectations
  prior(exponential(1), class = "sigma")#,  # Residual SD
  #prior(exponential(1), class = "lagr")   # Phylogenetic smoothing parameter
)

# Define MCMC settings
CHAINS <- 4
ITER <- 10000
WARMUP <- ITER/10
CORES <- 4
SEED <- 4472

# Prior predictive check
model_phylo_ppc <- brm(
  day5_nasal_mean ~ isolation_year_std + (1|gr(strain, cov = A)),
  data = data_std,
  data2 = list(A = A),
  family = lognormal(),
  prior = prior_spec,
  sample_prior = 'only',
  chains = CHAINS,
  iter = ITER,
  warmup = WARMUP,
  cores = CORES,
  seed = SEED,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr")

pp_check(model_phylo_ppc)



# Model 1: Phylogenetic regression (trait ~ year)
# The (1|gr(species, cov = A)) term models phylogenetic correlation
# br() = brms residual variance parameterization

model_phylo <- brm(
  day5_lung_mean_std ~ isolation_year_std + (1|gr(strain, cov = A)),
  data = data_std,
  data2 = list(A = A),
  family = gaussian(),
  prior = prior_spec,
  chains = 4,
  iter = 5000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  save_all_pars = TRUE,
  backend = "cmdstanr"  # Use cmdstanr for faster sampling (optional)
  # backend = "rstan"   # Alternative: use rstan (default, slower)
)

# Model 2: Null model (trait ~ 1, no year effect)
# For model comparison
model_null <- brm(
  day5_nasal_mean_std ~ 1 + (1|gr(strain, cov = A)),
  data = data_std,
  data2 = list(A = A),
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4,
  iter = 5000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  backend = "cmdstanr",
  save_all_pars = TRUE
)


# Prior Predictive

# Fitted Models


# Diagnositics - MCMC

# Diagnositcs - Model

# Predictions

# The entire dataset is replicated once for each unique combination of variables, and predictions are made
posterior_prediction <- predictions(model_phylo,
            by ='isolation_year_std',
            conf_level = 0.9,
            type = "response", # posterior draws of the expected value 
            newdata = NULL # Unit-level predictions for each observed value in the dataset (empirical distribution)
) %>% 
  mutate(isolation_year = isolation_year_std * sd(data$isolation_year) + mean(data$isolation_year),
         across(c(estimate, conf.low, conf.high), .fns = ~.x * sd(data$day5_lung_mean) + mean(data$day5_lung_mean))) %>% 
  ggplot(aes(x = isolation_year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.7, fill = 'grey') + 
  geom_line(aes(y = estimate)) +
  geom_point(data = data, aes(x = isolation_year, y = day5_lung_mean)) +
  theme_bw() + 
  ylab('Day 5 Lung Titre') + 
  xlab('Isolation Year')


conditional_effect <- predictions(model_phylo,
            by ='isolation_year_std',
            conf_level = 0.9,

            re_formula = NA,
            newdata = NULL # Unit-level predictions for each observed value in the dataset (empirical distribution)
) %>% 
  mutate(isolation_year = isolation_year_std * sd(data$isolation_year) + mean(data$isolation_year),
         across(c(estimate, conf.low, conf.high), .fns = ~.x * sd(data$day5_lung_mean) + mean(data$day5_lung_mean)))%>% 
  ggplot(aes(x = isolation_year)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.7, fill = 'grey') + 
  geom_line(aes(y = estimate)) +
  geom_point(data = data, aes(x = isolation_year, y = day5_lung_mean)) +
  theme_bw() + 
  ylab('Day 5 Lung Titre') + 
  xlab('Isolation Year')

cowplot::plot_grid(posterior_prediction, conditional_effect, nrow = 1, align = 'h', axis = 'tb')





################################### OUTPUT #####################################
# Save output files, plots, or results

#################################### END #######################################
################################################################################  

# Bayesian Phylogenetic Regression using brms
# =============================================
# Evaluate whether a continuous trait has changed through time while 
# accounting for phylogenetic non-independence using Bayesian methods

# Install required packages (uncomment to run)
# install.packages(c("brms", "ape", "tidyverse", "tidybayes", "bayesplot", "ggplot2"))
# install.packages("cmdstanr")  # For faster sampling



# =============================================================================
# 1. PREPARE DATA
# =============================================================================

# Create example data or load your own
# Required: species names as row names, trait and year columns



# Load phylogenetic tree
# Example: tree <- read.tree("path/to/tree.nwk")



# =============================================================================
# 2. SPECIFY PHYLOGENETIC COVARIANCE MATRIX
# =============================================================================

# Create the covariance matrix from the phylogenetic tree
A <- ape::vcv.phylo(tree)  # Phylogenetic covariance matrix
# Alternative: A <- ape::vcv(tree, model = "Brownian")

# Prepare for brms (species must be in exact order as data)
data$species <- rownames(data)
rownames(A) <- tree$tip.label



# =============================================================================
# 3. SPECIFY PRIORS
# =============================================================================



# View default priors (uncomment to see what brms suggests)
# get_prior(trait_std ~ year_std + (1|gr(species, cov = A)),
#           data = data, family = gaussian(),
#           data2 = list(A = A))

# =============================================================================
# 4. FIT BAYESIAN PHYLOGENETIC REGRESSION MODELS
# =============================================================================



# =============================================================================
# 5. EXAMINE DIAGNOSTICS
# =============================================================================

# Summary of posterior distributions
summary(model_phylo)

# Check for convergence (Rhat < 1.01 is good, < 1.05 is acceptable)
cat("\nRhat values (convergence diagnostic):\n")
print(rhat(model_phylo))

# Effective sample size ratio
cat("\nEffective sample size ratio:\n")
print(neff_ratio(model_phylo))

# Trace plots - visualize MCMC chains
plot(model_phylo)

# More detailed diagnostic plots
par(mfrow = c(2, 2))

# Trace plots for specific parameters
posterior_samples <- as_draws_df(model_phylo)

# Density plots
bayesplot::mcmc_areas(model_phylo, prob = 0.95) + 
  labs(title = "Posterior Distribution of Parameters")

# Trace plots
bayesplot::mcmc_trace(model_phylo,
                      pars = c("b_year_std", "sigma"),
                      n_warmup = 1000)

par(mfrow = c(1, 1))

# =============================================================================
# 6. POSTERIOR PREDICTIONS & UNCERTAINTY
# =============================================================================



# Extract posterior samples
posterior_samples <- as_draws_df(model_phylo)

# Summarize key parameter (slope for year effect)
year_effect <- posterior_samples %>%
  pull(b_isolation_year_std) %>%
  summarise(
    Mean = mean(.),
    Median = median(.),
    SD = sd(.),
    Q2.5 = quantile(., 0.025),
    Q97.5 = quantile(., 0.975)
  )

cat("\n===== POSTERIOR SUMMARY: Year Effect on Trait =====\n")
print(year_effect)

# Probability that the effect is positive
prob_positive <- mean(posterior_samples$b_isolation_year_std > 0)
cat("\nProbability that year effect is positive:", prob_positive, "\n")

# Get predictions for new data
new_data <- expand_grid(
  year_std = seq(min(data$isolation_year_std), max(data$isolation_year_std), length.out = 20),
  species = "Species_A"  # Population-level prediction
)

# Conditional predictions (account for phylogenetic structure)
pred_conditional <- posterior_epred(
  model_phylo,
  newdata = new_data,
  re_formula = NA  # Population-level prediction (ignoring species intercepts)
)

# Add predictions to new_data
new_data$pred_mean <- apply(pred_conditional, 2, mean)
new_data$pred_lower <- apply(pred_conditional, 2, quantile, prob = 0.025)
new_data$pred_upper <- apply(pred_conditional, 2, quantile, prob = 0.975)

# Back-transform from standardized scale (optional)
new_data <- new_data %>%
  mutate(
    year_original = year_std * sd(data$isolation_year) + mean(data$isolation_year),
    trait_pred = pred_mean * sd(data$day3_nasal_mean) + mean(data$day3_nasal_mean),
    trait_lower = pred_lower * sd(data$day3_nasal_mean) + mean(data$day3_nasal_mean),
    trait_upper = pred_upper * sd(data$day3_nasal_mean) + mean(data$day3_nasal_mean)
  )

# =============================================================================
# 7. VISUALIZATION
# =============================================================================

# Plot observed data with posterior predictions
plot_data <- data %>%
  mutate(
    trait_original = day3_nasal_mean,
    year_original = isolation_year
  )

p1 <- ggplot(plot_data, aes(x = year_original, y = trait_original)) +
  geom_point(size = 4, alpha = 0.6, color = "darkblue") +
  geom_ribbon(data = new_data, 
              aes(x = year_original, ymin = trait_lower, ymax = trait_upper),
              alpha = 0.3, fill = "blue", inherit.aes = FALSE) +
  geom_line(data = new_data,
            aes(x = year_original, y = trait_pred),
            color = "blue", size = 1, inherit.aes = FALSE) +
  labs(x = "Year", y = "Trait Value",
       title = "Trait Change Through Time (Bayesian Phylogenetic Regression)",
       subtitle = "Blue line = posterior mean, ribbon = 95% credible interval") +
  theme_minimal() +
  theme(text = element_text(size = 12))

print(p1)

# Posterior distribution of the year effect
p2 <- posterior_samples %>%
  ggplot(aes(x = b_year_std)) +
  geom_density(fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(x = "Effect of Year (standardized)",
       y = "Density",
       title = "Posterior Distribution of Year Effect") +
  theme_minimal()

print(p2)

# =============================================================================
# 8. MODEL COMPARISON (LOO-IC)
# =============================================================================

# Add LOO criterion to models (for model comparison)
model_phylo <- add_criterion(model_phylo, "loo")
model_null <- add_criterion(model_null, "loo")

# Compare models
loo_compare(model_phylo, model_null)

# Interpretation:
# - Negative elpd_diff favors the first model
# - If |elpd_diff| > 4, there's meaningful difference
# - If |elpd_diff| < 4, models are similar

# =============================================================================
# 9. CHECK PRIOR SENSITIVITY (Optional but Recommended)
# =============================================================================

# Fit model with weakly informative priors to compare
prior_weak <- c(
  prior(normal(0, 10), class = "Intercept"),
  prior(normal(0, 5), class = "b")
)

model_phylo_weak <- brm(
  trait_std ~ year_std + (1|gr(species, cov = A)),
  data = data,
  data2 = list(A = A),
  family = gaussian(),
  prior = prior_weak,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  seed = 123,
  refresh = 0  # Suppress iteration printing
)

# Compare posterior estimates
cat("\n===== PRIOR SENSITIVITY CHECK =====\n")
cat("Strong prior - Year effect:\n")
print(posterior_samples %>% pull(b_year_std) %>% 
        summarise(Mean = mean(.), SD = sd(.)))

cat("\nWeak prior - Year effect:\n")
posterior_weak <- as_draws_df(model_phylo_weak)
print(posterior_weak %>% pull(b_year_std) %>% 
        summarise(Mean = mean(.), SD = sd(.)))

# =============================================================================
# 10. EXTRACT AND VISUALIZE PHYLOGENETIC EFFECTS
# =============================================================================

# Get species-specific intercepts (random effects)
ranef_species <- ranef(model_phylo)$species[, , "Intercept"]

species_effects <- data.frame(
  species = rownames(ranef_species),
  intercept = ranef_species[, "Estimate"],
  lower = ranef_species[, "Q2.5"],
  upper = ranef_species[, "Q97.5"]
)

# Plot phylogenetic effects (trait values adjusted for phylogeny)
p3 <- ggplot(species_effects, aes(x = reorder(species, intercept), y = intercept)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(x = "Species", y = "Phylogenetic Intercept",
       title = "Species-specific Deviations (Accounting for Phylogeny)") +
  theme_minimal()

print(p3)

# =============================================================================
# 11. POSTERIOR PREDICTIVE CHECKS
# =============================================================================

# Check if model predictions match observed data
pp_check(model_phylo, ndraws = 100)

# More detailed posterior predictive checks
pp_check(model_phylo, type = "stat", stat = "mean")
pp_check(model_phylo, type = "error_binned")
pp_check(model_phylo, type = "scatter")

# =============================================================================
# 12. RESULTS SUMMARY
# =============================================================================

results_summary <- data.frame(
  Parameter = "Year Effect (standardized)",
  Posterior_Mean = mean(posterior_samples$b_year_std),
  Posterior_SD = sd(posterior_samples$b_year_std),
  CI_Lower = quantile(posterior_samples$b_year_std, 0.025),
  CI_Upper = quantile(posterior_samples$b_year_std, 0.975),
  Prob_Positive = prob_positive,
  Model_LOO = model_phylo$criteria$loo$estimates["looic", "Estimate"]
)

print(results_summary)

# Save results
# write.csv(results_summary, "brms_phylo_results.csv", row.names = FALSE)
data_std %>%
  mutate(tropism_specificity_day5 =( day5_lung_mean - day5_nasal_mean)/(day5_nasal_mean+day5_lung_mean)) %>%
  ggplot(aes(x = isolation_year, y = tropism_specificity_day5)) +
  geom_point()

data_std %>%
  ggplot(aes(x = isolation_year, y = max_body_weight_loss_percent)) +
  geom_point()
# Save model for later use
# saveRDS(model_phylo, "model_phylo.rds")
# model_phylo <- readRDS("model_phylo.rds")

# Do both traits show directional change over time?
tropism_specificity_day5 = (lung_mean_5days - nasal_mean_5days) / 
  (lung_mean_5days + nasal_mean_5days)

trend_tropism <- brm(
  tropism_specificity_day5 ~ year_std + (1|gr(strain_id, cov = A)), 
  data = analysis_data, data2 = list(A = A), family = gaussian()
)

trend_pathogenicity <- brm(
  weight_loss ~ year_std + (1|gr(strain_id, cov = A)), 
  data = analysis_data, data2 = list(A = A), family = gaussian()
)

# Plot predictions over time with phylogenetic uncertainty

# Model A: Year alone explains weight loss
m_time_only <- brm(
  weight_loss ~ year_std + (1|gr(strain_id, cov = A)), 
  data = analysis_data, data2 = list(A = A), family = gaussian()
)

# Model B: Tropism explains weight loss (independent of year)
m_tropism_only <- brm(
  weight_loss ~ scale(tropism_specificity_day5) + (1|gr(strain_id, cov = A)), 
  data = analysis_data, data2 = list(A = A), family = gaussian()
)

# Model C: Both (does tropism add explanatory power after accounting for time trend?)
m_both <- brm(
  weight_loss ~ scale(tropism_specificity_day5) + year_std + (1|gr(strain_id, cov = A)), 
  data = analysis_data, data2 = list(A = A), family = gaussian()
)

# Compare
loo_compare(
  add_criterion(m_time_only, "loo"),
  add_criterion(m_tropism_only, "loo"),
  add_criterion(m_both, "loo")
)


data_std %>% mutate(tropism_specificity_day5 =( day5_lung_mean - day5_nasal_mean)/(day5_nasal_mean+day5_lung_mean)) %>% ggplot(aes(x = isolation_year, y = tropism_specificity_day5)) + geom_point()