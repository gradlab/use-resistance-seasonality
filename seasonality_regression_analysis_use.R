#This script runs the non-linear regressions to fit seasonal antibiotic use data to a sinusoidal model with a linear correction (Model A)

#Load libraries
library(tidyverse)
library(magrittr)
library(plotrix)

# ##################################################
# Inputs
# ##################################################

#Load functions
source("seasonality_regression_functions.R")

#Load use dataset
data.use = read_csv("raw_data/antibiotic_use_data.csv")

# ##################################################
# Run regressions
# ##################################################

#Run use regressions for each drug class using both a 12-month and 6-month period 
results = data.use %>%
  #convert years to characters
  mutate(year = as.character(year)) %>%
  #define the outcome variable (y)
  mutate(y = mean_daily_claims_per_10000ppl) %>%
  #nest dataframe by drug class
  nest(-drug_class) %>%
  #define 2 models (12-month and 6-month period) to run for each drug class 
  left_join(crossing(drug_class = c("Macrolides", "Nitrofurans", "Penicillins", "Quinolones", "Tetracyclines"), period = c(6, 12)), by = c("drug_class")) %>%
  mutate(omega = 2*pi/period) %>%
  #create a model matrix for each model
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "year"))) %>%
  #run the regression
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, on="year"), .f = model_A_func)) %>%
  #calculate the AIC for each model 
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  #extract model parameters into a table
  mutate(model_summary = pmap(.l = list(model = model, model_matrix = model_matrix, on = "year", num_covariates = 0), .f = get_model_summary_func)) 

# ##################################################
# Make parameters table
# ##################################################

#Edit use and resistance model parameters so that all amplitudes are positive
#and phases are between 0 and 12 (or 6, depending on the period) months
#Make tables of raw model parameters
model.params.raw = results %>%
  select(drug_class, period, omega, AIC, model_summary) %>%
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
model.params.edit = model.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-year) %>%
  gather(variable, value, -(c("drug_class", "term", "period", "omega", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("drug_class", "period", "omega", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
model.params.full = model.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(model.params.edit) %>%
  select(drug_class, period, omega, AIC, term, year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  #apply Benjamini-Hochberg corrections to p-values
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  arrange(drug_class, period, term)

# ##################################################
# Make seasonal deviates table
# ##################################################

#Calculate seasonal deviates
deviates = results %>%
  mutate(deviates = pmap(
    .l = list(data = data, model_summary = model_summary),
    .f = function(data, model_summary) {
      #make table of slopes and intercepts
      slopes_intercepts = model_summary %>%
        filter(term %in% c("slope", "intercept")) %>%
        select(term, estimate, year) %>%
        spread(term,estimate)
      
      #make table of seasonal deviates
      deviates = data %>%
        #add slopes and intercepts
        left_join(slopes_intercepts, by = c("year")) %>%
        #calculate de-trended seasonal deviate 
        mutate(deviate = y - (slope*month + intercept)) %>%
        #calculate mean seasonal deviates by month
        group_by(month) %>%
        summarize(seasonal_deviate = mean(deviate), sem = std.error(deviate)) %>%
        ungroup()
      return(deviates)
    }
  )) 

#Print deviates table
deviates.table = deviates %>%
  select(drug_class, omega, period, AIC, deviates) %>%
  unnest(deviates) %>%
  arrange(drug_class, period, month)

# ##################################################
# Save outputs
# ##################################################

write_csv(model.params.full, "tables/use_model_values.csv")
write_csv(deviates.table, "tables/use_seasonal_deviates.csv")
