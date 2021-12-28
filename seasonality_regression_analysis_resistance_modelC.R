#This script runs the non-linear regressions to fit seasonal antibiotic resistance data
#to a sinusoidal model with a linear correction and adjustment for age sex and site of
#infection (Model C)

#Load libraries
library(tidyverse)
library(magrittr)
library(plotrix)

# ##################################################
# Inputs
# ##################################################

#Load functions
source("seasonality_regression_functions.R")

#Load raw resistance data
data.SA = read_csv("raw_data/Saureus_antibiotic_resistance_data.csv")
data.EC = read_csv("raw_data/Ecoli_antibiotic_resistance_data.csv")
data.KP = read_csv("raw_data/Kpneumoniae_antibiotic_resistance_data.csv")

#Combine raw data and edit demographics columns
data.res = bind_rows(data.SA, data.EC, data.KP) %>%
  mutate(is_male = ifelse(sex == "M", 1, 0)) %>%
  mutate(is_SST = ifelse(site_of_infection == "skin_softtissue", 1, 0)) %>%
  mutate(is_BL = ifelse(site_of_infection == "blood", 1, 0)) %>%
  mutate(is_RT = ifelse(site_of_infection == "respiratory_tract", 1, 0)) %>%
  mutate(is_AB = ifelse(site_of_infection == "abscess_or_fluid_nos", 1, 0)) %>%
  select(-sex, -site_of_infection)

# ##################################################
# Run regressions
# ##################################################

#Run resistance seasonality regressions for each org/drug with 6 and 12m periods using Model C
results = data.res %>%
  #make hospital/year column
  mutate(hos_year = paste(hospital, as.character(year), sep = "_")) %>%
  #define the outcome variable (y): log2(MIC)
  mutate(y = log2(MIC)) %>%
  #nest dataframe by organism and drug
  nest(-organism, -drug_code, -drug_name, -drug_class) %>%
  #define 2 models (12-month and 6-month period) to run for each organism/drug combination
  left_join(
    crossing(drug_code = c("AMC", "AMP", "CIP", "ERY", "NIT", "OXA", "PEN", "TET"), period = c(6, 12)),
    by = c("drug_code")
  ) %>%
  mutate(omega = 2*pi/period) %>%
  #create a model matrix for each model
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  #run the regression
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, on="hos_year"), .f = model_C_func)) %>%
  #calculate the AIC for each model
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  #extract model parameters into a table
  mutate(model_summary = pmap(.l = list(model = model, model_matrix = model_matrix, on = "hos_year", num_covariates = 6), .f = get_model_summary_func)) 

# ##################################################
# Make parameters table
# ##################################################

#Edit use and resistance model parameters so that all sinusoid amplitudes are positive
#and phases are between 0 and 12 (or 6, depending on the period) months
#Make tables of raw model parameters
model.params.raw = results %>%
  select(organism, drug_code, drug_name, drug_class, period, omega, AIC, model_summary) %>% 
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
model.params.edit = model.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-hos_year) %>%
  gather(variable, value, -(c("organism", "drug_code", "drug_name", "drug_class", "term", "omega", "period", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("organism", "drug_code", "drug_name", "drug_class", "omega", "period", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
model.params.full = model.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(model.params.edit) %>%
  select(organism, drug_code, drug_name, drug_class, period, omega, AIC, term, hos_year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "beta_age", "beta_sex", "beta_BL", "beta_RT", "beta_SST", "beta_AB", "slope", "intercept"))) %>%
  #apply Benjamini-Hochberg corrections to amplitude p-values
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  arrange(organism, drug_code, period, term)

# ##################################################
# Make seasonal deviates table
# ##################################################

#Calculate seasonal deviates
deviates = results %>%
  mutate(deviates = pmap(
    .l = list(data = data, model_summary = model_summary),
    .f = function(data, model_summary) {
      
      #get beta demographics
      B_age = model_summary %>%
        filter(term == "beta_age") %>%
        pull(estimate)
      
      B_sex = model_summary %>%
        filter(term == "beta_sex")%>%
        pull(estimate)
      
      B_SST = model_summary %>%
        filter(term == "beta_SST") %>%
        pull(estimate)
      
      B_AB = model_summary %>%
        filter(term == "beta_AB")%>%
        pull(estimate)
      
      B_BL = model_summary %>%
        filter(term == "beta_RT") %>%
        pull(estimate)
      
      B_RT = model_summary %>%
        filter(term == "beta_BL")%>%
        pull(estimate)
      
      #make table of slopes and intercepts
      slopes_intercepts = model_summary %>%
        filter(term %in% c("slope", "intercept")) %>%
        select(term, estimate, hos_year) %>%
        spread(term,estimate)
      
      #make table of seasonal deviates
      deviates = data %>%
        #add slopes and intercepts
        left_join(slopes_intercepts, by = c("hos_year")) %>%
        #add beta_age and beta_sex
        mutate(beta_age = B_age) %>%
        mutate(beta_sex = B_sex) %>%
        mutate(beta_SST = B_SST) %>%
        mutate(beta_AB = B_AB) %>%
        mutate(beta_BL = B_BL) %>%
        mutate(beta_RT = B_RT) %>%
        #calculate de-trended seasonal deviate 
        mutate(deviate = y - (slope*month + intercept + beta_age*age + beta_sex*is_male + beta_SST*is_SST + beta_AB*is_AB + beta_BL*is_BL + beta_RT*is_RT)) %>%
        #calculate mean seasonal deviates by month
        group_by(month) %>%
        summarize(seasonal_deviate = mean(deviate), sem = std.error(deviate)) %>%
        ungroup()
      
      return(deviates)
    }
  )) 

#Print deviates table
deviates.table = deviates %>%
  select(organism, drug_code, drug_name, drug_class, omega, period, AIC, deviates) %>%
  unnest(deviates) %>%
  arrange(organism, drug_code, period, month)

# ##################################################
# Save outputs
# ##################################################

write_csv(model.params.full, "tables/resistance_modelC_values.csv")
write_csv(deviates.table, "tables/resistance_modelC_seasonal_deviates.csv")