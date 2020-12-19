#This script runs the non-linear regressions to fit the seasonal antibiotic use and resistance data to sinusoidal models. 
#Due to data sharing restrictions, we cannot provide the input files needed to run this code. 
#However, we do provide the output files from this code in the "data/" folder, which can be used to reproduce the figures in the publication.

#Load libraries 
library(tidyverse)
library(magrittr)
library(plotrix)

# ######################################################
# Inputs
# ######################################################

#Load raw use dataset (unpublished)
load("data_use_BOS.RData") #data.use.BOS

#Load raw resistance datasets for each organism (unpublished)
load("data_SA.RData") #data.SA
load("data_EC.RData") #data.EC
load("data_KP.RData") #data.KP

# ######################################################
# Regression functions
# ######################################################

#Function to create a model matrix: where rows represent data rows and columns represent
#each year or clinic/year. The model matrix consists of 0s and 1s, where a 1 means that this 
#data row is from the corresponding year or clinic/year. The model matrix is used in the 
#nls regression to allow for a separate slope and intercept term to be estimated for each
#year or clinic/year.
make_model_matrix_func = function(data, on) {
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #create model matrix - same number of rows as data rows, 0/1s, 1 if this row in clinic year
  model_matrix = model.matrix(formula(paste("~ 0 +", on)), data) %>%
    set_colnames(str_remove(colnames(.), on))
  
  #check that num rows in model matrix same as data
  stopifnot(dim(model_matrix) == c(nrow(data), n_cys))
  
  return(model_matrix)
}

#Function to run the nls regression: sinusoidal + linear fit function
run_regression_func = function(data, model_matrix, omega, A_init, on) {
  
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #get first guess for clinic per year intercepts: use mean values
  start_intercepts = data %>%
    group_by(!!sym(on)) %>%
    summarize(y = mean(y)) %>%
    arrange(!!sym(on)) %T>%
    #check that order matches model matrix columns
    {stopifnot(all(.[[on]] == colnames(model_matrix)))} %>%
    pull(y)
  
  #get first guess for clinic per year slopes: set as 0
  start_slopes = rep(0, n_cys)
  
  #fit the model
  model = nls(y ~ amplitude * cos(omega * (month - phase)) + drop(model_matrix %*% slope) * month + drop(model_matrix %*% intercept),
              start = list(amplitude = A_init, phase = 0, slope = start_slopes, intercept = start_intercepts),
              data = data) 
  
  return(model)
}


#Function to make table of regression output params
get_model_summary_func = function(model, model_matrix, on) {
  
  names = colnames(model_matrix)
  
  model_values1 = summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    set_colnames(c("term","estimate","std.error","statistic","p.value")) 
  
  model_values2 = confint.default(model) %>%
    set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")
  
  model_values = left_join(model_values1, model_values2, by="term") %>%
    #add hospital/year to intercept/slope terms
    mutate(!!on := c("","",names,names)) %>%
    mutate(term = case_when(str_detect(term, "slope") ~ "slope",
                            str_detect(term, "intercept") ~ "intercept",
                            TRUE ~ term))
  
  return(model_values)
}


#Function to make seasonal deviates table 
make_deviates_table_func = function(data, model_summary, on) {
  
  #make table of slopes and intercepts
  slopes_intercepts = model_summary %>%
    filter(term %in% c("slope", "intercept")) %>%
    select(term, estimate, !!sym(on)) %>%
    spread(term,estimate)
  
  #make table of seasonal deviates
  deviates = data %>%
    #add slopes and intercepts
    left_join(slopes_intercepts, by = c(on)) %>%
    #calculate detrended seasonal deviate from slope and intercept
    mutate(deviate = y - (slope*month + intercept)) %>%
    #calculate mean seasonal deviates by month
    group_by(month) %>%
    summarize(seasonal_deviate = mean(deviate), sem = std.error(deviate)) %>%
    ungroup()
  
  return(deviates)
}

# ######################################################
# Run regressions
# ######################################################

#Run use regressions for each drug class using both a 12-month and 6-month period
results.use = data.use.BOS %>%
  #convert years to characters
  mutate(year = as.character(year)) %>%
  #define the outcome variable (y): claims/1000ppl
  mutate(y = rxpkp) %>%
  #nest dataframe by drug class
  nest(-drug_class) %>%
  #define 2 models (12-month and 6-month period) to run for each drug class 
  left_join(crossing(drug_class = c("Macrolides", "Nitrofurans", "Penicillins", "Quinolones", "Tetracyclines"), omega = c(2*pi/6, 2*pi/12)), by = c("drug_class")) %>%
  mutate(model_type = paste0(as.character(2*pi/omega), "m_period")) %>%
  mutate(period = 2*pi/omega) %>%
  #create a model matrix for each model
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "year"))) %>%
  #run the regression
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init=0.1, on="year"), .f = run_regression_func)) %>%
  #calculate the AIC for each model 
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  #extract model parameters into a table
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "year"))) %>%
  #calculate monthly seasonal deviates from each model
  mutate(deviates = pmap(.l = list(data = data, model_summary = model_summary, on = "year"), .f = make_deviates_table_func))

#Run resistance regressions of S. aureus
results.SA = data.SA %>%
  #define the outcome variable (y): log2(MIC)
  mutate(y = as.double(log_value)) %>%
  #nest dataframe by organism and drug
  nest(-organism, -org, -drug_code, -drug_name, -drug_class) %>%
  #define 2 models (12-month and 6-month period) to run for each organism/drug combination 
  left_join(crossing(drug_code = c("CIP", "ERY", "NIT", "OXA", "PEN", "TET"), omega = c(2*pi/6, 2*pi/12)), by = c("drug_code")) %>%
  mutate(model_type = paste0("period_", as.character(2*pi/omega), "m")) %>%
  mutate(period = 2*pi/omega) %>%
  #create a model matrix for each model
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  #run the regression
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init=0.1, on="hos_year"), .f = run_regression_func)) %>%
  #calculate the AIC for each model
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  #extract model parameters into a table
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "hos_year"))) %>%
  #calculate monthly seasonal deviates from each model
  mutate(deviates = pmap(.l = list(data = data, model = model_summary, on = "hos_year"), .f = make_deviates_table_func))

#Run resistance regressions of E. coli 
results.EC = data.EC %>%
  mutate(y = as.double(log_value)) %>%
  nest(-organism, -org, -drug_code, -drug_name, -drug_class) %>%
  left_join(crossing(drug_code = c("AMC", "AMP", "CIP", "NIT", "TET"), omega = c(2*pi/6, 2*pi/12)), by = c("drug_code")) %>%
  mutate(model_type = paste0("period_", as.character(2*pi/omega), "m")) %>%
  mutate(period = 2*pi/omega) %>%
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init = 0.1, on = "hos_year"), .f = run_regression_func)) %>%
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "hos_year"))) %>%
  mutate(deviates = pmap(.l = list(data = data, model_summary = model_summary, on = "hos_year"), .f = make_deviates_table_func))

#Run resistance regressions of K. pneumo 
results.KP = data.KP %>%
  mutate(y = as.double(log_value)) %>%
  nest(-organism, -org, -drug_code, -drug_name, -drug_class) %>%
  left_join(crossing(drug_code = c("AMC", "CIP", "NIT", "TET"), omega = c(2*pi/6, 2*pi/12)), by = c("drug_code")) %>%
  mutate(model_type = paste0("period_", as.character(2*pi/omega), "m")) %>%
  mutate(period = 2*pi/omega) %>%
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init = 0.1, on = "hos_year"), .f = run_regression_func)) %>%
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "hos_year"))) %>%
  mutate(deviates = pmap(.l = list(data = data, model_summary = model_summary, on = "hos_year"), .f = make_deviates_table_func))

# ######################################################
# Make AIC comparison tables (Table S4 and S5)
# ######################################################

#Make use model comparison table (S4)
table.S4 = results.use %>%
  select(drug_class, model_type, AIC) %>%
  group_by(drug_class) %>%
  mutate(min = min(AIC)) %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(drug_class, model_type, AIC_2) %>%
  spread(model_type, AIC_2)

#Make resistance model comparison table (S5)
table.S5 = bind_rows(results.SA, results.EC, results.KP) %>%
  select(org, drug_name, model_type, AIC) %>%
  group_by(org, drug_name) %>%
  mutate(min = min(AIC)) %>%
  ungroup() %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(org, drug_name, model_type, AIC_2) %>%
  spread(model_type, AIC_2)

#Save AIC tables to file
write_csv(table.S4, "figures/Table_S4.csv")
write_csv(table.S5, "figures/Table_S5.csv")

# ######################################################
# Edit then print use and resistance model values to file
# ######################################################

#Function to convert negative amplitudes to positive and phases to range from 0 months to model period. 
#A sinusoid with amplitude -A is the same as a sinusoid with amplitude A and phase + 0.5*period.  
#Similarly, a sinusoid with phase p is the same a sinuoid with phase p-period or p+period. 
convert_a_phases_func = function(a_estimate, a_ci.lower, a_ci.upper, phase_estimate, phase_ci.lower, phase_ci.upper, period) {
  
  #if the amplitude is negative, add (-2*a) to a, a_ci.lower, a_ci.upper; add 1/2*period to phase, phase_ci.lower, phase_ci.upper
  if (a_estimate < 0) {
    a_ci.lower = a_ci.lower + -2*a_estimate
    a_ci.upper = a_ci.upper + -2*a_estimate
    a_estimate = a_estimate + -2*a_estimate
    
    phase_estimate = phase_estimate + period/2
    phase_ci.lower = phase_ci.lower + period/2
    phase_ci.upper = phase_ci.upper + period/2
  }
  
  #while phase is not between 1-period, add or subtract period 
  while (!(phase_estimate >= 0 & phase_estimate <= period)) {
    if(phase_estimate < 0) {
      phase_estimate = phase_estimate + period
      phase_ci.lower = phase_ci.lower + period
      phase_ci.upper = phase_ci.upper + period
    }
    
    if(phase_estimate > period) {
      phase_estimate = phase_estimate - period
      phase_ci.lower = phase_ci.lower - period
      phase_ci.upper = phase_ci.upper - period
    }
  }
  
  dat = data.frame(amplitude_estimate = a_estimate, amplitude_ci.lower = a_ci.lower,
                   amplitude_ci.upper = a_ci.upper, phase_estimate = phase_estimate,
                   phase_ci.lower = phase_ci.lower, phase_ci.upper = phase_ci.upper)
  return(dat)
}

#Make tables of raw model parameters
use.model.params.raw = results.use %>%
  select(drug_class, period, omega, AIC, model_summary) %>%
  unnest(model_summary)

res.model.params.raw = bind_rows(results.SA, results.EC, results.KP) %>%
  select(org, drug_code, drug_name, drug_class, period, omega, AIC, model_summary) %>% 
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
use.model.params.edit = use.model.params.raw %>%
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

res.model.params.edit = res.model.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-hos_year) %>%
  gather(variable, value, -(c("org", "drug_code", "drug_name", "drug_class", "term", "omega", "period", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(period = 2*pi/omega) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("org", "drug_code", "drug_name", "drug_class", "omega", "period", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
use.model.params.full = use.model.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(use.model.params.edit) %>%
  select(drug_class, period, omega, AIC, term, year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  arrange(drug_class, term)

res.model.params.full = res.model.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(res.model.params.edit) %>%
  select(org, drug_code, drug_name, drug_class, period, omega, AIC, term, hos_year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  arrange(org, drug_code, term)

#Filter models to only keep models with the lower AIC for each drug class or org/drug combination
use.model.params = use.model.params.full %>%
  group_by(drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) 

res.model.params = res.model.params.full %>%
  group_by(org, drug_code, drug_name, drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) 

#Make separate table of model params for E. coli/AMC and E. coli/AMP resistance with 12m period
EC.12m.params = res.model.params.full %>%
  filter(org == "E. coli" & drug_code %in% c("AMC", "AMP") & period == 12)

#Save model values to file
write_csv(use.model.params, "data/model_values_use.csv")
write_csv(res.model.params, "data/model_values_resistance.csv")
write_csv(EC.12m.params, "data/Ecoli_AMC_AMP_12m_model_values.csv")

# ######################################################
# Multiple testing correction on amplitude estimates
# ######################################################

#Apply Benjamini-Hochberg corrections to amplitude estimates
use.amplitudes.p_adj = use.model.params %>%
  filter(term == "amplitude") %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  select(-year)

res.amplitudes.p_adj = res.model.params %>%
  filter(term == "amplitude") %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  mutate(significant = p.value.BH < 0.05) %>%
  select(-hos_year) %>%
  arrange(org, drug_code)

#Save amplitude p-values to file
write_csv(use.amplitudes.p_adj, "data/model_amplitude_pvalues_use.csv")
write_csv(res.amplitudes.p_adj, "data/model_amplitude_pvalues_resistance.csv")

# ######################################################
# Make seasonal deviates tables and save to file
# ######################################################

#Make combined table of raw model parameters. Keep only the model with the lower AIC for each org/drug combination.
deviates.use = results.use %>%
  select(drug_class, period, AIC, deviates) %>%
  #filter to just model with lower AIC for each org/drug
  group_by(drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) %>%
  unnest(deviates)

deviates.res = bind_rows(results.SA, results.EC, results.KP) %>%
  select(org, drug_code, drug_name, drug_class, period, AIC, deviates) %>%
  #filter to just model with lower AIC for each org/drug
  group_by(org, drug_code, drug_name, drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) %>%
  unnest(deviates)

deviates.EC.12m = results.EC %>%
  filter(org == "E. coli" & drug_code %in% c("AMC", "AMP") & period == 12) %>%
  select(org, drug_code, drug_name, drug_class, period, deviates) %>%
  unnest(deviates)

#Save seasonal deviates tables to file
write_csv(deviates.use, "data/seasonal_deviates_use.csv")
write_csv(deviates.res, "data/seasonal_deviates_resistance.csv")
write_csv(deviates.EC.12m, "data/Ecoli_AMC_AMP_12m_seasonal_deviates.csv")

