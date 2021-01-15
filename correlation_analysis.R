#This script calculates Spearman correlations between antibiotic use and resistance seasonal deviates.

#Load libraries
library(tidyverse)
library(magrittr)
library(DescTools)

# ######################################################
# Inputs
# ######################################################

#Load use and resistance model values
use.model.params = read_csv("tables/model_values_use.csv")
res.model.params = read_csv("tables/model_values_resistance.csv")

#Load use and resistance seasonal deviates
deviates.use = read_csv("tables/seasonal_deviates_use.csv")
deviates.res = read_csv("tables/seasonal_deviates_resistance.csv")

# ######################################################
# Calculate correlations
# ######################################################

#Make combined use-resistance seasonal deviates table
use.res.deviates = left_join(
  deviates.res %>%
    #manually filter out org/drug combinations for which there was no seasonality in resistance
    filter(!(organism == "S. aureus" & drug_code %in% c("PEN", "TET"))) %>%
    filter(!(organism == "E. coli" & drug_code == "TET")) %>%
    filter(!(organism == "K. pneumoniae" & drug_code %in% c("AMC", "TET"))) %>%
    select(organism, month, drug_code, drug_name, drug_class, res = seasonal_deviate),
  
  deviates.use %>%
    select(use_class = drug_class, month, use = seasonal_deviate),
  
  by = c("month")
)

#Calculate spearman correlations between use and resistance seasonal deviates
correlations = use.res.deviates %>%
  nest(-organism, -drug_code, -drug_name, -drug_class, -use_class) %>%
  mutate(spearman_test = map(data, ~ cor.test(.$res, .$use, method = "spearman", use = "complete.obs", exact = TRUE))) %>%
  mutate(rho = map_dbl(spearman_test, ~ .$estimate),
         p.value = map_dbl(spearman_test, ~ .$p.value)
  ) %>%
  mutate(n_samples = map_dbl(data, ~nrow(.))) %>%
  mutate(ci = map2(rho, n_samples, ~ CorCI(.x, .y, alternative = "two.sided")),
         ci.lower = map_dbl(ci, ~ .["lwr.ci"]),
         ci.upper = map_dbl(ci, ~ .["upr.ci"])
  ) %>%
  mutate(p.value.BH = p.adjust(p.value, method = "BH")) %>%
  mutate(use_class_short = substr(use_class, 1,3)) %>%
  mutate(use_res = paste(drug_code, use_class_short, sep=" / ")) %>%
  select(organism, res_drug_code = drug_code, res_drug_name = drug_name, res_drug_class = drug_class, use_class,
         rho, ci.lower, ci.upper, p.value, p.value.BH)

#Save correlations table
write_csv(correlations, "tables/correlations.csv")
