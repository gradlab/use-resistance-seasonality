#This script calculates spearman correlations between use and resistance seasonal deviates

#Load libraries
library(tidyverse)
library(magrittr)
library(DescTools)

# ##################################################
# Inputs
# ##################################################

#Load use and resistance seasonal deviates
deviates.use = read_csv("tables/use_seasonal_deviates.csv")
deviates.res = read_csv("tables/resistance_modelA_seasonal_deviates.csv")

# ##################################################
# Calculate correlations
# ##################################################

#Filter use and resitsance deviates to best-fitting model by AIC between 6 and 12 month period
deviates.use.fil = deviates.use %>%
  group_by(drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank)

deviates.res.fil = deviates.res %>%
  group_by(organism, drug_code, drug_name, drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank)

#Make combined use-resistance seasonal deviates table
use.res.deviates = left_join(
  bind_rows(
    deviates.res.fil %>%
      mutate(join_by = month) %>%
      mutate(lag = "0 months"), 
    
    deviates.res.fil %>%
      mutate(join_by = month - 1) %>%
      mutate(lag = "1 month"),
    
    deviates.res.fil %>%
      mutate(join_by = month - 2) %>%
      mutate(lag = "2 months"),
    
    deviates.res.fil %>%
      mutate(join_by = month - 3) %>%
      mutate(lag = "3 months")
  ) %>%
    mutate(join_by = ifelse(join_by <= 0, join_by+12, join_by)) %>%
    filter(!(organism == "S. aureus" & drug_code %in% c("PEN", "TET"))) %>%
    filter(!(organism == "E. coli" & drug_code %in% c("AMC", "TET"))) %>%
    filter(!(organism == "K. pneumoniae" & drug_code %in% c("AMC", "TET"))) %>%
    
    select(organism, join_by, res_month = month, drug_code, drug_name, drug_class, lag, res = seasonal_deviate),
  
  deviates.use.fil %>%
    select(use_class = drug_class, join_by = month, use_month = month, use = seasonal_deviate),
  
  by = c("join_by")
)

#Calculate spearman correlations between use and resistance seasonal deviates
correlations = use.res.deviates %>%
  nest(-organism, -drug_code, -drug_name, -drug_class, -use_class, -lag) %>%
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
  select(organism, res_drug_code = drug_code, res_drug_name = drug_name, res_drug_class = drug_class, use_class, lag,
         rho, ci.lower, ci.upper, p.value, p.value.BH)

# ##################################################
# Save output
# ##################################################

write_csv(correlations, "tables/correlations.csv")