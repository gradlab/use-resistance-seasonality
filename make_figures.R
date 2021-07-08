#This script makes all the figures included in the publication.

#Load libraries
library(tidyverse)
library(magrittr)
library(cowplot)
library(ggpubr)

# ######################################################
# Inputs
# ######################################################

#Load antibiotic use data
data.use = read_csv("raw_data/antibiotic_use_data.csv")

#Load model values
use.model.params = read_csv("tables/model_values_use.csv")
res.model.params = read_csv("tables/model_values_resistance.csv")
EC.12m.model.params = read_csv("tables/Ecoli_AMC_AMP_12m_model_values.csv")
res.model.params.under65 = read_csv("tables/model_values_resistance_under65.csv")

#Load multiple testing results for amplitude estimates
use.amplitude.p = read_csv("tables/model_amplitude_pvalues_use.csv")
res.amplitude.p = read_csv("tables/model_amplitude_pvalues_resistance.csv")

#Load use and resistance seasonal deviates
use.deviates = read_csv("tables/seasonal_deviates_use.csv")
res.deviates = read_csv("tables/seasonal_deviates_resistance.csv")
EC.12m.deviates = read_csv("tables/Ecoli_AMC_AMP_12m_seasonal_deviates.csv")

#Load use-resistance correlations results
correlations = read_csv("tables/correlations.csv")

#Define colors to be used across figures
colors = setNames( c("#220050", "#b30059","#0091a8","#359023", "#ffa500"), 
                   c("Macrolides", "Nitrofurans", "Penicillins", "Quinolones", "Tetracyclines") )

# ######################################################
# Make Figure 1
# ######################################################

#Make figure 1a
labels = data.frame(drug_class = c("Penicillins", "Macrolides", "Quinolones", "Tetracyclines", "Nitrofurans"),
                    x.pos = c(4, 4, 4, 4, 4),
                    y.pos = c(6.6, 3.3, 2.4, 1.5, 0.7))

f1a = data.use %>%
  mutate(year_month = paste(as.character(year), as.character(str_pad(month, 2, pad = "0")), sep = "-")) %>%
  mutate(x = dense_rank(year_month)) %>%
  ggplot(aes(x=x, y=claims_per_10000ppl_per_day, group=drug_class, color=drug_class)) +
  geom_point(size = 0.7) +
  geom_smooth(aes(fill=drug_class), span = 0.2, size = 0.7) +
  geom_text(data = labels, aes(x=x.pos, y=y.pos, label=drug_class), hjust = 0, size = 4.5) +
  ylab("Mean daily claims/10,000 people") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(1, 13, 25, 37, 49), labels = c("Jan '11", "Jan '12", "Jan '13", "Jan '14", "Jan '15")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11)) 


#Make figure 1b
#Cosine function
cos_func = function(month, amplitude, phase, omega, intercept) {
  amplitude * cos(omega *(month - phase)) + intercept
} 

#Function to plot seasonal use model
plot_use_model_func = function(deviates, class, amplitude, phase, omega, a_lower, a_upper, sig) {
  
  col = colors[class]
  
  if (sig) { 
    title = paste(class, "*")
  } else {
    title = class
  }
  
  ci = data.frame(month=seq(1,12,0.01)) %>%
    mutate(lower_ci = map_dbl(month, ~cos_func(., a_lower, phase, omega, 0))) %>%
    mutate(upper_ci = map_dbl(month, ~cos_func(., a_upper, phase, omega, 0))) 
  
  p = ggplot(data = deviates, aes(x = month)) +
    geom_point(aes(x = month, y = seasonal_deviate), color = col, size = 1) +
    geom_errorbar(aes(x = month, ymin = seasonal_deviate - sem, ymax = seasonal_deviate + sem), width = 0.5, color = col) + 
    stat_function(fun = cos_func, args = list(a = amplitude, phase = phase, omega = omega, intercept = 0), size = 0.7, color = col) +
    geom_ribbon(data = ci, aes(x = month, ymin = lower_ci, ymax = upper_ci), fill = col, alpha = 0.3) +
    scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11)) +
    ggtitle(title) +
    xlab("Month") +
    theme_classic() +
    theme(legend.position="none",
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank()
    )
  
  return(p)
}

f1b_plots = use.model.params %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(drug_class, omega, term, estimate, ci.lower, ci.upper) %>%
  gather(variable, value, -(c("drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  left_join(use.amplitude.p %>% select(drug_class, amplitude_p.value.BH = p.value.BH), by = c("drug_class")) %>%
  mutate(sig = amplitude_p.value.BH < 0.05) %>%
  #add seasonal deviates table
  left_join(
    use.deviates %>%
      nest(-drug_class) %>%
      rename(deviates_table = data),
    by = c("drug_class")
  ) %>%
  #make plots
  mutate(plot = pmap(.l = list(deviates = deviates_table, class = drug_class, amplitude = amplitude_estimate,
                               phase = phase_estimate, omega = omega, a_lower = amplitude_ci.lower,
                               a_upper = amplitude_ci.upper, sig = sig),
                     .f = plot_use_model_func)) %>%
  pull(plot)

f1b = do.call(plot_grid, c(f1b_plots, nrow = 3, ncol = 2, align = "hv")) %>%
  annotate_figure(left = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 90)) 


#Combine figures 1a and 1b
f1 = plot_grid(plot_grid(f1a, nrow = 2, ncol = 1, rel_heights = c(2, 1)),
               f1b, nrow = 1, ncol = 2, rel_widths = c(1.2, 1), labels = c("A", "B"))


# ######################################################
# Make Figures 2, S1, S2
# ######################################################

#Function to plot use and resistance models together
plot_use_resistance_func = function(deviates, drug, class, u_a, u_p, u_o, u_low, u_up, r_a, r_p, r_o, r_low, r_up, ratio, r_sig) {
  col = colors[class]
  
  if (r_sig) { 
    title = paste(drug, "*")
  } else {
    title = drug
  }
  
  regressions = data.frame(month=seq(1,12,0.01)) %>%
    mutate(r_actual = map_dbl(month, ~cos_func(., r_a, r_p, r_o, 0))) %>%
    mutate(u_actual = map_dbl(month, ~cos_func(., u_a/ratio, u_p, u_o, 0))) %>%
    gather(type, value, -month) %>%
    mutate(leg = case_when(type == "r_actual" ~ paste(drug, "Resistance"), type == "u_actual" ~ paste(class, "Use"))) %>%
    mutate(leg = factor(leg, levels = c(paste(drug, "Resistance"), paste(class, "Use"))))
  
  ci = data.frame(month=seq(1,12,0.01)) %>%
    mutate(r_lower = map_dbl(month, ~cos_func(., r_low, r_p, r_o, 0))) %>%
    mutate(r_upper = map_dbl(month, ~cos_func(., r_up, r_p, r_o, 0))) %>%
    mutate(u_lower = map_dbl(month, ~cos_func(., u_low, u_p, u_o, 0))) %>%
    mutate(u_upper = map_dbl(month, ~cos_func(., u_up, u_p, u_o, 0)))
  
  p = ggplot(data = regressions) +
    geom_point(data = deviates, aes(x = month, y = seasonal_deviate), color = col, size = 1) +
    geom_errorbar(data = deviates, aes(x = month, ymin = seasonal_deviate - sem, ymax = seasonal_deviate + sem), width = 0.5, color = col) +
    geom_line(aes(x = month, y = value, color = leg, linetype = leg), size = 0.7) +
    geom_ribbon(data = ci, aes(x = month, ymin = r_lower, ymax = r_upper), fill = col, alpha = 0.3) +
    geom_ribbon(data = ci, aes(x = month, ymin = u_upper/ratio, ymax = u_lower/ratio), fill = "grey20", alpha = 0.3) +
    scale_color_manual(values = c(col, "grey20")) +
    scale_y_continuous(sec.axis = sec_axis(~. * ratio), limits = c(-.165, .165)) +
    scale_x_continuous(breaks=c(1, 3, 5, 7, 9, 11)) +
    ggtitle(title) +
    theme_classic() +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10),
          axis.title = element_blank()
    ) 
  
  return(p)
} 


use_res_plots = res.model.params %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(organism, drug_name, drug_class, omega, term, estimate, ci.lower, ci.upper, p.value) %>%
  gather(variable, value, -(c("organism","drug_name","drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  left_join(res.amplitude.p %>% select(organism, drug_name, amplitude_p.value.BH = p.value.BH), by = c("organism", "drug_name")) %>%
  mutate(res_sig = amplitude_p.value.BH < 0.05) %>%
  #add resistance seasonal deviates table
  left_join(
    res.deviates %>%
      nest(-organism, -drug_name, -drug_class) %>%
      rename(deviates_table = data),
    by = c("organism","drug_name","drug_class")
  ) %>%
  select(organism, drug_name, drug_class, res_amplitude = amplitude_estimate, res_phase = phase_estimate,
         res_omega = omega, res_upper = amplitude_ci.upper, res_lower = amplitude_ci.lower, res_sig, deviates_table) %>%
  #add use model params
  left_join(
    use.model.params %>%
      filter(term %in% c("amplitude", "phase")) %>%
      select(drug_class, omega, term, estimate, ci.lower, ci.upper) %>%
      gather(variable, value, -(c("drug_class", "term", "omega"))) %>%
      unite(temp, term, variable) %>%
      spread(temp, value) %>%
      select(drug_class, use_amplitude = amplitude_estimate, use_phase = phase_estimate, use_omega = omega,
             use_upper = amplitude_ci.upper, use_lower = amplitude_ci.lower),
    by = c("drug_class")
  ) %>%
  #get use-resistance scaling ratio
  mutate(u.r_ratio = abs(use_amplitude/res_amplitude)) %>%
  group_by(drug_class) %>%
  mutate(ratio = min(u.r_ratio)) %>%
  ungroup() %>%
  #make plots
  mutate(drug_name = ifelse(drug_name == "Amoxicillin/Clavulanate", "Amox/Clav", drug_name)) %>%
  mutate(plot = pmap(.l = list(deviates = deviates_table, drug = drug_name, class = drug_class, u_a = use_amplitude,
                               u_p = use_phase, u_o = use_omega, u_low = use_lower, u_up = use_upper,
                               r_a = res_amplitude, r_p = res_phase, r_o = res_omega, r_low = res_lower,
                               r_up = res_upper, ratio = ratio, r_sig = res_sig),
                     .f = plot_use_resistance_func)) 

#Make figure 2
f2_plots = use_res_plots %>%
  filter(organism == "S. aureus") %>%
  arrange(drug_class, drug_name) %>%
  pull(plot)

f2 = do.call(plot_grid, c(f2_plots, nrow = 2, ncol = 3)) %>%
  annotate_figure(left = text_grob("Seasonal deviates in resistance,\nadjusted for age and sex (log2(MIC))", size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270)) %>%
  annotate_figure(bottom = text_grob("Month", size = 11))

#Make figure S1
fS1_plots = use_res_plots %>%
  filter(organism == "E. coli") %>%
  arrange(drug_class, drug_name) %>%
  pull(plot)

fS1 = do.call(plot_grid, c(fS1_plots, nrow = 2, ncol = 3)) %>%
  annotate_figure(left = text_grob("Seasonal deviates in resistance,\nadjusted for age and sex (log2(MIC))", size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270)) %>%
  annotate_figure(bottom = text_grob("Month", size = 11))

#Make figure S2
fS2_plots = use_res_plots %>%
  filter(organism == "K. pneumoniae") %>%
  arrange(drug_class, drug_name) %>%
  pull(plot)

fS2 = do.call(plot_grid, c(fS2_plots, nrow = 2, ncol = 2)) %>%
  annotate_figure(left = text_grob("Seasonal deviates in resistance,\nadjusted for age and sex (log2(MIC))", size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270)) %>%
  annotate_figure(bottom = text_grob("Month", size = 11))

# ######################################################
# Make figure S3 
# ######################################################

#Plot E. coli AMC and AMP resistance with a 12-month period
fS3_plots = EC.12m.model.params %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(organism, drug_name, drug_class, omega, term, estimate, ci.lower, ci.upper, p.value) %>%
  gather(variable, value, -(c("organism","drug_name","drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(res_sig = amplitude_p.value < 0.05) %>%
  #add resistance seasonal deviates table
  left_join(
    EC.12m.deviates %>%
      nest(-organism, -drug_name, -drug_class) %>%
      rename(deviates_table = data),
    by = c("organism","drug_name","drug_class")
  ) %>%
  select(organism, drug_name, drug_class, res_amplitude = amplitude_estimate, res_phase = phase_estimate,
         res_omega = omega, res_upper = amplitude_ci.upper, res_lower = amplitude_ci.lower, deviates_table, res_sig) %>%
  #add use model params
  left_join(
    use.model.params %>%
      filter(term %in% c("amplitude", "phase")) %>%
      select(drug_class, omega, term, estimate, ci.lower, ci.upper) %>%
      gather(variable, value, -(c("drug_class", "term", "omega"))) %>%
      unite(temp, term, variable) %>%
      spread(temp, value) %>%
      select(drug_class, use_amplitude = amplitude_estimate, use_phase = phase_estimate, use_omega = omega,
             use_upper = amplitude_ci.upper, use_lower = amplitude_ci.lower),
    by = c("drug_class")
  ) %>%
  #get use-resistance scaling ratio
  mutate(u.r_ratio = abs(use_amplitude/res_amplitude)) %>%
  group_by(drug_class) %>%
  mutate(ratio = min(u.r_ratio)) %>%
  ungroup() %>%
  #make plots
  mutate(drug_name = ifelse(drug_name == "Amoxicillin/Clavulanate", "Amox/Clav", drug_name)) %>%
  mutate(plot = pmap(.l = list(deviates = deviates_table, drug = drug_name, class = drug_class, u_a = use_amplitude,
                               u_p = use_phase, u_o = use_omega, u_low = use_lower, u_up = use_upper,
                               r_a = res_amplitude, r_p = res_phase, r_o = res_omega, r_low = res_lower,
                               r_up = res_upper, ratio = ratio, r_sig = res_sig),
                     .f = plot_use_resistance_func)) %>%
  pull(plot)

#Make figure s3
fS3 = do.call(plot_grid, c(fS3_plots, nrow = 1, ncol = 2)) %>%
  annotate_figure(left = text_grob("Seasonal deviates in resistance,\nadjusted for age and sex (log2(MIC))", size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use\n(mean daily claims/10,000 people)", size = 11, rot = 270)) %>%
  annotate_figure(bottom = text_grob("Month", size = 11))

# ######################################################
# Make figure 3
# ######################################################

#Make amplitude plots function
plot_amplitudes_func = function(dat, title) {
  p = ggplot(data = dat, aes(x = reorder(drug_code, estimate), y = estimate, color = drug_class)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.4) +
    geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
    scale_y_continuous(breaks = c(0, 0.05, 0.1), limits = c(-0.02, 0.1)) + #CHANGED lower limits from -0.01 to -0.02 for under 65 analysis
    scale_color_manual(values = colors, name = "Antibiotic Class") +
    coord_flip() +
    ggtitle(title) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 12),
          axis.title = element_blank(),
          axis.text = element_text(size = 11),
          legend.title = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 11)
    )
  
  return(p)
}

f3_SA_plot = res.model.params %>%
  filter(term == "amplitude") %>%
  filter(organism == "S. aureus") %>%
  plot_amplitudes_func(., "S. aureus")

f3_EC_plot = res.model.params %>%
  filter(term == "amplitude") %>%
  filter(organism == "E. coli") %>%
  plot_amplitudes_func(., "E. coli")

f3_KP_plot = res.model.params %>%
  filter(term == "amplitude") %>%
  filter(organism == "K. pneumoniae") %>%
  plot_amplitudes_func(., "K. pneumoniae")

#Make figure 3
f3 = ggarrange(f3_SA_plot, f3_EC_plot, f3_KP_plot, nrow = 1, ncol = 3, common.legend = TRUE, legend = "right") %>%
  annotate_figure(bottom = text_grob("Amplitude of seasonality (log2(MIC) deviates)", size = 11))

# ######################################################
# Make figure 4
# ######################################################

#Make resistance phases table
res.phases = bind_rows(
  res.model.params %>%
    filter(term == "phase"),
  #add second peak for org/drugs with 6 month period
  res.model.params %>%
    filter(term == "phase") %>%
    filter(period == 6) %>%
    mutate(estimate = estimate + 6, ci.lower = ci.lower + 6, ci.upper = ci.upper + 6)
) %>%
  #filter out org/drug combinations that did not meet criteria for seasonality
  filter(!(organism == "S. aureus" & drug_code %in% c("PEN", "TET"))) %>%
  filter(!(organism == "E. coli" & drug_code == "TET")) %>%
  filter(!(organism == "K. pneumoniae" & drug_code %in% c("AMC", "TET"))) %>%
  #manually edit phase for K.pneumo/NIT (subtract 12m) for aesthetic purposes for figure
  mutate(estimate = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", estimate - 12, estimate),
         ci.lower = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", ci.lower - 12, ci.lower),
         ci.upper = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", ci.upper - 12, ci.upper)
  ) %>%
  mutate(org_drug = paste(organism, drug_code, sep = "/")) 

#Make use phases table
use.phases = bind_rows(
  use.model.params %>%
    filter(term == "phase"),
  #add second peak for org/drugs with 6 month period
  use.model.params %>%
    filter(term == "phase") %>%
    filter(period == 6) %>%
    mutate(estimate = estimate + 6, ci.lower = ci.lower + 6, ci.upper = ci.upper + 6)
) %>%
  filter(drug_class != "Tetracyclines")

#Define plot labels
labs = c("S. aureus/CIP" = expression(paste(italic("S. aureus"), " / CIP")),
         "S. aureus/ERY" = expression(paste(italic("S. aureus"), " / ERY")),
         "S. aureus/NIT" = expression(paste(italic("S. aureus"), " / NIT")),
         "S. aureus/OXA" = expression(paste(italic("S. aureus"), " / OXA")),
         "E. coli/AMC" = expression(paste(italic("E. coli"), " / AMC")),
         "E. coli/AMP" = expression(paste(italic("E. coli"), " / AMP")),
         "E. coli/CIP" = expression(paste(italic("E. coli"), " / CIP")),
         "E. coli/NIT" = expression(paste(italic("E. coli"), " / NIT")),
         "K. pneumoniae/CIP" = expression(paste(italic("K. pneumoniae"), " / CIP")),
         "K. pneumoniae/NIT" = expression(paste(italic("K. pneumoniae"), " / NIT"))
)


#Make figure 4
f4 = ggplot(data = res.phases) +
  facet_wrap(~drug_class, scales = "free", nrow = 5) +
  geom_point(aes(x = reorder(org_drug, estimate), y = estimate, color=drug_class), size = 2) +
  geom_errorbar(aes(x = reorder(org_drug, estimate), ymin = ci.lower, ymax = ci.upper, color=drug_class), width = 0.4) +
  geom_hline(data = use.phases, aes(yintercept = estimate, color = drug_class), size = 1) +
  geom_rect(data = use.phases, aes(ymin = ci.lower, ymax = ci.upper, xmin = -Inf, xmax = Inf, fill = drug_class), alpha = 0.3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = labs) +
  scale_y_continuous(breaks = seq(-2 ,11, 1), limits = c(-2.5, 11), #Changed -2.2 to -2.5
                     labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov")) +
  coord_flip() +
  ylab("Peak month(s) of seasonality") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank()
  )

# ######################################################
# Make figures 5, S4
# ######################################################

#Function to plot use-resistance correlations for each species 
#with lags of 0-3 months
plot_corr_func = function(dat, o) {
  d0 = dat %>%
    filter(organism == o) %>%
    filter(lag == "0 months") %>%
    mutate(use_class_short = substr(use_class, 1,3)) %>%
    mutate(use_res = paste(res_drug_code, use_class_short, sep = "/"))
  
  d1 = dat %>%
    filter(organism == o) %>%
    filter(lag == "1 month") %>%
    mutate(use_class_short = substr(use_class, 1,3)) %>%
    mutate(use_res = paste(res_drug_code, use_class_short, sep = "/"))
  
  d2 = dat %>%
    filter(organism == o) %>%
    filter(lag == "2 months") %>%
    mutate(use_class_short = substr(use_class, 1,3)) %>%
    mutate(use_res = paste(res_drug_code, use_class_short, sep = "/"))
  
  d3 = dat %>%
    filter(organism == o) %>%
    filter(lag == "3 months") %>%
    mutate(use_class_short = substr(use_class, 1,3)) %>%
    mutate(use_res = paste(res_drug_code, use_class_short, sep = "/"))
  
  y_ticks = d0 %>%
    mutate(res_label = paste0(res_drug_code, " /")) %>%
    arrange(rho) %>%
    pull(res_label)
  
  order = d0 %>%
    arrange(rho) %>%
    pull(use_res)
  
  p0 = d0 %>%
    mutate(use_res = factor(use_res, order)) %>%
    ggplot(aes(x = use_res, y = rho)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ylab(expression(paste("Spearman's ", rho))) +
    labs(color = "") +
    ylim(-1, 1) +
    ggtitle("0m lag") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  
  p1= d1 %>%
    mutate(use_res = factor(use_res, order)) %>%
    ggplot(aes(x = use_res, y = rho)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ylab(expression(paste("Spearman's ", rho))) +
    labs(color = "") +
    ylim(-1, 1) +
    ggtitle("1m lag") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  
  p2 = d2 %>%
    mutate(use_res = factor(use_res, order)) %>%
    ggplot(aes(x = use_res, y = rho)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ylab(expression(paste("Spearman's ", rho))) +
    labs(color = "") +
    ylim(-1, 1) +
    ggtitle("2m lag") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  
  p3 = d3 %>%
    mutate(use_res = factor(use_res, order)) %>%
    ggplot(aes(x = use_res, y = rho)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ylab(expression(paste("Spearman's ", rho))) +
    labs(color = "") +
    ylim(-1, 1) +
    ggtitle("3m lag") +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5)
    )
  
  
  labels = ggplot(data = d0, aes(x = reorder(use_res, rho), y = 1, fill = use_class)) +
    geom_tile(color = "white", alpha = 1) +
    geom_text(aes(label = use_class_short), color = "white") + 
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = y_ticks) +
    coord_flip() +
    ggtitle("Resistance / Use") +
    theme_minimal() + 
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.6),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 11),
          axis.title = element_blank()
    )
  
  p = plot_grid(labels, p0, p1, p2, p3, nrow = 1, ncol = 5, rel_widths = c(1, 1.5, 1.5, 1.5, 1.5), align = "h") 
  
  return(p)
}

f5_SA = plot_corr_func(correlations, "S. aureus")
f5_EC = plot_corr_func(correlations, "E. coli") %>% annotate_figure(top = "")
f5_KP = plot_corr_func(correlations, "K. pneumoniae") %>% annotate_figure(top = "")

#Make figure 5
f5 = f5_SA

#Make figure s4
fS4 = plot_grid(f5_EC, f5_KP, ncol = 1, nrow = 2, rel_heights = c(2, 1.3), labels = c("A", "B"))

# ######################################################
# Make table s4
# ######################################################

#Make amplitudes table for under 65 analysis
a_table = res.model.params.under65 %>%
  filter(term == "amplitude") %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ round(., 3)) %>%
  mutate(`amplitude (95% CI)` = paste0(as.character(estimate), " (", as.character(ci.lower), " to ", as.character(ci.upper), ")" )) %>%
  mutate(`amplitude p-value` = signif(p.value, 1)) %>%
  mutate(period = paste(period, "months")) %>%
  select(organism, drug_name, period, `amplitude (95% CI)`, `amplitude p-value`)

#Make phases table for under 65 analysis
p_table = res.model.params %>%
  filter(term == "phase") %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ round(., 1)) %>%
  mutate(`phase (95% CI)` = paste0(as.character(estimate), " (", as.character(ci.lower), " to ", as.character(ci.upper), ")" )) %>%
  select(organism, drug_name, `phase (95% CI)`)

#Make table s4
tableS4 = left_join(a_table, p_table, by = c("organism", "drug_name"))

# ######################################################
# Save all figures 
# ######################################################

ggsave(f1, filename = "figures/Figure1.pdf", width = 7.5, height = 5.5)
ggsave(f2, filename = "figures/Figure2.pdf", width = 7.5, height = 5.5)
ggsave(f3, filename = "figures/Figure3.pdf", width = 6.5, height = 4)
ggsave(f4, filename = "figures/Figure4.pdf", width = 6, height = 5.5)
ggsave(f5, filename = "figures/Figure5.pdf", width = 7.5, height = 5)

ggsave(fS1, filename = "figures/SupplementaryFigure1.pdf", width = 7.5, height = 5.5)
ggsave(fS2, filename = "figures/SupplementaryFigure2.pdf", width = 5, height = 5.5)
ggsave(fS3, filename = "figures/SupplementaryFigure3.pdf", width = 5, height = 2.75)
ggsave(fS4, filename = "figures/SupplementaryFigure4.pdf", width = 7.5, height = 8)
write_csv(tableS4, "figures/SupplementaryTable4.csv")