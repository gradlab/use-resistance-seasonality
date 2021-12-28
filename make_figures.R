#This script makes all the figures and tables included in the publication.

#Load libraries
library(tidyverse)
library(magrittr)
library(ggpubr)

# ##################################################
# Inputs
# ##################################################

#Load antibiotic use data
data.use = read_csv("raw_data/antibiotic_use_data.csv")

#Load antibiotic resistance data
data.SA = read_csv("raw_data/Saureus_antibiotic_resistance_data.csv")
data.EC = read_csv("raw_data/Ecoli_antibiotic_resistance_data.csv")
data.KP = read_csv("raw_data/Kpneumoniae_antibiotic_resistance_data.csv")

#Load model values
use.model.params = read_csv("tables/use_model_values.csv")
res.modelA.params = read_csv("tables/resistance_modelA_values.csv")
res.modelA.params.under65 = read_csv("tables/resistance_OP_under65_modelA_values.csv")
res.modelB.params = read_csv("tables/resistance_modelB_values.csv")
res.modelC.params = read_csv("tables/resistance_modelC_values.csv")

#Load use and resistance (model A) seasonal deviates
use.deviates = read_csv("tables/use_seasonal_deviates.csv")
res.modelA.deviates = read_csv("tables/resistance_modelA_seasonal_deviates.csv")

#Load use-resistance correlations results
correlations = read_csv("tables/correlations.csv")

#Define colors to be used across figures
colors = setNames( c("#220050", "#b30059","#0091a8","#359023", "#ffa500"), 
                   c("Macrolides", "Nitrofurans", "Penicillins", "Quinolones", "Tetracyclines") )

#Filter model parameter tables to the best-fitting model (by AIC)
#between a 6 and 12 month period for each drug class or bug/drug
filter_models_AIC_func = function(table, group_cols) {
  
  table.fil = table %>%
    group_by_at(vars(all_of(group_cols))) %>%
    mutate(rank = dense_rank(AIC)) %>%
    ungroup() %>%
    filter(rank == 1) %>%
    select(-AIC, -rank)
  
  return(table.fil)
}

use.model.params.fil = filter_models_AIC_func(use.model.params, c("drug_class"))
res.modelA.params.fil = filter_models_AIC_func(res.modelA.params, c("organism", "drug_name"))
res.modelA.params.under65.fil = filter_models_AIC_func(res.modelA.params.under65, c("organism", "drug_name"))

# ##################################################
# Make Figure 1
# ##################################################

#Make figure 1a
labels = data.frame(drug_class = c("Penicillins", "Macrolides", "Quinolones", "Tetracyclines", "Nitrofurans"),
                    x.pos = c(4, 4, 4, 4, 4),
                    y.pos = c(6.6, 3.3, 2.4, 1.5, 0.7))

f1a = data.use %>%
  mutate(year_month = paste(as.character(year), as.character(str_pad(month, 2, pad = "0")), sep = "-")) %>%
  mutate(x = dense_rank(year_month)) %>%
  ggplot(aes(x=x, y=mean_daily_claims_per_10000ppl, group=drug_class, color=drug_class)) +
  geom_point(size = 0.7) +
  geom_smooth(aes(fill=drug_class), span = 0.2, size = 0.7) +
  geom_text(data = labels, aes(x=x.pos, y=y.pos, label=drug_class), hjust = 0, size = 4.2) +
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

#Function to plot use model
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
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank()
    )
  
  return(p)
}

f1b_plots = use.model.params.fil %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(drug_class, omega, term, estimate, ci.lower, ci.upper, p.value.BH) %>%
  gather(variable, value, -(c("drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  # left_join(use.amplitude.p %>% select(drug_class, omega, amplitude_p.value.BH = p.value.BH), by = c("drug_class", "omega")) %>%
  mutate(sig = amplitude_p.value.BH < 0.05) %>%
  #add seasonal deviates table
  left_join(
    use.deviates %>%
      nest(-drug_class, -omega) %>%
      rename(deviates_table = data),
    by = c("drug_class", "omega")
  ) %>%
  
  #make plots
  mutate(plot = pmap(.l = list(deviates = deviates_table, class = drug_class, amplitude = amplitude_estimate,
                               phase = phase_estimate, omega = omega, a_lower = amplitude_ci.lower,
                               a_upper = amplitude_ci.upper, sig = sig),
                     .f = plot_use_model_func)) %>%
  pull(plot)

f1b = do.call(ggarrange, c(f1b_plots, nrow = 3, ncol = 2, align = "hv")) %>%
  annotate_figure(left = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 90)) 

#Combine figures 1a and 1b
f1 = ggarrange(ggarrange(f1a, nrow = 2, ncol = 1, heights = c(2, 1)),
               f1b, nrow = 1, ncol = 2, widths = c(1.2, 1), labels = c("A.", "B."))

# ##################################################
# Make Figures 2, S1, S2, S3
# ##################################################

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
    xlab("Month") +
    theme_classic() +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 10),
          axis.title.y = element_blank()
    ) 
  
  return(p)
} 

use_res_plots = res.modelA.params.fil %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(organism, drug_name, drug_class, omega, term, estimate, ci.lower, ci.upper, p.value.BH) %>%
  gather(variable, value, -(c("organism", "drug_name", "drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(res_sig = amplitude_p.value.BH < 0.05) %>%
  #add resistance seasonal deviates table
  left_join(
    res.modelA.deviates %>%
      nest(-organism, -drug_name, -drug_class, -omega) %>%
      rename(deviates_table = data),
    by = c("organism", "drug_name", "drug_class", "omega")
  ) %>%
  select(organism, drug_name, drug_class, res_amplitude = amplitude_estimate, res_phase = phase_estimate,
         res_omega = omega, res_upper = amplitude_ci.upper, res_lower = amplitude_ci.lower, res_sig, deviates_table) %>%
  #add use model params
  left_join(
    use.model.params.fil %>%
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

f2 = do.call(ggarrange, c(f2_plots, nrow = 2, ncol = 3)) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in resistance ("*log["2"]*"(MIC))"), size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270)) 

#Make figure S1
fS1_plots = use_res_plots %>%
  filter(organism == "E. coli") %>%
  arrange(drug_class, drug_name) %>%
  pull(plot)

fS1 = do.call(ggarrange, c(fS1_plots, nrow = 2, ncol = 3)) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in resistance ("*log["2"]*"(MIC))"), size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270))

#Make figure S2
fS2_plots = use_res_plots %>%
  filter(organism == "K. pneumoniae") %>%
  arrange(drug_class, drug_name) %>%
  pull(plot)

fS2 = do.call(ggarrange, c(fS2_plots, nrow = 2, ncol = 2)) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in resistance ("*log["2"]*"(MIC))"), size = 11, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (mean daily claims/10,000 people)", size = 11, rot = 270)) 

#Make figure S3
fS3_plots = res.modelA.params %>%
  filter(organism == "E. coli" & drug_code == "AMP" & period == 12) %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(organism, drug_name, drug_class, omega, term, estimate, ci.lower, ci.upper, p.value.BH) %>%
  gather(variable, value, -(c("organism","drug_name","drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(res_sig = amplitude_p.value.BH < 0.05) %>%
  #add resistance seasonal deviates table
  left_join(
    res.modelA.deviates %>%
      filter(organism == "E. coli" & drug_code == "AMP" & period == 12) %>%
      nest(-organism, -drug_name, -drug_class) %>%
      rename(deviates_table = data),
    by = c("organism","drug_name","drug_class")
  ) %>%
  select(organism, drug_name, drug_class, res_amplitude = amplitude_estimate, res_phase = phase_estimate,
         res_omega = omega, res_upper = amplitude_ci.upper, res_lower = amplitude_ci.lower, deviates_table, res_sig) %>%
  #add use model params
  left_join(
    use.model.params.fil %>%
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
  mutate(plot = pmap(.l = list(deviates = deviates_table, drug = drug_name, class = drug_class, u_a = use_amplitude,
                               u_p = use_phase, u_o = use_omega, u_low = use_lower, u_up = use_upper,
                               r_a = res_amplitude, r_p = res_phase, r_o = res_omega, r_low = res_lower,
                               r_up = res_upper, ratio = ratio, r_sig = res_sig),
                     .f = plot_use_resistance_func)) %>%
  pull(plot)

fS3 = do.call(ggarrange, c(fS3_plots, nrow = 1, ncol = 1)) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in resistance ("*log["2"]*"(MIC))"), size = 10, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use\n(mean daily claims/10,000 people)", size = 10, rot = 270))

# ##################################################
# Make Figure 3
# ##################################################

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
    theme(plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 11),
          axis.title = element_blank(),
          axis.text = element_text(size = 11),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10)
    )
  
  return(p)
}

f3_SA_plot = res.modelA.params.fil %>%
  filter(term == "amplitude") %>%
  filter(organism == "S. aureus") %>%
  plot_amplitudes_func(., "S. aureus")

f3_EC_plot = res.modelA.params.fil %>%
  filter(term == "amplitude") %>%
  filter(organism == "E. coli") %>%
  plot_amplitudes_func(., "E. coli")

f3_KP_plot = res.modelA.params.fil %>%
  filter(term == "amplitude") %>%
  filter(organism == "K. pneumoniae") %>%
  plot_amplitudes_func(., "K. pneumoniae")

f3 = ggarrange(f3_SA_plot, f3_EC_plot, f3_KP_plot, nrow = 1, ncol = 3, common.legend = TRUE, legend = "right") %>%
  annotate_figure(bottom = text_grob(expression("Amplitude of seasonality ("*log["2"]*"(MIC) deviates)"), size = 11))

# ##################################################
# Make Figure 4
# ##################################################

#Make resistance phases table
res.phases = bind_rows(
  res.modelA.params.fil %>%
    filter(term == "phase"),
  #add second peak for org/drugs with 6 month period
  res.modelA.params.fil %>%
    filter(term == "phase") %>%
    filter(period == 6) %>%
    mutate(estimate = estimate + 6, ci.lower = ci.lower + 6, ci.upper = ci.upper + 6)
) %>%
  #filter out org/drug combinations that did not meet criteria for seasonality
  filter(!(organism == "S. aureus" & drug_code %in% c("PEN", "TET"))) %>%
  filter(!(organism == "E. coli" & drug_code %in% c("AMC", "TET"))) %>%
  filter(!(organism == "K. pneumoniae" & drug_code %in% c("AMC", "TET"))) %>%
  #manually edit phase for K.pneumo/NIT (subtract 12m) for aesthetic purposes for figure
  mutate(estimate = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", estimate - 12, estimate),
         ci.lower = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", ci.lower - 12, ci.lower),
         ci.upper = ifelse(organism == "K. pneumoniae" & drug_code == "NIT", ci.upper - 12, ci.upper)
  ) %>%
  mutate(org_drug = paste(organism, drug_code, sep = "/")) 

#Make use phases table
use.phases = bind_rows(
  use.model.params.fil %>%
    filter(term == "phase"),
  #add second peak for org/drugs with 6 month period
  use.model.params.fil %>%
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
         "E. coli/AMP" = expression(paste(italic("E. coli"), " / AMP")),
         "E. coli/CIP" = expression(paste(italic("E. coli"), " / CIP")),
         "E. coli/NIT" = expression(paste(italic("E. coli"), " / NIT")),
         "K. pneumoniae/CIP" = expression(paste(italic("K. pneumoniae"), " / CIP")),
         "K. pneumoniae/NIT" = expression(paste(italic("K. pneumoniae"), " / NIT"))
)


#Plot figure 4
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
        strip.text.x = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_blank()
  )


# ##################################################
# Make Figures 5, S4
# ##################################################

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
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 11, hjust = 0.5)
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
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 11, hjust = 0.5)
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
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 11, hjust = 0.5)
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
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 11, hjust = 0.5)
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
  
  p = ggarrange(labels, p0, p1, p2, p3, nrow = 1, ncol = 5, widths = c(1, 1.5, 1.5, 1.5, 1.5), align = "h") 
  
  return(p)
}

f5_SA = plot_corr_func(correlations, "S. aureus")
f5_EC = plot_corr_func(correlations, "E. coli") %>% annotate_figure(top = "")
f5_KP = plot_corr_func(correlations, "K. pneumoniae") %>% annotate_figure(top = "")

f5 = f5_SA
fS4 = ggarrange(f5_EC, f5_KP, ncol = 1, nrow = 2, heights = c(2, 1.3), labels = c("A.", "B."))

# ##################################################
# Make Table 1
# ##################################################

#Make amplitude, phase tables for each model
amp_model1 = res.modelA.params %>%
  filter(term == "amplitude") %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ ifelse(abs(.) < 0.01, formatC(., format = "e", digits = 1), signif(., 2))) %>%
  mutate(tmp = paste0(estimate, " (", ci.lower, ", ", ci.upper, ")")) %>%
  mutate(tmp = ifelse(term == "amplitude" & p.value.BH < 0.05, paste(tmp, "*"), tmp)) %>%
  select(organism, drug_code, period, AIC, term, tmp) %>%
  spread(term, tmp) %>%
  rename(`Model A amplitude (95% CI)` = amplitude)

amp_model2 = res.modelB.params %>%
  filter(term == "amplitude") %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ ifelse(abs(.) < 0.01, formatC(., format = "e", digits = 1), signif(., 2))) %>%
  mutate(tmp = paste0(estimate, " (", ci.lower, ", ", ci.upper, ")")) %>%
  mutate(tmp = ifelse(term == "amplitude" & p.value.BH < 0.05, paste(tmp, "*"), tmp)) %>%
  select(organism, drug_code, period, AIC, term, tmp) %>%
  spread(term, tmp) %>%
  rename(`Model B amplitude (95% CI)` = amplitude)

amp_model3 = res.modelC.params %>%
  filter(term == "amplitude") %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ ifelse(abs(.) < 0.01, formatC(., format = "e", digits = 1), signif(., 2))) %>%
  mutate(tmp = paste0(estimate, " (", ci.lower, ", ", ci.upper, ")")) %>%
  mutate(tmp = ifelse(term == "amplitude" & p.value.BH < 0.05, paste(tmp, "*"), tmp)) %>%
  select(organism, drug_code, period, AIC, term, tmp) %>%
  spread(term, tmp) %>%
  rename(`Model C amplitude (95% CI)` = amplitude)

#Make model comparison table
table1 = amp_model1 %>%
  #filter to period with lowest AIC
  group_by(organism, drug_code) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  #join with models 2 and 3
  left_join(amp_model2 %>% select(-AIC), by = c("organism", "drug_code", "period")) %>%
  left_join(amp_model3 %>% select(-AIC), by = c("organism", "drug_code", "period")) %>%
  select(organism, drug_code, period, `Model A amplitude (95% CI)`, `Model B amplitude (95% CI)`, `Model C amplitude (95% CI)`)


# ##################################################
# Make Figure S5
# ##################################################

#Set colors
demographic_colors = setNames(c("#d7b5d8", "#fbb4b9", "#f768a1", "#c51b8a", "#7a0177", "#a1d99b", "#31a354", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"),
                              c("skin/soft tissue", "abscess or fluid NOS", "blood", "respiratory tract", "urinary tract", "F", "M", "00-19", "20-39", "40-64", "65+"))

#Plot monthly num of isolates by age group
fS5_a = bind_rows(data.SA, data.EC, data.KP) %>%
  mutate(age_group = case_when(age <= 19 ~ "00-19",
                               age > 19 & age <= 39 ~ "20-39",
                               age > 39 & age <= 64 ~ "40-64",
                               age > 64 ~ "65+",
                               TRUE ~ "NA"
  )) %>%
  group_by(organism, age_group, month) %>%
  summarise(num_isolates = n_distinct(isolate_ID)) %>%
  ungroup() %>%
  
  ggplot(aes(x = month, y = num_isolates, fill = age_group)) +
  facet_wrap(~organism, scales = "free_y") +
  geom_bar(stat = "identity") +
  ggtitle("Monthly num. of isolates by age group") +
  scale_fill_manual(values = demographic_colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11)) +
  ylab("num. of isolates") +
  theme_minimal() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 10, face = "italic"),
        legend.position = "bottom") 

#Plot monthly num of isolates by sex
fS5_b = bind_rows(data.SA, data.EC, data.KP) %>%
  group_by(organism, sex, month) %>%
  summarise(num_isolates = n_distinct(isolate_ID)) %>%
  ungroup() %>%
  
  ggplot(aes(x = month, y = num_isolates, fill = sex)) +
  facet_wrap(~organism, scales = "free_y") +
  geom_bar(stat = "identity") +
  ggtitle("Monthly num. of isolates by sex") +
  scale_fill_manual(values = demographic_colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11)) +
  ylab("num. of isolates") +
  theme_minimal() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 10, face = "italic"),
        legend.position = "bottom") 

#Plot monthly num of isolates by site of infection
fS5_c = bind_rows(data.SA, data.EC, data.KP) %>%
  group_by(organism, site_of_infection, month) %>%
  summarise(num_isolates = n_distinct(isolate_ID)) %>%
  ungroup() %>%
  mutate(site_of_infection = case_when(
    site_of_infection == "skin_softtissue" ~ "skin/soft tissue",
    site_of_infection == "abscess_or_fluid_nos" ~ "abscess or fluid NOS",
    site_of_infection == "blood" ~ "blood",
    site_of_infection == "respiratory_tract" ~ "respiratory tract",
    site_of_infection == "urinary_tract" ~ "urinary tract"
  )) %>%
  
  ggplot(aes(x = month, y = num_isolates, fill = site_of_infection)) +
  facet_wrap(~organism, scales = "free_y") +
  geom_bar(stat = "identity") +
  ggtitle("Monthly num. of isolates by site of infection") +
  scale_fill_manual(values = demographic_colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11)) +
  ylab("num. of isolates") +
  theme_minimal() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        strip.text.x = element_text(size = 10, face = "italic"),
        legend.position = "bottom") 

fS5 = ggarrange(fS5_a, fS5_b, fS5_c, nrow = 3, ncol = 1, labels = c("A.", "B.", "C."), align = "hv")

# ##################################################
# Make Table S1
# ##################################################

tableS1 = res.modelA.params.under65.fil %>%
  filter(term %in% c("amplitude", "phase")) %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ ifelse(abs(.) < 0.01, formatC(., format = "e", digits = 1), signif(., 2))) %>%
  mutate(tmp = paste0(estimate, " (", ci.lower, ", ", ci.upper, ")")) %>%
  mutate(tmp = ifelse(term == "amplitude" & p.value.BH < 0.05, paste(tmp, "*"), tmp)) %>%
  select(organism, drug_code, period, term, tmp) %>%
  spread(term, tmp) %>%
  rename(`Amplitude (95% CI)` = amplitude, `Phase (95% CI)` = phase)

# ##################################################
# Make Table S2
# ##################################################

#Filter model C params to same model periods as in model A
res.modelC.params.fil = res.modelA.params.fil %>%
  count(organism, drug_code, drug_name, drug_class, period) %>%
  select(-n) %>%
  left_join(res.modelC.params, by = c("organism", "drug_code", "drug_name", "drug_class", "period"))

#Make table of beta demographics for model C
tableS2 = res.modelC.params.fil %>%
  filter(term %in% c("beta_age", "beta_sex", "beta_BL", "beta_RT", "beta_SST", "beta_AB")) %>%
  mutate_at(vars(estimate, ci.lower, ci.upper), ~ ifelse(abs(.) < 0.01, formatC(., format = "e", digits = 1), signif(., 2))) %>%
  mutate(tmp = paste0(estimate, " (", ci.lower, ", ", ci.upper, ")")) %>%
  mutate(tmp = ifelse(p.value.BH < 0.05, paste(tmp, "*"), tmp)) %>%
  select(organism, drug_code, period, term, tmp) %>%
  spread(term, tmp) %>%
  select(organism, drug_code, period, beta_age, beta_sex, beta_BL, beta_RT, beta_SST, beta_AB)


# ##################################################
# Make Table S4
# ##################################################

#Make organisms combined dataset, add age group column
data.res = bind_rows(data.SA, data.EC, data.KP) %>%
  mutate(age_group = case_when(age <= 19 ~ "00-19",
                               age > 19 & age <= 39 ~ "20-39",
                               age > 39 & age <= 64 ~ "40-64",
                               age > 64 ~ "65+",
                               TRUE ~ "NA"
  ))

#Get total number of isolates per organism
totals_byOrg = data.res %>%
  count(organism, isolate_ID) %>%
  select(-n) %>%
  count(organism) %>%
  rename(total = n)

#Make Table S4: Total number of isolates by demographics
tableS4 = bind_rows(
  totals_byOrg %>%
    mutate(total = as.character(total)) %>%
    mutate(type = "Total") %>%
    spread(organism, total),
  
  lapply(c("hospital", "patient_type", "site_of_infection", "age_group", "sex"), function(t) {
    data.res %>%
      count(organism, isolate_ID, !!sym(t)) %>%
      select(-n) %>%
      count(organism, !!sym(t)) %>%
      left_join(totals_byOrg, by = c("organism")) %>%
      mutate(p = round(n/total*100, 1)) %>%
      mutate(x = paste0(as.character(n), " (", as.character(p), "%)")) %>%
      select(organism, type = !!sym(t), x) %>%
      spread(organism, x)
  }) %>%
    bind_rows()
)

# ##################################################
# Make Table S5
# ##################################################

#Make Table S5: Percent resistance by hospital
tableS5 = data.res %>%
  count(organism, hospital, drug_name, phenotype) %>%
  left_join(data.res %>%
              count(organism, hospital, drug_name) %>%
              rename(total = n),
            by = c("organism", "hospital", "drug_name")) %>%
  filter(phenotype == "NS") %>%
  mutate(p = round(n/total*100, 1)) %>%
  mutate(percent_resistance = paste0(as.character(p), "%")) %>%
  select(organism, hospital, drug = drug_name, percent_resistance) %>%
  spread(hospital, percent_resistance)

# ##################################################
# Make Tables S6, S7
# ##################################################

#Make use AIC table comparing 6 and 12-month period models
tableS6 = use.model.params %>%
  filter(term == "amplitude") %>%
  select(drug_class, period, AIC) %>%
  group_by(drug_class) %>%
  mutate(min = min(AIC)) %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(drug_class, period, AIC_2) %>%
  spread(period, AIC_2) %>%
  rename(`12-month period` = `12`, `6-month period` = `6`)

#Make resistance AIC table comparing 6 and 12-month period models
tableS7 = res.modelA.params %>%
  filter(term == "amplitude") %>%
  select(organism, drug_name, period, AIC) %>%
  group_by(organism, drug_name) %>%
  mutate(min = min(AIC)) %>%
  ungroup() %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(organism, drug_name, period, AIC_2) %>%
  spread(period, AIC_2)

# ##################################################
# Save figures and tables
# ##################################################

#Save figures as tiff
ggsave(f1, filename = "figures/Fig1.tiff", width = 7.5, height = 5, dpi = 300, units = "in")
ggsave(f2, filename = "figures/Fig2.tiff", width = 6.5, height = 5.5, dpi = 300, units = "in")
ggsave(f3, filename = "figures/Fig3.tiff", width = 6.5, height = 4, dpi = 300, units = "in")
ggsave(f4, filename = "figures/Fig4.tiff", width = 6, height = 5.5, dpi = 300, units = "in")
ggsave(f5, filename = "figures/Fig5.tiff", width = 7.5, height = 5, dpi = 300, units = "in")

ggsave(fS1, filename = "figures/S1_Fig.tiff", width = 6.5, height = 5.5, dpi = 300, units = "in")
ggsave(fS2, filename = "figures/S2_Fig.tiff", width = 5, height = 5.5, dpi = 300, units = "in")
ggsave(fS3, filename = "figures/S3_Fig.tiff", width = 3, height = 3, dpi = 300, units = "in")
ggsave(fS4, filename = "figures/S4_Fig.tiff", width = 7.5, height = 8, dpi = 300, units = "in")
ggsave(fS5, filename = "figures/S5_Fig.tiff", width = 6.5, height = 7, dpi = 300, units = "in")

#Save tables
write_csv(table1, "figures/Table1.csv")
write_csv(tableS1, "figures/S1_Table.csv")
write_csv(tableS2, "figures/S2_Table.csv")
write_csv(tableS4, "figures/S4_Table.csv")
write_csv(tableS5, "figures/S5_Table.csv")
write_csv(tableS6, "figures/S6_Table.csv")
write_csv(tableS7, "figures/S7_Table.csv")

#Save figures as pdf
ggsave(f1, filename = "figures/Fig1.pdf", width = 7.5, height = 5)
ggsave(f2, filename = "figures/Fig2.pdf", width = 6.5, height = 5.5)
ggsave(f3, filename = "figures/Fig3.pdf", width = 6.5, height = 4)
ggsave(f4, filename = "figures/Fig4.pdf", width = 6, height = 5.5)
ggsave(f5, filename = "figures/Fig5.pdf", width = 7.5, height = 5)

ggsave(fS1, filename = "figures/S1_Fig.pdf", width = 6.5, height = 5.5)
ggsave(fS2, filename = "figures/S2_Fig.pdf", width = 5, height = 5.5)
ggsave(fS3, filename = "figures/S3_Fig.pdf", width = 3, height = 3)
ggsave(fS4, filename = "figures/S4_Fig.pdf", width = 7.5, height = 8)
ggsave(fS5, filename = "figures/S5_Fig.pdf", width = 6.5, height = 7)
