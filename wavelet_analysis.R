#This script performs a wavelet analysis on the antibiotic use data

# ##################################################
# Inputs
# ##################################################

#Load antibiotic use dataset
data.use = read_csv("raw_data/antibiotic_use_data.csv")

#Edit input table
data.use = data.use %>%
  mutate(month = month - 1) %>%
  mutate(year_month = year + month/12) %>%
  mutate(t = (year - min(year)) * 12 + month) 

# ##################################################
# Compute wavelets
# ##################################################

#Compute use wavelet for each antibiotic class (uses WaveletComp package)
fence = function(x, lower, upper) pmax(lower, pmin(upper, x))

wt_plot_data = function(w) {
  with(w, {
    crossing(i = seq_along(axis.1), j = seq_along(axis.2)) %>%
      mutate(x = axis.1[i],
             y = axis.2[j],
             scale = Scale[j],
             amplitude = map2_dbl(j, i, ~ Ampl[.x, .y]),
             power = map2_dbl(j, i, ~ Power[.x, .y]),
             p.value = map2_dbl(j, i, ~ Power.pval[.x, .y]),
             coi_raw = approx(coi.1, coi.2, x)$y,
             coi = fence(coi_raw, min(y), max(y)),
             sig = as.double(p.value < 0.05))
  })
}

wavelet.results = data.use %>%
  nest(-drug_class) %>%
  mutate(wavelet = map(data, ~ WaveletComp::analyze.wavelet(., "mean_daily_claims_per_10000ppl", verbose = FALSE))) %>%
  mutate(wavelet_plot_data = map(wavelet, wt_plot_data)) 

#Make figure data table
fig_data = wavelet.results %>%
  select(drug_class, wavelet_plot_data) %>%
  unnest()

#Make wavelet figure (Fig S6)
wavelet_plot = fig_data %>%
  ggplot(aes(x-1, y)) +
  geom_tile(aes(fill = amplitude)) +
  geom_ribbon(aes(ymin = coi, ymax = max(y)), alpha = 0.25) +
  geom_contour(aes(z = as.double(p.value < 0.05)), color = 'black', size = 0.5, bins = 1) +
  facet_wrap(~drug_class) +
  scale_fill_distiller(palette = 'RdYlBu') +
  geom_hline(yintercept = 3.65, linetype = 2) +
  geom_hline(yintercept = 2.58, linetype = 2) +
  labs(x = "time (months since Jan 2011)", y = expression(log["2"]*"(period) (months)")) +
  theme_minimal() +
  theme(text = element_text(size = 12))

# ##################################################
# Save figure
# ##################################################

write_csv(fig_data, "figure_data/S6_Fig/S6_Fig_data.csv")
ggsave(wavelet_plot, filename = "figures/S6_Fig.tiff", width = 7.5, height = 5, dpi = 300, units = "in")
ggsave(wavelet_plot, filename = "figures/S6_Fig.pdf", width = 7.5, height = 5)
