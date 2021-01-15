# Code associated with "Large variation in the association between seasonal antibiotic use and resistance across multiple bacterial species and antibiotic classes"

### Daphne S. Sun, Stephen M. Kissler, Sanjat Kanjilal, Scott W. Olesen, Marc Lipsitch, Yonatan H. Grad

Included are all code and raw data need to reproduce all analyses associated with this preprint. Antibiotic use data was obtained from the Massachusetts All Payers Claims Database (https://www.chiamass.gov/ma-apcd/). Antibiotic resistance data for 3 organisms, *Staphylococcus aureus*, *Escherichia coli*, and *Klebsiella pneumoniae*, was obtained from Brigham and Women's Hospital and Massachusetts General Hospital. 

# Files 
- `README.md`: this file

- `regression_analysis.R`: code to run non-linear regressions to fit seasonal use and resistance data to sinusodal models
- `correlation_analysis.R`: code to calculate correlations between use and resistance seasonal deviates 
- `make_figures.R`: code to make all figures included in the publication from outputs of `regression_analysis.R` and `correlation_analysis.R`

- `raw_data/`: contains raw antibiotic use and resistance data
  - `antibiotic_use_data.csv`: raw antibiotic claims data for Boston, MA, aggregated by antibiotic class, year, and month
  - `Ecoli_antibiotic_resistance_data.csv`: raw *E. coli* antibiotic resistance data, where each row represents one antibiotic test on one isolate
  - `Kpneumoniae_antibiotic_resistance_data.csv`: raw *K. pneumoniae* antibiotic resistance data, where each row represents one antibiotic test on one isolate 
  - `Saureus_antibiotic_resistance_data`: raw *S. aureus* antibiotic resitsance data, where each row represents one antibiotic test on one isolate

- `tables/`: contains intermediate files output during analysis
  - `model_values_use.csv`: table of sinusoidal model-fitted parameters for seasonal antibiotic use of each antibiotic class (Output by `regression_analysis.R`)
  - `model_values_resistance.csv`: table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic combination (Output by `regression_analysis.R`)
  
  - `model_amplitude_pvalues_use.csv`: table of model-fitted amplitudes for seasonal antibiotic use (same as in `model_values_use.csv`), with p-values before and after multiple testing correction (Output by `regression_analysis.R`)
  - `model_amplitude_pvalues_resistance.csv`: table of model-fitted amplitudes for seasonal antibiotic resistance (same as in `model_values_resistance.csv`), with p-values before and after multiple testing correction (Output by `regression_analysis.R`)
  
  - `seasonal_deviates_use.csv`: table of seasonal deviates of use averaged by month for each antibiotic class (Output by `regression_analysis.R`)
  - `seasonal_deviates_resistance.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination (Output by `regression_analysis.R`)
  
  - `Ecoli_AMC_AMP_12m_model_values.csv`: table of model-fitted parameters for seasonal antibiotic resistance to Amoxicillin/Clavulanate and Ampicillin in *E. coli* using a 12-month period sinuosoidal model, as opposed to the 6-month period model, which is reported in the `model_values_resistance.csv` file. These data are used to generate Supplementary Figure 3. (Output by `regression_analysis.R`)
  - `Ecoli_AMC_AMP_12m_seasonal_deviates.csv`: table of seasonal deviates in resistance to Amoxicillin/Clavulanate and Ampicillin in *E. coli* calculated from the output of the model fits to 12-month period sinuosoidal model, as opposed to the 6-month period model, which is reported in the `seasonal_deviates_resistance.csv` file. These data are used to generate Supplementary Figure 3. (Output by `regression_analysis.R`)
  
  - `correlations.csv`: table of calculated Spearman correlation coefficients for seasonal use and resistance (Output by `correlation_analysis.R`)

- `figures/`: contains all figures associated with this publication
  - `Fig1.pdf`: "Seasonal patterns of antibiotic use by class" (Output by `make_figures.R`)
  - `Fig2.pdf`: "Seasonality of antibiotic use and resistance by class in *Staphylococcus aureus*" (Output by `make_figures.R`)
  - `Fig3.pdf`: "Amplitudes of seasonality of resistance by species and antibiotic class" (Output by `make_figures.R`)
  - `Fig4.pdf`: "Phases of seasonality of use and resistance by species and antibiotic class" (Output by `make_figures.R`)
  - `Fig5.pdf`: "Seasonal resistance to multiple antibiotics is positively correlated with seasonal use of penicillins and macrolides" (Output by `make_figures.R`)
  - `FigS1.pdf`: "Seasonality of antibiotic use and resistance by class in *Escherichia coli*" (Output by `make_figures.R`)
  - `FigS2.pdf`: "Seasonality of antibiotic use and resistance by class in *Klebsiella pneumoniae*" (Output by `make_figures.R`)
  - `FigS3.pdf`: "Seasonality of use and resistance for penicillins in *Escherichia coli* with a 12-month period model" (Output by `make_figures.R`)
  - `Table_S1.csv`: "Percent of claims by individual antibiotics within each class"
  - `Table_S2.csv`: "Total number of isolates by demographics"
  - `Table_S3.csv`: "Antibiotics included in analysis and percent resistance by hospital"
  - `Table_S4.csv`: "Comparison of the Akaike information criterion (AIC) values across two sinusoidal models for antibiotic use" (Output by `regression_analysis.R`)
  - `Table_S5.csv`: "Comparison of the Akaike information criterion (AIC) values across two sinusoidal models for antibiotic resistance" (Output by `regression_analysis.R`)
  

Correspondence: Daphne Sun <dssun@g.harvard.edu>