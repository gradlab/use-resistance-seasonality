# Code associated with "Large variation in the association between seasonal antibiotic use and resistance across multiple bacterial species and antibiotic classes"

### Daphne S. Sun, Stephen M. Kissler, Sanjat Kanjilal, Scott W. Olesen, Marc Lipsitch, Yonatan H. Grad

Included are all code and raw data need to reproduce all analyses associated with this preprint. Antibiotic use data was obtained from the Massachusetts All Payers Claims Database (https://www.chiamass.gov/ma-apcd/). Antibiotic resistance data for 3 organisms, *Staphylococcus aureus*, *Escherichia coli*, and *Klebsiella pneumoniae*, was obtained from Brigham and Women's Hospital and Massachusetts General Hospital. 

# Files 
- `README.md`: this file

- `seasonality_regression_functions.R`: functions to run non-linear sinusoidal regressions for use and resistance
- `seasonality_regression_analysis_use.R`: code to run non-linear regressions to fit seasonal use data to sinusoidal models
- `seasonality_regression_analysis_resistance_modelA.R`: code to run non-linear regressions to fit seasonal resistance data to sinusoidal models, unadjusted for demographics
- `seasonality_regression_analysis_resistance_modelB.R`: code to run non-linear regressions to fit seasonal resistance data to sinusoidal models, adjusted for patient age and sex
- `seasonality_regression_analysis_resistance_modelC.R`: code to run non-linear regressions to fit seasonal resistance data to sinusoidal models, adjusted for patient age, sex, and site of infection
- `correlation_analysis.R`: code to calculate correlations between use and resistance seasonal deviates 
- `wavelet_analysis.R`: code to perform wavelet analysis for antibiotic use data
- `make_figures.R`: code to make all figures included in the publication 

- `raw_data/`: contains raw antibiotic use and resistance data
  - `antibiotic_use_data.csv`: raw antibiotic claims data for Boston, MA, aggregated by antibiotic class, year, and month
  - `Ecoli_antibiotic_resistance_data.csv`: raw *E. coli* antibiotic resistance data, where each row represents one antibiotic test on one isolate
  - `Kpneumoniae_antibiotic_resistance_data.csv`: raw *K. pneumoniae* antibiotic resistance data, where each row represents one antibiotic test on one isolate 
  - `Saureus_antibiotic_resistance_data.csv`: raw *S. aureus* antibiotic resitsance data, where each row represents one antibiotic test on one isolate

- `tables/`: contains intermediate files output during analysis
  - `use_model_values.csv`: full table of sinusoidal model-fitted parameters for seasonal antibiotic use of each antibiotic class (Output by `seasonality_regression_analysis_use.R`)
  - `use_seasonal_deviates.csv`: table of seasonal deviates of use averaged by month for each antibiotic class (Output by `seasonality_regression_analysis_use.R`)
  
  - `resistance_modelA_values.csv`: full table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic, using a model unadjusted for demographics (Output by `seasonality_regression_analysis_resistance_modelA.R`)
  - `resistance_modelA_seasonal_deviates.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination, using a model unadjusted for demographics (Output by `seasonality_regression_analysis_resistance_modelA.R`)
  
  - `resistance_OP_under65_modelA_values.csv`: full table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic, using a model unadjusted for demographics, after restricting to isolates from outpatients under 65 years old (Output by `seasonality_regression_analysis_resistance_modelA.R`)
  - `resistance_OP_under65_modelA_seasonal_deviates.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination, using a model unadjusted for demographics, after restricting to isolates from outpatients under 65 years old (Output by `seasonality_regression_analysis_resistance_modelA.R`)
  
  - `resistance_modelB_values.csv`: full table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic, using a model adjusted for age and sex (Output by `seasonality_regression_analysis_resistance_modelB.R`)
  - `resistance_modelB_seasonal_deviates.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination, using a model adjusted for age and sex (Output by `seasonality_regression_analysis_resistance_modelB.R`)
 
  - `resistance_modelC_values.csv`: full table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic, using a model adjusted for age, sex, and site of infection (Output by `seasonality_regression_analysis_resistance_modelC.R`)
  - `resistance_modelC_seasonal_deviates.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination, using a model adjusted for age, sex, and site of infection (Output by `seasonality_regression_analysis_resistance_modelC.R`)
  
  - `correlations.csv`: table of calculated Spearman correlation coefficients for seasonal use and resistance (Output by `correlation_analysis.R`)

- `figures/`: contains all figures associated with this publication
  - `Fig1.pdf`: "Seasonal patterns of antibiotic use by class" (Output by `make_figures.R`)
  - `Fig2.pdf`: "Seasonality of antibiotic use and resistance by class in *Staphylococcus aureus*" (Output by `make_figures.R`)
  - `Fig3.pdf`: "Amplitudes of seasonality of resistance by species and antibiotic class" (Output by `make_figures.R`)
  - `Fig4.pdf`: "Phases of seasonality of use and resistance by species and antibiotic class" (Output by `make_figures.R`)
  - `Fig5.pdf`: "Seasonal resistance to multiple antibiotics is positively correlated with seasonal use of penicillins and macrolides" (Output by `make_figures.R`)
  - `Table1.csv`: "Comparison of estimated amplitudes of seasonality across three sinusoidal models for resistance" (Output by `make_figures.R`)
  
  - `S1_Fig.pdf`: "Seasonality of antibiotic use and resistance by class in *Escherichia coli*" (Output by `make_figures.R`)
  - `S2_Fig.pdf`: "Seasonality of antibiotic use and resistance by class in *Klebsiella pneumoniae*" (Output by `make_figures.R`)
  - `S3_Fig.pdf`: "Seasonality of use and resistance for penicillins in *Escherichia coli* with a 12-month period model" (Output by `make_figures.R`)
  - `S4_Fig.pdf`: "Spearman correlations between seasonal use and resistance with 0-3 months lag in *E. coli* and *K. pneumoniae*" (Output by `make_figures.R`)
  - `S5_Fig.pdf`: "Seasonal incidence of infection by demographic group or site of infection" (Output by `make_figures.R`)
  - `S6_Fig.pdf`: "Wavelet analysis of antibiotic use by class" (Output by `wavelet_analysis.R`)
  
  - `S1_Table.csv`: "Amplitudes and phases of seasonality of resistance in patients under 65 years old" (Output by `make_figures.R`)
  - `S2_Table.csv`: "Age, sex, and site of infection coefficients for adjusted sinusoidal model of seasonal resistance" (Output by `make_figures.R`)
  - `S3_Table.csv`: "Percent of claims by individual antibiotics within each class"
  - `S4_Table.csv`: "Total number of isolates by demographics" (Output by `make_figures.R`)
  - `S5_Table.csv`: "Antibiotics included in analysis and percent resistance by hospital" (Output by `make_figures.R`)
  - `S6_Table.csv`: "Comparison of the Akaike information criterion (AIC) values between 6- and 12-month period models for antibiotic use" (Output by `make_figures.R`)
  - `S7_Table.csv`: "Comparison of the AIC values between 6- and 12-month period models for antibiotic resistance" (Output by `make_figures.R`)


Correspondence: Daphne Sun <dssun@g.harvard.edu>