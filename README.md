# Code associated with "Large variation in the association between seasonal antibiotic use and resistance across multiple bacterial species and antibiotic classes"

### Daphne S. Sun, Stephen M. Kissler, Sanjat Kanjilal, Scott W. Olesen, Marc Lipsitch, Yonatan H. Grad

Due to data sharing restrictions, we are unable to publish the raw antibiotic use and resistance datasets used in this analysis. Because these raw datasets are used as inputs in `regression_analysis.R`, this script cannot be run. However, we have included the outputs of this code, including the sinusoidal model values and seasonal deviates tables, which can be used to run the `correlation_analysis.R` and `make_figures.R` scripts. These code can reproduce all analyses and figures associated with this publication except for Supplementary Tables 1, 2, and 3, which are summary tables of the raw use and resistance data.

# Files 
- `README.md`: this file

- `regression_analysis.R`: code to run non-linear regressions to fit seasonal use and resistance data to sinusodal models
- `correlation_analysis.R`: code to calculate correlations between use and resistance seasonal deviates 
- `make_figures.R`: code to make all figures included in the publication from outputs of `regression_analysis.R` and `correlation_analysis.R`

- `data/`
  - `antibiotic_use_data.csv`: table of total antibiotic claims per 1000 ppl by antibiotic class, year, and month
  
  - `model_values_use.csv`: table of sinusoidal model-fitted parameters for seasonal antibiotic use of each antibiotic class (Output by `regression_analysis.R`)
  - `model_values_resistance.csv`: table of sinusoidal model-fitted parameters for seasonal antibiotic resistance in each species-antibiotic combination (Output by `regression_analysis.R`)
  
  - `model_amplitude_pvalues_use.csv`: table of model-fitted amplitudes for seasonal antibiotic use (same as in `model_values_use.csv`), with p-values before and after multiple testing correction (Output by `regression_analysis.R`)
  - `model_amplitude_pvalues_resistance.csv`: table of model-fitted amplitudes for seasonal antibiotic resistance (same as in `model_values_resistance.csv`), with p-values before and after multiple testing correction (Output by `regression_analysis.R`)
  
  - `seasonal_deviates_use.csv`: table of seasonal deviates of use averaged by month for each antibiotic class (Output by `regression_analysis.R`)
  - `seasonal_deviates_resistance.csv`: table of seasonal deviates in resistance averaged by month for each species-antibiotic combination (Output by `regression_analysis.R`)
  
  - `Ecoli_AMC_AMP_12m_model_values.csv`: table of model-fitted parameters for seasonal antibiotic resistance to Amoxicillin/Clavulanate and Ampicillin in *E. coli* using a 12-month period sinuosoidal model, as opposed to the 6-month period model, which is reported in the `model_values_resistance.csv` file. These data are used to generate Supplementary Figure 3. (Output by `regression_analysis.R`)
  - `Ecoli_AMC_AMP_12m_seasonal_deviates.csv`: table of seasonal deviates in resistance to Amoxicillin/Clavulanate and Ampicillin in *E. coli* calculated from the output of the model fits to 12-month period sinuosoidal model, as opposed to the 6-month period model, which is reported in the `seasonal_deviates_resistance.csv` file. These data are used to generate Supplementary Figure 3. (Output by `regression_analysis.R`)
  
  - `correlations.csv`: table of calculated Spearman correlation coefficients for seasonal use and resistance (Output by `correlation_analysis.R`)

- `figures/`
  - `Fig1.pdf`: (Output by `make_figures.R`)
  - `Fig2.pdf`: (Output by `make_figures.R`)
  - `Fig3.pdf`: (Output by `make_figures.R`)
  - `Fig4.pdf`: (Output by `make_figures.R`)
  - `Fig5.pdf`: (Output by `make_figures.R`)
  - `FigS1.pdf`: (Output by `make_figures.R`)
  - `FigS2.pdf`: (Output by `make_figures.R`)
  - `FigS3.pdf`: (Output by `make_figures.R`)
  - `Table_S1.csv`
  - `Table_S2.csv`
  - `Table_S3.csv`
  - `Table_S4.csv`: (Output by `regression_analysis.R`)
  - `Table_S5.csv`: (Output by `regression_analysis.R`)
  

Correspondence: Daphne Sun <dssun@g.harvard.edu>
