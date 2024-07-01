# Salmonella: Conditional Incicence (CI) UK - NL
Ojective of the analysis: The code does look at how the risk of salmonellosis in humans depends on weather variables from historical data from England and Wales. It applies the same methodology in the Netherlands for comparison.

## Contents:
1. Conditional_probability_analysis_and_plots_UK.R: conditional probability calculated from UK data and applied to UK data
2. Conditional_probability_analysis_and_plots_UKtoNL.R: conditional probability calculated from UK data and applied to NL data
3. Conditional_probability_analysis_and_plots_NLtoNL.R: conditional probability calculated from NL data and applied to NL data

Each code includes the CI as a stratification of disease and categorized weather data with a quantile distribution (i.e. each weather bin contains the same number of observations), the subsequent application of the conditional probability calculated in England and Wales, a validation by comparing empirical/historical cases with modelled cases (i.e., reconstruction), and the visualization of the reconstruction and the conditional incidence of 3 simultaneous weather factors.

## Input files and adaptations required to use the code:
Diseases incidence: Simulated_Salmonella_environment_"width_char".csv. Data assumes a time lag applied of "width" duration in days. 
Weather data:  Simulated_Laboratory_"width_char".csv. Data assumes a time lag applied of "width" duration in days.
Resident information: Sum_ByLab_1987_2016.csv. Yearly number of residents at postcode resolution.

### Settings (lines of the code No.1, as an example)
- Select time lag (line 49, variable "width"). Here 7 days by default
- Select years of interest: "first_year" (line 54) and "last_year" (line 55). Here 2000 to 2016 by default for England and Wales and 2015 to 2019 in the Netherlands.
- Select weather variables of interest ("variable_1", line 394; "variable_2", line 395; and "variable_3", line 396)
- Set the size of weather bins (variable "delta_var3", lines 413 to 442; "delta_var1", lines 452 to 476; and "delta_var2", lines 485 to 525)

Built in R Version 4.2.2
