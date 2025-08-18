
Code for the paper "**Antibiotic resistance selection and deselection in municipal wastewater from 47 countries**"

The code for the core statistical analysis of the paper is available in *core_statistical_analysis.R*.



The raw data used for this statistical analysis will be provided in *data/CFU_data.xlsx* after the paper has been published.



Additionally, a script for calculating 1% mics for ecoli based on data from EUCAST is provided in *mic_1percent_ecoli_script.R*. 
Data used as input for the script is in *data/example_eucast_ecoli_mic_data.xlsx* and expected output in *output/mic_1percent_ecoli.xlsx*.

1. System requirements

The scripts require R, and have been tested using R version 4.4.0 on a Windows 11 x64 laptop with 32 GB RAM.
The expected run time on such a laptop is short, less than one minute.
All Required R packages can be installed from CRAN (https://cran.r-project.org/), and the versions in use are as follows:

| Package   | Version   | Purpose                                         |
|-----------|-----------|-------------------------------------------------|
| tidyverse | 2.0.0     | Data wrangling                                  |
| MASS      | 7.3.60.2  | Negative binomial modeling                      |
| emmeans   | 1.10.3    | Estimated marginal means and contrasts          |
| broom     | 1.0.6     | Tidying model outputs                           |
| magrittr  | 2.0.3     | Pipe operators (%<>%, %>%)                       |
| multcomp  | 1.4.25    | Multiple comparison adjustments                 |
| readxl    | 1.4.3     | Read from Excel files                                |
| writexl   | 1.5.0     | Write to Excel files                               |

2. Installation guide
The packages and versions in use can be installed with install.packages("<<packagename>>","version") for all packages, e.g install.packages("MASS",version="7.3.60.2"). 
Expected install time of the packages on a "normal" desktop computer is about 30 minutes, with the majority of time spent on installing tidyverse.  

3. Demo/Instructions for use

Data must be downloaded along with the scripts, and the path to the data updated as appropriate in the R script prior to running. 
Similarly, the path to where outputs are to be stored can be adapted in the end of the script, if need be.
The expected outputs with the provided input data will be provided under outputs/<<output_files>> after the paper has been published.


The functionality of the code is highly specific to the input format of the provided data, and it was not intended to be more generally applicable. 
If you would want to apply the statistical analysis framework presented on your own data, a suggestion is to look at the format of the data going into the function "counts_modeller" and adapt your input data accordingly.

