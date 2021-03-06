# Residency, recruitment, and staging duration of hatch-year roseate terns (*Sterna dougallii*) during the pre-migratory staging period 

### Davis, K.L., S.M. Karpanty, J.A. Spendelow, J.B. Cohen, M.A. Althouse, K.C. Parsons, C.F. Luttazi, Daniel H. Catlin, Daniel Gibson

### Code/Data DOI: <a href="https://zenodo.org/badge/latestdoi/151756110"><img src="https://zenodo.org/badge/151756110.svg" alt="DOI"></a>

### Please contact first author for questions about data or code: Kayla Davis (davisk93@msu.edu)
_________________________________________________________________________________________________________________________________________
## Abstract
Avian migratory stopover and staging sites represent important energetic bottlenecks and may influence population dynamics as much as breeding or wintering periods. Roseate terns (Sterna dougallii) are an ideal species to examine staging demography because >70% of the entire endangered northwest Atlantic population stages at accessible locations around Cape Cod, MA before southward migration. We quantified hatch-year tern weekly residency, weekly recruitment rate into the staging population, and derived weekly staging population growth rate during two post-breeding, pre-migratory staging seasons (2014 and 2015) at Cape Cod National Seashore. We also estimated hatch-year tern staging duration at Cape Cod staging grounds. Tern residency probability at Cape Cod National Seashore during 2014 and 2015 was nearly 1 during the first weeks of the season and decreased steadily over the last four weeks to ~0.5 in the final week of the study. Recruitment rates into the staging population, representing the weekly per capita increase in hatch-year terns during the staging season, indicated that most terns arrived on the staging grounds during the first 3–5 weeks of the staging season (16 July–19 August). We also identified differences in staging duration between birds from the two breeding regions. Hatch-year terns from the southernmost region spent less time staging at Cape Cod National Seashore than their northern counterparts in both 2014 and 2015. These differences may indicate alternative staging strategies for individuals originating in different regions and possibly reveal differences in conditions between these areas; for example, in the availability of ephemeral prey fish. 
## Data
[HYResights_AllProofed](https://github.com/davisk93/ROST-Staging-Demography/blob/master/HYResights_AllProofed.csv): csv file containing all HY ROST resighting data

[ch_ROST](https://github.com/davisk93/ROST-Staging-Demography/blob/master/ch_ROST.csv): csv file of HY ROST capture histories

[SODA files](https://github.com/davisk93/ROST-Staging-Demography/tree/master/SODA%20files): This subdirectory contains files used for the program SODA stopover duration analysis. See subdirectory README.md for descriptions of the files included for this analysis. 
## Code
[DemographyAnalysis_ROST](https://github.com/davisk93/ROST-Staging-Demography/blob/master/DemographyAnalysis_ROST.R): R script for running RMARK analyses of HY ROST residency, recruitment, and stopover population growth during 2014-2015 staging seasons at Cape Cod National Seashore, MA. The script also contains code to create figures of parameters from the RMARK analysis and create input files for the SODA analysis.

[SODA_figures_ROST](https://github.com/davisk93/ROST-Staging-Demography/blob/master/SODA_figures_ROST.r): R script for creating stopover duration figures using output data from program SODA analysis.
