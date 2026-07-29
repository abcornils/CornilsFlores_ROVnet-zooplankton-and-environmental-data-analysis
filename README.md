# CornilsFlores_ROVnet-zooplankton-and-environmental-data-analysis

# PS122 ROV data analysis — code supplement

This repository contains the code and data underlying the figures and
tables of the manuscript "Under-ice habitat enables overwintering survival 
for zooplankton in the Central Arctic Ocean" by 
Astrid Cornils, Hauke Flores, Kim Julie Wermeyer, Philipp Anhaus, 
Carin J. Ashjian, Robert G. Campbell, Giulia Castellani, Celia Gelfman, 
Nicole Hildebrandt, Clara J. M. Hoppe, Christian Katlein, Nadine Knüppel, 
Serdar Sakinan, Fokje L. Schaafsma, Katrin Schmidt, Katyanne M. Shoemaker, 
Martina Vortkamp, Barbara Niehoff 
published in Nature Ecology & Evolution (2026). 

## `CornilsFlores_ROVanalysis_essential_data.RData`

Consolidated source data file required to run `02_run_analysis.R`. Contains the
following data frames:

| Object(s) | Content | Source |
| `meta_schulz`, `meta_schulz1` | Daily averages of SST, SSS, MLD, bottom depth along the MOSAiC drift | Schulz K et al. (2023) The Eurasian Arctic Ocean along the MOSAiC drift (2019–2020): Core hydrographic parameters. Arctic Data Center. doi:10.18739/a21j9790b |
| `combined_data_0m`, `combined_data_10m` | PAR / under-ice light | Anhaus P et al. (2022) Spectral solar radiation over and under sea ice during the MOSAiC expedition 2019/20. PANGAEA. doi:10.1594/pangaea.979856 |
| `summary_ice_0m`, `summary_ice_10m` | Sea-ice draft | Anhaus P et al. (2022) Single-beam sea-ice draft from remotely operated vehicle (ROV) surveys during the MOSAiC expedition 2019/20. PANGAEA. doi:10.1594/pangaea.952801 |
| `chla_integrated` | Depth-integrated chlorophyll a, water column (0–20 m) | Hoppe CJM et al. (2023) Water column Chlorophyll a concentrations during the MOSAiC expedition (PS122) in the Central Arctic Ocean 2019–2020. PANGAEA. doi:10.1594/pangaea.963277 |
| `rov_chla_ice_sel` | Depth-integrated chlorophyll a, bottom sea ice | Hoppe CJM et al. (2026) Chlorophyll a concentrations in first year sea ice during the MOSAiC expedition (PS122) in the Central Arctic Ocean 2019–2020. PANGAEA. doi:10.1594/pangaea.990053 |
| `poc_integrated` | Depth-integrated particulate organic carbon (POC), water column | Hoppe CJM et al. (2025) Concentrations and stable isotopic composition of particulate organic carbon and nitrogen from water column samples collected in the central Arctic Ocean in 2019/2020 during the MOSAiC campaign. PANGAEA. doi:10.1594/pangaea.980518 |
| `merged_poc`, `merged_final` | ROV net station metadata, merged with all of the environmental parameters above (± 7 days) | ROV station metadata: see. Cornils et al. (2024) |
| `ps122_ROV_abundance`, `ps122_ROV_biovolume`, `ps122_ROV_biomass_final`, `ps122_ROV_carbon_final`, `ROV_biovolume_final`, `ps122_ROV_all_final` | Cornils A et al. (2024) Zooplankton abundance, biovolume, biomass and respiration from ROVnet samples during the MOSAiC expedition (PS122) in the Central Arctic Ocean. PANGAEA. doi:10.1594/PANGAEA.968897

## `CornilsFlores_ROVanalysis.R`

Loads the RData file above and reproduces every step of the analysis:
Table 1, Figures 1–5, and Supplementary Tables and Figures.

To run: place both files in the same folder, adjust the `setwd()` path at
the top of the script, and run top to bottom in a fresh R session. A fixed
random seed (`set.seed(42)`) is set at the start so that all stochastic
steps (NMDS, permutation tests, bootstrapping) give identical results on
every run. Final multi-panel figures were assembled from the individual
panels manually in Inkscape.

## Use of AI assistance

The analysis code was written by Astrid Cornils. An AI assistant (Claude)
was used to review, debug, and refactor the final version of the analysis
script prior to publication — specifically: separating data preparation
from the analysis script, fixing a package-masking bug (`plyr` vs.
`dplyr`), adding random seeds for reproducibility, and adding section
headers linking code blocks to the corresponding manuscript Tables/Figures.
All underlying data processing and statistical analysis choices are the
authors' own.

## Contact

Astrid Cornils, 
Alfred-Wegener-Institut, 
Helmholtz-Zentrum für Polar- und Meeresforschung
astrid.cornils@awi.de

