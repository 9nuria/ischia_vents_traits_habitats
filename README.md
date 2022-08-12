# Functional biodiversity shifts within and across benthic habitats under ocean acidification
----------------------

![GitHub](https://img.shields.io/badge/github-%23181717.svg?&style=for-the-badge&logo=github&logoColor=white")
![Gitlab](https://img.shields.io/badge/gitlab-%23FCA121.svg?&style=for-the-badge&logo=gitlab&logoColor=black")
![R](https://img.shields.io/badge/R-v4.0.3-276DC3?style=for-the-badge&logo=r&logoColor=white")
[![GitHub commits](https://badgen.net/github/commits/9nuria/ischia_vents_habitats.js)](https://github.com/9nuria/ischia_vents_habitats.js/commit/)
[![GitHub latest commit](https://badgen.net/github/last-commit/9nuria/ischia_vents_habitats.js)](https://github.com/9nuria/ischia_vents_habitats.js/commit/)

Scripts and data for functional diversity under ocean acidification across different habitat types

Data architecture> Scripts and data are from the folder named "seb" within the master branch valeriano/nuria_acid. Dowloaded from github on 2021-12-29. 

Folders
folder data> data for the analysis, including raw data 
folder FD> functional diversity files
folder plot> plots obtained

Scripts
A dataset for FD: preparing datasets for analysis using raw data
B FD: computing FD with the new package
C mean beta: Correlation of beta-diversity between taxonomic and functional facets as boxplot and scatter plot. Figures for SM 
D trait dominance: category abundance for each trait across habitats and pH conditions. Remove it? 
E mds permdisp: taxonomic and functional MDS based on hill numbers with stats
F functional hull: new figure FD with Ambient and Low pH conditions
G diversity indices: diversity indices across habitats and 2 pH conditions witth stats
H change trait abund: Change in trait abundance with pH across habitats

```{Session Info, echo = T}
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  12.2.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lme4_1.1-28     brms_2.16.3     Rcpp_1.0.8.3    Matrix_1.4-1    devtools_2.4.3  usethis_2.1.5   MASS_7.3-55    
 [8] vegan_2.5-7     lattice_0.20-45 permute_0.9-7   patchwork_1.1.1 ape_5.6-2       geometry_0.4.5  reshape2_1.4.4 
[15] mFD_1.0.1       forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
[22] tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1
```