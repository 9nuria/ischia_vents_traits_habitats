# Functional changes across marine habitats due to ocean acidification

![GitHub](https://img.shields.io/badge/GitHub-39457E?style=for-the-badge&logo=github&logoColor=white)
![Gitlab](https://img.shields.io/badge/GitLab-FFA500?style=for-the-badge&logo=gitlab&logoColor=white)
![RStudio](https://img.shields.io/badge/RStudio-75AADB?style=for-the-badge&logo=RStudio&logoColor=white)

Scripts and data for Functional changes across marine habitats due to ocean acidification   
This Github repository is structured as follows:

- :file_folder: [``R``](https://github.com/9nuria/ischia_vents_habitats/tree/main/R) contains the current script and R project related to this study
- :file_folder: [``data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data) contains the following sub-folders
  - :file_folder: [``1 – raw data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/1%20–%20raw%20data) which hosts original .csv files
  - :file_folder: [``2 – elaborated data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/2%20–%20data%20generated) which hosts the data built with this script
  - :file_folder: [``3 – model``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/3%20–%20model) which hosts the models built with this script
- :file_folder: [``outputs``](https://github.com/9nuria/ischia_vents_habitats/tree/main/outputs) contains the following sub-folders
  - :file_folder: [``Figures``](https://github.com/9nuria/ischia_vents_habitats/tree/main/outputs/Figures) which hosts the figures built with this script
  - :file_folder: [``Tables``](https://github.com/9nuria/ischia_vents_habitats/tree/main/outputs/Tables) which hosts the tables from this script

This analyze has been launched with the following machine parameters

```{Session Info, echo = T}
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.6-4      lattice_0.20-45  permute_0.9-7    forcats_0.5.2    stringr_1.4.1   
 [6] dplyr_1.0.10     purrr_0.3.5      readr_2.1.3      tidyr_1.2.1      tibble_3.1.8    
[11] ggplot2_3.3.6    tidyverse_1.3.2  reshape2_1.4.4   patchwork_1.1.2  mFD_1.0.1       
[16] MASS_7.3-58.1    lme4_1.1-30      Matrix_1.5-1     geometry_0.4.6.1 devtools_2.4.5  
[21] usethis_2.1.6    brms_2.18.0      Rcpp_1.0.9       ape_5.6-2   
```