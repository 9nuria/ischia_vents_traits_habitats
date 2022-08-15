# Functional biodiversity shifts within and across benthic habitats under ocean acidification

![GitHub](https://img.shields.io/badge/github-%23181717.svg?&style=for-the-badge&logo=github&logoColor=white")
![Gitlab](https://img.shields.io/badge/gitlab-%23FCA121.svg?&style=for-the-badge&logo=gitlab&logoColor=black")
![R](https://img.shields.io/badge/R-v4.0.3-276DC3?style=for-the-badge&logo=r&logoColor=white")

Scripts and data for functional diversity under ocean acidification across different habitat types    
This Github repository is structured as follows:

- :file_folder: [``R``](https://github.com/9nuria/ischia_vents_habitats/tree/main/R) contains the current script and R project related to this study
- :file_folder: [``data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data) contains the following sub-folders
  - :file_folder: [``1 – raw data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/1%20–%20raw%20data) which hosts original .csv files
  - :file_folder: [``2 – elaborated data``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/2%20–%20data%20generated) which hosts the data built with this script
  - :file_folder: [``3 – model``](https://github.com/9nuria/ischia_vents_habitats/tree/main/data/3%20–%20model) which hosts the models built with this script
- :file_folder: [``outputs``](https://github.com/9nuria/ischia_vents_habitats/tree/main/outputs) contains the following sub-folders
  - :file_folder: [``plot``](https://github.com/9nuria/ischia_vents_habitats/tree/main/outputs/plot) which hosts the plots built with this script

This analyze has been launched with the following machine parameters

```{Session Info, echo = T}
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  12.2.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vegan_2.5-7     lattice_0.20-45 permute_0.9-7   forcats_0.5.1   stringr_1.4.0  
 [6] dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.7   
[11] ggplot2_3.3.6   tidyverse_1.3.1 reshape2_1.4.4  patchwork_1.1.1 mFD_1.0.1      
[16] MASS_7.3-55     lme4_1.1-28     Matrix_1.4-1    geometry_0.4.5  devtools_2.4.3 
[21] usethis_2.1.5   brms_2.16.3     Rcpp_1.0.8.3    ape_5.6-2  
```