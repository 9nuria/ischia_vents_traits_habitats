# ischia_vents_habitats
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

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lme4_1.1-28     brms_2.16.3     Rcpp_1.0.8.3    Matrix_1.4-1    devtools_2.4.3  usethis_2.1.5   MASS_7.3-55    
 [8] vegan_2.5-7     lattice_0.20-45 permute_0.9-7   patchwork_1.1.1 ape_5.6-2       geometry_0.4.5  reshape2_1.4.4 
[15] mFD_1.0.1       forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
[22] tibble_3.1.7    ggplot2_3.3.6   tidyverse_1.3.1

loaded via a namespace (and not attached):
  [1] readxl_1.3.1         backports_1.4.1      plyr_1.8.6           igraph_1.2.11        splines_4.0.3       
  [6] crosstalk_1.2.0      rstantools_2.1.1     inline_0.3.19        digest_0.6.29        foreach_1.5.2       
 [11] htmltools_0.5.2      viridis_0.6.2        fansi_1.0.3          magrittr_2.0.3       checkmate_2.0.0     
 [16] memoise_2.0.1        cluster_2.1.2        tzdb_0.2.0           remotes_2.4.2        modelr_0.1.8        
 [21] RcppParallel_5.1.5   matrixStats_0.61.0   xts_0.12.1           prettyunits_1.1.1    colorspace_2.0-3    
 [26] rvest_1.0.2          xfun_0.30            haven_2.4.3          callr_3.7.0          crayon_1.5.1        
 [31] jsonlite_1.8.0       iterators_1.0.14     zoo_1.8-9            glue_1.6.2           gtable_0.3.0        
 [36] emmeans_1.7.2        distributional_0.3.0 pkgbuild_1.3.1       rstan_2.21.3         abind_1.4-5         
 [41] scales_1.2.0         mvtnorm_1.1-3        DBI_1.1.2            miniUI_0.1.1.1       viridisLite_0.4.0   
 [46] xtable_1.8-4         magic_1.6-0          stats4_4.0.3         StanHeaders_2.21.0-7 DT_0.21             
 [51] htmlwidgets_1.5.4    httr_1.4.2           threejs_0.3.3        posterior_1.2.1      ellipsis_0.3.2      
 [56] pkgconfig_2.0.3      loo_2.5.0            farver_2.1.0         dbplyr_2.1.1         utf8_1.2.2          
 [61] labeling_0.4.2       tidyselect_1.1.2     rlang_1.0.2          later_1.3.0          munsell_0.5.0       
 [66] cellranger_1.1.0     tools_4.0.3          cachem_1.0.6         cli_3.3.0            generics_0.1.2      
 [71] ade4_1.7-18          broom_0.7.12         ggridges_0.5.3       FD_1.0-12            fastmap_1.1.0       
 [76] knitr_1.39           processx_3.5.2       fs_1.5.2             dendextend_1.15.2    nlme_3.1-155        
 [81] mime_0.12            projpred_2.0.2       GA_3.2.2             xml2_1.3.3           brio_1.1.3          
 [86] compiler_4.0.3       bayesplot_1.9.0      shinythemes_1.2.0    rstudioapi_0.13      gamm4_0.2-6         
 [91] testthat_3.1.2       reprex_2.0.1         stringi_1.7.6        ps_1.6.0             desc_1.4.1          
 [96] Brobdingnag_1.2-7    nloptr_2.0.0         markdown_1.1         shinyjs_2.1.0        tensorA_0.36.2      
[101] vctrs_0.4.1          pillar_1.7.0         lifecycle_1.0.1      bridgesampling_1.1-2 estimability_1.3    
[106] httpuv_1.6.5         R6_2.5.1             promises_1.2.0.1     gridExtra_2.3        gawdis_0.1.3        
[111] sessioninfo_1.2.2    codetools_0.2-18     boot_1.3-28          colourpicker_1.1.1   gtools_3.9.2        
[116] assertthat_0.2.1     pkgload_1.2.4        rprojroot_2.0.2      withr_2.5.0          shinystan_2.6.0     
[121] mgcv_1.8-39          parallel_4.0.3       hms_1.1.1            grid_4.0.3           minqa_1.2.4         
[126] coda_0.19-4          cmdstanr_0.3.0       shiny_1.7.1          lubridate_1.8.0      base64enc_0.1-3     
[131] dygraphs_1.1.1.6    
```