

library(tidyverse)

# computing average and standard deviation and sample size of FE biomass in the 2 conditions

# quadrats tidy style
tab <- quadrats_fe_cover %>%
  as.data.frame %>%
  rownames_to_column("Quadrats") %>%
  left_join( select(sites_quadrats_info, Quadrats, pH, habitat )  ) %>%
  mutate(pH_hab = paste( substr(pH,1,3), habitat ,sep="_"), .after = Quadrats ) %>%
  select(- pH, - habitat)
  
# function to compute standard error
std_err <- function(x) {
  sd(x)/sqrt(length(x))
}

# average and standard error
pH_hab_mean <- tab %>%
  select(-Quadrats) %>%
  group_by(pH_hab) %>%
  summarise( across(.cols = everything(), .fns = mean, .names=NULL ) )
  
pH_hab_se <- tab %>%
  select(-Quadrats) %>%
  group_by(pH_hab) %>%
  summarise( across(.cols = everything(), .fns = std_err, .names=NULL ) )

  head(pH_hab_se)
  

  