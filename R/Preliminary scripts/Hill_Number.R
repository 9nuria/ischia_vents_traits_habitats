##### RAW VALUES ----
# Hill Number with Functionality
# SRA_tax = 0.057 ± 0.007     # SRL_tax = 0.039 ± 0.001
# CA_ tax = 0.076 ± 0.019     # CL_ tax = 0.119 ± 0.030
# RA_ tax = 0.036 ± 0.009     # RL_ tax = 0.019 ± 0.000
# CoA_tax = 0.052 ± 0.008     # CoL_tax = 0.047 ± 0.010

# Hill Number with Taxonomy
# SRA_tax = 0.397 ± 0.025     # SRL_tax = 0.206 ± 0.000
# CA_ tax = 0.579 ± 0.045     # CL_ tax = 0.623 ± 0.029
# RA_ tax = 0.316 ± 0.023     # RL_ tax = 0.199 ± 0.000
# CoA_tax = 0.409 ± 0.029     # CoL_tax = 0.218 ± 0.017

##### PERMANOVA W/ HILL NIMBERS ----
#### TAXONOMY 
### Deep Reefs
#           Df  Sum Sq  Mean Sq       F N.Perm Pr(>F)    
# Groups     1 0.14564 0.145644 14.3830    999  0.001 **
### Reefs
#           Df  Sum Sq  Mean Sq       F N.Perm Pr(>F)    
# Groups     1 0.054124 0.054124 8.2957    999  0.009 **
### Cave
#           Df  Sum Sq  Mean Sq       F N.Perm Pr(>F)    
# Groups     1 0.00249 0.0024937 0.1139    999   0.72
### Shallow Reefs
#           Df  Sum Sq  Mean Sq       F N.Perm Pr(>F)    
# Groups     1 0.048716 0.048716 9.6049    999  0.011 *

#### FUNCTIONALITY 
### Deep Reefs
#           Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0002816 0.00028157 0.2010    999  0.683
### Reefs
#           Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0027230 0.00272302 2.7305    999  0.094 .
### Cave
#           Df    Sum Sq   Mean Sq       F N.Perm Pr(>F) 
# Groups     1 0.0042020 0.0042024 0.42060    999  0.535
### Shallow Reefs
#           Df    Sum Sq   Mean Sq       F N.Perm Pr(>F) 
# Groups     1 0.0010199 0.0010199 2.54590    999  0.137

##### Code ----
mean_f1 = list()
for (Q in 1:n) { 
  beta_df[[Q]]        <- mFD::dist.to.df(quadrats_beta_hill[[Q]])
  beta_df.taxo[[Q]]   <- beta_df[[Q]] %>% dplyr::select(-funct_q1)
  habitat             <- c("CoL")
  quad                <- info[info$condition %in% habitat,]$Quadrats 
  divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,]         # select quadrat–habitat (i.e. 1s1)
  divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
  mean_f1[[Q]]        <- mean(divtaxo[[Q]]$taxo_q1)                                 # select quadrat SR for asb.2
}
round(mean(abind::abind(mean_f1)),3) ; round(sd(abind::abind(mean_f1)),3)
for (Q in 1:n) { 
  beta_df[[Q]]        <- mFD::dist.to.df(quadrats_beta_hill[[Q]])
  beta_df.taxo[[Q]]   <- beta_df[[Q]] %>% dplyr::select(-taxo_q1)
  habitat             <- c("CoL", "CoA")
  quad                <- info[info$condition %in% habitat,]$Quadrats 
  divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,]         # select quadrat–habitat (i.e. 1s1)
  divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
  beta_long_fun[[Q]]  <- long_to_wide_distance(divtaxo[[Q]])
  groupf[[Q]]         <- sapply(labels(beta_long_fun[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
  modf[[Q]]           <- betadisper(beta_long_fun[[Q]], groupf[[Q]], bias.adjust = TRUE)
  anova(modf[[Q]]) ; permutest(modf[[Q]], permutations = 999, pairwise = T) # ; plot(modf[[Q]]) ; boxplot(modf[[Q]])
}