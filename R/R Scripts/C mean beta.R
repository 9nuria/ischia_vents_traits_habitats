# mean of taxo and functional beta within / between sites & conditions
# Correlation of beta-diversity values between taxonomic and functional facets as boxplot and scatter plot. Figures for SM 


## preparing #####
rm(list=ls())

# libraries ---
library(tidyverse)
library(mFD)


# folder with ready to use data ----
dir_data<-"./data"
# folder with FD values ----
dir_results<-"./FD"
# folder with plot  ----
dir_plot<-"./plot"

# loading metadata built in script A ----
load(file.path(dir_data,"sites_quadrats_info.Rdata"))

# loading beta diversity values computed in script B ----
load(file=file.path(dir_results,"quadrats_beta_hill.Rdata") )


# building a dataframe with quadrats metadata and beta-diversity values ####

# from distance matrices to dataframe
beta_df<-mFD::dist.to.df(quadrats_beta_hill)

# habitats and pH level of quadrats, recoding pH into ambient vs low
quadrats_habph<-select(sites_quadrats_info, Quadrats, pH, habitat) %>%
  mutate(ph=substr(pH,1,3)) %>%
  select(-pH)

# left_joining first time for quadrat x1 and renaming
# same with x2
quadrats_beta_df<-left_join(beta_df, quadrats_habph, 
                by = c("x1" = "Quadrats")
                ) %>%
  rename(id_x1=x1, ph_x1=ph, hab_x1=habitat) %>%
  left_join(quadrats_habph, 
            by = c("x2" = "Quadrats")
  ) %>%
  rename(id_x2=x2, ph_x2=ph, hab_x2=habitat)
  
head(quadrats_beta_df)

# adding new variable with type of pair
quadrats_beta_df <- mutate(quadrats_beta_df,  
                           pair_type= case_when(
              hab_x1==hab_x2 & ph_x1==ph_x2 ~ "within_hab*ph",
              hab_x1=="cave" & hab_x2==hab_x1 & ph_x1!=ph_x2 ~ "interpH_cave",
              hab_x1=="reef" & hab_x2==hab_x1 & ph_x1!=ph_x2 ~ "interpH_reef",
              hab_x1=="shallow_reef" & hab_x2==hab_x1 & ph_x1!=ph_x2 ~ "interpH_shallow",
              hab_x1=="deep_reef" & hab_x2==hab_x1 & ph_x1!=ph_x2 ~ "interpH_deep",
              
              ph_x1=="low" & ph_x2=="low" & hab_x1!=hab_x2 ~ "interhab_lowpH",
              ph_x1=="amb" & ph_x2=="amb" & hab_x1!=hab_x2 ~ "interhab_ambpH",
              TRUE   ~ "_interhab_interpH"
                           )
)

head(quadrats_beta_df)

summary(as.factor(quadrats_beta_df$pair_type))


# boxplot for taxo and functional beta across type of pairs (Figure for SM) #####

beta_taxo_plot<-ggplot(data = quadrats_beta_df, 
                       aes(x = pair_type, y = taxo_q1) ) +
  geom_boxplot() +
  xlab("Habitat and pH of pair of quadrats") +
  ylab("Taxo beta (q=1)")

beta_funct_plot<-ggplot(data = quadrats_beta_df, 
                       aes(x = pair_type, y = funct_q1) ) +
  geom_boxplot() +
  xlab("Habitat and pH of pair of quadrats") +
  ylab("Func beta (q=1)")

boxplot_beta<-beta_taxo_plot / beta_funct_plot
ggsave(boxplot_beta, filename="boxplot_beta_taxo_funct.png", path = dir_plot, width = 9, height=5)

# plot of taxo vs functional beta within each type of pair #####
beta_cor_plot<-ggplot(data = quadrats_beta_df, 
                        aes(x = taxo_q1, y = funct_q1) ) +
  geom_point(size=0.05) +
  facet_wrap(vars(pair_type))
ggsave(beta_cor_plot, filename="beta_taxo_funct.png", path = dir_plot, width = 4)

