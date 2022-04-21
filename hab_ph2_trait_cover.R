# Change abundance traits

## preparing #####
rm(list=ls())

# libraries ---
library(tidyverse)
library(mFD)
library(patchwork)

# folder with ready to use data ----
dir_data<-"./data"

# loading dataset built in script A ----
load(file.path(dir_data, "sites_quadrats_info.Rdata"))
load(file=file.path(dir_data,"quadrats_species_cover.Rdata") )
load(file=file.path(dir_data,"sp_tr.Rdata") )

summary(sp_tr)
dim(quadrats_species_cover)


# function to compute standard error
std_err <- function(x) {
  sd(x)/sqrt(length(x))
}

# species cover in quadrats in longer style
quad_sp_long <- quadrats_species_cover %>%
  as.data.frame %>%
  rownames_to_column("Quadrats") %>%
  pivot_longer(cols = !"Quadrats", names_to="Species", values_to = "cover") %>%
  left_join( select(sites_quadrats_info, Quadrats, pH, habitat )  ) %>%
  mutate(pH_hab = paste( substr(pH,1,3), habitat ,sep="_"), .before = Quadrats ) %>%
  select(- pH, - habitat)
# head(quad_sp_long)

# trait values in logn format  
sp_tr_long <- sp_tr %>%
  rownames_to_column("Species") %>%
  mutate_if(is.factor, as.character) %>%
  pivot_longer( cols = !"Species", names_to="Trait", values_to = "Modality" )
# head(sp_tr_long)  
  
# merging then grouping by quadrats and modality to compute total cover
quad_sp_tr <- left_join(quad_sp_long, sp_tr_long, by="Species") %>%
  group_by(pH_hab, Quadrats, Trait, Modality) %>%
  summarize(cover=sum(cover))

head(quad_sp_tr)

# average and standard error
pH_hab_mean_se <- quad_sp_tr %>%
  select(-Quadrats) %>%
  group_by(pH_hab, Trait, Modality) %>%
  summarise( across(.cols = cover, 
                    .fns = list(mean = mean, se = std_err), 
                    .names="{.col}_{.fn}"
                    )
             )

head(pH_hab_mean_se)


