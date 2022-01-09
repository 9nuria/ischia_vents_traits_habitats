# Seb computing FD with new mFD package
# link (not to share at this point, software paper submitted to Ecography)
# https://cmlmagneville.github.io/mFD/

# category abundance for each trait across habitats and pH conditions. This figure and subplots are not considered for the article so far. 

## preparing #####
rm(list=ls())

# libraries ---
library(tidyverse)
library(mFD)


# folder with ready to use data ----
dir_data<-"./data"

# loading dataset built in script A ----
load(file.path(dir_data, "fe_tr.Rdata"))
load(file.path(dir_data, "traits_cat.Rdata"))
load(file.path(dir_data, "quadrats_species_cover.Rdata"))
load(file.path(dir_data, "quadrats_fe_cover.Rdata"))
load(file.path(dir_data,"sites_quadrats_info.Rdata"))

# combinations of hab*pH  ----

# names
hab_ph<-c("shallow_reef_ambient1","shallow_reef_ambient2","shallow_reef_low",
          "cave_ambient1","cave_ambient2","cave_low",
          "reef_ambient1","reef_ambient2", "reef_low",
          "deep_reef_ambient1", "deep_reef_ambient2", "deep_reef_low")

# color code for habitat and pH from Nuria
vcolors<-c("#93a1fa", "#93a1fa","#f7d305",
           "#6478f5","#6478f5", "#f5a511",
           "#3953f7","#3953f7","#f78e0c",
           "#0219ad", "#0219ad", "#f7560c")
names(vcolors)<-hab_ph


# folder to save plot as png ----
dir_plot_trait<-"./plot/traits"


## for each trait computing average cover of each modality in each hab*ph

# list to store results
habph_tr_moddom<-list()

# loop on trait 
for (t in names(fe_tr) ) {
  
  # modalities of trait t
  t_mod<-levels(fe_tr[,t])
  
  # matrix to store results
  t_mod_pcover<-matrix(NA, length(hab_ph) , length(t_mod), 
                       dimnames=list(hab_ph, t_mod) )

  
  # percentage of cover of each modality through double loop
  for (h in hab_ph) {
    
    # quadrats from h
    q_h<-row.names(sites_quadrats_info)[which(sites_quadrats_info[,"habitat_ph"]==h)] 
  
    for (m in t_mod) {
      
      # fe with modality m
      fe_m<-row.names(fe_tr)[which(fe_tr[,t]==m)]
      
      # average cover of m in h
      if( length(fe_m)==1) {  
        t_mod_pcover[h,m]<-mean(quadrats_fe_cover[q_h,fe_m])
      } else {
        t_mod_pcover[h,m]<-mean(apply(quadrats_fe_cover[q_h,fe_m],1,sum))
        }
      
      
    }# end of m
  }# end of h

    # storing
  habph_tr_moddom[[t]]<-t_mod_pcover
  
  # illustrating as stacked barplot after pivoting
  plot_t_mod_pcover<-t_mod_pcover %>% 
    as_tibble(t_mod_pcover, rownames="habitat_pH") %>%
    pivot_longer(cols = t_mod, names_to = "modality", values_to = "cover") %>%
    ggplot( aes(fill=modality, y=cover, x=habitat_pH)) + 
    geom_bar(position="stack", stat="identity") +
    theme(legend.position="top")
  
  # saving
    ggsave(plot_t_mod_pcover, filename=paste0(t,"_habpH_cover.png"), path = dir_plot_trait, width = 14, height=5)
  

}# end of t

##