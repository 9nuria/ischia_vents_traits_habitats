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
load(file.path(dir_data, "fe_tr.Rdata"))
load(file.path(dir_data, "traits_cat.Rdata"))
load(file.path(dir_data, "quadrats_species_cover.Rdata"))
load(file.path(dir_data, "quadrats_fe_cover.Rdata"))
load(file.path(dir_data,"sites_quadrats_info.Rdata"))

load(file.path(dir_data, "habph_fe_cover.Rdata"))


# computing mean cover of FEs in the 2 pH levels: Ambient and Low pH
habph2_fe_cover<-rbind( 
  cave_low = habph_fe_cover["cave_low",],
  cave_amb = apply(  habph_fe_cover[c("cave_ambient1", "cave_ambient2"),], 2, mean),
  
  deep_reef_low = habph_fe_cover["deep_reef_low",],
  deep_reef_amb = apply(  habph_fe_cover[c("deep_reef_ambient1", "deep_reef_ambient2"),], 2, mean),
  
  reef_low = habph_fe_cover["reef_low",],
  reef_amb = apply(  habph_fe_cover[c("reef_ambient1", "reef_ambient2"),], 2, mean),
  
  shallow_reef_low = habph_fe_cover["shallow_reef_low",],
  shallow_reef_amb = apply(  habph_fe_cover[c("shallow_reef_ambient1", "shallow_reef_ambient2"),], 2, mean)
  
)

#
sites_quadrats_info<-sites_quadrats_info %>%
  mutate(habitat_ph2= recode_factor(habitat_ph,
                            "cave_ambient1" = "cave_amb", "cave_ambient2"= "cave_amb" , "cave_low"="cave_low",
                            "deep_reef_ambient1"="deep_reef_amb", "deep_reef_ambient2"="deep_reef_amb", "deep_reef_low"="deep_reef_low", 
                            "reef_ambient1"= "reef_amb", "reef_ambient2"="reef_amb", "reef_low"="reef_low",  
                            "shallow_reef_ambient1"="shallow_reef_amb", "shallow_reef_ambient2"="shallow_reef_amb", "shallow_reef_low"= "shallow_reef_low"))





str(sites_quadrats_info$habitat_ph2)

# setting parameters for plot ####

# folder to save plot as png ----
dir_plot<-"./plot"
root_dir<-getwd()
setwd(dir_plot)

# color code for habitat and pH from Nuria ----

hab_ph2<-c("shallow_reef_amb","shallow_reef_low",
           "cave_amb", "cave_low",
           "reef_amb", "reef_low",
           "deep_reef_amb", "deep_reef_low")

vcolors<-c("#93a1fa", "#f7d305",
           "#6478f5", "#f5a511",
           "#3953f7", "#f78e0c",
           "#0219ad", "#f7560c")
names(vcolors)<-hab_ph2

## for each trait computing average cover of each modality in each hab*ph
###see script D trait dominance


# list to store results
habph2_tr_moddom<-list()

# loop on trait 
for (t in names(fe_tr) ) {
  
  # modalities of trait t
  t_mod<-levels(fe_tr[,t])
  
  # matrix to store results
  t_mod_pcover<-matrix(NA, length(hab_ph2) , length(t_mod), 
                       dimnames=list(hab_ph2, t_mod) )
  
  
  # percentage of cover of each modality through double loop
  for (h in hab_ph2) {
    
    # quadrats from h
    q_h<-row.names(sites_quadrats_info)[which(sites_quadrats_info[,"habitat_ph2"]==h)] 
    
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
  habph2_tr_moddom[[t]]<-t_mod_pcover
}  

#list as dataframe

df<- as.data.frame(habph2_tr_moddom)

# computing delta (response change) as Ambient- Low pH

delta.df<-rbind(
  shallow_reef=apply(df[c("shallow_reef_low", "shallow_reef_amb"),], 2, diff),
  cave = apply(df[c("cave_low", "cave_amb"),], 2, diff),
  reef = apply(df[c("reef_low", "reef_amb"),], 2, diff),
  deep_reef= apply( df[c("deep_reef_low", "deep_reef_amb"),], 2, diff))

delta.df<-as.data.frame(delta.df)
delta.df$habitat <- rownames(delta.df)


delta.dflong<-delta.df%>%
  tidyr::gather(trait, delta, form.a:chem.2)

delta.dflong$habitat<-factor(delta.dflong$habitat, ordered=TRUE, levels= c("shallow_reef", "cave", "reef", "deep_reef"))

#Figure


figall<-ggplot(delta.dflong, aes(trait, delta, fill = trait)) +
  geom_bar(stat = "identity", position=position_dodge(width=0.8, preserve = "single"), alpha=0.5) +
  # ylab("Cover trait change")+
  geom_hline(yintercept = 0) +
  scale_fill_discrete()+
  #  scale_y_continuous(breaks = 0:nlevels(DF$variable)) +
  theme_light() +
  theme(axis.text.y= element_text (face="bold", size= 14),
        axis.title.y = element_blank(),
        axis.text.x=element_blank(), 
        axis.title.x = element_blank(),
        legend.position="bottom",
        legend.text=element_text(size=14, face="bold"))+
  facet_grid(habitat~., 
             labeller=labeller(habitat=c("shallow_reef"= "Shallow Reef", "cave"="Cave", "reef"="Reef", "deep_reef"= "Deep Reef")))+
  theme(strip.text = element_text(size=20, color="white", face="bold"))
figall

ggsave("fig5_trait_change.png", plot= figall, path = dir_plot, device="png", height=35, width=35,  dpi=300)

#https://ggplot2-book.org/facet.html


