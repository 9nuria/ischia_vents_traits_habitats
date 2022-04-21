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

dim(quadrats_fe_cover)
dim(habph_fe_cover)

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

## by Seb ####

# computing relative change as (Low pH - Ambient) / Ambient = (Low pH /Ambient) - 1

# complete data frame

df_plot <- df %>%
  mutate( habitat= factor( c( rep("shallow_reef",2), rep("cave",2), rep("reef",2), rep("deep_reef",2) ) ) ) %>%
  mutate(ph=factor( rep(c("amb", "low"),4) ) )

df_plot

# trait: form
k="form"
df_plot_k <- df_plot %>% 
  select( habitat, ph, starts_with(k)) %>%
  pivot_longer(cols=starts_with(k), names_to = "categ", values_to = "cover") %>%
  mutate(categ=str_replace_all(string=categ, pattern=paste0(k,"."), replacement = "") )

#  stacked barplot for pH level with facets as habitat  
plot_k <- ggplot(df_plot_k,
                 aes(x=ph,y=cover,fill=categ))+
  geom_bar(stat = "identity",color="black")+
  facet_wrap(~habitat,nrow=1)

plot_k

## by Nuria##
# colors for categories
colors6<-c("#1b9e77", "#e7298a", "#7570b3", "#d95f02", "#0c2c84", "#c4c0c2")
colors4b<-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
colors3<-c("#1b9e77", "#7570b3", "#e7298a")
colors2<-c("#1b9e77", "#e7298a")
colorsph<-c("#bcd7e6", "#e31a1c")

# trait: form
k="form"
df_plot_k <- df_plot %>% 
  select( habitat, ph, starts_with(k)) %>%
  pivot_longer(cols=starts_with(k), names_to = "categ", values_to = "cover") %>%
  mutate(categ=str_replace_all(string=categ, pattern=paste0(k,"."), replacement = "") )

#order habitats by depth

df_plot_k$habitat<-(factor(df_plot_k$habitat, levels=c("shallow_reef","cave","reef","deep_reef"), labels=c("Shallow reef","Cave","Reef","Deep reef")))
df_plot_k$ph<-(factor(df_plot_k$ph, levels=c("amb","low"), labels=c("ambient","low")))

##### Plot 1 Morphological form
plot1 <- ggplot(df_plot_k,aes(x=ph,y=cover,color=categ, group=categ))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="ambient")-0.6,
                xmax=which(levels(as.factor(ph))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="low")-0.5,
                xmax=which(levels(as.factor(ph))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors4b, name="", labels=c("Encrusting", "Filaments", "Massive", "Tree"))+
  scale_y_continuous(name ="Cover (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(name="", labels= c("", ""))+
  labs(title = "Morphological form")+
  theme_bw()+
  theme (plot.title = element_text(size=14, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~habitat,nrow=1)+(theme(strip.text=element_text(size=14)))

plot1

# trait: feeding
f="feeding"

df_plot_f <- df_plot %>% 
  select( habitat, ph, starts_with(f)) %>%
  pivot_longer(cols=starts_with(f), names_to = "categ", values_to = "cover") %>%
  mutate(categ=str_replace_all(string=categ, pattern=paste0(f,"."), replacement = "") )

#order habitats by depth

df_plot_f$habitat<-(factor(df_plot_f$habitat, levels=c("shallow_reef","cave","reef","deep_reef"), labels=c("Shallow reef","Cave","Reef","Deep reef")))
df_plot_f$ph<-(factor(df_plot_f$ph, levels=c("amb","low"), labels=c("ambient","low")))

plot2 <- ggplot(df_plot_f,aes(x=ph,y=cover,color=categ, group=categ))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="ambient")-0.6,
                xmax=which(levels(as.factor(ph))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="low")-0.5,
                xmax=which(levels(as.factor(ph))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors6, name="", labels=c("Autotroph", "Filter feeder", "Herbivor/Grazer", "Carnivor", "Detritivor", "Parasite"))+
  scale_y_continuous(name ="Cover (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(name="", labels= c("", ""))+
  labs(title = "Feeding")+
  theme_bw()+
  guides(color=guide_legend(ncol=2))+
  theme (plot.title = element_text(size=14, hjust=0.5),
         axis.title.x=element_blank())+
  facet_wrap(~habitat,nrow=1)+(theme(strip.background = element_blank(),
                                     strip.text.x = element_blank()))

plot2
# trait: growth
g="growth"

df_plot_g <- df_plot %>% 
  select( habitat, ph, starts_with(g)) %>%
  pivot_longer(cols=starts_with(g), names_to = "categ", values_to = "cover") %>%
  mutate(categ=str_replace_all(string=categ, pattern=paste0(g,"."), replacement = "") )

#order habitats by depth

df_plot_g$habitat<-(factor(df_plot_g$habitat, levels=c("shallow_reef","cave","reef","deep_reef"), labels=c("Shallow reef","Cave","Reef","Deep reef")))
df_plot_g$ph<-(factor(df_plot_g$ph, levels=c("amb","low"), labels=c("ambient","low")))

plot3 <- ggplot(df_plot_g,aes(x=ph,y=cover,color=categ, group=categ))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="ambient")-0.6,
                xmax=which(levels(as.factor(ph))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="low")-0.5,
                xmax=which(levels(as.factor(ph))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors3, name="", labels=c("Low", "Moderate", "High"))+
  scale_y_continuous(name ="Cover (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(name="", labels= c("", ""))+
  labs(title = "Growth")+
  theme_bw()+
  theme (plot.title = element_text(size=14, hjust=0.5),
         axis.title.x=element_blank())+
  facet_wrap(~habitat,nrow=1)+(theme(strip.background = element_blank(),
                                     strip.text.x = element_blank()))

plot3

# trait: growth
c="calcification"

df_plot_c <- df_plot %>% 
  select( habitat, ph, starts_with(c)) %>%
  pivot_longer(cols=starts_with(c), names_to = "categ", values_to = "cover") %>%
  mutate(categ=str_replace_all(string=categ, pattern=paste0(c,"."), replacement = "") )

#order habitats by depth

df_plot_c$habitat<-(factor(df_plot_c$habitat, levels=c("shallow_reef","cave","reef","deep_reef"), labels=c("Shallow reef","Cave","Reef","Deep reef")))
df_plot_c$ph<-(factor(df_plot_c$ph, levels=c("amb","low"), labels=c("ambient","low")))

plot4 <- ggplot(df_plot_c,aes(x=ph,y=cover,color=categ, group=categ))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="ambient")-0.6,
                xmax=which(levels(as.factor(ph))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(ph))=="low")-0.5,
                xmax=which(levels(as.factor(ph))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="", labels=c("No", "Yes"))+
  scale_y_continuous(name ="Cover (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  labs(title = "Calcification")+
  theme_bw()+
  theme (plot.title = element_text(size=14, hjust=0.5),
         axis.title.x=element_blank())+
  facet_wrap(~habitat,nrow=1)+(theme(strip.background = element_blank(),
                                     strip.text.x = element_blank()))

plot4

#layout 4 plots

all4<-(plot1/plot2/plot3/plot4)
all4

getwd()

ggsave("fig5_functions_cover_4plots.png", plot= all4, device="png", height=25, width=20, units="cm", dpi=300)
