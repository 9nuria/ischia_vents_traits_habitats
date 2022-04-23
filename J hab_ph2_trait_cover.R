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

#include ph and habitat ggplot . Nuria
quad_sp_tr<-quad_sp_tr%>%
  mutate(ph= recode_factor(pH_hab,
                           "amb_cave"= "amb", "amb_deep_reef"="amb", "amb_reef"="amb", "amb_shallow_reef"="amb", 
                           "low_cave"= "low","low_deep_reef"="low", "low_reef" = "low", "low_shallow_reef"="low"),
         habitat=recode_factor(pH_hab,
                               "amb_cave"= "cave", "amb_deep_reef"="deep_reef", "amb_reef"="reef", "amb_shallow_reef"="shallow_reef", 
                               "low_cave"= "cave","low_deep_reef"="deep_reef", "low_reef" = "reef", "low_shallow_reef"="shallow_reef"))


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

# colors for categories

colors6<-c("#1b9e77", "#e7298a", "#7570b3", "#d95f02", "#0c2c84", "#c4c0c2")
colors4b<-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
colors3<-c("#1b9e77", "#7570b3", "#e7298a")
colors2<-c("#1b9e77", "#e7298a")
colorsph<-c("#bcd7e6", "#e31a1c")

#order habitats by depth

pH_hab_mean_se$habitat<-(factor(pH_hab_mean_se$habitat, levels=c("shallow_reef","cave","reef","deep_reef"), labels=c("Shallow reef","Cave","Reef","Deep reef")))
pH_hab_mean_se$ph<-(factor(pH_hab_mean_se$ph, levels=c("amb","low"), labels=c("ambient","low")))

# trait: form

form.pH_hab_mean_se<- pH_hab_mean_se %>% 
  filter(Trait=="form")

plot1 <- ggplot(form.pH_hab_mean_se,aes(x=ph,y=cover_mean,color=Modality, group=Modality))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=cover_mean-cover_se, ymax=cover_mean+cover_se), width=.2,
                position=position_dodge(0.05))+
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

feed.pH_hab_mean_se<- pH_hab_mean_se %>% 
  filter(Trait=="feeding")

plot2 <- ggplot(feed.pH_hab_mean_se,aes(x=ph,y=cover_mean,color=Modality, group=Modality))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=cover_mean-cover_se, ymax=cover_mean+cover_se), width=.2,
                position=position_dodge(0.05))+
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
  facet_wrap(~habitat,nrow=1)+(theme(strip.text=element_text(size=14)))

plot2

# trait: growth
growth.pH_hab_mean_se<- pH_hab_mean_se %>% 
  filter(Trait=="growth")

plot3 <- ggplot(growth.pH_hab_mean_se,aes(x=ph,y=cover_mean,color=Modality, group=Modality))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=cover_mean-cover_se, ymax=cover_mean+cover_se), width=.2,
                position=position_dodge(0.05))+
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
  facet_wrap(~habitat,nrow=1)+(theme(strip.text=element_text(size=14)))

plot3

#trait calcification
cal.pH_hab_mean_se<- pH_hab_mean_se %>% 
  filter(Trait=="calcification")

plot4 <- ggplot(cal.pH_hab_mean_se,aes(x=ph,y=cover_mean,color=Modality, group=Modality))+
  geom_line(linetype="dashed")+
  geom_point(size=5)+
  geom_errorbar(aes(ymin=cover_mean-cover_se, ymax=cover_mean+cover_se), width=.2,
                position=position_dodge(0.05))+
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
  facet_wrap(~habitat,nrow=1)+(theme(strip.text=element_text(size=14)))

plot4

#layout 4 plots

all4<-(plot1/plot2/plot3/plot4)
all4

# setting parameters for plot ####

# folder to save plot as png ----
dir_plot<-"./plot"
root_dir<-getwd()
setwd(dir_plot)
getwd()

ggsave("fig5_functions_cover_se_4plots.png", plot= all4, device="png", height=25, width=20, units="cm", dpi=300)
