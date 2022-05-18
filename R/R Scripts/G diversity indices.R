#diversity indices across habitats and 2 pH conditions

## preparing #####
rm(list=ls())

# libraries ---
library(tidyverse)
library(mFD)
library(patchwork)


# folder with ready to use data ----
dir_data<-"./data"
dir_results<-"./FD"
# folder to save plot as png ----
dir_plot<-"plot"

# loading dataset built in script A ----
load(file.path(dir_results, "quadrats_biodiv.Rdata"))

# setting parameters for plot ####

quadrats_biodiv$condition<-factor(quadrats_biodiv$condition,levels=c("SRA","SRL","CL","CA","RL","RA","CoL","CoA"))#nu
quadrats_biodiv$condition<- recode_factor (quadrats_biodiv$condition, SRA= "Shallow Reef Ambient pH", SRL="Shallow Reef Low pH", 
                                           CA= "Cave Ambient pH", CL= "Cave Low pH", 
                                           RA= "Reef Ambient pH", RL= "Reef Low pH", 
                                           CoA= "Deep Reef Ambient pH", CoL= "Deep Reef Low pH")

vcolors<-c("#93a1fa", "#f7d305",
           "#6478f5", "#f5a511",
           "#3953f7","#f78e0c",
           "#0219ad", "#f7560c")
levels(quadrats_biodiv$condition)

hab_ph2<-c("Shallow Reef Ambient pH","Shallow Reef Low pH",
           "Cave Ambient pH", "Cave Low pH",
           "Reef Ambient pH", "Reef Low pH",
           "Deep Reef Ambient pH", "Deep Reef Low pH")

hab_ph3<-c("Shallow Reef\nAmbient pH","Shallow Reef\nLow pH",
           "Cave\nAmbient pH", "Cave\nLow pH",
           "Reef\nAmbient pH", "Reef\nLow pH",
           "Deep Reef\nAmbient pH", "Deep Reef\nLow pH")


names(vcolors)<-hab_ph2

########Plots 

# species richness
nbsp_plot<-ggplot(data = quadrats_biodiv, aes(x =condition, y = Nb_sp, color= condition)) +
  geom_boxplot(color="black", outlier.shape=NA, alpha=0.2, fill=vcolors)+
  geom_jitter(aes(color = condition), size=2, shape=16, alpha=0.5, position=position_jitter(0.3)) +
  scale_color_manual(values=vcolors)+
  ylab("Species richness") +
  ggtitle("A) Species richness")+
  scale_y_continuous(limits= c(0, 30), breaks=c(0, 10,20,30))+
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold"),
        #axis.text.x= element_text(size= 12, angle=15),
        axis.text.x=element_blank(),
        axis.text.y= element_text (face="plain", size= 14),
        axis.title.y = element_text (face="bold", size= 14),
        legend.position = "none",
        axis.title.x = element_blank())
nbsp_plot

# FE richness
fe_plot<-ggplot(data = quadrats_biodiv, aes(x =condition, y = FE_richness, color= condition)) +
  geom_boxplot(color="black", outlier.shape=NA, alpha=0.2, fill=vcolors)+
  geom_jitter(aes(color = condition), size=2, shape=16, alpha=0.5, position=position_jitter(0.3)) +
  scale_color_manual(values=vcolors)+
  ylab("FE richness\n(# groups)") +
  ggtitle("B) Functional entities richness")+
  theme_bw() +
  scale_y_continuous(limits= c(0, 30), breaks=c(0, 10,20,30))+
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold"),
        #axis.text.x= element_text(size= 12, angle=15),
        axis.text.x=element_blank(),
        axis.text.y= element_text (face="plain", size= 14),
        axis.title.y = element_text (face="bold", size= 14),
        legend.position = "none",
        axis.title.x = element_blank())
fe_plot

# multidim funct dispersion

fdisp_plot<-ggplot(data = quadrats_biodiv, aes(x =condition, y = fdis, color= condition)) +
  geom_boxplot(color="black", outlier.shape=NA, alpha=0.2, fill=vcolors)+
  geom_jitter(aes(color = condition), size=2, shape=16, alpha=0.5, position=position_jitter(0.3)) +
  scale_color_manual(values=vcolors)+
  ylab("FDis") +
  ggtitle("C) Functional dispersion")+
  theme_bw() +
  scale_y_continuous(limits= c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8))+
  scale_x_discrete(labels=hab_ph3)+
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold"),
        axis.text.x= element_text(size= 14,  vjust=0.5, face="bold"),
        axis.text.y= element_text (face="plain", size= 14),
        legend.position = "none",
        axis.title.y = element_text (face="bold", size= 14),
        axis.title.x = element_blank())
fdisp_plot

# assembling 
boxplot<-nbsp_plot / fe_plot / fdisp_plot 
boxplot

ggsave("fig4_indicesdiv.png", plot= boxplot, path = dir_plot, device="png", height=35, width=35, units="cm", dpi=300)

############stats glms

library(lme4)

quadrats_biodiv<-quadrats_biodiv %>%
  mutate(ph2= recode_factor(pH, 
                            "ambient1" = "ambient", "ambient2"= "ambient" , "low"="low"))

quadrats_biodiv$ph2<-factor(quadrats_biodiv$ph2, ordered=FALSE, levels= c("ambient", "low"))
quadrats_biodiv$habitat<-factor(quadrats_biodiv$habitat, ordered=FALSE, levels= c("shallow_reef", "cave", "reef", "deep_reef"))

#design: habitat 4 levels, pH conditions 2 levels
str(quadrats_biodiv)

m1<-glm(Nb_sp ~ habitat*ph2, family=poisson, quadrats_biodiv)
m1

summary(m1)

m2<-glm(FE_richness ~ habitat*ph2, family=poisson, quadrats_biodiv)
m2

summary(m2)

m3<-glm(FE_richness ~ habitat*ph2, family=gaussian(), quadrats_biodiv)
m3
summary(m3)



