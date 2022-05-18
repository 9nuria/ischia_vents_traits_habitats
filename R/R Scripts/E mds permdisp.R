
# Script to plot MDS on beta taxonomic and functional diversity with the Hill number framework 
# Script to calculate multivariaye homogeneity of group variances for funct and taxonomic beta- diversity 

#library load
library('mFD')
library('reshape2')
library('geometry')
library('ape')
library('ggplot2')
library('patchwork')
library('vegan')
library('MASS')
library('tidyverse')
library('devtools')


# folder with ready to use data ----
dir_data<-"./data"
dir_plot<-"./plot"

# loading dataset built in script B ----
load(file.path(dir_data, "quadrats_beta_hill.Rdata"))
load(file.path(dir_data, "sites_quadrats_info.Rdata"))

####
#Pcoa functional

pcoa_fun <- as.data.frame(cmdscale(quadrats_beta_hill$funct_q1, 4))# number axes
#cmdscale Classical multidimensional scaling (MDS) of a data matrix

info<-sites_quadrats_info%>%
  dplyr::select(Quadrats,condition,Description.condition, pH, habitat)

pcoa_fun<- pcoa_fun %>%
  bind_cols(info)

#Pcoa taxonomic


pcoa_taxo <- as.data.frame(cmdscale(quadrats_beta_hill$taxo_q1, 4))
str(pcoa_taxo)

pcoa_taxo<- pcoa_taxo %>%
  bind_cols(info)

#Eigenvalues

x.func <- cmdscale(quadrats_beta_hill$funct_q1, eig=TRUE) # coordinates points, eigenvalues= eig

cumsum(x.func$eig[x.func$eig>=0]) / sum(x.func$eig[x.func$eig>0])

# 0.6242943 0.7218448 0.7750432 0.8124273

x.tax<- cmdscale(quadrats_beta_hill$taxo_q1,eig=TRUE)
cumsum(x.tax$eig[x.tax$eig>=0]) / sum(x.tax$eig[x.tax$eig>0])
#0.2394259 0.3809141 0.4459513 0.5033202
###


########### MDS figure  1 legend


str(pcoa_fun$condition)
pcoa_fun$condition<-as.factor(pcoa_fun$condition)

pcoa_fun$condition<-factor(pcoa_fun$condition,levels=c("SRA","SRL","CL","CA","RL","RA","CoL","CoA"))

pcoa_fun$condition<- recode_factor (pcoa_fun$condition, SRA= "Shallow Reef Ambient pH", SRL="Shallow Reef Low pH", CL= "Cave Low pH", CA= "Cave Ambient pH", RL= "Reef Low pH", RA= "Reef Ambient pH",
                                    CoL= "Deep Reef Low pH", CoA= "Deep Reef Ambient pH")

levels(pcoa_fun$condition)


#1 legend-functional 

fig.fun.v1v2<-ggplot(pcoa_fun, aes(x=V1, y=V2, group=condition, fill=condition, color=condition, shape=condition)) + 
  geom_point(size=6, aes(color=condition, shape=condition)) +
  scale_colour_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_fill_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_shape_manual(values= c(21,21,24,24,22,22,23,23))+
  labs(title="B) Functional diversity", x="V1", y="V2") +
  scale_x_continuous(limits=c(-0.55, 0.6))+
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  theme_light(base_size=20)+
  theme(legend.title=element_blank(),
        legend.position="bottom")

fig.fun.v1v2

fig.fun.v3v4<-ggplot(pcoa_fun, aes(x=V3, y=V4, group=condition, fill=condition, color=condition)) + 
  geom_point(size=6, aes(color=condition, shape=condition)) +
  scale_colour_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_fill_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_shape_manual(values= c(21,21,24,24,22,22,23,23))+
  labs(title="B) Functional diversity", x="V3", y="V4") +
  scale_x_continuous(limits=c(-0.55, 0.6))+
  scale_y_continuous(limits=c(-0.55, 0.55)) +
  theme_light(base_size=20)+
  theme(legend.title=element_blank(),
        legend.position="bottom")

fig.fun.v3v4



#ggsave("figures/fig5_mds_funct_subset_traits.png", plot= fig.fun2.def, device="png", height=25, width=35, units="cm", dpi=300)



###############
#taxonomic

pcoa_taxo$condition<-factor(pcoa_taxo$condition,levels=c("SRA","SRL","CL","CA","RL","RA","CoL","CoA"))

str(pcoa_taxo$condition)

pcoa_taxo$condition<- recode_factor (pcoa_taxo$condition, SRA= "Shallow Reef Ambient pH", SRL="Shallow Reef Low pH", CL= "Cave Low pH", CA= "Cave Ambient pH", RL= "Reef Low pH", RA= "Reef Ambient pH",
                                    CoL= "Deep Reef Low pH", CoA= "Deep Reef Ambient pH")

levels(pcoa_taxo$condition)


fig.tax.v1v2<-ggplot(pcoa_taxo, aes(x=V1, y=V2, group=condition, fill=condition, color=condition, shape=condition)) + 
  geom_point(size=6) +
  scale_colour_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_fill_manual( values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_shape_manual(values= c(21,21,24,24,22,22,23,23))+
  #scale_fill_discrete(name= "Habitat/pH condition", labels= c("SRA","SRL","CL","CA","RL","RA","CoL","CoA"))+
  labs(title="A) Taxonomic diversity", x="V1", y="V2") +
  scale_x_continuous(limits=c(-0.55, 0.6))+
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  theme_light(base_size=20)+
  theme(legend.position='none')

fig.tax.v1v2

fig.tax.v3v4<-ggplot(pcoa_taxo, aes(x=V3, y=V4, group=condition, fill=condition, color=condition, shape=condition)) + 
  geom_point(size=6) +
  scale_colour_manual(values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_fill_manual( values = c("#93a1fa", "#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad"))+
  scale_shape_manual(values= c(21,21,24,24,22,22,23,23))+
  #scale_fill_discrete(name= "Habitat/pH condition", labels= c("SRA","SRL","CL","CA","RL","RA","CoL","CoA"))+
  labs(title="A) Taxonomic diversity", x="V3", y="V4") +
  scale_x_continuous(limits=c(-0.55, 0.6))+
  scale_y_continuous(limits=c(-0.55, 0.55)) +
  theme_light(base_size=20)+
  theme(legend.position='none')

fig.tax.v3v4


####figure MDS combined: taxonomic and functional

mdsv1v2<-(fig.tax.v1v2/fig.fun.v1v2)
mdsv1v2

mdsv3v4<-(fig.tax.v3v4/fig.fun.v3v4)
mdsv3v4

ggsave("fig3_mds_funct_taxo_v1v2.png", plot= mdsv1v2, path = dir_plot, device="png", height=35, width=35, units="cm", dpi=300)

ggsave("fig_mds_funct_taxo_sm_v3v4.png", plot= mdsv3v4, path = dir_plot, device="png", height=35, width=35, units="cm", dpi=300)


### 4 (habitats) Permdisp, 4 permanova, and pairwise comparisons for funct and taxo

##taxonomic

#1. to calculate the multivariate dispersion (variances; average distance to centroids)
#habitat shallow reefs



# shallow reefs, functional

beta_df<-mFD::dist.to.df(quadrats_beta_hill)
beta_df.fun<-beta_df%>%
  dplyr::select(-taxo_q1)

habitat <- c("SRA",  "SRL")
quad <- info[info$condition %in% habitat,]$Quadrats #select the quadrats within the habitat vector 1s1, 2s1 etc...
divf <- beta_df.fun[beta_df.fun$x1 %in% quad,] #selection of quadrats Shallow reefs for asb.1
divf <- divf[divf$x2 %in% quad,]#selection of quadrats Shallow reefs for asb.2

source(file = "long_to_wide_dist.R")
beta_long_fun <- long_to_wide_distance(divf)
str(divf)


groupf <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})


#https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/permutest.betadisper
#1) To calculate Permidsp (multivariate dispersion, variance)>permidsp
modf <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE) #bias.adjust= TRUE for unequal number of samples

#2) To test if the dispersions (variance) are different> Anova or use a permutation test to generate 
#a permutation distribution of F under the Null hypothesis of no difference in dispersion between groups> permutest.betadisper

anova(modf)

## Permutation test for F
permutest(modf, permutations = 999)

plot(modf)
boxplot(modf)


# shallow reefs, taxonomic
beta_df.taxo<-beta_df%>%
  dplyr::select(-funct_q1)

divtaxo <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] #selection of quadrats Shallow reefs for asb.1
divtaxo <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})

modtaxo <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo)
permutest(modtaxo, permutations = 999)

plot(modtaxo)
boxplot(modtaxo)
####caves, funct
habitat <- c("CA",  "CL")
quad <- info[info$condition %in% habitat,]$Quadrats #select the quadrats within the habitat vector 1s1, 2s1 etc...
divf <- beta_df.fun[beta_df.fun$x1 %in% quad,] #selection of quadrats Shallow reefs for x1
divf <- divf[divf$x2 %in% quad,]#selection of quadrats Shallow reefs for x2

beta_long_fun <- long_to_wide_distance(divf)
str(divf)

groupf <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})

modf <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf)

permutest(modf, permutations = 999, pairwise=TRUE)

plot(modf)
boxplot(modf)

####caves, taxo

divtaxo <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] #selection of quadrats for x1
divtaxo <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})

modtaxo <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo)
permutest(modtaxo, permutations = 999)

plot(modtaxo)
boxplot(modtaxo)

####reefs, funct

habitat <- c("RA",  "RL")
quad <- info[info$condition %in% habitat,]$Quadrats #select the quadrats within the habitat vector 1s1, 2s1 etc...
divf <- beta_df.fun[beta_df.fun$x1 %in% quad,] #selection of quadrats for x1
divf <- divf[divf$x2 %in% quad,]#selection of quadrats  for x2

beta_long_fun <- long_to_wide_distance(divf)
str(divf)

groupf <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})

modf <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf)

permutest(modf, permutations = 999, pairwise=TRUE)

plot(modf)
boxplot(modf)

####reefs, taxo 

divtaxo <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] #selection of quadrats for x1
divtaxo <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})

modtaxo <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo)
permutest(modtaxo, permutations = 999)

plot(modtaxo)
boxplot(modtaxo)

#### deep reefs, funct

habitat <- c("CoA",  "CoL")
quad <- info[info$condition %in% habitat,]$Quadrats #select the quadrats within the habitat vector 1s1, 2s1 etc...
divf <- beta_df.fun[beta_df.fun$x1 %in% quad,] #selection of quadrats for x1
divf <- divf[divf$x2 %in% quad,]#selection of quadrats  for x2

beta_long_fun <- long_to_wide_distance(divf)
str(divf)

groupf <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})

modf <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf)

permutest(modf, permutations = 999, pairwise=TRUE)

plot(modf)
boxplot(modf)

#### deep reefs, taxo 

divtaxo <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] #selection of quadrats for x1
divtaxo <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})

modtaxo <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo)
permutest(modtaxo, permutations = 999)

plot(modtaxo)
boxplot(modtaxo)


###################################

