# Seb computing FD with new mFD package
# link (not to share at this point, software paper submitted to Ecography)
# https://cmlmagneville.github.io/mFD/


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
load(file.path(dir_data, "habph_fe_cover.Rdata"))
load(file.path(dir_data,"sites_quadrats_info.Rdata"))


# total cover
quadrat_covertot<-apply(quadrats_species_cover, 1, sum)
summary(quadrat_covertot)

# number of sp per quadrat
quadrat_nbsp<-apply(quadrats_species_cover, 1, function(x) { length(which(x>0) ) } )
summary(quadrat_nbsp) # from 2 to 27

# number of FE per quadrat
quadrat_nbFE<-apply(quadrats_fe_cover, 1, function(x) { length(which(x>0) ) } )
summary(quadrat_nbFE) # from 2 to 23

# quadrats with at least 5 Fes to be able to compute FRic
quadrat_sup4FE<-names(which(quadrat_nbFE>4))
length(quadrat_sup4FE) # 279 out of the 294 quadrats

## computing Gower distance between FEs ####
fe_fe_dist<-funct.dist(fe_tr, tr_cat = traits_cat, metric = "gower")
range(fe_fe_dist)

## multidimensional space for FD indices ####

# building functional spaces from PCoA ----
fe_fspaces<-quality.fspaces(sp_dist = fe_fe_dist )

# comparing their quality with mean absolute deviation ----
round(fe_fspaces$quality_fspaces,3)
# => lowest mAD is or 4D

# FE coordinates in the 4D funct space ----
fe_4D_coord<-fe_fspaces$details_fspaces$sp_pc_coord[,1:4]


# compute FDis and FIde for all quadrats ---
quadrats_multidimFD<-alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, 
                                       asb_sp_w = quadrats_fe_cover,
                                       ind_vect = c("fdis", "fide", "fspe", "fori"), 
                                       scaling = TRUE, details_returned = TRUE
                                       )

# compute FRic, FDis and FIde for all habitats_pH ---
habph_multidimFD<-alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, 
                                       asb_sp_w = habph_fe_cover,
                                       ind_vect = c("fric", "fdis", "fide", "fdiv"), 
                                       scaling = TRUE, details_returned = TRUE
                                       )


## alpha Hill numbers taxo and FD indices #####

# compute taxonomic Hill numbers on FEs (q=0 for number, q= 1 for Shannon) ----
quadrats_taxo_hill<-alpha.fd.hill(asb_sp_w = quadrats_fe_cover,
              sp_dist = fe_fe_dist,
              q = c(0,1) , tau = "min" ,
              details_returned = FALSE )
colnames(quadrats_taxo_hill)<-c("FE_richness", "FE_shannon")
head(quadrats_taxo_hill)

# compute functional Hill numbers on FEs (q=1, Shannon-like) ----
quadrats_funct_hill<-alpha.fd.hill(asb_sp_w = quadrats_fe_cover,
                                  sp_dist = fe_fe_dist,
                                  q = 1 , tau = "mean" ,
                                  details_returned = FALSE )
head(quadrats_funct_hill)

## merging all diversity indices with quadrats info
quadrats_biodiv<-data.frame(sites_quadrats_info,
                            Total_cover=quadrat_covertot,
                            Nb_sp=quadrat_nbsp,
                            quadrats_taxo_hill, 
                            quadrats_funct_hill,
                            quadrats_multidimFD$functional_diversity_indices[,-1]
                            )
summary(quadrats_biodiv)


## beta Hill numbers taxo and FD indices #####

# compute taxonomic beta with Hill numbers on FEs (q= 1) ----
quadrats_beta_taxo_hill<-beta.fd.hill(asb_sp_w = quadrats_fe_cover,
             sp_dist = fe_fe_dist,
             q = 1 , tau = "min" ,
             details_returned = FALSE )


# compute functional beta with Hill numbers on FEs (q=1) ----
quadrats_beta_funct_hill<-beta.fd.hill(asb_sp_w = quadrats_fe_cover,
                                      sp_dist = fe_fe_dist,
                                      q = 1 , tau = "mean" ,
                                      details_returned = FALSE )

# list with both beta ---
quadrats_beta_hill<-list (taxo_q1= quadrats_beta_taxo_hill$q1 ,
                     funct_q1 = quadrats_beta_funct_hill$q1)

## saving all tables for FD analyses ######

# folder with FD values ----
dir_results<-"./FD"

# saving ----
save(fe_4D_coord, file=file.path(dir_results,"fe_4D_coord.Rdata") )#data for new figure func space with only 2pH levels
save(quadrats_biodiv, file=file.path(dir_results,"quadrats_biodiv.Rdata") )#species richness added as a biodiv metric 
save(quadrats_multidimFD, file=file.path(dir_results,"quadrats_multidimFD.Rdata") )# species richness added as a biodiv metric 
save(quadrats_beta_hill, file=file.path(dir_results,"quadrats_beta_hill.Rdata") )#new plot of trait and cover
save(habph_multidimFD, file=file.path(dir_results,"habph_multidimFD.Rdata") )#species richness added as a biodiv metric




## end of script #####