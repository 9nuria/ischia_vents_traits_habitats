## This script is to prepare datasets for analyzes using data sent by Nuria on 2021-05-05 #####

## preparing #####
rm(list=ls())
library(tidyverse)
library(mFD)

# folder to load raw data 
dir_raw_data<-"./data/from_nuria"

# folder to save ready to use data
dir_data<-"./data/"


## information about sites and quadrats ####
sites_quadrats_info<-read.csv(file.path(dir_raw_data,"Sites_Quadrats.csv") )
head(sites_quadrats_info)
row.names(sites_quadrats_info)<-sites_quadrats_info$Quadrats

# new variable merging habitat type with pH level
sites_quadrats_info$habitat_ph<-paste(sites_quadrats_info$habitat, 
                                      sites_quadrats_info$pH, 
                                      sep="_")
# names of these combinations
habph<-unique(sites_quadrats_info$habitat_ph) # 12 levels
habph

# number of quadrats per habitat_pH
habph_nbquadrats<-sites_quadrats_info %>%
  group_by(habitat_ph) %>%
  summarize(Ntot=n() )
habph_nbquadrats # highly variable from 12 to 54


## matrices of species cover in sites ####
# from each matrix to a long table
dat_cast<-read.csv(file.path(dir_raw_data,"Data_Cover_tCastello.csv") ) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")
head(dat_cast)

dat_cora<-read.csv(file.path(dir_raw_data,"Data_Cover_tCoralligenous.csv") )%>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")

dat_chia<-read.csv(file.path(dir_raw_data,"Data_Cover_tChiane.csv") ) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")

dat_cave<-read.csv(file.path(dir_raw_data,"Data_Cover_tCaves.csv") )%>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")

# row binding tables and going back to wide table then matrix and replacing NA with 0
quadrats_species_cover0<-rbind(dat_cast, dat_cora, dat_chia, dat_cave) %>%
                    pivot_wider(names_from = species, values_from = cover) 

quadrats_species_cover<-as.matrix( select(quadrats_species_cover0,  !X) )
row.names(quadrats_species_cover)<-quadrats_species_cover0$X
quadrats_species_cover[is.na(quadrats_species_cover)]<-0

summary( apply(quadrats_species_cover,1,sum) ) # OK
dim(quadrats_species_cover) # 294 quadrats and 232 species

# checking all species present at least once
species_sumcover<-apply (quadrats_species_cover, 2, sum)
range(species_sumcover)
species_sumcover[which(species_sumcover==0)] # 7 species absent from all quadrats

# removing those absent species and reordering quadrats as in sites_quadrats
quadrats_species_cover<-quadrats_species_cover[sites_quadrats_info$Quadrats,
                                               names(which(species_sumcover>0) )]
dim(quadrats_species_cover) # 294 quadrats and 225 species



# species trait values and FEs ####

# loading data from csv ----
sp_tr<-read.csv(file.path(dir_raw_data,"species_traits.csv") )
nrow(sp_tr) # 233 species

# species in quadrats with no trait values ---  
species_notraits<-colnames(quadrats_species_cover)[ ! (colnames(quadrats_species_cover) %in% sp_tr$Species) ]

species_notraits # => 9 taxa + bare substrata

# removing them from quadrat*species cover matrix
quadrats_species_cover<-quadrats_species_cover[, 
                        colnames(quadrats_species_cover) %in% sp_tr$Species]
dim(quadrats_species_cover)

# names of species
sp_nm<-colnames(quadrats_species_cover )
length(sp_nm) # 215 species


# filtering absent species from trait dataframe
sp_tr<-sp_tr[sp_tr$Species %in% sp_nm,]

# table with trait coding (Q = categorical , O = ordinal) ----
traits_cat<-data.frame( trait_name= c("form", "feeding", "growth",
                                      "calcification", "mobility", 
                                      "agerepromaturity", "chem") ,
                        trait_type= c("N", "N", "O",
                                      "N", "O", 
                                      "O", "O")
)


# species traits dataframe with correct coding of variables ----
row.names(sp_tr)<-sp_tr$Species
sp_tr<-sp_tr[,traits_cat$trait_name] 


sp_tr$form<-as.factor(sp_tr$form)
levels(sp_tr$form) # both "c" and "c " => correcting
sp_tr$form<-as.factor(gsub(as.character(sp_tr$form),pattern = " ", replacement = "" ) )
levels(sp_tr$form) # OK

sp_tr$feeding<-as.factor(sp_tr$feeding)
levels(sp_tr$feeding) # OK

sp_tr$calcification<-as.factor(sp_tr$calcification)
levels(sp_tr$calcification) # OK

sp_tr$growth<-as.ordered(sp_tr$growth)
levels(sp_tr$growth) # OK

sp_tr$mobility<-as.ordered(sp_tr$mobility)
levels(sp_tr$mobility) # OK

sp_tr$agerepromaturity<-as.ordered(sp_tr$agerepromaturity)
levels(sp_tr$agerepromaturity) # OK

sp_tr$chem<-as.ordered(sp_tr$chem)
levels(sp_tr$chem) # OK


# clustering species into FE ----
sp_to_fe<-mFD::sp.to.fe(sp_tr = sp_tr, tr_cat = traits_cat)

# looking at number of species per FE
summary(sp_to_fe$fe_nb_sp) # => from 1 to 14 (median= 2)

# names and number of FEs
fe_nm<-unique(sp_to_fe$fe_nm)
length(fe_nm) # 74 FE

# list of species in each FE
fe_sp<-list()
for (k in fe_nm) {
  fe_sp[[k]]<-names( sp_to_fe$sp_fe[which(sp_to_fe$sp_fe==k)] )
}# end of k

# trait values of FE
fe_tr<-sp_to_fe$fe_tr
head(fe_tr)

# FE trait values and cover in quadrats ----
quadrats_fe_cover<-matrix(0, 
                            nrow = nrow(quadrats_species_cover), 
                            ncol=length(fe_nm),
                            dimnames=list( row.names(quadrats_species_cover), fe_nm ) )

for (k in fe_nm)
{
  sp_k<-fe_sp[[k]]
  
  if (length(sp_k)==1) {
    quadrats_fe_cover[,k]<-quadrats_species_cover[, sp_k ]
    
  } else {
    quadrats_fe_cover[,k]<-apply( quadrats_species_cover[, sp_k ], 1, sum)
  }
    
}# end of k


## average cover of species and FE in each habitat_ph ####
habph_species_cover<-matrix(0, 
                            nrow = length(habph), 
                            ncol=ncol(quadrats_species_cover),
                            dimnames=list( habph, colnames(quadrats_species_cover) ) )
habph_fe_cover<-matrix(0, 
                            nrow = length(habph), 
                            ncol=ncol(quadrats_fe_cover),
                            dimnames=list( habph, colnames(quadrats_fe_cover) ) )

for (k in habph)
{
  quad_k<-sites_quadrats_info[which(sites_quadrats_info$habitat_ph==k),"Quadrats"]
  habph_species_cover[k,]<-apply( quadrats_species_cover[quad_k,], 2, mean)
  habph_fe_cover[k,]<-apply( quadrats_fe_cover[quad_k,], 2, mean)
}# end of k

dim(habph_species_cover)

## saving all tables for FD analyses ######
save(sites_quadrats_info, file=file.path(dir_data,"sites_quadrats_info.Rdata") )
save(quadrats_species_cover, file=file.path(dir_data,"quadrats_species_cover.Rdata") )
save(habph_species_cover, file=file.path(dir_data,"habph_species_cover.Rdata") )
save(sp_tr, file=file.path(dir_data,"sp_tr.Rdata") )
save(traits_cat, file=file.path(dir_data,"traits_cat.Rdata") )
save(sp_to_fe, file=file.path(dir_data,"sp_to_fe.Rdata") )
save(fe_tr, file=file.path(dir_data,"fe_tr.Rdata") )
save(fe_sp, file=file.path(dir_data,"fe_sp.Rdata") )
save(quadrats_fe_cover, file=file.path(dir_data,"quadrats_fe_cover.Rdata") )
save(habph_fe_cover, file=file.path(dir_data,"habph_fe_cover.Rdata") )

## end of script #####


