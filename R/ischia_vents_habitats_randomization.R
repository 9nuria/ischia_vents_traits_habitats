# Functional biodiversity shifts within and across benthic habitats under ocean acidification ----------------------
## Núria Teixidó*, Enric Ballesteros, Samir Alliouane, Maria Cristina Gambi, Jean-Pierre Gattuso, Kristy 
## Kroeker, Fiorenza Micheli, Alice Mirasole, Sebastien Villéger, Cinzia De Vittor, Valeriano Parravacini 
## *corresponding author. Email: nuria.teixido@imev-mer.fr; nuria.teixido@szn.it 

rm(list=ls()) ; options(mc.cores = parallel::detectCores()) ; setwd("..")

## Loading packages and data ---------------------------------------------------------------------------------------

# Packages
library(ape)
library(brms)
library(devtools)
library(geometry)
library(lme4)
library(MASS)
library(Matrix) 
library(mFD)
library(patchwork)
library(Rcpp)
library(reshape2)
library(tidyverse)
library(vegan)

# Directories
dir_raw_data   <- "./data/1 – raw data"              # folder to load raw data 
dir_data       <- "./data/2 – data generated"        # folder to save ready to use data
dir_plot       <- "./outputs/plot"                   # folder with plot 
dir_plot_trait <- "./outputs/plot/traits"            # folder to save plot as png
dir_model      <- "./data/3 – model"                 # folder to save models

# Generic functions
std_err <- function(x) {sd(x)/sqrt(length(x))}       # function to compute standard error
long_to_wide_distance <- function (df) {             # function to convert distances
  matNames <- sort(unique(as.character(unlist(df[1:2])))) ; colnames(df)[3] <- "value"
  myMat <- matrix(0, length(matNames), length(matNames), dimnames = list(matNames, matNames))
  myMat[as.matrix(df[c(1,2)])] <- df$value ; myMat[as.matrix(df[c(2,1)])] <- df$value
  return(as.dist(myMat))}

# Raw data
sites_quadrats_info <- read.csv(file.path(dir_raw_data,"Sites_Quadrats.csv"))
sp_tr               <- read.csv(file.path(dir_raw_data,"species_traits.csv"))
dat_cast            <- read.csv(file.path(dir_raw_data,"Data_Cover_tCastello.csv")) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")
dat_cora            <- read.csv(file.path(dir_raw_data,"Data_Cover_tCoralligenous.csv")) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")
dat_chia            <- read.csv(file.path(dir_raw_data,"Data_Cover_tChiane.csv")) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")
dat_cave            <- read.csv(file.path(dir_raw_data,"Data_Cover_tCaves.csv")) %>% 
  pivot_longer(cols = ! X, names_to = "species" , values_to = "cover")

# Generated data
load(file.path(dir_data, "fe_tr.Rdata"))
load(file.path(dir_data, "traits_cat.Rdata"))
load(file.path(dir_data, "quadrats_species_cover.Rdata"))
load(file.path(dir_data, "quadrats_fe_cover.Rdata"))
load(file.path(dir_data, "quadrats_beta_hill.Rdata"))
load(file.path(dir_data, "habph_fe_cover.Rdata"))
load(file.path(dir_data, "sites_quadrats_info.Rdata"))
load(file.path(dir_data, "Predicted_values.RData")) 

# Model 
load(file.path(dir_model,"mn.RData"))

# Number of iterations you desire
n = 100

## Data preparation ------------------------------------------------------------------------------------------------
# SCRIPT A ---------------------------------------------------------------------------------------------------------

row.names(sites_quadrats_info) <- sites_quadrats_info$Quadrats

# New variable merging habitat type with pH level
sites_quadrats_info$habitat_ph <- paste(sites_quadrats_info$habitat, sites_quadrats_info$pH, sep="_")

# Names of these combinations
habph <- unique(sites_quadrats_info$habitat_ph) # 12 levels

# Number of quadrats per habitat_pH
(habph_nbquadrats <- sites_quadrats_info %>% group_by(habitat_ph) %>% summarize(Ntot=n())) # highly variable (12–54)

## Randomization ---------------------------------------------------------------------------------------------------
# SCRIPT M ---------------------------------------------------------------------------------------------------------

Infos_FEs = sites_quadrats_info %>% dplyr::select(., Quadrats, pH, habitat)
Infos_FEs$pH[which(Infos_FEs$pH == "ambient1")] = "ambient"
Infos_FEs$pH[which(Infos_FEs$pH == "ambient2")] = "ambient"
Infos_FEs_split = Infos_FEs %>% mutate(., condition = paste(habitat, pH, sep = "_")) %>% group_split(condition)

## Starting the randomization process
# Empty vectors of random sampling selection for each condition
Sample_low_cave      <- vector("list", length = n) ; Sample_amb_cave      <- vector("list", length = n) 
Sample_low_deep      <- vector("list", length = n) ; Sample_amb_deep      <- vector("list", length = n) 
Sample_low_reef      <- vector("list", length = n) ; Sample_amb_reef      <- vector("list", length = n) 
Sample_low_shallow   <- vector("list", length = n) ; Sample_amb_shallow   <- vector("list", length = n) 
# Empty vectors of output for each condition
cave_FEs_low         <- vector("list", length = n) ; cave_FEs_amb         <- vector("list", length = n)
deep_reef_FEs_low    <- vector("list", length = n) ; deep_reef_FEs_amb    <- vector("list", length = n)
reef_FEs_low         <- vector("list", length = n) ; reef_FEs_amb         <- vector("list", length = n)
shallow_reef_FEs_low <- vector("list", length = n) ; shallow_reef_FEs_amb <- vector("list", length = n)
# Datasets randomized
data_random          <- vector("list", length = n)

# Randomization process in order to get 12 quadrats per condition
for (Q in 1:n) {
  # Cave
  Sample_amb_cave[[Q]]      <- sample(Infos_FEs_split[[1]]$Quadrats, size = 12, replace = F) %>% data.frame()
  colnames(Sample_amb_cave[[Q]]) = "Quadrats" 
  cave_FEs_amb[[Q]]         <- merge(dat_cave, Sample_amb_cave[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  Sample_low_cave[[Q]]      <- sample(Infos_FEs_split[[2]]$Quadrats, size = 12, replace = F) %>% data.frame() 
  colnames(Sample_low_cave[[Q]]) = "Quadrats"
  cave_FEs_low[[Q]]         <- merge(dat_cave, Sample_low_cave[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  # Deep Reefs
  Sample_amb_deep[[Q]]      <- sample(Infos_FEs_split[[3]]$Quadrats, size = 12, replace = F) %>% data.frame()
  colnames(Sample_amb_deep[[Q]]) = "Quadrats" 
  deep_reef_FEs_amb[[Q]]    <- merge(dat_cora, Sample_amb_deep[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  Sample_low_deep[[Q]]      <- sample(Infos_FEs_split[[4]]$Quadrats, size = 12, replace = F) %>% data.frame() 
  colnames(Sample_low_deep[[Q]]) = "Quadrats"
  deep_reef_FEs_low[[Q]]    <- merge(dat_cora, Sample_low_deep[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  # Reefs
  Sample_amb_reef[[Q]]      <- sample(Infos_FEs_split[[5]]$Quadrats, size = 12, replace = F) %>% data.frame()
  colnames(Sample_amb_reef[[Q]]) = "Quadrats" 
  reef_FEs_amb[[Q]]         <- merge(dat_chia, Sample_amb_reef[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  Sample_low_reef[[Q]]      <- sample(Infos_FEs_split[[6]]$Quadrats, size = 12, replace = F) %>% data.frame() 
  colnames(Sample_low_reef[[Q]]) = "Quadrats"
  reef_FEs_low[[Q]]         <- merge(dat_chia, Sample_low_reef[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  # Shallow reefs
  Sample_amb_shallow[[Q]]   <- sample(Infos_FEs_split[[7]]$Quadrats, size = 12, replace = F) %>% data.frame()
  colnames(Sample_amb_shallow[[Q]]) = "Quadrats"
  shallow_reef_FEs_amb[[Q]] <- merge(dat_cast, Sample_amb_shallow[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  Sample_low_shallow[[Q]]   <- sample(Infos_FEs_split[[8]]$Quadrats, size = 12, replace = F) %>% data.frame() 
  colnames(Sample_low_shallow[[Q]]) = "Quadrats"
  shallow_reef_FEs_low[[Q]] <- merge(dat_cast, Sample_low_shallow[[Q]], by.x = "X", by.y = "Quadrats", all.y = T)
  # Randomized dataset
  data_random[[Q]]          <- rbind(cave_FEs_amb[[Q]], cave_FEs_low[[Q]], deep_reef_FEs_amb[[Q]], 
                                     deep_reef_FEs_low[[Q]], reef_FEs_amb[[Q]], reef_FEs_low[[Q]],
                                     shallow_reef_FEs_amb[[Q]], shallow_reef_FEs_low[[Q]]) %>% data.frame() }

# Row binding tables and going back to wide table then matrix and replacing NA with 0
quadrats_species_cover0 <- vector("list", length = n) ; quadrats_species_cover <- vector("list", length = n)
for (Q in 1:n) {
  quadrats_species_cover0[[Q]] <- data_random[[Q]] %>% pivot_wider(names_from = species, values_from = cover) 
  quadrats_species_cover[[Q]]  <- as.matrix(dplyr::select(quadrats_species_cover0[[Q]], !X))
  row.names(quadrats_species_cover[[Q]]) <- quadrats_species_cover0[[Q]]$X
  quadrats_species_cover[[Q]][is.na(quadrats_species_cover[[Q]])] = 0
  summary(apply(quadrats_species_cover[[Q]],1,sum))
  dim(quadrats_species_cover[[Q]]) # Use print() function – 96 quadrats and 232 species
  }

# Checking all species present at least once, remove those absent and reoder quadrats as in sites_quadrats
species_sumcover <- vector("list", length = n) 
for (Q in 1:n) {
  species_sumcover[[Q]] <- apply(quadrats_species_cover[[Q]], 2, sum)
  range(species_sumcover[[Q]])
  # print(length(species_sumcover[[Q]][which(species_sumcover[[Q]] == 0)])) 
  # 41–74 sp absent from all quadrats 
  quadrats_species_cover[[Q]] <- quadrats_species_cover[[Q]][, names(which(species_sumcover[[Q]]>0))]
  # print(dim(quadrats_species_cover[[Q]])) 
  # 96 quadrats and 161–195 sp
  }

# Identify Species in quadrats with no trait values and remove them
species_notraits <- vector("list", length = n) 
for (Q in 1:n) {
  species_notraits[[Q]] <- colnames(quadrats_species_cover[[Q]])[!(colnames(quadrats_species_cover[[Q]]) %in% 
                                                               sp_tr$Species)] 
  quadrats_species_cover[[Q]] <- quadrats_species_cover[[Q]][,colnames(quadrats_species_cover[[Q]]) %in% 
                                                               sp_tr$Species] 
  # print(dim(quadrats_species_cover[[Q]])) 
  # 96 quadrats and 156–186 sp
  }

# Names of species and Filtering absent species from trait dataframe
sp_nm <- vector("list", length = n) ; sp_tr_2 <- vector("list", length = n) 
for (Q in 1:n) {
  sp_nm[[Q]] <- colnames(quadrats_species_cover[[Q]]) 
  sp_tr_2[[Q]] <- sp_tr[sp_tr$Species %in% sp_nm[[Q]],]
  } ; sp_tr <- sp_tr_2 ; rm(sp_tr_2)

# Table with trait coding (Q = categorical , O = ordinal)
traits_cat <- data.frame(trait_name  = c("form", "feeding", "growth", "calcification", 
                                         "mobility", "agerepromaturity", "chem"),
                         trait_type  = c("N", "N", "O", "N", "O", "O", "O"))

# Species traits dataframe with correct coding of variables 
for (Q in 1:n) {
  row.names(sp_tr[[Q]])       <- sp_tr[[Q]]$Species
  sp_tr[[Q]]                  <- sp_tr[[Q]][,traits_cat$trait_name] 
  sp_tr[[Q]]$form             <- as.factor(gsub(as.character(sp_tr[[Q]]$form), pattern = " ", replacement = ""))
  sp_tr[[Q]]$feeding          <- as.factor(sp_tr[[Q]]$feeding)
  sp_tr[[Q]]$calcification    <- as.factor(sp_tr[[Q]]$calcification)
  sp_tr[[Q]]$growth           <- as.ordered(sp_tr[[Q]]$growth)
  sp_tr[[Q]]$mobility         <- as.ordered(sp_tr[[Q]]$mobility)
  sp_tr[[Q]]$agerepromaturity <- as.ordered(sp_tr[[Q]]$agerepromaturity)
  sp_tr[[Q]]$chem             <- as.ordered(sp_tr[[Q]]$chem) 
  }

# Clustering species into FE
sp_to_fe <- vector("list", length = n) ; fe_nm <- vector("list", length = n)
for (Q in 1:n) {
  sp_to_fe[[Q]] <- mFD::sp.to.fe(sp_tr = sp_tr[[Q]], tr_cat = traits_cat) 
  # looking at number of species per FE
  # print(summary(sp_to_fe[[Q]]$fe_nb_sp)) # From 1 to 14 (median = mostly 2)
  # Names and number of FEs
  fe_nm[[Q]] <- unique(sp_to_fe[[Q]]$fe_nm)
  # print(length(fe_nm[[Q]])) # 63–72 FE (before maximum without randomization = 74)
  }

# List of species in each FE
fe_sp <- vector("list", length = n) ; fe_tr <- vector("list", length = n)
for (Q in 1:n) {
  fe_sp[[Q]] <- vector("list", length = n)
  for (k in fe_nm[[Q]]) {
    fe_sp[[Q]][[k]] <- names(sp_to_fe[[Q]]$sp_fe[which(sp_to_fe[[Q]]$sp_fe == k)])
    }
  fe_sp[[Q]][sapply(fe_sp[[Q]], is.null)] <- NULL # Remove Null elements from the lists
  fe_tr[[Q]] <- sp_to_fe[[Q]]$fe_tr
  }

# FE trait values and cover in quadrats 
quadrats_fe_cover <- vector("list", length = n) ; sp_k <- vector("list", length = n)
for (Q in 1:n) {
quadrats_fe_cover[[Q]] <- matrix(0, nrow = nrow(quadrats_species_cover[[Q]]), ncol = length(fe_nm[[Q]]),
                            dimnames = list(row.names(quadrats_species_cover[[Q]]), fe_nm[[Q]]))
  for (k in fe_nm[[Q]]) {
    if (length(fe_sp[[Q]][[k]]) == 1) {
    quadrats_fe_cover[[Q]][,k] <- quadrats_species_cover[[Q]][, fe_sp[[Q]][[k]]]
  } else {
    quadrats_fe_cover[[Q]][,k] <- apply(quadrats_species_cover[[Q]][, fe_sp[[Q]][[k]]], 1, sum) 
    } 
    } 
  }

# Average cover of species and FE in each habitat – pH
habph_species_cover <- vector("list", length = n) ; habph_fe_cover <- vector("list", length = n)
for (k in habph) {
  quad_k <- sites_quadrats_info[which(sites_quadrats_info$habitat_ph == k), "Quadrats"] }

for (Q in 1:n) {
  habph_species_cover[[Q]]   <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_species_cover[[Q]]),
                                       dimnames = list(habph, colnames(quadrats_species_cover[[Q]])))
  habph_fe_cover[[Q]]        <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_fe_cover[[Q]]),
                                       dimnames = list(habph, colnames(quadrats_fe_cover[[Q]]))) 
  for (k in habph) {
    habph_species_cover[[Q]][k,] <- apply(quadrats_species_cover[[Q]][quad_k,], 2, mean)
    habph_fe_cover[[Q]][k,] <- apply(quadrats_fe_cover[[Q]][quad_k,], 2, mean) 
    }
  }
  
## Quick exploration and functional entities computation -----------------------------------------------------------
# SCRIPT B ---------------------------------------------------------------------------------------------------------
###### Functional indices ------------------------------------------------------------------------------------------

