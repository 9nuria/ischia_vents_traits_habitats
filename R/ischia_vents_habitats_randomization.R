# Functional biodiversity shifts within and across benthic habitats under ocean acidification ----------------------
## Núria Teixidó*, Enric Ballesteros, Samir Alliouane, Maria Cristina Gambi, Jean-Pierre Gattuso, Kristy 
## Kroeker, Fiorenza Micheli, Alice Mirasole, Sebastien Villéger, Cinzia De Vittor, Valeriano Parravacini 
## *corresponding author. Email: nuria.teixido@imev-mer.fr; nuria.teixido@szn.it 

rm(list=ls()) ; options(mc.cores = parallel::detectCores(), warn = - 1) ; setwd("..")

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
dir_scripts    <- "./R"                               # folder to load scripts

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
n = 5 ; source(file.path(dir_scripts,"Lists_and_vectors.R"))

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
for (Q in 1:n) {
  quadrats_species_cover0[[Q]] <- data_random[[Q]] %>% pivot_wider(names_from = species, values_from = cover) 
  quadrats_species_cover[[Q]]  <- as.matrix(dplyr::select(quadrats_species_cover0[[Q]], !X))
  row.names(quadrats_species_cover[[Q]]) <- quadrats_species_cover0[[Q]]$X
  quadrats_species_cover[[Q]][is.na(quadrats_species_cover[[Q]])] = 0
  summary(apply(quadrats_species_cover[[Q]],1,sum))
  dim(quadrats_species_cover[[Q]]) # Use print() function – 96 quadrats and 232 species
  }

# Checking all species present at least once, remove those absent and reoder quadrats as in sites_quadrats
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
for (Q in 1:n) {
  species_notraits[[Q]] <- colnames(quadrats_species_cover[[Q]])[!(colnames(quadrats_species_cover[[Q]]) %in% 
                                                               sp_tr$Species)] 
  quadrats_species_cover[[Q]] <- quadrats_species_cover[[Q]][,colnames(quadrats_species_cover[[Q]]) %in% 
                                                               sp_tr$Species] 
  # print(dim(quadrats_species_cover[[Q]])) 
  # 96 quadrats and 156–186 sp
  }

# Names of species and Filtering absent species from trait dataframe
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
for (Q in 1:n) {
  sp_to_fe[[Q]] <- mFD::sp.to.fe(sp_tr = sp_tr[[Q]], tr_cat = traits_cat) 
  # looking at number of species per FE
  # print(summary(sp_to_fe[[Q]]$fe_nb_sp)) # From 1 to 14 (median = mostly 2)
  # Names and number of FEs
  fe_nm[[Q]] <- unique(sp_to_fe[[Q]]$fe_nm)
  # print(length(fe_nm[[Q]])) # 63–72 FE (before maximum without randomization = 74)
  }

# List of species in each FE
for (Q in 1:n) {
  fe_sp[[Q]] <- vector("list", length = n)
  for (k in fe_nm[[Q]]) {
    fe_sp[[Q]][[k]] <- names(sp_to_fe[[Q]]$sp_fe[which(sp_to_fe[[Q]]$sp_fe == k)])
    }
  fe_sp[[Q]][sapply(fe_sp[[Q]], is.null)] <- NULL # Remove Null elements from the lists
  fe_tr[[Q]] <- sp_to_fe[[Q]]$fe_tr
  }

# FE trait values and cover in quadrats 
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
for (Q in 1:n) {
  habph_species_cover[[Q]]   <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_species_cover[[Q]]),
                                dimnames = list(habph, colnames(quadrats_species_cover[[Q]])))
  habph_fe_cover[[Q]]        <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_fe_cover[[Q]]),
                                dimnames = list(habph, colnames(quadrats_fe_cover[[Q]]))) 
  for (k in habph) {
    quad_k <- sites_quadrats_info[which(sites_quadrats_info$habitat_ph == k), "Quadrats"]
    quad_k <- quad_k[quad_k %in% row.names(quadrats_species_cover[[Q]])]
    habph_species_cover[[Q]][k,] <- apply(quadrats_species_cover[[Q]][which(quadrats_species_cover[[Q]] %in% 
                                                                              quad_k),], 2, mean)
    habph_fe_cover[[Q]][k,] <- apply(quadrats_fe_cover[[Q]][quad_k,], 2, mean) 
    } 
  }

## Quick exploration and functional entities computation -----------------------------------------------------------
# SCRIPT B ---------------------------------------------------------------------------------------------------------
###### Functional indices ------------------------------------------------------------------------------------------

# Basic statistics
for (Q in 1:n) { 
  # Total cover
  quadrat_covertot[[Q]] <- apply(quadrats_species_cover[[Q]], 1, sum) 
  # Number of sp per quadrat
  quadrat_nbsp[[Q]]     <- apply(quadrats_species_cover[[Q]], 1, function(x) {length(which(x>0))})
  # print(summary(quadrat_nbsp[[Q]]))                 # from 2 to 27
  # Number of FE per quadrat
  quadrat_nbFE[[Q]]     <- apply(quadrats_fe_cover[[Q]], 1, function(x) {length(which(x>0))})
  # print(summary(quadrat_nbFE[[Q]]))                 # from 2 to 23
  # Quadrats with at least 5 Fes to be able to compute FRic
  quadrat_sup4FE[[Q]]   <- names(which(quadrat_nbFE[[Q]] >= 5))
  # print(length(quadrat_sup4FE[[Q]]))                # 89-95 out of the 96 quadrats 
  }

# Functional statistics
for (Q in 1:n) { 
  # computing Gower distance between FEs
  fe_fe_dist[[Q]]          <- funct.dist(fe_tr[[Q]], tr_cat = traits_cat, metric = "gower")
  # Building functional spaces from PCoA
  fe_fspaces[[Q]]          <- quality.fspaces(sp_dist = fe_fe_dist[[Q]])
  # Comparing their quality with mean absolute deviation
  # print(round(fe_fspaces[[Q]]$quality_fspaces,3))   #  lowest mad is 4D
  # FE coordinates in the 4D funct space
  fe_4D_coord[[Q]]         <- fe_fspaces[[Q]]$details_fspaces$sp_pc_coord[,1:4]
  # Compute FDis and FIde for all quadrats
  quadrats_multidimFD[[Q]] <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord[[Q]], asb_sp_w = quadrats_fe_cover[[Q]], 
                                                ind_vect = c("fdis", "fide", "fspe", "fori"), 
                                                scaling = TRUE, details_returned = FALSE)
  # Compute FRic, FDis and FIde for all habitats_pH
  habph_multidimFD[[Q]]    <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord[[Q]], asb_sp_w = habph_fe_cover[[Q]],
                                                ind_vect = c("fric", "fdis", "fide", "fdiv"), 
                                                scaling = TRUE, details_returned = FALSE)
  # Compute taxonomic Hill numbers on FEs (q=0 for number, q= 1 for Shannon
  quadrats_taxo_hill[[Q]]  <- alpha.fd.hill(asb_sp_w = quadrats_fe_cover[[Q]], sp_dist = fe_fe_dist[[Q]], 
                                            q = c(0,1), tau = "min", details_returned = FALSE)
  colnames(quadrats_taxo_hill[[Q]]) <- c("FE_richness", "FE_shannon")
  # Compute functional Hill numbers on FEs (q=1, Shannon-like)
  quadrats_funct_hill[[Q]] <- alpha.fd.hill(asb_sp_w = quadrats_fe_cover[[Q]], sp_dist = fe_fe_dist[[Q]], 
                                            q = 1, tau = "mean", details_returned = FALSE)
  # Merging all diversity indices with quadrats info
  quadrats_biodiv[[Q]]     <- data.frame(Total_cover = quadrat_covertot[[Q]], Nb_sp = quadrat_nbsp[[Q]], 
                                         quadrats_taxo_hill[[Q]], quadrats_funct_hill[[Q]], 
                                         quadrats_multidimFD[[Q]]$functional_diversity_indices[,-1]) 
  quadrats_biodiv[[Q]]     <- merge(sites_quadrats_info, quadrats_biodiv[[Q]], by = "row.names", all.y = TRUE) %>% 
    remove_rownames %>% column_to_rownames(var = "Row.names")
  # Compute taxonomic and functional beta with Hill numbers on FEs (q=1)
  quadrats_betax_hill[[Q]] <- beta.fd.hill(asb_sp_w = quadrats_fe_cover[[Q]], sp_dist = fe_fe_dist[[Q]], 
                                           q = 1, tau = "min", details_returned = FALSE)
  quadrats_befun_hill[[Q]] <- beta.fd.hill(asb_sp_w = quadrats_fe_cover[[Q]], sp_dist = fe_fe_dist[[Q]], 
                                           q = 1, tau = "mean", details_returned = FALSE )
  quadrats_beta_hill[[Q]]  <- list(taxo_q1= quadrats_betax_hill[[Q]]$q1, funct_q1 = quadrats_befun_hill[[Q]]$q1)
  } 
  
# SCRIPT E ---------------------------------------------------------------------------------------------------------
###### PERMANOVA Exploration ---------------------------------------------------------------------------------------
# Script to plot MDS on beta taxonomic and functional diversity with the Hill number framework. 
# Script to calculate multivariate homogeneity of group variances for funct and taxonomic beta-diversity 

info <- sites_quadrats_info %>% dplyr::select(Quadrats, condition, Description.condition, pH, habitat)

# Shallow reefs, functional
for (Q in 1:n) { 
beta_df[[Q]]        <- mFD::dist.to.df(quadrats_beta_hill[[Q]])
beta_df.fun[[Q]]    <- beta_df[[Q]] %>% dplyr::select(-taxo_q1)
habitat             <- c("SRA",  "SRL")
quad                <- info[info$condition %in% habitat,]$Quadrats                # select quadrat–habitat (i.e. 1s1)
divf[[Q]]           <- beta_df.fun[[Q]][beta_df.fun[[Q]]$x1 %in% quad,]           # select quadrat SR for asb.1
divf[[Q]]           <- divf[[Q]][divf[[Q]]$x2 %in% quad,]                         # select quadrat SR for asb.2
beta_long_fun[[Q]]  <- long_to_wide_distance(divf[[Q]])
groupf[[Q]]         <- sapply(labels(beta_long_fun[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modf[[Q]]           <- betadisper(beta_long_fun[[Q]], groupf[[Q]], bias.adjust = T)
anova(modf[[Q]])                                                                  # To test if Var1 != Var2
permutest(modf[[Q]], permutations = 999) # ; plot(modf[[Q]]) ; boxplot(modf[[Q]]) # Permutation test for F
}

# Shallow reefs, taxonomic
for (Q in 1:n) { 
beta_df.taxo[[Q]]   <- beta_df[[Q]] %>% dplyr::select(-funct_q1)
divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,]         # select quadrat–habitat (i.e. 1s1)
divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
beta_long_taxo[[Q]] <- long_to_wide_distance(divtaxo[[Q]])
groupt[[Q]]         <- sapply(labels(beta_long_taxo[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modtaxo[[Q]]        <- betadisper(beta_long_taxo[[Q]], groupt[[Q]], bias.adjust = TRUE)
anova(modtaxo[[Q]]) ; permutest(modtaxo[[Q]], permutations = 999) # ; plot(modtaxo[[Q]]) ; boxplot(modtaxo[[Q]])
}

# Caves, functional
for (Q in 1:n) { 
habitat             <- c("CA",  "CL")
quad                <- info[info$condition %in% habitat,]$Quadrats 
divf[[Q]]           <- beta_df.fun[[Q]][beta_df.fun[[Q]]$x1 %in% quad,] 
divf[[Q]]           <- divf[[Q]][divf[[Q]]$x2 %in% quad,]
beta_long_fun[[Q]]  <- long_to_wide_distance(divf[[Q]])
groupf[[Q]]         <- sapply(labels(beta_long_fun[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modf[[Q]]           <- betadisper(beta_long_fun[[Q]], groupf[[Q]], bias.adjust = TRUE)
anova(modf[[Q]]) ; permutest(modf[[Q]], permutations = 999, pairwise = T) # ; plot(modf[[Q]]) ; boxplot(modf[[Q]])
}

# Caves, taxonomic
for (Q in 1:n) { 
divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,]
divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
beta_long_taxo[[Q]] <- long_to_wide_distance(divtaxo[[Q]])
groupt[[Q]]         <- sapply(labels(beta_long_taxo[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modtaxo[[Q]]        <- betadisper(beta_long_taxo[[Q]], groupt[[Q]], bias.adjust = TRUE)
anova(modtaxo[[Q]]) ; permutest(modtaxo[[Q]], permutations = 999) # ; plot(modtaxo[[Q]]) ; boxplot(modtaxo[[Q]])
}

# Reefs, functional
for (Q in 1:n) { 
habitat             <- c("RA",  "RL")
quad                <- info[info$condition %in% habitat,]$Quadrats 
divf[[Q]]           <- beta_df.fun[[Q]][beta_df.fun[[Q]]$x1 %in% quad,] 
divf[[Q]]           <- divf[[Q]][divf[[Q]]$x2 %in% quad,]
beta_long_fun[[Q]]  <- long_to_wide_distance(divf[[Q]])
groupf[[Q]]         <- sapply(labels(beta_long_fun[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modf[[Q]]           <- betadisper(beta_long_fun[[Q]], groupf[[Q]], bias.adjust = TRUE)
anova(modf[[Q]]) ; permutest(modf[[Q]], permutations = 999, pairwise = T) # ; plot(modf[[Q]]) ; boxplot(modf[[Q]])
}

# Reefs, taxonomic
for (Q in 1:n) { 
divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,] 
divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
beta_long_taxo[[Q]] <- long_to_wide_distance(divtaxo[[Q]])
groupt[[Q]]         <- sapply(labels(beta_long_taxo[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modtaxo[[Q]]        <- betadisper(beta_long_taxo[[Q]], groupt[[Q]], bias.adjust = TRUE)
anova(modtaxo[[Q]]) ; permutest(modtaxo[[Q]], permutations = 999) # ; plot(modtaxo[[Q]]) ; boxplot(modtaxo[[Q]])
}

# Deep reefs, functional
for (Q in 1:n) { 
habitat             <- c("CoA",  "CoL")
quad                <- info[info$condition %in% habitat,]$Quadrats 
divf[[Q]]           <- beta_df.fun[[Q]][beta_df.fun[[Q]]$x1 %in% quad,] 
divf[[Q]]           <- divf[[Q]][divf[[Q]]$x2 %in% quad,]
beta_long_fun[[Q]]  <- long_to_wide_distance(divf[[Q]])
groupf[[Q]]         <- sapply(labels(beta_long_fun[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modf[[Q]]           <- betadisper(beta_long_fun[[Q]], groupf[[Q]], bias.adjust = TRUE)
anova(modf[[Q]]) ; permutest(modf[[Q]], permutations = 999, pairwise = T) # ; plot(modf[[Q]]) ; boxplot(modf[[Q]])
}

# Deep reefs, taxonomic 
for (Q in 1:n) { 
divtaxo[[Q]]        <- beta_df.taxo[[Q]][beta_df.taxo[[Q]]$x1 %in% quad,] 
divtaxo[[Q]]        <- divtaxo[[Q]][divtaxo[[Q]]$x2 %in% quad,]
beta_long_taxo[[Q]] <- long_to_wide_distance(divtaxo[[Q]])
groupt[[Q]]         <- sapply(labels(beta_long_taxo[[Q]]), function(x) {info[info$Quadrats == x,]$condition})
modtaxo[[Q]]        <- betadisper(beta_long_taxo[[Q]], groupt[[Q]], bias.adjust = TRUE)
anova(modtaxo[[Q]]) ; permutest(modtaxo[[Q]], permutations = 999) # ; plot(modtaxo[[Q]]) ; boxplot(modtaxo[[Q]])
}

# SCRIPT I ---------------------------------------------------------------------------------------------------------
###### Functional statistics ---------------------------------------------------------------------------------------

for (Q in 1:n) { 
  tab[[Q]]         <- quadrats_fe_cover[[Q]] %>% as.data.frame %>% rownames_to_column("Quadrats") %>%
  left_join(dplyr::select(sites_quadrats_info, Quadrats, pH, habitat)) %>%
  mutate(pH_hab = paste(substr(pH,1,3), habitat ,sep = "_"), .after = Quadrats ) %>% dplyr::select(- pH, -habitat)
  # average and standard error
  pH_hab_mean[[Q]] <- tab[[Q]] %>% dplyr::select(-Quadrats) %>% group_by(pH_hab) %>%
    summarise(across(.cols = everything(), .fns = mean, .names = NULL))
  pH_hab_se[[Q]]   <- tab[[Q]] %>% dplyr::select(-Quadrats) %>% group_by(pH_hab) %>%
    summarise(across(.cols = everything(), .fns = std_err, .names=NULL))
  }

## Functional entities for each habitat regarding  pH conditions ---------------------------------------
# SCRIPT M ---------------------------------------------------------------------------------------------------------

# Cave
for (Q in 1:n) {
  sp_tr_FEs[[Q]] = sp_tr[[Q]] %>% mutate(., species = rownames(sp_tr[[Q]])) 
  # Ambient pH
  cave_sp_FEs_Amb[[Q]] = cave_FEs_amb[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(., species)
  sp_tr_FEs_Amb[[Q]] = merge(sp_tr_FEs[[Q]], cave_sp_FEs_Amb[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames() %>% column_to_rownames(var="species") %>% na.omit()
  FE_C_amb[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Amb[[Q]], tr_cat = traits_cat)$fe_nm)
  # Low pH
  cave_sp_FEs_Low[[Q]] = cave_FEs_low[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(species)
  sp_tr_FEs_Low[[Q]] = merge(sp_tr_FEs[[Q]], cave_sp_FEs_Low[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames %>% column_to_rownames(var="species") %>% na.omit()
  FE_C_low[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Low[[Q]], tr_cat = traits_cat)$fe_nm) }

# Deep reefs
for (Q in 1:n) {
  # Ambient pH
  deep_reef_sp_FEs_Amb[[Q]] = deep_reef_FEs_amb[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(., species)
  sp_tr_FEs_Amb[[Q]] = merge(sp_tr_FEs[[Q]], deep_reef_sp_FEs_Amb[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames() %>% column_to_rownames(var="species") %>% na.omit()
  FE_DR_amb[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Amb[[Q]], tr_cat = traits_cat)$fe_nm)
  # Low pH
  deep_reef_sp_FEs_Low[[Q]] = deep_reef_FEs_low[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(species)
  sp_tr_FEs_Low[[Q]] = merge(sp_tr_FEs[[Q]], deep_reef_sp_FEs_Low[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames %>% column_to_rownames(var="species") %>% na.omit()
  FE_DR_low[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Low[[Q]], tr_cat = traits_cat)$fe_nm) }

# Reefs
for (Q in 1:n) {
  # Ambient pH
  reef_sp_FEs_Amb[[Q]] = reef_FEs_amb[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(., species) 
  sp_tr_FEs_Amb[[Q]] = merge(sp_tr_FEs[[Q]], reef_sp_FEs_Amb[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames() %>% column_to_rownames(var="species") %>% na.omit()
  FE_R_amb[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Amb[[Q]], tr_cat = traits_cat)$fe_nm)
  # Low pH
  reef_sp_FEs_Low[[Q]] = reef_FEs_low[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(species)
  sp_tr_FEs_Low[[Q]] = merge(sp_tr_FEs[[Q]], reef_sp_FEs_Low[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames %>% column_to_rownames(var="species") %>% na.omit()
  FE_R_low[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Low[[Q]], tr_cat = traits_cat)$fe_nm) }

# Shallow reefs
for (Q in 1:n) {
  # Ambient pH
  shallow_reef_sp_FEs_Amb[[Q]] = shallow_reef_FEs_amb[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(., species)
  sp_tr_FEs_Amb[[Q]] = merge(sp_tr_FEs[[Q]], shallow_reef_sp_FEs_Amb[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames() %>% column_to_rownames(var="species") %>% na.omit()
  FE_SR_amb[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Amb[[Q]], tr_cat = traits_cat)$fe_nm)
  # Low pH
  shallow_reef_sp_FEs_Low[[Q]] = shallow_reef_FEs_low[[Q]] %>% dplyr::filter(., cover > 0) %>% dplyr::select(species)
  sp_tr_FEs_Low[[Q]] = merge(sp_tr_FEs[[Q]], shallow_reef_sp_FEs_Low[[Q]], by = "species", all.y = T) %>% distinct() %>% 
    remove_rownames %>% column_to_rownames(var="species") %>% na.omit()
  FE_SR_low[Q] <- length(mFD::sp.to.fe(sp_tr = sp_tr_FEs_Low[[Q]], tr_cat = traits_cat)$fe_nm) }

# FEs summary
data_summary_FEs_std = data.frame(Habitat = rep(c("Shallow Reefs", "Cave", "Reefs", "Deep Reefs"), each = 2),
                                  pH = rep(c("Ambient", "Low"), 4),
                                  FEs = c(round(mean(FE_SR_amb), 0), round(mean(FE_SR_low), 0), 
                                          round(mean(FE_C_amb),  0), round(mean(FE_C_low),  0), 
                                          round(mean(FE_R_amb),  0), round(mean(FE_R_low),  0), 
                                          round(mean(FE_DR_amb), 0), round(mean(FE_DR_low), 0)))
rm(FE_SR_amb, FE_SR_low, FE_C_amb, FE_C_low, FE_R_amb, FE_R_low, FE_DR_amb, FE_DR_low)

## Analyze and Vizualisation ---------------------------------------------------------------------------------------
# SCRIPT F ---------------------------------------------------------------------------------------------------------

# Computing mean cover of FEs in the 2 pH levels: Ambient and Low pH
for (Q in 1:n) {
  habph2_fe_cover[[Q]]    <- rbind(
    cave_amb         = apply(habph_fe_cover[[Q]][c("cave_ambient1", "cave_ambient2"),], 2, mean),
    cave_low         = habph_fe_cover[[Q]]["cave_low",],
    deep_reef_amb    = apply(habph_fe_cover[[Q]][c("deep_reef_ambient1", "deep_reef_ambient2"),], 2, mean),
    deep_reef_low    = habph_fe_cover[[Q]]["deep_reef_low",],
    reef_amb         = apply(habph_fe_cover[[Q]][c("reef_ambient1", "reef_ambient2"),], 2, mean),
    reef_low         = habph_fe_cover[[Q]]["reef_low",],
    shallow_reef_amb = apply(habph_fe_cover[[Q]][c("shallow_reef_ambient1", "shallow_reef_ambient2"),], 2, mean),
    shallow_reef_low = habph_fe_cover[[Q]]["shallow_reef_low",]) 
  # Computing FRic, FDis and FIde for all habitats * 2 levels of pH
  habph2_multidimFD[[Q]]  <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord[[Q]], asb_sp_w = habph2_fe_cover[[Q]], 
                                               ind_vect = c("fric", "fdis", "fide", "fdiv"), scaling = T, 
                                               details_returned = T)
  # Computing Specific richness
  RS[[Q]] <- merge(data_random[[Q]], Infos_FEs, by.x = "X", by.y = "Quadrats") %>% dplyr::filter(cover > 0) %>% 
    group_by(pH, habitat) %>% summarise(RS = round(length(unique(species)))) %>% arrange(habitat) %>% ungroup() %>% 
    select(RS) %>% mutate(Zone = seq(1,8,1))
  # Summary
  habph2_label <- rownames(habph2_fe_cover[[Q]])
  data_stat_FEs[[Q]] <- habph2_multidimFD[[Q]]$functional_diversity_indices}

RS = round(bind_rows(RS) %>% group_by(Zone) %>% summarise(RS = mean(RS))) %>% dplyr::select(-Zone)
data_stat_FEs  <- data_stat_FEs %>% bind_rows() %>% mutate(., Zone = rep(seq(1,8,1), n)) %>% 
  group_by(Zone) %>% summarise_all(mean) %>% mutate(., FE = round(sp_richn)) %>% 
  mutate(RS = RS$RS) %>% data.frame() %>% mutate(row = habph2_label) %>% column_to_rownames(var = "row") %>% 
  dplyr::select(RS, FE, fric, fdis, fdiv, fide_PC1, fide_PC2, fide_PC3, fide_PC4)

#### Making Figure 2 and Supplementary Figure 4 --------------------------------------------------------------------

## Plotting parameters
# Color palette
hab_ph2            <- c("shallow_reef_amb", "shallow_reef_low", "cave_amb", "cave_low", "reef_amb", "reef_low", 
                        "deep_reef_amb", "deep_reef_low")
vcolors            <- c("#93a1fa", "#f7d305", "#6478f5", "#f5a511", "#3953f7", "#f78e0c", "#0219ad", "#f7560c")
names(vcolors)     <- hab_ph2

# Coordinates of all species
for (Q in 1:n) { pool_coord[[Q]] <- habph2_multidimFD[[Q]]$details$sp_faxes_coord %>% data.frame() %>% 
  rownames_to_column(var = "FE") } ; pool_coord <- pool_coord %>% bind_rows() %>% group_by(FE) %>% 
  summarise_all(mean) %>% column_to_rownames(var = "FE") %>% as.matrix()

# vertices of all FEs in 4D
pool_vert_nm       <- rownames(pool_coord)
# range of axes
range_faxes_coord  <- range(pool_coord[,1:4])
range_axes         <- range_faxes_coord + c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1

# indices values
habph2_fd          <- data_stat_FEs

# habph2_multidimFD informations
for (Q in 1:n) { 
  species_avg_pst[[Q]] <- habph2_multidimFD[[Q]]$details$asb_sp_occ %>% data.frame() %>% 
    rownames_to_column(var = "FE") 
  asb_sp_relatw[[Q]] <- habph2_multidimFD[[Q]]$details$asb_sp_relatw %>% data.frame() %>% 
    rownames_to_column(var = "FE") } 
species_avg_pst <- species_avg_pst %>% bind_rows() %>% group_by(FE) %>% mutate_all(., ~replace_na(.,0)) %>% 
  summarise_all(mean) %>% column_to_rownames(var = "FE") 
asb_sp_relatw <- asb_sp_relatw %>% bind_rows() %>% group_by(FE) %>% mutate_all(., ~replace_na(.,0)) %>% 
  summarise_all(mean) %>% column_to_rownames(var = "FE") %>% as.matrix()

# Building the figure
pairs_axes         <- list(c(1,2), c(3,4)) ; FD_xy = list()
labels = c("Shallow Reef Ambient pH", "Shallow Reef Low pH", "Cave Ambient pH", "Cave Low pH",
           "Reef Ambient pH", "Reef Low pH", "Deep Reef Ambient pH", "Deep Reef Low pH") %>% data.frame()
rownames(labels) = hab_ph2
for (z in 1:length(pairs_axes)) {
  xy <- pairs_axes[[z]]                                                      # names of axes 
  ggplot_list <- list()                                                      # list to store ggplot
  for (v in hab_ph2) {
    col_v <- as.character(vcolors[v])                                        # color for habitat*pH levels
    sp_v  <- colnames(species_avg_pst)[(which(species_avg_pst[v,] >= 0.5))]  # species present in v
    # background with axes range set + title
    ggplot_v <- background.plot(range_faxes = range_axes, faxes_nm = paste0("PC", xy), color_bg = "grey95")
    ggplot_v <- ggplot_v + labs(subtitle=labels[v,])
    # convex hull of species pool
    ggplot_v <- pool.plot(ggplot_bg = ggplot_v, sp_coord2D = pool_coord[,xy], vertices_nD = pool_vert_nm, 
                          plot_pool = F, color_ch = NA, fill_ch = "white", alpha_ch = 1)
    # plot convex hull of assemblage but not species
    ggplot_v <- fric.plot(ggplot_bg = ggplot_v, asb_sp_coord2D = list(vv = pool_coord[sp_v,xy]),
                          asb_vertices_nD = list(vv = sp_v), plot_sp = F,
                          color_ch = c(vv = col_v), fill_ch = c(vv = col_v), alpha_ch = c(vv = 0.1))
    # plot species weights using plot.fide without showing mean value
    ggplot_v <- fide.plot(ggplot_bg = ggplot_v, asb_sp_coord2D = list(vv = pool_coord[sp_v,xy]),
                          asb_sp_relatw = list(vv = asb_sp_relatw[v,sp_v]),
                          asb_fide_coord2D = list(vv = habph2_fd[v, paste0("fide_PC", xy)]),
                          plot_sp = T, shape_sp = c(vv = 21), color_sp = c(vv = col_v), 
                          fill_sp = c(vv =paste0(col_v, "70")), shape_fide = c(vv = 23), size_fide = c(vv = 1),
                          color_fide = c(vv = col_v), fill_fide = c(vv = col_v), color_segment = c(vv = col_v),
                          width_segment = c(vv = 0.5), linetype_segment = c(vv = 1))
    ggplot_list[[v]] <- ggplot_v}                                             # ggplot_v storing in list
  
  # patchwork of plots : 4 habitats (rows) * 2 columns (pH)
  FD_xy[[z]] <- (ggplot_list[[1]] + ggplot_list[[2]]) / (ggplot_list[[3]] + ggplot_list[[4]]) / 
    (ggplot_list[[5]] + ggplot_list[[6]]) / (ggplot_list[[7]] + ggplot_list[[8]])}
