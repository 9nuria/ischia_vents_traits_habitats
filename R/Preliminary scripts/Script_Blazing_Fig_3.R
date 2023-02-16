rm(list=ls()) ; options(mc.cores = parallel::detectCores(), warn = - 1) ; #setwd("..")

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
dir_scripts    <- "./R"                              # folder to load scripts

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

# Graphics
theme_box <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(panel.border = element_rect(color = "black", fill = NA),
          axis.ticks = element_line(color = NA),
          axis.text = element_text(color = NA),
          panel.grid.major = element_line(color = NA)) }

# Number of iterations you desire
n = 1 ; source(file.path(dir_scripts,"Lists_and_vectors.R"))

## Data preparation ------------------------------------------------------------------------------------------------
# SCRIPT A ---------------------------------------------------------------------------------------------------------

row.names(sites_quadrats_info) <- sites_quadrats_info$Quadrats

# New variable merging habitat type with pH level
sites_quadrats_info$habitat_ph <- paste(sites_quadrats_info$habitat, sites_quadrats_info$pH, sep="_")

# Names of these combinations
habph <- unique(sites_quadrats_info$habitat_ph) # 12 levels

# Number of quadrats per habitat_pH
(habph_nbquadrats <- sites_quadrats_info %>% group_by(habitat_ph) %>% summarize(Ntot=n())) # highly variable (12–54)

## Testing the difference between ambient sites --------------------------------------------------------------------
# SCRIPT N ---------------------------------------------------------------------------------------------------------

# shallow reefs
data_permanova_cast <- vegdist(
  merge(sites_quadrats_info, dat_cast, by.x = "Quadrats", by.y = "X") %>% dplyr::filter(., pH != "low") %>% 
    dplyr::select(Quadrats, species, cover) %>% pivot_wider(., names_from = Quadrats, values_from = cover) %>% 
    column_to_rownames(., var = "species") %>% as.matrix() %>% t(), "bray") %>% 
  as.matrix() %>% melt(., varnames = c("row", "col")) %>% long_to_wide_distance()
Gr  <- sites_quadrats_info %>% dplyr::filter(habitat == "shallow_reef", pH %in% c("ambient1", "ambient2"))
(permanova_cast     <- adonis2(data_permanova_cast ~ pH, data = Gr, permutations = 999, method="bray")) # No diff.

# reefs
data_permanova_reef <- vegdist(
  merge(sites_quadrats_info, dat_chia, by.x = "Quadrats", by.y = "X") %>% dplyr::filter(., pH != "low") %>% 
    dplyr::select(Quadrats, species, cover) %>% pivot_wider(., names_from = Quadrats, values_from = cover) %>% 
    column_to_rownames(., var = "species") %>% as.matrix() %>% t(), "bray") %>% 
  as.matrix() %>% melt(., varnames = c("row", "col")) %>% long_to_wide_distance()
Gr  <- sites_quadrats_info %>% dplyr::filter(habitat == "reef", pH %in% c("ambient1", "ambient2"))
(permanova_reef     <- adonis2(data_permanova_reef ~ pH, data = Gr, permutations = 999, method="bray")) # No diff.

# deep reefs
data_permanova_deep <- vegdist(
  merge(sites_quadrats_info, dat_cora, by.x = "Quadrats", by.y = "X") %>% dplyr::filter(., pH != "low") %>% 
    dplyr::select(Quadrats, species, cover) %>% pivot_wider(., names_from = Quadrats, values_from = cover) %>% 
    column_to_rownames(., var = "species") %>% as.matrix() %>% t(), "bray") %>% 
  as.matrix() %>% melt(., varnames = c("row", "col")) %>% long_to_wide_distance()
Gr   <- sites_quadrats_info %>% dplyr::filter(habitat == "deep_reef", pH %in% c("ambient1", "ambient2"))
(permanova_deep     <- adonis2(data_permanova_deep ~ pH, data = Gr, permutations = 999, method="bray")) # No diff.

# cave – Different sampling effort in the cave according to line 81
# For obscure reason, it fails the whole script
# cave_amb1           <- sites_quadrats_info %>% dplyr::filter(habitat == "cave", pH == "ambient1")
# cave_amb2           <- sites_quadrats_info %>% dplyr::filter(habitat == "cave", pH == "ambient2")
# cave_amb1           <- sample(unique(cave_amb1$Quadrats), length(unique(cave_amb2$Quadrats)), replace = F)
# cave_amb2           <- unique(cave_amb2$Quadrats)
# Quadrats_cave       <- c(cave_amb1, cave_amb2) ; rm(cave_amb1, cave_amb2)
# data_permanova_cave <- vegdist(
#  merge(sites_quadrats_info, dat_cave, by.x = "Quadrats", by.y = "X") %>% 
#    dplyr::filter(., pH != "low", Quadrats %in% Quadrats_cave) %>% dplyr::select(Quadrats, species, cover) %>% 
#    pivot_wider(., names_from = Quadrats, values_from = cover) %>% 
#    column_to_rownames(., var = "species") %>% as.matrix() %>% t(), "bray") %>% 
#  as.matrix() %>% melt(., varnames = c("row", "col")) %>% long_to_wide_distance()
# Gr   <- sites_quadrats_info %>% dplyr::filter(Quadrats %in% Quadrats_cave)
# (permanova_cave     <- adonis2(data_permanova_cave ~ pH, data = Gr, permutations = 999, method="bray")) # No diff.

## Randomization ---------------------------------------------------------------------------------------------------
# SCRIPT M ---------------------------------------------------------------------------------------------------------

Infos_FEs       <- sites_quadrats_info %>% dplyr::select(., Quadrats, pH, habitat)
Infos_FEs$pH[which(Infos_FEs$pH == "ambient1")] = "ambient"
Infos_FEs$pH[which(Infos_FEs$pH == "ambient2")] = "ambient"
Infos_FEs_split <- Infos_FEs %>% mutate(., condition = paste(habitat, pH, sep = "_")) %>% group_split(condition)
sp_tr_overall   <- sp_tr

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

# Setting up the dataset
for (Q in 1:n) { 
  FE_tr_infos[[Q]]          <- fe_tr[[Q]] %>% rownames_to_column(var = "FE")
  FE_tr_infos[[Q]]          <- habph2_multidimFD[[Q]]$details$asb_sp_occ %>% t() %>% data.frame() %>% 
    rownames_to_column(var = "FE") %>% left_join(FE_tr_infos[[Q]])
  data_trait_occupancy[[Q]] <- habph2_multidimFD[[Q]]$details$sp_faxes_coord %>% data.frame() %>% 
    rownames_to_column(var = "FE") %>% left_join(FE_tr_infos[[Q]])} 

data_trait_occupancy <- data_trait_occupancy %>% bind_rows() %>% dplyr::select(., -FE) %>% 
  group_by(form, feeding, growth, calcification, mobility, agerepromaturity, chem) %>% 
  summarise_all(mean) %>% data.frame() %>% 
  mutate(Low_presence = (cave_low + deep_reef_low + reef_low + shallow_reef_low)/4,
         Amb_presence = (cave_amb + deep_reef_amb + reef_amb + shallow_reef_amb)/4) %>% 
  dplyr::select(-cave_amb, -cave_low, -deep_reef_amb, -deep_reef_low, -reef_amb, -reef_low, 
                -shallow_reef_amb, -shallow_reef_low) 

# Colors
col_form  <- data.frame(form             = levels(data_trait_occupancy$form), 
                        col_form         = c("#D8F3DC", "#FFCDD2", "#F2D388", "#CC9966"))
col_feed  <- data.frame(feeding          = levels(data_trait_occupancy$feeding), 
                        col_feed         = c("#D9ED92", "#FFCC99", "#9CCC65", "#CC9999", "#BCAAA4", "#F4EEFF"))
col_grow  <- data.frame(growth           = levels(data_trait_occupancy$growth), 
                        col_grow         = c("#3F51B5", "#9ADCFF", "#BDBDBD"))
col_calc  <- data.frame(calcification    = levels(data_trait_occupancy$calcification), 
                        col_calc         = c("#8DB596", "#BEDBBB"))
col_mobi  <- data.frame(mobility         = levels(data_trait_occupancy$mobility), 
                        col_mobi         = c("#F6BD60", "#FFF6BF"))
col_matu  <- data.frame(agerepromaturity = levels(data_trait_occupancy$agerepromaturity), 
                        col_matu         = c("#6886C5", "#C1EFFF"))
col_chem  <- data.frame(chem             = levels(data_trait_occupancy$chem), 
                        col_chem         = c("#FF8787", "#F8C4B4"))

# Add new variables
data_trait_occupancy <- data_trait_occupancy %>% left_join(col_form) %>% left_join(col_feed) %>% 
  left_join(col_grow) %>% left_join(col_calc) %>% left_join(col_mobi) %>% left_join(col_matu) %>% 
  left_join(col_chem)

## Start Plotting
## 1st row – Low pH Conditions
{hull_total               <- data_trait_occupancy %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy %>% dplyr::filter(Low_presence > 0.1)
  # Form
  hull_form_partial <- data_trait_occupancy_low %>% group_by(form) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(form)
  Form_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_form_partial, aes(fill = form), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_form)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_form)) +
    geom_point(aes(fill = form), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Feeding
  hull_feeding_partial <- data_trait_occupancy_low %>% group_by(feeding) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(feeding)
  feeding_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_feeding_partial, aes(fill = feeding), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_feed)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_feed)) +
    geom_point(aes(fill = feeding), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Growth
  hull_growth_partial <- data_trait_occupancy_low %>% group_by(growth) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(growth)
  growth_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_growth_partial, aes(fill = growth), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_grow)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_grow)) +
    geom_point(aes(fill = growth), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Calcification
  hull_calcification_partial <- data_trait_occupancy_low %>% group_by(calcification) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(calcification)
  calcification_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_calcification_partial, aes(fill = calcification), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_calc)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_calc)) +
    geom_point(aes(fill = calcification), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Mobility
  hull_mobility_partial <- data_trait_occupancy_low %>% group_by(mobility) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(mobility)
  mobility_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_mobility_partial, aes(fill = mobility), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_mobi)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_mobi)) +
    geom_point(aes(fill = mobility), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Age Repro Maturity
  hull_agerepromaturity_partial <- data_trait_occupancy_low %>% group_by(agerepromaturity) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(agerepromaturity)
  agerepromaturity_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_agerepromaturity_partial, aes(fill = agerepromaturity), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_matu)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_matu)) +
    geom_point(aes(fill = agerepromaturity), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")
  # Chem
  hull_chem_partial <- data_trait_occupancy_low %>% group_by(chem) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_low <- data_trait_occupancy_low %>% arrange(chem)
  Chem_Low = ggplot(data_trait_occupancy_low, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_chem_partial, aes(fill = chem), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_low$col_chem)) +
    scale_fill_manual(values = unique(data_trait_occupancy_low$col_chem)) +
    geom_point(aes(fill = chem), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    theme(legend.position = "none")}

# 2nd row – Ambient pH Conditions
{hull_total               <- data_trait_occupancy %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy %>% dplyr::filter(Amb_presence > 0.1)
  # Form
  hull_form_partial <- data_trait_occupancy_amb %>% group_by(form) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(form)
  Form_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_form_partial, aes(fill = form), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_form)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_form)) +
    geom_point(aes(fill = form), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Form") + theme(legend.position = "none")
  # Feeding
  hull_feeding_partial <- data_trait_occupancy_amb %>% group_by(feeding) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(feeding)
  feeding_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_feeding_partial, aes(fill = feeding), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_feed)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_feed)) +
    geom_point(aes(fill = feeding), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Feeding") + theme(legend.position = "none")
  # Growth
  hull_growth_partial <- data_trait_occupancy_amb %>% group_by(growth) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(growth)
  growth_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_growth_partial, aes(fill = growth), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_grow)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_grow)) +
    geom_point(aes(fill = growth), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Growth") + theme(legend.position = "none")
  # Calcification
  hull_calcification_partial <- data_trait_occupancy_amb %>% group_by(calcification) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(calcification)
  calcification_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_calcification_partial, aes(fill = calcification), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_calc)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_calc)) +
    geom_point(aes(fill = calcification), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Calcification") + theme(legend.position = "none")
  # Mobility
  hull_mobility_partial <- data_trait_occupancy_amb %>% group_by(mobility) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(mobility)
  mobility_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_mobility_partial, aes(fill = mobility), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_mobi)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_mobi)) +
    geom_point(aes(fill = mobility), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Mobility") + theme(legend.position = "none")
  # Age Repro Maturity
  hull_agerepromaturity_partial <- data_trait_occupancy_amb %>% group_by(agerepromaturity) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(agerepromaturity)
  agerepromaturity_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_agerepromaturity_partial, aes(fill = agerepromaturity), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_matu)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_matu)) +
    geom_point(aes(fill = agerepromaturity), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Age reproductive maturity") + theme(legend.position = "none")
  # Chem
  hull_chem_partial <- data_trait_occupancy_amb %>% group_by(chem) %>% slice(chull(PC1, PC2)) 
  data_trait_occupancy_amb <- data_trait_occupancy_amb %>% arrange(chem)
  Chem_Amb = ggplot(data_trait_occupancy_amb, aes(x = PC1, y = PC2)) + theme_box() +
    geom_polygon(data = hull_total, col = "black", fill = NA) + 
    geom_polygon(data = hull_chem_partial, aes(fill = chem), col = "darkgrey", alpha = 0.75) +
    scale_colour_manual(values = unique(data_trait_occupancy_amb$col_chem)) +
    scale_fill_manual(values = unique(data_trait_occupancy_amb$col_chem)) +
    geom_point(aes(fill = chem), color = "black", shape = 21, alpha = .5) +
    scale_x_continuous(name = "") + scale_y_continuous(name = "") + guides(size = "none") +
    ggtitle("Chemical defenses") + theme(legend.position = "none")}

## Statistics
VTot = cxhull::cxhull(data_trait_occupancy[,8:11] %>% as.matrix())$volume
# Low pH Conditions
{# Form
  df = data_trait_occupancy_low %>% dplyr::select(form, PC1, PC2, PC3, PC4) %>% group_split(form)
  for (i in 1:4) {
    if (nrow(df[[i]]) > 4) {volume_form_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_form_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Feeding
  df = data_trait_occupancy_low %>% dplyr::select(feeding, PC1, PC2, PC3, PC4) %>% group_split(feeding)
  for (i in 1:4) {
    if (nrow(df[[i]]) > 4) {volume_feed_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else if (nrow(df[[i]]) >= 3) {
      volume_feed_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}
    else {volume_feed_low[i] = 0 }}
  # Growth
  df = data_trait_occupancy_low %>% dplyr::select(growth, PC1, PC2, PC3, PC4) %>% group_split(growth)
  for (i in 1:3) {
    if (nrow(df[[i]]) > 4) {volume_grow_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_grow_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Calcification
  df = data_trait_occupancy_low %>% dplyr::select(calcification, PC1, PC2, PC3, PC4) %>% 
    group_split(calcification)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_calc_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_calc_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Mobility
  df = data_trait_occupancy_low %>% dplyr::select(mobility, PC1, PC2, PC3, PC4) %>% 
    group_split(mobility)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_mobi_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_mobi_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Age Repro Maturity
  df = data_trait_occupancy_low %>% dplyr::select(agerepromaturity, PC1, PC2, PC3, PC4) %>% 
    group_split(agerepromaturity)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_matu_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_matu_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Chem
  df = data_trait_occupancy_low %>% dplyr::select(chem, PC1, PC2, PC3, PC4) %>% group_split(chem)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_chem_low[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_chem_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}}
# Ambient pH Conditions
{# Form
  df = data_trait_occupancy_amb %>% dplyr::select(form, PC1, PC2, PC3, PC4) %>% group_split(form)
  for (i in 1:4) {
    if (nrow(df[[i]]) > 4) {volume_form_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_form_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Feeding
  df = data_trait_occupancy_amb %>% dplyr::select(feeding, PC1, PC2, PC3, PC4) %>% group_split(feeding)
  for (i in 1:5) {
    if (nrow(df[[i]]) > 4) {volume_feed_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else if (nrow(df[[i]]) >= 3) {
      volume_feed_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}
    else {volume_feed_amb[i] = 0 }}
  # Growth
  df = data_trait_occupancy_amb %>% dplyr::select(growth, PC1, PC2, PC3, PC4) %>% group_split(growth)
  for (i in 1:3) {
    if (nrow(df[[i]]) > 4) {volume_grow_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_grow_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Calcification
  df = data_trait_occupancy_amb %>% dplyr::select(calcification, PC1, PC2, PC3, PC4) %>% 
    group_split(calcification)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_calc_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_calc_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Mobility
  df = data_trait_occupancy_amb %>% dplyr::select(mobility, PC1, PC2, PC3, PC4) %>% 
    group_split(mobility)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_mobi_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_mobi_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Age Repro Maturity
  df = data_trait_occupancy_amb %>% dplyr::select(agerepromaturity, PC1, PC2, PC3, PC4) %>% 
    group_split(agerepromaturity)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_matu_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_matu_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}
  # Chem
  df = data_trait_occupancy_amb %>% dplyr::select(chem, PC1, PC2, PC3, PC4) %>% group_split(chem)
  for (i in 1:2) {
    if (nrow(df[[i]]) > 4) {volume_chem_amb[i] = (cxhull::cxhull(df[[i]][,2:5] %>% as.matrix()))$volume
    } else {volume_chem_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}}}

# Table volume
data_trait_occupancy <- data_trait_occupancy %>% drop_na()
data_volume = data.frame(Code     = c(unique(data_trait_occupancy$form), 
                                      unique(data_trait_occupancy$feeding), 
                                      unique(data_trait_occupancy$growth), 
                                      unique(data_trait_occupancy$calcification), 
                                      unique(data_trait_occupancy$mobility), 
                                      unique(data_trait_occupancy$agerepromaturity), 
                                      unique(data_trait_occupancy$chem)),
                         Category = c(rep("Form", 4), rep("Feeding", 6), rep("Growth", 3),
                                      rep("Calcification", 2), rep("Mobility", 2), rep("AgeRepro", 2),
                                      rep("Chem_Def", 2)),
                         Trait    = c("Encrusting", "Filaments", "Massive", "Tree", "Autotrophs",
                                      "Filters", "Grazers", "Carnivors", "Detritivors", "Parasites",
                                      "Low", "Moderate", "High", "Non-calcifyers", "calcifyers",
                                      "Sessile", "Mobile", "lower_1yr", "higher_1yr", "No", "Yes"),
                         VAmb     = round(c(volume_form_amb, volume_feed_amb, NA, volume_grow_amb, 
                                       volume_calc_amb, volume_mobi_amb, volume_matu_amb, 
                                       volume_chem_amb) / VTot,3),
                         VLow     = round(c(volume_form_low, volume_feed_low[c(1,2)], NA, 
                                       volume_feed_low[3], NA, volume_feed_low[4], volume_grow_low, 
                                       volume_calc_low, volume_mobi_low, volume_matu_low, 
                                       volume_chem_low) / VTot,3),
                         FE_Amb   = c(table(data_trait_occupancy_amb$form), 
                                      c(table(data_trait_occupancy_amb$feeding)[c(1:5)], 0),
                                      table(data_trait_occupancy_amb$growth), table(data_trait_occupancy_amb$calcification),
                                      table(data_trait_occupancy_amb$mobility), table(data_trait_occupancy_amb$agerepromaturity),
                                      table(data_trait_occupancy_amb$chem)),
                         FE_Low   = c(table(data_trait_occupancy_low$form), 
                                      c(table(data_trait_occupancy_low$feeding)[c(1:2)], 0, 
                                        table(data_trait_occupancy_low$feeding)[3], 0, table(data_trait_occupancy_low$feeding)[4]),
                                      table(data_trait_occupancy_low$growth), table(data_trait_occupancy_low$calcification),
                                      table(data_trait_occupancy_low$mobility), table(data_trait_occupancy_low$agerepromaturity),
                                      table(data_trait_occupancy_low$chem)))

# Final plot
# Low pH
Low_occupancy_trait = Form_Low + feeding_Low + growth_Low + calcification_Low + mobility_Low + 
  agerepromaturity_Low + Chem_Low + plot_layout(nrow = 1, widths = 7, heights = 1)
Low_occupancy_trait = 
  wrap_elements(grid::textGrob('Low pH',rot = 90, gp = grid::gpar(fontsize = 12, fontface = 'bold'))) +
  Low_occupancy_trait + plot_layout(nrow = 1, widths = c(1,30), heights = 1)
# Ambient pH
Amb_occupancy_trait = Form_Amb + feeding_Amb + growth_Amb + calcification_Amb + mobility_Amb + 
  agerepromaturity_Amb + Chem_Amb + plot_layout(ncol = 7, widths = 7, heights = 1)
Amb_occupancy_trait = 
  wrap_elements(grid::textGrob('Ambient pH',rot = 90, gp = grid::gpar(fontsize = 12, fontface = 'bold'))) +
  Amb_occupancy_trait + plot_layout(nrow = 1, widths = c(1,30), heights = 1)
# Plot combined
(Fig_Trait_occupancy = (Amb_occupancy_trait / Low_occupancy_trait))

data_Figure_3 <- list(
  list(data_volume),
  list(Form_Amb, feeding_Amb, growth_Amb, calcification_Amb, mobility_Amb, agerepromaturity_Amb, Chem_Amb),
  list(Form_Low, feeding_Low, growth_Low, calcification_Low, mobility_Low, agerepromaturity_Low, Chem_Low),
  list(fe_tr[[1]], fe_fspaces[[1]]),
  Fig_Trait_occupancy)

mFD::traits.faxes.cor(
  sp_tr          = fe_tr[[1]], 
  sp_faxes_coord = fe_fspaces[[1]]$details_fspaces$sp_pc_coord, 
  tr_nm          = NULL, 
  faxes_nm       = NULL,
  name_file      = NULL, 
  color_signif   = "red",
  color_non_signif = "gray80",
  plot = T)