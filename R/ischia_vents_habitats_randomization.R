# Functional biodiversity shifts within and across benthic habitats under ocean acidification ----------------------
## Núria Teixidó*, Enric Ballesteros, Samir Alliouane, Maria Cristina Gambi, Jean-Pierre Gattuso, Kristy 
## Kroeker, Fiorenza Micheli, Alice Mirasole, Sebastien Villéger, Cinzia De Vittor, Valeriano Parravacini 
## *corresponding author. Email: nuria.teixido@imev-mer.fr; nuria.teixido@szn.it 

## Code led by Jeremy Carlot, Sebastien Villeger, Valeriano Parravicini and Núria Teixidó

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
library(ecole)

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
sp_tr               <- read.csv(file.path(dir_raw_data,"species_traits.csv"), sep = ";")
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
n = 100 ; source(file.path(dir_scripts,"Lists_and_vectors.R"))

## Data preparation ------------------------------------------------------------------------------------------------
# SCRIPT A ---------------------------------------------------------------------------------------------------------

row.names(sites_quadrats_info) <- sites_quadrats_info$Quadrats

# New variable merging habitat type with pH level
sites_quadrats_info$habitat_ph <- paste(sites_quadrats_info$habitat, sites_quadrats_info$pH, sep="_")

# Names of these combinations
habph <- unique(sites_quadrats_info$habitat_ph) # 12 levels

# Number of quadrats per habitat_pH
(habph_nbquadrats <- sites_quadrats_info %>% group_by(habitat_ph) %>% summarize(Ntot=n())) # highly variable (12–54)

## Main raw documents ----------------------------------------------------------------------------------------------
# SCRIPT Q ---------------------------------------------------------------------------------------------------------

raw_data   <- rbind(dat_cast, dat_cora, dat_chia, dat_cave) %>% data.frame()
# Define number of species observed
raw_pst    <- raw_data %>% dplyr::filter(cover > 0)
sp_raw     <- unique(raw_pst$species) %>% data.frame() # 225 species observed
# Remove species without traits information
sp_to_keep <- sp_raw$.[sp_raw$. %in% sp_tr$Species]    # 215 species observed
# Number of observations
raw_pstobs <- raw_data %>% dplyr::filter(species %in% sp_tr$Species, cover > 0) %>% 
  rename(Quadrats = X) %>% right_join(sites_quadrats_info %>% dplyr::select(Quadrats, Description.condition))

### FOR NURIA – Only one xlsx file in data/2 – data generated
xlsx::write.xlsx(raw_pstobs, file = paste(dir_data, "raw_pstobs.xlsx", sep = "/"), row.names = F)

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

# SCRIPT P ---------------------------------------------------------------------------------------------------------
###### PERMANOVA Exploration #2 ------------------------------------------------------------------------------------
# Script to quantify species richness, FE & Fdis difference >> Figure S5

for (Q in 1:10) { # Block to 10 instead of n loops because it's taking a lot of memory use.
  # Shallow reefs
  habitat              <- c("SRA")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_SRA[[Q]] <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  habitat              <- c("SRL")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_SRL[[Q]] <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  # Cave
  habitat              <- c("CA")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_CA[[Q]]  <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  habitat              <- c("CL")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_CL[[Q]]  <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  # Reefs
  habitat              <- c("RA")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_RA[[Q]]  <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  habitat              <- c("RL")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_RL[[Q]]  <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  # Deep reefs
  habitat              <- c("CoA")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_CoA[[Q]] <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  habitat              <- c("CoL")
  quad                 <- info[info$condition %in% habitat,]$Quadrats   
  AOV_dataset_CoL[[Q]] <- quadrats_biodiv[[Q]][which (quadrats_biodiv[[Q]]$Quadrats %in% quad),]
  # General dataset
  AOV_dataset[[Q]]     <- rbind(AOV_dataset_SRA[[Q]], AOV_dataset_SRL[[Q]], AOV_dataset_CA[[Q]], 
                                AOV_dataset_CL[[Q]] , AOV_dataset_RA[[Q]] , AOV_dataset_RL[[Q]], 
                                AOV_dataset_CoA[[Q]], AOV_dataset_CoL[[Q]])
  euc_dist_sr[[Q]]     <- stats::dist(AOV_dataset[[Q]]$Nb_sp, method = "euclidean")
  pairwise_sr[[Q]]     <- ecole::permanova_pairwise(x = euc_dist_sr[[Q]], 
                                                    grp = AOV_dataset[[Q]]$condition, 
                                                    permutations = 999)
  pairwise_sr[[Q]]     <- pairwise_sr[[Q]][c(1,14,23,28),c(1,3,5)]
  euc_dist_fe[[Q]]     <- stats::dist(AOV_dataset[[Q]]$FE_richness, method = "euclidean")
  pairwise_fe[[Q]]     <- ecole::permanova_pairwise(x = euc_dist_fe[[Q]], 
                                                    grp = AOV_dataset[[Q]]$condition, 
                                                    permutations = 999)
  pairwise_fe[[Q]]     <- pairwise_fe[[Q]][c(1,14,23,28),c(1,3,5)]
  euc_dist_fdis[[Q]]   <- stats::dist(AOV_dataset[[Q]]$fdis, method = "euclidean")
  pairwise_fdis[[Q]]   <- ecole::permanova_pairwise(x = euc_dist_fdis[[Q]], 
                                                    grp = AOV_dataset[[Q]]$condition, 
                                                    permutations = 999)
  pairwise_fdis[[Q]]   <- pairwise_fdis[[Q]][c(1,14,23,28),c(1,3,5)]}

# Species richness stats
pairwise_sr            <- pairwise_sr %>% bind_rows() %>% group_by(pairs) %>% summarise_all(mean)
pairwise_sr_general    <- adonis2(AOV_dataset[[1]]$Nb_sp ~ AOV_dataset[[1]]$condition, AOV_dataset[[1]])

# Functional richness stats
pairwise_fe            <- pairwise_fe %>% bind_rows() %>% group_by(pairs) %>% summarise_all(mean)
pairwise_fe_general    <- adonis2(AOV_dataset[[1]]$FE_richness ~ AOV_dataset[[1]]$condition, AOV_dataset[[1]])

# Functional dispersion stats
pairwise_fdis          <- pairwise_fdis %>% bind_rows() %>% group_by(pairs) %>% summarise_all(mean)
pairwise_fdis_general  <- adonis2(AOV_dataset[[1]]$fdis ~ AOV_dataset[[1]]$condition, AOV_dataset[[1]])

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
for (Q in 1:n) { 
  pool_coord[[Q]] <- habph2_multidimFD[[Q]]$details$sp_faxes_coord %>% data.frame() %>% rownames_to_column("FE") %>% 
    left_join(., sp_to_fe[[Q]]$fe_tr %>% rownames_to_column("FE"), by = "FE") }
pool_coord <- pool_coord %>% bind_rows() %>% 
  group_by(form, feeding, growth, calcification, mobility, agerepromaturity, chem) %>% summarise_all(mean) 
pool_coord <- pool_coord %>% data.frame() %>% mutate(., FE = paste("fe_",seq(1, length(pool_coord$FE)), sep = "")) %>% 
  column_to_rownames("FE")

# Re-attribute the accurate FE for each iteration
for (Q in 1:n) {
  FE_within_iter[[Q]] <- merge(sp_to_fe[[Q]]$fe_tr, 
                         pool_coord[,1:7] %>% rownames_to_column("True_FE"), 
                         by = c("form", "feeding", "growth", "calcification", "mobility", "agerepromaturity", "chem"), 
                         all.x = T) %>% dplyr::select(True_FE)
  FE_within_iter[[Q]] <- FE_within_iter[[Q]]$True_FE %>% as.vector()}

# vertices of all FEs in 4D
pool_coord <- pool_coord %>% 
  dplyr::select(., -c(form, feeding, growth, calcification, mobility, agerepromaturity, chem)) %>% as.matrix()
pool_vert_nm       <- rownames(pool_coord)
# range of axes
range_faxes_coord  <- range(pool_coord[,1:4])
range_axes         <- range_faxes_coord + c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1

# indices values
habph2_fd          <- data_stat_FEs

# habph2_multidimFD informations
for (Q in 1:n) { 
  colnames(habph2_multidimFD[[Q]]$details$asb_sp_occ)    = FE_within_iter[[Q]]
  colnames(habph2_multidimFD[[Q]]$details$asb_sp_relatw) = FE_within_iter[[Q]]
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
    sp_v  <- colnames(species_avg_pst)[(which(species_avg_pst[v,] > 0.5))]   # species present in v
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

# SCRIPT E ---------------------------------------------------------------------------------------------------------
#### Making Figure SM X --------------------------------------------------------------------------------------------
# MDS plot  on beta taxonomic and functional diversity with the Hill number framework 

# PCOA functional and taxo
sites_info <- sites_quadrats_info %>% dplyr::select(Quadrats, condition, Description.condition, pH, habitat)
for (Q in 1:n) {
  pcoa_fun[[Q]] <- as.data.frame(cmdscale(quadrats_beta_hill[[Q]]$funct_q1, 4)) %>% data.frame() %>% 
    rownames_to_column(var = "Quadrats") 
  pcoa_taxo[[Q]] <- as.data.frame(cmdscale(quadrats_beta_hill[[Q]]$taxo_q1, 4)) %>% data.frame() %>% 
    rownames_to_column(var = "Quadrats")}
pcoa_fun <- pcoa_fun %>% bind_rows() %>% group_by(Quadrats) %>% mutate_all(., ~replace_na(.,0)) %>% 
  summarise_all(mean) %>% ungroup() %>% left_join(sites_info)
pcoa_taxo <- pcoa_taxo %>% bind_rows() %>% group_by(Quadrats) %>% mutate_all(., ~replace_na(.,0)) %>% 
  summarise_all(mean) %>% ungroup() %>% left_join(sites_info)

# Eigenvalues to get an overall idea
for (Q in 1:n) { 
  x.fun[[Q]] <- cmdscale(quadrats_beta_hill[[Q]]$funct_q1, eig=TRUE) 
  x.tax[[Q]] <- cmdscale(quadrats_beta_hill[[Q]]$taxo_q1,eig=TRUE)
  }
cumsum(x.fun[[1]]$eig[x.fun[[1]]$eig>=0]) / sum(x.fun[[1]]$eig[x.fun[[1]]$eig>0])       # 0.62 0.72 0.78 0.81
cumsum(x.tax[[1]]$eig[x.tax[[1]]$eig>=0]) / sum(x.tax[[1]]$eig[x.tax[[1]]$eig>0])       # 0.24 0.38 0.45 0.50

# Rename Functional PCOA dataset
pcoa_fun$condition  <- as.factor(pcoa_fun$condition)
pcoa_fun$condition  <- factor(pcoa_fun$condition, levels = c("SRA", "SRL", "CL", "CA", "RL", "RA", "CoL", "CoA"))
pcoa_fun$condition  <- recode_factor(pcoa_fun$condition, SRA = "Shallow Reef Ambient pH", 
                                     SRL = "Shallow Reef Low pH", CL = "Cave Low pH", CA = "Cave Ambient pH", 
                                     RL = "Reef Low pH", RA = "Reef Ambient pH", CoL = "Deep Reef Low pH", 
                                     CoA = "Deep Reef Ambient pH")

# Rename Taxonomic PCOA dataset
pcoa_taxo$condition <- factor(pcoa_taxo$condition, levels = c("SRA", "SRL", "CL", "CA", "RL", "RA", "CoL", "CoA"))
pcoa_taxo$condition <- recode_factor(pcoa_taxo$condition, SRA = "Shallow Reef Ambient pH", 
                                     SRL = "Shallow Reef Low pH", CL = "Cave Low pH", CA = "Cave Ambient pH", 
                                     RL = "Reef Low pH", RA = "Reef Ambient pH", CoL = "Deep Reef Low pH", 
                                     CoA = "Deep Reef Ambient pH")

# V1 and V2
fig.fun.v1v2 <- ggplot(pcoa_fun, aes(x = V1, y = V2, group = condition, fill = condition, color = condition, 
                                     shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "B) Functional diversity", x = "V1", y = "V2") + scale_x_continuous(limits = c(-0.6, 0.6)) + 
  scale_y_continuous(limits  = c(-0.4, 0.4)) + theme(legend.title = element_blank(), legend.position = "none")
fig.tax.v1v2 <- ggplot(pcoa_taxo, aes(x = V1, y = V2, group = condition, fill = condition, color = condition, 
                                      shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "A) Taxonomic diversity", x = "V1", y = "V2") + scale_x_continuous(limits = c(-0.6, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.position='none')

# SCRIPT G ---------------------------------------------------------------------------------------------------------
#### Making Figure 3 -----------------------------------------------------------------------------------------------
#### Functional diversity indices based on Hill numbers across habitats and among pH zones (Low and Ambient pH)

# setting parameters for plot
for (Q in 1:n) {
  quadrats_biodiv[[Q]]$condition <- factor(quadrats_biodiv[[Q]]$condition, 
                                      levels = c("SRA", "SRL", "CL", "CA", "RL", "RA", "CoL", "CoA"))
  quadrats_biodiv[[Q]]$condition <- recode_factor(quadrats_biodiv[[Q]]$condition, SRA = "Shallow Reef Ambient pH", 
                                             SRL ="Shallow Reef Low pH", CA = "Cave Ambient pH", CL = "Cave Low pH", 
                                             RA = "Reef Ambient pH", RL = "Reef Low pH", 
                                             CoA = "Deep Reef Ambient pH", CoL = "Deep Reef Low pH")
  }
hab_ph2        <- c("Shallow Reef Ambient pH", "Shallow Reef Low pH", "Cave Ambient pH", "Cave Low pH",
                    "Reef Ambient pH", "Reef Low pH", "Deep Reef Ambient pH", "Deep Reef Low pH")
hab_ph3        <- c("Shallow Reef\nAmbient pH", "Shallow Reef\nLow pH", "Cave\nAmbient pH", "Cave\nLow pH",
                    "Reef\nAmbient pH", "Reef\nLow pH", "Deep Reef\nAmbient pH", "Deep Reef\nLow pH")
names(vcolors) <- hab_ph2
quadrats_biodiv_avg = quadrats_biodiv %>% bind_rows() %>% group_by(Quadrats, site, condition, X) %>% 
  summarise(Total_cover = mean(Total_cover), Nb_sp = mean(Nb_sp), FE_richness = mean(FE_richness), 
            FE_shannon = mean(FE_shannon), FD_q1 = mean(FD_q1), fdis = mean(fdis), fspe = mean(fspe), 
            fori = mean(fori), fide_PC1 = mean(fide_PC1), fide_PC2 = mean(fide_PC2), fide_PC3 = mean(fide_PC3), 
            fide_PC4 = mean(fide_PC4)) %>% ungroup() %>% data.frame()

# Species richness
nbsp_plot <- ggplot(data = quadrats_biodiv_avg, aes(x = condition, y = Nb_sp, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("Species richness") + ggtitle("A) Species richness") + 
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_blank(),
        axis.text.y = element_text(face ="plain", size = 14), 
        axis.title.y = element_text(face ="bold", size = 14),
        legend.position = "none", axis.title.x = element_blank())

# FE richness
fe_plot <- ggplot(data = quadrats_biodiv_avg, aes(x = condition, y = FE_richness, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("FE richness\n(# groups)") + ggtitle("B) Functional richness") +
  scale_y_continuous(limits = c(0, 30), breaks=c(0, 10, 20, 30)) + theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_blank(),
        axis.text.y = element_text(face = "plain", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none", axis.title.x = element_blank())

# Functional dispersion
fdisp_plot <- ggplot(data = quadrats_biodiv_avg, aes(x = condition, y = fdis, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("FDis") + ggtitle("C) Functional dispersion") + theme_bw() +
  scale_y_continuous(limits= c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) + scale_x_discrete(labels = hab_ph3) +
  theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "none",
        axis.text.x = element_text(size = 14,  vjust = 0.5, face = "bold"), axis.title.x = element_blank(), 
        axis.text.y = element_text (face="plain", size=14), axis.title.y = element_text (face = "bold", size = 14))

# Assembling 
boxplot <- nbsp_plot / fe_plot / fdisp_plot 

## Data trait occupancy  -------------------------------------------------------------------------------------------
# SCRIPT O ---------------------------------------------------------------------------------------------------------

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
    } else if (nrow(df[[i]]) >= 3) {
      volume_form_low[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}
    else {volume_form_low[i] = 0 }}
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
    } else if (nrow(df[[i]]) >= 3) {
      volume_form_amb[i] = ((cxhull::cxhull(df[[i]][,2:3] %>% as.matrix()))$volume)^2}
    else {volume_form_amb[i] = 0 }}
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
                         Trait.   = c("Encrusting", "Filaments", "Massive", "Tree", "Autotrophs",
                                      "Filters", "Carnivors", "Grazers", "Detritivors", "Parasites",
                                      "Low", "Moderate", "High", "Non-calcifyers", "calcifyers",
                                      "Sessile", "Mobile", "lower_1yr", "higher_1yr", "No", "Yes"),
                         VAmb     = round(c(volume_form_amb, volume_feed_amb, NA, volume_grow_amb, 
                                            volume_calc_amb, volume_mobi_amb, volume_matu_amb, 
                                            volume_chem_amb) / VTot,3),
                         VLow     = round(c(volume_form_low, volume_feed_low[c(1,2)], NA, 
                                            volume_feed_low[3], NA, volume_feed_low[4], volume_grow_low, 
                                            volume_calc_low, volume_mobi_low, volume_matu_low, 
                                            volume_chem_low) / VTot,3))

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

# SCRIPT K ---------------------------------------------------------------------------------------------------------
#### Making Figure 4 -----------------------------------------------------------------------------------------------
# Bayesian prediction probability of the likelihood of benthic cover 

# species cover in quadrats in longer style
quad_sp_long <- vector("list")
for (Q in 1:n) {
  quad_sp_long[[Q]] <- quadrats_species_cover[[Q]] %>% as.data.frame %>% rownames_to_column("Quadrats") %>%
    pivot_longer(cols = !"Quadrats", names_to = "Species", values_to = "cover") %>%
    left_join(dplyr::select(sites_quadrats_info, Quadrats, pH, habitat)) %>%
    mutate(pH_hab = paste(substr(pH,1,3), habitat, sep = "_"), .before = Quadrats) %>% dplyr::select(- pH, -habitat)
}
quad_sp_long = quad_sp_long %>% bind_rows() %>% group_by(pH_hab, Quadrats, Species) %>% 
  summarise(cover = mean(cover)) %>% ungroup() %>% data.frame()

# trait values in long format  
sp_tr_long <- sp_tr_overall %>% mutate_if(is.integer, as.character) %>%
  pivot_longer(cols = !"Species", names_to="Trait", values_to = "Modality") 

# merging then grouping by quadrats and modality to compute total cover
pH_hab <- left_join(quad_sp_long, sp_tr_long, by = "Species") %>%
  group_by(pH_hab, Quadrats, Trait, Modality) %>% summarize(cover = sum(cover)) %>% ungroup() %>% 
  mutate(ph = recode_factor(pH_hab, "amb_cave" = "amb", "amb_deep_reef" = "amb", "amb_reef" = "amb", 
                            "amb_shallow_reef" = "amb", "low_cave" = "low", "low_deep_reef" = "low", 
                            "low_reef" = "low", "low_shallow_reef" = "low"),
         habitat = recode_factor(pH_hab, "amb_cave" = "cave", "amb_deep_reef" = "deep_reef", "amb_reef" = "reef", 
                                 "amb_shallow_reef" = "shallow_reef", "low_cave" = "cave", 
                                 "low_deep_reef" = "deep_reef", "low_reef" = "reef", 
                                 "low_shallow_reef" = "shallow_reef"))

# Looking at ambient vs low pH
info$pH[which(info$pH == "ambient1")] = "ambient" ; info$pH[which(info$pH == "ambient2")] = "ambient"

# Organizing the dataset
table(pH_hab$Trait, pH_hab$Modality) ; data_model <- vector("list", 7)
# Reproduction
Reproduction = pH_hab %>% dplyr::filter(., Trait == "agerepromaturity")
Reproduction_1 = Reproduction %>% dplyr::filter(., Modality == "1") 
Reproduction_2 = Reproduction %>% dplyr::filter(., Modality == "2")
matrice_reproduction = data.frame(quadrats = Reproduction_1$Quadrats, `1` = Reproduction_1$cover, 
                                  `2` = Reproduction_2$cover) %>% column_to_rownames('quadrats')
matrice_reproduction = list(matrice_reproduction) ; names(matrice_reproduction) = "tr" 
rm(Reproduction_1, Reproduction_2)
variables <- lapply(rownames(matrice_reproduction$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[1]] <- lapply(matrice_reproduction, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Calcification
Calcification = pH_hab %>% dplyr::filter(., Trait == "calcification")
Calcification_a = Calcification %>% dplyr::filter(., Modality == "a") 
Calcification_b = Calcification %>% dplyr::filter(., Modality == "b")
matrice_calcification = data.frame(quadrats = Calcification_a$Quadrats, a = Calcification_a$cover, 
                                   b = Calcification_b$cover) %>% column_to_rownames('quadrats')
matrice_calcification = list(matrice_calcification) ; names(matrice_calcification) = "tr" 
rm(Calcification_a, Calcification_b)
variables <- lapply(rownames(matrice_calcification$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[2]] <- lapply(matrice_calcification, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Chemistry
Chem = pH_hab %>% dplyr::filter(., Trait == "chem")
Chem_1 = Chem %>% dplyr::filter(., Modality == "1") ; Chem_2 = Chem %>% dplyr::filter(., Modality == "2")
matrice_chem = data.frame(quadrats = Chem_1$Quadrats, `1` = Chem_1$cover, `2` = Chem_2$cover) %>% 
  column_to_rownames('quadrats')
matrice_chem = list(matrice_chem) ; names(matrice_chem) = "tr" ; rm(Chem_1, Chem_2)
variables <- lapply(rownames(matrice_chem$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[3]] <- lapply(matrice_chem, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Feeding
Feeding = pH_hab %>% dplyr::filter(., Trait == "feeding")
Feeding_a = Feeding %>% dplyr::filter(., Modality == "a") ; Feeding_b = Feeding %>% dplyr::filter(., Modality == "b") 
Feeding_c = Feeding %>% dplyr::filter(., Modality == "c") ; Feeding_d = Feeding %>% dplyr::filter(., Modality == "d")
Feeding_e = Feeding %>% dplyr::filter(., Modality == "e") ; Feeding_h = Feeding %>% dplyr::filter(., Modality == "h")
matrice_feeding = data.frame(quadrats = Feeding_a$Quadrats, a = Feeding_a$cover, b = Feeding_b$cover, 
                             c = Feeding_c$cover, d = Feeding_d$cover, e = Feeding_e$cover, 
                             h = Feeding_h$cover) %>% column_to_rownames('quadrats')
matrice_feeding = list(matrice_feeding) ; names(matrice_feeding) = "tr" 
rm(Feeding_a, Feeding_b, Feeding_c, Feeding_d, Feeding_e, Feeding_h)
variables <- lapply(rownames(matrice_feeding$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[4]] <- lapply(matrice_feeding, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Form
Form = pH_hab %>% dplyr::filter(., Trait == "form")
Form_a = Form %>% dplyr::filter(., Modality == "a") ; Form_b = Form %>% dplyr::filter(., Modality == "b") 
Form_c = Form %>% dplyr::filter(., Modality == "c") ; Form_d = Form %>% dplyr::filter(., Modality == "d") 
matrice_form = data.frame(quadrats = Form_a$Quadrats, a = Form_a$cover, b = Form_b$cover, c = Form_c$cover, 
                          d = Form_d$cover) %>% column_to_rownames('quadrats')
matrice_form = list(matrice_form) ; names(matrice_form) = "tr" ; rm(Form_a, Form_b, Form_c, Form_d)
variables <- lapply(rownames(matrice_form$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[5]] <- lapply(matrice_form, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Growth
Growth = pH_hab %>% dplyr::filter(., Trait == "growth")
Growth_1 = Growth %>% dplyr::filter(., Modality == "1") ; Growth_2 = Growth %>% dplyr::filter(., Modality == "2") 
Growth_3 = Growth %>% dplyr::filter(., Modality == "3")
matrice_growth = data.frame(quadrats = Growth_1$Quadrats, `1` = Growth_1$cover, `2` = Growth_2$cover, 
                            `3` = Growth_3$cover) %>% column_to_rownames('quadrats')
matrice_growth = list(matrice_growth) ; names(matrice_growth) = "tr" ; rm(Growth_1, Growth_2, Growth_3)
variables <- lapply(rownames(matrice_growth$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[6]] <- lapply(matrice_growth, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Mobility
Mobility = pH_hab %>% dplyr::filter(., Trait == "mobility")
Mobility_1 = Mobility %>% dplyr::filter(., Modality == "1") 
Mobility_2 = Mobility %>% dplyr::filter(., Modality == "2")
matrice_mobility = data.frame(quadrats = Mobility_1$Quadrats, `1` = Mobility_1$cover, `2` = Mobility_2$cover) %>% 
  column_to_rownames('quadrats')
matrice_mobility = list(matrice_mobility) ; names(matrice_mobility) = "tr" ; rm(Mobility_1, Mobility_2)
variables <- lapply(rownames(matrice_mobility$tr), 
                    function(x) { info[info$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[7]] <- lapply(matrice_mobility, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})

# Remove additional datasets not useful anymore
rm(Calcification, Chem, Feeding, Form, Growth, Mobility, Reproduction,matrice_calcification, matrice_chem, 
   matrice_feeding, matrice_form, matrice_growth, matrice_mobility, matrice_reproduction)

# Bayesian Model
# The model lasts more than 4 hours with a 2,3 GHz Dual-Core Intel Core i5 and 16 GB memory
# mn = vector("list", 7) ; for (i in 1:7) {
# mn[[i]] <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = data_model[[i]]$tr, 
#                  family = multinomial(), control=list(adapt_delta=0.99, max_treedepth=15), 
#                  chains = 2, cores = 2, iter = 4000, warmup = 1000, backend = "cmdstanr")}

# Predictions
# The predictions take time as well
# Predicted_values = vector("list", 7) ; for (i in 1:7) {Predicted_values[[i]] = 
#                                                          predict(mn[[i]], data_model[[i]]$tr)}

data_predicted = vector("list", 7)
{data_predicted[[1]] = data.frame(Habitat = rep(data_model[[1]]$tr$habitat, 2), pH = rep(data_model[[1]]$tr$pH, 2), 
                                  Condition = c(rep("1",294), rep("2",294)),
                                  Cover = c(Predicted_values[[1]][1:294,1,1], Predicted_values[[1]][1:294,1,2]), 
                                  std.error = c(Predicted_values[[1]][1:294,2,1], Predicted_values[[1]][1:294,2,2]), 
                                  Q2.5 = c(Predicted_values[[1]][1:294,3,1], Predicted_values[[1]][1:294,3,2]), 
                                  Q97.5 = c(Predicted_values[[1]][1:294,4,1], Predicted_values[[1]][1:294,4,2]),
                                  Trait = rep("Reproduction", 294*1*2))
  data_predicted[[2]] = data.frame(Habitat = rep(data_model[[2]]$tr$habitat, 2), pH = rep(data_model[[2]]$tr$pH, 2), 
                                   Condition = c(rep("a",294), rep("b",294)),
                                   Cover = c(Predicted_values[[2]][1:294,1,1], Predicted_values[[2]][1:294,1,2]), 
                                   std.error = c(Predicted_values[[2]][1:294,2,1], Predicted_values[[2]][1:294,2,2]), 
                                   Q2.5 = c(Predicted_values[[2]][1:294,3,1], Predicted_values[[2]][1:294,3,2]), 
                                   Q97.5 = c(Predicted_values[[2]][1:294,4,1], Predicted_values[[2]][1:294,4,2]),
                                   Trait = rep("Calcification", 294*1*2))
  data_predicted[[3]] = data.frame(Habitat = rep(data_model[[3]]$tr$habitat, 2), pH = rep(data_model[[3]]$tr$pH, 2), 
                                   Condition = c(rep("1",294), rep("2",294)),
                                   Cover = c(Predicted_values[[3]][1:294,1,1], Predicted_values[[3]][1:294,1,2]), 
                                   std.error = c(Predicted_values[[3]][1:294,2,1], 
                                                 Predicted_values[[3]][1:294,2,2]), 
                                   Q2.5 = c(Predicted_values[[3]][1:294,3,1], Predicted_values[[3]][1:294,3,2]), 
                                   Q97.5 = c(Predicted_values[[3]][1:294,4,1], Predicted_values[[3]][1:294,4,2]),
                                   Trait = rep("Chem", 294*1*2))
  data_predicted[[4]] = data.frame(Habitat = rep(data_model[[4]]$tr$habitat, 6), pH = rep(data_model[[4]]$tr$pH, 6), 
                                   Condition = c(rep("a",294), rep("b",294), rep("c", 294), rep("d", 294), 
                                                 rep("e", 294), rep("h", 294)),
                                   Cover = c(Predicted_values[[4]][1:294,1,1], Predicted_values[[4]][1:294,1,2], 
                                             Predicted_values[[4]][1:294,1,3], Predicted_values[[4]][1:294,1,4], 
                                             Predicted_values[[4]][1:294,1,5], Predicted_values[[4]][1:294,1,6]), 
                                   std.error = c(Predicted_values[[4]][1:294,2,1], Predicted_values[[4]][1:294,2,2], 
                                                 Predicted_values[[4]][1:294,2,3], Predicted_values[[4]][1:294,2,4], 
                                                 Predicted_values[[4]][1:294,2,5], Predicted_values[[4]][1:294,2,6]), 
                                   Q2.5 = c(Predicted_values[[4]][1:294,3,1], Predicted_values[[4]][1:294,3,2], 
                                            Predicted_values[[4]][1:294,3,3], Predicted_values[[4]][1:294,3,4], 
                                            Predicted_values[[4]][1:294,3,5], Predicted_values[[4]][1:294,3,6]), 
                                   Q97.5 = c(Predicted_values[[4]][1:294,4,1], Predicted_values[[4]][1:294,4,2], 
                                             Predicted_values[[4]][1:294,4,3], Predicted_values[[4]][1:294,4,4], 
                                             Predicted_values[[4]][1:294,4,5], Predicted_values[[4]][1:294,4,5]),
                                   Trait = rep("Feeding", 294*6*2))
  data_predicted[[5]] = data.frame(Habitat = rep(data_model[[5]]$tr$habitat, 4), pH = rep(data_model[[5]]$tr$pH, 4), 
                                   Condition = c(rep("a",294), rep("b",294), rep("c",294), rep("d",294)),
                                   Cover = c(Predicted_values[[5]][1:294,1,1], Predicted_values[[5]][1:294,1,2], 
                                             Predicted_values[[5]][1:294,1,3], Predicted_values[[5]][1:294,1,4]), 
                                   std.error = c(Predicted_values[[5]][1:294,2,1], Predicted_values[[5]][1:294,2,2], 
                                                 Predicted_values[[5]][1:294,2,3], Predicted_values[[5]][1:294,2,3]), 
                                   Q2.5 = c(Predicted_values[[5]][1:294,3,1], Predicted_values[[5]][1:294,3,2], 
                                            Predicted_values[[5]][1:294,3,3], Predicted_values[[5]][1:294,3,4]), 
                                   Q97.5 = c(Predicted_values[[5]][1:294,4,1], Predicted_values[[5]][1:294,4,2], 
                                             Predicted_values[[5]][1:294,4,3], Predicted_values[[5]][1:294,4,3]),
                                   Trait = rep("Form", 294*4*2))
  data_predicted[[6]] = data.frame(Habitat = rep(data_model[[6]]$tr$habitat, 3), pH = rep(data_model[[6]]$tr$pH, 3), 
                                   Condition = c(rep("1",294), rep("2",294), rep("3",294)),
                                   Cover = c(Predicted_values[[6]][1:294,1,1], Predicted_values[[6]][1:294,1,2], 
                                             Predicted_values[[6]][1:294,1,3]), 
                                   std.error = c(Predicted_values[[6]][1:294,2,1], Predicted_values[[6]][1:294,2,2], 
                                                 Predicted_values[[6]][1:294,2,3]), 
                                   Q2.5 = c(Predicted_values[[6]][1:294,3,1], Predicted_values[[6]][1:294,3,2], 
                                            Predicted_values[[6]][1:294,3,3]), 
                                   Q97.5 = c(Predicted_values[[6]][1:294,4,1], Predicted_values[[6]][1:294,4,2], 
                                             Predicted_values[[6]][1:294,4,3]),
                                   Trait = rep("Growth", 294*3*2))
  data_predicted[[7]] = data.frame(Habitat = rep(data_model[[7]]$tr$habitat, 2), pH = rep(data_model[[7]]$tr$pH, 2), 
                                   Condition = c(rep("1",294), rep("2",294)),
                                   Cover = c(Predicted_values[[7]][1:294,1,1], Predicted_values[[7]][1:294,1,2]), 
                                   std.error = c(Predicted_values[[7]][1:294,2,1], Predicted_values[[7]][1:294,2,2]), 
                                   Q2.5 = c(Predicted_values[[7]][1:294,3,1], Predicted_values[[7]][1:294,3,2]), 
                                   Q97.5 = c(Predicted_values[[7]][1:294,4,1], Predicted_values[[7]][1:294,4,2]),
                                   Trait = rep("Mobility", 294*1*2))}
data_predicted = data_predicted %>% bind_rows()

# SCRIPT J ---------------------------------------------------------------------------------------------------------
# Vizualisation bayesian model
data_predicted_viz = data_predicted %>% group_by(Habitat, pH, Condition, Trait) %>% summarise_all(mean) 
data_predicted_viz$Habitat = factor(data_predicted_viz$Habitat, levels = c('shallow_reef','cave','reef','deep_reef'))

colors6  <-c("#1b9e77", "#e7298a", "#7570b3", "#d95f02", "#0c2c84", "#c4c0c2")
colors4b <-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
colors3  <-c("#1b9e77", "#7570b3", "#e7298a")
colors2  <-c("#1b9e77", "#e7298a")
colorsph <-c("#bcd7e6", "#e31a1c")

# Trait: Form
form_data_predicted_viz <- data_predicted_viz %>% dplyr::filter(Trait == "Form") 
plot1 <- ggplot(form_data_predicted_viz, aes(x = pH, y = Cover, color = Condition, group = Condition)) +
  geom_line(linetype = "dashed") + geom_point(size = 5) +
  geom_errorbar(aes(ymin = Cover - std.error, ymax = Cover + std.error), width=.2, position = position_dodge(0.05)) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "ambient") - 0.6,
                xmax = which(levels(as.factor(pH)) == "ambient") + 0.5, ymin = -Inf, ymax = Inf), 
            fill = "#deebf7", color = "NA", alpha = 0.05, inherit.aes = F) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "low") - 0.5,
                xmax = which(levels(as.factor(pH)) == "low") + 0.6, ymin = -Inf, ymax = Inf),
            fill = "#fdbb84", color = "NA", alpha = 0.05, inherit.aes = F) +
  scale_color_manual(values= colors4b, name="", labels = c("Encrusting", "Filament", "Massive", "Tree")) +
  scale_y_continuous(name = "Cover (%)", limits = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_discrete(name = "", labels = c("","")) + labs(title = "Morphological form") + theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title.x = element_blank()) +
  facet_wrap(~Habitat, nrow = 1, labeller = 
               labeller(Habitat = c("shallow_reef" = "Shallow Reef", "cave" = "Cave", "reef" = "Reef", 
                                    "deep_reef" = "Deep Reef"))) + theme(strip.text = element_text(size = 14))

# Trait: Feeding
feed_data_predicted_viz <- data_predicted_viz %>% dplyr::filter(Trait=="Feeding") 
plot2 <- ggplot(feed_data_predicted_viz, aes(x = pH, y = Cover, color = Condition, group = Condition)) +
  geom_line(linetype = "dashed") + geom_point(size = 5) +
  geom_errorbar(aes(ymin = Cover - std.error, ymax = Cover + std.error), width=.2, position = position_dodge(0.05)) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "ambient") - 0.6,
                xmax = which(levels(as.factor(pH)) == "ambient") + 0.5, ymin = -Inf, ymax = Inf), 
            fill = "#deebf7", color = "NA", alpha = 0.05, inherit.aes = F) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "low") - 0.5,
                xmax = which(levels(as.factor(pH)) == "low") + 0.6, ymin = -Inf, ymax = Inf),
            fill = "#fdbb84", color = "NA", alpha = 0.05, inherit.aes = F) +
  scale_color_manual(values = colors6, name = "", labels = c("Autotroph", "Filter feeder", "Herbivor/Grazer", 
                                                             "Carnivore", "Detritivor", "Parasite")) +
  scale_y_continuous(name = "Cover (%)", limits = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_discrete(name = "", labels = c("","")) + labs(title = "Feeding activity") + theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title.x = element_blank()) +
  facet_wrap(~Habitat, nrow = 1, labeller = 
               labeller(Habitat = c("shallow_reef" = "Shallow Reef", "cave" = "Cave", "reef" = "Reef", 
                                    "deep_reef" = "Deep Reef"))) + theme(strip.text = element_text(size = 14))

# Trait: Growth
growth_data_predicted_viz <- data_predicted_viz %>% dplyr::filter(Trait == "Growth") 
plot3 <- ggplot(growth_data_predicted_viz, aes(x = pH, y = Cover, color = Condition, group = Condition)) +
  geom_line(linetype = "dashed") + geom_point(size = 5) +
  geom_errorbar(aes(ymin = Cover - std.error, ymax = Cover + std.error), width=.2, position = position_dodge(0.05)) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "ambient") - 0.6,
                xmax = which(levels(as.factor(pH)) == "ambient") + 0.5, ymin = -Inf, ymax = Inf), 
            fill = "#deebf7", color = "NA", alpha = 0.05, inherit.aes = F) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "low") - 0.5,
                xmax = which(levels(as.factor(pH)) == "low") + 0.6, ymin = -Inf, ymax = Inf),
            fill = "#fdbb84", color = "NA", alpha = 0.05, inherit.aes = F) +
  scale_color_manual(values = colors3, name = "", labels = c("Low", "Moderate", "High")) +
  scale_y_continuous(name = "Cover (%)", limits = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_discrete(name = "", labels = c("","")) + labs(title = "Growth rate") + theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title.x = element_blank()) +
  facet_wrap(~Habitat, nrow = 1, labeller = 
               labeller(Habitat = c("shallow_reef" = "Shallow Reef", "cave" = "Cave", "reef" = "Reef", 
                                    "deep_reef" = "Deep Reef"))) + theme(strip.text = element_text(size = 14))

# Trait: Calcification
cal_data_predicted_viz <- data_predicted_viz %>% dplyr::filter(Trait == "Calcification") 
plot4 <- ggplot(cal_data_predicted_viz, aes(x = pH, y = Cover, color = Condition, group = Condition)) +
  geom_line(linetype = "dashed") + geom_point(size = 5) +
  geom_errorbar(aes(ymin = Cover - std.error, ymax = Cover + std.error), width=.2, position = position_dodge(0.05)) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "ambient") - 0.6,
                xmax = which(levels(as.factor(pH)) == "ambient") + 0.5, ymin = -Inf, ymax = Inf), 
            fill = "#deebf7", color = "NA", alpha = 0.05, inherit.aes = F) +
  geom_rect(aes(xmin = which(levels(as.factor(pH)) == "low") - 0.5,
                xmax = which(levels(as.factor(pH)) == "low") + 0.6, ymin = -Inf, ymax = Inf),
            fill = "#fdbb84", color = "NA", alpha = 0.05, inherit.aes = F) +
  scale_color_manual(values= colors2, name="", labels=c("No", "Yes"))+
  scale_y_continuous(name = "Cover (%)", limits = c(0,100), breaks = c(0, 25, 50, 75, 100)) +
  scale_x_discrete(name = "", labels= c("Ambient", "Low pH")) + labs(title = "Calcification") + theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title.x = element_blank()) +
  facet_wrap(~Habitat, nrow = 1, labeller = 
               labeller(Habitat = c("shallow_reef" = "Shallow Reef", "cave" = "Cave", "reef" = "Reef", 
                                    "deep_reef" = "Deep Reef"))) + theme(strip.text = element_text(size = 14))

# Everything combined
all4 <- (plot1/plot2/plot3/plot4)

# SCRIPT L ---------------------------------------------------------------------------------------------------------
#### Making Figure 5 -----------------------------------------------------------------------------------------------
# taxonomic and functional diversity and processes of ecosystem function across habitats and among pH zones

data_trait     <- data_predicted %>% group_by(., Habitat, pH, Condition, Trait) %>% summarise_all(mean) %>% 
  data.frame()
data_trait     <- data_trait %>% group_split(Trait)
# Complexity
Habitat        <- data_trait[[4]] %>% dplyr::filter(., Condition %in% c("c","d")) %>% 
  group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Primary production
Prim_Prod      <- data_trait[[3]] %>% dplyr::filter(., Condition == "a") %>% group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Herbivory
Herbivory      <- data_trait[[3]] %>% dplyr::filter(., Condition == "c") %>% group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Predation
Predation      <- data_trait[[3]] %>% dplyr::filter(., Condition == "d") %>% group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Filter Feeders
Filters        <- data_trait[[3]] %>% dplyr::filter(., Condition == "b") %>% group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Calcification
Calcification  <- data_trait[[1]] %>% dplyr::filter(., Condition == "b") %>% group_by(Habitat, pH, Trait) %>% 
  summarize(., Cover = sum(Cover), std.error = mean(std.error), Q2.5 = sum(Q2.5), Q97.5 = sum(Q97.5))
# Dataset functioning
data_functions <- rbind(Habitat, Prim_Prod, Herbivory, Filters, Calcification) %>% data.frame() %>% 
  mutate(., Function = c(rep("Canopy-forming species", 8), rep("Autotrophs", 8), rep("Herbivores/Grazers", 8), 
                         rep("Filter feeders", 8), rep("Calcifiers", 8))) %>% 
  dplyr::select(Function, Habitat, pH, Cover, std.error, Q2.5, Q97.5) 

# Quantify the difference between Low and ambient
data_functions       = data_functions %>% group_split(pH) 
change_dataset       = data.frame(Cover = seq(1, length(data_functions[[1]]$Q2.5), 1), Q2.5 = NA, Q97.5 = NA)
change_dataset$Cover = data_functions[[2]]$Cover - data_functions[[1]]$Cover
change_dataset$Q2.5  = data_functions[[2]]$Q2.5  - data_functions[[1]]$Q2.5
change_dataset$Q97.5 = data_functions[[2]]$Q97.5 - data_functions[[1]]$Q97.5
change_dataset$std.e = sqrt(data_functions[[2]]$std.error^2 + data_functions[[1]]$std.error^2)
Qmin_value = NA ; Qmax_value = NA ; for (i in 1:length(data_functions[[1]]$Q2.5)) {
  Qmin_value[i]      = min(change_dataset[i, c(2:3)])
  Qmax_value[i]      = max(change_dataset[i, c(2:3)])}
change_dataset$Q2.5  = Qmin_value ; change_dataset$Q97.5 = Qmax_value
Function_Change      = cbind(data_functions[[1]][,c(1:2)], change_dataset) %>% data.frame()

# Statistics from Figure 3
sub_data_change = quadrats_biodiv_avg %>% group_by(condition) %>% 
  summarise(SR = mean(Nb_sp), SR_std = std_err(Nb_sp), FER = mean(FE_richness), FER_std = std_err(FE_richness),
            FDis = mean(fdis), FDis_std = std_err(fdis)) %>% data.frame() %>% 
  mutate(Habitat = c(rep("shallow_reef",2), rep("cave",2), rep("reef",2), rep("deep_reef", 2)),
         pH = rep(c("Ambient pH", "Low pH"),4)) %>% 
  dplyr::select(Habitat, pH, SR, FER, FDis, SR_std, FER_std, FDis_std) %>% group_split(pH)

# Quantify the change for diversity indexes
change_subdataset          = data.frame(SR = rep(NA, 4), FER = rep(NA, 4), FDis = rep(NA, 4),
                                        SR_std = rep(NA, 4), FER_std = rep(NA, 4), FDis_std = rep(NA, 4))
change_subdataset$SR       = sub_data_change[[2]]$SR - sub_data_change[[1]]$SR
change_subdataset$FER      = sub_data_change[[2]]$FER - sub_data_change[[1]]$FER
change_subdataset$FDis     = sub_data_change[[2]]$FDis - sub_data_change[[1]]$FDis
change_subdataset$SR_std   = sqrt(sub_data_change[[2]]$SR_std^2 + sub_data_change[[1]]$SR_std^2)
change_subdataset$FER_std  = sqrt(sub_data_change[[2]]$FER_std^2 + sub_data_change[[1]]$FER_std^2)
change_subdataset$FDis_std = sqrt(sub_data_change[[2]]$FDis_std^2 + sub_data_change[[1]]$FDis_std^2)

# Build the statistics dataset
Cave = t((data_stat_FEs[2,c(1:2,4)] - data_stat_FEs[1,c(1:2,4)]) / data_stat_FEs[1,c(1:2,4)] * 100) %>% data.frame() %>% rownames_to_column("Index_label")
colnames(Cave) = c("Index_label", "Index")
Deep = t((data_stat_FEs[4,c(1:2,4)] - data_stat_FEs[3,c(1:2,4)]) / data_stat_FEs[3,c(1:2,4)] * 100) %>% data.frame() %>% rownames_to_column("Index_label")
colnames(Deep) = c("Index_label", "Index")
Reef = t((data_stat_FEs[6,c(1:2,4)] - data_stat_FEs[5,c(1:2,4)]) / data_stat_FEs[5,c(1:2,4)] * 100) %>% data.frame() %>% rownames_to_column("Index_label")
colnames(Reef) = c("Index_label", "Index")
Shal = t((data_stat_FEs[8,c(1:2,4)] - data_stat_FEs[7,c(1:2,4)]) / data_stat_FEs[7,c(1:2,4)] * 100) %>% data.frame() %>% rownames_to_column("Index_label")
colnames(Shal) = c("Index_label", "Index")
Stat_Change = data.frame(Habitat = rep(c('shallow_reef','cave','reef','deep_reef'), each = 3), 
                         rbind(Shal, Cave, Reef, Deep))

# Vizualisation
Stat_Change$Habitat = factor(Stat_Change$Habitat, levels = rev(c('shallow_reef','cave','reef','deep_reef')))
Function_Change$Habitat = factor(Function_Change$Habitat, levels = rev(c('shallow_reef','cave','reef','deep_reef')))
color_gradient <- colorRampPalette(c("red4", "brown3", "brown1", "skyblue2", "royalblue3"))
plot(rep(1,1000),col=color_gradient(1000),pch=19,cex=3) # Viz palette

Fig5sub1 = ggplot(Stat_Change) + geom_hline(yintercept = 0) + 
  facet_wrap(~Habitat, ncol = 4, labeller = labeller(Habitat = c("shallow_reef" = "Shallow Reef", "cave" = "Cave", 
                                                                 "reef" = "Reef", "deep_reef" = "Deep Reef"))) +
  geom_segment(aes(x = Index_label, xend = Index_label, y = 0, yend = Index, color = Index), 
               position = position_dodge(.7), size = 2, linetype = 1) +
  geom_point(aes(x = Index_label, y = Index, fill = Index), position = position_dodge(.7), 
             size = 7, shape = 21, color = "black") +
  coord_flip() + theme_bw() + labs(x = "", color = "") + 
  scale_y_continuous(name = "Change in biodiversity (%)", limits = c(-62, 62), breaks = seq(-60, 60, 30)) +
  scale_fill_gradientn(colours = color_gradient(10)) + scale_color_gradientn(colours = color_gradient(10)) +
  scale_x_discrete(labels = c("RS" = "Species richness", "FE" = "Functional entity\nrichness",
                              "fdis" = "Functional dispersion")) +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12), legend.text = element_text(size = 12), 
        axis.line.x = element_blank(), axis.ticks.x = element_line(), strip.text.x = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())    

color_gradient <- colorRampPalette(c("red4", "brown3", "brown1", "skyblue1", "skyblue2", "steelblue1", "royalblue3"))
  
Fig5sub2 = Function_Change %>%
  mutate(Function = fct_relevel(Function, "Calcifiers", "Canopy-forming species", "Herbivores/Grazers", "Filter feeders", "Autotrophs")) %>%
  ggplot() + geom_hline(yintercept = 0) + 
  facet_wrap(~Habitat, ncol = 4, labeller = labeller(Habitat = c("shallow_reef" = "", "cave" = "", 
                                                                 "reef" = "", "deep_reef" = ""))) +
  geom_segment(aes(x = Function, xend = Function, y = 0, yend = Cover, color = Cover), position = position_dodge(.7), size = 2, linetype = 1) +
  geom_point(aes(x = Function, y = Cover, fill = Cover), position = position_dodge(.7), size = 7, shape = 21, color = "black") +
  coord_flip() + theme_bw() + labs(x = "", color = "") + 
  scale_y_continuous(name = "Change in cover (%)", limits = c(-43, 43), breaks = c(-25, 0, 25)) +
  scale_fill_gradientn(colours = color_gradient(100)) + scale_color_gradientn(colours = color_gradient(100)) +
  theme(legend.position = "none", axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12), legend.text = element_text(size = 12), strip.background = element_blank(),
        axis.line.x = element_line(), axis.ticks.x = element_line(), strip.text.x = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())  

# Combine
Figure_5 = Fig5sub1 / Fig5sub2 + plot_layout(heights = c(1,3))

# SCRIPT E ---------------------------------------------------------------------------------------------------------
#### Making Supplementary Figure 6 ---------------------------------------------------------------------------------

# V3 and V4
fig.fun.v3v4 <- ggplot(pcoa_fun, aes(x = V3, y = V4, group = condition, fill = condition, color = condition, 
                                     shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "D) Functional diversity", x = "V3", y = "V4") + scale_x_continuous(limits = c(-0.6, 0.6)) + 
  scale_y_continuous(limits  = c(-0.4, 0.4)) + theme(legend.title = element_blank(), legend.position = "bottom")
fig.tax.v3v4 <- ggplot(pcoa_taxo, aes(x = V3, y = V4, group = condition, fill = condition, color = condition, 
                                      shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "C) Taxonomic diversity", x = "V3", y = "V4") + scale_x_continuous(limits = c(-0.6, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.position='none')

mds_top <- fig.tax.v1v2 + fig.tax.v3v4 
mds_bot <- fig.fun.v1v2 + fig.fun.v3v4 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
mds     <- mds_top / mds_bot + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

# SCRIPT C ---------------------------------------------------------------------------------------------------------
#### Making Supplementary Figures S7 and S8 ------------------------------------------------------------------------
# mean of taxo and functional beta within / between sites & conditions
# Correlation of beta-diversity values between taxonomic and functional facets as boxplot and scatter plot S. Methods 

# From distance matrices to dataframe
beta_df = beta_df %>% bind_rows() %>% group_by(x1, x2) %>% summarise_all(mean) %>% data.frame()

# Habitats and pH level of quadrats, recoding pH into ambient vs low
quadrats_habph <- dplyr::select(sites_quadrats_info, Quadrats, pH, habitat) %>% mutate(ph=substr(pH,1,3)) %>%
  dplyr::select(-pH)

# Left_joining first time for quadrat x1 and renaming | # same with x2
quadrats_beta_df <- left_join(beta_df, quadrats_habph,by = c("x1" = "Quadrats")) %>%
  rename(id_x1 = x1, ph_x1 = ph, hab_x1 = habitat) %>% left_join(quadrats_habph, by = c("x2" = "Quadrats")) %>% 
  rename(id_x2 = x2, ph_x2 = ph, hab_x2 = habitat)

# Adding new variable with type of pair
quadrats_beta_df <- mutate(quadrats_beta_df, pair_type = case_when(
  hab_x1 == hab_x2 & ph_x1  == ph_x2 ~ "within_hab*ph",
  hab_x1 == "cave" & hab_x2 == hab_x1 & ph_x1 != ph_x2 ~ "interpH_cave",
  hab_x1 == "reef" & hab_x2 == hab_x1 & ph_x1 != ph_x2 ~ "interpH_reef",
  hab_x1 == "shallow_reef" & hab_x2 == hab_x1 & ph_x1 != ph_x2 ~ "interpH_shallow",
  hab_x1 == "deep_reef" & hab_x2 == hab_x1 & ph_x1 != ph_x2 ~ "interpH_deep",
  ph_x1  == "low" & ph_x2 == "low" & hab_x1 != hab_x2 ~ "interhab_lowpH",
  ph_x1  == "amb" & ph_x2 == "amb" & hab_x1 != hab_x2 ~ "interhab_ambpH",
  TRUE   ~ "_interhab_interpH"))

# plot of taxo vs functional beta within each type of pair
beta_cor_plot <- ggplot(data = quadrats_beta_df, aes(x = taxo_q1, y = funct_q1) ) + geom_point(size=0.05) +
  facet_wrap(vars(pair_type))

# Boxplot for taxo and functional beta across type of pairs (Figure for SM)
beta_taxo_plot  <- ggplot(data = quadrats_beta_df, aes(x = pair_type, y = taxo_q1) ) + geom_boxplot() +
  xlab("Habitat and pH of pair of quadrats") + ylab("Taxo beta (q=1)")
beta_funct_plot <- ggplot(data = quadrats_beta_df, aes(x = pair_type, y = funct_q1) ) + geom_boxplot() +
  xlab("Habitat and pH of pair of quadrats") + ylab("Func beta (q=1)")
boxplot_beta <- beta_taxo_plot / beta_funct_plot

# SCRIPT D ---------------------------------------------------------------------------------------------------------
#### Making Supplementary Figures X --------------------------------------------------------------------------------

# Color palette
hab_ph  <- c("shallow_reef_ambient1", "shallow_reef_ambient2", "shallow_reef_low", "cave_ambient1", "cave_ambient2",
             "cave_low", "reef_ambient1", "reef_ambient2", "reef_low", "deep_reef_ambient1", "deep_reef_ambient2", 
             "deep_reef_low")
vcolors <- c("#93a1fa", "#93a1fa", "#f7d305", "#6478f5", "#6478f5", "#f5a511", "#3953f7", "#3953f7", "#f78e0c", 
             "#0219ad", "#0219ad", "#f7560c")
names(vcolors) <- hab_ph

# This figure is for pure illustration
fe_tr = fe_tr[[1]] ; quadrats_fe_cover = quadrats_fe_cover[[1]]

habph_tr_moddom <- list() ; plot_t_mod_pcover = list () ; for (t in names(fe_tr) ) {
  # Modalities of trait t
  t_mod <- levels(fe_tr[,t])
  # Matrix to store results
  t_mod_pcover <- matrix(NA, length(hab_ph) , length(t_mod), dimnames=list(hab_ph, t_mod))
  # Percentage of cover of each modality through double loop
  for (h in hab_ph) {
    # Quadrats from h
    q_h <- row.names(sites_quadrats_info)[which(sites_quadrats_info[,"habitat_ph"] == h)] 
    q_h <- q_h[q_h %in% row.names(quadrats_fe_cover)]
    for (m in all_of(t_mod)) {
      # FEs with modality m
      fe_m <- row.names(fe_tr)[which(fe_tr[,t] == m)]
      # Average cover of m in h
      if(length(fe_m) == 1) {  
        t_mod_pcover[h,m] <- mean(quadrats_fe_cover[q_h,fe_m])
      } else {
        t_mod_pcover[h,m] <- mean(apply(quadrats_fe_cover[q_h,fe_m],1,sum))}}}
  # storing
  habph_tr_moddom[[t]] <- t_mod_pcover
  # illustrating as stacked barplot after pivoting
  plot_t_mod_pcover[[t]] <- t_mod_pcover %>% as_tibble(t_mod_pcover, rownames = "habitat_pH") %>%
    pivot_longer(cols = t_mod, names_to = "modality", values_to = "cover") %>%
    ggplot(aes(fill = modality, y = cover, x = habitat_pH)) + geom_bar(position = "stack", stat = "identity") +
    theme(legend.position = "top") }

## Savings figures -------------------------------------------------------------------------------------------------

# Main Figures  
ggsave(FD_xy[[1]], filename = "Figure_2.png", path = dir_plot, device = "png", width = 6, height = 12,        # 2
       dpi = 300)              
ggsave(Fig_Trait_occupancy, filename = "Figure_3.png", path = dir_plot, device = "png", width = 20,           # 3
       height = 5, dpi = 300)   
ggsave(all4, filename = "Figure_4.png", path = dir_plot, device="png", height = 25, width = 20,               # 4
       units = "cm", dpi = 300)
ggsave(Figure_5, filename = "Figure_5.png",  path = dir_plot, device="png", height = 15, width = 20,          # 5
       units = "cm", dpi = 300)

# Supplementary figures
ggsave(FD_xy[[2]], filename = "Figure_S4.png", path = dir_plot, device = "png", width = 6,                    # S4
       height = 12)              
ggsave(mds, filename = "Figure_S6.png", path = dir_plot, device = "png", height = 35, width = 35,             # S6
       units = "cm", dpi = 300)
ggsave(beta_cor_plot, filename = "Figure_S7.png", path = dir_plot, width = 4, height = 7.5)                   # S7
ggsave(boxplot_beta, filename = "Figure_S8.png", path = dir_plot, width = 9, height = 5)                      # S8
ggsave(boxplot, filename = "Figure_S9.png", path = dir_plot, device = "png", height = 35, width = 35,         # S9
       units = "cm", dpi = 300)
for (t in names(fe_tr)) {  
  ggsave(plot_t_mod_pcover[[t]], filename = paste0("Figure_Sx_", t, "_.png"), path = dir_plot_trait,          # Sx
         width = 14, height = 5)}

## Saving outputs --------------------------------------------------------------------------------------------------

save(quadrats_species_cover, file = file.path(dir_data, "quadrats_species_cover.Rdata"))
save(quadrats_multidimFD, file = file.path(dir_data, "quadrats_multidimFD.Rdata"))
save(habph_species_cover, file = file.path(dir_data, "habph_species_cover.Rdata"))
save(quadrats_beta_hill, file = file.path(dir_data, "quadrats_beta_hill.Rdata"))
save(quadrats_fe_cover, file = file.path(dir_data, "quadrats_fe_cover.Rdata"))
save(habph_multidimFD, file = file.path(dir_data, "habph_multidimFD.Rdata"))
save(quadrats_biodiv, file = file.path(dir_data, "quadrats_biodiv.Rdata"))
save(habph_fe_cover, file = file.path(dir_data, "habph_fe_cover.Rdata"))
save(fe_4D_coord, file = file.path(dir_data, "fe_4D_coord.Rdata"))
save(data_random, file = file.path(dir_data, "data_random.Rdata"))
save(traits_cat, file = file.path(dir_data, "traits_cat.Rdata"))
save(sp_to_fe, file = file.path(dir_data, "sp_to_fe.Rdata"))
save(fe_tr, file = file.path(dir_data, "fe_tr.Rdata"))
save(fe_sp, file = file.path(dir_data, "fe_sp.Rdata"))