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

## Data preparation ------------------------------------------------------------------------------------------------
# SCRIPT A ---------------------------------------------------------------------------------------------------------

row.names(sites_quadrats_info) <- sites_quadrats_info$Quadrats

# New variable merging habitat type with pH level
sites_quadrats_info$habitat_ph <- paste(sites_quadrats_info$habitat, sites_quadrats_info$pH, sep="_")

# Names of these combinations
habph <- unique(sites_quadrats_info$habitat_ph) # 12 levels

# Number of quadrats per habitat_pH
(habph_nbquadrats <- sites_quadrats_info %>% group_by(habitat_ph) %>% summarize(Ntot=n())) # highly variable (12–54)

# Row binding tables and going back to wide table then matrix and replacing NA with 0
quadrats_species_cover0 <- rbind(dat_cast, dat_cora, dat_chia, dat_cave) %>%
  pivot_wider(names_from = species, values_from = cover) 
quadrats_species_cover <- as.matrix(dplyr::select(quadrats_species_cover0, !X))
row.names(quadrats_species_cover) <- quadrats_species_cover0$X
quadrats_species_cover[is.na(quadrats_species_cover)] <- 0
summary(apply(quadrats_species_cover,1,sum))
dim(quadrats_species_cover) # 294 quadrats and 232 species

# Checking all species present at least once
species_sumcover <- apply(quadrats_species_cover, 2, sum)
range(species_sumcover)
species_sumcover[which(species_sumcover == 0)] # 7 species absent from all quadrats

# Removing those absent species and reordering quadrats as in sites_quadrats
quadrats_species_cover <- quadrats_species_cover[sites_quadrats_info$Quadrats, names(which(species_sumcover>0))]
dim(quadrats_species_cover) # 294 quadrats and 225 species

# Species in quadrats with no trait values 
(species_notraits <- colnames(quadrats_species_cover)[!(colnames(quadrats_species_cover) %in% sp_tr$Species)])

# Removing them from quadrat*species cover matrix
quadrats_species_cover <- quadrats_species_cover[,colnames(quadrats_species_cover) %in% sp_tr$Species]

# Names of species
sp_nm <- colnames(quadrats_species_cover)

# Filtering absent species from trait dataframe
sp_tr <- sp_tr[sp_tr$Species %in% sp_nm,]

# Table with trait coding (Q = categorical , O = ordinal)
traits_cat <- data.frame(trait_name  = c("form", "feeding", "growth", "calcification", 
                                         "mobility", "agerepromaturity", "chem"),
                         trait_type  = c("N", "N", "O", "N", "O", "O", "O"))

# Species traits dataframe with correct coding of variables 
row.names(sp_tr)       <- sp_tr$Species
sp_tr                  <- sp_tr[,traits_cat$trait_name] 
sp_tr$form             <- as.factor(gsub(as.character(sp_tr$form), pattern = " ", replacement = ""))
sp_tr$feeding          <- as.factor(sp_tr$feeding)
sp_tr$calcification    <- as.factor(sp_tr$calcification)
sp_tr$growth           <- as.ordered(sp_tr$growth)
sp_tr$mobility         <- as.ordered(sp_tr$mobility)
sp_tr$agerepromaturity <- as.ordered(sp_tr$agerepromaturity)
sp_tr$chem             <- as.ordered(sp_tr$chem)

# Clustering species into FE
sp_to_fe <- mFD::sp.to.fe(sp_tr = sp_tr, tr_cat = traits_cat)

# Names and number of FEs
fe_nm <- unique(sp_to_fe$fe_nm)

# List of species in each FE
fe_sp <- list() ; for (k in fe_nm) {fe_sp[[k]]<-names( sp_to_fe$sp_fe[which(sp_to_fe$sp_fe==k)])}

# Trait values of FE
fe_tr <- sp_to_fe$fe_tr

# FE trait values and cover in quadrats 
quadrats_fe_cover <- matrix(0, nrow = nrow(quadrats_species_cover), ncol=length(fe_nm),
                            dimnames = list( row.names(quadrats_species_cover), fe_nm))
for (k in fe_nm) {
  sp_k<-fe_sp[[k]]
  if (length(sp_k)==1) {
    quadrats_fe_cover[,k]<-quadrats_species_cover[, sp_k ]
  } else {
    quadrats_fe_cover[,k]<-apply( quadrats_species_cover[, sp_k ], 1, sum) } }

# Average cover of species and FE in each habitat – pH
habph_species_cover   <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_species_cover),
                                dimnames = list(habph, colnames(quadrats_species_cover)))
habph_fe_cover        <- matrix(0, nrow  = length(habph), ncol = ncol(quadrats_fe_cover),
                                dimnames = list(habph, colnames(quadrats_fe_cover)))
for (k in habph) {
  quad_k <- sites_quadrats_info[which(sites_quadrats_info$habitat_ph == k), "Quadrats"]
  habph_species_cover[k,] <- apply(quadrats_species_cover[quad_k,], 2, mean)
  habph_fe_cover[k,] <- apply(quadrats_fe_cover[quad_k,], 2, mean) }

## Quick exploration and functionnal entities computation ----------------------------------------------------------
# SCRIPT B ---------------------------------------------------------------------------------------------------------
###### Functional indices ------------------------------------------------------------------------------------------

# Total cover
quadrat_covertot <- apply(quadrats_species_cover, 1, sum)

# Number of sp per quadrat
quadrat_nbsp <- apply(quadrats_species_cover, 1, function(x) {length(which(x>0))})

# Number of FE per quadrat
quadrat_nbFE <- apply(quadrats_fe_cover, 1, function(x) { length(which(x>0))})

# Quadrats with at least 5 Fes to be able to compute FRic
quadrat_sup4FE <- names(which(quadrat_nbFE >= 5))

# computing Gower distance between FEs
fe_fe_dist <- funct.dist(fe_tr, tr_cat = traits_cat, metric = "gower")

# Building functional spaces from PCoA
fe_fspaces <- quality.fspaces(sp_dist = fe_fe_dist )

# Comparing their quality with mean absolute deviation
round(fe_fspaces$quality_fspaces,3) # => lowest mAD is or 4D

# FE coordinates in the 4D funct space
fe_4D_coord <- fe_fspaces$details_fspaces$sp_pc_coord[,1:4]

# Compute FDis and FIde for all quadrats
quadrats_multidimFD <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, asb_sp_w = quadrats_fe_cover, 
                                         ind_vect = c("fdis", "fide", "fspe", "fori"), 
                                         scaling = TRUE, details_returned = TRUE)

# Compute FRic, FDis and FIde for all habitats_pH
habph_multidimFD <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, asb_sp_w = habph_fe_cover,
                                      ind_vect = c("fric", "fdis", "fide", "fdiv"), 
                                      scaling = TRUE, details_returned = TRUE)

# Compute taxonomic Hill numbers on FEs (q=0 for number, q= 1 for Shannon
quadrats_taxo_hill <- alpha.fd.hill(asb_sp_w = quadrats_fe_cover, sp_dist = fe_fe_dist, q = c(0,1), tau = "min",
                                    details_returned = FALSE )
colnames(quadrats_taxo_hill) <- c("FE_richness", "FE_shannon")

# Compute functional Hill numbers on FEs (q=1, Shannon-like)
quadrats_funct_hill <- alpha.fd.hill(asb_sp_w = quadrats_fe_cover, sp_dist = fe_fe_dist, q = 1, tau = "mean",
                                     details_returned = FALSE)

# Merging all diversity indices with quadrats info
quadrats_biodiv <- data.frame(sites_quadrats_info, Total_cover=quadrat_covertot, Nb_sp=quadrat_nbsp,
                              quadrats_taxo_hill, quadrats_funct_hill, 
                              quadrats_multidimFD$functional_diversity_indices[,-1])

# compute taxonomic and functional beta with Hill numbers on FEs (q=1)
quadrats_beta_taxo_hill   <- beta.fd.hill(asb_sp_w = quadrats_fe_cover, sp_dist = fe_fe_dist, q = 1, tau = "min",
                                          details_returned = FALSE )
quadrats_beta_funct_hill  <- beta.fd.hill(asb_sp_w = quadrats_fe_cover, sp_dist = fe_fe_dist, q = 1, tau = "mean",
                                          details_returned = FALSE )
quadrats_beta_hill        <- list(taxo_q1= quadrats_beta_taxo_hill$q1, funct_q1 = quadrats_beta_funct_hill$q1)

# SCRIPT E ---------------------------------------------------------------------------------------------------------
###### PERMANOVA Exploration ---------------------------------------------------------------------------------------

info <- sites_quadrats_info %>% dplyr::select(Quadrats, condition, Description.condition, pH, habitat)

# Shallow reefs, functional
beta_df        <- mFD::dist.to.df(quadrats_beta_hill)
beta_df.fun    <- beta_df %>% dplyr::select(-taxo_q1)
habitat        <- c("SRA",  "SRL")
quad           <- info[info$condition %in% habitat,]$Quadrats    # select the quadrat within the habitat (i.e. 1s1)
divf           <- beta_df.fun[beta_df.fun$x1 %in% quad,]         # selection of quadrats Shallow reefs for asb.1
divf           <- divf[divf$x2 %in% quad,]                       # selection of quadrats Shallow reefs for asb.2
beta_long_fun  <- long_to_wide_distance(divf)
groupf         <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})
modf           <- betadisper(beta_long_fun, groupf, bias.adjust = T)
anova(modf)                                                      # To test if dispersions (variance) are different
permutest(modf, permutations = 999) ; plot(modf) ; boxplot(modf) # Permutation test for F

# Shallow reefs, taxonomic
beta_df.taxo   <- beta_df %>% dplyr::select(-funct_q1)
divtaxo        <- beta_df.taxo[beta_df.taxo$x1 %in% quad,]       # selection of quadrats Shallow reefs for asb.1
divtaxo        <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt         <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})
modtaxo        <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo) ; permutest(modtaxo, permutations = 999) ; plot(modtaxo) ; boxplot(modtaxo)

# Caves, functional
habitat        <- c("CA",  "CL")
quad           <- info[info$condition %in% habitat,]$Quadrats 
divf           <- beta_df.fun[beta_df.fun$x1 %in% quad,] 
divf           <- divf[divf$x2 %in% quad,]
beta_long_fun  <- long_to_wide_distance(divf)
groupf         <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})
modf           <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf) ; permutest(modf, permutations = 999, pairwise = T) ; plot(modf) ; boxplot(modf)

# Caves, taxonomic
divtaxo        <- beta_df.taxo[beta_df.taxo$x1 %in% quad,]
divtaxo        <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt         <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})
modtaxo        <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo) ; permutest(modtaxo, permutations = 999) ; plot(modtaxo) ; boxplot(modtaxo)

# Reefs, functional
habitat        <- c("RA",  "RL")
quad           <- info[info$condition %in% habitat,]$Quadrats 
divf           <- beta_df.fun[beta_df.fun$x1 %in% quad,] 
divf           <- divf[divf$x2 %in% quad,]
beta_long_fun  <- long_to_wide_distance(divf)
groupf         <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})
modf           <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf) ; permutest(modf, permutations = 999, pairwise = T) ; plot(modf) ; boxplot(modf)

# Reefs, taxonomic
divtaxo        <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] 
divtaxo        <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt         <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})
modtaxo        <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo) ; permutest(modtaxo, permutations = 999) ; plot(modtaxo) ; boxplot(modtaxo)

# Deep reefs, functional
habitat        <- c("CoA",  "CoL")
quad           <- info[info$condition %in% habitat,]$Quadrats 
divf           <- beta_df.fun[beta_df.fun$x1 %in% quad,] 
divf           <- divf[divf$x2 %in% quad,]
beta_long_fun  <- long_to_wide_distance(divf)
groupf         <- sapply(labels(beta_long_fun), function(x) {info[info$Quadrats == x,]$condition})
modf           <- betadisper(beta_long_fun, groupf, bias.adjust = TRUE)
anova(modf) ; permutest(modf, permutations = 999, pairwise = T) ; plot(modf) ; boxplot(modf)

# Deep reefs, taxonomic 
divtaxo        <- beta_df.taxo[beta_df.taxo$x1 %in% quad,] 
divtaxo        <- divtaxo[divtaxo$x2 %in% quad,]
beta_long_taxo <- long_to_wide_distance(divtaxo)
groupt         <- sapply(labels(beta_long_taxo), function(x) {info[info$Quadrats == x,]$condition})
modtaxo        <- betadisper(beta_long_taxo, groupt, bias.adjust = TRUE)
anova(modtaxo) ; permutest(modtaxo, permutations = 999) ; plot(modtaxo) ; boxplot(modtaxo)

# SCRIPT I ---------------------------------------------------------------------------------------------------------
###### Functional statistics ---------------------------------------------------------------------------------------

tab <- quadrats_fe_cover %>% as.data.frame %>% rownames_to_column("Quadrats") %>%
  left_join(dplyr::select(sites_quadrats_info, Quadrats, pH, habitat)) %>%
  mutate(pH_hab = paste(substr(pH,1,3), habitat ,sep = "_"), .after = Quadrats ) %>% dplyr::select(- pH, -habitat)

# average and standard error
pH_hab_mean <- tab %>% dplyr::select(-Quadrats) %>% group_by(pH_hab) %>%
  summarise(across(.cols = everything(), .fns = mean, .names = NULL))
pH_hab_se <- tab %>% dplyr::select(-Quadrats) %>% group_by(pH_hab) %>%
  summarise(across(.cols = everything(), .fns = std_err, .names=NULL))

## Analyze and Vizualisation ---------------------------------------------------------------------------------------
# SCRIPT F ---------------------------------------------------------------------------------------------------------
#### Making Figure 2 and Supplementary Figure 4 --------------------------------------------------------------------

# Computing mean cover of FEs in the 2 pH levels: Ambient and Low pH
habph2_fe_cover <- rbind(
  cave_low         = habph_fe_cover["cave_low",],
  cave_amb         = apply(habph_fe_cover[c("cave_ambient1", "cave_ambient2"),], 2, mean),
  deep_reef_low    = habph_fe_cover["deep_reef_low",],
  deep_reef_amb    = apply(habph_fe_cover[c("deep_reef_ambient1", "deep_reef_ambient2"),], 2, mean),
  reef_low         = habph_fe_cover["reef_low",],
  reef_amb         = apply(habph_fe_cover[c("reef_ambient1", "reef_ambient2"),], 2, mean),
  shallow_reef_low = habph_fe_cover["shallow_reef_low",],
  shallow_reef_amb = apply(habph_fe_cover[c("shallow_reef_ambient1", "shallow_reef_ambient2"),], 2, mean))

# Computing FRic, FDis and FIde for all habitats * 2 levels of pH
habph2_multidimFD  <- alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, asb_sp_w = habph2_fe_cover, ind_vect = 
                                          c("fric", "fdis", "fide", "fdiv"), scaling = T, details_returned = T)

## Plotting parameters
# Color palette
hab_ph2            <- c("shallow_reef_amb", "shallow_reef_low", "cave_amb", "cave_low", "reef_amb", "reef_low", 
                        "deep_reef_amb", "deep_reef_low")
vcolors            <- c("#93a1fa", "#f7d305", "#6478f5", "#f5a511", "#3953f7", "#f78e0c", "#0219ad", "#f7560c")
names(vcolors)     <- hab_ph2

# Coordinates of all species
pool_coord         <- habph2_multidimFD$details$sp_faxes_coord
# vertices of all FEs in 4D
pool_vert_nm       <- habph2_multidimFD$details$pool_vert_nm
# range of axes
range_faxes_coord  <- range(pool_coord[,1:4])
range_axes         <- range_faxes_coord + c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1

# indices values
habph2_fd          <- habph2_multidimFD$functional_diversity_indices

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
    sp_v  <- names(which(habph2_multidimFD$details$asb_sp_occ[v,] == 1))     # species present in v
    # background with axes range set + title
    ggplot_v <- background.plot(range_faxes = range_axes, faxes_nm = paste0("PC", xy), color_bg = "grey95")
    ggplot_v <- ggplot_v + labs(subtitle=labels[v,])
    # convex hull of species pool
    ggplot_v <- pool.plot(ggplot_bg = ggplot_v, sp_coord2D = pool_coord[,xy], vertices_nD = pool_vert_nm, 
                          plot_pool = F, color_ch = NA, fill_ch = "white", alpha_ch = 1)
    # plot convex hull of assemblage but not species
    ggplot_v <- fric.plot(ggplot_bg = ggplot_v, asb_sp_coord2D = list(vv = pool_coord[sp_v,xy]),
                          asb_vertices_nD = list(vv = habph2_multidimFD$details$asb_vert_nm[[v]]), plot_sp = F,
                          color_ch = c(vv = col_v), fill_ch = c(vv = col_v), alpha_ch = c(vv = 0.1))
    # plot species weights using plot.fide without showing mean value
    ggplot_v <- fide.plot(ggplot_bg = ggplot_v, asb_sp_coord2D = list(vv = pool_coord[sp_v,xy]),
                          asb_sp_relatw = list(vv = habph2_multidimFD$details$asb_sp_relatw[v,sp_v]),
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
#### Making Figure 3 -----------------------------------------------------------------------------------------------

# PCOA functional and taxo
pcoa_fun <- as.data.frame(cmdscale(quadrats_beta_hill$funct_q1, 4)) 
pcoa_fun <- sites_quadrats_info %>% dplyr::select(Quadrats, condition, Description.condition, pH, habitat) %>% 
  bind_cols(pcoa_fun)
pcoa_taxo <- as.data.frame(cmdscale(quadrats_beta_hill$taxo_q1, 4))
pcoa_taxo <- sites_quadrats_info %>% dplyr::select(Quadrats, condition, Description.condition, pH, habitat) %>% 
  bind_cols(pcoa_taxo)

# Eigenvalues
x.func <- cmdscale(quadrats_beta_hill$funct_q1, eig=TRUE) 
cumsum(x.func$eig[x.func$eig>=0]) / sum(x.func$eig[x.func$eig>0]) # 0.6242943 0.7218448 0.7750432 0.8124273
x.tax <- cmdscale(quadrats_beta_hill$taxo_q1,eig=TRUE)
cumsum(x.tax$eig[x.tax$eig>=0]) / sum(x.tax$eig[x.tax$eig>0])     # 0.2394259 0.3809141 0.4459513 0.5033202

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
  labs(title = "B) Functional diversity", x = "V1", y = "V2") + scale_x_continuous(limits = c(-0.55, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.title = element_blank(), legend.position = "bottom")
fig.tax.v1v2 <- ggplot(pcoa_taxo, aes(x = V1, y = V2, group = condition, fill = condition, color = condition, 
                                      shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "A) Taxonomic diversity", x = "V1", y = "V2") + scale_x_continuous(limits = c(-0.55, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.position='none')
mdsv1v2 <- (fig.tax.v1v2/fig.fun.v1v2)

# SCRIPT G ---------------------------------------------------------------------------------------------------------
#### Making Figure 4 -----------------------------------------------------------------------------------------------

# setting parameters for plot
quadrats_biodiv$condition <- factor(quadrats_biodiv$condition, 
                                    levels = c("SRA", "SRL", "CL", "CA", "RL", "RA", "CoL", "CoA"))
quadrats_biodiv$condition <- recode_factor(quadrats_biodiv$condition, SRA = "Shallow Reef Ambient pH", 
                                           SRL ="Shallow Reef Low pH", CA = "Cave Ambient pH", CL = "Cave Low pH", 
                                           RA = "Reef Ambient pH", RL = "Reef Low pH", CoA = "Deep Reef Ambient pH", 
                                           CoL = "Deep Reef Low pH")
hab_ph2        <- c("Shallow Reef Ambient pH", "Shallow Reef Low pH", "Cave Ambient pH", "Cave Low pH",
                    "Reef Ambient pH", "Reef Low pH", "Deep Reef Ambient pH", "Deep Reef Low pH")
hab_ph3        <- c("Shallow Reef\nAmbient pH", "Shallow Reef\nLow pH", "Cave\nAmbient pH", "Cave\nLow pH",
                    "Reef\nAmbient pH", "Reef\nLow pH", "Deep Reef\nAmbient pH", "Deep Reef\nLow pH")
names(vcolors) <- hab_ph2

# Species richness
nbsp_plot <- ggplot(data = quadrats_biodiv, aes(x = condition, y = Nb_sp, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("Species richness") + ggtitle("A) Species richness") + 
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_blank(),
        axis.text.y = element_text(face ="plain", size = 14), 
        axis.title.y = element_text(face ="bold", size = 14),
        legend.position = "none", axis.title.x = element_blank())

# FE richness
fe_plot <- ggplot(data = quadrats_biodiv, aes(x = condition, y = FE_richness, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("FE richness\n(# groups)") + ggtitle("B) Functional richness") +
  scale_y_continuous(limits = c(0, 30), breaks=c(0, 10, 20, 30)) + theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"), axis.text.x = element_blank(),
        axis.text.y = element_text(face = "plain", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        legend.position = "none", axis.title.x = element_blank())

# Functional dispersion
fdisp_plot <- ggplot(data = quadrats_biodiv, aes(x = condition, y = fdis, color = condition)) +
  geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.2, fill = vcolors) +
  geom_jitter(aes(color = condition), size = 2, shape = 16, alpha = 0.5, position = position_jitter(0.3)) +
  scale_color_manual(values = vcolors) + ylab("FDis") + ggtitle("C) Functional dispersion") + theme_bw() +
  scale_y_continuous(limits= c(0, 0.8), breaks=c(0, 0.2, 0.4, 0.6, 0.8)) + scale_x_discrete(labels = hab_ph3) +
  theme(plot.title = element_text(size = 14, face = "bold"), legend.position = "none",
        axis.text.x = element_text(size = 14,  vjust = 0.5, face = "bold"), axis.title.x = element_blank(), 
        axis.text.y = element_text (face="plain", size=14), axis.title.y = element_text (face = "bold", size = 14))

# Assembling 
boxplot <- nbsp_plot / fe_plot / fdisp_plot 

# SCRIPT K ---------------------------------------------------------------------------------------------------------
#### Making Figure 5 -----------------------------------------------------------------------------------------------

# species cover in quadrats in longer style
quad_sp_long <- quadrats_species_cover %>% as.data.frame %>% rownames_to_column("Quadrats") %>%
  pivot_longer(cols = !"Quadrats", names_to = "Species", values_to = "cover") %>%
  left_join(dplyr::select(sites_quadrats_info, Quadrats, pH, habitat)) %>%
  mutate(pH_hab = paste(substr(pH,1,3), habitat, sep = "_"), .before = Quadrats) %>% dplyr::select(- pH, -habitat)

# trait values in logn format  
sp_tr_long <- sp_tr %>% rownames_to_column("Species") %>% mutate_if(is.factor, as.character) %>%
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
#  mn[[i]] <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = data_model[[i]]$tr, 
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
# Vizualisation
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
  scale_color_manual(values= colors4b, name="", labels = c("Encrusting", "Filaments", "Massive", "Tree")) +
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
                                                             "Carnivor", "Detritivor", "Parasite")) +
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

# SCRIPT E ---------------------------------------------------------------------------------------------------------
#### Making Supplementary Figure 6 ---------------------------------------------------------------------------------

# V3 and V4
fig.fun.v3v4 <- ggplot(pcoa_fun, aes(x = V3, y = V4, group = condition, fill = condition, color = condition, 
                                     shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "B) Functional diversity", x = "V3", y = "V4") + scale_x_continuous(limits = c(-0.55, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.title = element_blank(), legend.position = "bottom")
fig.tax.v3v4 <- ggplot(pcoa_taxo, aes(x = V3, y = V4, group = condition, fill = condition, color = condition, 
                                      shape = condition)) + theme_light(base_size = 20) + geom_point(size = 6) +
  scale_colour_manual(values = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_fill_manual(values   = c("#93a1fa","#f7d305","#f5a511","#6478f5","#f78e0c","#3953f7","#f7560c","#0219ad")) +
  scale_shape_manual(values  = c(rep(21, 2), rep(24, 2), rep(22, 2), rep(23, 2))) + 
  labs(title = "A) Taxonomic diversity", x = "V3", y = "V4") + scale_x_continuous(limits = c(-0.55, 0.6)) + 
  scale_y_continuous(limits  = c(-0.5, 0.5)) + theme(legend.position='none')
mdsv3v4 <- (fig.tax.v3v4/fig.fun.v3v4)

# SCRIPT C ---------------------------------------------------------------------------------------------------------
#### Making Supplementary Figures S7 and S8 ------------------------------------------------------------------------

# From distance matrices to dataframe
beta_df <- mFD::dist.to.df(quadrats_beta_hill)

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

habph_tr_moddom <- list() ; plot_t_mod_pcover = list () ; for (t in names(fe_tr) ) {
  # Modalities of trait t
  t_mod <- levels(fe_tr[,t])
  # Matrix to store results
  t_mod_pcover <- matrix(NA, length(hab_ph) , length(t_mod), dimnames=list(hab_ph, t_mod))
  # Percentage of cover of each modality through double loop
  for (h in hab_ph) {
    # Quadrats from h
    q_h <- row.names(sites_quadrats_info)[which(sites_quadrats_info[,"habitat_ph"] == h)] 
    for (m in t_mod) {
      # FEs with modality m
      fe_m <- row.names(fe_tr)[which(fe_tr[,t] == m)]
      # Average cover of m in h
      if( length(fe_m) == 1) {  
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
ggsave("Figure_3.png", plot = mdsv1v2, path = dir_plot, device = "png", height = 35, width = 35,              # 3
       units = "cm", dpi = 300)
ggsave("Figure_4.png", plot = boxplot, path = dir_plot, device = "png", height = 35, width = 35,              # 4
       units = "cm", dpi = 300)
ggsave("Figure_5.png", plot = all4, path = dir_plot, device="png", height = 25, width = 20,                   # 5
       units = "cm", dpi = 300)

# Supplementary figures
ggsave(FD_xy[[2]], filename = "Figure_S4.png", path = dir_plot, device = "png", width = 6,                    # S4
       height = 12)              
ggsave(mdsv3v4, filename = "Figure_S6.png", path = dir_plot, device = "png", height = 35,                     # S6
       width = 35, units = "cm", dpi = 300)
ggsave(beta_cor_plot, filename = "Figure_S7.png", path = dir_plot, width = 4)                                 # S7
ggsave(boxplot_beta, filename = "Figure_S8.png", path = dir_plot, width = 9, height = 5)                      # S8
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
save(traits_cat, file = file.path(dir_data, "traits_cat.Rdata"))
save(sp_to_fe, file = file.path(dir_data, "sp_to_fe.Rdata"))
save(fe_tr, file = file.path(dir_data, "fe_tr.Rdata"))
save(fe_sp, file = file.path(dir_data, "fe_sp.Rdata"))