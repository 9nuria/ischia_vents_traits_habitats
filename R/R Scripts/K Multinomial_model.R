## Setting parameters and loading packages ----
rm(list = ls()) ; options(mc.cores = parallel::detectCores()); library(Matrix) ; library(tidyverse) ; library(brms) ; library(patchwork)

# Annotations will start at the line 39
## Loading datasets ----
# folder with ready to use data 
dir_data<-"./data"

# loading dataset built in script A 
load(file.path(dir_data, "sites_quadrats_info.Rdata"))
load(file=file.path(dir_data,"quadrats_species_cover.Rdata") )
load(file=file.path(dir_data,"sp_tr.Rdata") )
infos <- read.csv("data/from_nuria/Sites_Quadrats.csv")
infos$pH[which(infos$pH == "ambient1")] = "ambient" ; infos$pH[which(infos$pH == "ambient2")] = "ambient"

## Cleaning the data the way we want ----
# species cover in quadrats in longer style
quad_sp_long <- quadrats_species_cover %>% as.data.frame %>% rownames_to_column("Quadrats") %>%
  pivot_longer(cols = !"Quadrats", names_to="Species", values_to = "cover") %>%
  left_join( select(sites_quadrats_info, Quadrats, pH, habitat )  ) %>%
  mutate(pH_hab = paste( substr(pH,1,3), habitat ,sep="_"), .before = Quadrats ) %>%
  select(- pH, - habitat)

# trait values in logn format  
sp_tr_long <- sp_tr %>% rownames_to_column("Species") %>% mutate_if(is.factor, as.character) %>%
  pivot_longer( cols = !"Species", names_to="Trait", values_to = "Modality" )

# merging then grouping by quadrats and modality to compute total cover
pH_hab <- left_join(quad_sp_long, sp_tr_long, by="Species") %>%
  group_by(pH_hab, Quadrats, Trait, Modality) %>%
  summarize(cover=sum(cover)) %>% ungroup() %>% 
  mutate(ph= recode_factor(pH_hab,
                           "amb_cave"= "amb", "amb_deep_reef"="amb", "amb_reef"="amb", "amb_shallow_reef"="amb", 
                           "low_cave"= "low","low_deep_reef"="low", "low_reef" = "low", "low_shallow_reef"="low"),
         habitat=recode_factor(pH_hab,
                               "amb_cave"= "cave", "amb_deep_reef"="deep_reef", "amb_reef"="reef", "amb_shallow_reef"="shallow_reef", 
                               "low_cave"= "cave","low_deep_reef"="deep_reef", "low_reef" = "reef", "low_shallow_reef"="shallow_reef"))

## Organizing the dataset ----
table(pH_hab$Trait, pH_hab$Modality) ; data_model <- vector("list", 7)
# Reproduction
Reproduction = pH_hab %>% filter(., Trait == "agerepromaturity")
Reproduction_1 = Reproduction %>% filter(., Modality == "1") ; Reproduction_2 = Reproduction %>% filter(., Modality == "2")
matrice_reproduction = data.frame(quadrats = Reproduction_1$Quadrats, `1` = Reproduction_1$cover, `2` = Reproduction_2$cover) %>% column_to_rownames('quadrats')
matrice_reproduction = list(matrice_reproduction) ; names(matrice_reproduction) = "tr" ; rm(Reproduction_1, Reproduction_2)
variables <- lapply(rownames(matrice_reproduction$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[1]] <- lapply(matrice_reproduction, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Calcification
Calcification = pH_hab %>% filter(., Trait == "calcification")
Calcification_a = Calcification %>% filter(., Modality == "a") ; Calcification_b = Calcification %>% filter(., Modality == "b")
matrice_calcification = data.frame(quadrats = Calcification_a$Quadrats, a = Calcification_a$cover, b = Calcification_b$cover) %>% column_to_rownames('quadrats')
matrice_calcification = list(matrice_calcification) ; names(matrice_calcification) = "tr" ; rm(Calcification_a, Calcification_b)
variables <- lapply(rownames(matrice_calcification$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[2]] <- lapply(matrice_calcification, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Chem (What is that ???)
Chem = pH_hab %>% filter(., Trait == "chem")
Chem_1 = Chem %>% filter(., Modality == "1") ; Chem_2 = Chem %>% filter(., Modality == "2")
matrice_chem = data.frame(quadrats = Chem_1$Quadrats, `1` = Chem_1$cover, `2` = Chem_2$cover) %>% column_to_rownames('quadrats')
matrice_chem = list(matrice_chem) ; names(matrice_chem) = "tr" ; rm(Chem_1, Chem_2)
variables <- lapply(rownames(matrice_chem$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[3]] <- lapply(matrice_chem, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Feeding
Feeding = pH_hab %>% filter(., Trait == "feeding")
Feeding_a = Feeding %>% filter(., Modality == "a") ; Feeding_b = Feeding %>% filter(., Modality == "b") ; Feeding_c = Feeding %>% filter(., Modality == "c")
Feeding_d = Feeding %>% filter(., Modality == "d") ; Feeding_e = Feeding %>% filter(., Modality == "e") ; Feeding_h = Feeding %>% filter(., Modality == "h")
matrice_feeding = data.frame(quadrats = Feeding_a$Quadrats, a = Feeding_a$cover, b = Feeding_b$cover, c = Feeding_c$cover,
                             d = Feeding_d$cover, e = Feeding_e$cover, h = Feeding_h$cover) %>% column_to_rownames('quadrats')
matrice_feeding = list(matrice_feeding) ; names(matrice_feeding) = "tr" ; rm(Feeding_a, Feeding_b, Feeding_c, Feeding_d, Feeding_e, Feeding_h)
variables <- lapply(rownames(matrice_feeding$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[4]] <- lapply(matrice_feeding, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Form
Form = pH_hab %>% filter(., Trait == "form")
Form_a = Form %>% filter(., Modality == "a") ; Form_b = Form %>% filter(., Modality == "b") 
Form_c = Form %>% filter(., Modality == "c") ; Form_d = Form %>% filter(., Modality == "d") 
matrice_form = data.frame(quadrats = Form_a$Quadrats, a = Form_a$cover, b = Form_b$cover, c = Form_c$cover, d = Form_d$cover) %>% column_to_rownames('quadrats')
matrice_form = list(matrice_form) ; names(matrice_form) = "tr" ; rm(Form_a, Form_b, Form_c, Form_d)
variables <- lapply(rownames(matrice_form$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[5]] <- lapply(matrice_form, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Growth
Growth = pH_hab %>% filter(., Trait == "growth")
Growth_1 = Growth %>% filter(., Modality == "1") ; Growth_2 = Growth %>% filter(., Modality == "2") ; Growth_3 = Growth %>% filter(., Modality == "3")
matrice_growth = data.frame(quadrats = Growth_1$Quadrats, `1` = Growth_1$cover, `2` = Growth_2$cover, `3` = Growth_3$cover) %>% column_to_rownames('quadrats')
matrice_growth = list(matrice_growth) ; names(matrice_growth) = "tr" ; rm(Growth_1, Growth_2, Growth_3)
variables <- lapply(rownames(matrice_growth$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[6]] <- lapply(matrice_growth, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})
# Mobility
Mobility = pH_hab %>% filter(., Trait == "mobility")
Mobility_1 = Mobility %>% filter(., Modality == "1") ; Mobility_2 = Mobility %>% filter(., Modality == "2")
matrice_mobility = data.frame(quadrats = Mobility_1$Quadrats, `1` = Mobility_1$cover, `2` = Mobility_2$cover) %>% column_to_rownames('quadrats')
matrice_mobility = list(matrice_mobility) ; names(matrice_mobility) = "tr" ; rm(Mobility_1, Mobility_2)
variables <- lapply(rownames(matrice_mobility$tr), function(x) { infos[infos$Quadrats==x,]}) ; variables <- do.call(rbind,variables)
data_model[[7]] <- lapply(matrice_mobility, function(d) {
  dat_m <- data.frame(variables) ; tr <- as.matrix(round(d,0)) ; dat_m$tr <- tr ; dat_m$s <- rowSums(tr) ; dat_m})

# Remove additional datasets not useful anymore
rm(Calcification, Chem, Feeding, Form, Growth, Mobility, Reproduction, 
   matrice_calcification, matrice_chem, matrice_feeding, matrice_form, matrice_growth, matrice_mobility, matrice_reproduction)

## Model ----
# The model lasts more than 4 hours with a 2,3 GHz Dual-Core Intel Core i5 and 16 GB memory
# mn = vector("list", 7) ; for (i in 1:7) {
#  mn[[i]] <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = data_model[[i]]$tr, family = multinomial(), control=list(adapt_delta=0.99, max_treedepth=15), 
#                  chains = 2, cores = 2, iter = 4000, warmup = 1000, backend = "cmdstanr")}
load(file=file.path(getwd(),"/model/mn.RData"))

## Predictions ----
# The predictions take time as well
# Predicted_values = vector("list", 7) ; for (i in 1:7) {Predicted_values[[i]] = predict(mn[[i]], data_model[[i]]$tr)}
load(file=file.path(getwd(),"/model/Predicted_values.RData")) ; data_predicted = vector("list", 7)
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
                                 std.error = c(Predicted_values[[3]][1:294,2,1], Predicted_values[[3]][1:294,2,2]), 
                                 Q2.5 = c(Predicted_values[[3]][1:294,3,1], Predicted_values[[3]][1:294,3,2]), 
                                 Q97.5 = c(Predicted_values[[3]][1:294,4,1], Predicted_values[[3]][1:294,4,2]),
                                 Trait = rep("Chem", 294*1*2))
data_predicted[[4]] = data.frame(Habitat = rep(data_model[[4]]$tr$habitat, 6), pH = rep(data_model[[4]]$tr$pH, 6), 
                                 Condition = c(rep("a",294), rep("b",294), rep("c", 294), rep("d", 294), rep("e", 294), rep("h", 294)),
                                 Cover = c(Predicted_values[[4]][1:294,1,1], Predicted_values[[4]][1:294,1,2], Predicted_values[[4]][1:294,1,3], 
                                           Predicted_values[[4]][1:294,1,4], Predicted_values[[4]][1:294,1,5], Predicted_values[[4]][1:294,1,6]), 
                                 std.error = c(Predicted_values[[4]][1:294,2,1], Predicted_values[[4]][1:294,2,2], Predicted_values[[4]][1:294,2,3],
                                               Predicted_values[[4]][1:294,2,4], Predicted_values[[4]][1:294,2,5], Predicted_values[[4]][1:294,2,6]), 
                                 Q2.5 = c(Predicted_values[[4]][1:294,3,1], Predicted_values[[4]][1:294,3,2], Predicted_values[[4]][1:294,3,3],
                                          Predicted_values[[4]][1:294,3,4], Predicted_values[[4]][1:294,3,5], Predicted_values[[4]][1:294,3,6]), 
                                 Q97.5 = c(Predicted_values[[4]][1:294,4,1], Predicted_values[[4]][1:294,4,2], Predicted_values[[4]][1:294,4,3],
                                           Predicted_values[[4]][1:294,4,4], Predicted_values[[4]][1:294,4,5], Predicted_values[[4]][1:294,4,5]),
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
                                 Cover = c(Predicted_values[[6]][1:294,1,1], Predicted_values[[6]][1:294,1,2], Predicted_values[[6]][1:294,1,3]), 
                                 std.error = c(Predicted_values[[6]][1:294,2,1], Predicted_values[[6]][1:294,2,2], Predicted_values[[6]][1:294,2,3]), 
                                 Q2.5 = c(Predicted_values[[6]][1:294,3,1], Predicted_values[[6]][1:294,3,2], Predicted_values[[6]][1:294,3,3]), 
                                 Q97.5 = c(Predicted_values[[6]][1:294,4,1], Predicted_values[[6]][1:294,4,2], Predicted_values[[6]][1:294,4,3]),
                                 Trait = rep("Growth", 294*3*2))
data_predicted[[7]] = data.frame(Habitat = rep(data_model[[7]]$tr$habitat, 2), pH = rep(data_model[[7]]$tr$pH, 2), 
                                 Condition = c(rep("1",294), rep("2",294)),
                                 Cover = c(Predicted_values[[7]][1:294,1,1], Predicted_values[[7]][1:294,1,2]), 
                                 std.error = c(Predicted_values[[7]][1:294,2,1], Predicted_values[[7]][1:294,2,2]), 
                                 Q2.5 = c(Predicted_values[[7]][1:294,3,1], Predicted_values[[7]][1:294,3,2]), 
                                 Q97.5 = c(Predicted_values[[7]][1:294,4,1], Predicted_values[[7]][1:294,4,2]),
                                 Trait = rep("Mobility", 294*1*2))}
data_predicted = data_predicted %>% bind_rows()

## Vizualisation ----
data_predicted_viz = data_predicted %>% group_by(Habitat, pH, Condition, Trait) %>% summarise_all(mean) 
data_predicted_viz$Habitat = factor(data_predicted_viz$Habitat, levels=c('shallow_reef','cave','reef','deep_reef'))

(Final_Plot <- ggplot(data = data_predicted_viz, aes(x = pH, y = Cover, color = Condition, group = Condition)) +
  geom_line(linetype="dashed") + geom_point(size=5)+
  geom_errorbar(aes(ymin = Cover - std.error, ymax = Cover + std.error), width=.2, position=position_dodge(0.05)) +
  geom_rect(aes(xmin=which(levels(as.factor(pH)) == "ambient")-0.6, xmax=which(levels(as.factor(pH)) == "ambient")+0.5, ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE ) +
  geom_rect(aes(xmin=which(levels(as.factor(pH)) == "low")-0.5, xmax=which(levels(as.factor(pH)) == "low")+0.6, ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE) +
  scale_y_continuous(name = "Cover (%)", limits=c(0,100), breaks = c(0,25,50,75,100)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH")) + theme_bw() +
  theme(plot.title = element_text(size=14, hjust=0.5), axis.title.x = element_blank()) +
  facet_grid(Trait~Habitat, labeller = labeller(Habitat = c("cave" = "Cave", "deep_reef" = "Deep Reef", "reef" = "Reef", "shallow_reef" = "Shallow Reef"))) + 
  theme(strip.text = element_text(size = 14)))