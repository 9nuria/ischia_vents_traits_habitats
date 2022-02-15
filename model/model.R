
library(reshape2)
library(brms)


load("data/DataLong.RData")
load("data/sp_tr.Rdata")

infos <- read.csv("data/from_nuria/Sites_Quadrats.csv")
infos$pH[which(infos$pH == "ambient1")] = "ambient"
infos$pH[which(infos$pH == "ambient2")] = "ambient"


dat <- dat[which(dat$species %in% rownames(sp_tr)),]

str(dat)
dat$quadrat = as.character(dat$quadrat)
dat$pH <- sapply(dat$quadrat, function(x) {infos[which(infos$Quadrats == x),]$pH})
dat$habitat <- rep(NA, nrow(dat))

#dat$habitat[which(dat$condition == "CA1")] = "cave"
#dat$habitat[which(dat$condition == "CA2")] = "cave"
#dat$habitat[which(dat$condition == "CL")] = "cave"
#dat$habitat[which(dat$condition  == "CoA1")] = "coralligenous"
#dat$habitat[which(dat$condition  == "CoA2")] = "coralligenous"
#dat$habitat[which(dat$condition  == "COL")] = "coralligeous"
#dat$habitat[which(dat$condition  == "RA")] = "reef"
#dat$habitat[which(dat$condition  == "RL")] = "reef"
#dat$habitat[which(dat$condition  == "SRA")] = "shallow_reef"
#dat$habitat[which(dat$condition  == "SRL")] = "shallow_reef"

tr <- sp_tr

matrices <- lapply(colnames(tr), function (x) {
  
  d <- dat  
  d$t <- sapply(d$species, function (y) {tr[y,x]})
  d <- d[which(!is.na(d$t)==T),]
  da <-  dcast(d, quadrat ~ t , sum, value.var = "cover")
  rownames(da) <- da$quadrat
  da <- da[,2:ncol(da)]
  da
})

variables <- lapply(rownames(matrices[[1]]), function(x) { infos[infos$Quadrats==x,]})
variables <- do.call(rbind,variables)

names(matrices) = colnames(tr)


data_models <- lapply(matrices, function(d) {
  
  dat_m <- data.frame(variables)
  tr <- as.matrix(round(d,0))
  dat_m$tr <- tr
  dat_m$s <- rowSums(tr)
  dat_m
  
})




models <- lapply(1:length(data_models), function(x) {
  
dat_m <- data_models[[x]]

mn <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = dat_m, family = multinomial(), control=list(adapt_delta=0.99, max_treedepth=15), 
                  chains = 2, cores = 2, iter = 4000, warmup = 1000,
                  inits = "0",
                  backend = "cmdstanr",
                  threads = 28)
  
  mn
  
})

save(models, file = "model/models_rand_2.RData")


load("model/models_rand_2.RData")


conditions <- make_conditions(models[[1]], "habitat")

conditional_effects(models[["form"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["feeding"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["growth"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["calcification"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["mobility"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["agerepromaturity"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["chem"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)

