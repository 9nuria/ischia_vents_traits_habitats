
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



priors_mod <- lapply(data_models, function(dat_m) {
  
  mn <- brms::get_prior(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = dat_m, family = multinomial())
  
  mn
  
})



models <- lapply(1:length(data_models), function(x) {
  
dat_m <- data_models[[x]]

p <- priors_mod[[x]]
  
mn <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = dat_m, family = multinomial(), control=list(adapt_delta=0.99, max_treedepth=15), 
                  chains = 2, cores = 2, iter = 4000, warmup = 1000,
                  priors = p,
                  inits = "0",
                  backend = "cmdstanr",
                  threads = 28)
  
  mn
  
})

save(models, file = "model/models_rand_2.RData")


load("model/models_rand_2.RData")

re_formula=“NULL”

summary(models[[1]])
bayes_R2(models[[1]])

conditions <- make_conditions(models[[2]], "habitat")
conditional_effects(models[[2]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)

conditional_effects(models[[1]] , categorical = TRUE,  effects ="habitat")
conditional_effects(models[[3]] , categorical = TRUE
conditional_effects(models[[4]] , categorical = TRUE)
conditional_effects(models[[5]] , categorical = TRUE)
conditional_effects(models[[6]] , categorical = TRUE)
conditional_effects(models[[7]] , categorical = TRUE)
conditional_effects(models[[8]] , categorical = TRUE)
conditional_effects(models[[9]] , categorical = TRUE)
conditional_effects(models[[10]] , categorical = TRUE)
conditional_effects(models[[11]] , categorical = TRUE)
conditional_effects(models[[12]] , categorical = TRUE)
conditional_effects(models[[13]] , categorical = TRUE)


