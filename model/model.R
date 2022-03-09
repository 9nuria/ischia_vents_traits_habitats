
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
  
################################ Bayesian model 
# Model formula. This takes a while, except if run from server. Otherwise, load directly the output of the model. 
# The formula uses pH as a fixed factor, pH as a fixed factor and habitat as a random factor to account for the differences among habitats. 

mn <- brms::brm(tr | trials(s) ~ 1 + pH + (1 + pH | habitat), data = dat_m, family = multinomial(), control=list(adapt_delta=0.99, max_treedepth=15), 
                  chains = 2, cores = 2, iter = 4000, warmup = 1000,
                  inits = "0",
                  backend = "cmdstanr",
                  threads = 28)
  
  mn
  
})

save(models, file = "model/models_rand_2.RData")

######################################## Ploting the model (list of 7 models/variables)

load("model/models_rand_2.RData")

# Check output of the model 
#print(models, pars = "b")
#summary (models)
#coef(models)

conditions <- make_conditions(models[[1]], "habitat")

#order habitats by depth
conditions<-conditions%>%
  arrange(factor(habitat, levels=c("shallow_reef","cave","reef","deep_reef")))

# colors for categories

colors6<-c("#1b9e77", "#e7298a", "#7570b3", "#d95f02", "#0c2c84", "#f783be")
colors4b<-c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
colors3<-c("#1b9e77", "#7570b3", "#e7298a")
colors2<-c("#1b9e77", "#e7298a")
colorsph<-c("#bcd7e6", "#e31a1c")##ff9d87

#the 7 models for each category
# condiitional effects
conditional_effects(models[["form"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["feeding"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["growth"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["calcification"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["mobility"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["agerepromaturity"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
conditional_effects(models[["chem"]] , effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)

#plots

#model 1, form, 4 levels
ce.form<-conditional_effects(models[["form"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.form.df<-as.data.frame(ce.form$`pH:cats__`)
ce.form.df$cond__<-factor(ce.form.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                          labels=c("Shallow reef","Cave","Reef","Deep reef"))

  
plot1<-ggplot(ce.form.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
               xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
               ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors4b, name="Form", labels=c("Encrusting", "Filaments", "Massive", "Tree"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot1   

### plot final layout

p1<-ggplot(ce.form.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors4b,  name="", labels=c("Encrusting", "Filaments", "Massive", "Tree"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("", ""))+
  labs(title = "Morphological form")+
  theme_bw()+
  theme(plot.title = element_text(size=12, hjust=0.5),
        axis.title.x= element_blank())+
   guides(color=guide_legend(ncol=2))+
   facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
p1   



#model 2, feeding, 6 levels
 
ce.feed<-conditional_effects(models[["feeding"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.feed.df<-as.data.frame(ce.feed$`pH:cats__`)
ce.feed.df$cond__<-factor(ce.feed.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                          labels=c("Shallow reef","Cave","Reef","Deep reef"))


plot2<-ggplot(ce.feed.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors6, name="Feeding", labels=c("Autotroph", "Filter feeder", "Herbivor/Grazer", "Carnivor", "Detritivor", "Parasite"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot2   


#plot2 final layout

p2<-ggplot(ce.feed.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors6, name="", labels=c("Autotroph", "Filter feeder", "Herbivor/Grazer", "Carnivor", "Detritivor", "Parasite"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("", ""))+
  theme_bw()+
  guides(color=guide_legend(ncol=2))+
  labs(title = "Feeding")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                      strip.text.x = element_blank()))
p2   


#model 3, growth, 3 levels
ce.growth<-conditional_effects(models[["growth"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.growth.df<-as.data.frame(ce.growth$`pH:cats__`)
ce.growth.df$cond__<-factor(ce.growth.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                          labels=c("Shallow reef","Cave","Reef","Deep reef"))


plot3<-ggplot(ce.growth.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors3, name="Growth", labels=c("Low", "Moderate", "High"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot3  

#plot 3, final layout

p3<-ggplot(ce.growth.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors3, name="", labels=c("Low", "Moderate", "High"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("", ""))+
  theme_bw()+
  labs(title = "Growth")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                            strip.text.x = element_blank()))
p3  

#model 4, calcification, 2 levels
ce.ca<-conditional_effects(models[["calcification"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.ca.df<-as.data.frame(ce.ca$`pH:cats__`)
ce.ca.df$cond__<-factor(ce.ca.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                            labels=c("Shallow reef","Cave","Reef","Deep reef"))


plot4<-ggplot(ce.ca.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="Calcification", labels=c("No", "Yes"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot4  

########final layout

p4<-ggplot(ce.ca.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="", labels=c("No", "Yes"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme_bw()+
  labs(title = "Calcification")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                      strip.text.x = element_blank()))
p4 

#model 5, mobility, 2 levels

ce.mob<-conditional_effects(models[["mobility"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.mob.df<-as.data.frame(ce.mob$`pH:cats__`)
ce.mob.df$cond__<-factor(ce.mob.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                        labels=c("Shallow reef","Cave","Reef","Deep reef"))

plot5<-ggplot(ce.mob.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="Mobility", labels=c("Sessile", "Vagile"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot5  

#layout, plot 5
p5<-ggplot(ce.mob.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="", labels=c("Sessile", "Vagile"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("", ""))+
  theme_bw()+
  labs(title = "Mobility")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                      strip.text.x = element_blank()))
p5  


#model 6, age reproduct maturity, 2 levels

ce.rep<-conditional_effects(models[["agerepromaturity"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.rep.df<-as.data.frame(ce.rep$`pH:cats__`)
ce.rep.df$cond__<-factor(ce.rep.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                         labels=c("Shallow reef","Cave","Reef","Deep reef"))

plot6<-ggplot(ce.rep.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="Age reproductive \n maturity", labels=c("< 1 year", "> 1 Year"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot6 


# model 6, final layout

p6<-ggplot(ce.rep.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="", labels=c("< 1 year", "> 1 Year"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("", ""))+
  theme_bw()+
  labs(title = "Age reproductive maturity")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                      strip.text.x = element_blank()))
p6 


#model 7, chemistry, 2 levels

ce.chem<-conditional_effects(models[["chem"]], effects ="pH", conditions = conditions, re_formula = NULL, categorical = TRUE)
ce.chem.df<-as.data.frame(ce.chem$`pH:cats__`)
ce.chem.df$cond__<-factor(ce.chem.df$cond__, levels=c("habitat = shallow_reef", "habitat = cave","habitat = reef", "habitat = deep_reef"),
                         labels=c("Shallow reef","Cave","Reef","Deep reef"))

plot7<-ggplot(ce.chem.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=5)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="Chemical defenses", labels=c("Non", "Yes"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme (
    axis.title.x=element_blank())+
  theme_bw()+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.text=element_text(size=14)))
plot7 

#plot 7, final layout

p7<-ggplot(ce.chem.df,aes(x=effect1__, y=estimate__, color=effect2__, group=effect2__)) + 
  geom_point(size=3)+ 
  geom_errorbar(aes(ymin = estimate__-se__, ymax = estimate__+se__),width=0.5, size = 0.5)+
  geom_line(linetype="dashed")+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="ambient")-0.6,
                xmax=which(levels(as.factor(effect1__))=="ambient")+0.5,
                ymin=-Inf, ymax=Inf), 
            fill="#deebf7", color="NA", alpha= 0.05,inherit.aes=FALSE )+
  geom_rect(aes(xmin=which(levels(as.factor(effect1__))=="low")-0.5,
                xmax=which(levels(as.factor(effect1__))=="low")+0.6,
                ymin=-Inf, ymax=Inf),
            fill="#fdbb84", color="NA", alpha= 0.05, inherit.aes=FALSE)+
  scale_color_manual(values= colors2, name="", labels=c("Non", "Yes"))+
  scale_y_continuous(name ="Probabiity", limits=c(0,1.00), breaks = c(0, 0.25,0.50,0.75,1.00)) +
  scale_x_discrete(name="", labels= c("Ambient", "Low pH"))+
  theme_bw()+
  labs(title = "Chemical defenses")+
  theme (
    plot.title = element_text(size=12, hjust=0.5),
    axis.title.x=element_blank())+
  facet_wrap(~cond__, nrow=1)+ (theme(strip.background = element_blank(),
                                      strip.text.x = element_blank()))
p7 




# layout 7 plots

all7<-(p1/p2/p3/p4/p5/p6/p7)+
  plot_layout(widths=1, heights=1)
all7

all4<-(p1/p2/p3/p4)
all4


### save

# folder to save plot as png ----
dir_plot<-"./plot"

ggsave("fig5_functions_7plots.png", plot= all7, device="png", height=25, width=20, units="cm", dpi=300)
ggsave("fig5_functions_4plots.png", plot= all4, device="png", height=25, width=20, units="cm", dpi=300)

