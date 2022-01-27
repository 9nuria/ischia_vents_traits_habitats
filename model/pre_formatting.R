#library load
library(reshape2)

#Following the discussion to have equal quadrats for all the sites at the Castello, the low condition includes only south side (S2). Then, Ambient 1 is Ambient south (S1), and Ambient 2 (N1) is Ambient north.
#Then, we need to remove the quadrats low north (N2): 1N2, 2N2, 3N2, 4N2, 5N2, 6N2, 7N2, 8N2, 9N2, 10N2,11N2 and 12N2.


# load files
cast <- read.csv("data/from_nuria/Data_Cover_tCastello.csv")
cave <- read.csv("data/from_nuria/Data_Cover_tCaves.csv")
chiane <- read.csv("data/from_nuria/Data_Cover_tChiane.csv")
cor <- read.csv("data/from_nuria/Data_Cover_tCoralligenous.csv")

infos <- read.csv("data/from_nuria/Sites_Quadrats.csv")

#change them into long format
cast.l <- melt(cast, id = c("X"))
cave.l <- melt(cave, id = c("X"))
chiane.l <- melt(chiane, id = c("X"))
cor.l <- melt(cor, id = c("X"))

colnames(cast.l) <- c("quadrat", "species", "cover")
colnames(cave.l) <- c("quadrat", "species", "cover")
colnames(chiane.l) <- c("quadrat", "species", "cover")
colnames(cor.l) <- c("quadrat", "species", "cover")


md.cast <- do.call(rbind, lapply(as.character(cast.l$quadrat), function(x) {
  
  infos[infos$Quadrats == x, 2:6]                       
  
}))

md.cave <- do.call(rbind, lapply(as.character(cave.l$quadrat), function(x) {
  
  infos[infos$Quadrats == x, 2:6]                       
  
}))

md.chiane <- do.call(rbind, lapply(as.character(chiane.l$quadrat), function(x) {
  
  infos[infos$Quadrats == x, 2:6]                       
  
}))

md.cor <- do.call(rbind, lapply(as.character(cor.l$quadrat), function(x) {
  
  infos[infos$Quadrats == x, 2:6]                       
  
}))


fl.cast <- cbind(md.cast, cast.l)
fl.cave <- cbind(md.cave, cave.l)
fl.chiane <- cbind(md.chiane, chiane.l)
fl.cor <- cbind(md.cor, cor.l)

dat <- rbind(fl.cast, fl.cave, fl.chiane, fl.cor) 

save(dat, file = "data/DataLong.RData")
