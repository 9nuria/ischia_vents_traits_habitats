library("hypervolume")

data_volume_low = list()
data_volume_amb = list()
levels_trait    = list()
chull.coords    = list(list(), list(), list(), list(), list(), list(), list())
hpts            = list(list(), list(), list(), list(), list(), list(), list())
pts             = list(list(), list(), list(), list(), list(), list(), list())
chull.area      = list(list(), list(), list(), list(), list(), list(), list())
chull.poly      = list(list(), list(), list(), list(), list(), list(), list())
xy.coords       = list(list(), list(), list(), list(), list(), list(), list())

# Define the functional Volume
VTot                   = hypervolume(data_trait_occupancy[,8:10])@Volume
ATot                   = sp::Polygon(data_trait_occupancy[,8:9][c(chull(data_trait_occupancy[,8:9]), 
                                     chull(data_trait_occupancy[,8:9])[1]),], hole = F)@area
Trait                  = colnames(data_trait_occupancy)[1:7]
hypervolume            = vector("list", 2)
hypervolume[[1]]       = matrix(nrow = 6, ncol = 7)
hypervolume[[2]]       = matrix(nrow = 6, ncol = 7)

# Ambient pH Conditions
for (i in 1:7) {
  data_volume_amb[[i]] = data_trait_occupancy_amb %>% dplyr::select(1:11) %>% group_split(Trait[i])
  data_volume_amb[[i]] = data_volume_amb[[i]][[1]] %>% data.frame()
  levels_trait[[i]]    = data_trait_occupancy %>% dplyr::select(., Trait[i]) %>% group_split(Trait[i])
  levels_trait[[i]]    = levels_trait[[i]][[1]] %>% data.frame() 
  levels_trait[[i]]    = unique(levels_trait[[i]][,1]) %>% as.vector() 
  for (j in 1:length(levels_trait[[i]])) {
    if (nrow(data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]),]) < 1) {
      hypervolume[[1]][j,i]  <- NA
    } else if (nrow(data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]),]) <= 2) {
      hypervolume[[1]][j,i]  <- 0 
    } else {
      hpts[[i]][[j]]         <- chull(x = data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]), 8], 
                                      y = data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]), 9])
      hpts[[i]][[j]]         <- c(hpts[[i]][[j]], hpts[[i]][[j]][1]) 
      xy.coords[[i]][[j]]    <- cbind(data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]), 8],
                                      data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == levels_trait[[i]][j]), 9])
      chull.coords[[i]][[j]] <- xy.coords[[i]][[j]][hpts[[i]][[j]],]
      chull.poly[[i]][[j]]   <- sp::Polygon(chull.coords[[i]][[j]], hole = F) 
      hypervolume[[1]][j,i]  <- round((chull.poly[[i]][[j]]@area / ATot)*100, 2) 
      # hypervolume[[1]][j,i] = round((((hypervolume(data_volume_amb[[i]][which(data_volume_amb[[i]][,i] == 
      #                                levels_trait[[i]][j]), c(8:10)])@Volume) / VTot) * 100),2) 
    }}}

# Low pH Conditions
for (i in 1:7) {
  data_volume_low[[i]] = data_trait_occupancy_low %>% dplyr::select(1:11) %>% group_split(Trait[i])
  data_volume_low[[i]] = data_volume_low[[i]][[1]] %>% data.frame()
  levels_trait[[i]]    = data_trait_occupancy %>% dplyr::select(., Trait[i]) %>% group_split(Trait[i])
  levels_trait[[i]]    = levels_trait[[i]][[1]] %>% data.frame() 
  levels_trait[[i]]    = unique(levels_trait[[i]][,1]) %>% as.vector() 
  for (j in 1:length(levels_trait[[i]])) {
    if (nrow(data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]),]) < 1) {
        hypervolume[[2]][j,i]  <- NA
      } else if (nrow(data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]),]) <= 2) {
        hypervolume[[2]][j,i]  <- 0 
      } else {
        hpts[[i]][[j]]         <- chull(x = data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]), 8], 
                                        y = data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]), 9])
        hpts[[i]][[j]]         <- c(hpts[[i]][[j]], hpts[[i]][[j]][1]) 
        xy.coords[[i]][[j]]    <- cbind(data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]), 8],
                                        data_volume_low[[i]][which(data_volume_low[[i]][,i] == levels_trait[[i]][j]), 9])
        chull.coords[[i]][[j]] <- xy.coords[[i]][[j]][hpts[[i]][[j]],]
        chull.poly[[i]][[j]]   <- sp::Polygon(chull.coords[[i]][[j]], hole = F) 
        hypervolume[[2]][j,i]  <- round((chull.poly[[i]][[j]]@area / ATot)*100, 2) 
      # hypervolume[[2]][j,i] = round((((hypervolume(data_volume_low[[i]][which(data_volume_low[[i]][,i] == 
      #                               levels_trait[[i]][j]), c(8:11)], method = "svm")@Volume) / VTot) * 100),2)
    }}}
