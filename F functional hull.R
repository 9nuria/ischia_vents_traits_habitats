#new figure func space with only 2pH levels


## preparing #####
rm(list=ls())

# libraries ---
library(tidyverse)
library(mFD)


# folder with ready to use data ----
dir_data<-"./data"
dir_results<-"./FD"


# loading dataset built in script A& B ----
load(file.path(dir_data, "habph_fe_cover.Rdata"))
load(file.path(dir_results,"fe_4D_coord.Rdata"))

# computing mean cover of FEs in the 2 pH levels: Ambient and Low pH
habph2_fe_cover<-rbind( 
  cave_low = habph_fe_cover["cave_low",],
  cave_amb = apply(  habph_fe_cover[c("cave_ambient1", "cave_ambient2"),], 2, mean),
  
  deep_reef_low = habph_fe_cover["deep_reef_low",],
  deep_reef_amb = apply(  habph_fe_cover[c("deep_reef_ambient1", "deep_reef_ambient2"),], 2, mean),
  
  reef_low = habph_fe_cover["reef_low",],
  reef_amb = apply(  habph_fe_cover[c("reef_ambient1", "reef_ambient2"),], 2, mean),
  
  shallow_reef_low = habph_fe_cover["shallow_reef_low",],
  shallow_reef_amb = apply(  habph_fe_cover[c("shallow_reef_ambient1", "shallow_reef_ambient2"),], 2, mean)
  
)
  

# compute FRic, FDis and FIde for all habitats * 2 levels of pH ---
habph2_multidimFD<-alpha.fd.multidim(sp_faxes_coord = fe_4D_coord, 
                                    asb_sp_w = habph2_fe_cover,
                                    ind_vect = c("fric", "fdis", "fide", "fdiv"), 
                                    scaling = TRUE, details_returned = TRUE
)



# setting parameters for plot ####

# folder to save plot as png ----
dir_plot<-"./plot"
root_dir<-getwd()
setwd(dir_plot)

# color code for habitat and pH from Nuria ----

hab_ph2<-c("shallow_reef_amb","shallow_reef_low",
            "cave_amb", "cave_low",
            "reef_amb", "reef_low",
            "deep_reef_amb", "deep_reef_low")

vcolors<-c("#93a1fa", "#f7d305",
           "#6478f5", "#f5a511",
           "#3953f7", "#f78e0c",
           "#0219ad", "#f7560c")
names(vcolors)<-hab_ph2

# coordinates of all species
pool_coord<-habph2_multidimFD$details$sp_faxes_coord

# vertices of all fe in 4D ----
pool_vert_nm<-habph2_multidimFD$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(pool_coord[,1:4])
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1

# indices values
habph2_fd<-habph2_multidimFD$functional_diversity_indices


## plotting along pairs of axes ####

pairs_axes<-list( c(1,2), c(3,4) )

for ( z in 1:length(pairs_axes) )
  {
  
  # names of axes   
  xy<-pairs_axes[[z]]
  
# list to store ggplot
  ggplot_list<-list()
  
  
  for (v in hab_ph2 ) {
    # v="deep_reef_low"
    
    # color for habitat*pH levels
    col_v<-as.character(vcolors[v])
    
    # species present in v
    sp_v<-names(which(habph2_multidimFD$details$asb_sp_occ[v,]==1))
    
    # background with axes range set + title
    ggplot_v<-background.plot(range_faxes=range_axes,
                              faxes_nm=paste0("PC", xy), 
                              color_bg="grey95")
    ggplot_v<-ggplot_v + labs(subtitle=v)
    
    # convex hull of species pool
    ggplot_v<-pool.plot(ggplot_bg=ggplot_v,
                        sp_coord2D=pool_coord[,xy],
                        vertices_nD=pool_vert_nm,
                        plot_pool=FALSE,
                        color_ch=NA, fill_ch="white", alpha_ch=1
    )
    
    # plot convex hull of assemblage but not species
    ggplot_v<-fric.plot( ggplot_bg=ggplot_v, 
                         asb_sp_coord2D=list(vv=pool_coord[sp_v,xy]),
                         asb_vertices_nD=list(vv=habph2_multidimFD$details$asb_vert_nm[[v]]),
                         plot_sp = FALSE,
                         color_ch=c(vv=col_v),
                         fill_ch=c(vv=col_v),
                         alpha_ch=c(vv=0.1)
    )
    
    
    # plot species weights using plot.fide without showing mean value
    ggplot_v<- fide.plot(ggplot_bg=ggplot_v,
                         asb_sp_coord2D=list(vv=pool_coord[sp_v,xy]),
                         asb_sp_relatw=list(vv=habph2_multidimFD$details$asb_sp_relatw[v,sp_v]),
                         asb_fide_coord2D=list(vv=habph2_fd[v,paste0("fide_PC", xy)]),
                         plot_sp = TRUE,
                         shape_sp = c(vv=21),
                         color_sp =c(vv=col_v),
                         fill_sp = c(vv=paste0(col_v,"70") ),
                         shape_fide = c(vv=23),
                         size_fide = c(vv=1),
                         color_fide = c(vv=col_v),
                         fill_fide = c(vv=col_v),
                         color_segment = c(vv=col_v),
                         width_segment = c(vv=0.5),
                         linetype_segment = c(vv=1)
    )
    
    
    # ggplot_v storing in list
    ggplot_list[[v]]<-ggplot_v
    
  }# end of v
  
  
  # patchwork of plots : 4 habitats (rows) * 2 columns (pH)
  FD_xy<- (ggplot_list[[1]] + ggplot_list[[2]] ) / 
    (ggplot_list[[3]] + ggplot_list[[4]] ) / 
    (ggplot_list[[5]] + ggplot_list[[6]] ) / 
    (ggplot_list[[7]] + ggplot_list[[8]] )
  
  
  # saving as png ----
  sz<-3
  ggsave(FD_xy, filename = paste0("FD_pc",xy[1],"vs",xy[2],"low_vs_amb.png"),
         device = "png",width=2*sz, height=4*sz)
  
  
  
}# end of axes_id



setwd(root_dir)

