library(lidR)
library(sf)
library(sp)
library(tidyverse)
library(sfheaders)
library(data.table)
library(purrr)
 
library(RStoolbox)
library(tidyterra)
library(car)

library(rasterVis)
library(grid)
library('latticeExtra')



las_list <- read_rds('las_list.rds')  
rastimagelist <- read_rds('rastimagelist.rds')  



##### all trees 
data2 <- read.csv('E:/OneDrive - caa.columbia.edu/songzhaoying/matou/matouGWAS/fam location all trees.csv')
data2$fam <- as.factor(data2$fam)
data1 <- data2  
shp <- sfheaders::sf_point(data1, x = "x", y = "y",z='z',keep = TRUE )


library(terra)
shp <- sf::st_set_crs(shp,32650)

extr_matou <- function(las_list,rastimagelist){
  lapply(las_list, function(x){
    expr <- tryCatch({
      library("lidR")
      library("rgdal")
      library(sfheaders)
      library(stars)
      library(raster)
      library(tidyverse)
      library(sf)
      library(data.table)
      tictoc:: tic("processing las file")
      opt_output_files(x) <-paste0(tempdir(),rnorm(1) ,"/{ORIGINALFILENAME}_{ID}")
      opt_chunk_size(x) <- 0
      opt_chunk_buffer(x) <- 20
      classified_ctg <- classify_ground(x, csf())
      dtm <- rasterize_terrain(classified_ctg, 1, tin())
      tictoc:: toc()
      tictoc:: tic("processing normalize_h")
      ctg_norm2 <- normalize_height(classified_ctg, dtm)
      opt_select(ctg_norm2) <- "xyz"
      opt_filter(ctg_norm2) <- "-drop_z_below 0"
      crff  <- shp$fam
      library(data.table)
      crown_se  <- lapply(crff, function(fdx){
        expr <- tryCatch({
          sf2  <-  shp %>% dplyr:: filter(fam == fdx)
          opt_output_files(ctg_norm2) <-''
          subset3 <- clip_roi(ctg_norm2, sf2, radius=3)
          # las3 = readLAS(subset3)
          lasrff <-  filter_poi(subset3, Z>=2)
          plot_crowns <- plot_metrics(lasrff, func = .stdtreemetrics,sf2,radius =15)
          message(paste0(fdx,':area=',plot_crowns$convhull_area)) 
          #  las3 = subset3
          # chm <- rasterize_canopy(las3, 0.2, p2r(0.15))
          # ttops2 <- locate_trees(las3, lmf(ws= 6 , hmin = 2))
          # algo2 <- dalponte2016(chm, ttops2)
          # ctg_segm <- segment_trees(las3, algo2)
          # crown_polo = crown_metrics(ctg_segm, func = .stdtreemetrics, geom = "convex")
          fer <- payload(subset3)  %>% dplyr::mutate(treeID=fdx)
          fer$treeID <- as.factor(fer$treeID)
          names(fer) <- c('x','y','z','treeID')
          fer <- as.data.frame(fer)
          message(paste0('project',fdx))
          
          fd1 <- rastimagelist[names(rastimagelist) == x$type]
          
          lapply(fd1, function(ralis){
            expr <- tryCatch({
              dsta <- terra::extract(ralis,fer[,c('x','y')],xy=T ) %>% mutate(treeID=fer$treeID,
                                                                              Z=fer$z,
                                                                              tree_H =plot_crowns$Z ,
                                                                              area=plot_crowns$convhull_area
              ) 
              message(paste0(fdx,':area=',plot_crowns$convhull_area)) 
              return(dsta)
            })})
          
        },error = function(e){
          NULL
        })
      }) %>% invoke(rbind,.)
      return(crown_se)
    },error = function(e) {
      message('Caught an error!')
      cat("ERROR :", conditionMessage(e), "/n")
      print(e)},
    print(paste("processing Loop_", x$type, sep = "_")))
  })
}


extr_matou(las_list, rastimagelist = rastimagelist ) %>%
  saveRDS('data_21.rds')


