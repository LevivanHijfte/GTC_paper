#### Library ####

library(sf)
library(jsonify)
library(tidyverse)
library(parallel)
library(fmesher)
library(sfheaders)
library(patchwork)

#### Functions ####

#merge objects in close vicinity
funtry1 <- function(var, distance = 30, tolerance = 2) {
  
  simplified_var <- st_simplify(var, dTolerance = tolerance)
  buffered_var   <- st_buffer(simplified_var, 30)
  overlap_list   <- st_overlaps(buffered_var)
  
  filtered_overlaps = seq_along(overlap_list) |> 
    lapply(function(x) c(x,overlap_list[[x]])) |> 
    lapply(function(x) BBB$LAMA$id[x])
  
  return(filtered_overlaps)
}

#
funtry2 <- function(var) {
  
  overlap_matrix <- sapply(var, function(x) sapply(var, function(y) any(x %in% y)))
  groups         <- unique(apply(overlap_matrix, 1, function(x) paste(which(x), collapse = ",")))
  result         <- lapply(groups, function(g) unique(unlist(var[as.numeric(unlist(strsplit(g, ",")))])))
  return(result)
}

#
funtry3 = function(object){
  
  f = list()
  
  for (i in 1:length(y)) {
    f[[i]] = lapply(1:length(y),function(x) any(y[[i]] %in% y[[x]])) |> 
      unlist() |> 
      which()
  }
  return(f)
}

#
funx2   = function(object = object, value_name = value_name){
  
  require("parallel")
  require("jsonify")
  
  if(!"sf" %in% class(object)){
    return("object is not of type sf")
  }
  
  proc =  mclapply(object$measurements, function(x) 
    if(is.na(x)){
      return(NA)
      
    } else{ 
      
      return(jsonify::from_json(x)[[value_name]])
      
    }
  )
  proc = mclapply(proc, function(x) ifelse(is.null(x), 0, x))
  
  return(unlist(proc))
}

####data####

dir          = "dir_GeoJSON_files_CD31/CD3/LAMA2_stainings"
files        = list.files(dir)
excludeLAMA2 = files[c(5,13)]
excludeCD31  = files[c(5,13)]

####lists for preprocessing####

windows = list()
cells   = list()
vessels = list()

####extract cell and ROI info####

for (i in seq_along(files)) {
  
  tmp                  = st_read(paste(dir, files[2], sep = "/"), crs = NA)
  tmp$annotationsimple = sapply(tmp$classification, function(x) ifelse(is.na(x), NA, from_json(x)$name))
  
  ann = tmp[tmp$objectType == "annotation",]
  ann = ann[st_area(ann$geometry)>100,]
  
  cellobj        = tmp[tmp$objectType == "cell",] |> 
    mutate(centroid = st_centroid(geometry))
  
  cellobj2 = data.frame(id    =  cellobj$id,
                        CD3   =  funx2(object = cellobj, value_name = "Red610: Nucleus: Mean")) |> 
    mutate(class_bio = case_when(CD3  >  300   ~ "T_cell", 
                                 .default = "Cell")) |> 
    right_join(cellobj, by = "id")
  i=2
  windows[[i]]    = tmp[is.na(tmp$classification)&tmp$objectType == "annotation",]
  cells[[i]]      = cellobj2
  vessels[[i]]    = list(LAMA  = ann[ann$annotationsimple == "LAMA2",] |> drop_na(id),
                         CD31  = ann[ann$annotationsimple == "CD31",]  |> drop_na(id))
  gc()
}

#-------------------------------------------
#### LAMA2 preprocessing ####
#-------------------------------------------

nonconvex = list()

for(i in seq_along(vessels)) {
  print(names(vessels)[i])
  if (names(vessels)[i] %in% excludeLAMA2) {
    
    next  
    
  }else{
    
    BBB = vessels[[i]]
    inb = funtry1(var = BBB$LAMA, distance = 30)
    x   = inb
    
    repeat {
      
      y = funtry2(var = x)
      
      if (all(unlist(lapply(funtry3(object = y), function(x) length(x) == 1)))) {
        break
      } 
      
      x = y
      print("it")
    }
    
    names(y) = 1:length(y)
    
    final = data.frame(group = unlist(lapply(1:length(y), function(x) rep(names(y)[x], length(y[[x]])))),
                       id    = unlist(y))
    
    HULL = BBB$LAMA |> 
      inner_join(final, by = "id") |>
      group_by(group) |> 
      summarise(do_union = TRUE) |> 
      mutate(nonconvex = mclapply(geometry, function(x) fm_nonconvex_hull(x, convex = -0.001, concave = 100) |>
                                    sf_remove_holes(), 
                                  mc.cores = 20, mc.cleanup = T) |> 
               unlist(recursive = FALSE) |> 
               st_sfc())
    
    parts     = st_cast(st_union(HULL$nonconvex),"POLYGON")
    cluster   = st_intersects(HULL$nonconvex, parts)
    
    if (any(which(sapply(cluster, function(x) length(x)>1)))) {
      
      double  = which(sapply(cluster, function(x) length(x)>1))
      b       = lapply(double, function(x) st_union(parts[cluster[[x]],])) |> unlist(recursive = FALSE) |> st_sfc()
      c       = st_sfc(c(parts, b))
      parts2  = c[-unlist(cluster[unlist(double)])]
      cluster = st_intersects(HULL$nonconvex, parts2)
      
    } else {print("no doublets")}
    
    nonconvex[[i]] = cbind(st_sf(nonconvex = HULL$nonconvex), cluster = unlist(cluster)) |> 
      group_by(cluster) |> 
      summarize()
    gc()
  }
}

names(windows)   = files
names(cells)     = files
names(vessels)   = files
names(nonconvex) = files
nonconvex        = nonconvex[!names(nonconvex) %in% excludeLAMA2]

#-------------------------------------------
#### save LAMA2 preprocessing ####
#-------------------------------------------

save(windows, cells, vessels, nonconvex, file = "dir_nonconvex_LAMA2_data.rdata")
