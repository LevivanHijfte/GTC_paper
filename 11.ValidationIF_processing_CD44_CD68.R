####library####

library(furrr)
library(spatstat)
library(tidyverse)
library(fpc)
library(sf)
library(flux)
library(jsonify)
library(progressr)
library(readxl)
library(stars)
library(parallel)

####functions####
#this function moves xy coordinates for a ppp objects 
funx = function (x, xcorrection, ycorrection, centre = c(0,0), angle = 0, slot = "geometry"){
  dat        = x
  correct    = x[[slot]]
  corrected  = correct
  
  for (j in 1:length(correct)){
    for (i in 1:length(correct[[j]])) {
      for (k in 1:length(correct[[j]][[i]])) {
        corrected[[j]][[i]][[k]][,1] = correct[[j]][[i]][[k]][,1] + xcorrection
        corrected[[j]][[i]][[k]][,2] = correct[[j]][[i]][[k]][,2] + ycorrection
      }
    }
  }
  
  dat[[slot]] = corrected
  
  return(dat)
}

funx2 = function(object, value_name){
  
  require("parallel")
  
  if(!"sf" %in% class(object)){
    return("object is not of type sf")
  }
  
  proc =  mclapply(object$measurements, function(x) 
    if(is.na(x)){
      return(NA)
    } else{ 
      return(jsonify::from_json(x)[[value_name]])
    },
    mc.cores = 15, mc.cleanup = T)
  
#  if(any(lapply(proc, function(x) is.null(x)))) {warning("Intensity values are missing")}
  
  proc = mclapply(proc, function(x) ifelse(is.null(x), 0, x),mc.cores = 15, mc.cleanup = T)
  return(unlist(proc))
}

####data####

dir        = "dir_GeoJSON_files"
metadat   = readxl::read_xlsx("Table_S7", skip = 2)
panel2samp = metadat$Panel_2[!metadat$Panel_2 == "NA"]
panel1samp = metadat$Panel_1[!metadat$Panel_1 == "NA"]

####Panel_2####

ppp_objects_CD44 = list()
ppp_border_CD44  = list()

for(i in seq_along(panel2samp)){

  CD3value      = "FITC: Nucleus: Mean"
  CD44value     = "Red610: Cytoplasm: Mean"
  CD3threshold  = metadat[metadat$Panel_2 == panel2samp[i],]$Panel_2_CD3_nucleus_mean
  CD44threshold = metadat[metadat$Panel_2 == panel2samp[i],]$Panel_2_CD44_cytoplasm_mean
  
  obj   = st_read(paste(dir, paste0(panel2samp[i], ".geojson"), sep = ""), crs = NA)
  cells = data.frame(id =   obj$id,
                     CD44 = funx2(object = obj, value_name = CD44value),
                     CD3  = funx2(object = obj, value_name = CD3value))
  
  # Create annotation for all cells in the GeoJSON object
  
  cells = cells |> 
    mutate(CD44 = as.numeric(CD44),
           CD3  = as.numeric(CD3),
           CD44call = CD44 >= CD44threshold, 
           CD3call  = CD3  >= CD3threshold, 
           classification = case_when(CD44call ==T & CD3call == T ~ "T_cell",
                                      CD44call ==T & CD3call == F ~ "Gemistocyte",
                                      CD44call ==F & CD3call == T ~ "T_cell",
                                      CD44call ==F & CD3call == F ~ "Cell", 
                                      .default = "Cell"))
  
  geom    = obj$geometry
  dat     = as.data.frame(st_coordinates(st_centroid(geom)))
  dat$id  = obj$id
  
  names(dat) = c("x", "y", "id")
  windowid = obj$id[1]
  
  cells  = cells[cells$id != windowid,]
  dat    = dat[dat$id     != windowid,]
  window = as.owin(obj[[6]][[1]])
  inside = inside.owin(x = dat$x, y = dat$y, w = window)
  dat2   = dat[inside,]
  cells2 = cells[inside,]
  stopifnot(all(dat2$id == cells2$id))
  
  ppp_objects_CD44[[i]]      = ppp(x      = dat2$x, y = dat2$y, 
                                   marks  = as.factor(cells2$classification), 
                                   window = window)
  
  ppp_border_CD44[[i]] = obj$geometry[1]

}

####Panel_1####

ppp_objects_CD68    = list()
ppp_border_CD68     = list()

for(i in seq_along(panel1samp)){

  CD3value      = "Cy5: Nucleus: Mean"
  CD68value     = "Aqua: Cytoplasm: Max"
  CD3threshold  = metadat[metadat$Panel_1 == panel1samp[i],]$Panel_1_CD3_nucleus_mean
  CD68threshold = metadat[metadat$Panel_1 == panel1samp[i],]$Panel_1_CD68_cytoplasm_max
  
  obj = st_read(paste(dir, paste0(panel1samp[i], ".geojson"), sep = "/"), crs = NA)
  
  if (metadat$x[i] != "NA") {
    obj = funx(x = obj, xcorrection = as.numeric(metadat[metadat$Panel_1 == panel1samp[i],]$x), ycorrection = as.numeric(metadat[metadat$Panel_1 == panel1samp[i],]$y)) 
  } else {print("No XY manipulation")}  
  cells = data.frame(id =   obj$id,
                     CD68 = funx2(object = obj, value_name = CD68value),
                     CD3  = funx2(object = obj,  value_name = CD3value))
  
  # Create annotation for all cells in the GeoJSON object
  cells = cells |> 
    mutate(CD68 = as.numeric(CD68),
           CD3  = as.numeric(CD3),
           CD68call = CD68 >= CD68threshold, 
           CD3call  = CD3  >= CD3threshold, 
           classification = case_when(CD68call ==T & CD3call == T ~ "T_cell",
                                      CD68call ==T & CD3call == F ~ "TAM",
                                      CD68call ==F & CD3call == T ~ "T_cell",
                                      CD68call ==F & CD3call == F ~ "Cell", 
                                     .default = "Cell"))
  
  geom   = obj$geometry
  dat    = as.data.frame(st_coordinates(st_centroid(geom)))
  dat$id = obj$id

  names(dat) = c("x", "y", "id")
  windowid   = obj$id[1]
  
  cells  = cells[cells$id != windowid,]
  dat    = dat[dat$id     != windowid,]
  window = as.owin(obj[[6]][[1]])
  inside = inside.owin(x = dat$x, y = dat$y, w = window)
  dat2   = dat[inside,]
  cells2 = cells[inside,]
  stopifnot(all(dat2$id == cells2$id))
  x          = ppp(x      = dat2$x, y = dat2$y, 
                   marks  = as.factor(cells2$classification), 
                   window = window)
  
  if(metadat$Centerx[i] != "NA"){
    x = spatstat.geom::rotate.ppp(X      =  x, 
                                  centre =  c(as.numeric(metadat$Centerx[i]),as.numeric(metadat$Centery[i])), 
                                  angle  =  as.numeric(metadat$angle[i]), 
                                  rescue =  T)
  }else{print("No rotation")}
  
  ppp_objects_CD68[[i]] = x
  ppp_border_CD68[[i]]  = obj$geometry[1]
  gc()
}
####processing panel 2####

density_val_CD44_T_cells  = list()
density_CD44              = list()
for(i in 1:length(ppp_objects_CD44)){
  tmp1 = st_as_stars(density(subset(ppp_objects_CD44[[i]], marks == "Gemistocyte"), sigma = 500, dimyx = 500))
  tmp2 = st_as_sf(tmp1) 
  
  tmp3 = st_as_stars(density(ppp_objects_CD44[[i]], sigma = 500, dimyx = 500))
  tmp4 = st_as_sf(tmp3) 
  
  density_CD44[[i]] = data.frame(tmp2$geometry, v = tmp2$v/tmp4$v)
  
  density_val_CD44_T_cells[[i]] = data.frame(x    = ppp_objects_CD44[[i]]$x, 
                                             y    = ppp_objects_CD44[[i]]$y, 
                                             mark = ppp_objects_CD44[[i]]$marks)
} 

q = list()
for (j in 1:length(density_CD44)) {
  mat = matrix(nrow = length(density_CD44[[j]][[2]]), ncol = 3)
  colnames(mat) = c("xmean", "ymean", "v")
  for (i in 1:length(density_CD44[[j]][[2]])) {
    mat[i, "xmean"] = mean(density_CD44[[j]][[1]][[i]][[1]][,1])
    mat[i, "ymean"] = mean(density_CD44[[j]][[1]][[i]][[1]][,2])
    mat[,"v"]      = density_CD44[[j]][[2]]
  }
  q[[j]] = as.data.frame(mat)
}

for (i in 1:length(density_CD44)) {
  density_val_CD44_T_cells[[i]]$v = NA
  for(j in 1:nrow(density_val_CD44_T_cells[[i]])){
    x = density_CD44[[i]][[2]][which(q[[i]]$xmean == q[[i]]$xmean[which.min(abs(q[[i]]$xmean - density_val_CD44_T_cells[[i]]$x[j]))] &
                                     q[[i]]$ymean == q[[i]]$ymean[which.min(abs(q[[i]]$ymean - density_val_CD44_T_cells[[i]]$y[j]))])]
    if (length(x) == 0) {
           density_val_CD44_T_cells[[i]]$v[j] = NA  
    } else{density_val_CD44_T_cells[[i]]$v[j] = x
    }
  }
}

####processing panel 1####

density_val_CD44_TAMs    = list()
density_CD44_cons_slides = list()
TAMppplist               = list()

for (i in seq_along(datacd44$Panel_1)) {
  obj1 = ppp_objects_CD44[[i]]
  obj2 = ppp_objects_CD68[[i]]
  if(length(obj1) <1) {print("false")}
  else{
    new_window =  intersect.owin(as.owin(obj2), as.owin(obj1))
    
    TAMppplist[[i]]  =  ppp(x      = obj2$x, 
                            y      = obj2$y, 
                            marks  = obj2$marks, 
                            window = new_window)
  }
}

for (i in seq_along(datacd44$Panel_1)) {
  obj1 = ppp_objects_CD44[[i]]
  obj2 = ppp_objects_CD68[[i]]
  if(length(obj1) <1) {print("false")
  }else{
    new_window =  intersect.owin(as.owin(obj2), as.owin(obj1))
    
    gemppp  =  ppp(x      = obj1$x, 
                   y      = obj1$y, 
                   marks  = obj1$marks, 
                   window = new_window)
    
    TAMppp  =  TAMppplist[[i]]
    
    tmp1 = st_as_stars(density(subset(gemppp, marks == "Gemistocyte"), sigma = 500, dimyx = 500))
    tmp2 = st_as_sf(tmp1) 
    
    tmp3 = st_as_stars(density(gemppp, sigma = 500, dimyx = 500))
    tmp4 = st_as_sf(tmp3) 
    
    density_try2 = data.frame(tmp2$geometry, v = tmp2$v/tmp4$v)
    
    u      = data.frame(x    = TAMppp$x, 
                        y    = TAMppp$y, 
                        mark = TAMppp$marks)
    
    mat = matrix(nrow = length(density_try2[[1]]), ncol = 3)
    colnames(mat) = c("xmean", "ymean", "v")
    
    for (j in 1:length(density_try2[[1]])) {
      mat[j, "xmean"] = mean(density_try2[[1]][[j]][[1]][,1])
      mat[j, "ymean"] = mean(density_try2[[1]][[j]][[1]][,2])
      mat[,"v"]       = density_try2[[2]]
    }
    v = as.data.frame(mat)
    
    u$v = NA
    for(j in 1:nrow(u)){
      x = density_try2[[2]][which(v$xmean == v$xmean[which.min(abs(v$xmean - u$x[j]))] &
                                  v$ymean == v$ymean[which.min(abs(v$ymean - u$y[j]))])]
      if (length(x) == 0) {
        u$v[j] = NA  
      } else{u$v[j] = x
      }
      
    }
    density_CD44_cons_slides[[i]] = density_try2
    density_val_CD44_TAMs[[i]]    = u
  }
  print("true")
}

names(ppp_objects_CD44)         = metadat$Names_data
names(ppp_border_CD44)          = metadat$Names_data
names(ppp_objects_CD68)         = metadat$Names_data
names(ppp_border_CD68)          = metadat$Names_data
names(density_val_CD44_T_cells) = metadat$Names_data
names(density_CD44)             = metadat$Names_data
names(density_CD44_cons_slides) = metadat$Names_data
names(density_val_CD44_TAMs)    = metadat$Names_data

saveRDS(ppp_objects_CD44,         file = "dir_ppp_objects_CD44.rds")
saveRDS(ppp_border_CD44,          file = "dir_ppp_border_CD44.rds")
saveRDS(ppp_objects_CD68,         file = "dir_ppp_objects_CD68.rds")
saveRDS(ppp_border_CD68,          file = "dir_ppp_border_CD68.rds")
saveRDS(density_CD44,             file = "dir_density_CD44.rds")
saveRDS(density_val_CD44_T_cells, file = "dir_density_val_CD44_T_cells.rds")
saveRDS(density_CD44_cons_slides, file = "dir_density_CD44_cons_slides.RDS")
saveRDS(density_val_CD44_TAMs,    file = "dir_density_val_CD44_TAMs.rds")

