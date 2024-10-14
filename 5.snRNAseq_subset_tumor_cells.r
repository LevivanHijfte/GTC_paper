####data####

tumors     = readRDS("dir_annotated_tumors.RDS")

annotation = data.frame()
for (i in seq_along(tumors)) {
  a = data.frame(cellid = rownames(tumors[[i]]@meta.data), state = tumors[[i]]@meta.data$state, bioidentsintegrated = tumors[[i]]@meta.data$bioidents)
  annotation = rbind(annotation,a)
}

####functions####

load_seuratobject <- function(location, cellnames, project, min.cells = 3, min.features = 200) {
  object = Read10X(data.dir = location)
  object = object[, colnames(object) %in% cellnames]
  object = CreateSeuratObject(counts = object, 
                              min.cells = min.cells, 
                              min.features = min.features, 
                              project = project)
  mito.features_object = grep(pattern = "^MT-", x=rownames(x=object), value=T)
  percent.mito_object = Matrix::colSums(x = GetAssayData(object = object, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object, slot = "counts"))
  object[["percent.mito"]] = percent.mito_object
  return(object)
}

normalizeobject <- function(object, state, nfeatures = 3000) {
  object = SCTransform(object, vst.flavor = "v2", verbose = FALSE)|> 
    RunPCA(npcs = 30, verbose = F)
  object[["state"]] = state
  return(object)
}
####data preprocessing####

allstates = annotation |> 
            mutate(cellid = substr(cellid, 1,nchar(cellid)-2)) |> 
            filter(bioidentsintegrated %in% c("Oligo-like tumor state", 
                                              "stem-like tumor state",
                                              "Astro-like tumor state", 
                                              "G1/S/G2/M"))


location1 = "path_to_matrix_1"
object_1  = load_seuratobject(location = location1, cellnames = allstates[allstates$state == "1", ]$cellid, project = "glioma")
object_1  = normalizeobject(object_1, state = "1", nfeatures = 3000)

location2 =  "path_to_matrix_2"
object_2  = load_seuratobject(location = location2, cellnames = allstates[allstates$state == "2", ]$cellid, project = "glioma")
object_2  = normalizeobject(object_2, state = "2", nfeatures = 3000)

location3 = "path_to_matrix_3"
object_3  = load_seuratobject(location = location3, cellnames = allstates[allstates$state == "3", ]$cellid, project = "glioma")
object_3  = normalizeobject(object_3, state = "3", nfeatures = 3000)

location4 = "path_to_matrix_4"
object_4  = load_seuratobject(location = location4, cellnames = allstates[allstates$state == "4", ]$cellid, project = "glioma")
object_4  = normalizeobject(object_4, state = "4", nfeatures = 3000)

location5 ="path_to_matrix_5"
object_5  = load_seuratobject(location = location5, cellnames = allstates[allstates$state == "5_1", ]$cellid, project = "glioma")
object_5  = normalizeobject(object_5, state = "5_1", nfeatures = 3000)

location6 = "path_to_matrix_6"
object_6 = load_seuratobject(location = location6, cellnames = allstates[allstates$state == "6", ]$cellid, project = "glioma")
object_6 = normalizeobject(object_6, state = "6", nfeatures = 3000)

location7 = "path_to_matrix_7"
object_7 = load_seuratobject(location = location7, cellnames = allstates[allstates$state == "5_2", ]$cellid, project = "glioma")
object_7 = normalizeobject(object_7, state = "5_2", nfeatures = 3000)


####integrate data####

reference.list    <- c(object_1, object_2, object_3, object_4, object_5, object_6, object_7) # Makes a list of the objects you want to merge
features          <- SelectIntegrationFeatures(object.list = reference.list, nfeatures = 3000)
ifnb.list         <- PrepSCTIntegration(object.list = reference.list, anchor.features = features)
glioma.anchors    <- FindIntegrationAnchors(object.list = ifnb.list, 
                                            anchor.features = features, 
                                            normalization.method = "SCT", reduction = "rpca") # Finds the common sources of variation
glioma.integrated <- IntegrateData(anchorset = glioma.anchors, normalization.method = "SCT") # Integrates the two datasets

# Post-processing merged data
glioma.integrated <- RunPCA(object = glioma.integrated, verbose = F)
ElbowPlot(glioma.integrated, ndims = 40)
glioma.integrated <- RunUMAP(object = glioma.integrated, reduction = "pca", dims = 1:40) # Again, you can play around with nrs of dims
glioma.integrated <- FindNeighbors(object = glioma.integrated, dims = 1:40)
glioma.integrated <- FindClusters(object = glioma.integrated, resolution = .6)

DefaultAssay(object = glioma.integrated) = "integrated"
saveRDS(glioma.tumor, "dir_tumor_cells.RDS")
