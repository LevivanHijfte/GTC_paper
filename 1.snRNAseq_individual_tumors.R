####libraries####
library(Seurat)
library(sctransform)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(patchwork)
library(corrplot)
library(SingleCellExperiment)
library(scDblFinder)
library(cowplot)

####functions####

doubletcheck_preprocessing  <- function(location, project)
{
      doubletcheck = Read10X(data.dir = location) |> 
                     CreateSeuratObject(min.cells = 3, min.features = 200, project = project) |> 
                     SCTransform(vst.flavor = "v2", verbose = FALSE) |> 
                     RunPCA(npcs = 30, verbose = F) |> 
                     RunUMAP(dims = 1:30, verbose = F) |> 
                     FindNeighbors(dims = 1:30, verbose = F) |> 
                     FindClusters(resolution = 1, verbose = F)
                     
      return(doubletcheck)
}

find_doublets = function(object, seed = 123456) { 
      dbl      = Seurat::as.SingleCellExperiment(object)
      top.mam  = scran::getTopHVGs(dbl, prop=0.1)
      dbl.dens = scDblFinder::computeDoubletDensity(dbl,
                                                     subset.row=top.mam, 
                                                     d=ncol(reducedDim(dbl)))
      dbl$DoubletScore = dbl.dens
      stopifnot(colnames(dbl) == colnames(object))
      object$DoubletScore = dbl$DoubletScore
      rm(dbl, dbl.dens, top.mam)
      gc()
      object$DoubletScoreLog = log1p(object$DoubletScore)
      return(object)
}

visualize_doublets <- function(object, intercept = 0){
      a = FeaturePlot(object, "DoubletScore")
      b = ggplot(data.frame("doublets" = object$DoubletScoreLog), aes(doublets)) +
          geom_density() +
          theme_bw() +
          geom_vline(xintercept = intercept)
      return(a+b+plot_layout(ncol = 2))
}

load_seuratobject <- function(location, doubletnames, project, min.cells = 3, min.features = 200) {
      object = Read10X(data.dir = location)
      object = object[, !colnames(object) %in% doubletnames]
      object = CreateSeuratObject(counts = object, 
                                  min.cells = min.cells, 
                                  min.features = min.features, 
                                  project = project)
      mito.features_object = grep(pattern = "^MT-", x=rownames(x=object), value=T)
      percent.mito_object  = Matrix::colSums(x = GetAssayData(object = object, slot="counts")[mito.features_object,]) / Matrix::colSums(x = GetAssayData(object = object, slot = "counts"))
      object[["percent.mito"]] = percent.mito_object
      return(object)
}

filterandnormalizeobject <- function(object, state, minnFeature_RNA, minnCount_RNA, maxpercent.mito = 0.1, nfeatures = 2000) {
  object = subset(x = object, subset = nFeature_RNA > minnFeature_RNA & nCount_RNA >  minnCount_RNA & percent.mito <maxpercent.mito)
  object = SCTransform(object, vst.flavor = "v2", verbose = FALSE) |> 
           RunPCA(npcs = 30, verbose = F) |> 
           RunUMAP(dims = 1:30, verbose = F) |> 
           FindNeighbors(dims = 1:30, verbose = F) |> 
           FindClusters(resolution = 1, verbose = F)
  object[["state"]] = state
  return(object)
}

####markers####

standardmarkers = c("P2RY12", "GFAP", "EGFR", "RBFOX3", "CD2", "MOG", "CD163", "CD34", "PDGFRB", "CD83", "CCL3", "CCL4")

####loading datasets####

#####1#####
# Loading and pre-processing dataset 1
location1    = "path_to_matrix_1"
dbltch       = doubletcheck_preprocessing(location = location1, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 1.25)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.5])

object_1     = load_seuratobject(location     = location1,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_1, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(1200))
b = VlnPlot(object = object_1, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(2000))
c = VlnPlot(object = object_1, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_1 = filterandnormalizeobject(object          = object_1, 
                                    state           = "1", 
                                    minnFeature_RNA = 1200, 
                                    minnCount_RNA   = 2000)

#####2#####
# Loading and pre-processing dataset 2
location2    = "path_to_matrix_2"
dbltch       = doubletcheck_preprocessing(location = location2, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 1.5)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.5])

object_2     = load_seuratobject(location     = location2,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_2, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(1200))
b = VlnPlot(object = object_2, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(1500))
c = VlnPlot(object = object_2, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_2 = filterandnormalizeobject(object          = object_2, 
                                    state           = "2", 
                                    minnFeature_RNA = 1200, 
                                    minnCount_RNA   = 1500)

#####3#####
# Loading and pre-processing dataset 3
location3    = "path_to_matrix_3"
dbltch       = doubletcheck_preprocessing(location = location3, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 1.75)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.75])
  
object_3     = load_seuratobject(location     = location3,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_3, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(800))
b = VlnPlot(object = object_3, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(1000))
c = VlnPlot(object = object_3, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_3 = filterandnormalizeobject(object          = object_3, 
                                    state           = "3", 
                                    minnFeature_RNA = 800, 
                                    minnCount_RNA   = 1000)

#####4#####
# Loading and pre-processing dataset 4
location4    = "path_to_matrix_4"
dbltch       = doubletcheck_preprocessing(location = location6, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 1.5)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.5])

object_4     = load_seuratobject(location     = location4,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_4, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(800))
b = VlnPlot(object = object_4, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(800))
c = VlnPlot(object = object_4, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_4 = filterandnormalizeobject(object          = object_4, 
                                    state           = "4", 
                                    minnFeature_RNA = 800, 
                                    minnCount_RNA   = 800)

#####5#####
# Loading and pre-processing dataset 5
location5    = "path_to_matrix_5"
dbltch       = doubletcheck_preprocessing(location = location5, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 1.25)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.25])

object_5     = load_seuratobject(location     = location5,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_5, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(600))
b = VlnPlot(object = object_5, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(600))
c = VlnPlot(object = object_5, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_5 = filterandnormalizeobject(object          = object_5, 
                                    state           = "5_1", 
                                    minnFeature_RNA = 600, 
                                    minnCount_RNA   = 600)

#####6####
# Loading and pre-processing dataset 6
location6    = "path_to_matrix_6"
dbltch       = doubletcheck_preprocessing(location = location6, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept =1.6)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 1.6])

object_6     = load_seuratobject(location     = location6,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_6, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(800))
b = VlnPlot(object = object_6, features = c("nCount_RNA"),   pt.size = 0.5) + scale_y_log10() + geom_hline(yintercept = c(1000))
c = VlnPlot(object = object_6, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_6 = filterandnormalizeobject(object          = object_6, 
                                    state           = "6", 
                                    minnFeature_RNA = 1000, 
                                    minnCount_RNA   = 1300)

#####7####
# Loading and pre-processing dataset 7
location7    = "path_to_matrix_7"
dbltch       = doubletcheck_preprocessing(location = location7, project = "glioma")
dbltch       = find_doublets(dbltch)
visualize_doublets(dbltch, intercept = 2)
doubletnames = names(dbltch$DoubletScoreLog[dbltch$DoubletScoreLog > 2])

object_7     = load_seuratobject(location     = location7,
                                 doubletnames = doubletnames, 
                                 project      = "glioma")

a = VlnPlot(object = object_7, features = c("nFeature_RNA"), pt.size = 0.5) + geom_hline(yintercept = c(1000))
b = VlnPlot(object = object_7, features = c("nCount_RNA"),   pt.size = 0.5) + geom_hline(yintercept = c(1500))
c = VlnPlot(object = object_7, features = c("percent.mito"), pt.size = 0.5)

a+b+c+plot_layout(ncol = 3)

object_7 = filterandnormalizeobject(object          = object_7, 
                                    state           = "5_2", 
                                    minnFeature_RNA = 1400, 
                                    minnCount_RNA   = 1700)


states = c(object_1[["state"]][[1]][1],
           object_2[["state"]][[1]][1],
           object_3[["state"]][[1]][1],
           object_4[["state"]][[1]][1],
           object_5[["state"]][[1]][1],
           object_6[["state"]][[1]][1],
           object_7[["state"]][[1]][1])

saveRDS(file = "/mnt/data/users/levi/snRNAseqGLASS/primary_tumors.RDS",
                         object = c(object_1,
                                    object_2,
                                    object_3,
                                    object_4, 
                                    object_5,
                                    object_6,
                                    object_7))
