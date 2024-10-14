####data####

tumors       = readRDS("/mnt/data/users/levi/snRNAseqGLASS/primary_tumors.RDS")
cell_marker  = read_xlsx("d:/Documents/Results/Gemistocytes/manuscript/Submission/Nature Communications/7.Supplementary_tables.xlsx", sheet = "Table S10", skip = 2)

####functions####

cellscore <- function(object, markers, nbins = 30, assay = "SCT", slot = "data") {
  for (l in 1:length(unique(markers$Celltype))) {
    
    name      = unique(markers$Celltype)[l]    
    marker    = markers[markers$Celltype == name,]$gene
    
    data      = GetAssayData(object, assay = assay, slot = slot)
    genes     = rownames(data)
    markersub = marker[marker %in% genes]
    a         = colMeans(data[markersub,]) 
    c         = data.frame(gene = rownames(data), value = rowMeans(data)) |> 
      remove_rownames() |> 
      mutate(bins = ntile(value, n = nbins))
    
    genesa    = list()
    
    for(i in 1:length(markersub)){
      bin         = c[c$gene %in% markersub[i],]$bin
      genesa[[i]] = c |> 
        dplyr::filter(bins == bin) |> 
        sample_n(size = 100) |> 
        pull(gene)
    }
    
    f = mclapply(genesa, function(a) colMeans(data[a,]))
    f = as.data.frame(do.call(cbind, f))
    g = data.frame(sample = rownames(f), test = a, control = rowMeans(f)) |> 
      remove_rownames() |> 
      mutate(result = test - control)
    
    stopifnot(all(Cells(object) == g$sample))
    object[[name]] =  g$result
    
  }
  return(object)
}

####Calculate enrichment score####

for (i in 1:length(tumors)) {
tumors[[i]] = cellscore(object = tumors[[i]], markers = c) 
}

for (i in 1:length(tumors)) {
  a = cowplot::plot_grid(DimPlot(tumors[[i]], label = T) + 
                           ggtitle(tumors[[i]]$state[1]), 
                         DotPlot(tumors[[i]],features = celltypes, col.min = 0) + 
                           theme(axis.text.x = element_text(angle = 60, hjust=1)) &
                           scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))),
                         ncol = 1)
  
  b = cowplot::plot_grid(a, FeaturePlot(tumors[[i]], features = celltypes, order = T, ncol = 3) &
                           scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))), ncol = 2, rel_widths = c(1,2))
  
  ggsave(plot = b, filename = paste0("dir_visualization",paste(tumors[[i]]$state[1], "pdf", sep = ".")), 
         height = 12, 
         width = 20)
}

ot = "Oligo-like tumor state"
at = "Astro-like tumor state"
st = "stem-like tumor state"
ut = "undetermined tumor state"
m  = "TAMs"
o  = "Oligodendrocytes"
a  = "Astrocytes"
e  = "Endothelial/Pericytes"
d  = "G1/S/G2/M"
n  = "Neurons"
l  = "T cells"

bioidents = rbind(

#1
data.frame(tumor     = "1", 
           seurat_clusters = 0:(length(unique(tumors[[1]]$seurat_clusters))-1), 
           bioidents = c(ot, at, m,  ot, 
                         at, at, at, ot, 
                         at, at, o,  ot,
                         at, ot, at, d, 
                         a,  e)),

#2
data.frame(tumor     = "2", 
           seurat_clusters = 0:(length(unique(tumors[[2]]$seurat_clusters))-1), 
           bioidents = c(ot, at, at, at, #3
                         m,  at, ot, ot, #7
                         at, n,  o,  m,  #11
                         ot, a,  n,  n,  #15
                         n,  n,  st, n)),  

#3
data.frame(tumor     = "3", 
           seurat_clusters = 0:(length(unique(tumors[[3]]$seurat_clusters))-1), 
           bioidents = c(ot, at, at, m,  #3
                         at, m,  o,  at, #7
                         ot, o,  at, a,  #11
                         at, at, st, e)),

#4
data.frame(tumor     = "4", 
           seurat_clusters = 0:(length(unique(tumors[[6]]$seurat_clusters))-1), 
           bioidents = c(at, m,  at, ot, #3
                         at, at, ot, at, #7
                         ot, ot, at, at, #11
                         at, at, ut, o, #15
                         n,  d,  a,  ot, #19
                         e)),

#5_1
data.frame(tumor     = "5_1", 
           seurat_clusters = 0:(length(unique(tumors[[7]]$seurat_clusters))-1), 
           bioidents = c(m,  ot, o,  m, #3
                         o,  ot, m,  o, #7
                         o,  at, at, at, #11
                         m,  ot, m,  at,
                         l,  st, n)),

#6
data.frame(tumor     = "6", 
           seurat_clusters = 0:(length(unique(tumors[[8]]$seurat_clusters))-1), 
           bioidents = c(m,  at, at, ot, #3
                         ot, ot, at, at, #7
                         at, m,  n,  o, #11
                         n,  a,  ot, m, #15
                         n,  n)),

#5_2
data.frame(tumor     = "5_2", 
           seurat_clusters = 0:(length(unique(tumors[[9]]$seurat_clusters))-1), 
           bioidents = c(o,  o,  o,  o, #3
                         m,  o,  o,  m, #7
                         n,  at, ot, o, #11
                         n,  o,  ot, n, #15
                         e,  ot))
)

for (i in 1:length(tumors)) {

tumorname     = tumors[[i]]$state[1]
clusteridents = bioidents[bioidents$tumor == tumorname,-1] |> 
  mutate(seurat_clusters = as.factor(seurat_clusters))
  
tumors[[i]]@meta.data = tumors[[i]]@meta.data |> 
  rownames_to_column("cellid") |>  
  left_join(clusteridents, by = "seurat_clusters") |> 
  column_to_rownames("cellid")
}
for (i in 1:length(tumors)) {
 Idents(tumors[[i]]) = "bioidents" 
}

saveRDS(tumors, "Dir_annotated_tumors.RDS")
