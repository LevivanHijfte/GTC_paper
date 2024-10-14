library(infercnv)
library(ensembldb)
library(Matrix)
library(tidyverse)
library(Seurat)

tumors         = readRDS("dir_annotated_tumor_files.RDS")
normalcells    = c("TAMs","Oligodendrocytes","Astrocytes","Endothelial/Pericytes","Neurons","T cells")

####data processing####

for(i in seq_along(tumors)){
dat = tumors[[i]]@assays$RNA@counts

####cell annitation####

cellannotation = tumors[[i]]@meta.data %>% 
                 rownames_to_column("cellid") %>% 
                 dplyr::select(cellid, bioidents)

cellannotation = cellannotation[match(cellannotation$cellid, colnames(dat)),]
colnames(cellannotation) = NULL
write.table(cellannotation, file = "dir_to_cellannotation.txt", sep = "\t",row.names = F)

####run inferCNV####

refgroups = normalcells[normalcells %in% unique(tumors[[i]]@meta.data$bioidents)]
outdir    = paste("dir_to_store_cnv_inferences", tumors[[i]]@meta.data$state[1], sep = "")

if (!dir.exists(outdir)) {
  dir.create(outdir)
} 

obj = CreateInfercnvObject(raw_counts_matrix = dat,
                           annotations_file  =  "dir_to_cellannotation.txt",
                           delim             = "\t",
                           gene_order_file   = "import_gene_order_file",
                           ref_group_names   = refgroups)

infercnv_obj = infercnv::run(obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=outdir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             plot_steps = F,
                             HMM=T, 
                             num_threads = 5)
}
