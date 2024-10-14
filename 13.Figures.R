####library####
library(VennDiagram)
library(ComplexHeatmap)
library(RColorBrewer)
library(sctransform)
library(tidyverse)
library(patchwork)
library(parallel)
library(circlize)
library(corrplot)
library(ggridges)
library(ggrepel)
library(viridis)
library(rstatix)
library(Seurat)
library(readxl)
library(DESeq2)
library(ggpubr)
library(ggbeeswarm)
library(infercnv)
library(WGCNA)
library(CellChat)
library(Seurat)
library(ggalluvial)
library(NMF)
library(furrr)
library(spatstat)
library(tidyverse)
library(fpc)
library(sf)
library(flux)
library(jsonify)
library(stars)
library(data.table)

#### Orders ####
.reorder_celltypes   = c("Total", "Tumor_cells", "TAMs", "CD4_T_cells", "CD8_T_cells", "B_cells") 
.reorder_celltypes2  = c("Total_cells", "Tumor_cells", "Macrophages", "T_cells", "CD4_T_cells", "CD8_T_cells", "B_cells") 
.reorder_celltypes3  = c("B_cells", "CD8_T_cells", "CD4_T_cells", "TAMs", "Tumor_cells", "Total") 
.reorder_celltypes4  = c("Total", "Tumor_cells","TAMs", "T_cells", "CD4_T_cells", "CD8_T_cells", "B_cells")
.reorder_GTC         = c("-", "+", "++", "+++")
.comparisons_GTC     = list(c("-","+"), c("-","++"), c("-","+++"), c("+++","+"), c("++","+"), c("++","+++")) 

#### Colors ####
.col_ptid = c("102" = '#f7fcfd', "103" = '#66c2a4', "105" = '#00441b', "108" = '#fff7ec',
              "110" = '#fc8d59', "111" = '#7f0000', "113" = '#fff7fb', "115" = '#74a9cf', 
              "117" = '#023858', "118" = '#fff7f3', "121" = '#f768a1', "122" = '#49006a',
              "124" = '#ffffff', "125" = '#ccece6', "126" = '#238b45', "129" = '#fdd49e',
              "131" = '#d7301f', "133" = '#d0d1e6', "136" = '#0570b0', "137" = '#fcc5c0',
              "139" = '#ae017e', "141" = '#252525', "143" = '#e5f5f9', "145" = '#41ae76',
              "146" = '#fee8c8', "147" = '#ef6548', "148" = '#ece7f2', "149" = '#3690c0',
              "150" = '#fde0dd', "156" = '#dd3497', "157" = '#737373', "158" = '#99d8c9',
              "160" = '#006d2c', "162" = '#fdbb84', "166" = '#b30000', "172" = '#a6bddb',
              "174" = '#045a8d', "177" = '#fa9fb5', "178" = '#7a0177')

.col_celltype  = c("Total_cells"       = "#b30000", "Tumor_cells"       = "#fc8d59", 
                   "Macrophages"       = "#045a8d", "T_cells"           = "#00441b", 
                   "Thelper_cells"     = "#238b45", "Cytotoxic_T_cells" = "#66c2a4", 
                   "B_cells"           = "#fed976")

.col_cuff_GTC   = c("-"   = '#e5f5f9', "+"   = '#99d8c9',
                   "++"  = '#238b45',"+++" = '#00441b')

.col_cuff3     = c("Cuff: -"   = '#e5f5f9',"Cuff: +"   = '#99d8c9',
                   "Cuff: ++"  = '#238b45',"Cuff: +++" = '#00441b',
                   "Cuff: NA"  = "white")

#### data ####
##### VECTRA ####
data        = read_xlsx("Table S9", skip = 2)
metadata    = read_xlsx("Table S1", skip = 2)

##### Bulk RNA-seq ####
DE_test     = read_rds("output_DE_test_script.RDS")
depathways  = read_xlsx("Table S13", skip = 2)

##### snRNA-seq data ####
glioma.integrated     = readRDS("output_integrated_seurat_object_all_cells")
glioma.tumor          = readRDS("output_integrated_seurat_object_tumor_cells")
TAMs                  = readRDS("output_integrated_seurat_object_TAMs")
cellChat              = readRDS("output_integrated_seurat_CellChat_analysis")

##### Spatial protein data ####
dat.norm.neg          = readRDS("Output_normalization_NS_protein_data")
dat.norm.quant        = readRDS("Output_normalization_NS_RNA_data")
annotation_NS_protein = read_xlsx("Table S5",  skip = 2)
annotation_NS_RNA     = read_xlsx("Table S3",  skip = 2)
modulesspatial        = read_xlsx("Table S14", skip = 2)

#### Validation stainings ####
cd44data               = readRDS("dir_ppp_objects_CD44.rds")
cd68data               = readRDS("dir_ppp_objects_CD68.rds")
cd44border             = readRDS("dir_ppp_border_CD44.rds")
density_try            = readRDS("dir_density_CD44.rds")
s                      = readRDS("dir_density_val_CD44_T_cells.rds")
z                      = readRDS("dir_density_val_CD44_TAMs.rds")
CRYAB_data             = read.csv("dir_CRYAB_data_export") 
annotation_CD44_CD68   = read_xlsx("Table S7", skip = 2) |> 
                         filter(Panel_1 != "NA")
load("dir_nonconvex_LAMA2_data.rdata") # data LAMA2_CD31_CD3 staining

##### Markers ####
markers     = read_xlsx("Table S10", skip = 2)
reastdf     = read_xlsx("Table S11", skip = 2)[,1:2]

#### Functions ####
# Calculate enrichment scores for single cell RNA-sequencing data
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

#calculates z-score for inter-cell distance values in PPP objects
funtest = function(object, k = 1, not_interested = NA, n = 1, reps = 100, cores = 20){
  require(tidyverse)
  require(parallel)
  require(data.table)
  require(spatstat)
  #test random
  exclude         = names(which(table(object$marks)<n))
  marks(object)   = ifelse(object$marks %in% c(exclude, not_interested), "exclude", as.character(object$marks))
  
  samp            = replicate(n = reps, object$marks[sample(seq_along(object$marks))], simplify = F)
  
  distrand        = mclapply(1:length(samp), function(x) as.data.frame(nndist(X = object, by = as.factor(samp[[x]]), k = k)), mc.cores = cores)
  
  rand_mean_sd    = rbindlist(distrand, use.names = T, idcol = "rep") |> 
    data.frame(marks = unlist(samp)) |> 
    filter(marks != "exclude") |>
    pivot_longer(cols = !c(marks,rep), names_to = "to", values_to = "distance_rand") |> 
    mutate(marks_to = paste(marks, to, sep = "-")) |> 
    select(rep, marks_to, distance_rand) |> 
    group_by(rep, marks_to) |> 
    summarize(dist_rand_mean = mean(distance_rand),
              dist_rand_sd   = sd(distance_rand)) |> 
    group_by(marks_to) |> 
    summarize(sd_rand   = mean(dist_rand_sd), 
              mean_rand = mean(dist_rand_mean))
  
#test nn
  distinf = nndist(object, by = as.factor(object$marks), k = k) |> 
    data.frame(marks = object$marks) |> 
    filter(marks != "exclude") |>
    pivot_longer(cols = ! marks,names_to = "to", values_to = "distance") |> 
    mutate(marks_to = paste(marks, to, sep = "-")) |> 
    select(marks_to, distance) |> 
    left_join(rand_mean_sd, by = "marks_to") |> 
    group_by(marks_to) |> 
    mutate(z      = (distance - mean_rand)/sd_rand, 
           median = median(z))
  
  return(distinf)
}

#### Figures ####
#### 1 ####
##### 1c ####
#count matrix heatmap
Counts_hm = data[, c("ROI_nr","res","Study_ID", "CD8_T_cells","CD4_T_cells","B_cells","Macrophages","Tumor_cells")]

# select Study_ID's with a median of 0 tumor cells
notum = Counts_hm |> 
        dplyr::select(!c(ROI_nr, res)) |>
        group_by(Study_ID) |>
        summarise_all(median) |> 
        filter(Tumor_cells == 0) |>
        pull(Study_ID) 

#filter samples without tumor cells
Counts_hm = Counts_hm[!Counts_hm$Study_ID %in% notum,] 

#logtransform and scale data
hmdat_hm1 = Counts_hm |>
            distinct(ROI_nr, .keep_all = T) |> 
            select(!c(res, Study_ID)) |> 
            column_to_rownames("ROI_nr") |>
            log1p() |>
            scale() |>
            t() 

#create metadata table
metadatahm_hm1 = Counts_hm |> 
                 mutate(T_cells = CD4_T_cells + CD8_T_cells) |> 
                 left_join(metadata, by = c("ROI_nr", "Study_ID")) |>
                 distinct(.keep_all = T) |>
                 mutate("log1p(T_cells)" = log1p(T_cells), 
                        Patient = substr(Study_ID, 1,3)) |> 
                 dplyr::filter(cuff_score != "NA") 

#order metadata table
metadatahm_hm1 = metadatahm_hm1[match(colnames(hmdat_hm1), metadatahm_hm1$ROI_nr),]

#sort on Total T cells
reorderded_T = order(metadatahm_hm1$`log1p(T_cells)`, decreasing = T)

# color gradient for T cells
col_fun = circlize::colorRamp2(c(0, 8), c("#e5f5f9", "#00441b"))

## Annotation heatmap
an_hm1 = ComplexHeatmap::columnAnnotation(
         "log1p(T_cells)" = metadatahm_hm1$`log1p(T_cells)`,
         Cuff             = metadatahm_hm1$cuff_score,
         Study_ID         = metadatahm_hm1$Patient,
         WHO_2021         = metadatahm_hm1$WHO_2021,
         Resection        = metadatahm_hm1$res, 
         Sex              = metadatahm_hm1$Sample_Sex, 
         CNV              = metadatahm_hm1$CNV_load, 
         na_col           = c("white"),
         col = list("log1p(T_cells)" = col_fun,
                    Study_ID         = ..col_ptid,
                    WHO_2021         = c("A2" = "#deebf7", "A3" = "#6baed6", "A4" = "#08306b", "NA" = "white"),
                    Resection        = c("1" = "#deebf7", "2" = "#6baed6", "3" = "#08306b", "4" = "grey60"),
                    Sex       = c("F" = "#fdd49e", "M" = "#fff7ec"),
                    CNV         = c("high" = "#fdd49e", "low" = "#fff7ec","NA" = 'white'),
                    Cuff       = c("NA" = 'white', "-" = '#e5f5f9',"+" = '#99d8c9',"++" = '#238b45',"+++" = '#00441b')))

#plot heatmap counts
hm1 = ComplexHeatmap::Heatmap(hmdat_hm1, 
                              column_order       = reorderded_T,
                              name               = "Counts (scaled)",
                              top_annotation     = an_hm1, 
                              column_names_rot   = 90, 
                              height             = unit(5, "cm"),
                              column_names_gp    = grid::gpar(fontsize = 2),
                              row_title          = "Count data",
                              column_dend_height = unit(0, "pt"),
                              show_column_names  = F)

colnames_hm1 = colnames(hmdat_hm1)
colorder_hm1 = hm1 |> ComplexHeatmap::column_order() 

#data heatmap inter-cell distances
hmdat_2 = data |>
          filter(ROI_nr %in% colnames_hm1) |> 
          select(c(ROI_nr,
                   distance_score_from_B_cells_to_CD4_T_cells, 
                   distance_score_from_B_cells_to_B_cells,
                   distance_score_from_B_cells_to_Tumor_cells,
                   distance_score_from_CD4_T_cells_to_CD4_T_cells,
                   distance_score_from_CD4_T_cells_to_B_cells,
                   distance_score_from_CD4_T_cells_to_Tumor_cells)) |>
          rename_with(~ sub("distance_score_from_", "", .x), .cols = everything()) |> # remove "distance_score_from_" from every column name
          distinct(ROI_nr, .keep_all = T) |> 
          column_to_rownames("ROI_nr") |>
          t()  

#plot heatmap inter-cell distances
hmdat_2f           = apply(hmdat_2, 2, as.numeric)  
rownames(hmdat_2f) = rownames(hmdat_2)
hm2 = ComplexHeatmap::Heatmap(hmdat_2f,
                              column_names_rot  = 90, 
                              height            = unit(5, "cm"),
                              column_names_gp   = grid::gpar(fontsize = 2),
                              column_order      = colorder_hm1,
                              na_col            = "grey90",
                              row_title         = "Distance data",
                              show_column_names = F) 

#check col order with hm1
colorder_hm2 = hm2 |> 
  ComplexHeatmap::column_order() |>
  as.vector()

#integrate and plot heatmaps
ht_list = hm1 %v% hm2

pdf(file = "1c.pdf", height = 10, width = 18)
ComplexHeatmap::draw(ht_list) 
dev.off()

##### 1f ####
#exclude failed stainings
cellsub     = cells[-c(5,13)] 
vesselsub   = vessels[-c(5,13)]

#calculate inter-cell distances
tmp31       = mclapply(seq_along(cellsub), function(x) st_nearest_feature(cellsub[[x]]$centroid, vesselsub[[x]]$CD31$geometry),mc.cleanup = T, mc.cores = 15)
cells31     = lapply(seq_along(cellsub), function(x) cbind(cellsub[[x]], tmp31[[x]]))

#integrate and subset data all samples
result = list()
rows   = list()

for (i in seq_along(cellsub)) {
  set.seed(23)
  
  a = cells31[[i]]
  b = vesselsub[[i]]
  r = which(a$class_bio == "T_cell")
  
  if (length(r)>2000) {
    r = sample(r, size = 2000)
  }
  
  s = sample(which(a$class_bio == "Cell"), size = length(r))
  rows[[i]] = c(r,s)
  c = a[c(r,s),]
  
  result[[i]] =   mclapply(1:length(c$centroid), function(x) st_distance(c$centroid[x], b$CD31$geometry[c$`tmp31[[x]]`][x]), mc.cleanup = T, mc.cores = 15) |> 
                  unlist()
  
  print("next")
  gc()
}

cells31sub = lapply(seq_along(cells31),function(x) cells31[[x]][rows[[x]],])
try        = lapply(seq_along(cells31sub), function(x) cbind(cells31sub[[x]],result[[x]]))

for (i in 1:length(try)) {
  colnames(try[[i]]) = c("id","CD3","class_bio","objectType",
                         "classification","isLocked","measurements","geometry",
                         "annotationsimple","centroid","tmp31", "result")
}

for (i in seq_along(try)) {
  try[[i]] = try[[i]][,c("id","CD3","class_bio","objectType","annotationsimple","tmp31", "result")]
}

try2 = data.table::rbindlist(try, use.names = T, id = "sample")

try3 = try2 |> 
       filter(result<600) |> 
       mutate(result = result*0.3441)

#test differences between data distributions
ks.test(try2[try2$class_bio == "Cell",]$result, try2[try2$class_bio == "T_cell",]$result)

#plot density distance to CD31+ object
pdf("1f.pdf", width = 5, height = 3)
ggplot() +
  geom_rect(data = data.frame(a=1),xmin = -Inf, xmax = 17.84183,   ymin = -2, ymax = 2,   fill = "red", alpha = .2) +
  geom_density(data = try3, aes(result, col = class_bio),linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = c("#e41a1c", "#377eb8"))
dev.off()

##### 1g ####
#filter samples that failed QC
cells2              = cells[-c(5,13)]

#test which cells are contained by which cuff
whichcuff           = mclapply(c(1:13), function(x) st_within(cells2[[x]]$centroid, nonconvexsub[[x]], sparse = T) |>  
                                 lapply(function(x) ifelse(length(x)>0,x,0)) |> 
                                 unlist(), mc.cores = 15, mc.cleanup = T)
cells3              = lapply(seq_along(cells2), function(x) cbind(cells2[[x]], whichcuff[[x]])) |> 
                      lapply(function(x) x |> dplyr::select(id, class_bio, centroid, `whichcuff[[x]]`))
cellscomb           = data.table::rbindlist(cells3, use.names = T, idcol = "samp")  
colnames(cellscomb) = c("samp","id", "class_bio", "centroid", "whichcuff")

data = cellscomb |>
       mutate(cuff = ifelse(whichcuff == 0, "No", "Yes")) |> 
       group_by(samp, class_bio, cuff) |> 
       summarise(sum = n()) |> 
       pivot_wider(names_from = cuff, values_from = sum) |> 
       mutate(percentage = Yes /(No+Yes)*100)

#calculate areas tumor stroma and t cell cuffs
area_cuffs      = sapply(nonconvexsub[c(1:13)], function(x) sum(st_area(x)))
area_all        = sapply(windows[c(1:4,6:12,14,15)], function(x) sum(st_area(x)))
ara_cuff_sample = data.frame(names = names(area_all), samp = 1:13, area_window = area_all-area_cuffs, cuffs = area_cuffs)
cons            = 1000000

#test and plot difference in T cell count adjusted for tissue area
ptest = data |> 
        filter(class_bio == "T_cell") |> 
        left_join(ara_cuff_sample) |> 
        mutate(Cuff = Yes/cuffs*cons, 
               Stroma = No/area_window*cons) |>
        ungroup() |> 
        select(Cuff, Stroma) |> 
        pivot_longer(cols = everything(), names_to = "Compartment", values_to = "Area adjusted T cell quantity") |> 
        rstatix::wilcox_test(`Area adjusted T cell quantity`~Compartment, paired = T) |> 
        rstatix::add_significance()

pdf("1g.pdf", height = 2.5, width = 3.5)
data |> 
  filter(class_bio == "T_cell") |> 
  left_join(ara_cuff_sample) |> 
  mutate(Cuff = Yes/cuffs*cons, 
         Stroma = No/area_window*cons) |>
  ungroup() |> 
  select(Cuff, Stroma) |> 
  pivot_longer(cols = everything(), names_to = "Compartment", values_to = "Area adjusted T cell quantity") |> 
  ggplot(aes(Compartment, `Area adjusted T cell quantity`, col = Compartment)) +
  geom_boxplot(size = 1) +  
  theme_classic() +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  ggpubr::stat_pvalue_manual(data = ptest, label = "p.signif", y.position = 60) +
  theme(legend.position = "none")
dev.off()

##### 1h ####
#subset cell counts
cell_counts = data[,c("ROI_nr","res","Study_ID","CD8_T_cells", "CD4_T_cells","B_cells", "Macrophages", "Tumor_cells", "number.of.All.cells")] |> 
              rename("number.of.All.cells" = "Total_cells") |>   
              left_join(metadata, by = c("ROI_nr", "Study_ID"))

samenstelling_TME = cell_counts |> 
                    dplyr::select(c(ROI_nr,Study_ID, res, cuff_score, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells)) |> 
                    dplyr::filter(cuff_score != "NA") |> 
                    mutate(T_cells    = CD8_T_cells + CD4_T_cells, 
                           Annotation = ifelse(cuff_score %in% c("+", "++", "+++"), "T_cells_in_PS", "T_cell_absent"),
                           Annotation = ifelse(Annotation == "T_cell_absent" & T_cells >= 3, "T_cells_in_TS", Annotation), 
                           Total_cells = as.numeric(Total_cells))|> 
                    dplyr::select(-cuff_score, -T_cells) |> 
                    pivot_longer(cols = !c("Study_ID", "res", "Annotation", "ROI_nr"))

#Test and plot differences in cell counts between ROI annotations
wilcox_1.D = samenstelling_TME |> 
             dplyr::group_by(name) |> 
             rstatix::wilcox_test(formula = value ~ Annotation,  
                                  paired  = F) |> 
             rstatix::adjust_pvalue(method = "fdr") |> 
             mutate()

pdf("1h.pdf", height = 6, width = 3.5)
ggplot(samenstelling_TME, aes(factor(Annotation, level = c("T_cells_in_PS", "T_cells_in_TS", "T_cell_absent")), log1p(value))) +
  geom_boxplot(aes(fill = Annotation)) +
  facet_wrap(~factor(name, level = c("Total_cells", "Tumor_cells", "Macrophages", "CD4_T_cells", "CD8_T_cells", "B_cells")), scales = "fixed", nrow = 3)+
  labs(x = "Region annotation", y = "log(value + 1)") + 
  theme_bw() +
  theme(axis.text.x     = element_text(angle = 60, hjust=1), 
        legend.position = "none") +
  scale_fill_manual(values = c("#99d8c9", "#00441b", "#238b45")) +
  stat_pvalue_manual(wilcox_1.D, 
                     label = "p.adj.signif", 
                     y.position = c(10,8,9,10,8,9,10,9,10,8,9,10,9,10), 
                     tip.length = 0,  
                     hide.ns = T, 
                     size = 3)+
  coord_cartesian(ylim = c(0,10.5)) +
  labs(x = "Log(count + 1)")
dev.off()

##### 1i ####
# create data frame containing relevant inter-inter cell distance z-scores
data_distance_test = data |> 
  left_join(samenstelling_TME[, c("ROI_nr", "Annotation")], by="ROI_nr") |> 
  filter(Annotation %in% c("T_cells_in_PS", "T_cells_in_TS")) |>
  dplyr::select(distance_score_from_CD8_T_cells_to_CD8_T_cells,
                distance_score_from_CD8_T_cells_to_CD4_T_cells,
                distance_score_from_CD4_T_cells_to_CD4_T_cells, 
                distance_score_from_B_cells_to_B_cells,
                distance_score_from_Macrophage_cells_to_Macrophage_cells,
                distance_score_from_Macrophage_cells_to_CD4_T_cells,
                distance_score_from_Macrophage_cells_to_B_cells,
                distance_score_from_Macrophage_cells_to_CD8_T_cells,
                distance_score_from_Tumor_cells_to_Tumor_cells,
                distance_score_from_B_cells_to_CD4_T_cells,
                distance_score_from_CD4_T_cells_to_Macrophage_cells,
                distance_score_from_CD8_T_cells_to_Macrophage_cells,
                distance_score_from_CD4_T_cells_to_CD8_T_cells,
                distance_score_from_B_cells_to_Macrophage_cells,
                distance_score_from_CD4_T_cells_to_B_cells,
                distance_score_from_B_cells_to_CD8_T_cells,
                distance_score_from_CD8_T_cells_to_B_cells,
                ROI_nr, 
                Annotation) |>
  pivot_longer(cols = c(distance_score_from_CD8_T_cells_to_CD8_T_cells, 
                        distance_score_from_CD4_T_cells_to_CD4_T_cells, 
                        distance_score_from_B_cells_to_B_cells, 
                        distance_score_from_Macrophage_cells_to_Macrophage_cells,
                        distance_score_from_Tumor_cells_to_Tumor_cells,
                        distance_score_from_B_cells_to_CD4_T_cells,
                        distance_score_from_CD8_T_cells_to_Macrophage_cells,
                        distance_score_from_CD4_T_cells_to_Macrophage_cells,
                        distance_score_from_CD4_T_cells_to_CD8_T_cells,
                        distance_score_from_B_cells_to_Macrophage_cells,
                        distance_score_from_CD4_T_cells_to_B_cells,
                        distance_score_from_CD8_T_cells_to_CD4_T_cells,
                        distance_score_from_B_cells_to_CD8_T_cells,
                        distance_score_from_CD8_T_cells_to_B_cells,        
                        distance_score_from_Macrophage_cells_to_CD4_T_cells,
                        distance_score_from_Macrophage_cells_to_B_cells,
                        distance_score_from_Macrophage_cells_to_CD8_T_cells)) |> 
  dplyr::filter(value != "NA") |>
  mutate(value = as.numeric(value)) |> 
  mutate(name = case_when(name == "distance_score_from_CD8_T_cells_to_CD8_T_cells"             ~ "CD8_T_cells-CD8_T_cells",
                          name == "distance_score_from_B_cells_to_B_cells"                     ~ "B_cell-B_cell",
                          name == "distance_score_from_CD4_T_cells_to_CD4_T_cells"             ~ "CD4_T_cells-CD4_T_cells",
                          name == "distance_score_from_Macrophage_cells_to_Macrophage_cells"   ~ "TAM-TAM",
                          name == "distance_score_from_Macrophage_cells_to_CD4_T_cells"        ~ "TAM-CD4_T_cells",
                          name == "distance_score_from_Macrophage_cells_to_B_cells"            ~ "TAM-B_cell",
                          name == "distance_score_from_Macrophage_cells_to_CD8_T_cells"        ~ "TAM-CD8_T_cells",
                          name == "distance_score_from_Tumor_cells_to_Tumor_cells"             ~ "Tumor_cell-tumor_cell",
                          name == "distance_score_from_B_cells_to_CD4_T_cells"                 ~ "B_cell-CD4_T_cells",
                          name == "distance_score_from_CD4_T_cells_to_Macrophage_cells"        ~ "CD4_T_cells-TAM",
                          name == "distance_score_from_CD8_T_cells_to_Macrophage_cells"        ~ "CD8_T_cells-TAM",
                          name == "distance_score_from_CD4_T_cells_to_CD8_T_cells"             ~ "CD4_T_cells-CD8_T_cells",
                          name == "distance_score_from_B_cells_to_Macrophage_cells"            ~ "B_cell-TAM",
                          name == "distance_score_from_CD4_T_cells_to_B_cells"                 ~ "CD4_T_cells-B_cell",
                          name == "distance_score_from_CD8_T_cells_to_CD4_T_cells"             ~ "CD8_T_cells-CD4_T_cells",
                          name == "distance_score_from_B_cells_to_CD8_T_cells"                 ~ "B_cell-CD8_T_cells",
                          name == "distance_score_from_CD8_T_cells_to_B_cells"                 ~ "CD8_T_cells-B_cell"),
         name_to     = sapply(strsplit(name, "-"),"[",2),
         name_from   = sapply(strsplit(name, "-"),"[",1), 
         order = case_when(name_to == "B_cell"     ~ 1,
                           name_to == "CD4_T_cell" ~ 2,
                           name_to == "CD8_T_cell" ~ 3,
                           name_to == "TAM"        ~ 4))

# test and plot differences of inter-inter cell distance z-scores between ROI annotations
wilcox_1.D = data_distance_test |> 
             dplyr::group_by(name) |> 
             rstatix::wilcox_test(formula     = value ~ Annotation, 
                                  paired = F) |> 
             rstatix::adjust_pvalue(method = "fdr") |> 
             rstatix::add_significance(p.col = "p.adj") 

a = ggplot(data_distance_test |> filter(name %in% wilcox_1.D[1:4,]$name),aes(reorder(name, order), value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA,aes(fill = Annotation)) +
    theme_bw() +
    theme(strip.text      = element_text(size=.1),
          axis.text.x     = element_text(angle = 60, hjust=1), 
          strip.text.x    = element_text(size = 6),
          axis.title.x    = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#238b45", "#99d8c9")) +
    labs(y = "z-score inter-cell distance", title = "B_cell") +
    coord_cartesian(ylim = c(-6,14)) +
    stat_pvalue_manual(wilcox_1.D[1:4,], 
                       x = "name",
                       label = "p.adj.signif", 
                       y.position = 13, tip.length = 0,  hide.ns = F, size = 3) 

b = ggplot(data_distance_test |> filter(name %in% wilcox_1.D[5:8,]$name),aes(reorder(name, order), value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA,aes(fill = Annotation)) +
    theme_bw() +
    theme(strip.text      = element_text(size=.1),
          axis.text.x     = element_text(angle = 60, hjust=1), 
          strip.text.x    = element_text(size = 6),
          axis.title.x    = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#238b45", "#99d8c9")) +
    labs(y = "z-score inter-cell distance", title = "CD4_T_cells") +
    coord_cartesian(ylim = c(-6,14)) +
    stat_pvalue_manual(wilcox_1.D[5:8,], 
                       x = "name",
                       label = "p.adj.signif", 
                       y.position = 13, tip.length = 0,  hide.ns = F, size = 3) 

c = ggplot(data_distance_test |> filter(name %in% wilcox_1.D[13:16,]$name),aes(reorder(name, order), value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA,aes(fill = Annotation)) +
    theme_bw() +
    theme(strip.text      = element_text(size=.1),
          axis.text.x     = element_text(angle = 60, hjust=1), 
          strip.text.x    = element_text(size = 6),
          axis.title.x    = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#238b45", "#99d8c9")) +
    labs(y = "z-score inter-cell distance", title = "TAM") +
    coord_cartesian(ylim = c(-6,14)) +
    stat_pvalue_manual(wilcox_1.D[13:16,], 
                       x = "name",
                       label = "p.adj.signif", 
                       y.position = 13, tip.length = 0,  hide.ns = F, size = 3) 

d = ggplot(data_distance_test |> filter(name %in% wilcox_1.D[9:12,]$name),aes(reorder(name, order), value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA,aes(fill = Annotation)) +
    theme_bw() +
    theme(strip.text      = element_text(size=.1),
          axis.text.x     = element_text(angle = 60, hjust=1), 
          strip.text.x    = element_text(size = 6),
          axis.title.x    = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("#238b45", "#99d8c9")) +
    labs(y = "z-score inter-cell distance", title = "CD8_T_cell") +
    coord_cartesian(ylim = c(-6,14)) +
    stat_pvalue_manual(wilcox_1.D[9:12,], 
                       x = "name",
                       label = "p.adj.signif", 
                       y.position = 13, hide.ns = F, size = 3)

pdf("1i.pdf", height = 7, width = 9)
(a | b) / (c | d) + patchwork::plot_layout(guides = "collect")
dev.off()

#### 2 ####
##### 2c ####
#test and plot differences between GTC scores for ROIs 
data_C6 = cell_counts |> 
          mutate(Total_cells = as.numeric(Total_cells)) |> 
          tidyr::pivot_longer(cols = !c(Study_ID, res, ROI_nr, cuff_score, GTC_score, IDH_mutation, WHO_2021, Sample_Sex, CNV_load, Resection_nr)) |> 
          dplyr::filter(GTC_score != "NA") |> 
          dplyr::mutate(GTC_score = forcats::fct_relevel(GTC_score, "-", "+", "++", "+++")) |> 
          dplyr::mutate(value = log1p(value)) |> 
          dplyr::filter(!name %in% c("T_cells"))   
        
wilcox_C6 = data_C6 |> 
            dplyr::group_by(name) |> 
            rstatix::wilcox_test(formula = value ~ GTC_score, 
                                 comparisons = .comparisons_GTC, 
                                 paired = F) |> 
            rstatix::adjust_pvalue(method = "fdr")

pdf("2c.pdf", width = 9, height = 3)
ggplot(data_C6, aes(x=GTC_score, y=value)) + 
  geom_boxplot(aes(fill = GTC_score)) +
  facet_grid(~factor(name, level = .reorder_celltypes2), scales = "fixed") +
  scale_fill_manual(values = .col_cuff_GTC) +
  theme_bw() +
  theme(legend.position="none") +
  coord_cartesian(ylim = c(0,12.5)) +
  ggpubr::stat_pvalue_manual(data = wilcox_C6, y.position = c(     
    9.3, 10, 10.7, 11.4,
    8.6, 9.3, 10, 10.7, 11.4,
    8.6, 9.3, 10, 10.7, 11.4,  
    8.6, 9.3, 10, 10.7, 11.4, 12.1,
    8.6, 9.3, 10, 8.6, 9.3, 10), tip.length = 0,  hide.ns = T, label = "{p.adj.signif}", size = 3) + #give y.position as many numbers as you have significant comparisons, max 3
  labs(x = "Gemistocytes", y = "Log(counts + 1)") 
dev.off()

##### 2d ####
#test and plot distribution differences of GTC-scores depending on ROI annotation
data_C4 = cell_counts |>
          mutate(T_cells = CD4_T_cells + CD8_T_cells) |> 
          dplyr::filter(GTC_score != "NA")

datac4test = data_C4 |> 
  mutate(T_cell_reg = case_when(T_cells > 3 & cuff_score == "-" ~ "T_cells in TS",
                                T_cells > 3 & cuff_score != "-" ~ "T_cells in PS",
                                .default = "T_cells absent")) 

chisq.test(as.matrix(table(datac4test$T_cell_reg,datac4test$GTC_score)))
chisq.test(as.matrix(table(datac4test$T_cell_reg,datac4test$GTC_score))[1:2,])
chisq.test(as.matrix(table(datac4test$T_cell_reg,datac4test$GTC_score))[2:3,])
chisq.test(as.matrix(table(datac4test$T_cell_reg,datac4test$GTC_score))[1:3,])

pdf("2d.pdf", height = 3.5, width = 4.5)
data_C4 |> 
  mutate(T_cell_reg = case_when(T_cells > 3 & cuff_score == "-" ~ "T_cells in TS",
                                T_cells > 3 & cuff_score != "-" ~ "T_cells in PS",
                                .default = "T_cells absent")) |>
  ggplot(aes(factor(T_cell_reg, levels = c("T_cells absent", "T_cells in TS", "T_cells in PS")), fill = GTC_score)) +
  geom_bar(position="fill") +
  scale_fill_manual("Gemistocyte", values = c(.col_cuff_GTC)) +
  theme_bw() +
  labs(x = "ROI annotation", y = "% Stamps")
dev.off()  

##### 2e ####
#test and plot distribution differences of GTC-scores depending on ROI cuff annotation
data_C4_table = table(data_C4$GTC_score, data_C4$cuff_score)  
chisq.test(data_C4_table[2:3,])

pdf("2e.pdf", height = 3, width = 4)
data_C4 |> 
  dplyr::mutate(GTC_score = factor(GTC_score, levels = .reorder_GTC)) |> 
  ggplot(aes(GTC_score, fill = cuff_score)) +
  geom_bar(position="fill") +
  scale_fill_manual("Cuffs", values = c(.col_cuff_GTC)) +
  theme_bw() +
  labs(x = "Gemistocytes", y = "% Stamps")
dev.off()

##### 2f ####
#DEG between bulk GTC high and low samples
bdt = DE_test |> 
      rownames_to_column("gene") |> 
      left_join(markers |> 
                  dplyr::filter(Celltype %in% c("Tcell", "mac", "mic")), by = "X") |> 
      mutate(Celltype = ifelse(padj < 0.05 &  log2FoldChange < -.5 | padj < 0.05 & log2FoldChange > .5, Celltype, NA),
             Celltypeorder = ifelse(is.na(Celltype), 1,2),
             Celltype = ifelse(is.na(Celltype), NA, Celltype), 
             Celltype = case_when(Celltype == "mic" ~ "TAM",
                                  Celltype == "mac" ~ "TAM",
                                  .default = Celltype), 
             Name  = ifelse(padj < 0.05 &  log2FoldChange < -.5 | padj < 0.05 & log2FoldChange > .5, X, NA),
             size  = ifelse(is.na(Celltype), 1,1.1)) |>
      arrange(Celltypeorder)

pdf("2f.pdf", height = 4, width = 7)
ggplot(bdt, aes(log2FoldChange, -log10(padj),label = Name, col = Celltype, size = size)) +
  geom_point() +
  scale_size(range = c(.5,1.5)) +
  theme_bw() + 
  coord_cartesian(xlim = c(-3.1,3.1)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  geom_vline(xintercept = c(-.5,.5), linetype = 2)  +
  scale_color_manual(values = c("#f46d43","#74add1"), na.value = "grey")
dev.off()

##### 2g ####
#Pathway enrichments for differentially expressed genes between GTC high and low bulk RNA measurements
pdf("2g.pdf", width = 8, height = 4)
depathways |> 
  mutate(z.score        = as.numeric(z.score), 
         X.log.p.value. = as.numeric(X.log.p.value.)) |>
  slice_max(X.log.p.value., n = 30) |>
  drop_na(z.score) |> 
  ggplot(aes(z.score, reorder(Ingenuity.Canonical.Pathways, z.score), fill = X.log.p.value.)) +
  geom_col() +
  theme_bw() +
  ylab("Pathways") +
  scale_fill_viridis_c(option = "plasma")
dev.off()

#### 3 ####
##### 3a ####
#umap of integrated glioma dataset
pdf("3a.pdf", height = 4, width = 6)
DimPlot(glioma.integrated, group.by = "bioidentsintegrated", pt.size = .0005,
        cols = c(
          "#ff5733",  # Vivid Red-Orange
          "#4caf50",  # Vivid Green
          "#795548",  # Vivid Brown
          "#f44336",  # Vivid Red
          "#9c27b0",  # Vivid Purple
          "#ffc300",  # Vivid Yellow
          "#e91e63",  # Vivid Pink
          "#607d8b",  # Vivid Gray
          "#ff9800",  # Vivid Orange
          "#2196f3"   # Vivid Blue
        ))
dev.off()

##### 3b ####
#umap of integrated tumor cells
pdf("3b.pdf", height = 4, width = 5)
DimPlot(glioma.tumor, group.by = "seurat_clusters",
        cols = c(
          "#ff6e54",  # Coral
          "#ffcc67",  # Yellow
          "#6eeb83",  # Mint Green
          "#51b2e5",  # Sky Blue
          "#e772d3",  # Bright Pink
          "#a25de6",  # Purple
          "#ffb65e",  # Orange
          "#7ddc7d",  # Lime Green
          "#4fa5c9",  # Blue
          "#ff8b88",  # Peach
          "#98e4b5",  # Light Green
          "#ef9a9a",  # Light Red
          "#bd93d8",  # Lavender
          "#ffa38b",  # Light Orange
          "#80deea",  # Light Blue
          "#ff9fff",  # Light Pink
          "#4b4b4b"   # Gray
        ))
dev.off()

##### 3c ####
#filter bulk GTC genes 
gemmark      = DE_test |>
               rownames_to_column("X") |> 
               dplyr::filter( padj < 0.001 & log2FoldChange > .5) |>
               pull(X)

#calculate enrichment of bulk GTC gene signature in tumor cell populations 
feat         = head(VariableFeatures(glioma.tumor, assay = "integrated"), 2000)
marktest     = data.frame(Celltype = "gemistocyte", gene = feat[feat %in% gemmark])
glioma.tumor = cellscore(glioma.tumor, marktest)
col_fun      = circlize::colorRamp2(c(-1, 0, 2.5), viridis(20)[c(1,5,20)])
tmp          = AverageExpression(object = glioma.tumor, group.by = "seurat_clusters", features = marktest$gene)

colann       = glioma.tumor@meta.data[,c("seurat_clusters", "bioidents")] |> 
               distinct(seurat_clusters, .keep_all = T) |> 
               arrange(seurat_clusters)

#plot heatmap of bulk GTC signature in tumor cell populations
column_ha = HeatmapAnnotation(annotation = colann$bioidents, 
                              col        = list(annotation = setNames(c("#ca0020", "#f4a582", "#92c5de", "#0571b0"), 
                                                                      unique(colann$bioidents))))

pdf("3c.pdf", height = 6, width = 8)
Heatmap(matrix          = t(scale(t(tmp[["SCT"]]))), 
        col             = col_fun,
        top_annotation  = column_ha,
        row_names_gp    = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 10))
dev.off()

##### 3d ####
#featureplot of bulk GTC and Astro-like enrichment scores in tumor cell populations
l = FeaturePlot(object = glioma.tumor,  features = "gemistocyte",      order = T, max.cutoff = .4) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

m = FeaturePlot(object = glioma.tumor,  features = "astro.like.tumor", order = T, max.cutoff = .8) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))

pdf("3d.pdf", height = 8, width = 6)
l+m
dev.off()

##### 3e ####
# barplot of differential reactive astrocyte marker gene expression bulk RNA seq data
pdf(file = "3e.pdf", height = 5, width = 12)
DE_test |>
  rownames_to_column("gene") |> 
  dplyr::filter(gene %in% reastdf$gene) |> 
  ggplot(aes(reorder(gene,log2FoldChange),log2FoldChange, fill = -log10(padj))) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

##### 3f ####
# dotplot of differential reactive astrocyte marker gene expression bulk RNA seq data
tmp2 = unique(reastdf$gene[reastdf$gene %in% VariableFeatures(glioma.tumor, assay = "integrated")])

glioma.tumor@meta.data = glioma.tumor@meta.data |>
                         mutate(bioidents3 = ifelse(gemistocyte > 0.2, "Gemistocyte tumor state", bioidents))

Idents(glioma.tumor) = "bioidents3"
dat = DotPlot(glioma.tumor, features = tmp2)
dat = dat$data

Order = c("Oligo-like tumor state","G1/S/G2/M","stem-like tumor state"
          ,"Astro-like tumor state","Gemistocyte tumor state")

Order2 = dat |> filter(id == "Gemistocyte tumor state") |> 
  arrange(desc(avg.exp.scaled)) |> pull(features.plot)

pdf("3f.pdf", width = 8, height = 2.5)
ggplot(dat, aes(factor(features.plot,level = Order2), factor(id, level = Order))) +
  geom_point(aes(col = avg.exp.scaled, size = pct.exp)) +
  theme_bw()+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  theme(legend.position = "top")
dev.off()

##### 3g ####
# filter NanoString spatial protein assay data for normalization and housekeeping markers and outliers
housekeepers_ns_protein = c("S6", "Histone H3", "GAPDH")
Negcontrol_ns_protein   = c("Rb IgG", "Ms IgG1", "Ms IgG2a")

dat.norm.neg = dat.norm.neg[!rownames(dat.norm.neg) %in% c(housekeepers_ns_protein, Negcontrol_ns_protein, low_targets, "PanCk", "HYB-POS"),
                            !colnames(dat.norm.neg) %in% c("2_11", "3_12", "6_07")]

#remove 2_11 because it is an outlier
annotation_vis = annotation_NS_protein[!annotation_NS_protein$Segment_tags %in% "Control" & !annotation_NS_protein$name %in% "2_11",] 
stopifnot(all(colnames(dat.norm.neg) %in% annotation_vis$name))

geman      = annotation_vis |> 
             dplyr::filter(GTC_score != "NA") |> 
             select(name, GTC_score2) |> 
             mutate(a = "a") |> 
             column_to_rownames("name")

#test and plot differential protein presence in NanoString ROIs
datdeseq = dat.norm.neg[,colnames(dat.norm.neg) %in% rownames(geman)]

geman = geman[colnames(datdeseq),]
stopifnot(all(colnames(datdeseq) == rownames(geman)))
dds    = DESeqDataSetFromMatrix(countData = round(datdeseq),
                                colData   = geman, 
                                design    = ~ GTC_score2)

sizeFactors(dds) = 1
dds              = DESeq(dds)
res              = results(dds, name = "GTC_score2_2High_vs_1Low") |> 
                   as.data.frame()

pdf("3g.pdf", height = 4.5, width = 5.5)
res |> 
  rownames_to_column("Protein") |> 
  mutate(label = ifelse(Protein %in% c("Fibronectin", "SMA", "OX40L", "PTEN", "CD44", "CD127"), Protein, NA)) |> 
  ggplot(aes(log2FoldChange, -log10(padj), label = label)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-.5,.5),    linetype = 2) +
  geom_label_repel()
dev.off()

#### 4 ####
##### 4b ####
#relate NanoString spatial RNA ROI traits to gene modules from WGCNA analysis
moduleColors = modulesspatial$moduleColor[match(colnames(datExpr), modulesspatial$gene)]
nGenes       = ncol(datExpr)
nSamples     = nrow(datExpr)

datTraits    = annotation_NS_RNA |> 
               filter(qc == "yes") |> 
               dplyr::select(name, GTC_score2, T_cell_count_area, Resection, ROI_tag) |> 
               mutate(T_cells_in_perivascular_space = ifelse(ROI_tag   =="T_cells_in_perivascular_space", 1,0), 
                      T_cells_in_tumor_stroma        = ifelse(ROI_tag   =="T_cells_in_Tumor_stroma",       1,0),
                      T_cells_absent                = ifelse(ROI_tag   =="T_cells_absent",                1,0), 
                      Resection                     = ifelse(Resection == "Primary",                      1,2),
                      T_cell_count_area             = as.numeric(T_cell_count_area)) |> 
               select(!ROI_tag) |> 
               column_to_rownames("name")

MEs               = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs               = MEs[,names(MEs) != "MEgrey"]
datTraits         = datTraits[match(paste(substr(rownames(MEs),4,4),substr(rownames(MEs),7,8), sep = "_"),rownames(datTraits)),]
moduleTraitCor0   = cor(MEs, datTraits, use = "p")
moduleTraitCor    = moduleTraitCor0[order(moduleTraitCor0[,3],decreasing = T),]
moduleTraitPvalue = data.frame(corPvalueStudent(moduleTraitCor, nSamples))

MEs = MEs[,match(rownames(moduleTraitCor),names(MEs))]

textMatrix = matrix(nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))

for (i in 1: nrow(moduleTraitCor)) {
  for (j in 1:ncol(moduleTraitCor)) {
    textMatrix[i,j] =   paste("R = ", signif(moduleTraitCor[i,j], 2), "\n", "p = ",
                              signif(moduleTraitPvalue[i,j], 1), sep = "")
  }  }

#plot gene modules relevant for spatial tissue phenotypes
pdf("4b.pdf", width = 5,height = 6)
par(mar = c(10, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor[c("MEblack", "MEturquoise", "MEtan", "MEpink", "MEyellow", "MEred", "MEgreen"),
                                       c("T_cells_in_tumor_stroma","T_cells_in_perivascular_space","T_cells_absent")],
               xLabels = colnames(moduleTraitCor)[c(3,4,5)],
               yLabels = paste("module", 1:7),
               ySymbols = rownames(moduleTraitCor)[c(1,2,5,6,8,10,12)],
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = .9,
               cex.lab.y = 1.2,
               cex.lab.x =1.2,
               zlim = c(-1,1))
dev.off()

##### 4c ####
#  plot correlation of gene scores for T cells in PS with gene scores for GTC quantity 
pdf("4c.pdf", height = 3.5, width = 5)
modulesspatial |> 
  mutate(colors = ifelse(moduleColor %in% c("tan", "yellow", "pink"), moduleColor, NA), 
         order  = ifelse(moduleColor %in% c("tan", "yellow", "pink"),2,1),
         colors = case_when(colors == "yellow" ~ "Module 4",
                            colors == "tan"    ~ "Module 3",
                            colors == "pink"   ~ "Module 5")) |>
  select(GS.Gemistocytes, GS.T.cells.in.perivascular.space, colors, order) |>
  arrange(order) |> 
  ggplot(aes(GS.Gemistocytes, GS.T.cells.in.perivascular.space, col = colors, size = order)) +
  geom_point() +
  theme_bw() +
  scale_size(range = c(.5,1)) +
  scale_color_manual(values =  c("#f46d43", "#43a2f4", "#f4d742"), na.value = "grey") +
  stat_cor(method = "pearson",cor.coef.name = "r", inherit.aes = F, aes(GS.Gemistocytes, GS.T.cells.in.perivascular.space)) +
  xlab("Gene score gemistocytes") + ylab("Gene score T cells in perivascular space")
dev.off()

##### 4d ####
#calculate enrichment scores for spatial gene modules in 3b in single nucleus RNA-sequencing data
gemistocyte_cells                      = rownames(glioma.tumor@meta.data)[glioma.tumor$gemistocyte>0.2]
glioma.integrated$cellid               = rownames(glioma.integrated@meta.data)
glioma.integrated$bioidentsintegrated2 = ifelse(glioma.integrated$cellid %in% gemistocyte_cells, "Gemistocyte tumor state", glioma.integrated$bioidentsintegrated)
Idents(glioma.integrated)              = "bioidentsintegrated2"

modulessub = modulesspatial |> 
             select(gene, moduleColor) |> 
             dplyr::rename("Celltype" = moduleColor) |> 
             filter(Celltype %in% c("tan", "turquoise", "black", "pink", "yellow", "green", "red")) |> 
             mutate(Celltype2 = case_when(Celltype == "green"     ~ 1, 
                                          Celltype == "red"       ~ 2,
                                          Celltype == "pink"      ~ 3,
                                          Celltype == "yellow"    ~ 4,
                                          Celltype == "black"     ~ 5,
                                          Celltype == "tan"       ~ 6,
                                          Celltype == "turquoise" ~ 7),
                    Celltype = case_when(Celltype == "green"      ~ "Module_6", 
                                         Celltype == "red"        ~ "Module_7",
                                         Celltype == "pink"       ~ "Module_4",
                                         Celltype == "yellow"     ~ "Module_5",
                                         Celltype == "black"      ~ "Module_1",
                                         Celltype == "tan"        ~ "Module_3",
                                         Celltype == "turquoise"  ~ "Module_2"))

glioma.integrated = cellscore(glioma.integrated, modulessub)

tmp = DotPlot(glioma.integrated, features = unique(modulessub$Celltype))$data |> 
      mutate(order   = case_when(features.plot == "Module_1" ~ 1, 
                                 features.plot == "Module_5" ~ 6,
                                 features.plot == "Module_6" ~ 2,
                                 features.plot == "Module_4" ~ 4,
                                 features.plot == "Module_7" ~ 3,
                                 features.plot == "Module_3" ~ 5,
                                 features.plot == "Module_2" ~ 7)) |> 
      mutate(order2  = case_when(id  == "Astro-like tumor state"   ~ 9, 
                                 id  == "Oligo-like tumor state"   ~ 8,
                                 id  == "Oligodendrocytes"         ~ 2,
                                 id  == "TAMs"                     ~ 12,
                                 id  == "G1/S/G2/M"                ~ 7,
                                 id  == "Astrocytes"               ~ 3,
                                 id  == "stem-like tumor state"    ~ 5, 
                                 id  == "undetermined"             ~ 11,
                                 id  == "Endothelial cells"        ~ 6,
                                 id  == "Neurons"                  ~ 1,
                                 id  == "T cells"                  ~ 4,
                                 id  == "Gemistocyte tumor state"  ~ 10))

#plot dotplot gene enrichment score of cell populations
pdf("4d.pdf", height = 3.5, width = 5)
ggplot(tmp, aes(reorder(features.plot, order), reorder(id, order2), size = pct.exp, col= avg.exp.scaled)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1), 
        axis.title  = element_blank()) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

##### 4e ####
#calculate macrophage/microglia enrichment scores in TAMs of snRNA-seq data
#calculate Module-3 enrichment score in TAMs of snRNA-seq data
progrmac =  markers |> 
            dplyr::filter(Celltype %in% c("mic", "mac")) |> 
            dplyr::select(!Ref)

markersspatial = modulesspatial |> 
                 dplyr::filter(moduleColor %in% c("tan", "turquoise", "yellow")) |> 
                 dplyr::mutate(Celltype = case_when(moduleColor == "tan"       ~ "Module_3",
                                                    moduleColor == "turquoise" ~ "Module_2",
                                                    moduleColor == "yellow"    ~ "Module_5")) |> 
                 dplyr::select(Celltype, gene)

tammarkers         = rbind(progrmac, markersspatial)
feat               = head(VariableFeatures(TAMs, assay = "integrated"), 2000)
marktest           = tammarkers[tammarkers$gene %in% feat,]

TAMs = cellscore(TAMs, marktest)

dat = data.frame(
      "Microglia"  = TAMs$mic,
      "Macrophage" = TAMs$mac,
      "Module_3"   = TAMs$Module_3) |>
      mutate(score = case_when(Microglia > 0 & Macrophage < 0 & Module_3 < 0 ~ "Microglia",
                               Microglia < 0 & Macrophage > 0 & Module_3 < 0 ~ "Macrophage",
                               Module_3  > 0 ~ "Module_3", 
                               .default =  "Undetermined")) |> 
      rownames_to_column('cell_id')

TAMs@meta.data = TAMs@meta.data |>
                 rownames_to_column("cell_id") |> 
                 left_join(dat, by = "cell_id") |> 
                 column_to_rownames("cell_id")

TAMs@meta.data          = TAMs@meta.data |>
                          mutate(trysum = Macrophage - Microglia)
TAMs$trysum2            = ifelse(TAMs$trysum >0.5, 0.5, TAMs$trysum)
TAMs$trysum3            = ifelse(TAMs$trysum2 < -0.5, -0.5, TAMs$trysum2)
TAMs@meta.data$Module_3 = TAMs$Module_3.x

#Plot UMAP representations of Macrophage/microglia scale and Module 3 score
a = FeaturePlot(TAMs, features = "trysum3",  pt.size = .2) +
    scale_colour_gradientn(colours = c("#b2182b","#f4a582","#fddbc7","#f7f7f7", "#92c5de","#2166ac","#053061")) +
    ggtitle("Macrophage-microglia scale")

b = FeaturePlot(TAMs, features = "Module_3",pt.size = .3, order = T) +
    scale_colour_gradientn(colours = c("grey90","grey90","grey90","#92c5de","#2166ac","#053061","#053061")) +
    ggtitle("Module_3 score")

pdf("4e.pdf", height = 6, width = 5)
a/b
dev.off()

##### 4f ####
#violin plot of enrichment scores for MDMs and microglia scores in Module-3 positive TAMs
pdf("4f.pdf", height = 2, width = 3)
TAMs@meta.data |> 
  dplyr::filter(score == "Module_3") |> 
  dplyr::select(Microglia, Macrophage) |> 
  pivot_longer(cols = everything()) |> 
  ggplot(aes(name, value, col = name)) +
  geom_violin(linewidth = 1) +
  theme_classic()+
  scale_color_manual(values = c("sienna3", "steelblue3")) +
  theme(legend.position = "none") +
  xlab("Score")
dev.off()

##### 4g ####
#plot showing which Module 3 genes are enriched in Module 3 positive TAMs
module3genes                  = markersspatial[markersspatial$Celltype == "Module_3",]$gene
TAMs@meta.data$module_3_score = ifelse(TAMs@meta.data$Module_3.x > 0.0, "True", "False")
Idents(TAMs)                  = "module_3_score"
avgtams                       = AverageExpression(TAMs, assays = "SCT", slot = "data")

pdf("4g.pdf", height = 4, width = 6)
as.data.frame(log1p(avgtams$SCT)) |> 
  rownames_to_column("gene") |> 
  mutate(Module_3 = gene %in% module3genes, 
         label    = ifelse(gene %in% c("CD83", "CCL3", "CCL4", "CCL4L2", "CCL3L3", "SPP1", "PDK4", "FOS", "SRGN", "SAT1", "IL1B", "OLR1"), gene, "")) |>
  arrange(Module_3) |> 
  filter(gene != "MALAT1") |> 
  ggplot(aes(True, False, col = Module_3, label = label)) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("grey", "sienna3")) +
  geom_text_repel(max.overlaps = 30, inherit.aes = F, aes(True, False, label = label), size = 3) +
  xlab("Module_3 positive") + ylab("Module_3 negative")
dev.off()  

##### 4h ####
#dotplot showing expression of reactive microglia markers in TAMs
pdf("4h.pdf", height = 4, width = 5.5)
DotPlot(glioma.integrated, features = c("SPP1", "TNF", "IL1A","IL1B", "C1QA","C1QB")) +
  scale_colour_gradientn(colours = c("grey90","grey90","#92c5de","#2166ac","#053061"))  +
  theme(axis.text.x  = element_text(angle = 60, hjust=1))
dev.off()

##### 4i ####
#cellchat analysis estimation of SPP1-CD44 pathway communication in snRNA-seq data
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")
pdf(file = "/mnt/data/users/levi/gemistocytes/plots/SPP1signaling.pdf", height = 4,width = 6)
netAnalysis_signalingRole_network(cellChat, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)
dev.off()

#### 5 ####
##### 5a ####
# calculate average expression of snRNA-seq determined markers for GTCs and Module 3 TAMs in bulk and spatial data
M3TAMs = rownames(TAMs@meta.data)[TAMs@meta.data$module_3_score == "True"]
glioma.integrated@meta.data$bioidentsintegrated3 = ifelse(rownames(glioma.integrated@meta.data) %in% M3TAMs, "Module_3_TAM", glioma.integrated@meta.data$bioidentsintegrated2)

glioma.integrated = PrepSCTFindMarkers(glioma.integrated)

gemistomarkers    = FindMarkers(glioma.integrated, 
                                ident.1 = "Gemistocyte tumor state", 
                                assay = "SCT",
                                group.by = "bioidentsintegrated3")
M3TAMmarkers      = FindMarkers(glioma.integrated, 
                                ident.1 = "Module_3_TAM", 
                                assay = "SCT",
                                group.by = "bioidentsintegrated3")

gemistomarkers_filtered = gemistomarkers |> 
                          filter(avg_log2FC>1,                                     
                                 p_val_adj<0.001) |> 
                          rownames_to_column("gene") |> 
                          mutate(celltype = "GTC")

M3TAMmarkers_filtered   = M3TAMmarkers |> 
                          filter(avg_log2FC>1, 
                                 p_val_adj<0.001) |> 
                          rownames_to_column("gene")|> 
                          mutate(celltype = "M3_microglia")

lis = sapply(strsplit(rownames(expressbulk), "_"), "[",2)

#plot correlation of GTC and module 3 TAM expression signatures in bulk and spatial expression data
a = data.frame("Module_3_microglia"   = expressbulk[lis %in% M3TAMmarkers_filtered$gene,]            |> apply(2, mean), 
               "Gemistocytes"         = expressbulk[lis %in% gemistomarkers_filtered$gene,]          |> apply(2, mean)) |>
    ggplot(aes(Gemistocytes, Module_3_microglia)) +
    geom_point() +
    theme_classic() +
    stat_cor(method = "spearman", cor.coef.name = "rho") +
    geom_smooth(method = "lm") + 
    ggtitle("Bulk")

b = data.frame("Module_3_microglia"   = dat_ns_RNA[rownames(dat_ns_RNA) %in% M3TAMmarkers_filtered$gene,]   |> log1p() |> apply(2, mean), 
               "Gemistocytes"         = dat_ns_RNA[rownames(dat_ns_RNA) %in% gemistomarkers_filtered$gene,] |> log1p() |> apply(2, mean)) |> 
    ggplot(aes(Gemistocytes, Module_3_microglia)) +
    geom_point() +
    theme_classic() +
    stat_cor(method = "spearman", cor.coef.name = "rho") +
    geom_smooth(method = "lm") + 
    ggtitle("Spatial")

pdf("5a.pdf", width = 7, height = 3.5)
a + b
dev.off()

##### 5d ####
#bar plot of cell proportions in whole slide validation scans
tmp4        = lapply(cd44data, function(x) as.data.frame(table(x$marks)))
names(tmp4) = names(cd44data)

a = data.table::rbindlist(tmp4, use.names = T, idcol = "Names_data")  |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |>
  filter(`Gemistocyte annotation` == "Gemistocyte high") |>
  mutate(Var1 = case_when(Var1 == "Cell" ~ "1Negative",
                          Var1 == "T_cell" ~ "CD3+",
                          Var1 == "Gemistocyte" ~ "CD44+")) |> 
  ggplot(aes(as.factor(Names_data), Freq, fill = Var1)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#fddbc7", "#d6604d", "#67001f")) +
  theme_minimal() 

b = data.table::rbindlist(tmp4, use.names = T, idcol = "Names_data")  |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |>
  filter(gemistocyte == "Gemistocyte low") |> 
  mutate(Var1 = case_when(Var1 == "Cell" ~ "1Negative",
                          Var1 == "T_cell" ~ "CD3+",
                          Var1 == "Gemistocyte" ~ "CD44+")) |> 
  ggplot(aes(as.factor(Names_data), Freq, fill = Var1)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#fddbc7", "#d6604d", "#67001f")) +
  theme_minimal() 

tmp5 = lapply(cd68data, function(x) as.data.frame(table(x$marks)))
names(tmp5) = datacd44$Panel_1

f = data.table::rbindlist(tmp5, use.names = T, idcol = "Names_data")  |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |>
  filter(gemistocyte == "Gemistocyte high") |> 
  mutate(Var1 = case_when(Var1 == "Cell" ~ "1Negative",
                          Var1 == "T_cell" ~ "CD3+",
                          Var1 == "TAM" ~ "CD68+")) |> 
  ggplot(aes(as.factor(Names_data), Freq, fill = Var1)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#e0e0e0", "#878787", "#1a1a1a")) +
  theme_minimal() 

g = data.table::rbindlist(tmp5, use.names = T, idcol = "Names_data")  |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |>
  filter(gemistocyte == "Gemistocyte low") |> 
  mutate(Var1 = case_when(Var1 == "Cell" ~ "1Negative",
                          Var1 == "T_cell" ~ "CD3+",
                          Var1 == "TAM" ~ "CD68+")) |> 
  ggplot(aes(as.factor(Names_data), Freq, fill = Var1)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#e0e0e0", "#878787", "#1a1a1a")) +
  theme_minimal() 

pdf("5d1.pdf", height = 3, width = 8)
(a + b) + plot_layout(widths = c(2,1.5), guides = "collect")
dev.off()

pdf("5d2", height = 3, width = 8)
(f + g) + plot_layout(widths = c(2,1.5), guides = "collect")
dev.off()

##### 5e ####
# Test and plot CRYAB+ cell quantities in GTChigh and low tumor samples
pval = CRYAB_data |> 
  select(Study_ID2,Area.µm.2, Gemistocyte, FITC..Cytoplasm..Median) |> 
  filter(FITC..Cytoplasm..Median > 200) |>
  group_by(Study_ID2,Area.µm.2, Gemistocyte) |> 
  summarize(sum = n()) |>
  ungroup() |> 
  mutate(rat = sum/Area.µm.2) |> 
  wilcox_test(rat~Gemistocyte) |> 
  add_y_position() |> 
  add_significance()

pdf("5e.pdf", height = 3, width = 4.5)
CRYAB_data |> 
  select(Study_ID2,Area.µm.2, Gemistocyte, FITC..Cytoplasm..Median) |> 
  filter(FITC..Cytoplasm..Median > 200) |>
  group_by(Study_ID2,Area.µm.2, Gemistocyte) |> 
  summarize(sum = n()) |>
  ungroup() |> 
  mutate(rat = sum/Area.µm.2) |> 
  ggplot(aes(Gemistocyte, rat, colour =  Gemistocyte)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() + 
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  stat_pvalue_manual(pval, y.position = 0.006, label = "p.signif")
dev.off()

##### 5f ####
#Nearest neighbor analysis showing z scores between Tcells, TAMs and GTCs
tmp  = lapply(cd44data, function(x) funtest(object = x, cores = 20))
tmp3 = lapply(cd68data, function(x) funtest(object = x, cores = 25))

names(tmp)  = names(cd44data)
names(tmp3) = names(cd68data)

tmp2 = data.table::rbindlist(tmp, use.names = T, idcol = "Names_data") |> 
  group_by(Panel_2, marks_to) |> 
  summarise(zmean = median(z)) |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |> 
  mutate(marks_from = sapply(strsplit(marks_to, "-"), "[", 1))|> 
  filter(!grepl("Cell",marks_to))

pvalc = tmp2 |> 
  group_by(marks_to) |> 
  rstatix::wilcox_test(zmean~gemistocyte) |> 
  rstatix::adjust_pvalue() |> 
  rstatix::add_significance() |> 
  rstatix::add_xy_position(x = "marks_to")

c = ggplot(tmp2, aes(marks_to, zmean, col = gemistocyte))  +
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~marks_from, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  stat_pvalue_manual(pvalc, label = "p.adj.signif", xmin = "marks_to", xmax = NULL, y = 1.5)

tmp4        = data.table::rbindlist(tmp3, use.names = T, idcol = "Names_data") |> 
  group_by(Panel_1, marks_to) |> 
  summarise(zmean = median(z)) |> 
  left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")], by = c("Names_data")) |> 
  mutate(marks_from = sapply(strsplit(marks_to, "-"), "[", 1)) |> 
  filter(!grepl("Cell",marks_to))

pvald = tmp4 |> 
  group_by(marks_to) |> 
  rstatix::wilcox_test(zmean~`Gemistocyte annotation`) |> 
  rstatix::adjust_pvalue() |> 
  rstatix::add_significance() |> 
  rstatix::add_xy_position(x = "marks_to")

d = ggplot(tmp4, aes(marks_to, zmean, col = `Gemistocyte annotation`))  +
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~marks_from, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  coord_cartesian(ylim = c(-2,2)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  stat_pvalue_manual(pvald, label = "p.adj.signif", xmin = "marks_to", xmax = NULL, y = 1.5)

pdf("5f.pdf", width = 8, height = 4)
c+d + plot_layout(guides = "collect")
dev.off()

##### 5h,i,j ####
#devide cells from whole slide multiplex image scans in 20 bins depending on CD44 density.
tmp = lapply(seq_along(s), function(i) {data.frame(x = s[[i]][s[[i]]$mark == "T_cell",]$x,
                                                   y = s[[i]][s[[i]]$mark == "T_cell",]$y,
                                                   v = s[[i]][s[[i]]$mark == "T_cell",]$v,
                                                   cluster =  dbscan(s[[i]][s[[i]]$mark == "T_cell",c("x","y")], 
                                                                     MinPts = 3, 
                                                                     scale  = FALSE, 
                                                                     eps    = 100,
                                                                     method = c("hybrid", "raw", "dist"))$cluster)})

binrange = data.table::rbindlist(s, use.names = T, idcol = "names") |> 
  drop_na(v) |> 
  mutate(bin = ntile(v, n= 20)) |> 
  group_by(bin) |> 
  dplyr::summarize(min = min(v), 
                   max = max(v))

tmp2 = rbindlist(tmp, use.names = T, idcol = "name") |> 
  drop_na(v) |> 
  group_by(name,cluster) |> 
  dplyr::summarize(mean = mean(v)) |> filter(cluster != 0) 

tmp3 = numeric()

for(i in 1:nrow(tmp2)){
  tmp3[i] = binrange[which(tmp2$mean[i]>binrange$min & tmp2$mean[i]<binrange$max),]$bin
}

tmp2$bin = tmp3

for (i in 1:length(s)) {
  colnames(s[[i]]) = c("x", "y", "mark", "v")
}

for (i in 1:length(density_try)) {
  colnames(density_try[[i]]) = c("geometry", "CD44_density")
}

#plot association between CD44 density and immune cell presence
bins = 20
aa = data.table::rbindlist(s, use.names = T, idcol = "names") |> 
    drop_na(v) |> 
    mutate(bin = ntile(v, n= bins)) |>
    left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                dplyr::rename("names" = Names_data), by = "names") |> 
    ggplot(aes(bin, fill = as.factor(names))) +
    geom_bar(col = "black") +
    scale_fill_manual(values = c(c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d"),
                                 c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))) +
    theme_bw() +
    labs(x = "CD44 density bins", y = "Cells", title = "Cell proportions per bin slide 1") +
    theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) 

ab = data.table::rbindlist(s, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin = ntile(v, n= bins)) |>
     group_by(bin) |> 
     dplyr::summarize(mean = mean(v)) |> 
     ggplot(aes(bin, mean)) + 
     geom_point() +
     geom_line(alpha = .2) +
     labs(x = "CD44 density bins", y = "CD44 density") +
     theme_bw() +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

ac = data.table::rbindlist(s, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin = ntile(v, n= bins), 
            mark2 = ifelse(mark == "T_cell", "T_cell", "Cell")) |>
     group_by(names, bin, mark2) |> 
     dplyr::summarize(sum = n()) |>
     ungroup() |> 
     pivot_wider(names_from = mark2, values_from = sum, values_fill = 0) |>
     mutate(fraction = T_cell/(Cell+T_cell)*100, 
            names = as.factor(names), 
            bin = as.numeric(bin)) |>
     left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                 dplyr::rename("names" = Names_data), by = "names")|>    
     filter(`Gemistocyte annotation` == "High") |> 
     mutate(`Gemistocyte high` = as.factor(names)) |> 
     ggplot(aes(bin,fraction)) +  
     geom_line(aes(group = `Gemistocyte high`)) +
     geom_point(aes(fill = `Gemistocyte high`), color='black', shape=21, size = 2) +
     theme_bw() +
     scale_fill_manual(values = c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d")) +
     labs(x = "CD44 density bins", y = "T cell fraction")+
     coord_cartesian(ylim = c(0,17))  +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) 

ad = data.table::rbindlist(s, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin = ntile(v, n= bins), 
            mark2 = ifelse(mark == "T_cell", "T_cell", "Cell")) |>
     group_by(names, bin, mark2) |> 
     dplyr::summarize(sum = n()) |>
     ungroup() |> 
     pivot_wider(names_from = mark2, values_from = sum, values_fill = 0) |> 
     mutate(fraction = T_cell/(Cell+T_cell)*100, 
            names = as.factor(names), 
            bin = as.numeric(bin)) |>
     left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                 dplyr::rename("names" = Names_data), by = "names") |> 
     filter(`Gemistocyte annotation` == "Low") |>  
     mutate(`Gemistocyte low` = as.factor(names)) |> 
     ggplot(aes(bin,fraction)) +  
     geom_line(aes(group = `Gemistocyte low`)) +
     geom_point(aes(fill = `Gemistocyte low`), color='black', shape=21, size = 2) +
     theme_bw() +
     scale_fill_manual(values = c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c")) +
     labs(x = "CD44 density bins", y = "T cell fraction") +
     coord_cartesian(ylim = c(0,17))

ba = data.table::rbindlist(z, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin = ntile(v, n= bins)) |>
     left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                 dplyr::rename("names" = Names_data), by = "names")|> 
     ggplot(aes(bin, fill = as.factor(names))) +
     geom_bar(col = "black") +
     scale_fill_manual(values = c(c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d"),
                                  c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))) +
     theme_bw() +
     labs(x = "CD44 density bins", y = "Cells", title = "Cell proportions per bin slide 2") +
     theme(legend.position = "none",
           axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())

bb = data.table::rbindlist(z, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin = ntile(v, n= bins)) |>
     group_by(bin) |> 
     dplyr::summarize(mean = mean(v)) |> 
     ggplot(aes(bin, mean)) + 
     geom_point() +
     geom_line(alpha = .2) +
     labs(x = "CD44 density bins", y = "CD44 density") +
     theme_bw() +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) 

bc = data.table::rbindlist(z, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin   = ntile(v, n= bins),
            mark2 = ifelse(mark == "TAM", "TAM", "Cell")) |>
     group_by(names, bin, mark2) |> 
     dplyr::summarize(sum = n()) |>
     ungroup() |> 
     pivot_wider(names_from = mark2, values_from = sum, values_fill = 0) |> 
     mutate(fraction = TAM/(Cell+TAM)*100, 
            names = as.factor(names), 
            bin = as.numeric(bin)) |>
     left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                 dplyr::rename("names" = Names_data), by = "names")|> 
     filter(`Gemistocyte annotation` == "High") |>
     mutate(`Gemistocyte high` = as.factor(names)) |> 
     ggplot(aes(bin,fraction)) +  
     geom_line(aes(group = `Gemistocyte high`)) +
     geom_point(aes(fill = `Gemistocyte high`), color='black', shape=21, size = 2) +
     theme_bw() +
     scale_fill_manual(values = c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#99000d")) +
     labs(x = "CD44 density bins", y = "TAM fraction") +
     theme(axis.title.x=element_blank(),
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank()) 

bd = data.table::rbindlist(z, use.names = T, idcol = "names") |> 
     drop_na(v) |> 
     mutate(bin   = ntile(v, n= bins),
            mark2 = ifelse(mark == "TAM", "TAM", "Cell")) |>
     group_by(names, bin, mark2) |> 
     dplyr::summarize(sum = n()) |>
     ungroup() |> 
     pivot_wider(names_from = mark2, values_from = sum, values_fill = 0) |> 
     mutate(fraction = TAM/(Cell+TAM)*100, 
            names = as.factor(names), 
            bin = as.numeric(bin)) |>
     left_join(annotation_CD44_CD68[,c("Names_data", "Gemistocyte annotation")] |>  
                 dplyr::rename("names" = Names_data), by = "names") |> 
     filter(`Gemistocyte annotation` == "Low") |>
     mutate(`Gemistocyte low` = as.factor(names)) |> 
     ggplot(aes(bin,fraction)) +  
     geom_line(aes(group = `Gemistocyte low`)) +
     geom_point(aes(fill = `Gemistocyte low`), color='black', shape=21, size = 2) +
     theme_bw() +
     scale_fill_manual(values = c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c")) +
     labs(x = "CD44 density bins", y = "TAM fraction") +
     coord_cartesian(ylim = c(0,55))

pdf("5hij.pdf", height = 9.5, width = 7)
((aa/ab/ac/ad) + plot_layout(heights = c(1.5,1.5,2.5,2.5)) | (ba/bb/bc/bd) + plot_layout(heights = c(1.5,1.5,2.5,2.5))) + plot_layout(guides = "collect") & theme(legend.position="bottom")
dev.off()

#### Supplementary figures ####
#### S1 ####
##### s1e ####
# Create cell count matrix
counts_per_stamp = cell_counts[,c(3,1,4:9)] |> 
                   dplyr::mutate(ptid        = substr(Study_ID, 1, nchar(Study_ID) -3), 
                                 Total_cells = as.numeric(Total_cells)) |> 
                   tidyr::pivot_longer(cols = !c(Study_ID, ROI_nr, ptid)) |> 
                   dplyr::left_join(metadata[,c("Study_ID", "Initial_Recurrent")] |> distinct(.keep_all = T), by = "Study_ID") |> 
                   dplyr::mutate(Initial_Recurrent = as.factor(Initial_Recurrent)) |> 
                   mutate(name = case_when(name == "Macrophages"       ~ "TAMs", 
                                           name == "Total_cells"       ~ "Total", 
                                           name == "Tumor_cells"       ~ "Tumor_cells", 
                                           .default = name))

#plot cell quantities per resection
pdf("s1e.pdf", height = 6, width = 12)
ggplot(counts_per_stamp, aes(x = Study_ID, y = log1p(value), fill = ptid)) +
  geom_boxplot(size = .2, outlier.size = .2) +
  scale_fill_manual(values = ..col_ptid, aes(alpha = 0.5)) +
  facet_wrap(~factor(name, level = .reorder_celltypes), ncol = 1, as.table = T, ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position="none", 
        strip.text = element_text(size=10)) +
  labs(y = "log(Count + 1)", x = "Sample ID") +
  geom_vline(xintercept=c(2.5,4.5,5.5,6.5,8.5,
                          9.5,11.5,13.5,15.5,17.5,
                          19.5,21.5,23.5,25.5,27.5,
                          28.5,31.5,33.5,34.5,35.5,
                          37.5,39.5,41.5,43.5,45.5,
                          48.5,50.5,51.5,53.5,55.5,
                          57.5,59.5,61.5,63.5,66.5,
                          69.5,70.5,73.5), linetype = "dotted", size = .2)
dev.off()

##### s1f ####
# Plot cell quantities averaged per sample
pdf("s1f.pdf", height = 5, width = 3)
cell_counts |>   
  select(Study_ID, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells) |> 
  mutate(Total_cells = as.numeric(Total_cells),
         CD8_T_cells = as.numeric(CD8_T_cells),
         CD4_T_cells = as.numeric(CD4_T_cells),
         B_cells = as.numeric(B_cells),
         Macrophages = as.numeric(Macrophages),
         Tumor_cells = as.numeric(Tumor_cells)) |>
  drop_na(Total_cells) |>
  group_by(Study_ID) |> 
  filter(n()>3) |> 
  filter(!Study_ID %in% notum) |> 
  group_by(Study_ID) |> 
  summarize(across(everything(), median)) |> 
  mutate_if(is.numeric, round, 0) |> 
  ungroup() |> 
  pivot_longer(cols = !Study_ID) |>
  mutate(value = log1p(value)) |> 
  mutate(name = case_when(name == "CD8_T_cells"       ~ "CD8_T_cells", 
                          name == "CD4_T_cells"       ~ "CD4_T_cells", 
                          name == "Macrophages"       ~ "TAMs", 
                          name == "Total_cells"       ~ "Total", 
                          name == "B_cells"           ~ "B_cells",
                          name == "Tumor_cells"       ~ "Tumor_cells")) |> 
  ggplot(aes(x = factor(name, level = .reorder_celltypes3), y = value)) +
  geom_boxplot(aes(color = name), outlier.shape = NA) +
  geom_quasirandom(varwidth = T, size = .6) +
  theme_bw() +
  scale_color_manual(values = .col_celltype) +
  labs(y = "log(count + 1)", x = "Cell type", color = "Cell type", fill = "Cell type") +
  theme(legend.position = "none") +
  coord_flip()
dev.off()

##### s1g ####
#plot cell quantities split for tumor high and low ROIs
data_3.2 = cell_counts |>   
  select(Study_ID, WHO_2021, Resection_nr, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells) |> 
  mutate(Total_cells = as.numeric(Total_cells),
         CD8_T_cells = as.numeric(CD8_T_cells),
         CD4_T_cells = as.numeric(CD4_T_cells),
         B_cells = as.numeric(B_cells),
         Macrophages = as.numeric(Macrophages),
         Tumor_cells = as.numeric(Tumor_cells), 
         T_cells     = CD8_T_cells+CD4_T_cells,
         class       = ifelse(median(Tumor_cells) > Tumor_cells, "Low", "High")) |>
  drop_na(Total_cells) |>
  group_by(Study_ID, WHO_2021, Resection_nr, class) |> 
  filter(n()>3) |> 
  filter(!Study_ID %in% notum) |> 
  summarize(across(everything(), median)) |> 
  mutate_if(is.numeric, round, 0) |> 
  ungroup() |> 
  pivot_longer(cols = !c(Study_ID, WHO_2021, Resection_nr, class)) |>
  mutate(value = log1p(value)) |> 
  mutate(name = case_when(name == "T_cells"           ~ "T_cells", 
                          name == "CD8_T_cells"       ~ "CD8_T_cells", 
                          name == "CD4_T_cells"       ~ "CD4_T_cells", 
                          name == "Macrophages"       ~ "TAMs", 
                          name == "Total_cells"       ~ "Total", 
                          name == "B_cells"           ~ "B_cells",
                          name == "Tumor_cells"       ~ "Tumor_cells"))  |> 
           dplyr::filter(!name %in% c("Total", "Tumor_cells"))

wilcox_3.2 = data_3.2 |> 
             dplyr::group_by(name) |> 
             rstatix::wilcox_test(formula = value ~ class, 
                                  comparisons = c("High", "Low"), 
                                  paired = F) |> 
             rstatix::adjust_pvalue(method = "fdr")

a = ggplot(data_3.2, aes(x = class, y = value)) +
    geom_boxplot(aes(fill= class), size = .2) +
    facet_grid(~factor(name, level = .reorder_celltypes4), as.table = T, scales="free") +
    scale_fill_manual(values = c("High" = "#238b45","Low" = "#99d8c9")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") +
    coord_cartesian(ylim = c(0,7)) +
    labs(y = "Log(count + 1)", x = "Tumor class") +
    ggpubr::stat_pvalue_manual(data = wilcox_3.2, y.position = c(6.5), hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}", size = 3.5) 

##### s1h ####
#paired analysis between cell quantities of primary and last recurrence resection cell quantities.
#cell quantities are averaged per resection
data_3.1 = cell_counts |>   
  select(Study_ID, Initial_Recurrent, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells) |> 
  mutate(Total_cells = as.numeric(Total_cells),
         CD8_T_cells = as.numeric(CD8_T_cells),
         CD4_T_cells = as.numeric(CD4_T_cells),
         T_cells      = CD8_T_cells + CD4_T_cells,
         B_cells = as.numeric(B_cells),
         Macrophages = as.numeric(Macrophages),
         Tumor_cells = as.numeric(Tumor_cells)) |>
  drop_na(Total_cells) |>
  group_by(Study_ID, Initial_Recurrent) |> 
  filter(n()>3) |> 
  filter(!Study_ID %in% notum) |> 
  group_by(Study_ID, Initial_Recurrent) |> 
  summarize(across(everything(), median)) |> 
  mutate_if(is.numeric, round, 0) |> 
  ungroup() |> 
  pivot_longer(cols = !c(Study_ID, Initial_Recurrent)) |>
  mutate(value = log1p(value)) |> 
  mutate(name = case_when(name == "CD8_T_cells"       ~ "CD8_T_cells", 
                          name == "CD4_T_cells"       ~ "CD4_T_cells", 
                          name == "Macrophages"       ~ "TAMs", 
                          name == "T_cells"           ~ "T_cells",
                          name == "Total_cells"       ~ "Total", 
                          name == "B_cells"           ~ "B_cells",
                          name == "Tumor_cells"       ~ "Tumor_cells")) |> 
  filter(Initial_Recurrent != "NA")

wilcox_3.1 = data_3.1 |> 
             dplyr::group_by(name) |> 
             rstatix::wilcox_test(formula = value ~ Initial_Recurrent, 
                                  comparisons = c("I","R"), 
                                  paired = F) |> 
             rstatix::adjust_pvalue(method = "fdr") 

b = ggplot(data_3.1, aes(x = Initial_Recurrent, y = value)) +
    geom_boxplot(aes(fill = Initial_Recurrent), size = .2) +
    facet_grid(~factor(name, level = .reorder_celltypes4), as.table = T, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") +
    coord_cartesian(ylim = c(0,8)) +
    scale_fill_manual(values = c("I" = "#fee8c8", "R" = "#99d8c9")) +
    labs(y = "Log(count + 1)", x = "Sample Type") +
    scale_x_discrete(labels = c("I" = "Initial", "R" = "Recurrent")) +
    ggpubr::stat_pvalue_manual(data = wilcox_3.1, y.position = c(6.5), hide.ns = TRUE, label = "p = {scales::pvalue(p.adj)}", size = 3.5)

##### s1i ####
#comparison of cell quantities between WHO grades
#cell quantities were averaged per resection
data_X4 = cell_counts |>   
  select(Study_ID, WHO_2021, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells) |> 
  mutate(Total_cells = as.numeric(Total_cells),
         CD8_T_cells = as.numeric(CD8_T_cells),
         CD4_T_cells = as.numeric(CD4_T_cells),
         T_cells      = CD8_T_cells + CD4_T_cells,
         B_cells = as.numeric(B_cells),
         Macrophages = as.numeric(Macrophages),
         Tumor_cells = as.numeric(Tumor_cells)) |>
  drop_na(Total_cells) |>
  group_by(Study_ID, WHO_2021) |> 
  filter(n()>3) |> 
  filter(!Study_ID %in% notum) |> 
  group_by(Study_ID, WHO_2021) |> 
  summarize(across(everything(), median)) |> 
  mutate_if(is.numeric, round, 0) |> 
  ungroup() |> 
  pivot_longer(cols = !c(Study_ID, WHO_2021)) |>
  mutate(value = log1p(value)) |> 
  mutate(name = case_when(name == "CD8_T_cells"       ~ "CD8_T_cells", 
                          name == "CD4_T_cells"       ~ "CD4_T_cells", 
                          name == "Macrophages"       ~ "TAMs", 
                          name == "T_cells"           ~ "T_cells",
                          name == "Total_cells"       ~ "Total", 
                          name == "B_cells"           ~ "B_cells",
                          name == "Tumor_cells"       ~ "Tumor_cells")) |> 
  filter(WHO_2021 != "NA") |> 
  filter(name != "B_cells")

wilcox_X4 = data_X4 |> 
            dplyr::group_by(name) |> 
            mutate(value = as.numeric(value)) |> 
            rstatix::wilcox_test(value ~ WHO_2021) |> 
            rstatix::adjust_pvalue(method = "fdr")

c = ggplot(data_X4, aes(x = WHO_2021, y = value)) +
    geom_boxplot(aes(fill = WHO_2021), size = .2) +
    facet_grid(~factor(name, level = .reorder_celltypes4), as.table = T, scales="free") +
    theme_bw() +
    theme(legend.position="none") +
    coord_cartesian(ylim = c(0,8)) +
    labs(y = "Log(counts + 1)", x = "WHO_2021") +
    scale_fill_manual(values = c("A2" = "#ffffe5","A3" = "#ffffbf","A4" = "#fdae61")) +
    ggpubr::stat_pvalue_manual(data = wilcox_X4, y.position = c(5.5), hide.ns = T, label = "p = {scales::pvalue(p.adj)}", size = 3.5)

pdf("s1ghi.pdf", height = 7, width = 8)
a/b/c
dev.off()

#### S2 ####
##### s2b ####
#density plot of cell distance to CD31+ structures per sample
pdf("s2b.pdf", width = 12, height = 3)
try2 |> 
  filter(result<600) |> 
  mutate(result = result*0.3441) |>
  ggplot(aes(result, col = class_bio)) +
  geom_density(linewidth = 1) +
  theme_classic() +
  scale_color_manual(values = c("#e41a1c", "#377eb8")) +
  facet_wrap(~sample, nrow = 2, scales = "free")
dev.off()

##### s2c ####
cells2        = cells[-c(5,13)]
whichcuff     = mclapply(c(1:13), function(x) st_within(cells2[[x]]$centroid, nonconvexsub[[x]], sparse = T) |>  
                           lapply(function(x) ifelse(length(x)>0,x,0)) |> 
                           unlist(), mc.cores = 15, mc.cleanup = T)
cells3        = lapply(seq_along(cells2), function(x) cbind(cells2[[x]], whichcuff[[x]])) |> 
                lapply(function(x) x |> dplyr::select(id, class_bio, centroid, `whichcuff[[x]]`))
cellscomb           = data.table::rbindlist(cells3, use.names = T, idcol = "samp")  
colnames(cellscomb) = c("samp","id", "class_bio", "centroid", "whichcuff")

barplotd = cellscomb |> 
           filter(class_bio != "Cell") |>
           group_by(samp,class_bio, whichcuff) |> 
           summarize(sum = as.numeric(n())) |>
           filter(sum >2) |> 
           mutate(rank = rank(sum))

#barplots individual samples T cell cuffs
barplotd1 = barplotd  |> filter(whichcuff != 0) |> mutate(in_cuff = ifelse(whichcuff == 0, "No", "Yes"))

plotlist = list()
for(i in 1:length(cells3)){
  plotlist[[i]] =   ggplot(data = barplotd1[barplotd1$samp == i,]) +
    geom_col(aes(sum, fct_reorder(as.character(whichcuff), rank), fill = sum)) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "top") +
    scale_fill_viridis_c(option = "magma", trans = "log", limits = c(1,400),begin = 0, end = .85) +
    coord_cartesian(xlim = c(1,320))
  
}

#barplots individual samples T cell stroma
barplotd2 = barplotd |> filter(whichcuff == 0)

plotlist2 = list()
for(i in 1:length(cells3)){
  plotlist2[[i]] =   ggplot(data = barplotd2[barplotd2$samp == i,]) +
    geom_col(aes(sum, fct_reorder(as.character(whichcuff), rank), fill = sum)) +
    scale_x_log10() +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "top") +
    scale_fill_viridis_c(option = "magma", trans = "log", limits = c(1,90000),begin = 0, end = .85)+
    coord_cartesian(xlim = c(1,11300))
  
}

order = barplotd |> 
  group_by(samp) |> 
  summarise(n = n()) |> 
  arrange(desc(n)) |> pull(samp)

# compbined barplots
a = (plotlist[[order[1]]] | plotlist[[order[2]]]  | plotlist[[order[3]]]  | plotlist[[order[4]]]  | plotlist[[order[5]]] | plotlist[[order[6]]] | plotlist[[order[7]]] | plotlist[[order[8]]] | 
       plotlist[[order[9]]] | plotlist[[order[10]]] | plotlist[[order[11]]] | plotlist[[order[12]]] | plotlist[[order[13]]]) + plot_layout(guides = "collect")

b = (plotlist2[[order[1]]] | plotlist2[[order[2]]]  | plotlist2[[order[3]]]  | plotlist2[[order[4]]]  | plotlist2[[order[5]]] | plotlist2[[order[6]]] | plotlist2[[order[7]]] | plotlist2[[order[8]]] | 
       plotlist2[[order[9]]] | plotlist2[[order[10]]] | plotlist2[[order[11]]] | plotlist2[[order[12]]] | plotlist2[[order[13]]]) + plot_layout(guides = "collect")

pdf("s2c.pdf", height = 3.5, width = 25)
a/b + plot_layout(heights = c(40,1))
dev.off()

##### s2d ####
#Test difference in representation of spatial phenotypes between initial and recurrent tumors 
samenstelling_TME = cell_counts |> 
                    dplyr::select(c(ROI_nr,Study_ID, res, cuff_score, CD8_T_cells, CD4_T_cells, B_cells, Macrophages, Tumor_cells, Total_cells)) |> 
                    dplyr::filter(cuff_score != "Unknown") |> 
                    mutate(T_cells    = CD8_T_cells + CD4_T_cells, 
                           Annotation = ifelse(cuff_score %in% c("+", "++", "+++"), "T_cells_in_PS", "T_cell_absent"),
                           Annotation = ifelse(Annotation == "T_cell_absent" & T_cells >= 3, "T_cells_in_S", Annotation), 
                           Total_cells = as.numeric(Total_cells))|> 
                    select(-cuff_score, -T_cells) |> 
                    pivot_longer(cols = !c("Study_ID", "res", "Annotation", "ROI_nr"))

tmp44 = samenstelling_TME |> 
        group_by(Study_ID, Annotation) |> 
        summarise(sum = n())

tmp45 = data.frame(sum = rowSums(table(tmp44[,c(1,2)])))  |> 
        rownames_to_column("Study_ID") |> 
        left_join(pivot_counts_classindependent[,c("Study_ID", "Sample_Type")], by = "Study_ID") |> 
        distinct(Study_ID,.keep_all = T) |> 
        drop_na(Sample_Type) 

tmp45 |>  
  group_by(sum, Sample_Type) |> 
  summarize(n()) |> 
  pivot_wider(names_from = "Sample_Type", values_from = `n()`) |> 
  column_to_rownames("sum") |> 
  fisher_test()

pdf("s2d.pdf", height = 3.5, width = 5)
ggplot(tmp45, aes(Sample_Type, sum, fill = as.factor(sum))) +
  geom_col(position = "fill") +
  theme_minimal() +
  scale_fill_manual(values = c("#99d8c9","#238b45", "#00441b"))
dev.off()

##### s2e ####
#Test difference in immune cell quantities between Spatial phenotypes corrected for tumor cell quantity
pval_corr_tumorcells = samenstelling_TME |> 
                       pivot_wider(names_from = name,
                                   values_from = value) |> 
                       pivot_longer(cols = c(5,6,7,8,10)) |> 
                       mutate(pers = value/(Tumor_cells+value)) |>
                       filter(name != "Total_cells") |> 
                       group_by(name) |> 
                       rstatix::wilcox_test(formula     = pers ~ Annotation,  
                                            paired = F) |> 
                       rstatix::adjust_pvalue(method = "fdr")

pdf("s2e.pdf", height = 5, width = 5)
samenstelling_TME |> 
  pivot_wider(names_from = name,
              values_from = value) |> 
  pivot_longer(cols = c(5,6,7,8,10)) |> 
  mutate(pers = value/(Tumor_cells+value)) |> 
  filter(name != "Total_cells") |> 
  ggplot(aes(factor(Annotation, level = c("T_cells_in_PS", "T_cells_in_S", "T_cell_absent")), pers)) +
  geom_boxplot(aes(fill = Annotation)) +
  facet_wrap(~factor(name, level = c("Macrophages", "CD4_T_cells", "CD8_T_cells", "B_cells")), scales = "fixed", nrow = 3)+
  labs(x = "Region annotation", y = "log(value + 1)") + 
  theme_bw() +
  theme(axis.text.x     = element_text(angle = 60, hjust=1), 
        legend.position = "none") +
  scale_fill_manual(values = c("#99d8c9", "#00441b", "#238b45")) +
  stat_pvalue_manual(pval_corr_tumorcells, 
                     label = "p.adj.signif", 
                     y.position = c(1.1,1.2,1.3,1.1,1.2,1.3,1.1,1.2,1.3,1.1,1.2), 
                     tip.length = 0,  
                     hide.ns = T, 
                     size = 3)+
  coord_cartesian(ylim = c(0,1.4)) +
  labs(x = "Log(count + 1)")
dev.off()

#### S3 ####
##### s3a ####
#test difference T cell cuff presence
#comparisons between Initial and Recurrent samples, and between WHO grades
tmp44 |> 
  left_join(metadata[,c("Study_ID", "Initial_Recurrent", "WHO_2021")] |> 
              distinct(Study_ID, .keep_all = T), by = "Study_ID") |>
  pivot_wider(names_from = Annotation, values_from = sum, values_fill = 0) |> 
  mutate(T_cell_cuffs_present = T_cells_in_PS > 0) |>
  drop_na(Initial_Recurrent) |>
  group_by(T_cell_cuffs_present, Initial_Recurrent) |> 
  summarize(Propotrions_of_samples = n()) |>
  ungroup() |> 
  pivot_wider(names_from = Initial_Recurrent, values_from = Propotrions_of_samples) |> 
  column_to_rownames("T_cell_cuffs_present") |> 
  fisher_test()

tmp44 |> 
  left_join(metadata[,c("Study_ID", "Initial_Recurrent", "WHO_2021")] |> 
              distinct(Study_ID, .keep_all = T), by = "Study_ID") |>
  pivot_wider(names_from = Annotation, values_from = sum, values_fill = 0) |> 
  mutate(T_cell_cuffs_present = T_cells_in_PS > 0) |>
  filter(Initial_Recurrent != "NA") |>
  filter(WHO_2021 != "NA") |>
  group_by(T_cell_cuffs_present, WHO_2021) |> 
  summarize(Propotrions_of_samples = n())|>
  ungroup() |> 
  pivot_wider(names_from = WHO_2021, values_from = Propotrions_of_samples) |> 
  column_to_rownames("T_cell_cuffs_present") |> 
  fisher_test()

a = tmp44 |> 
    left_join(metadata[,c("Study_ID", "Initial_Recurrent", "WHO_2021")] |> 
                distinct(Study_ID, .keep_all = T), by = "Study_ID") |>
    pivot_wider(names_from = Annotation, values_from = sum, values_fill = 0) |> 
    mutate(T_cell_cuffs_present = T_cells_in_PS > 0) |>
    filter(Initial_Recurrent != "NA") |>
    filter(WHO_2021 != "NA") |>
    group_by(T_cell_cuffs_present, Initial_Recurrent) |> 
    summarize(Propotrions_of_samples = n()) |> 
    ggplot(aes(Initial_Recurrent, Propotrions_of_samples, fill = T_cell_cuffs_present)) +
    geom_col(position = "fill") +
    theme_bw() +
    scale_fill_manual(values = c("#99d8c9","#238b45"))

b = tmp44 |> 
    left_join(metadata[,c("Study_ID", "Initial_Recurrent", "WHO_2021")] |> 
               distinct(Study_ID, .keep_all = T), by = "Study_ID") |>
    pivot_wider(names_from = Annotation, values_from = sum, values_fill = 0) |> 
    mutate(T_cell_cuffs_present = T_cells_in_PS > 0) |>
    filter(Initial_Recurrent != "NA") |>
    filter(WHO_2021 != "NA") |>
    group_by(T_cell_cuffs_present, WHO_2021) |> 
    summarize(Propotrions_of_samples = n()) |> 
    ggplot(aes(WHO_2021, Propotrions_of_samples, fill = T_cell_cuffs_present)) +
    geom_col(position = "fill") + 
    theme_bw() +
    scale_fill_manual(values = c("#99d8c9","#238b45"))

pdf(file = "s3a.pdf", height = 3, width = 6)
a + b + patchwork::plot_layout(guides = "collect")
dev.off()

##### s3c ####
#Visualize T cell cuff presence per resection
percentage_cuffs2 = cell_counts |> 
                    dplyr::select(c(Study_ID, res, cuff_score, GTC_score)) |> 
                    dplyr::group_by(Study_ID, res) |> 
                    dplyr::count(Study_ID, cuff_score) |> 
                    tidyr::pivot_wider(names_from = cuff_score, values_from = n, values_fill = 0) |> 
                    dplyr::mutate("Total"     = (`-`+`+`+`++`+`+++`+`NA`),
                                  "Cuff: -"   = (`-`  *100/`Total`),
                                  "Cuff: +"   = (`+`  *100/`Total`),
                                  "Cuff: ++"  = (`++` *100/`Total`),
                                  "Cuff: +++" = (`+++`*100/`Total`),
                                  "Cuff: NA"  = (`NA` *100/`Total`)) |> 
                    dplyr::mutate(ptid = substr(Study_ID, 1, nchar(Study_ID) -3)) |> 
                    dplyr::select(!c("-", "+", "++", "+++", "NA", "Total")) |> 
                    tidyr::pivot_longer(cols = !c("Study_ID", "res", "ptid"), names_to = "name", values_to = "Percentage") 

pdf(file = "s3c.pdf", width = 15,height = 5)
ggplot(percentage_cuffs2, aes(x = "", y = Percentage, fill = name)) +
  geom_bar(stat = "identity", width = 0.5, linewidth = 0.001,color="grey90") +
  coord_polar("y", start = 0, direction = 1) +
  facet_grid(rows = vars(as.factor(res)), 
             cols = vars(as.factor(ptid)), 
             switch = "y") +
  theme_minimal() +
  labs(x = "Resection", y = "Patient ID") +
  scale_fill_manual(values = c(.col_cuff3)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.spacing = unit(0.01, "lines"),
        strip.text.y.left = element_text(angle=0))
dev.off()

##### s3d ####
#Test the difference of T cell cuff presence in ROIs between Initial and Recurrent samples 
pdf(file = "s3d.pdf", height = 4, width = 4)
cell_counts |> 
  filter(Initial_Recurrent != "NA") |> 
  filter(cuff_score != "NA") |> 
  ggplot(aes(Initial_Recurrent, fill = cuff_score)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = .col_cuff_GTC ) +
  theme_bw()
dev.off()

dattest = cell_counts |> 
          filter(Initial_Recurrent != "NA") |> 
          filter(cuff_score != "NA")

dattestchisq = t(table(dattest$Initial_Recurrent, dattest$cuff_score))
fisher.test(dattestchisq)

##### s3e ####
#Boxplot of ROI cell quantities split for Cuff size
samenstelling_cuffs2 = cell_counts |> 
                       dplyr::select(c(Study_ID, res, cuff_score, CD8_T_cells, CD4_T_cells, 
                                       B_cells, Macrophages, Tumor_cells, Total_cells)) |> 
                       dplyr::group_by(cuff_score, res, Study_ID) |> 
                       mutate(Total_cells = as.numeric(Total_cells)) |> 
                       tidyr::pivot_longer(cols = !c("cuff_score", "res", "Study_ID"), names_to = "name", values_to = "values") |> 
                       mutate(name = case_when(name == "CD8_T_cells" ~ "CD8_T_cells",
                                               name == "CD4_T_cells" ~ "CD4_T_cells",
                                               name == "Macrophages" ~ "TAMs",
                                               name == "Total_cells" ~ "Total",
                                               .default = name)) |> 
                       filter(cuff_score != "NA")

wilcox_1.D = samenstelling_cuffs2 |> 
             dplyr::group_by(name) |> 
             rstatix::wilcox_test(formula     = values ~ cuff_score, 
                                  comparisons = .comparisons_GTC, 
                                  paired = F) |> 
             rstatix::adjust_pvalue(method = "fdr")

pdf("s3e.pdf", height = 5, width = 7)
ggplot(samenstelling_cuffs2, aes(x = cuff_score, y = log1p(values))) +
  geom_boxplot(aes(fill = cuff_score)) +
  facet_wrap(~factor(name, level = .reorder_celltypes4), ncol = 3) +
  scale_fill_manual(values = .col_cuff_GTC) +
  labs(x = "Cuff size", y = "log(value + 1)") + 
  theme_bw() +
  theme(legend.position = "none") + # hier plot je 1 waarde (median) per pt-res-cuff2 groep
  coord_cartesian(ylim = c(0,10.5)) +
  ggpubr::stat_pvalue_manual(data = wilcox_1.D, y.position = c(5, 6, 7, 8, 9, 10, 
                                                               5, 6, 7, 8, 9, 10,
                                                               9, 10, 
                                                               5, 6, 7, 8, 9, 10,
                                                               10), tip.length = 0,  hide.ns = T, label = "{p.adj.signif}", size = 3)
dev.off()

##### s3f ####
#Boxplot of inter-cell distances split for ROI cuff size
data_distance_test = data |>
                     left_join(metadata[,c("ROI_nr", "cuff_score")], by = "ROI_nr") |> 
                     filter(cuff_score %in% c("-","+", "++", "+++")) |>
                     select(distance_score_from_CD8_T_cells_to_CD8_T_cells,
                            distance_score_from_CD8_T_cells_to_CD4_T_cells,
                            distance_score_from_CD4_T_cells_to_CD4_T_cells, 
                            distance_score_from_B_cells_to_B_cells,
                            distance_score_from_Macrophage_cells_to_Macrophage_cells,
                            distance_score_from_Macrophage_cells_to_CD4_T_cells,
                            distance_score_from_Macrophage_cells_to_B_cells,
                            distance_score_from_Macrophage_cells_to_CD8_T_cells,
                            distance_score_from_B_cells_to_CD4_T_cells,
                            distance_score_from_CD4_T_cells_to_Macrophage_cells,
                            distance_score_from_CD8_T_cells_to_Macrophage_cells,
                            distance_score_from_CD4_T_cells_to_CD8_T_cells,
                            distance_score_from_B_cells_to_Macrophage_cells,
                            distance_score_from_CD4_T_cells_to_B_cells,
                            distance_score_from_B_cells_to_CD8_T_cells,
                            distance_score_from_CD8_T_cells_to_B_cells,
                            ROI_nr, 
                            cuff_score) |>
                     pivot_longer(cols = c(distance_score_from_CD8_T_cells_to_CD8_T_cells, 
                                           distance_score_from_CD4_T_cells_to_CD4_T_cells, 
                                           distance_score_from_B_cells_to_B_cells, 
                                           distance_score_from_Macrophage_cells_to_Macrophage_cells,
                                           distance_score_from_B_cells_to_CD4_T_cells,
                                           distance_score_from_CD8_T_cells_to_Macrophage_cells,
                                           distance_score_from_CD4_T_cells_to_Macrophage_cells,
                                           distance_score_from_CD4_T_cells_to_CD8_T_cells,
                                           distance_score_from_B_cells_to_Macrophage_cells,
                                           distance_score_from_CD4_T_cells_to_B_cells,
                                           distance_score_from_CD8_T_cells_to_CD4_T_cells,
                                           distance_score_from_B_cells_to_CD8_T_cells,
                                           distance_score_from_CD8_T_cells_to_B_cells,        
                                           distance_score_from_Macrophage_cells_to_CD4_T_cells,
                                           distance_score_from_Macrophage_cells_to_B_cells,
                                           distance_score_from_Macrophage_cells_to_CD8_T_cells)) |> 
                     filter(value != "NA") |>
                     mutate(value = as.numeric(value)) |> 
                     mutate(name = case_when(name == "distance_score_from_CD8_T_cells_to_CD8_T_cells"             ~ "CD8_T_cells-CD8_T_cells",
                                             name == "distance_score_from_B_cells_to_B_cells"                     ~ "B_cell-B_cell",
                                             name == "distance_score_from_CD4_T_cells_to_CD4_T_cells"             ~ "CD4_T_cells-CD4_T_cells",
                                             name == "distance_score_from_Macrophage_cells_to_Macrophage_cells"   ~ "TAM-TAM",
                                             name == "distance_score_from_Macrophage_cells_to_CD4_T_cells"        ~ "TAM-CD4_T_cells",
                                             name == "distance_score_from_Macrophage_cells_to_B_cells"            ~ "TAM-B_cell",
                                             name == "distance_score_from_Macrophage_cells_to_CD8_T_cells"        ~ "TAM-CD8_T_cells",
                                             name == "distance_score_from_B_cells_to_CD4_T_cells"                 ~ "B_cell-CD4_T_cells",
                                             name == "distance_score_from_CD4_T_cells_to_Macrophage_cells"        ~ "CD4_T_cells-TAM",
                                             name == "distance_score_from_CD8_T_cells_to_Macrophage_cells"        ~ "CD8_T_cells-TAM",
                                             name == "distance_score_from_CD4_T_cells_to_CD8_T_cells"             ~ "CD4_T_cells-CD8_T_cells",
                                             name == "distance_score_from_B_cells_to_Macrophage_cells"            ~ "B_cell-TAM",
                                             name == "distance_score_from_CD4_T_cells_to_B_cells"                 ~ "CD4_T_cells-B_cell",
                                             name == "distance_score_from_CD8_T_cells_to_CD4_T_cells"             ~ "CD8_T_cells-CD4_T_cells",
                                             name == "distance_score_from_B_cells_to_CD8_T_cells"                 ~ "B_cell-CD8_T_cells",
                                             name == "distance_score_from_CD8_T_cells_to_B_cells"                 ~ "CD8_T_cells-B_cell"),
                            name_to     = sapply(strsplit(name, "-"),"[",2),
                            name_from   = sapply(strsplit(name, "-"),"[",1), 
                            order  = case_when(name_to    == "B_cell"      ~ 10,
                                               name_to    == "CD4_T_cells" ~ 20,
                                               name_to    == "CD8_T_cells" ~ 30,
                                               name_to    == "TAM"         ~ 40), 
                            order2 = case_when(name_from  == "B_cell"      ~ 1,
                                               name_from  == "CD4_T_cells" ~ 2,
                                               name_from  == "CD8_T_cells" ~ 3,
                                               name_from  == "TAM"         ~ 4), 
                            order3 = order + order2)



wilcox_sup.2.f = data_distance_test |> 
                 arrange(order3) |> 
                 dplyr::group_by(name) |> 
                 rstatix::wilcox_test(formula     = value ~ cuff_score,
                                      paired = F) |> 
                 rstatix::adjust_pvalue(method = "fdr") |> 
                 rstatix::add_significance(p.col = "p.adj") |> 
                 add_xy_position(x = "name", group = "cuff_score", dodge = 0.75, step.increase = .05)


a = ggplot(data_distance_test|> filter(name %in% wilcox_sup.2.f$name[1:24]),aes(reorder(name, order3), value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.shape = NA,aes(fill = cuff_score)) +
    theme_bw() +
    theme(strip.text      = element_text(size=.1),
          axis.text.x     = element_text(angle = 60, hjust=1), 
          strip.text.x    = element_text(size = 6),
          axis.title.x    = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = .col_cuff_GTC) +
    labs(y = "z-score inter-cell distance", title = "B_cell") +
    coord_cartesian(ylim = c(-5.5,30)) +
    stat_pvalue_manual(wilcox_sup.2.f[1:24,], hide.ns = T, tip.length = 0,
                       y.position = c(22,24,28,30,
                                      22,24,30,
                                      20,22,24,26,28,30)) 

b = ggplot(data_distance_test|> filter(name %in% wilcox_sup.2.f$name[25:48]),aes(reorder(name, order3), value)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA,aes(fill = cuff_score)) +
  theme_bw() +
  theme(strip.text      = element_text(size=.1),
        axis.text.x     = element_text(angle = 60, hjust=1), 
        strip.text.x    = element_text(size = 6),
        axis.title.x    = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = .col_cuff_GTC) +
  labs(y = "z-score inter-cell distance", title = "CD8_T_cell") +
  coord_cartesian(ylim = c(-5.5,30)) +
  stat_pvalue_manual(wilcox_sup.2.f[25:48,]|> mutate(xmin = xmin - 4, xmax = xmax - 4), hide.ns = T, tip.length = 0,
                     y.position = c(22,24,26,28,
                                    22,24,
                                    
                                    20))


c = ggplot(data_distance_test|> filter(name %in% wilcox_sup.2.f$name[73:96]),aes(reorder(name, order3), value)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA,aes(fill = cuff_score)) +
  theme_bw() +
  theme(strip.text      = element_text(size=.1),
        axis.text.x     = element_text(angle = 60, hjust=1), 
        strip.text.x    = element_text(size = 6),
        axis.title.x    = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = .col_cuff_GTC) +
  labs(y = "z-score inter-cell distance", title = "TAM") +
  coord_cartesian(ylim = c(-5.5,30)) +
  stat_pvalue_manual(wilcox_sup.2.f[49:72,]|> mutate(xmin = xmin - 12, xmax = xmax - 12), hide.ns = T, tip.length = 0,
                     y.position = c(      24,30,
                                          24,28,30,
                                          22,24,29,
                                          22))

d = ggplot(data_distance_test|> filter(name %in% wilcox_sup.2.f$name[49:72]),aes(reorder(name, order3), value)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.shape = NA,aes(fill = cuff_score)) +
  theme_bw() +
  theme(strip.text      = element_text(size=.1),
        axis.text.x     = element_text(angle = 60, hjust=1), 
        strip.text.x    = element_text(size = 6),
        axis.title.x    = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = .col_cuff_GTC) +
  labs(y = "z-score inter-cell distance", title = "CD4_T_cell") +
  coord_cartesian(ylim = c(-5.5,30)) +
  stat_pvalue_manual(wilcox_sup.2.f[73:96,]|> mutate(xmin = xmin - 8, xmax = xmax - 8), hide.ns = T, tip.length = 0,
                     y.position = c(24,28,
                                    20,22,28,30,30,30,30,30
                     )) 

pdf("s3f.pdf", height = 11, width = 15)
(a | b) / (c | d) + plot_layout(guides = "collect")
dev.off()

##### s3g ####
#boxplots of CD8 T cell fraction depending on Cuff size 
testfigures3e = cell_counts |> 
                mutate(T_cells = CD8_T_cells + CD4_T_cells) |> 
                filter(T_cells > 0) |> 
                dplyr::select(c(Study_ID, res, cuff_score, CD8_T_cells, CD4_T_cells)) |> 
                dplyr::filter(cuff_score != "NA")  |> 
                mutate(sum = CD8_T_cells+CD4_T_cells, 
                       CD8frac = CD8_T_cells/sum*100,
                       CD4frac = CD4_T_cells/sum*100) |> 
                wilcox_test(CD8frac~cuff_score) |> 
                adjust_pvalue() |> 
                add_y_position()

pdf("s3g.pdf", height = 5,width = 3.5)
cell_counts |> 
  mutate(T_cells = CD8_T_cells + CD4_T_cells) |> 
  filter(T_cells > 0) |> 
  dplyr::select(c(Study_ID, res, cuff_score, CD8_T_cells, CD4_T_cells)) |> 
  dplyr::filter(cuff_score != "NA")  |> 
  mutate(sum     = CD8_T_cells+CD4_T_cells, 
         CD8frac = CD8_T_cells/sum*100,
         CD4frac = CD4_T_cells/sum*100) |>
  ggplot(aes(cuff_score, CD8frac)) +
  geom_boxplot(aes(fill = cuff_score)) +
  scale_fill_manual(values = .col_cuff_GTC) +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(y = "%CD8+ of CD3+", x = "Cuff score") +
  stat_pvalue_manual(testfigures3e, hide.ns = T, tip.length = 0, y.position = c(102, 106))
dev.off()

#### S5 ####
##### s5b ####
#Boxplots of cell quantities split for GTC score and for spatial tissue phenotypes
wilcox_test = data_C6 |> 
              left_join(samenstelling_TME[,c("ROI_nr", "Annotation")] |> 
                          distinct(ROI_nr, .keep_all = T), by = "ROI_nr") |> 
              group_by(Annotation, name) |> 
              wilcox_test(value~GTC_score)|> 
              rstatix::adjust_pvalue(method = "fdr")

pdf("s5b.pdf", height = 7,width = 10)
data_C6 |> 
  left_join(samenstelling_TME[,c("ROI_nr", "Annotation")] |> 
              distinct(ROI_nr, .keep_all = T), by = "ROI_nr") |> 
  ggplot(aes(GTC_score, value)) +
  geom_boxplot(aes(fill = GTC_score)) +
  facet_grid(factor(Annotation, level = c("T_cell_absent", "T_cells_in_S", "T_cells_in_PS"))~factor(name, level = .reorder_celltypes2)) +
  scale_fill_manual(values = c(.col_cuff_GTC), na.value = "white") +
  theme_bw() +
  theme(legend.position="top") +
  labs(x = "Gemistocyte score", y = "log(count + 1)") +
  ggpubr::stat_pvalue_manual(data = wilcox_test, tip.length = 0, y.position = c(9,8.5,8,
                                                                                9.5,
                                                                                9.5,
                                                                                9.5,9,8.5,
                                                                                9.5,9,8.5,
                                                                                9,8.5,
                                                                                8.5,8,
                                                                                8,
                                                                                9,8.5,7.5,8,
                                                                                8,
                                                                                9,8.5,7.5,8),  
                             hide.ns = T, label = "{p.adj.signif}", size = 3) +
  coord_cartesian(ylim = c(0,10))
dev.off()

##### s5c ####
#Visualization of GTC score distributions per sample 
percentage_gem = cell_counts |> 
                 dplyr::select(c(Study_ID, res, GTC_score)) |> 
                 dplyr::group_by(Study_ID, res) |> 
                 dplyr::count(Study_ID, GTC_score) |> 
                 tidyr::pivot_wider(names_from = GTC_score, values_from = n, values_fill = 0) |>
                 dplyr::mutate("Total" = (`-`+`+`+`++`+`+++`+`NA`),
                               "-"     = (`-`*100/`Total`),
                               "+"     = (`+`*100/`Total`),
                               "++"    = (`++`*100/`Total`),
                               "+++"   = (`+++`*100/`Total`)) |>
                 dplyr::mutate(ptid = substr(Study_ID, 1, nchar(Study_ID) -3)) |> 
                 dplyr::select(!c("Total","NA")) |> 
                 tidyr::pivot_longer(cols = !c("Study_ID", "res", "ptid"), names_to = "name", values_to = "Percentage") 

pdf("s5c.pdf", height = 2, width = 12)
ggplot(percentage_gem, aes(x = "", y = Percentage, fill = name)) +
  geom_bar(stat = "identity", width = 0.5, linewidth = 0.001, color="grey90") +
  coord_polar("y", start = 0, direction = 1) +
  facet_grid(rows = vars(as.factor(res)), 
             cols = vars(as.factor(ptid)), 
             switch = "y") +
  theme_minimal() +
  labs(x = "Resection", y = "Patient ID") +
  scale_fill_manual(values = c(.col_cuff_GTC), na.value = "white") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.spacing = unit(0.01, "lines"),
        strip.text.y.left = element_text(angle=0)) 
dev.off()

##### s5d ####
#correlations of averaged cell quantities and GTC score in tumor resections
cell_counts |> 
  filter(Initial_Recurrent != "NA") |> 
  filter(GTC_score != "NA") |>
  group_by(Initial_Recurrent, GTC_score) |> 
  summarise(num = n()) |> 
  pivot_wider(names_from = "Initial_Recurrent", values_from = num) |> 
  column_to_rownames("GTC_score") |> 
  chisq_test()

pdf("s5d.pdf", height = 3.5, width = 5)
cell_counts |> 
  filter(Initial_Recurrent != "NA") |> 
  filter(GTC_score != "NA") |>
  ggplot(aes(Initial_Recurrent, fill = GTC_score)) +
  geom_bar(position = "fill") +
  theme_bw() +
  scale_fill_manual(values = c(.col_cuff_GTC))
dev.off()

##### s5e ####
cell_counts |> 
  filter(WHO_2021 != "NA") |> 
  filter(GTC_score != "NA") |>
  group_by(GTC_score, WHO_2021) |> 
  summarise(num = n()) |> 
  pivot_wider(names_from = "WHO_2021", values_from = num) |> 
  column_to_rownames("GTC_score") |> 
  chisq_test()

pdf("s5e.pdf", height = 3.5, width = 5)
cell_counts |> 
  filter(WHO_2021 != "NA") |> 
  filter(GTC_score != "NA") |>
  ggplot(aes(WHO_2021, fill = GTC_score)) +
  geom_bar(position = "fill") +
  theme_bw() +
  scale_fill_manual(values = c(.col_cuff_GTC))
dev.off()

##### s5f ####
pdf("s5f.pdf", height = 5, width = 9)
data_C6 |> 
  mutate(gemperstamp_value = case_when(GTC_score == "-" ~ 0,
                                       GTC_score == "+" ~ 1,
                                       GTC_score == "++" ~ 2,
                                       GTC_score == "+++" ~ 3)) |> 
  group_by(Study_ID, name) |> 
  summarize(nROI    = n(),
            Gemsum  = mean(gemperstamp_value),
            cellsum = mean(value)) |>
  ggplot(aes(Gemsum, cellsum)) +
  geom_point(size = 1.5, alpha = .5) +
  theme_bw() +
  stat_cor(cor.coef.name = "r", method = "pearson") +
  geom_smooth(method = "lm") +
  facet_wrap(~factor(name, levels = c("Total_cells", "Tumor_cells", "Macrophages", "CD4_T_cells", "CD8_T_cells", "B_cells")), scales = "free")
dev.off()

##### s5g ####
#DE test of NanoString spatial RNA assay for GTC-high vs GTC-low ROIs
gemann = annotation_NS_RNA |> 
         column_to_rownames("name") |> 
         filter(GTC_score != "NA") 

gemann = gemann[colnames(dat.norm.quant),]

dds    = DESeqDataSetFromMatrix(countData = round(dat.norm.quant),
                                colData   = gemann, 
                                design    = ~ annotation)

sizeFactors(dds) = 1
dds              = DESeq(dds)
res              = results(dds, name="annotation_low_vs_high") |> 
  as.data.frame()

p = res |> 
    rownames_to_column("gene") |> 
    left_join(markers |>  filter(Celltype %in% c("Tcell", "mac", "mic")), by = "gene") |> 
    mutate(Celltype = ifelse(padj < 0.05 &  log2FoldChange < -.5 | padj < 0.05 & log2FoldChange > .5, Celltype, NA),
           Celltypeorder = ifelse(is.na(Celltype), 1,2),
           Celltype = ifelse(is.na(Celltype), NA, Celltype), 
           Celltype = case_when(Celltype == "mic" ~ "TAM",
                                Celltype == "mac" ~ "TAM",
                                .default = Celltype), 
           Name  = ifelse(padj < 0.05 &  log2FoldChange < -.5 | padj < 0.05 & log2FoldChange > .5, gene, NA),
           size  = ifelse(is.na(Celltype), 1,1.1)) |>
    arrange(Celltypeorder) |> 
    distinct(gene, .keep_all = T)

pdf("s5g.pdf",height = 3, width = 5)
ggplot(p, aes(log2FoldChange*-1, -log10(padj),  col = Celltype, size = size)) +
  geom_point() +
  theme_bw()    +
  scale_size(range = c(1,2)) +
  scale_color_manual(values = c("#f46d43"), na.value = "grey")+
  coord_cartesian(xlim = c(-2,2)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  geom_vline(xintercept = c(-.5,.5), linetype = 2)+
  xlab("Log2FoldChange")
dev.off()

##### s5h ####
#correlation of z scores for GTC high/low DE test of bulk and spatial RNA readouts
spatial_bulk_comparison = bulkdetest2 |>
                          rownames_to_column("X") |> 
                          dplyr::rename("name" = X, 
                                        "statbulk" = stat) |> 
                          dplyr::select(name, statbulk) |> 
                          inner_join(res |> dplyr::rename("statspatial" = stat) |>
                                       rownames_to_column("name") |> 
                                       dplyr::select(name, statspatial), by = "name")

pdf("s5h.pdf", width = 4, height = 3.5)
ggplot(spatial_bulk_comparison, aes(statbulk, statspatial*-1)) +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_cor(method = "pearson", cor.coef.name = "r") +
  xlab("Bulk_stat") +
  ylab("Spatial_stat")
dev.off()

##### s5i ####
#UMAP representations of snRNAseq data showing colors for the 7 different samples 
cols = c("#fdbf6f","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")

Idents(glioma.integrated) = "state"
a = DimPlot(glioma.integrated, shuffle = T) +
  scale_color_manual(values = cols)

Idents(glioma.tumor) = "state"
b = DimPlot(glioma.tumor, shuffle = T) +
  scale_color_manual(values = cols)

Idents(TAMs) = "state"
c = DimPlot(TAMs, shuffle = T) +
  scale_color_manual(values = cols)

png(file = "s5i.png", width = 12, height = 3, res = 300, units = "in")
a + b + c
dev.off()

#### S6 ####
##### s6a ####
#violin plots of single cell feature counts for 7 snRNA-seq samples
vln = VlnPlot(glioma.integrated, features = "nFeature_RNA", group.by = "state")$data

pdf("s6a.pdf", width = 3, height = 2.5)
ggplot(vln, aes(ident, nFeature_RNA)) +
  geom_violin() +
  theme_bw() +
  labs(x = "Sample", y = "Feature count")
dev.off()

##### s6b ####
#violin plot showing enrichment scores for marker genes of standard cell types in the central nervous system in snRNAseq data
glioma.integrated@meta.data$bioidentsintegrated4 = ifelse(glioma.integrated@meta.data$bioidentsintegrated %in% c("Oligo-like tumor state", "Astro-like tumor state", "stem-like tumor state", 
                                                                                                                 "G1/S/G2/M"), "Undefined", glioma.integrated@meta.data$bioidentsintegrated) 

Idents(glioma.integrated) = "bioidentsintegrated4"
cells = unique(markers$Celltype)
dat   = VlnPlot(glioma.integrated, features = cells[c(1:6)], stack = T)$data

.celltype_singlecell = c("Undefined","T cells","Endothelial cells","Neurons","TAMs","Oligodendrocytes","Astrocytes")

pdf("s6b.pdf", height = 4, width = 8)
ggplot(dat, aes(factor(ident, level = .celltype_singlecell),expression, fill = ident, col = ident)) +
  geom_hline(yintercept = 0, linetype = 2)+
  geom_violin(alpha = .2) +
  theme_classic() +
  facet_wrap(~factor(feature, level = c("ast", "oli", "mic", "neu", "end", "Tcell")),ncol = 6) +
  labs(y = "Enrichment score", x = "Celltypes") +
  theme(legend.position = "none", axis.text.x=element_text(size = 10,angle = 60, hjust=1))+
  coord_flip()
dev.off()

##### s6c ####
#violin plot showing enrichment scores for marker genes of known single tumor cell types in snRNAseq data
markerssub   = markers |> 
               dplyr::filter(Celltype %in%  c("G2/M", "G1/S")) |> 
               mutate(Celltype = "G1.S.G2.M") |>  
               drop_na(gene)

glioma.tumor = cellscore(glioma.tumor, markers = markerssub)

datvln = data.frame("ident"            = VlnPlot(glioma.tumor, group.by = "seurat_clusters", 
                                                 features = c("stem.like.tumor"))$data$ident,
                    "stem-like tumor"  = VlnPlot(glioma.tumor, group.by = "seurat_clusters", 
                                                 features = c("stem.like.tumor"))$data$stem.like.tumor,
                    "astro-like tumor" = VlnPlot(glioma.tumor, group.by = "seurat_clusters", 
                                                 features = c("astro.like.tumor"))$data$astro.like.tumor,
                    "oligo-like tumor" = VlnPlot(glioma.tumor, group.by = "seurat_clusters", 
                                                 features = c("oligo.like.tumor"))$data$oligo.like.tumor,
                    "G1/S/G2/M"        = VlnPlot(glioma.tumor, group.by = "seurat_clusters", 
                                                 features = c("G1.S.G2.M"))$data$G1.S.G2.M) |> 
  group_by(ident) |> 
  mutate(order = mean(astro.like.tumor)) |> 
  pivot_longer(cols = !c(ident, order))

pdf("s6c.pdf", height = 5, width = 6)
ggplot(datvln, aes(reorder(ident, order), value))+
  geom_hline(yintercept = 0, col = "red", size = .2)  +
  geom_violin(size = .2) +
  facet_wrap(~factor(name, levels = c("astro.like.tumor", "oligo.like.tumor", "stem.like.tumor", "G1.S.G2.M")), ncol =  1, scales = "free_y", ) +
  theme_bw() +
  labs(x = "Identity", y = "Enrichment score") +
  theme(strip.text = element_text(size = 20))
dev.off()

##### s6d ####
#CNV estimates based on snRNA-seq data 
setwd("Dir_to_CNV_files")

allsamp = list()
for (i in seq_along(list.files(getwd()))) {
  tmp            = readRDS(paste(getwd(), list.files(getwd())[i],  sep = "/"))
  tumorcellnames = names(tmp@observation_grouped_cell_indices)
  tmp2           = list()
  
  for (j in 1:length(tumorcellnames)) {
    namelist = colnames(tmp@expr.data[,tmp@observation_grouped_cell_indices[[tumorcellnames[j]]]])
    tmp2[[j]] = tmp@expr.data[,namelist] |>  
      as.data.frame()               |> 
      rownames_to_column("gene")    |> 
      pivot_longer(cols = !gene)    |> 
      left_join(tmp@gene_order                                       |> 
                  rownames_to_column('gene')                         |> 
                  mutate(mid = start + (stop - start)), by = "gene") |> 
      group_by(gene, mid, chr)      |> 
      summarize(x = median(value))
    colnames(tmp2[[j]]) = c("gene", "mid", "chr", tumorcellnames[j])
  }
  
  names(tmp2) = tumorcellnames
  stopifnot(all(tmp2$`Astro-like tumor state`$gene == tmp2$`G1/S/G2/M`$gene) && all(tmp2$`Astro-like tumor state`$gene == tmp2$`Oligo-like tumor state`$gene))
  
  df = tmp2[[1]]
  for(l in 2:length(tumorcellnames)){
    df[,tumorcellnames[l]] = tmp2[[tumorcellnames[l]]][,tumorcellnames[l]]
  }
  df$sample = list.files(getwd())[i]
  
  allsamp[[i]] = df
}

try = bind_rows(allsamp)

pdf("s6d.pdf", height = 4.5, width = 8)
try |> 
  pivot_longer(cols = c("Astro-like tumor state", "Oligo-like tumor state", "G1/S/G2/M","stem-like tumor state","undetermined tumor state")) |> 
  filter(name != "undetermined tumor state") |> 
  mutate(Chromosome = substr(chr, 4,5), 
         Chromosome = substr(paste0(0,Chromosome),nchar(paste0(0,Chromosome))-1,nchar(paste0(0,Chromosome)))) |>
  ggplot(aes(mid,value,col = name)) +
  geom_hline(yintercept = 1, size = .5) +
  geom_smooth(span = .6, size = .5)+
  facet_grid(sample~Chromosome, scales = "free_x",space = "free_x") +
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  coord_cartesian(ylim = c(0.90, 1.12)) +
  scale_color_manual(values = c("#381D2A", "#FF6F59", "#3E6990", "#06A77D"))
dev.off()

##### s6e ####
#CNV estimates based on bulk methylation data
ids = data.frame(met_ID       = c("203519500059_R06C01", 
                                  "203519500059_R04C01", 
                                  "203175700014_R05C01",
                                  "203519500055_R04C01",
                                  "203175830107_R06C01",
                                  "203430580025_R05C01",
                                  "203430580025_R04C01"),
                 sampleID     = c("1","6","2","3","4","5_1","5_2"))


pdf("s6e.pdf", height = 4, width = 8)
cnv_data_bulk  |> 
  dplyr::filter(type %in% c("loc.start", "loc.end")) |> 
  mutate(ID = case_when(ID == "203519500059_R06C01"~"1", 
                        ID == "203519500059_R04C01"~"6", 
                        ID == "203175700014_R05C01"~"2",
                        ID == "203519500055_R04C01"~"3",
                        ID == "203175830107_R06C01"~"4",
                        ID == "203430580025_R05C01"~"5_1",
                        ID == "203430580025_R04C01"~"5_2")) |> 
  ggplot(aes(x = pos, y = seg.mean.l2fc, group = id, col = after_stat(y))) +
  geom_hline(yintercept = 0, linetype = 1, size = .5) +
  facet_grid(ID~chrom, scales = "free", space = "free") +
  geom_line(lwd = 1) +
  geom_line(data =  cnv_data_bulk |>  dplyr::filter(type == "vline"),lty=2,col="black") +
  theme_minimal() +
  theme(axis.title         = element_text(face = "bold", size = rel(1)),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x        = element_text(angle = 90, vjust = 0.45, hjust = 1),
        panel.spacing      = unit(0.1, "lines"), 
        legend.position    = "none") +
  coord_cartesian(ylim = c(-.5, .5)) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  scale_colour_gradientn(colours = c("#053061", "#053061", "#2166ac", "#2166ac", "#4393c3", "#92c5de",
                                     "darkgrey","#f4a582", "#d6604d", "#b2182b", "#67001f", "#67001f"),)
dev.off()

##### s6g ####
#Correlation between bulk and snRNA-seq CNV profiles
alldat  = rbind(bulkcnv)
alldat2 = cnv_data_bulk |> 
          left_join(annotationbulk[,c("methylation.sid", "Sample_Name")], by = c("ID" = "methylation.sid")) |> 
          mutate(join = paste(Sample_Name, name, sep = "_"))

snrna_samples = data.frame(sample2 = c("1","2","3","4","5_1", "5_2", "6"),
                           Study_ID = c("110_R1", "115_R1", "129_R1", "131_R3", "147_R1", "147_R2", "113_R1"))

cnv_data_bulk$pos = cnv_data_bulk$pos*1000000
alldat3           = alldat2[,c(1,2,6,8,9,10,14)] |> 
                    pivot_wider(names_from = type, values_from = pos)
alldat3$unique    = seq_along(alldat3$ID)

tmp100 = lapply(seq_along(try2$Study_ID),function(x) alldat3[which(alldat3$Sample_Name %in% try2[x,]$Study_ID & 
                                                                   alldat3$chrom       %in% try2[x,]$chr &
                                                                   alldat3$loc.start    <   try2[x,]$mid & 
                                                                   alldat3$loc.end      >   try2[x,]$mid),]$unique)

try2$unique = sapply(tmp100, function(x) ifelse(is_empty(x), NA, x))

tmp300 = try2 |> 
         group_by(unique) |> 
         summarize(value = mean(value))

pdf("s6g.pdf", width = 3.5, height = 3)
alldat3 |> 
  left_join(tmp300[,c("unique", "value")], by = "unique") |>
  mutate(value = value - 1) |> 
  ggplot(aes(`seg.mean.l2fc`, value)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_smooth(method = "lm") +
  stat_cor()
dev.off()

##### s6h ####
#Expression of GTC-high bulk DE genes in snRNAseq cell populations
feat     = head(VariableFeatures(glioma.integrated, assay = "integrated"),2000)
marktest = data.frame(Celltype = "gemistocyte", gene = feat[feat %in% gemmark])
col_fun  = circlize::colorRamp2(c(-.5, 0, 2.5), viridis(20)[c(1,5, 20)])
tmp      = AverageExpression(object   = glioma.integrated, 
                             group.by = "bioidentsintegrated", 
                             features = marktest$gene)

pdf("s6h.pdf", height = 10, width = 3)
Heatmap(t(scale(t(tmp$SCT))), 
        col             = col_fun,
        row_names_gp    = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 8))
dev.off()

##### s6i ####
#density plot of enrichment scores for GTC DE genes in snRNAseq tumor cell populations
pdf("s6i.pdf", height = 2, width = 2.5)
ggplot(glioma.tumor@meta.data, aes(gemistocyte)) + 
  geom_density() + 
  geom_vline(xintercept = .2, linetype = 2) +
  theme_classic() +
  labs(x= "Gemistocyte score")
dev.off()

##### s6j ####
#cell numbers and proportions in the single cell GTC cluster
glioma.tumor@meta.data$gemistocyte = glioma.tumor$gemistocyte

pdf("s6j.pdf", height = 3, width = 7)
glioma.tumor@meta.data |>
  mutate(gemstate = ifelse(gemistocyte > 0.2, "gemistocyte", "non_gemistocyte")) |> 
  group_by(state, gemstate) |>
  summarise(sum = n()) |> 
  pivot_wider(names_from = gemstate, values_from = sum) |>  
  mutate(a_Percentage = (gemistocyte/(non_gemistocyte + gemistocyte))*100, 
         b_Count = gemistocyte/38) |>
  pivot_longer(col = c(b_Count, a_Percentage)) |>
  ggplot(aes(state, value, fill = name)) +
  geom_col(position  = "dodge") +
  coord_cartesian(ylim = c(0,5)) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 38, name = "Count")) +
  theme_classic() +
  ylab("Percentage") +
  scale_fill_manual(values = c("#f46d43","#74add1"))
dev.off()

#### S7 ####
##### s7b ####
#requires Figure 5D
#heatmap of NanoString spatial protein assay data
col_fun         = colorRamp2(c(0, 200), c("white", "red"))

column_ha       = HeatmapAnnotation(T_cells     = as.numeric(annotation_vis$Tcellnumber), 
                                    Region      = annotation_vis$Segment_tags, 
                                    Gemistocyte = annotation_vis$GTC_score, 
                                    col = list(T_cells     = col_fun, 
                                               Region      = c("T_cells_present_stroma" = "#d9f0d3",
                                                               "T_cells_present_vessel" = "#7fbf7b", 
                                                               "T_cells_absent"         = "#1b7837"),
                                               Gemistocyte = c("No" = "white", "Low" = "#d9f0d3", "Medium" = "#7fbf7b", 
                                                               "High" = "#1b7837", "NA" = "lightgrey")))

dat = t(scale(t(log(dat.norm.neg))))
dat = dat[,annotation_vis$name]

pdf("s7b.pdf",height = 9, width = 12)
Heatmap(dat,
        top_annotation = column_ha,
        row_km = 2, 
        column_km = 2,
        column_labels = annotation_vis$name, 
        show_parent_dend_line = FALSE)
dev.off()

##### s7c ####
#boxplot showing near and far edge distance of ROIs annotated as "T cells in perivascular space"
tmp = ROIdistance |> 
      drop_na(distance_vessel_close) |> 
      pivot_longer(cols = 4:5)

pdf("s7c.pdf", height = 3, width = 3)
ggplot(tmp, aes(name, as.numeric(value))) +
  geom_boxplot() +
  theme_classic()
dev.off()

##### s7d ####
#Requires figure 4b
#Association of all WGCNA gene modules with ROI characteristics in NanoString spatial RNA assay data
pdf("s7d.pdf", width = 5,height = 6)
par(mar = c(10, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor[,c(-2,-6)],
               xLabels = colnames(moduleTraitCor)[c(-2,-6)],
               yLabels = paste("module", c(1,2,8,9,3,4,10,5,11,6,12,7,13)),
               ySymbols = rownames(moduleTraitCor),
               colors = blueWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = .9,
               cex.lab.y = 1.2,
               cex.lab.x =1.2,
               zlim = c(-1,1))
dev.off()

##### s7e ####
#Dotplot showing enrichment scores of all WGCNA gene modules in snRNAseq cell populations
gemistocyte_cells = rownames(glioma.tumor@meta.data)[glioma.tumor$gemistocyte>0.2]
glioma.integrated$cellid = rownames(glioma.integrated@meta.data)
glioma.integrated$bioidentsintegrated2 = ifelse(glioma.integrated$cellid %in% gemistocyte_cells, "Gemistocyte tumor state", glioma.integrated$bioidentsintegrated)
Idents(glioma.integrated) = "bioidentsintegrated2"

modulessub = modulesspatial |> 
  select(substanceBXH, moduleColor) |> 
  dplyr::rename("gene" = substanceBXH,
                "Celltype" = moduleColor)

glioma.integrated = cellscore(glioma.integrated, modulessub)

tmp = DotPlot(glioma.integrated, features = unique(modulessub$Celltype))$data |> 
  filter(features.plot != "grey") |> 
  mutate(features.plot = case_when(features.plot == "black"       ~ "Module 1",
                                   features.plot == "turquoise"   ~ "Module 2",
                                   features.plot == "salmon"      ~ "Module 8",
                                   features.plot == "blue"        ~ "Module 9",
                                   features.plot == "tan"         ~ "Module 3",
                                   features.plot == "pink"        ~ "Module 4",
                                   features.plot == "purple"      ~ "Module 10",
                                   features.plot == "yellow"      ~ "Module 5",
                                   features.plot == "brown"       ~ "Module 11",
                                   features.plot == "green"       ~ "Module 6",
                                   features.plot == "greenyellow" ~ "Module 12",
                                   features.plot == "red"         ~ "Module 7",
                                   features.plot == "magenta"     ~ "Module 13"))

order_try = c("Astrocytes", "T cells",'Oligodendrocytes', "Neurons", "Endothelial cells", "TAMs",  
              "Gemistocyte tumor state","Astro-like tumor state", "Oligo-like tumor state", "stem-like tumor state", "G1/S/G2/M")

order_try2 = paste("Module", c(10, 11, 12, 8, 1, 2, 5, 3,4, 13, 7, 6, 9))

tmp$features.plot = factor(tmp$features.plot, levels = order_try2)
tmp$id            = factor(tmp$id, levels = order_try)

pdf("s7e.pdf", height = 4.5, width = 7)
ggplot(tmp, aes(features.plot, id, size = pct.exp, col= avg.exp.scaled)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust=1), 
        axis.title = element_blank()) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()

##### s7f ####
#UMAP representation of TAM snRNAseq data showing enrichment for WGCNA gene modules
pdf("s7f.pdf", height = 3, width = 12)
FeaturePlot(TAMs, features = c("Module_2", "Module_3", "Module_5"),pt.size = .3, order = T,ncol = 3) &
  scale_colour_gradientn(colours = c("grey90","grey90","grey90","#92c5de","#2166ac","#053061","#053061")) 
dev.off()

##### s7g ####
#Density plot of Module-3 enrichment scores in snRNAseq TAM population
pdf("s7g.pdf", height = 2, width = 3)
TAMs@meta.data |> 
  ggplot(aes(Module_3))+
  geom_density() + 
  geom_vline(xintercept = 0, linetype = 2) +
  theme_classic()
dev.off()

##### s7h ####
#Net analysis in cellchat 
ht1 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "outgoing", width = 4, height = 17, font.size = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellChat, pattern = "incoming", width = 4, height = 17, font.size = 10)
pdf(file = "s7h.pdf", height = 10,width = 8)
ht1+ht2
dev.off()

#### s8 ####
##### s8d ####
#Images of GeoJSON exports validation multiplex IF stainings
size   = .5
a = ggplot() +
  geom_sf(data = density_try[[1]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[1]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[1]][s[[1]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void() +
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
b = ggplot() +
  geom_sf(data = density_try[[2]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[2]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[2]][s[[2]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
c = ggplot() +
  geom_sf(data = density_try[[3]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[3]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[3]][s[[3]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
d = ggplot() +
  geom_sf(data = density_try[[4]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[4]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[4]][s[[4]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
e = ggplot() +
  geom_sf(data = density_try[[5]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[5]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[5]][s[[5]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
f = ggplot() +
  geom_sf(data = density_try[[6]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[6]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[6]][s[[6]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
g = ggplot() +
  geom_sf(data = density_try[[7]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[7]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[7]][s[[7]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
h = ggplot() +
  geom_sf(data = density_try[[8]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[8]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[8]][s[[8]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
i = ggplot() +
  geom_sf(data = density_try[[9]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[9]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[9]][s[[9]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
j = ggplot() +
  geom_sf(data = density_try[[10]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[10]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[10]][s[[10]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
k = ggplot() +
  geom_sf(data = density_try[[11]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[11]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[11]][s[[11]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
l = ggplot() +
  geom_sf(data = density_try[[12]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[12]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[12]][s[[12]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
m = ggplot() +
  geom_sf(data = density_try[[13]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[13]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[13]][s[[13]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
n = ggplot() +
  geom_sf(data = density_try[[14]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[14]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[14]][s[[14]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)
o = ggplot() +
  geom_sf(data = density_try[[15]], aes(geometry= geometry, fill = CD44_density), color = NA) +
  geom_sf(data =  cd44border[[15]][[1]], fill = NA, color = "black", linewidth = 1) +
  geom_point(data = s[[15]][s[[15]]$mark == "T_cell",], aes(x,y, col = mark), size = size) +
  theme_void()+
  scale_color_manual(values = "darkorange1") +
  scale_fill_viridis(option="mako") +
  theme(legend.position = "none") +
  geom_line(data = data.frame(x = c(1000, 6797.1), y = c(1000,1000)), aes(x=x,y=y), linewidth = 2)

pdf("s8d.pdf", height = 15, width = 30)
((a|b|c|d|e)/(f|h|i|j|l))/ 
  (g|k|m|n) + plot_layout(guides = "collect")
dev.off()

##### s8e ####
#density plot of CRYAB intensity signal in 30 stainings of whole slide tumor samples
pdf("s8e.pdf", height = 2, width = 4)
CRYAB_data |>
  dplyr::filter(FITC..Cytoplasm..Median<2000) |> 
  ggplot(aes(x = FITC..Cytoplasm..Median)) +
  geom_density(alpha = .3) +
  geom_vline(xintercept = 200, linetype = 2, col = "darkred") +
  theme_classic()
dev.off()