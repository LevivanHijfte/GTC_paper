####Library####
library(readxl)
library(WGCNA)
library(svglite)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(readxl)

# data

annotation    = read_xlsx("dir_publication_table_S3",  skip = 2)
load("result_preprocessing_WGCNA")

# Define numbers of genes and samples
nGenes       = ncol(datExpr)
nSamples     = nrow(datExpr)

datTraits    = annotation |> 
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

###########individual genes to trait in individual modules###########

# names (colors) of the modules
modNames                    = substring(names(MEs), 3)
geneModuleMembership        = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue                    = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue)             = paste("p.MM", modNames, sep="")

namestraits           = colnames(datTraits)
geneTraitSignificance = as.data.frame(matrix(nrow = ncol(datExpr), ncol = length(namestraits)))
GSPvalue              = as.data.frame(matrix(nrow = ncol(datExpr), ncol = length(namestraits)))

for (trait in 1:length(namestraits)) {
  geneTraitSignificance[,trait]   = cor(datExpr, as.data.frame(datTraits[,trait]), use = "p")
  GSPvalue[,trait]                = corPvalueStudent(as.matrix(geneTraitSignificance[,trait]), nSamples)
}

colnames(geneTraitSignificance) = paste("GS.", namestraits, sep="")
names(GSPvalue)              = paste("p.GS.", namestraits, sep="")

##################sumarize gene info#############################

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = names(datExpr),
                       moduleColor  = moduleColors,
                       geneTraitSignificance,
                       as.data.frame(GSPvalue))

# Order modules by their significance for tcells
modOrder = order(-abs(cor(MEs, tcell_per_area, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.T_cell_count))
geneInfo  = geneInfo0[geneOrder, ]

saveRDS(geneInfo, "dir_spatialgenemodules.RDS")
