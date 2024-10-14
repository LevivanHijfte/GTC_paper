################
#
#based on vignette: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
#
###############
####Library####

library(CellChat)
library(Seurat)
library(ggalluvial)
library(NMF)

###Data####

glioma.integrated = readRDS("integrated_glioma_data")
load("CellChat_data_standard_workflow")

####processing####

labels         <- Idents(glioma.integrated)
meta           <- data.frame(group = labels, row.names = names(labels)) 
cellChat       <- createCellChat(object = glioma.integrated, group.by = "ident", assay = "SCT")
CellChatDB     <- CellChatDB.human
showDatabaseCategory(CellChatDB.human)
CellChatDB.use <- subsetDB(CellChatDB.human, search = c("Cell-Cell Contact", "Secreted Signaling", "ECM-Receptor"))
future::plan("multisession", workers = 15)
options(future.globals.maxSize= 600000000)
cellChat@DB    <- CellChatDB.use
cellChat       <- subsetData(cellChat)
cellChat       <- identifyOverExpressedGenes(cellChat)
cellChat       <- identifyOverExpressedInteractions(cellChat)
cellChat       <- computeCommunProb(cellChat, )
cellChat       <- filterCommunication(cellChat, min.cells = 9)
cellChat       <- computeCommunProbPathway(cellChat)
cellChat       <- aggregateNet(cellChat)
saveRDS(cellChat, file = "dir_file_export_CellChat")
