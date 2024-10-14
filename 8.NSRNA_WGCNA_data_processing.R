##############library###########################

library(svglite)
library(tidyverse)
library(WGCNA)
library(readxl)
library(gridExtra)
library(VennDetail)

###############options##########################

options(stringsAsFactors = FALSE)

#############load data##########################

nsdat    = readRDS("dir_normalized_NS_RNA_data") #according to van Hijfte et al. 2024; https://doi.org/10.1016/j.isci.2022.105760
nsdat    = log2(nsdat)

############process data#########################

#filter data and metadata
nsdat     = nsdat[, traitData$Resection != "Control"]

#prepare data for analysis
datExpr0  = as.data.frame(t(nsdat))

#flag bad quality genes
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#remove bad genes and samples
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#############sample filtering######################

sampleTree = hclust(dist(datExpr0), method = "average")

# Plot the sample tree
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
abline(h = 25, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 25, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust == 1)
datExpr     = datExpr0[keepSamples, ]
nGenes      = ncol(datExpr)
nSamples    = nrow(datExpr)

samples             = rownames(datExpr)

#collect garbage
collectGarbage()

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

sft = pickSoftThreshold(datExpr, 
                        powerVector = powers, 
                        verbose = 5,
                        networkType = "signed")

# Plot the results

par(mfrow = c(1, 2))

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2], labels = powers, cex=.9 ,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))

text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = .9, col = "red")

dev.off()

#determine SpftPower
softPower = 16

adjacency = adjacency(datExpr, power = softPower, type = "signed")

########calculate Topological Overlap Matrix############

TOM      = TOMsimilarity(adjacency, TOMType = "signed")
dissTOM  = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro            = geneTree, 
                            distM             = dissTOM,
                            deepSplit         = 3, 
                            pamRespectsDendro = FALSE,
                            minClusterSize    = 30)

# !!!!!!label 0 is for genes without strong associations!!!!!!!
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#############merging of similar modules###############

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs    = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
dev.off()

# merge similar clusters
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

# Plot new clusters
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#################save data #######################

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder   = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs          = mergedMEs

# Save module colors and labels for use in subsequent parts
save(datExpr, MEs, moduleLabels, moduleColors, geneTree, dissTOM, 
     file = "dir_rdata_file.rdata")
