######library######

library(tidyverse)
library(readxl)
library(ggcorrplot)
library(RColorBrewer)

######functions######

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

######Data######

data           = read_xlsx("d:/Documents/Results/Gemistocytes/manuscript/Submission/Nature Communications/7.Supplementary_tables.xlsx", sheet = "Table S8",  skip = 2) |> 
                 column_to_rownames("ROI")
proteinan      = read_xlsx("d:/Documents/Results/Gemistocytes/manuscript/Submission/Nature Communications/7.Supplementary_tables.xlsx", sheet = "Table S5",  skip = 2)
colors         = brewer.pal(4,"Reds")
colors2        = brewer.pal(6,"Blues")

# Control markers
housekeepers         = c("S6", "Histone H3", "GAPDH")
Negcontrol           = c("Rb IgG", "Ms IgG1", "Ms IgG2a")

####visualization raw data####

#check for referential clustering based on AOI Type or Tissue ID

#Housekeeper geomean: this captures signal strength.
#IgG geomean: this captures background (negative controls), but in most experiments also reflects signal strength, 
#as AOIs with more on-target signal also have more background.

#Negative control geomean
igg        = data[grepl("IgG", rownames(data)),]

iggcontrol = data.frame(gm_mean = sort(apply(data[grepl("IgG", rownames(data)),], 2, gm_mean))) %>%
             rownames_to_column(var = "name") %>%
             inner_join(proteinan[,c(1,3,8,6)], by = "name")

for (i in 1:nrow(iggcontrol)) {
  if (iggcontrol$Segment_tags[i] == "T_cells_present_vessel") 
  { iggcontrol$coltag[i] = colors[1] } 
  else if (iggcontrol$Segment_tags[i] == "T_cells_present_stroma")
  { iggcontrol$coltag[i] = colors[2] } 
  else if (iggcontrol$Segment_tags[i] == "Control")
  { iggcontrol$coltag[i] = colors[3] } 
  else 
  { iggcontrol$coltag[i] = colors[4] }
}

for (i in 1:nrow(iggcontrol)) {
  if (iggcontrol$sample[i] == 1)
  { iggcontrol$colslide[i] = colors2[1] }
  else if (iggcontrol$sample[i] == 2)
  { iggcontrol$colslide[i] = colors2[2] }
  else if (iggcontrol$sample[i] == 3)
  { iggcontrol$colslide[i] = colors2[3] }
  else if (iggcontrol$sample[i] == 4)
  { iggcontrol$colslide[i] = colors2[4] }
  else if (iggcontrol$sample[i] == 5)
  { iggcontrol$colslide[i] = colors2[5] }
  else
  { iggcontrol$colslide[i] = colors2[6] }
}

par(mfrow = c(3,1), oma=c(0, 0, 0, 10))
 hist(log2(iggcontrol$gm_mean), main = "Histogram IgG control")
 barplot(sort(log2(iggcontrol$gm_mean)), col = iggcontrol$coltag, las = 2, main = "IgG control segment tag")
       legend(par('usr')[2], par('usr')[4], 
       legend = unique(iggcontrol$Segment_tags), fill = unique(iggcontrol$coltag), bty='n', xpd=NA)
 barplot(sort(log2(iggcontrol$gm_mean)),  col = iggcontrol$colslide,  las = 2, main = "IgG control slide number")
       legend(par('usr')[2], par('usr')[4], 
       legend = unique(iggcontrol$sample), fill = unique(iggcontrol$colslide), bty='n', xpd=NA)
dev.off()

######QC######

qc.roi.log = data.frame("FOV_registration"     = as.numeric(proteinan$Fov_counted)/as.numeric(proteinan$Fov_count)>=.75,
                        "Binding_density"      = as.numeric(proteinan$BindingDensity)>=0.1 & as.numeric(proteinan$BindingDensity)<=2.25, #for nCounter MAX or FLEX: 0.1-2.25, for nCounter SPRINT 0.1-1.8
                        "Positive_control"     = t(mean(as.numeric(data[grepl("POS", rownames(data)),]))/data[grepl("POS", rownames(data)),]>=0.3&
                                                   mean(as.numeric(data[grepl("POS", rownames(data)),]))/data[grepl("POS", rownames(data)),]<=3), 
                        "Minimum_surface_area" = as.numeric(proteinan$AOI_surface_area)>10000,
                        "Minimum-nuclei_count" = as.numeric(proteinan$AOI_nuclei_count)>50)

qc.roi     = data.frame("FOV_registration"     = as.numeric(proteinan$Fov_counted)/as.numeric(proteinan$Fov_count),
                        "Binding_density"      = as.numeric(proteinan$BindingDensity), 
                        "Positive_control"     = t(mean(as.numeric(data[grepl("POS", rownames(data)),]))/data[grepl("POS", rownames(data)),]),
                        "Minimum_surface_area" = as.numeric(proteinan$AOI_surface_area),
                        "Minimum_nuclei_count" = as.numeric(proteinan$AOI_nuclei_count))

#redo QC without areas with normalization factor area ~900 as this is an indication of a hybridization issue

roi.fail       = qc.roi %>% filter(HYB.POS>800) %>% rownames_to_column(var = "names") %>% pull(names)

datacounts_fil = data[, !colnames(data) %in% roi.fail]
annotation_fil = proteinan[!proteinan$sample %in% roi.fail,]

qc.roi.fil     = data.frame("FOV_registration"     = as.numeric(annotation_fil$Fov_counted)/as.numeric(annotation_fil$Fov_count),
                            "Binding_density"      = as.numeric(annotation_fil$BindingDensity), 
                            "Positive_control"     = t(mean(as.numeric(datacounts_fil[grepl("POS", rownames(datacounts_fil)),]))/datacounts_fil[grepl("POS", rownames(datacounts_fil)),]),
                            "Minimum_surface_area" = as.numeric(annotation_fil$AOI_surface_area),
                            "Minimum_nuclei_count" = as.numeric(annotation_fil$AOI_nuclei_count))

####choose normalizing strategy####

Controls = datacounts_fil[c(Negcontrol, housekeepers),] %>% rbind("Area"         = as.numeric(annotation_fil$AOI_surface_area), 
                                                                  "nuclei_count" = as.numeric(annotation_fil$AOI_nuclei_count), 
                                                                  "Countsums"    = as.numeric(colSums(datacounts_fil)))
cor      = cor(t(Controls))

#corplot controls

p = ggcorrplot(cor, outline.col = "white")
p + scale_fill_gradient2(limit = c(round((min(cor)-.1),1),1), 
                         low = "blue", high =  "red", mid = "white", 
                         midpoint = 1-(.5*(1-round((min(cor)-.1),1))))

####scale to positive scaling factor####

dat.fil.pos = t(t(datacounts_fil)*qc.roi.fil$HYB.POS)

####scale to nuclei count or area####

#scaling factors
scaling.facors = data.frame("nuclei"   = gm_mean(as.numeric(annotation_fil$AOI_nuclei_count))/as.numeric(annotation_fil$AOI_nuclei_count),
                             "area"    = gm_mean(as.numeric(annotation_fil$AOI_surface_area))/as.numeric(annotation_fil$AOI_surface_area), 
                             row.names = row.names(qc.roi.fil))

#scale nuclei

dat.fil.nuc  = t(t(dat.fil.pos)*scaling.facors$nuclei)

#scale area

dat.fil.area = t(t(dat.fil.pos)*scaling.facors$area)

#visualization

boxplot(log(dat.fil.pos), cex = .2, pch = 19, las = 2, xaxt = "n", main = "Raw data", xlab = "Regions of Interest")
boxplot(log(dat.fil.nuc), cex = .2, pch = 19, las = 2, xaxt = "n", main = "Scaled for nuclei count", xlab = "Regions of Interest")
boxplot(log(dat.fil.area), cex = .2, pch = 19, las = 2, xaxt = "n", main = "Scaled for area", xlab = "Regions of Interest")

######normalization######
#####negative control plot####

neg.plot = data.frame(t(log(dat.fil.nuc[rownames(dat.fil.nuc) %in% c("Rb IgG", "Ms IgG1", "Ms IgG2a"),])))

ggplot(neg.plot, aes(x = 1:72)) + 
  geom_line(aes(y = Rb.IgG,   colour = "Rb.IgG"),   size = 1) + 
  geom_line(aes(y = Ms.IgG1,  colour = "Ms.IgG1"),  size = 1) + 
  geom_line(aes(y = Ms.IgG2a, colour = "Ms.IgG2a"), size = 1) + 
  ylab("count") + xlab("ROIs") + ggtitle("IgG controls") + labs(colour = "Genes") +
  theme_minimal()
dev.off()

####normalization factors####

neg.norm         = as.data.frame(dat.fil.pos[Negcontrol,])
neg.norm.fact    = mean(apply(neg.norm,  2, gm_mean)) / apply(neg.norm,  2, gm_mean)

####Spike-in normalization####

dat.norm.neg  = t(t(dat.fil.pos)  * neg.norm.fact)

# visualize normalization

boxplot(dat.norm.neg, log = "y", cex = .2, pch = 19, las = 2, xaxt = "n", main = "IgG control normalization", xlab = "Regions of Interest")

####determining background####

bgsd       = apply(dat.fil.pos[Negcontrol,], 2, mean)
cbg        = t(t(dat.fil.pos)/bgsd)

boxplot(cbg, log = "y", cex = .2, pch = 19, las = 2, xaxt = "n", main = "Signal to noise ratio", xlab = "Regions of Interest")

pl = cbg[order(apply(cbg, 1, median)),]

low_targets = data.frame("threshold" = apply(pl, 1, max)>3) %>% 
              rownames_to_column("probes") %>% 
              filter(threshold == F) %>% 
              pull(probes)

myColors <- ifelse(rownames(pl) %in% low_targets, "red", "white")
boxplot(t(pl), log = "y", las = 2, main = "Signal to Noise Ratio", col = myColors, cex = .2, pch = 19, las = 2, cex.axis = .8)
abline(h = c(1,3), col = c("blue","Red"))

cbgfil = cbg[!rownames(cbg) %in% low_targets,]

#Technical replicate

techrep = dat.norm.neg[, grepl("Control", colnames(dat.norm.neg))]
techrep = as.data.frame(log(techrep))

plot(C203_protein_012_Control ~ F206_protein_007_Control, data = techrep, ylab = "Control 1", xlab = "control 2", main = "Technical replicate")
abline(fit <- lm(C203_protein_012_Control ~ F206_protein_007_Control, data = techrep), col='red')
legend("topleft", bty="n", legend=paste("R2 is", format(summary(fit)$adj.r.squared, digits=4)))

saveRDS(dat.norm.neg, file = "dir_processed_NS_protein_data.RDS")
