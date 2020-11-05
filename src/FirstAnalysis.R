library(Seurat)
library(ggplot2)
library(openxlsx)

outdir <- "./output/"
outP <- paste0(outdir, "plots/")
dir.create(outP, recursive = TRUE, showWarnings = FALSE)
outT <-  paste0( outdir, "tables/")
dir.create(outT, recursive = TRUE,showWarnings = FALSE )


uncultured <- readRDS("./data/seurat.uncultured.20200903.rds")
Idents(uncultured) <- uncultured@meta.data$shortName


cultured <- readRDS("./data/all_cells.cultured.annotated.20200826.rds")
Idents(cultured) <- cultured@meta.data$AssignedIdentities
# unique(cultured@meta.data$Group) Guives the 2D and 3D groups
# groups <- c("hVM 3D d14" ,"hVM 3D d30", "hVM 2D d14", "hVM 2D d30")


# day to culture analysis -----

d143D <- subset(cultured, Group=="hVM 3D d14")
d142D <- subset(cultured, Group=="hVM 2D d14")
d303D <- subset(cultured, Group=="hVM 3D d30")
d302D <- subset(cultured, Group=="hVM 2D d30")



markersUNC <- FindMarkers(uncultured, ident.1 = "Glut/Sero N", ident.2 = NULL, min.diff.pct = 0.1 )
markersUNC <- markersUNC[order(markersUNC$pct.1, decreasing =TRUE),]
markersCUL <- FindMarkers(cultured, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL, min.diff.pct = 0.1)
markersCUL <- markersCUL[order(markersCUL$pct.1, decreasing =TRUE),]
markersd143D <- FindMarkers(d143D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL, min.diff.pct = 0.1)
markersd143D <- markersd143D[order(markersd143D$pct.1, decreasing =TRUE),]
markersd142D <- FindMarkers(d142D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL, min.diff.pct = 0.1)
markersd142D <- markersd142D[order(markersd142D$pct.1, decreasing =TRUE),]
markersd303D <- FindMarkers(d303D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL, min.diff.pct = 0.1)
markersd303D <- markersd303D[order(markersd303D$pct.1, decreasing =TRUE),]
markersd302D <- FindMarkers(d302D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL, min.diff.pct = 0.1)
markersd302D <- markersd302D[order(markersd302D$pct.1, decreasing =TRUE),]


write.xlsx(markersUNC,paste0(outT,"Uncultured.xlsx"))
write.xlsx(markersCUL, paste0(outT, "Cultured.xlsx"))
write.xlsx(markersd142D, paste0(outT, "d142D"))
write.xlsx(markersd143D, paste0(outT, "d143D"))
write.xlsx(markersd302D, paste0(outT, "d302D"))
write.xlsx(markersd303D, paste0(outT, "d303D"))





# pltos ====
common <- subset(markersUNC, rownames(markersCUL) %in% rownames(markersUNC))
expressed <- subset(common, common$p_val_adj<0.001 & common$avg_logFC>0)

genes <- rownames(expressed)
pUNC <- FeaturePlot(uncultured, features = genes, reduction = "umap", cols = c("#ababab","#ff0000"))
ggsave("FPlotUncultured.pdf", plot = pUNC, path = outP, dpi = 600, height = 25, units = "cm")
pCUL <- FeaturePlot(cultured, features = genes, reduction = "umap", cols = c("#ababab","#ff0000"))
ggsave("FplotCultured.pdf", plot = pCUL, path = outP, dpi = 600, height = 25, units = "cm")

vUNC <- VlnPlot(uncultured, features = genes,pt.size = 0)
ggsave("VPlotUncultured.pdf", plot = vUNC, path = outP, dpi = 600, height = 25, units = "cm")
vCUL <- VlnPlot(cultured, features = genes,pt.size = 0)
ggsave("VPlotCultured.pdf", plot = vCUL, path = outP, dpi = 600, height = 25, units = "cm")

dotUNC <- DotPlot(uncultured, features = genes)+coord_flip()
ggsave("DotUncultured.pdf", plot = dotUNC, path = outP, dpi = 600)
dotCUL <- DotPlot(cultured, features = genes)+coord_flip()
ggsave("DotCultured.pdf", plot = dotCUL, path = outP, dpi = 600)

# rerunning UMAP----

UncSer <- subset(uncultured, subset = shortName== "Glut/Sero N")
# UncNoSer <- subset(uncultured, subset = shortName== "Glut/Sero N", invert=TRUE)
CultSer <- subset(cultured,  idents = "Glutamatergic/Serotoninergic")
# CulNoSer <- subset(cultured,  idents = "Glutamatergic/Serotoninergic", invert=TRUE)


