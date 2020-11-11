
{
source("./src/functions.R")
outdir <- "./output/"
outP <- paste0(outdir, "plots/")
dir.create(outP, recursive = TRUE, showWarnings = FALSE)
outT <-  paste0( outdir, "tables/")
dir.create(outT, recursive = TRUE,showWarnings = FALSE )
}
markers <- GetGeneList("./data/markers/")

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



markersUNC <- FindMarkers(uncultured, ident.1 = "Glut/Sero N", ident.2 = NULL)
# markersUNC <- markersUNC[order(markersUNC$pct.1, decreasing =TRUE),]
markersCUL <- FindMarkers(cultured, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
# markersCUL <- markersCUL[order(markersCUL$pct.1, decreasing =TRUE),]
markersd143D <- FindMarkers(d143D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
# markersd143D <- markersd143D[order(markersd143D$pct.1, decreasing =TRUE),]
markersd142D <- FindMarkers(d142D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
# markersd142D <- markersd142D[order(markersd142D$pct.1, decreasing =TRUE),]
markersd303D <- FindMarkers(d303D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
# markersd303D <- markersd303D[order(markersd303D$pct.1, decreasing =TRUE),]
markersd302D <- FindMarkers(d302D, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
# markersd302D <- markersd302D[order(markersd302D$pct.1, decreasing =TRUE),]


write.xlsx(markersUNC,paste0(outT,"Uncultured.xlsx"),row.names=T)
write.xlsx(markersCUL, paste0(outT, "Cultured.xlsx"),row.names=T)
write.xlsx(markersd142D, paste0(outT, "d142D.xlsx"),row.names=T)
write.xlsx(markersd143D, paste0(outT, "d143D.xlsx"),row.names=T)
write.xlsx(markersd302D, paste0(outT, "d302D.xlsx"),row.names=T)
write.xlsx(markersd303D, paste0(outT, "d303D.xlsx"),row.names=T)



# pltos ====
common <- subset(markersUNC, rownames(markersCUL) %in% rownames(markersUNC))
write.xlsx(common, paste0(outT, "common.xlsx"),row.names=T)
expressed <- subset(common, common$p_val_adj<0.001 & common$avg_logFC>0)
write.xlsx(expressed, paste0(outT, "Common_Expressed.xlsx"),row.names=T)

genes <- rownames(expressed)
pUNC <- FeaturePlot(uncultured, features = genes, reduction = "umap", cols = c("grey", "red", "black"))
ggsave("FPlotUncultured.pdf", plot = pUNC, path = outP, dpi = 600, width = 40, height = 30, units = "cm")
pCUL <- FeaturePlot(cultured, features = genes, reduction = "umap", cols = c("grey", "red", "black"))
ggsave("FplotCultured.pdf", plot = pCUL, path = outP, dpi = 600,  width = 40, height = 30, units = "cm")

vUNC <- VlnPlot(uncultured, features = genes,pt.size = 0)
ggsave("VPlotUncultured.pdf", plot = vUNC, path = outP, dpi = 600, width = 40, height = 30, units = "cm")
vCUL <- VlnPlot(cultured, features = genes,pt.size = 0)
ggsave("VPlotCultured.pdf", plot = vCUL, path = outP, dpi = 600, width = 40, height = 30, units = "cm")

dotUNC <- DotPlot(uncultured, features = genes)+coord_flip()+ theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("DotUncultured.pdf", plot = dotUNC, path = outP, dpi = 600)
dotCUL <- DotPlot(cultured, features = genes)+coord_flip()+ theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("DotCultured.pdf", plot = dotCUL, path = outP, dpi = 600)


#Plot ser/glut markers
for(name in names(markers)){
FeaturePlot(uncultured, features = markers[[name]], reduction = "umap", cols = c("grey", "red", "black"))
ggsave(paste0("FeaturePlot_Uncultured_",name, ".pdf"), path = outP, dpi = 600, width = 50, height = 50, units = "cm")
FeaturePlot(cultured, features = markers[[name]], reduction = "umap", cols = c("grey", "red", "black"))
ggsave(paste0("FeaturePlot_Cultured_",name, ".pdf"), path = outP, dpi = 600,  width = 50, height = 50, units = "cm")
rm(name)
}

# integrating clusters ----
UncSer <- subset(uncultured, subset = shortName== "Glut/Sero N")
DefaultAssay(UncSer) <- "RNA"
UncSer[["State"]] <- rep("Uncultured", nrow(UncSer@meta.data))

CultSer <- subset(cultured,  idents = "Glutamatergic/Serotoninergic")
DefaultAssay(CultSer) <- "RNA"
CultSer[["State"]] <- rep("Cultured", nrow(CultSer@meta.data))

Ser.List <-  list(Uncultured=UncSer, Cultured=CultSer)
for (i in 1:length(Ser.List)){
  Ser.List[[i]] <- NormalizeData(Ser.List[[i]], verbose = FALSE)
  Ser.List[[i]] <- FindVariableFeatures(Ser.List[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

reference.list <- Ser.List[c("Uncultured","Cultured")]
Ser.Anchors <- FindIntegrationAnchors(object.list = Ser.List, dims = 1:30)
Ser.Integrated <- IntegrateData(anchorset = Ser.Anchors, dims = 1:30)
rm(Ser.List,reference.list,Ser.Anchors)

Ser.Integrated@meta.data <- Ser.Integrated@meta.data[,colSums(is.na(Ser.Integrated@meta.data))==0]



Ser.Integrated <- FindVariableFeatures(object = Ser.Integrated)
Ser.Integrated <- ScaleData(Ser.Integrated)
Ser.Integrated <- RunPCA(Ser.Integrated, npcs = 30)
Ser.Integrated <- FindNeighbors(Ser.Integrated, reduction = "pca", dims = 1:30)
Ser.Integrated <- FindClusters(Ser.Integrated, resolution = 0.4)
Ser.Integrated <- RunUMAP(Ser.Integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(Ser.Integrated, reduction = "umap", group.by = "State")
p2 <- DimPlot(Ser.Integrated, reduction = "umap")
p1/p2
ggsave("IntegratedUMAP0.4.pdf",path = outP, device = "pdf", height = 25, width = 20, dpi = 600, units = "cm")
saveRDS(Ser.Integrated,"./data/Glut_Ser_Integrated0.4res.rds")

markerSer <- FindMarkers(Ser.Integrated,ident.1 = 0, ident.2 = NULL )
markerSer <- markerSer[order(markerSer$pct.1, decreasing =TRUE),]
markerSer1 <- FindMarkers(Ser.Integrated,ident.1 = 1, ident.2 = NULL )
markerSer1 <- markerSer1[order(markerSer$pct.1, decreasing =TRUE),]
markerSer2 <- FindMarkers(Ser.Integrated,ident.1 = 2, ident.2 = NULL )
markerSer2 <- markerSer[order(markerSer$pct.1, decreasing =TRUE),]

FeaturePlot(Ser.Integrated, features = "FEV", reduction = "umap")

# rerunning UMAP----

SerClust <- readRDS("./data/Glut_Ser_integrated.rds")

silhouette <- compute_silhouette_scores(SerClust)

#seems like the optimum is 0.1 or 0.3

optRes=c(0.1,0.3)
for (res in optRes){
SerClust <- FindNeighbors(SerClust, reduction = "pca", dims = 1:30)
SerClust <- FindClusters(SerClust, resolution = res)
SerClust <- RunUMAP(SerClust, reduction = "pca", dims = 1:30)
p1 <- DimPlot(SerClust, reduction = "umap", group.by = "State")
p2 <- DimPlot(SerClust, reduction = "umap")
p1/p2
ggsave(paste0("IntegratedOPTUMAP",res,".pdf"),path = outP, device = "pdf", height = 25, width = 20, dpi = 600, units = "cm")
saveRDS(SerClust,paste0("./data/Glut_Ser_Integrated",res,"res.rds"))
}
rm(optRes, res, p1, p2 )

## Plot specific glut/ser markers
SerClust <- readRDS("./data/Glut_Ser_integrated0.3res.rds")
for(name in names(markers)){
 FeaturePlot(SerClust, features = markers[[name]], reduction = "umap", cols = c("grey", "red", "black") ) 
 ggsave(paste0("FeaturePlot_Ser-Glut_",name,".pdf"), path = outP, dpi = 600,  width = 50, height = 50, units = "cm")
rm(name)
 }

markers0 <- FindMarkers(SerClust, ident.1 = 0, ident.2 = NULL)
markers0 <- subset(markers0, markers0$avg_logFC>0)
write.xlsx(markers0,paste0(outT,"Markers_Clust0.xlsx"),row.names=T)
markers2 <- FindMarkers(SerClust, ident.1 = 2, ident.2 = NULL)
markers2 <- subset(markers2, markers2$avg_logFC>0)
write.xlsx(markers2,paste0(outT,"Markers_Clust2.xlsx"),row.names=T)
markers4 <- FindMarkers(SerClust, ident.1 = 4, ident.2 = NULL)
markers4 <- subset(markers4, markers4$avg_logFC>0)
write.xlsx(markers4,paste0(outT,"Markers_Clust4.xlsx"),row.names=T)



# lineage analysis ----
sds <- slingshot(Embeddings(SerClust, "umap"), clusterLabels = SerClust$seurat_clusters, 
                 start.clus = NULL, stretch = 2)


cell_colors_clust <- cell_pal(SerClust$seurat_clusters, hue_pal())

pdf(paste0(outP, "lineage.pdf"))
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()
