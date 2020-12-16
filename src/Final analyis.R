#To Do
# Take the ucnultured data an select the markers that are expressed ONLY in Ser
# check how many are in bulk RNA data
# Check how many of them are in later timepoints
{
source("./src/functions.R")
outdir <- "./output/"
outP <- paste0(outdir, "plots/")
dir.create(outP, recursive = TRUE, showWarnings = FALSE)
outT <-  paste0( outdir, "tables/")
dir.create(outT, recursive = TRUE,showWarnings = FALSE )
pval=0.01
}

# markers <- GetGeneList("./data/markers/")

#read the data
uncultured <- readRDS("./data/seurat.uncultured.20200903.rds")
Idents(uncultured) <- uncultured@meta.data$shortName
cultured <- readRDS("./data/all_cells.cultured.annotated.20200826.rds")
Idents(cultured) <- cultured@meta.data$AssignedIdentities

#Diff exp-----


#Get the markers for Ser cluster, Dop, and progenitors
markersUNCSer <- FindMarkers(uncultured, ident.1 = "Glut/Sero N", ident.2 = NULL)
markersUNCSer <- subset(markersUNCSer, markersUNCSer$avg_logFC>0 & markersUNCSer$p_val_adj<pval)
markersUNCDop <- FindMarkers(uncultured, ident.1 = "Dopa N", ident.2 = NULL)
markersUNCDop <- subset(markersUNCDop, markersUNCDop$avg_logFC>0 & markersUNCDop$p_val_adj<pval)
markersUNCPrec <- FindMarkers(uncultured, ident.1  ="Prog", ident.2 =NULL)
markersUNCPrec <- subset(markersUNCPrec, markersUNCPrec$avg_logFC>0 & markersUNCPrec$p_val_adj<pval)


#Select the ones only on Ser
SerMarkers <- setdiff(x= markersUNCSer, y=markersUNCDop, all.x = T,by=0)
SerMarkers <- setdiff(x= SerMarkers, y=markersUNCPrec, all.x = T,by=0)
SerMarkers <- rownames_to_column(SerMarkers, var = "Gene")
write.xlsx(SerMarkers, paste0(outT,"Table_1_SerMarkersUnc.xlsx"),quote=F,rownames=TRUE)
rm(markersUNCDop,markersUNCPrec)

#Bulk counts ----
#Check the expression in the kirkeby data
equivalents <- read.csv("./data/BulkRNA/mart_export.txt", sep = "\t", header = T)
equivalents <- equivalents %>% distinct(Gene.name, .keep_all = T)
equivalents <- column_to_rownames(equivalents, var = "Gene.name")

RPKM1 <- read.csv("./data/BulkRNA/ES_published_samples.txt", sep = "\t", header = T)
RPKM2 <- read.csv("./data/BulkRNA/ES_sequenced_samples.txt", sep = "\t", header = T)
RPKM <- cbind(RPKM1,RPKM2)
rm(RPKM1,RPKM2)

#check on published
datadf <- data.frame(matrix(ncol = 16))
colnames(datadf) <- c("id",colnames(RPKM[,1:15]))
for (id in SerMarkers$Gene){
  row=data.frame(cbind(id, RPKM[equivalents[id,1],1:15]))#Select the published data (first15)
  colnames(row) <- colnames(datadf)
  datadf <- rbind(datadf,row)
}
datadf <- na.omit(datadf)
write.xlsx(datadf, paste0(outT,"Table_2_Counts_RPKM_Published.xlsx"),quotes=FALSE)
rownames(datadf) <- NULL

datadf <- column_to_rownames(datadf, var="id")
datadf$Sum <- rowSums(datadf)
sums <- data.frame(id=rownames(datadf),Total_RPKM=datadf[,"Sum"])
write.xlsx(sums, paste0(outT,"Table_2_Counts_RPKM_Published_Total.xlsx"),quotes=FALSE)


#check on all Kirkeby data
datadf <- data.frame(matrix(ncol = 58))
colnames(datadf) <- c("id",colnames(RPKM))
for (id in SerMarkers$Gene){
  row=data.frame(cbind(id, RPKM[equivalents[id,1],]))
  colnames(row) <- colnames(datadf)
  datadf <- rbind(datadf,row)
}
datadf <- na.omit(datadf)
write.xlsx(datadf, paste0(outT,"Table_2.1_Counts_RPKM_All.xlsx"),quotes=FALSE)
rownames(datadf) <- NULL

datadf <- column_to_rownames(datadf, var="id")
datadf$Sum <- rowSums(datadf)
sums <- data.frame(id=rownames(datadf),Total_RPKM=datadf[,"Sum"])
write.xlsx(sums, paste0(outT,"Table_2.1_Counts_RPKM_All_Total.xlsx"),quotes=FALSE)

rm(row,RPKM,equivalents,datadf, sums)

#Cultured data----

#Check cultured data
markersCULSer <- FindMarkers(cultured, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
markersCULSer <- subset(markersCULSer, markersCULSer$avg_logFC>0 & markersCULSer$p_val_adj<pval)
markersCULDop <- FindMarkers(cultured, ident.1 = "Dopamine Neurons", ident.2 = NULL)
markersCULDop <- subset(markersCULDop, markersCULDop$avg_logFC>0 & markersCULDop$p_val_adj<pval)
markersCULPrec <- FindMarkers(cultured, ident.1  ="Neural Precursor", ident.2 =NULL)
markersCULPrec <- subset(markersCULPrec, markersCULPrec$avg_logFC>0 & markersCULPrec$p_val_adj<pval)


#Select the ones only on Ser
SerMarkersCUL <- setdiff(x= markersCULSer, y=markersCULDop, all.x = T,by=0)
SerMarkersCUL <- setdiff(x= SerMarkersCUL, y=markersCULPrec, all.x = T,by=0)
SerMarkersCUL <- rownames_to_column(SerMarkersCUL, var = "Gene")
write.xlsx(SerMarkersCUL, paste0(outT,"SerMarkersCUL.xlsx"),quote=F,rownames=TRUE)

common <- inner_join(SerMarkers,SerMarkersCUL, by="Gene",suffix=c("_UNC","_CUL"))
write.xlsx(common,paste0(outT, "Table_3_Ser_Markers_UNC_CUL.xlsx"))

rm(markersCULDop, markersCULSer,markersCULPrec, SerMarkersCUL, common,SerMarkers)

#Across timepoint ----

#check the data by timepoints from uncultured to day 14 and 30
d14 <- c("hVM 3D d14","hVM 2D d14")
d30 <- c("hVM 3D d30" , "hVM 2D d30")
cultd14 <- subset(cultured, subset = Group == d14)
cultd30 <- subset(cultured, subset = Group == d30)

rm(d14,d30)

names(markersUNCSer) <- paste0(names(markersUNCSer),"_D0")
markersUNCSer <- rownames_to_column(markersUNCSer, var = "ID")

markersCULSer14 <- FindMarkers(cultd14, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
markersCULSer14 <- subset(markersCULSer14, markersCULSer14$avg_logFC>0 & markersCULSer14$p_val_adj<pval)
names(markersCULSer14) <- paste0(names(markersCULSer14),"_D14")
markersCULSer14 <- rownames_to_column(markersCULSer14, var = "ID")
markersCULSer30 <- FindMarkers(cultd30, ident.1 = "Glutamatergic/Serotoninergic", ident.2 = NULL)
markersCULSer30 <- subset(markersCULSer30, markersCULSer30$avg_logFC>0 & markersCULSer30$p_val_adj<pval)
names(markersCULSer30) <- paste0(names(markersCULSer30),"_D30")
markersCULSer30 <- rownames_to_column(markersCULSer30, var = "ID")



exp_time <- full_join(markersUNCSer,markersCULSer14, by="ID") %>% 
  full_join(., markersCULSer30, by="ID")
write.xlsx(exp_time, paste0(outT,"Table_3.1_markersUNC-TIME_all.xlsx"))

exp_time <- na.omit(exp_time)
write.xlsx(exp_time, paste0(outT,"Table_3.2_markersUNC-TIME_continous.xlsx"))


rm(list=setdiff(ls(), c("outdir","outP","outT")))

#Plotting----
timepoint <- read.xlsx(paste0(outT,"Table_3.1_markersUNC-TIME_all.xlsx"))
rpkm <- read.xlsx(paste0(outT,"Table_2.1_Counts_RPKM_All.xlsx"))
