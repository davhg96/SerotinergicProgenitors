#TO DO
# check diff exp at d0, select specific for Serotoninergic -DONE
# featplot across timepoint -DONE
# Bluk RNA barplot of these selected markers
# lfc line plot across timepoint
{
  source("./src/functions.R")
  outdir <- "./output/"
  outP <- paste0(outdir, "plots/")
  dir.create(outP, recursive = TRUE, showWarnings = FALSE)
  outP2 <- paste0(outP, "FeatPlots/")
  dir.create(outP2, recursive = TRUE, showWarnings = FALSE)
  outT <-  paste0( outdir, "tables/")
  dir.create(outT, recursive = TRUE,showWarnings = FALSE )
  pval=0.01
  rpkm_mean_thrshld <- 10
}


uncultured <- readRDS("./data/seurat.uncultured.20200903.rds")
Idents(uncultured) <- uncultured@meta.data$shortName

cultured <- readRDS("./data/all_cells.cultured.annotated.20200826.rds")
Idents(cultured) <- cultured@meta.data$AssignedIdentities

d14 <- c("hVM 3D d14","hVM 2D d14")
d30 <- c("hVM 3D d30" , "hVM 2D d30")
cultd14 <- subset(cultured, subset = Group == d14)
cultd30 <- subset(cultured, subset = Group == d30)

rm(cultured,d14,d30)

#FeatPlot diff exp D0

timepoint <- read.xlsx(paste0(outT,"Table_3.1_markersUNC-TIME_all.xlsx"))
rownames(timepoint) <- timepoint$ID
rpkm <- read.xlsx(paste0(outT,"Table_2.1_Counts_RPKM_All.xlsx"))
rownames(rpkm) <- rpkm$id
rpkm$mean <- rowMeans(rpkm[,-1]) 
rpkm <- rpkm[order(-rpkm$mean),,drop=FALSE]
rpkm <- rpkm[-c(1),]



for (gene in rpkm$id){
  FeaturePlot(uncultured, features = gene, cols = c("grey","red","black"))
  ggsave(paste0("FeatPlot_",gene,".pdf"), 
         path = outP2,device = "pdf")
}

#manual visualization of plots
marker_genes <- c("CRYBA2","FEV","GCH1","SLC6A4","SLC17A8","TPH2")

outPD14 <- paste0(outP2, "day14/")
dir.create(outPD14)
outPD30 <- paste0(outP2, "day30/")
dir.create(outPD30)

for (gene in marker_genes){
  FeaturePlot(cultd14, features = gene, cols = c("grey","red","black"))
  ggsave(paste0("FeatPlot_D14_",gene,".pdf"), 
         path = outPD14,device = "pdf")
  FeaturePlot(cultd30, features = gene, cols = c("grey","red","black"))
  ggsave(paste0("FeatPlot_D30_",gene,".pdf"), 
         path = outPD30,device = "pdf")
}
 

# Barplot

#All samples
plotdf <-  data.frame(ID=rpkm$id,
                      mean=apply(rpkm[,2:ncol(rpkm)],1,mean),
                      SD=apply(rpkm[,2:ncol(rpkm)],1,sd))
plotdf <- subset(plotdf, rownames(plotdf) %in% marker_genes)

ggplot(data = plotdf, aes(x=ID, y=mean))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.2,
                position=position_dodge(.9)) +
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 8))

ggsave("Barplot_expression_D0_all_samples-FSTL.pdf", path = outP,device = "pdf", height = 15, width = 30, units = "cm" )  

#LFC over time

datatp <- timepoint[marker_genes,c(1,3,8,13)]
datatp <- gather(datatp, key = "day", value = "avg_LFC", avg_logFC_D0,avg_logFC_D14, avg_logFC_D30)
datatp$day <- sapply(strsplit(datatp$day, "_"), tail, 1)
datatp$day <- sapply(substr(datatp$day,2,4),tail,1)

plotdf <- datatp
plotdf$ID <- as.factor(plotdf$ID)
plotdf$day <- as.numeric(plotdf$day)
labels <- data.frame(ID=levels(plotdf$ID), 
                     day=rep(3,length(levels(plotdf$ID))),
                     avg_LFC=rep(4,length(levels(plotdf$ID)))) 
rownames(labels) <- labels$ID
for (gene in labels$ID){
  info <- subset(plotdf, plotdf$ID==gene) %>% na.omit()
  labels[gene, "day"] <- tail(info$day,1)
  labels[gene, "avg_LFC"] <- tail(info$avg_LFC,1)
}

ggplot(plotdf, aes(x=day, y=avg_LFC, group=ID, color=ID))+
  geom_line() +
  geom_point(aes(fill=ID),shape=21, size=2) +
  ggrepel::geom_label_repel(data=labels, 
                            aes(label=ID,
                                x=day,
                                y=avg_LFC),
                            size=5,
                            box.padding   = 0.5, 
                            point.padding = 0.5,
                            segment.color = 'grey50')+
  theme_minimal() +
  theme(legend.position = "none")+
  xlim(0,35)+
  ggtitle("Average logFold change over time")
ggsave(paste0("AVG_LFC_overTime_Manual_pick", ".pdf"), 
       path = outP,device = "pdf", height = 15, width = 30, units = "cm" )
