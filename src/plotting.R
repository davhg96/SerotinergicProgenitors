

{
  source("./src/functions.R")
  outdir <- "./output/"
  outP2 <- paste0(outdir, "p_markers/")
  dir.create(outP2, recursive = TRUE, showWarnings = FALSE)
  outT <- paste0(outdir, "tables/")
  dir.create(outT, recursive = TRUE, showWarnings = FALSE)
}

#data dfs -----
uncultured <- read.xlsx("./output/tables/uncultured.xlsx", rowNames = T)
uncultured <- subset(uncultured, uncultured$avg_logFC>0)
uncultured <- rownames_to_column(uncultured, var = "ID")
write.xlsx(uncultured, file = "./output/tables/uncultured_expressed.xlsx", row.names=T)

cultured <-  read.xlsx("./output/tables/cultured.xlsx", rowNames = T)
cultured <- subset(cultured, cultured$avg_logFC>0)
cultured <- rownames_to_column(cultured, var = "ID")
write.xlsx(cultured, file= "./output/tables/cultured_expressed.xlsx", row.names=T)

common <- merge(uncultured, cultured, by="ID")
common <- subset(common, common$p_val_adj.x<0.001 & common$avg_logFC.x>0)

all.df <- data.frame(
  ID = c(common$ID, cultured$ID, uncultured$ID),
  ESM_ID = rep(3, nrow(common) + nrow(cultured) + nrow(uncultured)),
  data_from = c(
    rep("common", nrow(common)),
    rep("cultured", nrow(cultured)),
    rep("uncultured", nrow(uncultured))
  )
)

write.xlsx(all.df, file = "./data/samples/all_features.xlsx")

#Create the DF -----
features <- read.xlsx("./data/samples/all_features.xlsx")
RPKM <- read.csv("./data/samples/ES_published_samples.txt", sep = "\t", header = T)
equivalents <- read.csv("./data/samples/mart_export.txt", sep = "\t", header = T)
equivalents <- equivalents %>% distinct(Gene.name, .keep_all = T)
equivalents <- column_to_rownames(equivalents, var = "Gene.name")


for (row in 1:nrow(features)) {
  features[row, 2] <- equivalents[features[row, 1], 1]
}
rm(row)

samples <- read.xlsx("./data/samples/samples.xlsx")

# build the dataframe

plotdf <- data.frame(
  TH_group = c(
    rep("TH_high", nrow(features) * colSums(!is.na(samples))[1]),
    rep("TH_low", nrow(features) * colSums(!is.na(samples))[2])),
  sample_name = c(
    rep(samples[!is.na(samples$TH_high), 1], each = nrow(features)),
    rep(samples[!is.na(samples$TH_low), 2], each = nrow(features))),
  gene = rep(features$ID, sum(colSums(!is.na(samples)))),
  esmID = rep(features$ESM_ID, sum(colSums(!is.na(samples)))),
  data_from = rep(features$data_from, sum(colSums(!is.na(samples)))),
  rpkm = rep(3, nrow(features) * sum(colSums(!is.na(samples))))
)


# populate rpkm

for (n in 1:nrow(plotdf)){
  plotdf[n,"rpkm"] <- RPKM[plotdf[n,"esmID"],plotdf[n,"sample_name"]]
}
rm(n)

write.xlsx(plotdf, "./data/samples/plotDF_all.xlsx" )






#plot ----

datadf <- read.xlsx("./data/samples/plotDF_all.xlsx")

plotdf <- 
  datadf %>% filter(!is.na(rpkm)) %>%
  group_by(TH_group,gene,data_from) %>%
  summarize(meanrpkm=mean(rpkm),
            sd = sd(rpkm))
write.xlsx(plotdf, file = paste0(outP2,"plotdf_all.xlsx"))


ggplot(data = plotdf, aes(x=gene, y=meanrpkm, fill=TH_group))+
  facet_grid(rows = vars(data_from))+
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=meanrpkm-sd, ymax=meanrpkm+sd), width=.2,
                position=position_dodge(.9)) +
  theme(text = element_text(size=20),
    axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 25))
  
ggsave("barplot_expression_all.pdf", path = outP2,device = "pdf", height = 30, width = 75, units = "cm" )
  

plotdf$data_from <- as.factor(plotdf$data_from)
for(group in levels(plotdf$data_from)){
  data <- subset(plotdf, plotdf$data_from == group)
  ggplot(data = data, aes(x=gene, y=meanrpkm, fill=TH_group))+
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=meanrpkm-sd, ymax=meanrpkm+sd), width=.2,
                  position=position_dodge(.9)) +
    theme(text = element_text(size=20),
      axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5,size = 25 ))
  figname <- paste0("barplot_expression_",group,".pdf")
  ggsave(figname, path = outP2,device = "pdf", height = 15, width = 75, units = "cm" )
  
}



