

{
  source("./src/functions.R")
  outdir <- "./output/"
  outP <- paste0(outdir, "p_markers/")
  dir.create(outP, recursive = TRUE, showWarnings = FALSE)
  outT <- paste0(outdir, "tables/")
  dir.create(outT, recursive = TRUE, showWarnings = FALSE)
}


#Create the DF -----
features <- read.xlsx("./data/samples/markers.xlsx")
RPKM <- read.csv("./data/samplesES_published_samples.txt", sep = "\t", header = T)
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

write.xlsx(plotdf, "./data/samples/plotDF.xlsx" )

#plot ----

df <- data_summary(plotdf, varname ="rpkm", groupnames = c("TH_group","sample_name","data_from"))

ggplot(data = plotdf)+
  