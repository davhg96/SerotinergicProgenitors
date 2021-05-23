library(Seurat)
library(tidyverse)
library(openxlsx)
library(slingshot)
library(scales)



compute_silhouette_scores <- function(seurat,sil.from=.1, sil.to=1, by=.1,plot=T){
  library(cluster)
  library(ggplot2)
  library(cowplot)
  sil.results <- data.frame()
  for(resolution in seq(sil.from,sil.to,by = by)){
    seurat<-FindClusters(seurat,resolution=resolution)
    dists <- dist(Embeddings(seurat,reduction = "umap"))
    combined.fac <- factor(paste0(seurat$seurat_clusters))
    clusters<-as.integer(seurat$seurat_clusters)
    sil <- silhouette(clusters, dist = dists)
    if(plot){
      filename=paste0("silhouette.res=",resolution,".png")
      png(filename)
      plot(sil, border = NA)
      dev.off()
      cat("Plotted to: ",filename,"\n")
    }
    sil.results <- rbind(sil.results,c(resolution,mean(sil[,3]),length(levels(seurat$seurat_clusters))))
  }
  colnames(sil.results)<-c("Resolution","Average Si","nClusters")
  plot_grid(
    ggplot(sil.results,aes(x=factor(Resolution),y=`Average Si`))+geom_bar(stat="identity")+theme_cowplot()+xlab(""),
    ggplot(sil.results,aes(x=factor(Resolution),y=nClusters))+geom_bar(stat="identity")+theme_cowplot()+xlab("Resolution"),ncol=1)
  ggsave("silhouette_results.pdf",w=8,h=6)
  return(sil.results)
}




#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}




GetGeneList <- function(directory){
  FileList <- list.files(directory, pattern = ".txt") # ls the directory
  ResList <- list() #placeholder
  for (file in FileList){
    path <- paste0(directory,file)
    Filedata <- as.list(read.table(path, sep = "\n", header = TRUE)) #Extract the info
    ResList[[file]] <- Filedata[[1]] #Add the info to the result object
    
    
  }
  
  return(ResList)
}


# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
