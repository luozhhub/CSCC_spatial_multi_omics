library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)

### Fig2g
sc.metabolism.Seurat<-function (obj, method = "VISION", imputation = F, ncores = 2, 
                                metabolism.type = "KEGG") 
{
  countexp <- obj@assays$Spatial@counts
  countexp <- data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  if (metabolism.type == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (metabolism.type == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- row.names(countexp)
  }
  cat("Start quantify the metabolism activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  }
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  obj@assays$METABOLISM$score <- signature_exp
  obj
}

SM<-readRDS("merged_SM_8samples_merged_SM_1.rds")
DefaultAssay(SM) <- "SCT"
SM1 <- sc.metabolism.Seurat(obj = SM, 
                                   method = "VISION",     
                                   imputation = F, ncores = 2, 
                                   metabolism.type = "KEGG")

metabolism.matrix <- SM1@assays$METABOLISM$score
metabolism.matrix
SM1_score<-as.data.frame(t(metabolism.matrix))
rownames(SM1_score)<-colnames(SM1)
SM1<-AddMetaData(SM1,metadata =SM1_score)

input.pathway <- rownames(SM1@assays[["METABOLISM"]][["score"]])[41]
col = c('1'="#7fc97f", '2'="#beaed4", '3'="#fdc086", '4'="#386cb0", '5'="#e0abcf", '6'="#a34e3b", '7'="#efbd25", '8'="#1b9e77", '9'="#d95f02", '10'="#6c67a5", '11'="#d01b2a", '12'="#43acde", '13'="#666666")
cluster_order <- c(5, 6, 7, 10, 1, 2, 3, 4, 8, 9, 11, 12, 13)
SM1$cluster <- factor(SM1$integrated_snn_res.0.5, levels = as.character(cluster_order))
input.pathway <-input.pathway
input.parameter <-'cluster'
metadata <- SM1@meta.data
metabolism.matrix <- SM1@assays$METABOLISM$score

metadata[, input.parameter] <- as.character(metadata[, input.parameter])

metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, ])

gg_table <- c()
for (i in 1:length(input.pathway)) {
  gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], input.pathway[i], metabolism.matrix_sub[, i]))
}
gg_table <- data.frame(gg_table)
colnames(gg_table) <- c("cluster", "pathway", "score")
gg_table[, 3] <- as.numeric(as.character(gg_table[, 3]))
p <- ggplot(gg_table, aes(x = cluster, y = score, fill = cluster)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(limits = as.character(cluster_order)) +  # 按指定顺序显示cluster
  labs(x = "Cluster", y = "Metabolic Score") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank()   # 去除次要网格线
  ) +
  facet_wrap(~ pathway, scales = "free") +
  scale_fill_manual(values = col)

#### Fig2h
plot_gene = function (cluster){
  p1 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM21',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.4,images = 'SM25',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM27',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM32',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p5 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM34',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p6 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM35',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p7 = SpatialPlot(SM1, features  = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'N1',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p8 = SpatialPlot(SM1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'N5',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol  = 8,nrow = 1,common.legend = T,legend = 'right')
}

limit=c(-0.4, 0.4)
breaks = c(-1,1, 3)
p1<-plot_gene(cluster='Arginine and proline metabolism')
p1

#### Fig2i
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/N1-1-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/N5-3-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM25-2-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM32-2-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM34-1-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM35-2-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg

path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM-27-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
image(tiny_in, mz=136.04813, colorscale=jet.colors, smooth.image="gaussian")#neg



