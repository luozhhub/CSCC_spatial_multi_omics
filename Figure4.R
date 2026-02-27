library(Seurat)
library(magrittr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(patchwork)
library(ggtree)
library(BiocGenerics)
library(readr)
library(rtracklayer)
library(infercnv)
library(phylogram)
library(utils)
library(dendextend)
library(assertthat)
library(reticulate)
library(openxlsx)
library(scatterpie)
library(cowplot)
library(stats)
library(quadprog)
library(data.table)
library(Rfast)
library(ggrepel)
library(tibble)
library(clusterProfiler)
library(utils)
library(org.Hs.eg.db)
library(Cottrazm

#####  Fig 4a
#SM35
CSCC<-readRDS("/media/desk16/tly0202/TJ/all_sample_SME_seurat.rds.gz")
SM35 = CSCC[[10]]
SM35 <- SCTransform(SM35, assay = "Spatial", verbose = FALSE)
tumor_SM35 = subset(SM35,Morph_snn_res.1.5 %in% c(5,6,13,9))
DefaultAssay(tumor_SM35) <- "SCT"
Idents(tumor_SM35) = "Morph_snn_res.1.5"
SM35_4_colors = c("#2C9F2D","#F77F0E","#D62728","#2C77B4")
SM27_3_colors = c("#E04D50", "#4374A5", "#F08A21")
SpatialPlot(tumor_SM35,cols=SM35_4_colors)


#SM27
SM27 = CSCC[[5]]
SM27 <- SCTransform(SM27, assay = "Spatial", verbose = FALSE)
tumor_SM27 = subset(SM27,Morph_snn_res.1.5 %in% c(2,9,0))
DefaultAssay(tumor_SM27) <- "SCT"
Idents(tumor_SM27) = "Morph_snn_res.1.5"
SpatialPlot(tumor_SM27,cols=SM27_3_colors)


##### Fig 4b
#data select
gene_df = read.csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/Cottrazm/CSCC_gene_set.txt", sep="\t")
gene_df = gene_df %>% group_by(cluster) %>% top_n(200, avg_logFC)
seu = tumor_SM35
group="CIN"
levels = c(13, 5, 6, 9)
one_features = list(gene_df$Gene[gene_df$cluster == group])
one_scores = AddModuleScore(object = seu, features = one_features, name = group)
#head(one_scores@meta.data)
#plot
y_var = paste(group, "1", sep="")
df_one_scores = one_scores@meta.data[,c("Morph_snn_res.1.5", y_var)]
colnames(df_one_scores) = c("new_cluster", "scores")
df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
p = ggplot(df_one_scores, aes(x = new_cluster, y=scores, fill=new_cluster)) +
    geom_violin(trim=FALSE, color="darkred")+
    #stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="pointrange") +
    stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="crossbar", width=0.2) +
    theme_classic() +
    labs(title=paste("SM35", group, "gene set score"),x="Clusters", y = "Score") +
    theme(
      plot.margin = margin(t = 0.3, r = 0.3,b = 0.3,l = 0.3 ,unit = "cm"),
      #plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
      plot.title = element_text(size=16, hjust = 1.0),
      axis.title = element_text(size=16, colour = 'black'),
      axis.text.y = element_text(size = 12,colour = 'black'),
      axis.text.x = element_text(size = 14,colour = "black",angle = 0,hjust = 1),
      axis.ticks.length = unit(2,"mm")) +
    scale_fill_manual(values = SM35_4_colors) +
    NoLegend()

seu = tumor_SM27
group="CIN"
levels = c(2, 9, 0)
one_features = list(gene_df$Gene[gene_df$cluster == group])
one_scores = AddModuleScore(object = seu, features = one_features, name = group)
#head(one_scores@meta.data)
#plot
y_var = paste(group, "1", sep="")
df_one_scores = one_scores@meta.data[,c("Morph_snn_res.1.5", y_var)]
colnames(df_one_scores) = c("new_cluster", "scores")
df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
p = ggplot(df_one_scores, aes(x = new_cluster, y=scores, fill=new_cluster)) +
    geom_violin(trim=FALSE, color="darkred")+
    #stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="pointrange") +
    stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="crossbar", width=0.2) +
    theme_classic() +
    labs(title=paste("SM27", group, "gene set score"),x="Clusters", y = "Score") +
    theme(
      plot.margin = margin(t = 0.3, r = 0.3,b = 0.3,l = 0.3 ,unit = "cm"),
      #plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
      plot.title = element_text(size=16, hjust = 1.0),
      axis.title = element_text(size=16, colour = 'black'),
      axis.text.y = element_text(size = 12,colour = 'black'),
      axis.text.x = element_text(size = 14,colour = "black",angle = 0,hjust = 1),
      axis.ticks.length = unit(2,"mm")) +
    scale_fill_manual(values = SM27_3_colors) +
    NoLegend()
    
#### Fig4c
gene<-c('SPRR3','SPRR2A','SPRR2E','CNFN','S100A7','SPRR1B','SPRR2D','ANXA1','KRT6B','RHCG','EMP1','CD24','ZFP36','KRT16','FGFBP1','SFN','NDUFA4L2','NTS','SERPINB4','MT2A','MT1E','MT1G','KRT15' ,'MIR205HG','TNC','COL18A1','THBS2','COL7A1','PIK3R1','MMP11','STMN1','DPYSL3','IFITM3','GPC3','COL16A1','PTN')
DoHeatmap(object = tumor_SM35, features = gene)
DoHeatmap(object = tumor_SM27, features = gene)

    
#### Fig4e
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
DefaultAssay(tumor_SM35) <- "SCT"
Idents(tumor_SM35) = "layer"
data.markers <- FindAllMarkers(tumor_SM35, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers = data.markers %>% group_by(cluster) %>% arrange(cluster, avg_log2FC)
write.csv(data.markers, "SM35_DEG.csv", row.names = F, quote = F)
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
DefaultAssay(tumor_SM27) <- "SCT"
Idents(tumor_SM27) = "layer"
data.markers <- FindAllMarkers(tumor_SM27, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers = data.markers %>% group_by(cluster) %>% arrange(cluster, avg_log2FC)
write.csv(data.markers, "SM27_DEG.csv", row.names = F, quote = F)

##pathway generate
library(clusterProfiler)
library(org.Hs.eg.db)

data(geneList, package="DOSE")
data.markers = read.csv("SM27_DEG.csv")
GO_df_SM27 = map(c("Layer1", "Layer2", "Layer3"), function(x){
    data.markers_sub = data.markers[data.markers$cluster == x,]
    DEG.gene_up = as.character(data.markers_sub$gene) 
    DEG.entrez_id_up = mapIds(x = org.Hs.eg.db,
                          keys = DEG.gene_up,
                          keytype = "SYMBOL",
                          column = "ENTREZID")

    ego <- enrichGO(gene          = DEG.entrez_id_up,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
    ego@result$Layer = x
    ego@result$Sample = "SM27"
    return(ego@result)
}) %>% dplyr::bind_rows()

data.markers = read.csv("SM35_DEG.csv")#文件名在前面存的时候写错了
GO_df_SM35 = map(c("Layer1", "Layer2", "Layer3"), function(x){
    data.markers_sub = data.markers[data.markers$cluster == x,]
    DEG.gene_up = as.character(data.markers_sub$gene) 
    DEG.entrez_id_up = mapIds(x = org.Hs.eg.db,
                          keys = DEG.gene_up,
                          keytype = "SYMBOL",
                          column = "ENTREZID")

    ego <- enrichGO(gene          = DEG.entrez_id_up,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
    ego@result$Layer = x
    ego@result$Sample = "SM35"
    return(ego@result)
}) %>% dplyr::bind_rows()

GO_df = rbind(GO_df_SM27, GO_df_SM35)
write.csv(GO_df, "layer_GO.csv", row.names = F, quote = F)
readr::write_rds(GO_df, "layer_GO.rds.gz", compress = "gz")

GO_df <- readr::read_rds("layer_GO.rds.gz")
selected_pathways = c(
    "keratinocyte differentiation",
    "epidermal cell differentiation",
    "skin development",
    "epidermis development",
    "peptide cross-linking",

    "response to hypoxia",
    "response to decreased oxygen levels",
    "response to fatty acid",
    "response to copper ion",
    
    "oxidative phosphorylation",
    "aerobic respiration",
    "ATP synthesis coupled electron transport",
    "proton motive force-driven ATP synthesis",
    "electron transport chain",

    "integrin-mediated signaling pathway",
    "extracellular matrix organization",
    "transforming growth factor beta production",
    "basement membrane organization"
)
selected_pathways = selected_pathways

alldf = GO_df[GO_df$Description %in% selected_pathways,]
alldf$Layer_F <- apply( alldf[ , c("Layer", "Sample")] , 1 , paste , collapse = "_" )

alldf$Description = factor(alldf$Description, levels = c(selected_pathways))
p = ggplot(alldf,aes(Layer_F, Description))+
    theme_classic() + 
    geom_point(aes(size=Count,color=-log10(p.adjust),alpha=0.5))+
    scale_size(range = c(1, 6))+
   
    scale_color_gradient(low="blue",high="red", limits=c(0,5), oob = scales::squish)+
    labs(x="",y="", title = "GO BP pathways")+
    theme(
          legend.position = "right",
          legend.key.height = unit(8, "pt"),
          legend.key.width = unit(8, "pt"),
          axis.line=element_line(size=1),
          axis.text.y = element_text(size=12,face = "plain",colour = "#101010c2"),
          axis.text.x = element_text(size=12,face = "plain",colour = "#101010c2", angle = 45,hjust = 1),
          axis.title = element_text(size=20,face = "plain"),
          plot.margin = unit(c(1,1,1,1),"cm"))
pdf("./final_figures/Figure_4e.pdf", height=7, width=7)
p
dev.off()


##### Fig4f
library("cancersea")
data('available_pathways')
available_pathways
pathways <- list(angiogenesis=angiogenesis, apoptosis=apoptosis, cell_cycle=cell_cycle, differentiation=differentiation, dna_damage=dna_damage,
                dna_repair=dna_repair, emt=emt, hypoxia=hypoxia, inflammation=inflammation, invasion=invasion, metastasis=metastasis, proliferation=proliferation,
                quiescence=quiescence, stemness=stemness)
levels = c("Layer1", "Layer2", "Layer3")
all_groups = names(pathways)
all_groups = all_groups[c(3,6,10,8)]


p_list = purrr::map(1:length(all_groups), function(x){
            group = all_groups[x]
            
            #SM35
            df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
            seu = tumor_SM35
            one_features = list(pathways[[group]]$symbol)
            one_scores = AddModuleScore(object = seu, features = one_features, name = group, assay="SCT")
            y_var = paste(group, "1", sep="")
            df_one_scores = one_scores@meta.data[,c("orig.ident", "Morph_snn_res.1.5", y_var)]
            colnames(df_one_scores) = c("Sample", "new_cluster", "scores")
            df_one_scores$new_cluster = df_layer$Layer[match(df_one_scores$new_cluster, df_layer$Oringe_cluster)]
            df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
            df_one_scores_SM35 = df_one_scores
            
            #SM27
            df_layer_27 = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
            seu = tumor_SM27
            one_features = list(pathways[[group]]$symbol)
            one_scores = AddModuleScore(object = seu, features = one_features, name = group, assay="SCT")
            y_var = paste(group, "1", sep="")
            df_one_scores = one_scores@meta.data[,c("orig.ident", "Morph_snn_res.1.5", y_var)]
            colnames(df_one_scores) = c("Sample", "new_cluster", "scores")
            df_one_scores$new_cluster = df_layer_27$Layer[match(df_one_scores$new_cluster, df_layer_27$Oringe_cluster)]
            df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
            df_one_scores_SM37 = df_one_scores

            df_one_scores = rbind(df_one_scores_SM35, df_one_scores_SM37)

            p = ggplot(df_one_scores, aes(x = new_cluster, y=scores, fill=Sample)) +
                geom_violin(trim=TRUE)+
                stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="crossbar", width=0.2, position=position_dodge(0.9)) +
                theme_classic() +
                labs(title="",x="", y = paste(group, "set score")) +
                theme(
                  plot.margin = margin(t = 0.3, r = 0.3,b = 0.3,l = 0.3 ,unit = "cm"),
                  #plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
                  axis.title = element_text(size=22, colour = 'black'),
                  axis.text.y = element_text(size = 18,colour = 'black'),
                  axis.text.x = element_text(size = 20,colour = "black",angle = 45,hjust = 1),
                  axis.ticks.length = unit(2,"mm")
                  ) +
                scale_fill_manual(values = c("#E64B35", "#3C5488"))

            return(p)           
            
        })



library(cowplot)

legend <- get_legend(p_list[[1]])
all_niches <- plot_grid(p_list[[1]] + theme(legend.position="none"),
                        p_list[[2]]+ theme(legend.position="none"),
                        p_list[[3]]+ theme(legend.position="none"),
                        p_list[[4]]+ theme(legend.position="none"),
                        ncol = 4, align = "hv")
all_niches = plot_grid(all_niches, legend, rel_widths = c(4, .4))


pdf("./final_figures/Figure_4f.pdf", width=20, height=5) 
all_niches
dev.off()

#### EMT score
gene<-c("MT2A","NME2","IFITM3")
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
seu = tumor_SM35
one_scores = AddModuleScore(object = seu, features = list(gene), name = x)
y_var = paste(x, "1", sep="")
df_one_scores = one_scores@meta.data[,c("orig.ident", "Morph_snn_res.1.5", y_var)]
colnames(df_one_scores) = c("Sample", "new_cluster", "scores")
df_one_scores$new_cluster = df_layer$Layer[match(df_one_scores$new_cluster, df_layer$Oringe_cluster)]
df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
df_one_scores_SM35 = df_one_scores


df_layer_27 = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
seu = tumor_SM27
one_scores = AddModuleScore(object = seu, features = list(gene), name = x)
y_var = paste(x, "1", sep="")
df_one_scores = one_scores@meta.data[,c("orig.ident", "Morph_snn_res.1.5", y_var)]
colnames(df_one_scores) = c("Sample", "new_cluster", "scores")
df_one_scores$new_cluster = df_layer_27$Layer[match(df_one_scores$new_cluster, df_layer_27$Oringe_cluster)]
df_one_scores$new_cluster = factor(df_one_scores$new_cluster, levels = levels)
df_one_scores_SM37 = df_one_scores

df_one_scores = rbind(df_one_scores_SM35, df_one_scores_SM37)

p = ggplot(df_one_scores, aes(x = new_cluster, y=scores, fill=Sample)) +
  geom_violin(trim=TRUE) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult = 1), geom="crossbar", width=0.2, position=position_dodge(0.9)) +
  stat_compare_means(
    aes(group = new_cluster), # 按照 new_cluster 进行分组
    comparisons = list(c("Layer1", "Layer2"), c("Layer1", "Layer3"), c("Layer2", "Layer3")), # 设置要比较的组
    method = "t.test", # 你可以选择其他检验方法，如 "wilcox.test" (Mann-Whitney)
    label = "p.value", # 显示p值而不是显著性星号
    hide.ns = TRUE # 隐藏非显著性结果
  ) +
  theme_classic() +
  labs(title="", x="", y ="emt set score") +
  theme(
    plot.margin = margin(t = 0.3, r = 0.3, b = 0.3, l = 0.3 ,unit = "cm"),
    axis.title = element_text(size = 22, colour = 'black'),
    axis.text.y = element_text(size = 18, colour = 'black'),
    axis.text.x = element_text(size = 20, colour = "black", angle = 45, hjust = 1),
    axis.ticks.length = unit(2, "mm")
  ) +
  scale_fill_manual(values = c("#E64B35", "#3C5488"))


return(p)           

}


#####Fig4g
load("/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")
Idents(object = ST_list) <- "area"
#TJH37
SpatialPlot(ST_list[[5]],group.by = 'area',pt.size.factor = 1.8,label = T)
SpatialPlot(ST_list[[5]],group.by = 'seurat_clusters',pt.size.factor = 1.8,label = T)

#主要癌灶位于Cluster11
TJH_37_sub<-subset(ST_list[[5]], seurat_clusters %in% '11')
SpatialPlot(TJH_37_sub,group.by = 'area',pt.size.factor = 1.8,label = T)

##BGI数据的每个spot位置信息在rownames里
TJH_37_sub$CB<-rownames(TJH_37_sub@meta.data)
#第五个字符开始是位置信息
TJH_37_sub@meta.data$CB <- substring(TJH_37_sub@meta.data$CB, 5)


#先确定癌灶的范围
SpatialDimPlot(TJH_37_sub,pt.size.factor = 1,interactive = T,image.alpha = 1)

#癌灶大致位置范围在16500-21363
TJH_37_sub<-subset(TJH_37_sub, CB > '16500' &  CB < '21363')
SpatialPlot(TJH_37_sub,group.by = 'area',pt.size.factor = 1.8,label = T)

#右侧有3个离散点，影响可视化，删除
SpatialDimPlot(TJH_37_sub,pt.size.factor = 1,interactive = T,image.alpha = 1)

#位置分别是'20442','20443','21006'
TJH_37_sub$label<-'1'
TJH_37_sub$label[TJH_37_sub$CB %in% c('20442','20443','21006')]<-'0'
TJH_37_sub<-subset(TJH_37_sub, label %in% c('1'))

Layer1<-c('SPRR3','SPRR2A','SPRR2E','CNFN','S100A7','SPRR1B','SPRR2D','ANXA1','KRT6B','RHCG','EMP1','CD24')
Layer2<-c('ZFP36','KRT16','FGFBP1','SFN','NDUFA4L2','NTS','SERPINB4','MT2A','MT1E','MT1G','KRT15' ,'MIR205HG')
Layer3<-c('TNC','COL18A1','THBS2','COL7A1','PIK3R1','MMP11','STMN1','DPYSL3','IFITM3','GPC3','COL16A1','PTN')
gene<-c("MT2A","NME2","IFITM3")
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(Layer1), name = "Layer1",nbin = 10)
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(Layer2), name = "Layer2",nbin = 10)
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(Layer3), name = "Layer3",nbin = 10)
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(gene), name = "EMT",nbin = 10)
SpatialPlot(TJH_37_sub,features = c('Layer11'),pt.size.factor = 9) ###Layer1,Layer2,Layer3


TJH_37_sub <- SCTransform(TJH_37_sub, assay = "Spatial", verbose = FALSE)
DefaultAssay(TJH_37_sub) <- "SCT"


library("cancersea")
data('available_pathways')

pathways <- list(angiogenesis=angiogenesis, apoptosis=apoptosis, cell_cycle=cell_cycle, differentiation=differentiation, dna_damage=dna_damage,
                 dna_repair=dna_repair, emt=emt, hypoxia=hypoxia, inflammation=inflammation, invasion=invasion, metastasis=metastasis, proliferation=proliferation,
                 quiescence=quiescence, stemness=stemness)
                 
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(pathways[["dna_repair"]][["symbol"]]), name = "dna_repair",nbin = 10)
TJH_37_sub<- AddModuleScore(object = TJH_37_sub, features = list(pathways[["hypoxia"]][["symbol"]]), name = "hypoxia",nbin = 10)
TJH_37_sub <- AddModuleScore(object = TJH_37_sub, features = list(pathways[["proliferation"]][["symbol"]]), name = "proliferation",nbin = 10)
TJH_37_sub<- AddModuleScore(object = TJH_37_sub, features = list(pathways[["cell_cycle"]][["symbol"]]), name = "cell_cycle" ,nbin = 10)
TJH_37_sub<- AddModuleScore(object =TJH_37_sub, features = list(pathways[["invasion"]][["symbol"]]), name = "invasion" ,nbin = 10)
SpatialPlot(TJH_37_sub,features = c('cell_cycle'),pt.size.factor = 9) ### cell_cycle,dna_repair,hypoxia,EMT

##### Fig4h
plot_correlation <- function(correlation_result, x, y, x_label, y_label) {
  ggplot(TJH90_sub@meta.data, aes(x = x, y = y)) +
    geom_point(alpha = 1, col = '#4682B4') +  
    stat_smooth(method = "lm", col = "red") +  
    labs(
      title = paste("R = ", round(correlation_result$estimate, 2),  "\n", 
                    paste("P = ", formatC(correlation_result$p.value, format = "e", digits = 2))),
      x = x_label,
      y = y_label
    ) + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.title = element_text(size = 18, face = "bold", color = "black"),
      axis.text = element_text(size = 15, face = "bold", color = "black"),
      plot.title = element_text(size = 20, face = "bold", color = "black")
    )
}

p = cor.test(TJH_37_sub@meta.data$Layer11, TJH90_sub@meta.data$EMT1)
plot_correlation(p, TJH_37_sub@meta.data$Layer11, TJH90_sub@meta.data$EMT1, "layer1", "EMT") #### cell_cycle,dna_repair,hypoxia,EMT

p = cor.test(TJH_37_sub@meta.data$Layer31, TJH90_sub@meta.data$EMT1)
plot_correlation(p, TJH_37_sub@meta.data$Layer31, TJH90_sub@meta.data$EMT1, "layer3", "EMT")#### cell_cycle,dna_repair,hypoxia,EMT
saveRDS(TJH_37_sub,file="TJH_37_sub.rds")