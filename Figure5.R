library(Seurat)
library(tidyverse)
library(patchwork)
library(Cottrazm)
library(RColorBrewer)
library(patchwork)
library(ggtree)
library(BiocGenerics)
library(readr)
library(rtracklayer)
library(infercnv)
library(phylogram)
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
#read SME data
SME_objects = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/Cottrazm/all_sample_SME_seurat.rds.gz")
extract_subset_focal_spot = function(seu, clusters=NULL){
    Idents(object = seu) <- "Morph_snn_res.1.5"
    focal_spot = subset(x = seu, idents = clusters)
    return(focal_spot)
}

SM35_4_colors = c("#2C77B4","#2C9F2D","#F77F0E","#D62728") # 5，6，9，13
SM27_3_colors = c("#2C77B4", "#F77F0E", "#D62728") # 0，2，9 "#E04D50", "#4374A5", "#F08A21"
layer_colors =c("#d9cdab","#398E7B", "#F4C5A9", "#B4EDF4")# normal, layer1, layer2, layer3

trajection_color = c("#FAFE1B", "#48E5FC", "#E43DEC", "#F1655C", "#1A8B4C", "#27749F", "#7E0438", "#F0D6DF")

##Figure 5a, 5b,5c,5d
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

####SM35
SM35 = SME_objects[[10]]
SM35 <- SCTransform(SM35, assay = "Spatial", verbose = FALSE)
tumor_SM35 = extract_subset_focal_spot(SM35,c(11, 12, 13, 5, 6, 9))
tumor_SM35 = CellCycleScoring(tumor_SM35, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
df_layer = data.frame(Oringe_cluster = c(11,12, 13, 5, 6, 9), Layer = c("Normal", "Normal", "Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
tumor_SM35@meta.data$layer = factor(tumor_SM35@meta.data$layer, levels=c("Normal", "Layer1", "Layer2", "Layer3"))

#构建matrix
sample_ann <-  tumor_SM35@meta.data
gene_ann <- data.frame(gene_short_name = rownames(tumor_SM35@assays$Spatial), row.names =  rownames(tumor_SM35@assays$Spatial))
pd <- new("AnnotatedDataFrame",data=sample_ann)
fd <- new("AnnotatedDataFrame",data=gene_ann)
ct=as(as.matrix(tumor_SM35@assays$Spatial@counts), "sparseMatrix")#单细胞counts矩阵
sc_cds <- newCellDataSet(
  as(ct, "sparseMatrix"), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)

cds_1 <- sc_cds
cds_1 <- estimateSizeFactors(cds_1)
cds_1 <- estimateDispersions(cds_1)
cds_1 <- detectGenes(cds_1, min_expr = 0.1)
disp_table <- dispersionTable(cds_1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_1 <- setOrderingFilter(cds_1, unsup_clustering_genes$gene_id)
#plot_pc_variance_explained(cds_1, return_all = F) 
cds_1<- reduceDimension(cds_1, max_components = 2, num_dim = 8, reduction_method = 'tSNE', verbose = T)
cds_1<- clusterCells(cds_1)
#plot_cell_clusters(cds_1, 1, 2, color = "seurat_clusters")

##differential genes
expressed_genes <- row.names(subset(fData(cds_1),  num_cells_expressed >= 5))
clustering_DEG_genes <- differentialGeneTest(cds_1[expressed_genes,],fullModelFormulaStr = '~layer',cores = 10)


##trajactory
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1500]
cds_1 <- setOrderingFilter(cds_1,ordering_genes = ordering_genes)
cds_1 <- reduceDimension(cds_1, method = 'DDRTree')
cds_1 <- orderCells(cds_1, reverse = FALSE)

p1 = plot_cell_trajectory(cds_1, cell_size=2,color_by = "Pseudotime") + scale_color_gradient2(midpoint=mean(pData(cds_1)$Pseudotime), low = "#413FA0", mid = "#FDFFBE", high = "#A91435")
#ggsave("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/final_figures/Figure_5A_1.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
p2 = plot_cell_trajectory(cds_1, cell_size=2,color_by = "layer") + scale_color_manual(values=layer_colors) + guides(color = guide_legend(nrow = 2))
#ggsave("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/final_figures/Figure_5A_2.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
p3 = plot_cell_trajectory(cds_1, cell_size=2,color_by = "Phase") + scale_color_manual(values=c("#FDE723", "#22958B", "#430154")) + guides(color = guide_legend(nrow = 2))
p = p1 + p2 + p3

#SM27 
SM27 = SME_objects[[5]]
SM27 <- SCTransform(SM27, assay = "Spatial", verbose = FALSE)
tumor_SM27 = extract_subset_focal_spot(SM27,c(2, 0, 9))
tumor_SM27 = CellCycleScoring(tumor_SM27, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
tumor_SM27@meta.data$layer = factor(tumor_SM27@meta.data$layer, levels=c("Layer1", "Layer2", "Layer3"))

#构建matrix
sample_ann <-  tumor_SM27@meta.data
gene_ann <- data.frame(gene_short_name = rownames(tumor_SM27@assays$Spatial), row.names =  rownames(tumor_SM27@assays$Spatial))
pd <- new("AnnotatedDataFrame",data=sample_ann)
fd <- new("AnnotatedDataFrame",data=gene_ann)
ct=as(as.matrix(tumor_SM27@assays$Spatial@counts), "sparseMatrix")#单细胞counts矩阵
sc_cds <- newCellDataSet(
  as(ct, "sparseMatrix"), 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)

cds_2 <- sc_cds
cds_2 <- estimateSizeFactors(cds_2)
cds_2 <- estimateDispersions(cds_2)
cds_2 <- detectGenes(cds_2, min_expr = 0.1)
disp_table <- dispersionTable(cds_2)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds_2 <- setOrderingFilter(cds_2, unsup_clustering_genes$gene_id)
#plot_pc_variance_explained(cds_1, return_all = F) 
cds_2<- reduceDimension(cds_2, max_components = 2, num_dim = 8, reduction_method = 'tSNE', verbose = T)
cds_2<- clusterCells(cds_2)
#plot_cell_clusters(cds_1, 1, 2, color = "seurat_clusters")

##differential genes
expressed_genes <- row.names(subset(fData(cds_2),  num_cells_expressed >= 5))
clustering_DEG_genes <- differentialGeneTest(cds_2[expressed_genes,],fullModelFormulaStr = '~layer',cores = 10)


##trajactory
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1500]
cds_2 <- setOrderingFilter(cds_2,ordering_genes = ordering_genes)
cds_2 <- reduceDimension(cds_2, method = 'DDRTree')
cds_2 <- orderCells(cds_2, reverse = TRUE)

p1 = plot_cell_trajectory(cds_2, cell_size=2,color_by = "Pseudotime") + scale_color_gradient2(midpoint=mean(pData(cds_1)$Pseudotime), low = "#413FA0", mid = "#FDFFBE", high = "#A91435")
#ggsave("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/final_figures/Figure_5A_1.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
p2 = plot_cell_trajectory(cds_2, cell_size=2,color_by = "layer") + scale_color_manual(values=layer_colors[2:4]) + guides(color = guide_legend(nrow = 2))
#ggsave("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/final_figures/Figure_5A_2.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
p3 = plot_cell_trajectory(cds_2, cell_size=2,color_by = "Phase") + scale_color_manual(values=c("#FDE723", "#22958B", "#430154")) + guides(color = guide_legend(nrow = 2))
p = p1 + p2 + p3

readr::write_rds(cds_1, "/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_trajectory.rds.gz",compress = "gz")
readr::write_rds(cds_2, "/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM27_trajectory.rds.gz",compress = "gz")

df_group_SM35 = tumor_SM35@meta.data[,c("layer", "Phase")] %>% 
group_by(layer, Phase) %>% 
summarise(n = n())

p1 = ggplot(df_group_SM35, aes(x=layer, y=n, fill=Phase)) + 
  geom_bar(stat="identity", position="fill", color="black", lwd=0.25) + 
  labs(title="Cell cycle state of SM35",  y="Frequence") +
  theme(axis.title.x = element_blank())+ 
  guides(fill=F)+
  scale_fill_manual(values=c("#FDE723", "#22958B", "#430154"))


df_group_SM27 = tumor_SM27@meta.data[,c("layer", "Phase")] %>% 
group_by(layer, Phase) %>% 
summarise(n = n())

p2 = ggplot(df_group_SM27, aes(x=layer, y=n, fill=Phase)) + 
  geom_bar(stat="identity", position="fill", color="black", lwd=0.25) +
  labs(title="Cell cycle state of SM27", y="Frequence") +
  theme(axis.title.x = element_blank())+ 
  scale_fill_manual(values=c("#FDE723", "#22958B", "#430154"))

p = p1 + p2 + plot_layout(ncol=1)


#### Figure5e,5f
SM35_infercnv <-
  STCNVScore(
    TumorST = tumor_SM35,
    assay = "Spatial",
    Sample = 'SM35',
    OutDir = '/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/CNV/'
    )
cell_gene_state = c("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/CNV/SM35//HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat", 
                    "/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/CNV/SM27/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat")

cell_groupings = c("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/CNV/SM35/17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings", 
                    "/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/CNV/SM27//17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings")

lapply(cell_gene_state, function(x){
       cell_gene_state <- read.table(x,header = T)
}) %>% do.call(rbind, .) -> cell_gene_state_res

lapply(cell_groupings, function(x){
       cell_groupings <- read.table(x,header = T)
}) %>% do.call(rbind, .) -> cell_groupings_res

cell_gene_state_res %>%  mutate(cell_group_name = str_remove(cell_group_name, '.*\\.')) -> cell_gene_state_res
cell_groupings_res %>% mutate(cell_group_name = str_remove(cell_group_name, '.*\\.')) -> cell_groupings_res

cell_gene_state_res %>% mutate(cell_group_name = str_remove(cell_group_name, '.*\\.')) %>% 
    group_by(cell_group_name) %>%
    mutate(cnv_score = sum(abs(state - 3))) %>%
    distinct(cell_group_name, cnv_score) -> score_df

cell_groupings_res %>%
    left_join(score_df, by = 'cell_group_name') %>% select(bc = cell, cnv_score_ref_N1_N5 = cnv_score) -> score_df

##spot information
SM27 = SME_objects[[5]]
SM27 <- SCTransform(SM27, assay = "Spatial", verbose = FALSE)

SM35 = SME_objects[[10]]
SM35 <- SCTransform(SM35, assay = "Spatial", verbose = FALSE)

tumor_SM27 = extract_subset_focal_spot(SM27,c(0,2,9))
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]

tumor_SM35 = extract_subset_focal_spot(SM35,c(13, 5, 6, 9))
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]

#SM35
tumor_SM35_1 = RenameCells(object = tumor_SM35, add.cell.id = "SM35")
tumor_SM35_1@meta.data$CNV_score = score_df$cnv_score_ref_N1_N5[match(rownames(tumor_SM35_1@meta.data), score_df$bc)]
tumor_SM35_1 = tumor_SM35_1@meta.data[,c("layer", "CNV_score")]

tumor_SM35_1$layer = factor(tumor_SM35_1$layer, levels = c("Layer1", "Layer2", "Layer3"))
p = ggplot(tumor_SM35_1, aes(x = layer, y=CNV_score, fill=layer)) +
    geom_violin(trim=FALSE, color="darkred")+
    #stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="pointrange") +
    stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="crossbar", width=0.2) +
    theme_classic() +
    labs(title="",x="", y = "SM35 CNV Score") +
    theme(
      plot.margin = margin(t = 0.3, r = 0.3,b = 0.3,l = 0.3 ,unit = "cm"),
      #plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
      plot.title = element_text(size=16, hjust = 1.0),
      axis.title = element_text(size=16, colour = 'black'),
      axis.text.y = element_text(size = 12,colour = 'black'),
      axis.text.x = element_text(size = 14,colour = "black",angle = 0,hjust = 0.5),
      axis.ticks.length = unit(2,"mm")) +
    scale_fill_manual(values = c( "#3A8E7B", "#F4C4AA", "#B4EDF4")) +
    NoLegend()
p

tumor_SM27_1 = RenameCells(object = tumor_SM27, add.cell.id = "SM27")
tumor_SM27_1@meta.data$CNV_score = score_df$cnv_score_ref_N1_N5[match(rownames(tumor_SM27_1@meta.data), score_df$bc)]
tumor_SM27_1 = tumor_SM27_1@meta.data[, c("layer", "CNV_score")]

tumor_SM27_1$layer = factor(tumor_SM27_1$layer, levels = c("Layer1", "Layer2", "Layer3"))
p = ggplot(tumor_SM27_1, aes(x = layer, y=CNV_score, fill=layer)) +
    geom_violin(trim=FALSE, color="darkred")+
    #stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="pointrange") +
    stat_summary(fun.data="mean_sdl",fun.args = list(mult = 1), geom="crossbar", width=0.2) +
    theme_classic() +
    labs(title="",x="", y = "SM27 CNV Score") +
    theme(
      plot.margin = margin(t = 0.3, r = 0.3,b = 0.3,l = 0.3 ,unit = "cm"),
      #plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5),
      plot.title = element_text(size=16, hjust = 1.0),
      axis.title = element_text(size=16, colour = 'black'),
      axis.text.y = element_text(size = 12,colour = 'black'),
      axis.text.x = element_text(size = 14,colour = "black",angle = 0,hjust = 0.5),
      axis.ticks.length = unit(2,"mm")) +
    scale_fill_manual(values = c( "#3A8E7B", "#F4C4AA", "#B4EDF4")) +
    NoLegend()
p

tumor_SM27 = SM27
tumor_SM35 = SM35

#SM35
tumor_SM35_1 = RenameCells(object = tumor_SM35, add.cell.id = "SM35")
tumor_SM35_1@meta.data$CNV_score = score_df$cnv_score_ref_N1_N5[match(rownames(tumor_SM35_1@meta.data), score_df$bc)]
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
p = SpatialPlot(tumor_SM35_1, features = c("CNV_score"), stroke = 0.5, min.cutoff = 2000, max.cutoff = 5000,) + 
    theme(legend.position="top") + 
    scale_fill_gradientn(
        breaks = seq(2000, 5000, 1000), 
        limits = c(2000, 5000),
        name = 'CNV score',
        colours = SpatialColors(n = 100))
p

#SM27
tumor_SM27_1 = RenameCells(object = tumor_SM27, add.cell.id = "SM27")
tumor_SM27_1@meta.data$CNV_score = score_df$cnv_score_ref_N1_N5[match(rownames(tumor_SM27_1@meta.data), score_df$bc)]
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
p = SpatialPlot(tumor_SM27_1, features = c("CNV_score"), stroke = 0.5, min.cutoff = 2000, max.cutoff = 5000,) + 
    theme(legend.position="top") + 
    scale_fill_gradientn(
        breaks = seq(2000, 5000, 1000), 
        limits = c(2000, 5000),
        name = 'CNV score',
        colours = SpatialColors(n = 100))
p
