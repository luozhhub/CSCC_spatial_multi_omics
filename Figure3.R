.libPaths(c('~/seurat_v4/', .libPaths()))
library(Seurat, lib.loc = "/media/desk16/tly0202/seurat_v4/")
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
library(Cottrazm)

load("/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")
load("/media/desk16/tly02/Project/CESC/merged_SM_8samples_merged_SM_1.rdata")
###Fig3a
##空间连续性
#提取肿瘤样本的表皮细胞
#获取某个spot周围的6个spot
get_circle_n_spot = function(coor_df, spot){
  #spot = as.vector(spot)
  spot_top_left = c(row=spot["row"] + 1, col=spot["col"] - 1)
  spot_top_right = c(row=spot["row"] + 1, col=spot["col"] + 1)
  spot_right = c(row=spot["row"], col=spot["col"] + 2)
  spot_bottom_right = c(row=spot["row"] - 1, col=spot["col"] + 1)
  spot_bottom_left = c(row=spot["row"] - 1, col=spot["col"] - 1)
  spot_left = c(row=spot["row"], col=spot["col"] - 2)
  
  df_circle = data.frame(row=c(spot_top_left[1], spot_top_right[1],spot_right[1],spot_bottom_right[1],spot_bottom_left[1],spot_left[1]),
                         col=c(spot_top_left[2], spot_top_right[2],spot_right[2],spot_bottom_right[2],spot_bottom_left[2],spot_left[2])
  )
  
  df_circle$cor_name = paste0(df_circle$row, "x", df_circle$col)
  df_circle = coor_df[coor_df$cor_name %in% df_circle$cor_name,]
  return(df_circle)
}

#获取一些spot周围的spot
get_circle_all_spot = function(coor_df, spot_all){
  res = purrr::map(rownames(spot_all), function(x){
    df_one = spot_all[x,]
    #pinrt(x)
    df_one = get_circle_n_spot(coor_df=coor_df, spot=c("row"=spot_all[x, "row"], "col"=spot_all[x, "col"]))
    df_one = df_one %>% tibble::rownames_to_column(., "barcode")
    return(df_one)
  }) %>% bind_rows() %>% distinct()
  rownames(res) = res$barcode
  res=res[,!colnames(res) %in% c("barcode")]
  return(res)
}

####依次计算 "SM35", "SM27","SM25", "SM32", "SM34", "SM21"的连续性
sample_id = "SM21" #### 此处依次为 "SM35", "SM27","SM25", "SM32", "SM34", "SM21"
#提取坐标,所有spot
Idents(object = SM) <- "orig.ident"
whole_sample = subset(x = SM, idents = sample_id)
coor_w  = whole_sample@images[[sample_id]]@coordinates
# 转换坐标
xDiml <- max(coor_w$row)
yDiml <- max(coor_w$col)
coor_w$row = xDiml - coor_w$row
coor_w$cor_name = paste0(coor_w$row, "x", coor_w$col)


##肿瘤的spot
Idents(object = SM) <- "cellType"
focal_spot = subset(x =SM, idents = "Epi")
Idents(object = focal_spot) <- "orig.ident"
one_sample_spot = subset(x = focal_spot, idents = sample_id)

Idents(object = SM) <- "cellType"

#提取坐标,肿瘤spot
coor_1  = one_sample_spot@images[[sample_id]]@coordinates
# 转换坐标
xDiml <- max(coor_w$row)
yDiml <- max(coor_w$col)
coor_1$row = xDiml - coor_1$row
coor_1$cor_name = paste0(coor_1$row, "x", coor_1$col)


result_all_spots = get_circle_all_spot(coor_df = coor_w, spot_all = coor_1)
percent = nrow(coor_1)/nrow(result_all_spots)
cat("sample:", sample_id, "tumor:", nrow(coor_1), "all:", nrow(result_all_spots), "percent:", percent, "\n")

##sample: SM35 tumor: 1817 all: 2608 percent: 0.6967025
##sample: SM27 tumor: 1167 all: 1727 percent: 0.6757383
##sample: SM25 tumor: 1751 all: 3051 percent: 0.5739102
##sample: SM32 tumor: 1101 all: 1685 percent: 0.6534125
##sample: SM34 tumor: 1368 all: 2826 percent: 0.4840764
##sample: SM21 tumor: 1718 all: 3187 percent: 0.539065 

##纯度
DefaultAssay(focal_spot) = "RCTD"
Idents(object = focal_spot) <- "orig.ident"
subset_tumor_data <- subset(focal_spot, subset = orig.ident %in% c("SM35", "SM27","SM25", "SM32", "SM34", "SM21"))

data <- FetchData(subset_tumor_data, vars = c("EpithelialCells", "ident"))  # 替换 "ident" 为你的分组列名
median_values <- data %>%
  group_by(ident) %>%
  summarise(median_value = median(EpithelialCells))

###将上面计算每个样本的连续度值放在此处，注意样本名称排列顺序一致
median_values$continue = c(0.68,0.65,0.54,0.57,0.48,0.70)


# 绘制散点图
p = ggplot(median_values, aes(x = median_value, y = continue, color = ident)) +
  geom_point(size = 4) +  # 绘制点，设置点的大小
  geom_text(aes(label = ident), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +  # 添加标签
  labs(
    x = "Median Value",  # x 轴标签
    y = "Continue",      # y 轴标签
    title = "Scatter Plot of Median Value vs Continue"  # 图表标题
  ) +
  theme_minimal() +  # 使用简洁的主题
  theme(
    legend.position = "none",  # 隐藏图例
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(colour = "black", fill=NA, size=1)  # 添加边框
  ) +
  ylim(0.4, 0.75) +
  xlim(0.8,1)
p

###Fig3b1
DefaultAssay(SM) = "RCTD"
p = SpatialFeaturePlot(object = SM, features = "EpithelialCells", alpha = c(0.1, 1), images ="SM25")###"SM35", "SM27","SM25", "SM32", "SM34", "SM21"
p

# Basic density Fig3b2
df_SM27 = data[data$ident %in% c("SM27", "SM35", "SM21", "SM34"), ]
df_SM27$ident = factor(df_SM27$ident, levels = c("SM27", "SM35", "SM21", "SM34"))
p <- ggplot(df_SM27, aes(x=EpithelialCells)) + 
  geom_density(fill="lightblue") +
  theme_minimal() +  # 使用简洁的主题
  theme(
    legend.position = "none",  # 隐藏图例
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(colour = "black", fill=NA, size=1)  # 添加边框
  ) +
  xlab("Epi percentage")+
  geom_vline(
    aes(xintercept = 0.9),  # 关键修改：使用xintercept而非x
    color = "blue", 
    linetype = "dashed", 
    linewidth = 1  # 注意：ggplot2 3.4.0+推荐用linewidth替代size
  ) +
  facet_wrap(~ ident, ncol = 4)
p

#Figure 3e
#分组
df_group = data.frame(Oringe_sample = c("SM27", "SM35", "SM21", "SM34", "SM25", "SM32", "N1", "N5"), group = c("Continue", "Continue", "Infiltration", "Infiltration", "Intermediary", "Intermediary", "Normal", "Normal"))
focal_spot@meta.data$group = df_group$group[match(focal_spot@meta.data$orig.ident, df_group$Oringe_sample)]

DefaultAssay(focal_spot) <- "SCT"
Idents(focal_spot) = "group"
focal_spot <- PrepSCTFindMarkers(object = focal_spot)

VolcanoPlot=function(dif, log2FC=0.25, padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=1.2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*P.adj) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=16)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print(head(dd_text, n=3))
  print(tail(dd_text, n=3))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5,   #基因名字体大小
                       colour="black",alpha=1)
}


deg_genes = FindMarkers(focal_spot, ident.1 = "Continue", group.by = 'group', ident.2 = "Infiltration")

head(deg_genes)

dif1=data.frame(
  symbol=rownames(deg_genes),
  log2FoldChange=deg_genes$avg_log2FC,
  padj=deg_genes$p_val_adj)

#可视化


p = VolcanoPlot(dif1, padj=0.05, title="SM27/35 VS SM21/SM34",label.max = 10)
p

####Fig3c
##re annotation
cell_type = c("Mixed1", "Fibro", "B", "Epi", "Epi", "Epi", "Epi", "Immune", "Epi", "Endo", "Epi", "Fibro")
df_celltype = data.frame(Oringe_cluster = 0:11, type = cell_type)
ST_list[[1]]@meta.data$cellType = df_celltype$type[match(ST_list[[1]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("Epi", "Fibro", "Epi", "B", "Immune", "Epi", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:6, type = cell_type)
ST_list[[2]]@meta.data$cellType = df_celltype$type[match(ST_list[[2]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("Fibro", "Epi", "Endo", "B", "Epi", "Mixed", "Epi", "Fibrom", "Endo", "Mast")
df_celltype = data.frame(Oringe_cluster = 0:9, type = cell_type)
ST_list[[3]]@meta.data$cellType = df_celltype$type[match(ST_list[[3]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("Fibro", "Epi", "Immune", "Epi", "Epi", "Epi", "Mixed", "Fibro", "Fibro", "Fibro", "Epi", "Epi", "Endo")
df_celltype = data.frame(Oringe_cluster = 0:12, type = cell_type)
ST_list[[4]]@meta.data$cellType = df_celltype$type[match(ST_list[[4]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Mixed1", "Fibro", "B", "B", "Epi", "Epi", "Epi", "Epi", "Epi", "Mixed2", "Immune", "Epi", "Epi", "Immune", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:14, type = cell_type)
ST_list[[5]]@meta.data$cellType = df_celltype$type[match(ST_list[[5]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Immune", "Fibro", "Epi", "Epi", "Epi", "B", "T", "Fibro", "B", "Endo", "Myelo", "Mixed1", "Epi", "B", "Mast", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:15, type = cell_type)
ST_list[[6]]@meta.data$cellType = df_celltype$type[match(ST_list[[6]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Mixed1", "B", "Epi", "Fibro", "Mixed2", "Epi", "T", "Endo", "Mixed3", "Epi", "Myelo", "Mixed4", "Mast", "Fibro")
df_celltype = data.frame(Oringe_cluster = 0:13, type = cell_type)
ST_list[[7]]@meta.data$cellType = df_celltype$type[match(ST_list[[7]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Fibro", "B", "Mixed1", "Fibro", "Epi", "Mixed3", "Endo", "Epi", "Fibro", "Mixed2", "Fibro", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:11, type = cell_type)
ST_list[[8]]@meta.data$cellType = df_celltype$type[match(ST_list[[8]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Fibro", "Mixed1", "Epi", "Epi", "Epi", "Fibro", "Mixed2", "Epi", "Mixed3", "Epi", "B", "Epi", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:12, type = cell_type)
ST_list[[9]]@meta.data$cellType = df_celltype$type[match(ST_list[[9]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Fibro", "Fibro", "Epi", "Epi", "Epi", "B", "Epi", "T", "Epi", "B", "Epi", "Mixed1", "Epi", "Mixed2", "Epi", "Epi", "Immune", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:17, type = cell_type)
ST_list[[10]]@meta.data$cellType = df_celltype$type[match(ST_list[[10]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Fibro", "Fibro", "Epi", "Epi", "Epi", "Epi", "Mixed1", "B", "Fibro", "Epi", "Mixed2", "Goblet", "Epi", "Fibro", "Myeloi")
df_celltype = data.frame(Oringe_cluster = 0:14, type = cell_type)
ST_list[[11]]@meta.data$cellType = df_celltype$type[match(ST_list[[11]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("T", "Fibro", "Epi", "Fibro", "Fibro", "Epi", "Epi", "Fibro", "Mixed", "Fibro")
df_celltype = data.frame(Oringe_cluster = 0:9, type = cell_type)
ST_list[[12]]@meta.data$cellType = df_celltype$type[match(ST_list[[12]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]


cell_type = c("Fibro", "Immune", "Epi", "B", "Epi", "Goblet", "Fibro", "Immune", "Epi", "Fibro", "Endo", "Mixed")
df_celltype = data.frame(Oringe_cluster = 0:11, type = cell_type)
ST_list[[13]]@meta.data$cellType = df_celltype$type[match(ST_list[[13]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("Mixed1", "Mixed2", "Epi", "B", "Epi", "Fibro", "Fibro", "T", "Epi")
df_celltype = data.frame(Oringe_cluster = 0:8, type = cell_type)
ST_list[[14]]@meta.data$cellType = df_celltype$type[match(ST_list[[14]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

cell_type = c("Mixed1", "Mixed2", "Mixed3", "Fibro","Fibro", "Epi", "B", "B", "Fibro", "Fibro", "miexd4", "Fibro")
df_celltype = data.frame(Oringe_cluster = 0:11, type = cell_type)
ST_list[[15]]@meta.data$cellType = df_celltype$type[match(ST_list[[15]]@meta.data$SCT_snn_res.0.5, df_celltype$Oringe_cluster)]

save(ST_list, file = "/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")

##空间连接度
#提取肿瘤样本的表皮细胞
#获取某个spot周围的8个spot
get_circle_n_spot = function(coor_df, spot){
  #spot = as.vector(spot)
  spot_top_left = c(row=spot["row"] + 1, col=spot["col"] - 1)
  spot_top_middle= c(row=spot["row"] + 1, col=spot["col"])
  spot_top_right = c(row=spot["row"] + 1, col=spot["col"] + 1)
  spot_right = c(row=spot["row"], col=spot["col"] + 1)
  spot_bottom_right = c(row=spot["row"] - 1, col=spot["col"] + 1)
  spot_bottom_middle= c(row=spot["row"]-1, col=spot["col"])
  spot_bottom_left = c(row=spot["row"] - 1, col=spot["col"] - 1)
  spot_left = c(row=spot["row"], col=spot["col"] - 1)
  
  df_circle = data.frame(row=c(spot_top_left[1],spot_top_middle[1], spot_top_right[1],spot_right[1],spot_bottom_right[1],spot_bottom_middle[1],spot_bottom_left[1],spot_left[1]),
                         col=c(spot_top_left[2],spot_top_middle[2], spot_top_right[2],spot_right[2],spot_bottom_right[2],spot_bottom_middle[2],spot_bottom_left[2],spot_left[2])
  )
  
  df_circle$cor_name = paste0(df_circle$row, "x", df_circle$col)
  df_circle = coor_df[coor_df$cor_name %in% df_circle$cor_name,]
  return(df_circle)
}

#获取一些spot周围的spot
get_circle_all_spot = function(coor_df, spot_all){
  res = purrr::map(rownames(spot_all), function(x){
    df_one = spot_all[x,]
    #pinrt(x)
    df_one = get_circle_n_spot(coor_df=coor_df, spot=c("row"=spot_all[x, "row"], "col"=spot_all[x, "col"]))
    df_one = df_one %>% tibble::rownames_to_column(., "barcode")
    return(df_one)
  }) %>% bind_rows() %>% distinct()
  rownames(res) = res$barcode
  res=res[,!colnames(res) %in% c("barcode")]
  return(res)
}


##提取全部癌灶--
ST_list2<-list()
for (i in 1:15) {
  ST_list2[[i]] <- subset(ST_list[[i]], subset = cellType %in% c("Epi"))
}

#定义样本名称，ST15个样本
sample <- c('TJH08','TJH24','TJH26','TJH34','TJH37','TJH45','TJH46','TJH48','TJH66','TJH76','TJH90',"TJH14", "TJH17","TJH35","TJH55")

##修改样本信息，依次修改 orig.ident
for (i in 1:length(sample)) {
  ST_list2[[i]]$orig.ident <- sample[i]
}


##合并全部样本
ST_tumor<-merge(x=ST_list2[[1]],y=ST_list2[ -1 ])
table(ST_tumor$orig.ident)
SpatialPlot(ST_tumor,ncol = 4,group.by = 'cellType')


#提取坐标,肿瘤spot
###初始化一个空的数据框---
percent_results <- data.frame(sample = character(), percent = numeric(), stringsAsFactors = FALSE)

#计算空间连接度
# 遍历所有样本
for (i in 1:length(sample)) {
  # 提取单个样本的全部spot
  coor_w  = ST_list[[i]]@images[[sample[i]]]@coordinates
  # 转换坐标
  xDiml <- max(coor_w$row)
  yDiml <- max(coor_w$col)
  coor_w$row = xDiml - coor_w$row
  coor_w$cor_name = paste0(coor_w$row, "x", coor_w$col)
  
  # 提取单个样本的癌灶spot
  coor_1  = ST_tumor@images[[sample[i]]]@coordinates
  # 转换坐标
  coor_1$row = xDiml - coor_1$row
  coor_1$cor_name = paste0(coor_1$row, "x", coor_1$col)
  
  # 获取肿瘤spot的所有信息
  result_all_spots <- get_circle_all_spot(coor_df = coor_w, spot_all = coor_1)
  
  # 计算当前样本的百分比
  percent <- nrow(coor_1) / nrow(result_all_spots)
  
  # 将结果存储到数据框
  percent_results <- rbind(percent_results, data.frame(sample = sample[i], percent = percent))
}

table(percent_results$percent)


##切换至RCTD
DefaultAssay(ST_tumor)<-'RCTD'
Idents(ST_tumor)<-ST_tumor$orig.ident

##变量一为RCTD中上皮细胞成分
data <- FetchData(ST_tumor, vars = c("EpithelialCells", "ident"))  # 替换 "ident" 为你的分组列名
median_values <- data %>%
  group_by(ident) %>%
  summarise(median_value = median(EpithelialCells))

##continue应该是percent
median_values$continue = percent_results$percent



ggplot(median_values, aes(x = median_value, y = continue, color = ident)) +
  geom_point(size = 4) +  # 绘制点，设置点的大小
  geom_text(aes(label = ident), vjust = -0.5, hjust = 0.5, size = 4, color = "black") +  # 添加标签
  labs(
    x = "Median Value",  # x 轴标签
    y = "Continue",      # y 轴标签
    title = "Scatter Plot of Median Value vs Continue"  # 图表标题
  ) +
  theme_minimal() +  # 使用简洁的主题
  theme(
    legend.position = "none",  # 隐藏图例
    panel.grid = element_blank(),  # 去掉背景网格
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 给图片加上黑色边框
  )+
  ylim(0.45, 0.9) +
  xlim(0.65,0.95)




##figure 3d1
DefaultAssay(ST_list[[5]]) = "RCTD"
SpatialPlot(ST_list[[5]], features = "EpithelialCells",alpha = c(0.1, 1))

DefaultAssay(ST_list[[11]]) = "RCTD"
SpatialPlot(ST_list[[11]], features = "EpithelialCells",alpha = c(0.1, 1))



DefaultAssay(ST_list[[13]]) = "RCTD"
SpatialPlot(ST_list[[13]], features = "EpithelialCells",alpha = c(0.1, 1))


DefaultAssay(ST_list[[3]]) = "RCTD"
SpatialPlot(ST_list[[3]], features = "EpithelialCells",alpha = c(0.1, 1))


##figure 3d2
df_TJH37 = data[data$ident %in% c("TJH37", "TJH90", "TJH17", "TJH26"), ]
df_TJH37$ident = factor(df_TJH37$ident, levels = c("TJH37", "TJH90", "TJH17", "TJH26"))
p <- ggplot(df_TJH37, aes(x=EpithelialCells)) + 
  geom_density(fill="lightblue") +
  theme_minimal() +  # 使用简洁的主题
  theme(
    legend.position = "none",  # 隐藏图例
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(colour = "black", fill=NA, size=1)  # 添加边框
  ) +
  xlab("Epi percentage")+
  geom_vline(
    aes(xintercept = 0.9),  # 关键修改：使用xintercept而非x
    color = "blue", 
    linetype = "dashed", 
    linewidth = 1  # 注意：ggplot2 3.4.0+推荐用linewidth替代size
  ) +
  facet_wrap(~ ident, ncol = 4)
p


####Figure 3f
VolcanoPlot=function(dif, log2FC=0.25, padj=0.05, 
                     label.symbols=NULL, label.max=30,
                     cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=1.2) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*P.adj) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=16)
    ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print(head(dd_text, n=3))
  print(tail(dd_text, n=3))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                       #size=2.5,   #基因名字体大小
                       colour="black",alpha=1)
}

##取4个样本TJH37，90，17，26
sub_list3 <- list(ST_list2[[5]], ST_list2[[11]], ST_list2[[13]], ST_list2[[3]])
#定义样本信息
sub_list3[[1]]$orig.ident='TJH37'
sub_list3[[2]]$orig.ident='TJH90'
sub_list3[[3]]$orig.ident='TJH17'
sub_list3[[4]]$orig.ident='TJH26'


#合并样本
sub_merge<-merge(sub_list3[[1]], y = sub_list3[2:4])

##标准化
sub_merge<- SCTransform(sub_merge,assay = 'Spatial') %>% RunPCA(npcs = 30, verbose = T)

####分组
df_group = data.frame(Oringe_sample = c('TJH37','TJH90', 'TJH17','TJH26'), group = c("Continue", "Continue", "Infiltration", "Infiltration"))
sub_merge@meta.data$group = df_group$group[match(sub_merge@meta.data$orig.ident, df_group$Oringe_sample)]

Idents(sub_merge) = "group"
#进行差异分析并可视化
DefaultAssay(sub_merge) <- "SCT"
Idents(sub_merge)<-sub_merge$orig.ident
sub_merge <- PrepSCTFindMarkers(object = sub_merge)
sub_deg<-FindMarkers(sub_merge ,ident.1 = c('Continue'),ident.2 =  c('Infiltration'),min.pct = 0.25,logfc.threshold = 0.25)
dif2=data.frame(
  symbol=rownames(sub_deg),
  log2FoldChange=sub_deg$avg_log2FC,
  padj=sub_deg$p_val_adj)

p = VolcanoPlot(dif2, padj=0.05, title="TJH37/90 VS TJH17/26",label.max = 10)
p


##### Fig3g
dif1_gene<-dif1[dif1$padj< 0.05 & dif1$log2FoldChange > 0.25,]
dif1_gene<-as.data.frame(dif1_gene$symbol)
dif2_gene<-dif2[dif2$padj< 0.05 & dif2$log2FoldChange > 0.25,]
dif2_gene<-as.data.frame(dif2_gene$symbol)
sig_gene<-merge(dif1_gene,dif2_gene,by="symbol")


##富集分析----
library(clusterProfiler)
#SYMBOL转ENTREZID
gid <- bitr(unique(sig_gene$symbol), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(sig_gene, gid, by=c('symbol' = 'SYMBOL'))

#GO通路富集分析
GO = compareCluster(ENTREZID ~ group, data = markers, fun='enrichGO', OrgDb='org.Hs.eg.db',ont='BP')
p<-dotplot(GO, label_format=50,showCategory=25) + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_color_gradient(high="#4b5cc4",low="#FE8D3C")+ggtitle('GO-BP')
p

list=p[["plot_env"]][["df"]]$Description

list2=list[c(1:4,6:11,41:66)]
dotplot(GO, label_format=50,showCategory=list2) + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_color_gradient(high="#4b5cc4",low="#FE8D3C")+ggtitle('GO-BP')

saveRDS(ST_list,file="TJH.rds")
saveRDS(ST_list2,file="TJH-tumor.rds")
saveRDS(ST_tumor,file="ST_tumor.rds")




#### Fig3h
gene<-sig_gene$gene

SM<-AddModuleScore(object =SM, features = list(gene), name = "genescore",nbin = 20)
SM21 <- subset(SM, subset = orig.ident %in% c("SM21"))###"SM35", "SM27","SM25", "SM32", "SM34", "SM21"
epi_cell = subset(SM21, subset = cellType %in% c("Epi"))
plot_data <- data.frame(
  EpithelialCells = epi_cell@assays$RCTD@data["EpithelialCells", ],
  genescore1 = epi_cell@meta.data$genescore1
)

plot_data <- plot_data[order(plot_data$EpithelialCells, decreasing = TRUE), ]
library(segmented)
lm_model <- lm(genescore1 ~ EpithelialCells, data = plot_data)
seg_model <- segmented(lm_model, seg.Z = ~EpithelialCells, psi = median(plot_data$EpithelialCells))
ggplot(plot_data, aes(x = EpithelialCells, y = genescore1)) +
  geom_point(color = '#4682B4') +  # 设置散点颜色为 '#4682B4'
  geom_line(aes(y = predict(seg_model)), color = "red", linewidth = 1.5) +
  geom_vline(xintercept = seg_model$psi[2], linetype = "dashed") +
  labs(x = "Epithelial percentage", y = "Gene set score") +  # 设置坐标轴标题
  theme_bw() +  # 背景为白色
  theme(
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 14),  # 设置横坐标标题字体大小
    axis.title.y = element_text(size = 14),  # 设置纵坐标标题字体大小
    axis.text.x = element_text(size = 14),   # 设置横坐标刻度字体大小
    axis.text.y = element_text(size = 14)    # 设置纵坐标刻度字体大小
  )

for (i in 1:length(ST_list)) {
  ST_list[[i]] <- AddModuleScore(object = ST_list[[i]], features = list(gene), name = "genescore",nbin = 20)
}

epi_cell = subset(ST_list[[15]], subset = cellType %in% c("Epi"))
plot_data <- data.frame(
  EpithelialCells = epi_cell@assays$RCTD@data["EpithelialCells", ],
  genescore1 = epi_cell@meta.data$genescore1
)
plot_data <- plot_data[order(plot_data$EpithelialCells, decreasing = TRUE), ]
library(segmented)
lm_model <- lm(genescore1 ~ EpithelialCells, data = plot_data)
seg_model <- segmented(lm_model, seg.Z = ~EpithelialCells, psi = median(plot_data$EpithelialCells))
ggplot(plot_data, aes(x = EpithelialCells, y = genescore1)) +
  geom_point(color = '#4682B4') +  # 设置散点颜色为 '#4682B4'
  geom_line(aes(y = predict(seg_model)), color = "red", linewidth = 1.5) +
  geom_vline(xintercept = seg_model$psi[2], linetype = "dashed") +
  labs(x = "Epithelial percentage", y = "Gene set score") +  # 设置坐标轴标题
  theme_bw() +  # 背景为白色
  theme(
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 14),  # 设置横坐标标题字体大小
    axis.title.y = element_text(size = 14),  # 设置纵坐标标题字体大小
    axis.text.x = element_text(size = 14),   # 设置横坐标刻度字体大小
    axis.text.y = element_text(size = 14)    # 设置纵坐标刻度字体大小
  )