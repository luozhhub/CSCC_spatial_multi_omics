####FigS1
library(Seurat)
packageVersion("Seurat")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(tidydr)
library(RColorBrewer)
#Marker gene FigS1a-1h
#NUPR1,FABP1,COL1A1,VWF,CD68,CD79A,IGHG1,CD3E,KRT8,TOP2A


merged_SM_1<-readRDS("merged_SM_8samples_merged_SM_1.rds")
limit=c(-2, 3)
breaks = c(-2,0.5, 3)
colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100)
plot_gene = function (cluster){
  p1 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM21',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.4,images = 'SM25',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM27',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM32',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p5 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM34',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p6 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM35',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p7 = SpatialPlot(merged_SM_1, features  = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'N1',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p8 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'N5',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol  = 8,nrow = 1,common.legend = T,legend = 'right')
}

limit=c(-1, 3)
breaks = c(-1,1, 3)
p1<-plot_gene(cluster='KRT8')
p1


####FigS1i
limit=c(0, 2)
breaks = c(-2,2, 6)
p2<-plot_gene(cluster='TOP2A')
p2


###### FigS2
####### FigS2a 
col = ("#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00")
umap<-DimPlot(merged_SM_1,cols = col,group.by = 'orig.ident')
umap

####### FigS2b UMAP
library(Seurat)
library(data.table)
library(tidyverse)
BD_objects = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/BD_20samples/ColonData_Harmony.rds.gz")
umap<-DimPlot(BD_objects,group.by = 'AllTypes')
umap

####### FigS2c cellmarker dotplot
Idents(object = BD_objects <- "AllTypes"
list<-c('EPCAM','CDH1','KRT5','TP63', 'ITGAX','CSFIR','CD14','FCGR3A','CLDN5','VWF','CDH5','KDR','COl1A1','COL1A2','LUM','MZB1','CD79A','MS4A1','CD3D','CD3E','CD2','KIT','IL1RL1','MS4A2')
dot <- DotPlot(BD_objects, features = list,col.max = 1.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 10),                        
    axis.title.x = element_text(size = 12),                        
    axis.title.y = element_text(size = 12),                        
    legend.title = element_text(size = 10),                        
    legend.text = element_text(size = 10)
    )+ # 设置图例文本的字体大小
  scale_color_gradientn(colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100), breaks = c(1, 0, -1),                                          
                        labels = c("2", "0", "-1")))
dot

#####FigS2d
cell_prop = as.data.frame(t(merged_SM_1@assays[["RCTD"]]@counts))

# 将比例作为新的 meta.data 列
merged_SM_1<- AddMetaData(merged_SM_1, metadata =cell_prop)
DimPlot(merged_SM_1,cols = col,group.by = 'Bcells')
FeaturePlot(merged_SM_1,
            features = "BCells",
            reduction = "umap",
            pt.size = 1.5)
FeaturePlot(merged_SM_1,
            features = "EPCAM",
            reduction = "umap",
            pt.size = 1.5)            
            
###### FigS2e,S2f
library(spacexr)
library(Matrix)
library("tidyverse")
library("patchwork")
library("SingleR")
library("Seurat")
library("ggplot2")
library("ggsci")
library("msigdbr")
library('RColorBrewer')
library(tidyverse)
library(Seurat)
library(mistyR)
source("./misty_utilities.R")

sample_objects = readr::read_rds("sample_object.rds")
##calculating hallmarks score matrix
outdir = file.path(getwd(), "cell_type_data")
library("msigdbr")
hallmarks_score = function (seu){
    h_gene_sets = msigdbr(species = "human", category = "H")
    all_hallmarks = unique(h_gene_sets$gs_name)
    scores = purrr::map(all_hallmarks, function(x){
    one_features = list(h_gene_sets$human_gene_symbol[h_gene_sets$gs_name == x])
    seu = AddModuleScore(object = seu, features = one_features, name = x)
    return(seu[[paste(x, "1", sep="")]])
    }) %>% bind_cols()

    colnames(scores) = str_sub(colnames(scores), start = 10, end = -2)
    return(scores)
}

#save hallmarks matrix
for (i in 1:length(sample_objects)){
    one_sample = sample_objects[[i]]
    pathway_score = hallmarks_score(one_sample)
    sample_name = names(one_sample@images)
    print(sample_name)
    readr::write_rds(pathway_score, file=file.path(getwd(),"cell_type_data","hall_marks", paste(sample_name, "_hallmarks_score.rds.gz",sep="")), compress = "gz")
}

##add pathway score assay
add_pathway_assay <- function(visium_slide){
    #visium_slide
    #read pathway scores
    sample_name = names(visium_slide@images)
    print(sample_name)
    pathways = readr::read_rds(file.path(getwd(),"cell_type_data","hall_marks", paste(sample_name, "_hallmarks_score.rds.gz",sep="")))
    pathways = t(pathways)
    adt_assay = CreateAssayObject(counts = pathways[,colnames(visium_slide)])
    visium_slide[["HallMarks"]] = adt_assay  
    return(visium_slide)
  }
pathway_obj = purrr::map(cell_type_obj, add_pathway_assay)
outdir = file.path(getwd(), "cell_type_data")
readr::write_rds(pathway_obj, file=file.path(outdir, "pathway_cellType_linai_seu.rds.gz"),compress = "gz")


run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/tissue_structure/misty/cell_map/cm_") {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 5,
                      "para" = 15)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}

for (i in 1:8){
assay_label <- "RCTD"
slide = pathway_obj[[i]]
out_dir_1 = file.path(getwd(), "cell_type_data","misty_")
sample_name = names(slide@images)
    
assay <- assay_label
DefaultAssay(slide) <- assay
useful_features <- rownames(slide)
useful_features <- useful_features[! useful_features %in% "prolif"]

mout <- run_colocalization(slide = slide,
        useful_features = useful_features,
        out_label = sample_name,
        assay = assay,
        misty_out_alias = out_dir_1)
}

misty_outs <- list.dirs(file.path("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/", "cell_type_data"),full.names = T,recursive = F)
misty_outs
misty_outs <- misty_outs[grepl("RCTD", misty_outs)]
misty_res <- collect_results(misty_outs)
mistyR::plot_interaction_heatmap(misty_res, "intra", cutoff = 0)
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)


#####Fig S3
DefaultAssay(merged_SM_1)<-'RCTD'
SpatialPlot(merged_SM_1,features = 'EpithelialCells') #### BCells, EndothelialCells, EpithelialCells, FibroblastsCells, MastCells, MyeloidCells, TCells 


#####Fig S4
DefaultAssay(merged_SM_1)<-'RCTD'
SpatialPlot(merged_SM_1,features = 'EpithelialCells')

col = c('1'="#7fc97f", '2'="#beaed4", '3'="#fdc086", '4'="#386cb0", '5'="#e0abcf", '6'="#a34e3b", '7'="#efbd25", '8'="#1b9e77", '9'="#d95f02", '10'="#6c67a5", '11'="#d01b2a", '12'="#43acde", '13'="#666666")

#order  
My_levels <- c('1','4','9','2','3','8','12','13','11','5','6','7','10')
Idents(merged_SM_1) <- merged_SM_1@meta.data$integrated_snn_res.0.5
Idents(merged_SM_1) <- factor(Idents(merged_SM_1), levels= My_levels)

#subset sample
SM21<-subset(merged_SM_1,orig.ident %in% 'SM21')
SM25<-subset(merged_SM_1,orig.ident %in% 'SM25')
SM27<-subset(merged_SM_1,orig.ident %in% 'SM27')
SM32<-subset(merged_SM_1,orig.ident %in% 'SM32')
SM34<-subset(merged_SM_1,orig.ident %in% 'SM34')
SM35<-subset(merged_SM_1,orig.ident %in% 'SM35')

merged_SM_1@assays[["RCTD"]]@counts@Dimnames[[1]]


samples <- c("SM21", "SM25", "SM27", "SM32","SM34","SM35")
cell_types <- c("BCells", "TCells", "MastCells", "MyeloidCells", "FibroblastsCells","EndothelialCells", "EpithelialCells")


# Vlnplot

plot_cell_type <- function(cell_type) {
  plot_list <- list()
  
  for (sample in samples) {
    p <- VlnPlot(get(sample), features = cell_type, pt.size = 0, cols = col,y.max =1) +
      NoLegend() + xlab(NULL) + ylab('Cell Percentage') + ggtitle(paste(sample, cell_type))
    plot_list[[sample]] <- p
  }
  
  p <- ggarrange(plotlist = plot_list, ncol = 6, nrow = 1)
  
  
}
p

##### FigS6
###### FigS6a,S6c
load("/media/desk16/tly02/Project/CESC/merged_SM_8samples_merged_SM_1.rdata")
limit=c(-2, 3)
breaks = c(-2,0.5, 3)
colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100)
plot_gene = function (cluster){
  p1 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM21',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p2 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.4,images = 'SM25',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p3 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM27',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p4 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM32',image.alpha = 0)+ scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p5 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.5,images = 'SM34',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p6 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'SM35',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p7 = SpatialPlot(merged_SM_1, features  = cluster,stroke = 0.15, pt.size.factor =1.4,images = 'N1',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  p8 = SpatialPlot(merged_SM_1, features = cluster,stroke = 0.15,pt.size.factor = 1.6,images = 'N5',image.alpha = 0)+scale_fill_gradientn(colors = colours, limits = limit, breaks = breaks)+
    theme(legend.key.size = unit(0.7, "lines"))
  ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol  = 8,nrow = 1,common.legend = T,legend = 'right')
}

limit=c(-1, 3)
breaks = c(-1,1, 3)
p1<-plot_gene(cluster='SMS') ###SAT1
p1

###### FigS6b,6d
SM <- subset(merged_SM_1, subset = orig.ident %in% c("SM35")) 
SM@meta.data$group1<-NA
SM@meta.data$group1[SM$integrated_snn_res.0.5%in% c("5","6","7","10")] <- "Tumor"
SM@meta.data$group1[SM$integrated_snn_res.0.5%in% c("1","2","3","4",'8',"9","11","12","13")] <- "Other"
Idents(SM)<-"group1"
my_comparisons <- list(c("Other", "Tumor"))
VlnPlot(SM, features = 'SAT1', pt.size = 0) &
  theme_bw() &
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_text(color = 'black', face = "bold", size = 12), 
    axis.text.y = element_text(color = 'black', face = "bold",size=12),       
    axis.title.y = element_text(color = 'black', face = "bold", size = 15),     
    panel.grid.major = element_blank(),     
    panel.grid.minor = element_blank(),      
    panel.border = element_rect(color = "black", size = 1.2, linetype = "solid"),     
    panel.spacing = unit(0.12, "cm"),      
    plot.title = element_text(hjust = 0.5, face = "bold.italic"),    
    legend.position = 'none'
  ) &
  stat_compare_means(
    method = "t.test", hide.ns = F,   
    comparisons = my_comparisons,                   
    label = "p.format",                   
    bracket.size = 0.8,                  
    tip.length = 0, 
    size = 6
  ) &
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) &
  scale_fill_manual(values = c('#4682B4', "#BF1D2D"))  # 设置小提琴图的颜色为红色和蓝色




##### FigS7
load("/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")
sample <- c('TJH08','TJH24','TJH26','TJH34','TJH37','TJH45','TJH46','TJH48','TJH66','TJH76','TJH90',"TJH14", "TJH17","TJH35","TJH55")
c1<-SpatialPlot(TJH[[1]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.4,image.alpha =0 ,label = TRUE)
c1
c1<-SpatialPlot(TJH[[5]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.4,image.alpha =0 ,label = TRUE)
c1
#region defination
TJH[[5]]@meta.data$Region<-NA
TJH[[5]]@meta.data$Region[TJH[[5]]@meta.data$SCT_snn_res.0.5 %in% c("4","5","6","7","8","11","12","14")] <- "Tumor"
TJH[[5]]@meta.data$Region[TJH[[5]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","2","3","10","13","9")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[5]])<-TJH[[5]]@meta.data$Region
p1 = SpatialPlot(TJH[[5]], group.by = 'Region',stroke = 0.1, pt.size.factor =2.0,cols = col2)
p1


c2<-SpatialPlot(TJH[[1]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c2
table(TJH[[1]]$SCT_snn_res.0.5,TJH[[1]]$cellType)
#region defination
TJH[[1]]@meta.data$Region<-NA
TJH[[1]]@meta.data$Region[TJH[[1]]@meta.data$SCT_snn_res.0.5 %in% c("3","4","5","6","8","10")] <- "Tumor"
TJH[[1]]@meta.data$Region[TJH[[1]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","2","7","9","11")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[1]])<-TJH[[1]]@meta.data$Region
p2 = SpatialPlot(TJH[[1]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p2

c3<-SpatialPlot(TJH[[2]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c3
table(TJH[[2]]$SCT_snn_res.0.5,TJH[[2]]$cellType)
#region defination
TJH[[2]]@meta.data$Region<-NA
TJH[[2]]@meta.data$Region[TJH[[2]]@meta.data$SCT_snn_res.0.5 %in% c("3","4","5","6","8","10")] <- "Tumor"
TJH[[2]]@meta.data$Region[TJH[[2]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","2","7","9","11")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[2]])<-TJH[[2]]@meta.data$Region
p3 = SpatialPlot(TJH[[2]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.9,cols = col2)
p3

c4<-SpatialPlot(TJH[[3]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c4
table(TJH[[3]]$SCT_snn_res.0.5,TJH[[3]]$cellType)
#region defination
TJH[[3]]@meta.data$Region<-NA
TJH[[3]]@meta.data$Region[TJH[[3]]@meta.data$SCT_snn_res.0.5 %in% c("1","4","6")] <- "Tumor"
TJH[[3]]@meta.data$Region[TJH[[3]]@meta.data$SCT_snn_res.0.5 %in% c("0","2","3","5","7","8","9")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[3]])<-TJH[[3]]@meta.data$Region
p4 = SpatialPlot(TJH[[3]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.9,cols = col2)
p4


c5<-SpatialPlot(TJH[[4]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c5
table(TJH[[4]]$SCT_snn_res.0.5,TJH[[4]]$cellType)
#region defination
TJH[[4]]@meta.data$Region<-NA
TJH[[4]]@meta.data$Region[TJH[[4]]@meta.data$SCT_snn_res.0.5 %in% c("1","3","4","5","10","11")] <- "Tumor"
TJH[[4]]@meta.data$Region[TJH[[4]]@meta.data$SCT_snn_res.0.5 %in% c("0","2","6","7","8","9","12")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[4]])<-TJH[[4]]@meta.data$Region
p5 = SpatialPlot(TJH[[4]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p5

c6<-SpatialPlot(TJH[[6]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c6
table(TJH[[6]]$SCT_snn_res.0.5,TJH[[6]]$cellType)
#region defination
TJH[[6]]@meta.data$Region<-NA
TJH[[6]]@meta.data$Region[TJH[[6]]@meta.data$SCT_snn_res.0.5 %in% c("2","3","4","12","15")] <- "Tumor"
TJH[[6]]@meta.data$Region[TJH[[6]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","5","6","7","8","9","10","11","13","14")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')
Idents(TJH[[6]])<-TJH[[6]]@meta.data$Region
p6 = SpatialPlot(TJH[[6]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p6


c7<-SpatialPlot(TJH[[7]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c7
table(TJH[[7]]$SCT_snn_res.0.5,TJH[[7]]$cellType)
#region defination
TJH[[7]]@meta.data$Region<-NA
TJH[[7]]@meta.data$Region[TJH[[7]]@meta.data$SCT_snn_res.0.5 %in% c("2","5","9")] <- "Tumor"
TJH[[7]]@meta.data$Region[TJH[[7]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","3","4","6","7","8","10","11","12","13")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[7]])<-TJH[[7]]@meta.data$Region
p7 = SpatialPlot(TJH[[7]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.9,cols = col2)
p7

c8<-SpatialPlot(TJH[[8]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c8
table(TJH[[8]]$SCT_snn_res.0.5,TJH[[8]]$cellType)
#region defination
TJH[[8]]@meta.data$Region<-NA
TJH[[8]]@meta.data$Region[TJH[[8]]@meta.data$SCT_snn_res.0.5 %in% c("4","7","11")] <- "Tumor"
TJH[[8]]@meta.data$Region[TJH[[8]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","2","3","5","6","8","9","10")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[8]])<-TJH[[8]]@meta.data$Region
p8 = SpatialPlot(TJH[[8]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p8

c9<-SpatialPlot(TJH[[9]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c9
table(TJH[[9]]$SCT_snn_res.0.5,TJH[[9]]$cellType)
#region defination
TJH[[9]]@meta.data$Region<-NA
TJH[[9]]@meta.data$Region[TJH[[9]]@meta.data$SCT_snn_res.0.5 %in% c("2","3","4","7","9","11","12")] <- "Tumor"
TJH[[9]]@meta.data$Region[TJH[[9]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","5","6","8","10")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[9]])<-TJH[[9]]@meta.data$Region
p9 = SpatialPlot(TJH[[9]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p9

c10<-SpatialPlot(TJH[[10]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c10
table(TJH[[10]]$SCT_snn_res.0.5,TJH[[10]]$cellType)
#region defination
TJH[[10]]@meta.data$Region<-NA
TJH[[10]]@meta.data$Region[TJH[[10]]@meta.data$SCT_snn_res.0.5 %in% c("2","3","4","6","8","10","12","14","15","17")] <- "Tumor"
TJH[[10]]@meta.data$Region[TJH[[10]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","5","7","9","11","13","16")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[10]])<-TJH[[10]]@meta.data$Region
p10 = SpatialPlot(TJH[[10]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p10

c11<-SpatialPlot(TJH[[11]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c11
table(TJH[[11]]$SCT_snn_res.0.5,TJH[[11]]$cellType)
#region defination
TJH[[11]]@meta.data$Region<-NA
TJH[[11]]@meta.data$Region[TJH[[11]]@meta.data$SCT_snn_res.0.5 %in% c("2","3","4","5","9","12")] <- "Tumor"
TJH[[11]]@meta.data$Region[TJH[[11]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","6","7","8","10","11","13","14")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')
Idents(TJH[[11]])<-TJH[[11]]@meta.data$Region
p11 = SpatialPlot(TJH[[11]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p11

c12<-SpatialPlot(TJH[[12]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c12
table(TJH[[12]]$SCT_snn_res.0.5,TJH[[12]]$cellType)
#region defination
TJH[[12]]@meta.data$Region<-NA
TJH[[12]]@meta.data$Region[TJH[[12]]@meta.data$SCT_snn_res.0.5 %in% c("2","5","6")] <- "Tumor"
TJH[[12]]@meta.data$Region[TJH[[12]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","3","4","7","8","9")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[12]])<-TJH[[12]]@meta.data$Region
p12 = SpatialPlot(TJH[[12]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p12

c13<-SpatialPlot(TJH[[13]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c13
table(TJH[[13]]$SCT_snn_res.0.5,TJH[[13]]$cellType)
#region defination
TJH[[13]]@meta.data$Region<-NA
TJH[[13]]@meta.data$Region[TJH[[13]]@meta.data$SCT_snn_res.0.5 %in% c("2","4","8")] <- "Tumor"
TJH[[13]]@meta.data$Region[TJH[[13]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","3","5","6","7","9","10","11")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[13]])<-TJH[[13]]@meta.data$Region
p13 = SpatialPlot(TJH[[13]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p13




c14<-SpatialPlot(TJH[[14]], group.by = 'cellType',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c14
table(TJH[[14]]$SCT_snn_res.0.5,TJH[[14]]$cellType)
#region defination
TJH[[14]]@meta.data$Region<-NA
TJH[[14]]@meta.data$Region[TJH[[14]]@meta.data$SCT_snn_res.0.5 %in% c("2","4","8")] <- "Tumor"
TJH[[14]]@meta.data$Region[TJH[[14]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","3","5","6","7")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[14]])<-TJH[[14]]@meta.data$Region
p14 = SpatialPlot(TJH[[14]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p14

c15<-SpatialPlot(TJH[[15]], group.by = 'SCT_snn_res.0.5',stroke = 0.15, pt.size.factor =1.8,image.alpha =0 ,label = TRUE)
c15
table(TJH[[15]]$SCT_snn_res.0.5,TJH[[15]]$cellType)
#region defination
TJH[[15]]@meta.data$Region<-NA
TJH[[15]]@meta.data$Region[TJH[[15]]@meta.data$SCT_snn_res.0.5 %in% c("5")] <- "Tumor"
TJH[[15]]@meta.data$Region[TJH[[15]]@meta.data$SCT_snn_res.0.5 %in% c("0","1","2","3","4","6","7","8","9","10","11")] <- "Other"
col2<-c('Other'='#FDB462','Tumor'='#DC050C')

Idents(TJH[[15]])<-TJH[[15]]@meta.data$Region
p15 = SpatialPlot(TJH[[15]], group.by = 'Region',stroke = 0.1, pt.size.factor =1.8,cols = col2)
p15



##### FigS8
#### FigS8a
load("/media/desk16/tly02/Project/CESC/merged_SM_8samples_merged_SM_1.rdata")
load("/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")
gene<-read.csv("/media/desk16/tly0202/TJ/cor/gene.csv")
gene<-gene$gene

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

#### FigS8b
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









#### FigS9
##### FigS9a
library(reticulate)
use_python("/data2/projects/bioinfo/zhluo/software/miniconda3/envs/stlearn/bin/python", required=TRUE)
py_config()


STModiCluster_1 <- function(InDir = InDir,
                          Sample = Sample,
                          OutDir = NULL,
                          TumorST = TumorST,
                          res = 1.5) {
  if (is.null(OutDir) == TRUE) {
    OutDir <- paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }

  # Adjusted_expr_mtx
  #reticulate::use_condaenv("base", required = TRUE)
  reticulate::source_python(system.file("python/Rusedtile.py", package = "Cottrazm"))
  Adjusted_expr_mtx <- ME_normalize(inDir = InDir, outDir = OutDir, sample = Sample)

  # Create Morph seu.obj
  aa_try <- try(
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts),
    silent = T
  )
  if (is(aa_try, "try-error")) {
    library(Matrix)
    Adjusted_expr_mtx <- Matrix::readMM(paste(OutDir, Sample, "_raw_SME_normalizeA.mtx", sep = ""))
  } else {
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts)
    colnames(Adjusted_expr_mtx) <- rownames(TumorST@assays$Spatial@counts)
  }

  Adjusted_expr_mtxF <- t(as.matrix(Adjusted_expr_mtx))
  MorphMatirxSeurat <- CreateSeuratObject(counts = as.matrix(Adjusted_expr_mtxF))

  # Add morph feature as assay to TumorST
  MorphMatirxSeurat <- subset(MorphMatirxSeurat, cells = rownames(TumorST@meta.data))
  TumorST@assays$Morph <- MorphMatirxSeurat@assays$RNA

  # Use Morph as assay Cluster
  TumorST <- NormalizeData(TumorST, assay = "Morph")
  TumorST <- FindVariableFeatures(object = TumorST, mean.function = ExpMean, dispersion.function = LogVMR, assay = "Morph")
  TumorST <- ScaleData(object = TumorST, assay = "Morph") # ,vars.to.regress = c('Mito.percent','Ribo.percent'))
  TumorST <- RunPCA(object = TumorST, npcs = 50, verbose = FALSE, assay = "Morph")
  TumorST <- FindNeighbors(TumorST, reduction = "pca", dims = 1:50, assay = "Morph")
  TumorST <- RunUMAP(object = TumorST, dims = 1:50, assay = "Morph")
  TumorST <- FindClusters(TumorST, resolution = res, algorithm = 1, graph.name = "Morph_snn")

  TumorST@meta.data$seurat_clusters <- TumorST@meta.data[, paste("Morph_snn_res.", res, sep = "")]
  # plot cluster result
  .cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
  )
  pdf(paste(OutDir, Sample, "_Spatial_SeuratCluster.pdf", sep = ""), width = 7, height = 7)
  p <- SpatialDimPlot(TumorST, group.by = "seurat_clusters", cols = .cluster_cols, pt.size.factor = 1, alpha = 0.8) +
    scale_fill_manual(values = .cluster_cols)+
    labs(title = paste("Resolution = ", res, sep = ""))
  print(p)
  dev.off()

  pdf(paste(OutDir, Sample, "_UMAP_SeuratCluster.pdf", sep = ""), width = 7, height = 7)
  p <- DimPlot(TumorST, group.by = "seurat_clusters", cols = .cluster_cols) + labs(title = paste("Resolution = ", res, sep = "")) +
    scale_fill_manual(values = .cluster_cols)
  print(p)
  dev.off()

  # Add ImmuneScore
  Normalfeatures <- c("PTPRC","CD2","CD3D","CD3E","CD3G","CD5","CD7","CD79A",'MS4A1',"CD19")
  TumorST@meta.data$NormalScore <- apply(TumorST@assays$Morph@data[rownames(TumorST@assays$Morph@data) %in% Normalfeatures, ], 2, mean)

  pdf(paste(OutDir, Sample, "_NormalScore.pdf", sep = ""), width = 6, height = 4)
  p <- VlnPlot(TumorST, features = "NormalScore", pt.size = 0, group.by = "seurat_clusters", cols = .cluster_cols) +
    geom_boxplot() +
    geom_hline(yintercept = max(unlist(lapply(
      split(TumorST@meta.data[, c("seurat_clusters", "NormalScore")], TumorST@meta.data[, c("seurat_clusters", "NormalScore")]$seurat_clusters),
      function(test) median(test$NormalScore)
    ))), linetype = "dashed") +
    ggpubr::stat_compare_means() + NoLegend()
  print(p)
  dev.off()

  NormalCluster <- levels(TumorST$seurat_clusters)[order(unlist(lapply(
    split(
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")],
      TumorST@meta.data[, c("seurat_clusters", "NormalScore")]$seurat_clusters
    ),
    function(test) mean(test$NormalScore)
  )), decreasing = T)[1]]
  print(paste("NormalCluster = ", NormalCluster, sep = ""))

  # save CNV annotation file
  cellAnnotation <- data.frame(CellID = rownames(TumorST@meta.data), DefineTypes = TumorST@meta.data[, "seurat_clusters"])
  dir.create(paste(OutDir, "InferCNV", sep = ""))
  write.table(cellAnnotation, paste(OutDir, "InferCNV/CellAnnotation.txt", sep = ""), sep = "\t", row.names = F, col.names = F, quote = F)

  return(TumorST)
}

input_data_dir = "/data2/projects/bioinfo/zhluo/data/wupeng/spaceranger"
sample_list = list.files(input_data_dir)
sample_list = sample_list[1:10]
##file path
samples_dir = sample_list %>% file.path(input_data_dir, ., "outs/")
sample_names = c("AM4", "AMZHZ","N1", "N5", "SM27", "SM32", "SM21", "SM25", "SM34", "SM35")
SME_objects =  purrr::map(1:length(sample_list), function(x){
    InDir = samples_dir[x]
    Sample = sample_names[x]
    OutDir = paste(getwd(),"/",Sample,"/",sep = "")
    dir.create(OutDir)
    
    print(InDir)
    print(Sample)
    print(OutDir)
    
    TumorST <-STPreProcess(
                InDir = InDir,
                OutDir = OutDir,
                Sample = Sample)
    
    res = 1.5
    TumorST <-STModiCluster_1(
                    InDir = InDir,
                    Sample = Sample,
                    OutDir = OutDir,
                    TumorST = TumorST,
                    res = res)
    
    return(TumorST)
})
readr::write_rds(SME_objects, path=file.path(getwd(), "all_sample_SME_seurat.rds.gz"),compress = "gz")

##### FigS9b

CSCC<-readRDS("/media/desk16/tly0202/TJ/all_sample_SME_seurat.rds.gz")
SM35 = CSCC[[10]]
SM35 <- SCTransform(SM35, assay = "Spatial", verbose = FALSE)
tumor_SM35 = subset(SM35,Morph_snn_res.1.5 %in% c(5,6,13,9))

df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
DefaultAssay(tumor_SM35) <- "SCT"
Idents(tumor_SM35) = "layer"
SpatialPlot(tumor_SM35)
saveRDS(tumor_SM35,file="tumor_SM35.rds")


#SM27
SM27 = CSCC[[5]]
SM27 <- SCTransform(SM27, assay = "Spatial", verbose = FALSE)
tumor_SM27 = subset(SM27,Morph_snn_res.1.5 %in% c(2,9,0))

df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
DefaultAssay(tumor_SM27) <- "SCT"
Idents(tumor_SM27) = "layer"
SpatialPlot(tumor_SM27)
saveRDS(tumor_SM27,file="SM27-tumor.rds")

##### FigS9c

extract_subset_focal_spot = function(seu, clusters=NULL){
    Idents(object = seu) <- "Morph_snn_res.1.5"
    focal_spot = subset(x = seu, idents = clusters)
    return(focal_spot)
}

##get the outline of the area

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

##将spot分割为区域
get_spot_section = function(selected_df){
    spot_area = list()
    for (row_idx in rownames(selected_df)){
        #idx = 0
        if (row_idx %in% unlist(spot_area))
            {next}

        spot_1 = selected_df[row_idx, ]

        idx = length(spot_area) + 1      
        spot_area[[idx]] = c(row_idx)

        new_spot = get_circle_n_spot(selected_df, spot=c(row=spot_1$row[1], col=spot_1$col[1]))
        if(nrow(new_spot) == 0){next}
        else{
            new_area = c(spot_area[[idx]], rownames(new_spot))
            while(length(new_area) > length(spot_area[[idx]]))
                {
                spot_area[[idx]] = new_area
                new_area_df = get_circle_all_spot(coor_df = selected_df, spot_all=selected_df[spot_area[[idx]],])
                new_area = unique(c(spot_area[[idx]], rownames(new_area_df)))
            }
        }

    }
    spot_area = data.frame(barcode = unlist(spot_area), section = rep(seq(length(spot_area)), lengths(spot_area)))
    df_section = left_join(rownames_to_column(selected_df), spot_area, by = c('rowname'="barcode"))
    return(df_section)
}


##提取坐标
coor_1  = tumor_SM35@images$image@coordinates
#转换坐标
xDiml <- max(SM35@images$image@coordinates$row)
yDiml <- max(SM35@images$image@coordinates$col)
coor_1$row = xDiml - coor_1$row
coor_1$cor_name = paste0(coor_1$row, "x", coor_1$col)
#coor_1$group = ifelse(coor_1$col >37 & coor_1$col <100 & coor_1$row <77 & coor_1$row > 41, "selected", "all")

##remove area with spot less than 10
#selected_df = coor_1[coor_1$group == "selected",]
##get area of focal
spot_area = get_spot_section(coor_1) # 这里的参数可以改成:selected_df, coor_1

##remove focal less than 10 spot
spot_filter = spot_area %>% group_by(section) %>% filter(n()>10)

##get outer line of focal
for (i in spot_filter$rowname){
    spot = c("row"=spot_filter$row[spot_filter$rowname==i], "col"=spot_filter$col[spot_filter$rowname==i])
    suround = get_circle_n_spot(spot_filter, spot)
    if(nrow(suround) < 6) {
        spot_filter[spot_filter$rowname==i, "group"] = "Boundary"
    }
    else{
        spot_filter[spot_filter$rowname==i, "group"] = "Inner"
    }
}

##Figure S9c
library(ggthemes)
p = ggplot(spot_filter, aes(x=col, y=row, colour =group)) + geom_point() + theme_base() +
    theme(
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
    labs(x="", y="", title = "SM35 tumor boundary")

##Figure S9d
##calculating distance
outer_spot = spot_filter$rowname[spot_filter$group == "Boundary"]
spot_distance = function(spot_filter, outer_spot){
    result = purrr::map(spot_filter$rowname, function(i){
        dist_n = 0
        new_spot = NULL
        if (i %in% outer_spot) {dist_n = 1 }
        if(dist_n == 0){
            spot_1 = spot_filter[spot_filter$rowname==i,]
            new_spot = get_circle_n_spot(spot_filter, spot=c(row=spot_1$row[1], col=spot_1$col[1]))
            if(length(intersect(new_spot$rowname,outer_spot))){dist_n = 2}
            
        }
        m=2
        while(dist_n == 0)
            {
            m = m + 1
            new_spot = get_circle_all_spot(coor_df = spot_filter, spot_all=new_spot)
            if(length(intersect(new_spot$rowname,outer_spot))){dist_n = m}
        }
        dist_n
    }) %>% unlist()
    return(result)
}
spot_filter = as.data.frame(spot_filter)
rownames(spot_filter) = spot_filter$rowname
spot_filter$distance = spot_distance(spot_filter=spot_filter, outer_spot=outer_spot)
spot_filter$distance = factor(spot_filter$distance, levels = c(1,2,3))
spot_filter = spot_filter[spot_filter$section != 5,]
spot_filter$section = factor(spot_filter$section, levels = unique(spot_filter$section))

p = ggplot(spot_filter, aes(x=col, y=row, colour =distance)) + geom_point()+ theme_base() +
    theme(
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
    labs(x="", y="", title = "Distance to boundary")

subset_spot = spot_filter %>% select(group, section, distance)
tumor_SM35_sub = tumor_SM35[,rownames(subset_spot)]
subset_spot = subset_spot[colnames(tumor_SM35_sub),]
tumor_SM35_sub@meta.data = cbind(tumor_SM35_sub@meta.data, subset_spot)


##Figure S9e
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35_sub@meta.data$layer = df_layer$Layer[match(tumor_SM35_sub@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
df_group_SM35 = tumor_SM35_sub@meta.data[,c("layer", "distance")] %>% 
group_by(layer, distance) %>% 
summarise(n = n())
df_group_SM35

df_group_SM35$distance = factor(df_group_SM35$distance, levels = c(1,2,3))
p1 = ggplot(df_group_SM35, aes(x=layer, y=n, fill=distance)) + 
  geom_bar(stat="identity", position="fill", color="black", lwd=0.25) + 
  labs(title="SM35 layer spots distance to boundary",  y="Frequence") +
  theme(axis.title.x = element_blank(),
       plot.title = element_text(size=12, hjust = 0.3))+ 
  scale_fill_manual(values=c("#FDE723", "#22958B", "#430154"))




####Figure S9f
coor_1  = tumor_SM27@images$image@coordinates
xDiml <- max(SM27@images$image@coordinates$row)
yDiml <- max(SM27@images$image@coordinates$col)
coor_1$row = xDiml - coor_1$row
coor_1$cor_name = paste0(coor_1$row, "x", coor_1$col)
#coor_1$group = ifelse(coor_1$col >37 & coor_1$col <100 & coor_1$row <77 & coor_1$row > 41, "selected", "all")

##remove area with spot less than 10
#selected_df = coor_1[coor_1$group == "selected",]

##get area of focal
spot_area = get_spot_section(coor_1) # 这里的参数可以改成:selected_df, coor_1

##remove focal less than 10 spot
spot_filter = spot_area %>% group_by(section) %>% filter(n()>10)

##get outer line of focal
for (i in spot_filter$rowname){
    spot = c("row"=spot_filter$row[spot_filter$rowname==i], "col"=spot_filter$col[spot_filter$rowname==i])
    suround = get_circle_n_spot(spot_filter, spot)
    if(nrow(suround) < 6) {
        spot_filter[spot_filter$rowname==i, "group"] = "Boundary"
    }
    else{
        spot_filter[spot_filter$rowname==i, "group"] = "Inner"
    }
}

library(ggthemes)
p = ggplot(spot_filter, aes(x=col, y=row, colour =group)) + geom_point() + theme_base() +
    theme(
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
    labs(x="", y="", title = "SM27 tumor boundary")



####Figure S9g
outer_spot = spot_filter$rowname[spot_filter$group == "Boundary"]
spot_distance = function(spot_filter, outer_spot){
    result = purrr::map(spot_filter$rowname, function(i){
        dist_n = 0
        new_spot = NULL
        if (i %in% outer_spot) {dist_n = 1 }
        if(dist_n == 0){
            spot_1 = spot_filter[spot_filter$rowname==i,]
            new_spot = get_circle_n_spot(spot_filter, spot=c(row=spot_1$row[1], col=spot_1$col[1]))
            if(length(intersect(new_spot$rowname,outer_spot))){dist_n = 2}
            
        }
        m=2
        while(dist_n == 0)
            {
            m = m + 1
            new_spot = get_circle_all_spot(coor_df = spot_filter, spot_all=new_spot)
            if(length(intersect(new_spot$rowname,outer_spot))){dist_n = m}
        }
        dist_n
    }) %>% unlist()
    return(result)
}
spot_filter = as.data.frame(spot_filter)
rownames(spot_filter) = spot_filter$rowname
spot_filter$distance = spot_distance(spot_filter=spot_filter, outer_spot=outer_spot)
spot_filter$distance = factor(spot_filter$distance, levels = c(1,2,3,4))
spot_filter = spot_filter[!spot_filter$section %in% c(14,16,18,9),]
spot_filter$section = factor(spot_filter$section, levels = unique(spot_filter$section))

p = ggplot(spot_filter, aes(x=col, y=row, colour =distance)) + geom_point()+ theme_base() +
    theme(
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank()) +
    labs(x="", y="", title = "Distance to boundary")


#### Figure S9h
subset_spot = spot_filter %>% select(group, section, distance)
tumor_SM27_sub = tumor_SM27[,rownames(subset_spot)]
subset_spot = subset_spot[colnames(tumor_SM27_sub),]
tumor_SM27_sub@meta.data = cbind(tumor_SM27_sub@meta.data, subset_spot)


##SM27
df_layer_27 = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27_sub@meta.data$layer = df_layer_27$Layer[match(tumor_SM27_sub@meta.data$seurat_clusters, df_layer_27$Oringe_cluster)]
df_group_SM27 = tumor_SM27_sub@meta.data[,c("layer", "distance")] %>% 
group_by(layer, distance) %>% 
summarise(n = n())


df_group_SM27$distance = factor(df_group_SM27$distance, levels = c(1,2,3,4))
p1 = ggplot(df_group_SM27, aes(x=layer, y=n, fill=distance)) + 
  geom_bar(stat="identity", position="fill", color="black", lwd=0.25) + 
  labs(title="SM27 layer spots distance to boundary",  y="Frequence") +
  theme(axis.title.x = element_blank(),
       plot.title = element_text(size=12, hjust = 0.3))+ 
  scale_fill_manual(values=c("#FDE723", "#22958B", "#430154", "#7E6148"))

#### Figure S9i,S9j
df_markers_35 = read_csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_DEG.csv")
df_markers_27 = read_csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM37_DEG.csv")
layer1_35 = df_markers_35[df_markers_35$cluster == "Layer1",]
layer1_27 = df_markers_27[df_markers_27$cluster == "Layer1",]
shared_gene = intersect(layer1_35$gene,layer1_27$gene)

head(layer1_35[match(shared_gene, layer1_35$gene),])
head(layer1_27[match(shared_gene, layer1_27$gene),])

df_shared_fc = data.frame(gene=shared_gene, 
                          log2FC_SM35=layer1_35$avg_log2FC[match(shared_gene, layer1_35$gene)], 
                          log2FC_SM27=layer1_27$avg_log2FC[match(shared_gene, layer1_27$gene)] )

cor.test(df_shared_fc$log2FC_SM27, df_shared_fc$log2FC_SM35)
p4 = ggplot(df_shared_fc, aes(x=log2FC_SM35,y=log2FC_SM27)) +
  geom_point(size=0.8) +
  #geom_abline(slope=1,intercept=0,color="red",size=1)+
  scale_color_uchicago() +
  annotate("text", x = 2, y = 4, label = "r=0.55, p-value < 2.2e-16",
           color="#350E20FF",size = 4 )+
  #geom_smooth(method=lm,color="red")+
  labs(x="log2FC_SM35",y="log2FC_SM27",title="Layer 1")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12, color="black"),
        title = element_text(size = 12,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
       plot.margin = margin(1,1,1,0.2, "cm"))

layer3_35 = df_markers_35[df_markers_35$cluster == "Layer3",]
layer3_27 = df_markers_27[df_markers_27$cluster == "Layer3",]
shared_gene = intersect(layer3_35$gene,layer3_27$gene)

head(layer3_35[match(shared_gene, layer3_35$gene),])
head(layer3_27[match(shared_gene, layer3_27$gene),])

df_shared_fc = data.frame(gene=shared_gene, 
                          log2FC_SM35=layer3_35$avg_log2FC[match(shared_gene, layer3_35$gene)], 
                          log2FC_SM27=layer3_27$avg_log2FC[match(shared_gene, layer3_27$gene)] )

cor.test(df_shared_fc$log2FC_SM27, df_shared_fc$log2FC_SM35)
p5 = ggplot(df_shared_fc, aes(x=log2FC_SM35,y=log2FC_SM27)) +
  geom_point(size=0.8) +
  #geom_abline(slope=1,intercept=0,color="red",size=1)+
  scale_color_uchicago() +
  annotate("text", x = 1, y = 2, label = "r=0.39, p-value = 8.4e-13",
           color="#350E20FF",size = 4 )+
  #geom_smooth(method=lm,color="red")+
  labs(x="log2FC_SM35",y="log2FC_SM27",title="Layer 3")+
  theme_classic(base_line_size = 1) +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size=12, color="black"),
        title = element_text(size = 12,
                             color = "black",
                             hjust = 0.5),
        axis.text.x = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 12,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        legend.title = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(0.6,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
       plot.margin = margin(1,1,1,0.2, "cm"))





#### FigS10
SM27<-readRDS("/media/desk16/tly0202/TJ/SM27-tumor.rds")
SM35<-reasRDS("/media/desk16/tly0202/TJ/tumor_SM35.rds")
colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100)
SpatialPlot(SM27, features =c("SPRR3","SPRR2A,"SFN","SERPINB4","PIK3R1,"COL7A1"),pt.size.factor = 1.6,images = 'SM27')+ scale_fill_gradientn(colors = colours)
SpatialPlot(SM35, features =c("SPRR3","SPRR2A,"SFN","SERPINB4","PIK3R1,"COL7A1"),pt.size.factor = 1.6,images = 'SM35')+ scale_fill_gradientn(colors = colours)

####FigS11
##### FigS.11a
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)
library(Seurat, lib.loc = "/media/desk16/tly0202/seurat_v4/")
.libPaths(c('~/seurat_v4/', .libPaths()))
load("/home/zhluo/Project/CESC/script/new_script/new_annotation_all_sample.rdata")
Idents(object = ST_list) <- "area"
SpatialPlot(ST_list[[11]],group.by = 'area',pt.size.factor = 1.8,label = T)

###### TJH90_sub
SpatialPlot(ST_list[[11]],group.by = 'seurat_clusters',pt.size.factor = 1.8,label = T)

#主要癌灶位于Cluster9
TJH90_sub<-subset(ST_list[[11]], seurat_clusters %in% '9')
SpatialPlot(TJH90_sub,group.by = 'area',pt.size.factor = 1.8,label = T)
table(ST_list[[11]]$area)

#这个癌灶太零散，上面的办法无法实现，换方法
#三个脚本，进行spot分割
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



get_spot_section = function(selected_df){
  spot_area = list()
  for (row_idx in rownames(selected_df)){
    #idx = 0
    if (row_idx %in% unlist(spot_area))
    {next}
    
    spot_1 = selected_df[row_idx, ]
    
    idx = length(spot_area) + 1      
    spot_area[[idx]] = c(row_idx)
    
    new_spot = get_circle_n_spot(selected_df, spot=c(row=spot_1$row[1], col=spot_1$col[1]))
    if(nrow(new_spot) == 0){next}
    else{
      new_area = c(spot_area[[idx]], rownames(new_spot))
      while(length(new_area) > length(spot_area[[idx]]))
      {
        spot_area[[idx]] = new_area
        new_area_df = get_circle_all_spot(coor_df = selected_df, spot_all=selected_df[spot_area[[idx]],])
        new_area = unique(c(spot_area[[idx]], rownames(new_area_df)))
      }
    }
    
  }
  spot_area = data.frame(barcode = unlist(spot_area), section = rep(seq(length(spot_area)), lengths(spot_area)))
  df_section = left_join(rownames_to_column(selected_df), spot_area, by = c('rowname'="barcode"))
  return(df_section)
}


# 提取坐标
coordinates <- TJH90_sub@images$TJH90@coordinates
# 获取坐标的行列最大值
xDiml <- max(coordinates$row)
yDiml <- max(coordinates$col)

# 转换坐标
coordinates$row <- xDiml - coordinates$row
coordinates$cor_name <- paste0(coordinates$row, "x", coordinates$col)

# 将 spot 分割为区域
spot_area <- get_spot_section(coordinates)
# 统计每种类型的数量
section_counts <- table(spot_area$section)
print(section_counts)

# 移除少于50 个 spot 的区域
spot_filter <- spot_area %>% group_by(section) %>% filter(n() > 50 )
section_counts <- table(spot_filter$section)
print(section_counts)

ggplot(spot_filter, aes(x = col, y = row, colour = section)) + geom_point() 


spot_filter$label<-'1'
spot_filter<-spot_filter[,c('rowname','label')]
rownames(spot_filter)<-spot_filter$rowname

#把确定的spot信息加入到TJH90_sub@meta.data
TJH90_sub<-AddMetaData(TJH90_sub,metadata =spot_filter)

#取出保留的spot
TJH90_sub<-subset(TJH90_sub, label %in% '1')

Layer1<-c('SPRR3','SPRR2A','SPRR2E','CNFN','S100A7','SPRR1B','SPRR2D','ANXA1','KRT6B','RHCG','EMP1','CD24')
Layer2<-c('ZFP36','KRT16','FGFBP1','SFN','NDUFA4L2','NTS','SERPINB4','MT2A','MT1E','MT1G','KRT15' ,'MIR205HG')
Layer3<-c('TNC','COL18A1','THBS2','COL7A1','PIK3R1','MMP11','STMN1','DPYSL3','IFITM3','GPC3','COL16A1','PTN')
gene<-c("MT2A","NME2","IFITM3")
TJH90_sub <- AddModuleScore(object = TJH90_sub, features = list(Layer1), name = "Layer1",nbin = 10)
TJH90_sub <- AddModuleScore(object = TJH90_sub, features = list(Layer2), name = "Layer2",nbin = 10)
TJH90_sub <- AddModuleScore(object = TJH90_sub, features = list(Layer3), name = "Layer3",nbin = 10)
TJH90_sub  <- AddModuleScore(object = TJH90_sub , features = list(gene), name = "EMT",nbin = 10)
SpatialPlot(TJH90_sub,features = c('Layer11'),pt.size.factor = 9) ###Layer1,Layer2,Layer3


TJH90_sub <- SCTransform(TJH90_sub, assay = "Spatial", verbose = FALSE)
DefaultAssay(TJH90_sub) <- "SCT"


library("cancersea")
data('available_pathways')

pathways <- list(angiogenesis=angiogenesis, apoptosis=apoptosis, cell_cycle=cell_cycle, differentiation=differentiation, dna_damage=dna_damage,
                 dna_repair=dna_repair, emt=emt, hypoxia=hypoxia, inflammation=inflammation, invasion=invasion, metastasis=metastasis, proliferation=proliferation,
                 quiescence=quiescence, stemness=stemness)
                 
TJH90_sub <- AddModuleScore(object = TJH90_sub, features = list(pathways[["dna_repair"]][["symbol"]]), name = "dna_repair",nbin = 10)
TJH90_sub<- AddModuleScore(object = TJH90_sub, features = list(pathways[["hypoxia"]][["symbol"]]), name = "hypoxia",nbin = 10)
TJH90_sub <- AddModuleScore(object = TJH90_sub, features = list(pathways[["proliferation"]][["symbol"]]), name = "proliferation",nbin = 10)
TJH90_sub<- AddModuleScore(object = TJH90_sub, features = list(pathways[["cell_cycle"]][["symbol"]]), name = "cell_cycle" ,nbin = 10)
TJH90_sub<- AddModuleScore(object = TJH90_sub, features = list(pathways[["invasion"]][["symbol"]]), name = "invasion" ,nbin = 10)
SpatialPlot(TJH90_sub,features = c('cell_cycle'),pt.size.factor = 9) ### cell_cycle,dna_repair,hypoxia,EMT

##### FigS11b
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

p = cor.test(TJH90_sub@meta.data$Layer11, TJH90_sub@meta.data$EMT1)
plot_correlation(p, TJH90_sub@meta.data$Layer11, TJH90_sub@meta.data$EMT1, "layer1", "EMT") #### cell_cycle,dna_repair,hypoxia,EMT

p = cor.test(TJH90_sub@meta.data$Layer31, TJH90_sub@meta.data$EMT1)
plot_correlation(p, TJH90_sub@meta.data$Layer31, TJH90_sub@meta.data$EMT1, "layer3", "EMT")#### cell_cycle,dna_repair,hypoxia,EMT
saveRDS(TJH90_sub,file="TJH90_sub.rds")

####  FigS12
##### FigS12a
tumor_SM27<-readRDS("/media/desk16/tly0202/TJ/SM27-tumor.rds")
tumor_SM35<-reasRDS("/media/desk16/tly0202/TJ/tumor_SM35.rds")
cds_1 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_trajectory.rds.gz")
cds_2 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM27_trajectory.rds.gz")
tumor_SM35@meta.data$Pseudotime = cds_1@phenoData@data$Pseudotime[match(rownames(tumor_SM35@meta.data), rownames(cds_1@phenoData@data))]
tumor_SM35@meta.data$Pseudotime = cds_1@phenoData@data$Pseudotime[match(rownames(tumor_SM35@meta.data), rownames(cds_1@phenoData@data))]
p1 = SpatialPlot(tumor_SM35, features = "Pseudotime", stroke = 0.5) + theme(legend.position="right")

tumor_SM27@meta.data$Pseudotime = cds_2@phenoData@data$Pseudotime[match(rownames(tumor_SM27@meta.data), rownames(cds_2@phenoData@data))]
tumor_SM27@meta.data$Pseudotime = cds_2@phenoData@data$Pseudotime[match(rownames(tumor_SM27@meta.data), rownames(cds_2@phenoData@data))]
p2 = SpatialPlot(tumor_SM27, features = "Pseudotime", stroke = 0.5) + theme(legend.position="right")
p = p1 + p2 + plot_layout(ncol=1)

gene_show = c("SPRR3", "SPRR2A", "S100A7", "VEGFA","TNC", "MMP11")
p1 = plot_genes_in_pseudotime(cds_1[gene_show,], color_by = "layer") + 
    scale_color_manual(values=layer_colors) + 
    labs(title="Genes trajection of SM35")
    
p1 = plot_genes_in_pseudotime(cds_2[gene_show,], color_by = "layer") + 
    scale_color_manual(values=layer_colors[2:4]) + 
    labs(title="Genes trajection of SM27")

#### Fig13
##### Fig13a
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
tumor_SM35<-reasRDS("/media/desk16/tly0202/TJ/tumor_SM35.rds")
tumor_SM27<-readRDS("/media/desk16/tly0202/TJ/SM27-tumor.rds")
sce_SCENIC <- open_loom("/data2/projects/bioinfo/zhluo/project/wupeng/scenic/output_epi/merged_epi_SCENIC.loom")

regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
##SM35
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("layer1", "layer2", "layer2", "layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
cellinfo <- tumor_SM35@meta.data[,c('orig.ident','seurat_clusters', "layer", "nFeature_Spatial","nCount_Spatial")]#细胞meta信息

cellTypes <-  as.data.frame(subset(cellinfo,select = c('layer', "seurat_clusters")))
cellTypes$layer =paste("SM35_", cellTypes$layer, sep="")
rownames(cellTypes) =paste("8_", rownames(cellTypes), sep="")


selectedResolution <- "layer"
sub_regulonAUC <- regulonAUC
cells = intersect(colnames(sub_regulonAUC), rownames(cellTypes))
cellTypes = cellTypes[cells,]
sub_regulonAUC_SM35 = sub_regulonAUC[,cells]

rss <- calcRSS(AUC=getAUC(sub_regulonAUC_SM35),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC_SM35), selectedResolution])
rss=na.omit(rss)

rssPlot <- 
  plotRSS(
  rss,
  zThreshold = 1.2,
  cluster_columns = TRUE,
  order_rows = TRUE,
  thr=0.1,
  varName = "layer",
  col.low = '#330066',
  col.mid = '#66CC66',
  col.high = '#FFCC33')


##SM27
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
cellinfo <- tumor_SM27@meta.data[,c('orig.ident','seurat_clusters', "layer", "nFeature_Spatial","nCount_Spatial")]#细胞meta信息

cellTypes <-  as.data.frame(subset(cellinfo,select = c('layer', "seurat_clusters")))
cellTypes$layer =paste("SM27_", cellTypes$layer, sep="")
rownames(cellTypes) =paste("3_", rownames(cellTypes), sep="")


selectedResolution <- "layer"
sub_regulonAUC <- regulonAUC
cells = intersect(colnames(sub_regulonAUC), rownames(cellTypes))
cellTypes = cellTypes[cells,]
sub_regulonAUC_SM27 = sub_regulonAUC[,cells]

rss <- calcRSS(AUC=getAUC(sub_regulonAUC_SM27),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC_SM35), selectedResolution])
rss=na.omit(rss)

rssPlot_27 <- 
  plotRSS(
  rss,
  zThreshold = 1.2,
  cluster_columns = TRUE,
  order_rows = TRUE,
  thr=0.1,
  varName = "layer",
  col.low = '#330066',
  col.mid = '#66CC66',
  col.high = '#FFCC33')

library(ggheatmap)
library(reshape2)

rss_data_35 <- rssPlot$plot$data
rss_data_35<-dcast(rss_data_35, 
                Topic~rss_data_35$layer,
                value.var = 'Z')
rownames(rss_data_35) <- rss_data_35[,1]
rss_data_35 <- rss_data_35[,-1]

rss_data_27 <- rssPlot_27$plot$data
rss_data_27<-dcast(rss_data_27, 
                Topic~rss_data_27$layer,
                value.var = 'Z')
rownames(rss_data_27) <- rss_data_27[,1]
rss_data_27 <- rss_data_27[,-1]

rss_data <- merge(rss_data_35, rss_data_27, by=0, all=FALSE)
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]

colnames(rss_data)
col_ann <- data.frame(Group= c("Layer 2", "Layer 3", "Layer 1", "Layer 3", "Layer 2", "Layer 1"),
                      Sample = c("SM35", "SM35", "SM35", "SM27", "SM27", "SM27"))#列注释
rownames(col_ann) <- colnames(rss_data)

groupcol <- c( "#3A8E7B", "#F4C4AA", "#B4EDF4")
names(groupcol) <- c("Layer 1","Layer 2", "Layer 3")

samplecol <- c("#98D352","#FF7F0E")
names(samplecol) <- c("SM35","SM27")
col <- list(Group=groupcol, Sample=samplecol)

text_columns <- sample(colnames(rss_data),0)#不显示列名

p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = T,
               show_cluster_rows = F, show_cluster_cols = F,
               #scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)

### FigS13b 
tf_list = read.csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/TF/Homo_sapiens_TF.txt", sep="\t")
tf_list = tf_list[,c("Species", "Symbol")]
tf_list = tf_list[!tf_list$Symbol == "", ] 

data_markers_SM27 = read.csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM37_DEG.csv")
data_markers_SM35 = read.csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_DEG.csv")
data_markers_SM27 = data_markers_SM27[data_markers_SM27$gene %in% tf_list$Symbol,]
data_markers_SM35 = data_markers_SM35[data_markers_SM35$gene %in% tf_list$Symbol,]

Cluster1_TF = intersect(data_markers_SM27$gene[data_markers_SM27$cluster == "Layer1"], data_markers_SM35$gene[data_markers_SM35$cluster == "Layer1"])
Cluster2_TF = intersect(data_markers_SM27$gene[data_markers_SM27$cluster == "Layer2"], data_markers_SM35$gene[data_markers_SM35$cluster == "Layer2"])
Cluster3_TF = intersect(data_markers_SM27$gene[data_markers_SM27$cluster == "Layer3"], data_markers_SM35$gene[data_markers_SM35$cluster == "Layer3"])
all_tf = c(Cluster1_TF, Cluster2_TF, Cluster3_TF)
#构建表达矩阵
sub_expr = tumor_SM35@assays$SCT@counts[all_tf,]
SM35_cell_layer1 = rownames(tumor_SM35@meta.data[tumor_SM35@meta.data$layer == "Layer1",])
SM35_cell_layer2 = rownames(tumor_SM35@meta.data[tumor_SM35@meta.data$layer == "Layer2",])
SM35_cell_layer3 = rownames(tumor_SM35@meta.data[tumor_SM35@meta.data$layer == "Layer3",])
SM35_res_1 = apply(sub_expr[,SM35_cell_layer1], 1, mean) %>% as.data.frame()
SM35_res_2 = apply(sub_expr[,SM35_cell_layer2], 1, mean) %>% as.data.frame()
SM35_res_3 = apply(sub_expr[,SM35_cell_layer3], 1, mean) %>% as.data.frame()
SM35_res = cbind(SM35_res_1, SM35_res_2, SM35_res_3)
colnames(SM35_res) = c("SM35_layer1", "SM35_layer2", "SM35_layer3")

sub_expr = tumor_SM27@assays$SCT@counts[all_tf,]
SM27_cell_layer1 = rownames(tumor_SM27@meta.data[tumor_SM27@meta.data$layer == "Layer1",])
SM27_cell_layer2 = rownames(tumor_SM27@meta.data[tumor_SM27@meta.data$layer == "Layer2",])
SM27_cell_layer3 = rownames(tumor_SM27@meta.data[tumor_SM27@meta.data$layer == "Layer3",])
SM27_res_1 = apply(sub_expr[,SM27_cell_layer1], 1, mean) %>% as.data.frame()
SM27_res_2 = apply(sub_expr[,SM27_cell_layer2], 1, mean) %>% as.data.frame()
SM27_res_3 = apply(sub_expr[,SM27_cell_layer3], 1, mean) %>% as.data.frame()
SM27_res = cbind(SM27_res_1, SM27_res_2, SM27_res_3)
colnames(SM27_res) = c("SM27_layer1", "SM27_layer2", "SM27_layer3")
tf_expre = cbind(SM35_res, SM27_res)
tf_expre = tf_expre[all_tf,]


col_ann <- data.frame(Group= c("Layer 1", "Layer 2", "Layer 3", "Layer 1", "Layer 2", "Layer 3"),
                      Sample = c("SM35", "SM35", "SM35", "SM27", "SM27", "SM27"))#列注释
rownames(col_ann) <- colnames(tf_expre)
groupcol <- c( "#3A8E7B", "#F4C4AA", "#B4EDF4")
names(groupcol) <- c("Layer 1","Layer 2", "Layer 3")
samplecol <- c("#98D352","#FF7F0E")
names(samplecol) <- c("SM35","SM27")
col <- list(Group=groupcol, Sample=samplecol)
text_columns <- sample(colnames(tf_expre),0)#不显示列名

p <- ggheatmap(tf_expre,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = T,cluster_cols = T,
               show_cluster_rows = F, show_cluster_cols = F,
               scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,
               legendName="Relative value",
               text_show_cols = text_columns)

#### Fig13c-g
cds_1 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_trajectory.rds.gz")
cds_2 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM27_trajectory.rds.gz")
pData(cds_1)$ELF3 = log2(exprs(cds_1)["ELF3",] + 1)
p7 = plot_cell_trajectory(cds_1, cell_size=2,color_by = "ELF3") +
     scale_colour_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24",
                       midpoint = 3)#, breaks=c(-0.05,0.05,0.15))

pData(cds_2)$ELF3 = log2(exprs(cds_2)["ELF3",] + 1)
p8 = plot_cell_trajectory(cds_2, cell_size=2,color_by = "ELF3") +
     scale_colour_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24",
                       midpoint = 3)#, breaks=c(-0.05,0.05,0.15))

gene= c("ELF3")
p9 = SpatialPlot(tumor_SM35, features = gene, stroke = 0.5) + theme(legend.position="top")
p10 = SpatialPlot(tumor_SM27, features = gene, stroke = 0.5) + theme(legend.position="top")

p = p7 + p8 + p9 + p10 + plot_layout(ncol=2)

#### FigS13h
TF_Genes_score= read.csv("All_TF_targets_Genes_score.csv")
sel_TFs = c("ELF3","GRHL1","PRDM1","SOX2","TEAD2")
sel_TF_Genes0 = TF_Genes_score[TF_Genes_score$TF%in%sel_TFs, ]
sel_TF_Genes = sel_TF_Genes0[sel_TF_Genes0$Score>1, ]
sel_TF_Genes = rbind(sel_TF_Genes, sel_TF_Genes0[sel_TF_Genes0$Gene%in%sel_TFs,])
sel_TF_Genes = unique(sel_TF_Genes)
head(sel_TF_Genes)
gr <- sel_TF_Genes %>% graph_from_data_frame(directed = T)
V(gr)$type= V(gr)$name
V(gr)$type[V(gr)$type%in%sel_TFs] ='TF'
V(gr)$type[V(gr)$type%in%sel_TF_Genes$Gene] ='Target'
V(gr)$size[V(gr)$type=="TF"] =2
V(gr)$size[V(gr)$type=="Target"] =1
V(gr)$color ='white'
p = ggraph(gr, layout ='sugiyama') + geom_edge_link(aes(color ="#D6404E"), show.legend = F) + geom_node_point(color ='white')#网络图
p
pData = p$data
pData = pData[rev(order(pData$type)),]
pData$color[1:length(tf_col)] <- tf_col
pData$color[pData$color=="white"] <-'grey'
p + geom_point(data=pData,aes(x,y,color=color,size=size,stroke=1), show.legend = F) +scale_color_manual(values=pData$color) + geom_text(data=subset(pData,type=='Target'),aes(x,y,label=name),size=3,fontface="italic", angle=45, hjust=1)+ geom_text(data=subset(pData,type=='TF'),aes(x,y,label=name), size=4,fontface="bold")+ theme_graph()+scale_y_discrete(expand=expansion(mult=c(0.5,0.05)))
p
p = ggraph(gr, layout ='fr') + geom_edge_link(aes(color ="#D6404E"), show.legend = F) + geom_node_point(color ='white')#网络图
pData = p$data
pData = pData[rev(order(pData$type)),]
pData$color[1:length(tf_col)] <- tf_col
pData$color[pData$color=="white"] <-'grey'

p = p + geom_point(data=pData,aes(x,y,color=color,size=size,stroke=1), show.legend = F) + scale_color_manual(values=pData$color) + geom_text(data=subset(pData,type=='Target'),aes(x,y,label=name),size=2.5,fontface="italic")+geom_text(data=subset(pData,type=='TF'),aes(x,y,label=name), size=3,fontface="bold")+ theme_graph()
p



#### FigS14
extract_subset_focal_spot = function(seu, clusters=NULL){
    Idents(object = seu) <- "Morph_snn_res.1.5"
    focal_spot = subset(x = seu, idents = clusters)
    return(focal_spot)
}
h_gene_sets = msigdbr(species = "human", category = "H")
all_hallmarks = unique(h_gene_sets$gs_name)
SME_objects = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/Cottrazm/all_sample_SME_seurat.rds.gz")
SM35 = SME_objects[[10]]
SM35 <- SCTransform(SM35, assay = "Spatial", verbose = FALSE)
tumor_SM35 = extract_subset_focal_spot(SM35,c(11, 12, 13, 5, 6, 9))
df_layer = data.frame(Oringe_cluster = c(11,12, 13, 5, 6, 9), Layer = c("Normal", "Normal", "Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
tumor_SM35@meta.data$layer = factor(tumor_SM35@meta.data$layer, levels=c("Normal", "Layer1", "Layer2", "Layer3"))
x ="ANDROGEN_RESPONSE" ####BILE_ACID_METABOLISM
one_features = list(h_gene_sets$gene_symbol[h_gene_sets$gs_name == x])
tumor_SM35 = AddModuleScore(object = tumor_SM35, features = one_features, name = x)
cds_1 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM35_trajectory.rds.gz")
module_score = tumor_SM35@meta.data[, paste0(x, "1")]
pData(cds_1)[[x]] = module_score
p4 = plot_cell_trajectory(cds_1, cell_size=2,color_by = "ANDROGEN_RESPONSE") +
     scale_colour_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24")
     #scale_colour_gradient2()
p4
     
SM27 = SME_objects[[5]]
SM27 <- SCTransform(SM27, assay = "Spatial", verbose = FALSE)
tumor_SM27 = extract_subset_focal_spot(SM27,c(2, 0, 9))
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
tumor_SM27@meta.data$layer = factor(tumor_SM27@meta.data$layer, levels=c("Layer1", "Layer2", "Layer3"))
x ="ANDROGEN_RESPONSE" ####BILE_ACID_METABOLISM
one_features = list(h_gene_sets$gene_symbol[h_gene_sets$gs_name == x])
tumor_SM35 = AddModuleScore(object =tumor_SM27, features = one_features, name = x)
cds_2 = readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/sample_compare/DEG_normal_tumor/SM27_trajectory.rds.gz")
module_score = tumor_SM27@meta.data[, paste0(x, "1")]
pData(cds_2)[[x]] = module_score
p4 = plot_cell_trajectory(cds_2, cell_size=2,color_by = "ANDROGEN_RESPONSE") +
     scale_colour_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24")
     #scale_colour_gradient2()
p4







#### FigS15
##### FigS15a-d
library(scMetabolism)
library(tidyverse)
library(rsvd)
#### load data
TJH90_sub <- readRDS("~/TJ/score/TJH90_sub.rds")
TJH37_sub <- readRDS("~/TJ/score/TJH37_sub.rds")
SM27<-readRDS("/media/desk16/tly0202/TJ/SM27-tumor.rds")
SM35<-reasRDS("/media/desk16/tly0202/TJ/tumor_SM35.rds")


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

#####  FigS15a,b,SM27 or SM35 
SM27 <- sc.metabolism.Seurat(obj = SM27, 
                                    method = "VISION",     
                                    imputation = F, ncores = 2, 
                                    metabolism.type = "KEGG")
metabolism.matrix <- SM27@assays$METABOLISM$score
metabolism.matrix
score<-as.data.frame(t(metabolism.matrix))
SM27<-AddMetaData(SM27,metadata =score)
SpatialPlot(SM27,features = c('Fatty acid biosynthesis'),pt.size.factor = 9) 


#### FigS15c,d TJH90_sub or TJH37_sub
TJH90_sub1 <- sc.metabolism.Seurat(obj = TJH90_sub, 
                                    method = "VISION",     
                                    imputation = F, ncores = 2, 
                                    metabolism.type = "KEGG")
metabolism.matrix <- TJH90_sub1@assays$METABOLISM$score
metabolism.matrix
TJH90_score<-as.data.frame(t(metabolism.matrix))
TJH90_sub1<-AddMetaData(TJH90_sub1,metadata =TJH90_score)
SpatialPlot(TJH90_sub1,features = c('Fatty acid biosynthesis'),pt.size.factor = 9) 





##### FigS15e-f
plot_correlation <- function(correlation_result, x, y, x_label, y_label) {
  ggplot(TJH90_sub1@meta.data, aes(x = x, y = y)) +
    geom_point(alpha = 1, col = '#4682B4') +  # 绘制散点，alpha 设置透明度
    stat_smooth(method = "lm", col = "red") +  # 添加回归线
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

p = cor.test(TJH90_sub1@meta.data$Layer11, TJH90_sub1@meta.data$`Fatty acid biosynthesis`)
plot_correlation(p, TJH90_sub1@meta.data$Layer11, TJH90_sub1@meta.data$`Fatty acid biosynthesis`, "layer1", "Fatty acid biosynthesis")

p = cor.test(TJH90_sub1@meta.data$Layer31, TJH90_sub1@meta.data$`Fatty acid biosynthesis`)
plot_correlation(p, TJH90_sub1@meta.data$Layer31, TJH90_sub1@meta.data$`Fatty acid biosynthesis`, "layer3", "Fatty acid biosynthesis")