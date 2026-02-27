#### laod  inneed package


library(Seurat)
packageVersion("Seurat")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(tidydr)
library(RColorBrewer)

#SpaceRanger-out
input_data_dir <-"/data2/projects/bioinfo/zhluo/data/wupeng/spaceranger"
sample_list <- list.files(input_data_dir)
sample_list

sample_list <- sample_list[c(3:10)]
sample_list

## file path
samples_dir<- sample_list %>% file.path(input_data_dir, .,"outs")
samples_dir
##slice id
sample_names <- c("N1", "N5", "SM27","SM32","SM21", "SM25", "SM34","SM35")

## read samples,creat seurat object
sample_objects <- purrr::map(1:length(sample_list), function(x) {
  ## read data
  one_dir <- samples_dir[x]
  sample_id <- sample_list[x]
  slice_id <- sample_names[x]
  sample_object <- Load10X_Spatial(
    data.dir = one_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = slice_id,
    filter.matrix = TRUE
  )
  sample_object@project.name <- sample_id
  sample_object@meta.data$orig.ident <- slice_id
  sample_object <- RenameCells(object = sample_object, add.cell.id = slice_id)
  
  ##quality filter
  sample_object = PercentageFeatureSet(sample_object, pattern="^MT-", col.name="percent.mt")
  sample_object = sample_object[, (sample_object$nFeature_Spatial > 200) & (sample_object$percent.mt < 25)]
  ##UMI >=3
  genes_obtained = names(rowSums(sample_object@assays$Spatial@counts)[rowSums(sample_object@assays$Spatial@counts) >= 3])
  ##remove mito genes and ribo genes
  mito.genes <- grep('^MT-', genes_obtained, value = TRUE)
  ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', genes_obtained, value = TRUE)
  genes_obtained = genes_obtained[!genes_obtained %in% c(mito.genes, ribo.genes)]
  sample_object <- sample_object[rownames(sample_object) %in% genes_obtained, ]
  ##normalization
  sample_object = sample_object %>% SCTransform(assay="Spatial", verbose = F) %>% RunPCA(assay = "SCT", verbose = F) %>% FindNeighbors(reduction = "pca", dims = 1:30, verbose = F) %>% FindClusters(verbose = F) %>% RunUMAP(dims = 1:30, verbose = F)
    return(sample_object)
})

####merge data，integration
options(future.globals.maxSize = 30000 * 1024^2) 
genes_name_list = lapply(sample_objects, rownames)
genes.common = Reduce(intersect, genes_name_list)
st.features = SelectIntegrationFeatures(sample_objects, nfeatures = 2000, verbose = FALSE, assay = rep("SCT", length(sample_objects)))
sample_objects <- PrepSCTIntegration(object.list = sample_objects, anchor.features = st.features,verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = sample_objects, normalization.method = "SCT",verbose = FALSE, anchor.features = st.features)
ST <- IntegrateData(anchorset = int.anchors, features.to.integrate = genes.common, normalization.method = "SCT", verbose = FALSE)


DefaultAssay(ST) <- "integrated"

#VariableFeatures duplicate removal
uniFeature <- unique(unlist(lapply(sample_objects, VariableFeatures)))
VariableFeatures(ST) <- uniFeature

#PCA
ST<-RunPCA(ST, features = VariableFeatures(object = ST), verbose = FALSE)

library(harmony)
###remove patch
ST <- RunHarmony(ST,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")

#FindClusters
ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5)
ST <- RunUMAP(ST, dims = 1:30, reduction = "harmony")


merged_SM_1<-ST
DefaultAssay(merged_SM_1) <- "integrated"
merged_SM_1@meta.data$integrated_snn_res.0.5 <- as.numeric(as.character(merged_SM_1@meta.data$integrated_snn_res.0.5)) + 1
merged_SM_1@meta.data$integrated_snn_res.0.5 <- as.factor(merged_SM_1@meta.data$integrated_snn_res.0.5)

data.markers <- FindAllMarkers(merged_SM_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#####Fig1B----
col = c('1'="#7fc97f", '2'="#beaed4", '3'="#fdc086", '4'="#386cb0", '5'="#e0abcf", '6'="#a34e3b", '7'="#efbd25", '8'="#1b9e77", '9'="#d95f02", '10'="#6c67a5", '11'="#d01b2a", '12'="#43acde", '13'="#666666")

umap<-DimPlot(merged_SM_1,cols = col,group.by = 'integrated_snn_res.0.5')+theme_dr(xlength=0.22,ylength=0.22, arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+theme(panel.grid=element_blank())
umap

#####Fig1C cellmarker dotplot
Idents(object = merged_SM_1) <- "integrated_snn_res.0.5"
list<-c('CD79A','IGLC2','IGHG4','FAP', 'COL16A1','GPC6','P4HA3','FAM83D','CEP85','ZNF891','ALOX12B','KRT17','VEGFA','TM4SF1','KRT15','MIR205HG','MT1X','KRT6A','KRT6B','KRT14','ADAM12','POSTN','MMP11','FN1','C1QB','CD163','CD68',"KRT7","ECM1","CNFN","CEACAM1","VWF","PECAM1","PLVAP","MYH11","ACTA2","MUSTN1",'VIM',"HBB","HBA1","HBA2")
dot <- DotPlot(merged_SM_1, features = list,col.max = 1.5) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  
    axis.text.y = element_text(size = 10),                        
    axis.title.x = element_text(size = 12),                        
    axis.title.y = element_text(size = 12),                        
    legend.title = element_text(size = 10),                        
    legend.text = element_text(size = 10)
    )+ 
 ylab('Clusters') +  
  scale_color_gradientn(colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100), breaks = c(1, 0, -1),                                           # 设置显示的刻度值
                        labels = c("2", "0", "-1")))
dot


#####Fig1D cellRatio script 
cellRatioPlot=function (object = NULL, sample.name = NULL, celltype.name = NULL, 
                        col.width = 0.7, flow.alpha = 0.25, flow.curve = 0, fill.col = NULL) 
{
  meta <- object@meta.data
  ratio.info <- meta %>% dplyr::group_by(.data[[sample.name]], 
                                         .data[[celltype.name]]) %>% dplyr::summarise(num = n()) %>% 
    dplyr::mutate(rel_num = num/sum(num))
  if (is.null(fill.col)) {
    fill.col <- jjAnno::useMyCol("paired", n = length(unique(meta[, 
                                                                  celltype.name])))
  }
  else {
    fill.col <- fill.col
  }
  p <- ggplot2::ggplot(ratio.info, ggplot2::aes_string(x = sample.name, 
                                                       y = "rel_num")) + ggplot2::geom_col(ggplot2::aes_string(fill = celltype.name), 
                                                                                           width = col.width) + ggalluvial::geom_flow(ggplot2::aes_string(stratum = celltype.name, 
                                                                                                                                                          alluvium = celltype.name, fill = celltype.name), width = 0.5, 
                                                                                                                                      alpha = flow.alpha, knot.pos = flow.curve) + ggplot2::theme_bw() + 
    ggplot2::coord_cartesian(expand = 0) + ggplot2::scale_y_continuous(labels = scales::label_percent()) + 
    ggplot2::scale_fill_manual(values = fill.col, name = "Clusters") + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_text(size = ggplot2::rel(1.8), 
                                                     color = "black"), axis.title = ggplot2::element_text(size = ggplot2::rel(2), 
                                                                                                          color = "black"), legend.text = ggplot2::element_text(size = ggplot2::rel(1.8), 
                                                                                                                                                                color = "black"), legend.title = ggplot2::element_text(size = ggplot2::rel(2), 
                                                                                                                                                                                                                       color = "black"), plot.margin = ggplot2::margin(t = 20, l = 30, r = 20, unit = "pt")) + 
    ggplot2::xlab("") + ggplot2::ylab("Cell percent ratio")
  return(p)
}

cellRatio<-cellRatioPlot(object = merged_SM_1,
                         sample.name = "orig.ident",
                         celltype.name = "integrated_snn_res.0.5",
                         flow.curve = 0.5,fill.col = col)+
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5))+xlab('Sample')
cellRatio 


#####Fig1E RTCD
###RCTD                         
#single cell
ColonData <- readr::read_rds("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/BD_20samples/ColonData_Harmony.rds.gz")
meta.data = read.csv("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/BD_20samples/meta.data.csv", row.names = 1)
cell_types = as.factor(meta.data$AllTypes)
names(cell_types) = rownames(meta.data)
reference <- Reference(ColonData@assays$RNA@counts, cell_types=cell_types)

input_data_dir = "/data2/projects/bioinfo/zhluo/data/wupeng/spaceranger"
sample_list = list.files(input_data_dir)
sample_list = sample_list[3:10]
##loop RCTD
purrr::map(1:length(sample_list), function(x){
samples_dir = sample_list %>% file.path(input_data_dir, ., "outs")
samples_dir[x]
puck = read.VisiumSpatialRNA(samples_dir[x])
#RCTD
#doublet
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
#full
full_myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
full_file_name = paste("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/RCTD/", sample_list[x], ".full.rds.gz",sep="")
readr::write_rds(full_myRCTD, file=full_file_name,compress = "gz")
})

add_cellProp_assay <- function(visium_slide){
    #visium_slide
    #read cell type prop
    files = list.files("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/RCTD/")
    full_files = grep('full', files, value=TRUE)
    full_files = full_files %>% file.path("/data2/projects/bioinfo/zhluo/project/wupeng/SM_sample/res_05/SpaCET/notebook/singleCell/RCTD/", .)
    cell_type_dict = c(1,2,3,4,5,6,7,8)
    names(cell_type_dict) = c("N1", "N5", "SM27","SM32","SM21", "SM25", "SM34","SM35")
    
    sample_name = names(visium_slide@images)
    file_path = full_files[cell_type_dict[sample_name]]
    full_myRCTD = readr::read_rds(file_path)
    barcodes <- colnames(full_myRCTD@spatialRNA@counts)
    weights <- full_myRCTD@results$weights
    norm_weights <- normalize_weights(weights) #col: cell trpe, row: barcode 
    norm_weights = t(norm_weights)

    #add assay
    adt_assay = CreateAssayObject(counts = norm_weights[,colnames(visium_slide)])
    visium_slide[["RCTD"]] = adt_assay  
    return(visium_slide)
  }
  
sample_objects= purrr::map(sample_objects, add_cellProp_assay)


add_combined_RCTD_assay <- function(integrated_obj, sample_objects_list) {

  all_rctd_data <- list()
  
  for (sample_obj in sample_objects_list) {
    sample_name <- unique(sample_obj$orig.ident)
    
    if ("RCTD" %in% Assays(sample_obj)) {

      rctd_data <- GetAssayData(sample_obj, assay = "RCTD", slot = "counts")
      

      colnames(rctd_data) <- paste0(sample_name, "_", colnames(rctd_data))
      
      all_rctd_data[[sample_name]] <- rctd_data
    }
  }
  

  if (length(all_rctd_data) > 0) {

    common_celltypes <- Reduce(intersect, lapply(all_rctd_data, rownames))
    
    if (length(common_celltypes) > 0) {
      combined_rctd <- do.call(cbind, lapply(all_rctd_data, function(x) {
        x[common_celltypes, , drop = FALSE]
      }))
      

      common_cells <- intersect(colnames(integrated_obj), colnames(combined_rctd))
      combined_rctd <- combined_rctd[, common_cells]
      
      rctd_assay <- CreateAssayObject(counts = combined_rctd)
      

      integrated_obj[["RCTD"]] <- rctd_assay
      
      cat("Added RCTD assay with", length(common_celltypes), 
          "cell types and", length(common_cells), "cells\n")
    }
  }
  
  return(integrated_obj)
}

merged_SM_1 <- add_combined_RCTD_assay(merged_SM_1, sample_objects)


cell_prop = as.data.frame(merged_SM_1@assays[["RCTD"]]@counts) %>% t() %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  pivot_longer(-spot_id)
cell_prop$id<-cell_prop$spot_id

#head(cell_prop)
niche_info = merged_SM_1@meta.data %>% as.data.frame() %>% 
  rownames_to_column("spot_id") %>%
  select_at(c("spot_id", "orig.ident", "integrated_snn_res.0.5")) %>%
  dplyr::rename("Niche" = "integrated_snn_res.0.5") %>%
  mutate(Niche = paste0("Cluster_", Niche))


cellprops_info =  cell_prop  %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit()

#head(cellprops_info)
cell_props_summary_CT_pat <- cellprops_info %>%
  group_by(name, Niche) %>%
  summarize(median_CT = median(value))
head(cell_props_summary_CT_pat)

# Check niches that are unique to some patients
cell_prop <-  cell_prop %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit() %>%
  dplyr::select(-spot_id) %>%
  group_by(name) %>%
  nest() %>%
  mutate(wres = map(data, function(dat) {
    
    niches <- dat$Niche %>%
      unique() %>%
      set_names()
    
    map(niches, function(g) {
      
      test_data <- dat %>%
        mutate(test_group = ifelse(.data[["Niche"]] == g,
                                   "target", "rest")) %>%
        mutate(test_group = factor(test_group,
                                   levels = c("target", "rest")))
      
      wilcox.test(value ~ test_group, 
                  data = test_data,
                  alternative = "greater") %>%
        broom::tidy()
    }) %>% enframe("Niche") %>%
      unnest()
    
  }))

wilcox_types <- cell_prop %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
  dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
  dplyr::mutate(sign = ifelse(adj_pval < 0.005, "*", ""))

#head(wilcox_types)
ct_median_desc <- cell_prop %>%
  dplyr::select(data) %>%
  unnest() %>%
  group_by(name, Niche) %>%
  summarise(median_prop = median(value)) %>%
  mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop))
dim(ct_median_desc)

niche_car_df <- left_join(wilcox_types, ct_median_desc) %>% na.omit()
dim(niche_car_df)

ct_median_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, scaled_median_prop) %>%
  pivot_wider(names_from = Niche, values_from = scaled_median_prop) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_sign_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, sign) %>%
  pivot_wider(names_from = Niche, values_from = sign) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_median_desc_mat
ct_sign_desc_mat


library("ComplexHeatmap")
Heatmap(ct_median_desc_mat, name = "scaled comp", 
        rect_gp = gpar(col = "black", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf(ct_sign_desc_mat[i, j]), x, y, gp = gpar(fontsize = 14))
        })


####Fig1F
###annotion
cell_type = c("Immune/Fibro", "Fibro", "Fibro", "Immune", "Epi", "Epi", "Epi", "Fibro", "Myelo", "Epi", "Endo", "Fibro", "Fibro")
df_celltype = data.frame(Oringe_cluster = 1:13, type = cell_type)
merged_SM_1@meta.data$cellType = df_celltype$type[match(merged_SM_1@meta.data$integrated_snn_res.0.5, df_celltype$Oringe_cluster)]
Idents(merged_SM_1) = "cellType"

col<-c("Immune/Fibro"="#7fc97f","Fibro"="#beaed4","Immune"="#e0abcf", "Endo" ="#fdd700","Epi"="#a34e3b","Myelo"="#cfe6ba")

#H&E + cellType
h1<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = "N1")+NoLegend()
c1<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.4,images ="N1",image.alpha = 0,cols = col)
h1+c1


h2<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'N5')+NoLegend()
c2<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.6,images = 'N5',image.alpha = 0,cols = col)
h2+c2

h3<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM21')+NoLegend()
c3<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.6,images = 'SM21',image.alpha = 0,cols = col)+NoLegend()
h3+c3

h4<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM25')+NoLegend()
c4<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.4,images = 'SM25',image.alpha = 0,cols = col)+NoLegend()
h4+c4

h5<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM27')+NoLegend()
c5<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.5,images = 'SM27',image.alpha = 0,cols = col)
h5+c5

h6<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM32')+NoLegend()
c6<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.6,images = 'SM32',image.alpha = 0,cols = col)+NoLegend()
h6+c6

h7<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM34')+NoLegend()
c7<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.5,images = 'SM34',image.alpha = 0,cols = col)+NoLegend()
h7+c7

h8<-SpatialPlot(merged_SM_1,  pt.size.factor =0,images = 'SM35')+NoLegend()
c8<-SpatialPlot(merged_SM_1, group.by = 'cellType',stroke = 0.15, pt.size.factor =1.6,images = 'SM35',image.alpha = 0,cols = col)
h8+c8
