library(tidyverse)
library(Seurat)
library(ggplot2)
library(spacexr)
library(Matrix)
library(ComplexHeatmap)
library(monocle)
library(RColorBrewer)
library(patchwork)

#### Fig6a
##hallmarkers

CSCC<-readRDS("/media/desk16/tly0202/TJ/all_sample_SME_seurat.rds.gz")
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


one_sample = CSCC[[10]]
SM35_all_hallmark<-hallmarks_score(one_sample)


SM35 = CSCC[[10]]
tumor_SM35 = subset(SM35,Morph_snn_res.1.5 %in% c(5,6,13,9))
df_layer = data.frame(Oringe_cluster = c(13, 5, 6, 9), Layer = c("Layer1", "Layer2", "Layer2", "Layer3"))
tumor_SM35@meta.data$layer = df_layer$Layer[match(tumor_SM35@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
SM35_sub_hallmark = SM35_all_hallmark[rownames(tumor_SM35@meta.data),]
nrow(tumor_SM35@meta.data)
nrow(SM35_sub_hallmark)
tumor_SM35@meta.data = cbind(tumor_SM35@meta.data, SM35_sub_hallmark)

hall_path = colnames(SM35_all_hallmark)

all_mean = map(hall_path, function(i){
    mean_path = tumor_SM35@meta.data %>% select(layer, i) %>% group_by(layer) %>% summarise_at(vars(i),list(name = mean))
    mean_path$Pathway = i
    return(mean_path)
}) %>% bind_rows()

head(all_mean)
hall_path

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

scaled_data <- all_mean %>% group_by(Pathway) %>% mutate(scale_ind = scale_this(name))
head(scaled_data)

Immune = c("ALLOGRAFT_REJECTION",
                    "COAGULATION",
                    "COMPLEMENT",
                    "IL6_JAK_STAT3_SIGNALING",
                    "IL2_STAT5_SIGNALING",
                    "INFLAMMATORY_RESPONSE",
                    "TNFA_SIGNALING_VIA_NFKB",
                    "INTERFERON_ALPHA_RESPONSE",
                    "INTERFERON_GAMMA_RESPONSE"
            )
Metabolism = c(
                        "CHOLESTEROL_HOMEOSTASIS",
                        "GLYCOLYSIS",
                        "HEME_METABOLISM",
                        "BILE_ACID_METABOLISM",
                        "FATTY_ACID_METABOLISM",
                        "OXIDATIVE_PHOSPHORYLATION",
                        "XENOBIOTIC_METABOLISM")
Signaling = c("ANDROGEN_RESPONSE",
              "ESTROGEN_RESPONSE_EARLY",
              "ESTROGEN_RESPONSE_LATE",
              "KRAS_SIGNALING_UP",
              "KRAS_SIGNALING_DN",
              "MTORC1_SIGNALING",
              "NOTCH_SIGNALING",
              "PI3K_AKT_MTOR_SIGNALING",
              "HEDGEHOG_SIGNALING",
              "TGF_BETA_SIGNALING",
              "WNT_BETA_CATENIN_SIGNALING",
              "APOPTOSIS",
              "HYPOXIA",
              "PROTEIN_SECRETION",
              "UNFOLDED_PROTEIN_RESPONSE",
              "REACTIVE_OXYGEN_SPECIES_PATHWAY"
              )

Proliferation = c("E2F_TARGETS",
                  "G2M_CHECKPOINT",
                  "MYC_TARGETS_V1",
                  "MYC_TARGETS_V2",
                  "P53_PATHWAY",
                  "MITOTIC_SPINDLE")

scaled_data = scaled_data[scaled_data$Pathway %in% c(Immune, Metabolism, Signaling, Proliferation),]
scaled_data$Pathway = factor(scaled_data$Pathway, levels = c(Immune, Metabolism, Signaling, Proliferation))

p = ggplot(scaled_data, aes(layer, Pathway, fill= scale_ind)) + 
  geom_tile() +
  labs(title = "SM35") +
  scale_fill_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24")+
  coord_fixed(ratio = 1/3)
  
##SM27
one_sample = CSCC[[5]]
SM27_all_hallmark<-hallmarks_score(one_sample)


SM27 = CSCC[[5]]
tumor_SM27 = subset(SM27,Morph_snn_res.1.5 %in% c(2, 9, 0))
df_layer = data.frame(Oringe_cluster = c(2, 9, 0), Layer = c("Layer1", "Layer2", "Layer3"))
tumor_SM27@meta.data$layer = df_layer$Layer[match(tumor_SM27@meta.data$seurat_clusters, df_layer$Oringe_cluster)]
SM27_sub_hallmark = SM27_all_hallmark[rownames(tumor_SM27@meta.data),]
nrow(tumor_SM27@meta.data)
nrow(SM27_sub_hallmark)
tumor_SM27@meta.data = cbind(tumor_SM27@meta.data, SM27_sub_hallmark)

hall_path = colnames(SM27_all_hallmark)

all_mean_27 = map(hall_path, function(i){
    mean_path = tumor_SM27@meta.data %>% select(layer, i) %>% group_by(layer) %>% summarise_at(vars(i),list(name = mean))
    mean_path$Pathway = i
    return(mean_path)
}) %>% bind_rows()

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

scaled_data <- all_mean_27 %>% group_by(Pathway) %>% mutate(scale_ind = scale_this(name))
head(scaled_data)

Immune = c("ALLOGRAFT_REJECTION",
                    "COAGULATION",
                    "COMPLEMENT",
                    "IL6_JAK_STAT3_SIGNALING",
                    "IL2_STAT5_SIGNALING",
                    "INFLAMMATORY_RESPONSE",
                    "TNFA_SIGNALING_VIA_NFKB",
                    "INTERFERON_ALPHA_RESPONSE",
                    "INTERFERON_GAMMA_RESPONSE"
            )
Metabolism = c(
                        "CHOLESTEROL_HOMEOSTASIS",
                        "GLYCOLYSIS",
                        "HEME_METABOLISM",
                        "BILE_ACID_METABOLISM",
                        "FATTY_ACID_METABOLISM",
                        "OXIDATIVE_PHOSPHORYLATION",
                        "XENOBIOTIC_METABOLISM")
Signaling = c("ANDROGEN_RESPONSE",
              "ESTROGEN_RESPONSE_EARLY",
              "ESTROGEN_RESPONSE_LATE",
              "KRAS_SIGNALING_UP",
              "KRAS_SIGNALING_DN",
              "MTORC1_SIGNALING",
              "NOTCH_SIGNALING",
              "PI3K_AKT_MTOR_SIGNALING",
              "HEDGEHOG_SIGNALING",
              "TGF_BETA_SIGNALING",
              "WNT_BETA_CATENIN_SIGNALING",
              "APOPTOSIS",
              "HYPOXIA",
              "PROTEIN_SECRETION",
              "UNFOLDED_PROTEIN_RESPONSE",
              "REACTIVE_OXYGEN_SPECIES_PATHWAY"
              )

Proliferation = c("E2F_TARGETS",
                  "G2M_CHECKPOINT",
                  "MYC_TARGETS_V1",
                  "MYC_TARGETS_V2",
                  "P53_PATHWAY",
                  "MITOTIC_SPINDLE")

scaled_data = scaled_data[scaled_data$Pathway %in% c(Immune, Metabolism, Signaling, Proliferation),]
scaled_data$Pathway = factor(scaled_data$Pathway, levels = c(Immune, Metabolism, Signaling, Proliferation))

p = ggplot(scaled_data, aes(layer, Pathway, fill= scale_ind)) + 
  geom_tile() +
  labs(title = "SM27") +
  scale_fill_gradient2(low = "#106ca7",
                       mid = "#dbddd8",
                       high = "#bb1b24")+
  coord_fixed(ratio = 1/3)


##### Fig6c

library(tidyverse)
library(Cardinal)
library(ggplot2)
library(ggthemes)


#############################################step1 select the region############################################
##read the region
focal_file = list.files('/home/zhluo/Project/CESC/selected_resions/SM27_neg/', recursive = T, full.names = T, 
                            pattern = 'txt$')
selected_posi = lapply(focal_file, function(x){
  focal_1 = read.csv(x, sep="\t", header = F)
  colnames(focal_1) = c("Num", "X", "Y", "State")
  focal_1 = focal_1[focal_1$State == 1,]
  #focal_1$pos = paste0(focal_1$X, "x", focal_1$Y)
  return(focal_1)
  }) %>% do.call(rbind, .)

##subset
focal_1_data = tiny_in[,pixels(tiny_in, coord=list(x=selected_posi$X, y=selected_posi$Y))]
##plot
image(focal_1_data, mz=333.2788, zlim=c(0,100), colorscale=jet.colors, smooth.image="gaussian")

###########################################step2 outline the boundary###########################################
##calculating distance
focal_region = selected_posi
##plot data
ymax = 89  #pos 89 neg 94
xmax = 115 #pos 115 neg 121
focal_region$Y = ymax - focal_region$Y
focal_region = focal_region[,c("X", "Y")]
colnames(focal_region) = c("col", "row")
#校正位置坐标
focal_region$cor_name = paste0(focal_region$col, "x", focal_region$row)

##section
get_circle_n_spot = function(coor_df, spot){
  
  #spot = as.vector(spot)
  spot_top = c(spot["row"] + 1, spot["col"])
  spot_top_left = c(spot["row"] + 1, spot["col"] - 1)
  spot_top_right = c(spot["row"] + 1, spot["col"] + 1)
  spot_right = c(spot["row"], spot["col"] + 1)
  spot_bottom_right = c(spot["row"] - 1, spot["col"] + 1)
  spot_bottom_left = c(spot["row"] - 1, spot["col"] - 1)
  spot_bottom = c(spot["row"] - 1, spot["col"])
  spot_left = c(spot["row"], spot["col"] - 1)
  df_circle = data.frame(row=c(spot_top[1], spot_top_left[1], spot_top_right[1],spot_right[1],spot_bottom_right[1],spot_bottom_left[1],spot_left[1], spot_bottom[1]),
                         col=c(spot_top[2], spot_top_left[2], spot_top_right[2],spot_right[2],spot_bottom_right[2],spot_bottom_left[2],spot_left[2], spot_bottom[2])
  )
  
  df_circle$cor_name = paste0(df_circle$col, "x", df_circle$row)
  df_circle = coor_df[coor_df$cor_name %in% df_circle$cor_name,]
  return(df_circle)
}

for (idx in rownames(focal_region)){
  spot = c("row"=focal_region$row[rownames(focal_region)==idx], "col"=focal_region$col[rownames(focal_region)==idx])
  suround = get_circle_n_spot(focal_region, spot)
  if(nrow(suround) < 8) {
    focal_region[rownames(focal_region)==idx, "group"] = "Boundary"
  }
  else{
    focal_region[rownames(focal_region)==idx, "group"] = "Inner"
  }
}



p = ggplot(focal_region, aes(x=col, y=row, colour =group)) + geom_point() + theme_base() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) +
  labs(x="", y="", title = "Selected region")
p


###########################################step3 the distance###########################################
##calculating distance
outer_spot = focal_region[focal_region$group == "Boundary",]

calSpotSpotDistance <- function(spot1, spot2)
{
  spot1_xy <- as.numeric(unlist(strsplit(spot1,"x")))
  spot2_xy <- as.numeric(unlist(strsplit(spot2,"x")))
  sqrt((spot1_xy[1]-spot2_xy[1])^2+(spot1_xy[2]-spot2_xy[2])^2)
}

calSpotVecDistance <- function(spot1, vec2)
{
  min(sapply(vec2,calSpotSpotDistance,spot1=spot1))
}

calVecVecDistance <- function(vec1,vec2)  # vec1 M2CAF, vec2 border
{
  sapply(vec1,calSpotVecDistance,vec2=vec2)
}

focal_region$distance = calVecVecDistance(focal_region$cor_name, outer_spot$cor_name)

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
p = ggplot(focal_region, aes(x=col, y=row, colour =distance)) + geom_point() + theme_base() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) +
  labs(x="", y="", title = "SM35 boundary distance") + 
  scale_color_gradientn(colors = jet.colors(7))

p

###########################################step3 calculate correlation###########################################
##MZ intensity
iData(focal_1_data, "intensity")
iData(focal_1_data, "intensity")
featureData(focal_1_data)
pixelData(focal_1_data)

image(focal_1_data, mz=333.2788, zlim=c(0,100), colorscale=jet.colors, smooth.image="gaussian")

pixelData(focal_1_data)
fid <- features(focal_1_data, 70 < mz & mz < 963) #neg 914 pos 963

iData(focal_1_data, "log2intensity") <- log2(iData(focal_1_data) + 1)
intensity_df = iData(focal_1_data, "log2intensity")




focal_region_1 = focal_region
#focal_region_1$row = ymax - focal_region_1$row
pixelData(focal_1_data)$pos = paste0(pixelData(focal_1_data)@coord$x, "x", ymax- pixelData(focal_1_data)@coord$y)

pixelData(focal_1_data)$distance = focal_region_1$distance[match(pixelData(focal_1_data)$pos, focal_region_1$cor_name)]
distance_pix = pixelData(focal_1_data)

lapply(1:length(mz(focal_1_data)), function(x){
  res = cor.test(intensity_df[x,], distance_pix$distance)
  mz_value = mz(focal_1_data)[x]
  df_res = DataFrame(mz=mz_value, cor = res$estimate, p_value = res$p.value)
  return(df_res)
}) %>% do.call(rbind,.) -> mz_result

write.csv(mz_result,file = "/home/zhluo/Project/CESC/metabolism_result/mz_correlation_pos_SM27.csv", quote = F, row.names = F)

image(focal_1_data, mz=465.3028, zlim=c(0,5000), colorscale=jet.colors, smooth.image="gaussian")


##########################################merge table##########################################################
##阴离子
SM27_neg = read.csv("/home/zhluo/Project/CESC/metabolism_result/mz_correlation_negtive_SM27.csv") %>% 
  rename(mz="Neg(m/z)", cor="Correlation_SM27", p_value="P-Value_SM27")  %>% 
  na.omit()

SM35_neg = read.csv("/home/zhluo/Project/CESC/metabolism_result/mz_correlation_negtive_SM35.csv") %>% 
  rename(mz="Neg(m/z)", cor="Correlation_SM35", p_value="P-Value_SM35")  %>% 
  na.omit()

neg_table = SM35_neg %>% left_join(SM27_neg, by=c("Neg(m/z)" = "Neg(m/z)")) %>% na.omit() %>% arrange(desc(Correlation_SM35))
plot(neg_table$Correlation_SM35, neg_table$Correlation_SM27)

neg_table$sig = ifelse(abs(neg_table$`P-Value_SM35`) <= 0.05 & abs(neg_table$`P-Value_SM27`) <= 0.05, "Sig", "Not Sig")

p2=ggplot( neg_table, aes(x = Correlation_SM35, y = Correlation_SM27, color = sig)) +
  geom_point(size=0.9)+geom_smooth(method=lm,color="black",size=0.7)+
  scale_colour_manual(values=c("#8696a7", "#CB9C7A"))+
  labs(x="SM35",y="SM27",title="Negtive ion ")+
  annotate("text", x = 0.3, y =-0.3, label = "r = 0.72, p-value <2e-16",
           color="#350E20FF",size = 5 )+
  #annotate("text", x = 0.4, y =-0.3, label = "",
  #         color="#350E20FF",size = 5 )+
  #geom_text(data=data_text, mappingaes(x=x,y=y,label=label,colour=NULL),nudge_x=0.1,nudge_y=0.1)+
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  theme(title = element_text(size = 14,
                             color = "black"),
        axis.title = element_text(size=10),
        #axis.text = element_text(size=12, color="black"),
        axis.text.x = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 14,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.x= element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.y= element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        
        legend.title = element_blank(),
        legend.text = element_text(size=14), legend.key.size = unit(0.9,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )+ guides(colour = guide_legend(override.aes = list(size=2.5)))#图例点大小
p2


##阳离子
SM27_pos = read.csv("/home/zhluo/Project/CESC/metabolism_result/mz_correlation_positive_SM27.csv") %>% 
  rename(mz="Pos(m/z)", cor="Correlation_SM27", p_value="P-Value_SM27")  %>% 
  na.omit()

SM35_pos = read.csv("/home/zhluo/Project/CESC/metabolism_result/mz_correlation_positive_SM35.csv") %>% 
  rename(mz="Pos(m/z)", cor="Correlation_SM35", p_value="P-Value_SM35")  %>% 
  na.omit()


pos_table = SM35_pos %>% left_join(SM27_pos, by=c("Pos(m/z)" = "Pos(m/z)")) %>% na.omit() %>% arrange(desc(Correlation_SM35))
plot(pos_table$Correlation_SM35, pos_table$Correlation_SM27)


pos_table$sig = ifelse(abs(pos_table$`P-Value_SM35`) <= 0.05 & abs(pos_table$`P-Value_SM27`) <= 0.05, "Sig", "Not Sig")
#write_xlsx(pos_table, "/home/zhluo/Project/CESC/metabolism_result/pos_merged_correlation.xlsx")

p2=ggplot( pos_table, aes(x = Correlation_SM35, y = Correlation_SM27, color = sig)) +
  geom_point(size=0.9)+geom_smooth(method=lm,color="black",size=0.7)+
  scale_colour_manual(values=c("#8696a7", "#CB9C7A"))+
  labs(x="SM35",y="SM27",title="Positive ion")+
  annotate("text", x = 0.3, y =-0.3, label = "r = 0.51, p-value <2e-16",
           color="#350E20FF",size = 5 )+
  #geom_text(data=data_text, mapping=aes(x=x,y=y,label=label,colour=NULL),nudge_x=0.1,nudge_y=0.1)+
  theme_minimal() +
  coord_fixed() +  
  geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + 
  theme(title = element_text(size = 14,
                             color = "black"),
        axis.title = element_text(size=10),
        #axis.text = element_text(size=12, color="black"),
        axis.text.x = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.text.y = element_text(size = 16,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.x= element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        axis.title.y= element_text(size = 20,colour = "black", margin=margin(2, 2, 2, 2, "mm")),
        
        legend.title = element_blank(),
        legend.text = element_text(size=14), legend.key.size = unit(0.9,"cm"),
        axis.ticks.length.x = unit(2,"mm"),
        axis.ticks.length.y = unit(2,"mm")
  )+ guides(colour = guide_legend(override.aes = list(size=2.5)))#图例点大小
p2

#### fig6d
options(digits=7)
library(readxl)
library(writexl)

##neg metabolite
neg_xl <- read_excel("/home/zhluo/Project/CESC/metabolism_result/Qualitative.xlsx", sheet = "neg") %>% 
  select(mz, Formula, Metabolites, Class, KEGG, 'Compound ID') %>%  mutate_at(c('mz'), as.numeric)
neg_xl$mz = round(neg_xl$mz, 4)
neg_table["Neg(m/z)"] = round(neg_table["Neg(m/z)"], 4)
neg_table = neg_table %>% left_join(neg_xl, c("Neg(m/z)" = "mz"))
write_xlsx(neg_table, "/home/zhluo/Project/CESC/metabolism_result/neg_merged_correlation.xlsx")

##pos metabolite
pos_xl <- read_excel("/home/zhluo/Project/CESC/metabolism_result/Qualitative.xlsx", sheet = "pos") %>% 
  select(mz, Formula, Metabolites, Class, KEGG, 'Compound ID') %>%  mutate_at(c('mz'), as.numeric)
pos_xl$mz = round(pos_xl$mz, 4)
pos_table["Pos(m/z)"] = round(pos_table["Pos(m/z)"], 4)
pos_table = pos_table %>% left_join(pos_xl, c("Pos(m/z)" = "mz"))
write_xlsx(pos_table, "/home/zhluo/Project/CESC/metabolism_result/pos_merged_correlation.xlsx")


##pie plot
neg_class = na.omit(neg_table$Class[neg_table$sig == "Sig" & neg_table$Correlation_SM35 > 0])
neg_summary = c()
for (one_mz in neg_class){
  re = str_split(one_mz, ";", n = 2, simplify = FALSE)
  neg_summary = c(neg_summary, re[[1]][1])
}
neg_summary = neg_summary%>% table() %>% as.data.frame()
colnames(neg_summary) = c("group", "freq")
neg_summary$group = as.character(neg_summary$group)
neg_summary_1 = neg_summary[ neg_summary$group %in% c("Fatty Acyls", "Glycerolipids", "Glycerophospholipids"),]
sum_others = sum(neg_summary$freq[! neg_summary$group %in% c("Fatty Acyls", "Glycerolipids", "Glycerophospholipids")])
neg_summary_1[nrow(neg_summary_1) + 1,] = c("Others", sum_others)
neg_summary_1$group = factor(neg_summary_1$group, levels = neg_summary_1$group)
neg_summary_1$freq = as.numeric(neg_summary_1$freq)
neg_summary_1$value = round(neg_summary_1$freq / sum(neg_summary_1$freq), digits = 2)
colnames(neg_summary_1) = c("group", "freq", "value")
#df <- data.frame(value = c(share_num/total_11, 1-share_num/total_11),
#                 group = c("Spots In niche_10", "Spots not in niche_10"))
df = neg_summary_1

df2 <- df %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))
#df2 = df2[df2$group=="Metabolite",]

library(ggrepel)
library(ggsci)
p = ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = "black") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste0(round(value,3) * 100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_brewer() +
  theme_void() + scale_fill_npg()
p
##enrichment analysis
#neg_meta = na.omit(neg_table$`Compound ID`[neg_table$sig == "Sig" & neg_table$Correlation_SM35 > 0])
#neg_meta = substr(neg_meta, 1,11)




####fig6e
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

library(scMetabolism)
library(tidyverse)
library(rsvd)
tumor_SM35 <- sc.metabolism.Seurat(obj = tumor_SM35, 
                                        method = "VISION",     
                                        imputation = F, ncores = 2, 
                                        metabolism.type = "KEGG")
metabolism.matrix <- tumor_SM35@assays$METABOLISM$scoremetabolism.matrix
metabolism.matrix
SM35_score<-as.data.frame(t(metabolism.matrix))
tumor_SM35<-AddMetaData(tumor_SM35,metadata =SM35_score)
library(ggpubr)

ggbarplot(tumor_SM35@meta.data, 
          x = "layer", 
          y = "fatty acid biosynthesis",
          fill = "layer",
          palette = "jco",
          add = "mean_se",  # 显示均值和标准误
          add.params = list(width = 0.3),
          order = c("Layer1", "Layer2", "Layer3")) +
  stat_compare_means(method = "anova", 
                     label.y = max(tumor_SM35@meta.data$fatty acid biosynthesis, na.rm = TRUE) * 1.1) +
  stat_compare_means(comparisons = list(c("Layer1", "Layer3"),
                                        c("Layer1", "Layer2"),
                                        c("Layer2", "Layer3")),
                     method = "t.test",
                     label = "p.signif") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))


####fig6g
path_in = "/home/bioinformatics/Data_download/zhaojw/Spatial_data/CESC/metabolism/imzml/SM-27-neg.imzML"
tiny_in <- readMSIData(path_in, attach.only=TRUE)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
image(tiny_in, mz=253.2166, colorscale=jet.colors, smooth.image="gaussian")#neg


####fig6h
colours = colorRampPalette(colors = rev(brewer.pal(11, name = "Spectral")))(100)
p1 = SpatialPlot(tumor_SM35, features ="CPT1A",stroke = 0.15,pt.size.factor = 1.6,images = 'SM35',image.alpha = 0)+ scale_fill_gradientn(colors = colours)
p2=Vlnplot(tumor_SM35,features="CPT1A",group.by="layer",pt.size=0)
# 添加两两t检验
p2 + stat_compare_means(
  comparisons = list(
    c("Layer1", "Layer2"),
    c("Layer1", "Layer3"),
    c("Layer2", "Layer3")
  ),
  method = "t.test",
  label = "p.signif",  # 显示显著性符号：* < 0.05, ** < 0.01, *** < 0.001
  tip.length = 0.01,
  step.increase = 0.1  # 调整线条间距
)
p3 = SpatialPlot(tumor_SM27, features ="CPT1A",stroke = 0.15,pt.size.factor = 1.6,images = 'SM35',image.alpha = 0)+ scale_fill_gradientn(colors = colours)
