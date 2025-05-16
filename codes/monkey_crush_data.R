#Load packages
library(tidyverse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)
library(openxlsx)
library(scDblFinder)
library(BiocParallel)
library(gridExtra)
library(ggpubr)
source("utils/utilFxns.R")
source("utils/plottingFxns.R")

## Preprocessing
# Load files and generate full object
mk_ctrl1 <- Read10X("Exp_Mtx/monkey os 7dpi snrna-seq ON&chiasm/cellranger/M2OOD/filtered_feature_bc_matrix/")
mk_7dpi1 <- Read10X("Exp_Mtx/monkey os 7dpi snrna-seq ON&chiasm/cellranger/M1OOS/filtered_feature_bc_matrix/")
mk_7dpi2 <- Read10X("Exp_Mtx/monkey os 7dpi snrna-seq ON&chiasm/cellranger/M2OOS/filtered_feature_bc_matrix/")

# add animal tags
colnames(mk_ctrl1) <- paste("ctrl1", colnames(mk_ctrl1), sep = "_")
colnames(mk_7dpi1) <- paste("7dpi1", colnames(mk_7dpi1), sep = "_")
colnames(mk_7dpi2) <- paste("7dpi2", colnames(mk_7dpi2), sep = "_")
on_mat <- cbind(mk_ctrl1, mk_7dpi1, mk_7dpi2)
mk_merge <- CreateSeuratObject(on_mat)

# add animal tag to metadata
mk_merge@meta.data[colnames(mk_ctrl1), 'injury_group'] = "control"
mk_merge@meta.data[colnames(mk_7dpi1), 'injury_group'] = "7dpi"
mk_merge@meta.data[colnames(mk_7dpi2), 'injury_group'] = "7dpi"
rm(on_mat, mk_ctrl1, mk_ctrl2, mk_7dpi1, mk_7dpi2);gc()
mk_merge@meta.data$injury_group <- factor(mk_merge@meta.data$injury_group, levels = c("control", "7dpi"))

# initial cluster
mk_merge <- ClusterSeurat(mk_merge,use.SCT=F)

# drop doublet 
library("scDblFinder")
library("BiocParallel")
sce <- as.SingleCellExperiment(mk_merge) 
bp <- MulticoreParam(3, RNGseed=1234)
bpstart(bp) 
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3)) 
bpstop(bp)
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE) 
table(sce$scDblFinder.class) 
mk_merge_sce <- as.Seurat(sce)
DimPlot(mk_merge_sce,group.by = "scDblFinder.class")
mk_merge <- mk_merge_sce
rm(mk_merge_sce,sce,bp);gc()
mk_merge <- subset(mk_merge, scDblFinder.class == "singlet")

# filter cells
DefaultAssay(mk_merge) <- "RNA"
mk_merge <- subset(mk_merge, subset = nCount_RNA < 8000 & nFeature_RNA >200)  

# re-cluster
mk_merge <- ClusterSeurat(mk_merge,use.SCT = F)

# find DEG
dir.create("monkey_crush_data/Markers/", recursive = T)
Idents(mk_merge) <- "seurat_clusters"
all_markers <- FindAllMarkers(mk_merge, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,assay = "RNA")
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "monkey_crush_data/Markers/mk_merge_res0.5_Markers.rds")
openxlsx::write.xlsx(all_markers,file = "monkey_crush_data/Markers/mk_merge_res0.5_Markers.xlsx",rowNames=T)

# annotation for clusters
Idents(mk_merge) <- "seurat_clusters"
mk_merge@meta.data$cell_class <- "Other"

mk_merge@meta.data[WhichCells(mk_merge, idents = c(6,13,21)),]$cell_class = "Astrocyte"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(0,1,12,19)),]$cell_class = "Oligodendrocyte"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(2,5,9,11,14,15)),]$cell_class = "MP"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(3,4,7,16,17)),]$cell_class = "Fibroblast"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(8)),]$cell_class = "OPC"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(10)),]$cell_class = "VEC"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(18)),]$cell_class = "Mural cell"
mk_merge@meta.data[WhichCells(mk_merge, idents = c(20)),]$cell_class = "T cell"

# further analysis of MP subsets
Idents(mk_merge) <- "cell_class"
mk_mp <- subset(mk_merge, ident ="MP")

# run clustering pipeline
DefaultAssay(mk_mp) <- "RNA"
mk_mp <- ClusterSeurat(mk_mp,cluster_resolution = 0.2) # res 0.2

# Find DE markers
mk_mp_markers <- FindAllMarkers(mk_mp,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
mk_mp_markers <- mk_mp_markers %>% arrange(desc(avg_log2FC))
saveRDS(mk_mp_markers, file = "monkey_crush_data/Markers/mk_mp_res0.2_Markers.rds")
openxlsx::write.xlsx(mk_mp_markers, file = "monkey_crush_data/Markers/mk_mp_res0.2_Markers.xlsx",rowNames=F)

# delete heterogeneous cell populations (cluster1 contains more hetero-markers and has a lower ncount, removing low quality cells)
mk_mp <- subset(mk_mp, ident="1",invert=T)

# run clustering pipeline
DefaultAssay(mk_mp) <- "RNA"
mk_mp <- ClusterSeurat(mk_mp,cluster_resolution = 0.2) # res 0.2

# Find DE markers
mk_mp_markers <- FindAllMarkers(mk_mp,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
mk_mp_markers <- mk_mp_markers %>% arrange(desc(avg_log2FC))
dir.create("monkey_crush_data/subtype/mp/markers/")
saveRDS(mk_mp_markers, file = "monkey_crush_data/Markers/mk_mp_res0.2_Markers.rds")
openxlsx::write.xlsx(mk_mp_markers, file = "monkey_crush_data/Markers/mk_mp_res0.2_Markers.xlsx",rowNames=F)

# delete heterogeneous cell populations
mk_mp <- subset(mk_mp, ident="7",invert=T)

# annotation
Idents(mk_mp) <- "seurat_clusters"
mk_mp$cell_subclass <- "others"
mk_mp@meta.data[WhichCells(mk_mp, idents = c(0,3,6)),]$cell_subclass = "Microglia"
mk_mp@meta.data[WhichCells(mk_mp, idents = c(1,2,4,5)),]$cell_subclass = "Macrophage"

# Re-import the metadata back into multi_model by exporting the metadata and merging it in excel
meta_mp <- mk_mp@meta.data
meta_all <- mk_merge@meta.data
write.xlsx(meta_mp,file = "monkey_crush_data/meta_mp.xlsx",rowNames=T)
write.xlsx(meta_all,file = "monkey_crush_data/meta_all.xlsx",rowNames=T)
# Processed through excel and imported back
meta_all <- read.xlsx("monkey_crush_data/meta_all.xlsx")
mk_merge$cell_subclass <- meta_all$cell_subclass

# plot for cell_subclass
Idents(mk_merge) <- "cell_subclass"
seurat_cols <- c("#479D88","#3C77AF","#D1352B","#BBDD78","#EE934E","#FFBC40","#B383B9","#F5CFE4","#F6948D")
DimPlot(mk_merge, group.by = "cell_subclass",label = T,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "monkey_crush_data/plots/umapplot/umapplot_cellclass.pdf",width = 6.5,height = 5)
DimPlot(mk_merge, group.by = "cell_subclass",label = T, split.by = "injury_group",ncol = 2,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "monkey_crush_data/plots/umapplot/umapplot_cellclass splitby_injury_group.pdf",width = 8,height = 4)
# mk_merge$orig.ident <- factor(mk_merge$orig.ident, levels = c("ctrl1","ctrl2","ctrl3","ctrl4","7dpi1","7dpi2","1mpi1","1mpi2","3mpi1","3mpi2"))
DimPlot(mk_merge, group.by = "cell_subclass",label = T, split.by = "orig.ident", ncol = 2,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "monkey_crush_data/plots/umapplot/umapplot_cellclass splitby_orig.ident.pdf",width = 9,height = 8)

# Find celltype DE markers
Idents(mk_mp) <- "cell_subclass"
mk_merge_markers <- FindAllMarkers(mk_merge,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
mk_merge_markers <- mk_merge_markers %>% arrange(desc(avg_log2FC))
saveRDS(mk_merge_markers, file = "monkey_crush_data/Markers/mk_mp_cell_subclass_Markers.rds")
openxlsx::write.xlsx(mk_merge_markers, file = "monkey_crush_data/Markers/mk_mp_cell_subclass_Markers.xlsx",rowNames=F)

# celltype marker dotplot
olig_markers= c("MBP","MOG")
astro_markers=c("GFAP","AQP4")
fibro_markers=c("FN1","COL1A2")
opc_markers=c("PDGFRA","BRINP3")
vec_markers=c("FLT1","ERG")
mural_markers=c("ACTA2","CALD1")
micro_markers=c("P2RY12", "SRGAP2")
macro_marker=c("CD36","MSR1")
tcell_marker=c("BCL11B","SKAP1")
dir.create("monkey_crush_data/plots/dotplot/")
DotPlot(mk_merge, features = c(astro_markers,fibro_markers,macro_marker,micro_markers,mural_markers,olig_markers,opc_markers,tcell_marker,vec_markers), group.by = "cell_subclass", assay = "RNA",) + RotatedAxis()
ggsave(filename = "monkey_crush_data/plots/dotplot/cell_class_marker_dotplot_coord_flip.pdf",height = 4.5,width = 8)

# Cell Ratio Heat Map
# set Idents
Idents(mk_merge) <- "cell_subclass"

meta_df <- mk_merge@meta.data

prop_df <- meta_df %>%
  group_by(orig.ident, cell_subclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

prop_mat <- dcast(prop_df, cell_subclass ~ orig.ident, value.var = "proportion", fill = 0)
rownames(prop_mat) <- prop_mat$cell_subclass
prop_mat <- as.matrix(prop_mat[, -1])  # 去掉第一列 cell_subclass
prop_mat <- prop_mat[rev(1:nrow(prop_mat)),] # 行进行数据颠倒，与dotplot进行匹配
prop_mat <- prop_mat[,c(1,3,2)] # 行进行数据颠倒，与dotplot进行匹配

col_fun_prop <- circlize::colorRamp2(c(0, max(prop_mat)), c("grey90", "blue"))

pdf(file = "monkey_crush_data/plots/heatmap/cellsubclass_proportion_heatmap_splitby_orig.ident.pdf", 
    width = 2.5, height = 5)

Heatmap(
  prop_mat,
  name = "Proportion",
  col = col_fun_prop,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = F,
  show_column_names = TRUE,
  column_names_rot = 90,
  border = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.5)
)

dev.off()

# log2fc heat map
# Define your marker gene and its corresponding type
marker_dict <- list(
  Astrocyte = c("GFAP", "AQP4"),
  Fibroblast = c("FN1", "COL1A2"),
  Macrophage = c("CD36", "MSR1"),
  Microglia = c("P2RY12", "SRGAP2"),
  `Mural cell` = c("ACTA2", "CALD1"),
  Oligodendrocyte = c("MBP", "MOG"),
  OPC = c("PDGFRA", "BRINP3"),
  `T cell` = c("BCL11B", "SKAP1"),
  VEC = c("FLT1", "ERG")
)

gene_to_type <- stack(marker_dict)
colnames(gene_to_type) <- c("gene", "celltype")
features <- gene_to_type$gene
types <- gene_to_type$celltype

samples <- unique(mk_merge$orig.ident)

pct_list <- list()
logfc_list <- list()

for (s in samples) {
  message("Processing sample: ", s)
  
  sample_obj <- subset(mk_merge, subset = orig.ident == s)
  Idents(sample_obj) <- "cell_subclass"
  
  deg <- FindAllMarkers(
    object = sample_obj,
    only.pos = FALSE,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1
  )
  
  deg_filtered <- merge(deg, gene_to_type,
                        by.x = c("gene", "cluster"),
                        by.y = c("gene", "celltype"))
  
  logfc_vec <- setNames(deg_filtered$avg_log2FC, deg_filtered$gene)
  pct1_vec  <- setNames(deg_filtered$pct.1, deg_filtered$gene)
  
  logfc_list[[s]] <- setNames(rep(0, length(features)), features)
  logfc_list[[s]][names(logfc_vec)] <- logfc_vec
  
  pct_list[[s]] <- setNames(rep(0, length(features)), features)
  pct_list[[s]][names(pct1_vec)] <- pct1_vec
}

logfc_mat <- do.call(cbind, logfc_list)
pct_mat   <- do.call(cbind, pct_list)

logfc_mat <- logfc_mat[rev(features), ]
pct_mat   <- pct_mat[rev(features), ]

col_fun_pct <- circlize::colorRamp2(c(0, max(logfc_mat)), c("grey90", "blue"))

dir.create("monkey_crush_data/plots/heatmap/")
pdf(file = "monkey_crush_data/plots/heatmap/marker_heatmap_splitby_orig.ident.pdf",width = 2.5,height = 5)
Heatmap(
  logfc_mat,
  name = "log2FC",
  col = col_fun_pct,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = F,
  show_column_names = TRUE,
  column_names_rot = 90,
  border = FALSE,  # 关闭整体边框
  rect_gp = gpar(col = "black", lwd = 0.5))  # 每个cell有黑框线，线宽可调
dev.off()

# paired_barplot
counts <- table(mk_merge@meta.data$cell_subclass, mk_merge@meta.data$injury_group)
counts <- as.data.frame(counts)
colnames(counts) <- c("cell_type", "injury_group", "Freq")

# 添加百分比信息（按模型组归一化）
counts$Percent <- counts$Freq / ave(counts$Freq, counts$injury_group, FUN = sum)

# 只保留目标细胞
target_cells <- c("Macrophage", "Microglia", "T cell")
counts <- counts %>% filter(cell_type %in% target_cells)

target_models <- c("control", "7dpi")

model_data <- counts %>%
  mutate(injury_group = factor(injury_group, levels = c("control","7dpi")))

p <- ggplot(model_data, aes(x = cell_type, y = Percent * 100, fill = injury_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("control" = "#1f77b4", "7dpi" = "#ff7f0e")) +
  labs(
    title = paste("Control vs 7dpi"),
    x = "", y = "Percentage (%)"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold"),
    legend.position = "top"
  )
p

ggsave(filename = "monkey_crush_data/plots/proportion plot/paired_barplot_immune_cell.pdf",width = 3, height = 5.5)

# extract macrophage for scoring
# skull or blood score
skull.genes <- c("IGKC","GM36486","PIM1","GM10076","MFGE8","EAR6","LARS2","LRP1","MT-ND3",
                 "FN1","MT-ND4L","SLC6A6","ATP13A2","C77080","MAN2B1","STAB1","PLD3","GPC1","MOBP","CTSD")
blood.genes <- c("NGP","S100A9","NRGN","CSF2","IFNG","TNFSF11","NKG7","RGS16","IL17A","LCK","GIMAP5",
                 "GZMA","CD8B1","CD4","PTPRCAP","CD3E","LAT","CTSW","TRAC","CD3G")
Idents(mk_merge) <- "cell_subclass"
mk_macro <- subset(mk_merge,ident="Macrophage")
mk_macro <- AddModuleScore(mk_macro,features = list(skull.genes) ,name = "EAE_Skull_derived_score" )
mk_macro <- AddModuleScore(mk_macro,features = list(blood.genes) ,name = "EAE_Blood_derived_score" )
mk_macro@meta.data$Skull_derived_score <- rescale(mk_macro@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
mk_macro@meta.data$Blood_derived_score <- rescale(mk_macro@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
pal <- viridis::viridis(n = 10, option = "H", direction = 1)

FeaturePlot(mk_macro,features = "Skull_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "monkey_crush_data/plots/geneset score plot/EAE_Skull_derived_score_order_true_macrophage.pdf",width = 6,height = 5.6)
FeaturePlot(mk_macro,features = "Blood_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "monkey_crush_data/plots/geneset score plot/EAE_Blood_derived_score_order_true_macrophage.pdf",width = 6,height = 5.6)
# vlnplot
meta_df <- mk_macro@meta.data
macrophage_df <- meta_df %>%
  filter(cell_subclass == "Macrophage") %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, injury_group)
macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")
p<-ggplot(macrophage_long, aes(x = Score_Type, y = Score, fill = Score_Type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",        # 显示“*”等符号
    comparisons = list(c("Skull_derived_score", "Blood_derived_score")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,              # 控制*符号显示位置
    size=10
  ) +
  scale_fill_manual(values = c("#fcae59", "#db498e")) +
  labs(x = "", y = "Score") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")+
  theme_classic2(base_size = 14) +  # 改变主题样式
  theme(legend.position = "none")+   # 最后关闭 legend
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20, face = "bold"),  # 横轴标题
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵轴标题
    axis.text.x = element_text(size = 16, face = "bold",angle = 25,vjust = 0.85,hjust = 0.75),                  # 横轴刻度
    axis.text.y = element_text(size = 16, face = "bold")                   # 纵轴刻度
  )+
  ylim(c(0,14))
ggsave(filename = "monkey_crush_data/plots/geneset score plot/Skull_vs_blood_score_vlnplot_macrophage.pdf",height = 6,width = 4.5)

# Extracting t cells and macrophage for scoring
# skull or blood score
Idents(mk_merge) <- "cell_subclass"
mk_t <- subset(mk_merge,ident=c("T cell","Macrophage"))
mk_t <- AddModuleScore(mk_t,features = list(skull.genes) ,name = "EAE_Skull_derived_score" )
mk_t <- AddModuleScore(mk_t,features = list(blood.genes) ,name = "EAE_Blood_derived_score" )
mk_t@meta.data$Skull_derived_score <- rescale(mk_t@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
mk_t@meta.data$Blood_derived_score <- rescale(mk_t@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
FeaturePlot(mk_t,features = "Skull_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "monkey_crush_data/plots/geneset score plot/EAE_Skull_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)
FeaturePlot(mk_t,features = "Blood_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "monkey_crush_data/plots/geneset score plot/EAE_Blood_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)
# vlnplot
meta_df <- mk_t@meta.data
t_df <- meta_df %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, injury_group,cell_subclass)
t_long <- t_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")
library(ggpubr)
p <- ggplot(t_long, aes(x = cell_subclass, y = Score, fill = cell_subclass)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",        # 显示“*”等符号
    comparisons = list(c("Macrophage", "T cell")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,              # 控制*符号显示位置
    size=10
  ) +
  scale_fill_manual(values = c("#fcae59", "#db498e")) +
  labs(x = "", y = "Score") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")+
  theme_classic2(base_size = 14) +  # 改变主题样式
  theme(legend.position = "none")+   # 最后关闭 legend
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20, face = "bold"),  # 横轴标题
    axis.title.y = element_text(size = 20, face = "bold"),  # 纵轴标题
    axis.text.x = element_text(size = 16, face = "bold",angle = 25,vjust = 0.85,hjust = 0.75),                  # 横轴刻度
    axis.text.y = element_text(size = 16, face = "bold")                   # 纵轴刻度
  )+
  facet_wrap(~Score_Type, ncol = 1) +  # 按样本分面，ncol 可改
  ylim(c(0,14))
ggsave(filename = "monkey_crush_data/plots/geneset score plot/Skull_vs_blood_score_vlnplot_tcell&macrophage.pdf",height = 8,width = 4)


