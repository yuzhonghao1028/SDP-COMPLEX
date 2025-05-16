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
on_sn_ctrl1 <- Read10X("exp_mtx/mice_control_ON/filtered_feature_bc_matrix/")
on_sn_ctrl2 <- Read10X("exp_mtx/Ms_ON_Control2/filtered_feature_bc_matrix/")
on_sn_ctrl3 <- Read10X("exp_mtx/Ms_ON_sheath_Control1/filtered_feature_bc_matrix/")
on_sn_ctrl4 <- Read10X("exp_mtx/Ms_ON_sheath_Control2/filtered_feature_bc_matrix/")
on_sn_3wpi_sohu1 <- Read10X("exp_mtx/Ms_SOHU_3wpi_ON1/filtered_feature_bc_matrix/")
on_sn_1wpi_crush1 <- Read10X("exp_mtx/Ms_ON_1wpi_Crush1/filtered_feature_bc_matrix/")
on_sn_1wpi_aion1 <- Read10X("exp_mtx/Ms-ON-Aion-7dpi/filtered_feature_bc_matrix/")
on_sn_1wpi_aion2 <- Read10X("exp_mtx/Ms_ON_1wpi_AION2/filtered_feature_bc_matrix/")
on_sn_3wpi_ms <- Read10X("exp_mtx/Ms-ON-3wpi-MS1/filtered_feature_bc_matrix/")

# add animal tags
colnames(on_sn_ctrl1) <- paste("ctrl1", colnames(on_sn_ctrl1), sep = "_")
colnames(on_sn_ctrl2) <- paste("ctrl2", colnames(on_sn_ctrl2), sep = "_")
colnames(on_sn_ctrl3) <- paste("ctrl3", colnames(on_sn_ctrl3), sep = "_")
colnames(on_sn_ctrl4) <- paste("ctrl4", colnames(on_sn_ctrl4), sep = "_")
colnames(on_sn_3wpi_sohu1) <- paste("sohu.3wpi", colnames(on_sn_3wpi_sohu1), sep = "_")
colnames(on_sn_1wpi_crush1) <- paste("crush1.1wpi", colnames(on_sn_1wpi_crush1), sep = "_")
colnames(on_sn_1wpi_aion1) <- paste("aion1.1wpi", colnames(on_sn_1wpi_aion1), sep = "_")
colnames(on_sn_1wpi_aion2) <- paste("aion2.1wpi", colnames(on_sn_1wpi_aion2), sep = "_")
colnames(on_sn_3wpi_ms) <- paste("ms.3wpi", colnames(on_sn_3wpi_ms), sep = "_")

# crate Seurat object
on_mat <- cbind(on_sn_ctrl1, on_sn_ctrl2, on_sn_ctrl3, on_sn_ctrl4, on_sn_3wpi_sohu1, on_sn_1wpi_crush1, on_sn_1wpi_aion1, on_sn_1wpi_aion2,on_sn_3wpi_ms)
multi_model <- CreateSeuratObject(on_mat, project = "multi_model", min.cells = 3, min.features = 200)
table(multi_model$orig.ident)
rm(on_mat, on_sn_ctrl1, on_sn_ctrl2, on_sn_ctrl3, on_sn_ctrl4, on_sn_3wpi_sohu1, on_sn_1wpi_crush1, on_sn_1wpi_aion1, on_sn_1wpi_aion2,on_sn_3wpi_ms);gc()

# annotation for model_type
multi_model@meta.data$model_type <- "Other"
multi_model@meta.data[multi_model$orig.ident %in% c("ctrl1","ctrl2","ctrl3","ctrl4"),]$model_type = "Control"
multi_model@meta.data[multi_model$orig.ident %in% c("aion1.1wpi","aion2.1wpi"),]$model_type = "AION"
multi_model@meta.data[multi_model$orig.ident %in% c("crush1.1wpi","crush2.1wpi"),]$model_type = "Crush"
multi_model@meta.data[multi_model$orig.ident %in% c("sohu.3wpi"),]$model_type = "SOHU"
multi_model@meta.data[multi_model$orig.ident %in% c("ms.3wpi"),]$model_type = "MS"
multi_model$model_type <- factor(multi_model$model_type,levels = c("Control","Crush","SOHU","AION","MS"))
multi_model$orig.ident <- factor(multi_model$orig.ident,levels = c("ctrl1","ctrl2","ctrl3","ctrl4","crush1.1wpi","crush2.1wpi","sohu.3wpi","aion1.1wpi","aion2.1wpi","ms.3wpi"))

# ClusterSeurat
multi_model <- ClusterSeurat(multi_model,use.SCT=T,vars.to.regress = "orig.ident",use.harmony = T,harmony.by = "orig.ident")

# drop doublet 
DefaultAssay(multi_model) <- "RNA"
multi_model <- NormalizeData(multi_model, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
multi_model <- FindVariableFeatures(multi_model, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# dbr dfault
sce <- as.SingleCellExperiment(multi_model) 
bp <- MulticoreParam(3, RNGseed=1234)
bpstart(bp) 
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3)) 
bpstop(bp)
# 可以考虑使用BPPARAM参数对其进行多线程处理（假设有足够的RAM）
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE) 
# plotDoubletMap(sce) 
table(sce$scDblFinder.class) 
multi_model_sce <- as.Seurat(sce)
DimPlot(multi_model_sce,group.by = "scDblFinder.class")
multi_model <- multi_model_sce
multi_model <- subset(multi_model, scDblFinder.class == "singlet")
rm(multi_model_sce,sce,bp);gc()

# filter cells
multi_model <- subset(multi_model, subset = nCount_RNA < 7000 & nFeature_RNA >200)

# re-cluster
DefaultAssay(multi_model) <- "RNA"
multi_model <- ClusterSeurat(multi_model,use.SCT=T,vars.to.regress = "orig.ident",use.harmony = T,harmony.by = "orig.ident")

# find DEG
dir.create("multi_model_mice/Markers/", recursive = T)
Idents(multi_model) <- "seurat_clusters"
all_markers <- FindAllMarkers(multi_model, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,verbose = T,assay = "RNA")
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "multi_model_mice/Markers/multi_model_res0.5_Markers.rds")
openxlsx::write.xlsx(all_markers,file = "multi_model_mice/Markers/multi_model_res0.5_Markers.xlsx",rowNames=F)
rm(all_markers)

# delete heterogeneous cell populations
multi_model <- subset(multi_model,ident="9",invert=T)
multi_model <- subset(multi_model,ident="15",invert=T)
multi_model <- subset(multi_model,ident="18",invert=T)

# annotation for clusters
multi_model@meta.data$cell_class <- "Other"
# celltype annotation
multi_model@meta.data[WhichCells(multi_model, idents = c(0,3,8,13)),]$cell_class = "Oligodendrocyte"
multi_model@meta.data[WhichCells(multi_model, idents = c(1,2,14)),]$cell_class = "Astrocyte"
multi_model@meta.data[WhichCells(multi_model, idents = c(4,5,10,12)),]$cell_class = "Fibroblast"
multi_model@meta.data[WhichCells(multi_model, idents = c(6)),]$cell_class = "MP"
multi_model@meta.data[WhichCells(multi_model, idents = c(7)),]$cell_class = "VEC"
multi_model@meta.data[WhichCells(multi_model, idents = c(11)),]$cell_class = "OPC"
multi_model@meta.data[WhichCells(multi_model, idents = c(16)),]$cell_class = "T cell"
multi_model@meta.data[WhichCells(multi_model, idents = c(17)),]$cell_class = "Mural cell"
multi_model@meta.data[WhichCells(multi_model, idents = c(19)),]$cell_class = "Proliferating cell"

# further analysis of MP subsets
Idents(multi_model) <- "cell_class"
multi_model_mp <- subset(multi_model, idents ="MP")

# run clustering pipeline
DefaultAssay(multi_model_mp) <- "RNA"
multi_model_mp <- ClusterSeurat(multi_model_mp,use.SCT=T,vars.to.regress = "orig.ident",use.harmony = T,harmony.by = "orig.ident",cluster_resolution = 0.3) 

# Find DE markers
Idents(multi_model_mp) <- "seurat_clusters"
multi_model_mp_markers <- FindAllMarkers(multi_model_mp,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
multi_model_mp_markers <- multi_model_mp_markers %>% arrange(desc(avg_log2FC))
saveRDS(multi_model_mp_markers, file = "multi_model_mice/Markers/multi_model_mp_res0.3_Markers.rds")
openxlsx::write.xlsx(multi_model_mp_markers, file = "multi_model_mice/Markers/multi_model_mp_res0.3_Markers.xlsx",rowNames=F)

# delete heterogeneous cell populations
multi_model_mp <- subset(multi_model_mp,ident="5",invert=T)
multi_model_mp <- subset(multi_model_mp,ident="2",invert=T)

# annotation
multi_model_mp$cell_subclass <- "others"
multi_model_mp@meta.data[WhichCells(multi_model_mp, idents = c(0)),]$cell_subclass = "Microglia"
multi_model_mp@meta.data[WhichCells(multi_model_mp, idents = c(1,3,4,6,7)),]$cell_subclass = "Macrophage"
multi_model_mp@meta.data[WhichCells(multi_model_mp, idents = c(8)),]$cell_subclass = "pDC"

# Re-import the metadata back into multi_model by exporting the metadata and merging it in excel
meta_mp <- multi_model_mp@meta.data
meta_all <- multi_model@meta.data
write.xlsx(meta_mp,file = "multi_model_mice/meta_mp.xlsx",rowNames=T)
write.xlsx(meta_all,file = "multi_model_mice/meta_all.xlsx",rowNames=T)
# Processed through excel and imported back
meta_all <- read.xlsx("multi_model_mice/meta_all.xlsx")
alldata_merge$cell_subclass <- meta_all$cell_subclass

# plot for cell_subclass
Idents(multi_model) <- "cell_subclass"
seurat_cols <- c("#479D88","#3C77AF","#D1352B","#BBDD78","#EE934E","#FFBC40","#B383B9","#AECDE1","#8FA4AE","#F5CFE4","#F6948D")
DimPlot(multi_model, group.by = "cell_subclass",label = T,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "multi_model_mice/plots/umapplot/umapplot_cellclass.pdf",width = 6.5,height = 5)
DimPlot(multi_model, group.by = "cell_subclass",label = T, split.by = "model_type",ncol = 2,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "multi_model_mice/plots/umapplot/umapplot_cellclass splitby_model_type.pdf",width = 9,height = 12)
# multi_model$orig.ident <- factor(multi_model$orig.ident, levels = c("ctrl1","ctrl2","ctrl3","ctrl4","7dpi1","7dpi2","1mpi1","1mpi2","3mpi1","3mpi2"))
DimPlot(multi_model, group.by = "cell_subclass",label = T, split.by = "orig.ident", ncol = 5,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "multi_model_mice/plots/umapplot/umapplot_cellclass splitby_orig.ident.pdf",width = 20,height = 9)

# Find celltype DE markers
Idents(multi_model) <- "cell_subclass"
all_markers <- FindAllMarkers(multi_model, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,verbose = T)
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "multi_model_mice/Markers/multi_model_Markers_cell_subclass.rds")
openxlsx::write.xlsx(all_markers,file = "multi_model_mice/Markers/multi_model_Markers_cell_subclass.xlsx",rowNames=F)
rm(all_markers)

# celltype marker dotplot
tcell_markers=c("Cd3e",'Bcl11b')
olig_markers= c("Mbp","Mobp")
astro_markers=c("Gfap","Aqp4")
fibro_markers=c("Col1a2","Ddr2")
opc_markers=c("Pdgfra","Sox10")
vec_markers=c("Flt1","Erg")
micro_markers=c("Hexb", "P2ry12")
macro_markers=c("F13a1", "Dab2")
pdc_markers=c("Siglech", "Irf8")
mural_markers=c("Pdgfrb","Rgs5")
proliferating_marker <- c("Top2a",'Cenpf')
dir.create("multi_model_mice/plots/dotplot/")
DotPlot(multi_model, features = c(astro_markers,fibro_markers,macro_markers,micro_markers,mural_markers,olig_markers,opc_markers,pdc_markers,proliferating_marker,tcell_markers,vec_markers), group.by = "cell_subclass", assay = "RNA") + RotatedAxis()
ggsave(filename = "multi_model_mice/plots/dotplot/cell_class_marker_dotplot_cellsubclass_coord_flip.pdf",height = 4.5,width = 9)

# Cell Ratio Heat Map
# set Idents
Idents(multi_model) <- "cell_subclass"
# metadata
meta_df <- multi_model@meta.data
# count the number & percentage of each cell type in each orig.ident
prop_df <- meta_df %>%
  group_by(orig.ident, cell_subclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
# to matrix format: rows are cell_subclass, columns are orig.ident
library(reshape2)
prop_mat <- dcast(prop_df, cell_subclass ~ orig.ident, value.var = "proportion", fill = 0)
rownames(prop_mat) <- prop_mat$cell_subclass
prop_mat <- as.matrix(prop_mat[, -1])  # 去掉第一列 cell_subclass
prop_mat <- prop_mat[rev(1:nrow(prop_mat)),] # 行进行数据颠倒，与dotplot进行匹配
# Setting the color gradient
col_fun_prop <- circlize::colorRamp2(c(0, max(prop_mat)), c("grey90", "blue"))
# heatmap
pdf(file = "multi_model_mice/plots/heatmap/cellsubclass_proportion_heatmap_splitby_orig.ident.pdf", 
    width = 5, height = 5)
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
  Astrocyte = c("Gfap","Aqp4"),
  Fibroblast = c("Col1a2","Ddr2"),
  Macrophage = c("F13a1", "Dab2"),
  Microglia = c("Hexb", "P2ry12"),
  `Mural cell` = c("Pdgfrb","Rgs5"),
  Oligodendrocyte = c("Mbp","Mobp"),
  OPC = c("Pdgfra","Sox10"),
  pDC = c("Siglech", "Irf8"),
  `Proliferating cell` = c("Top2a",'Cenpf'),
  `T cell` = c("Cd3e",'Bcl11b'),
  VEC = c("Flt1","Erg")
)
# Collect all genes and corresponding cell types
gene_to_type <- stack(marker_dict)
colnames(gene_to_type) <- c("gene", "celltype")
features <- gene_to_type$gene
types <- gene_to_type$celltype
# Get all samples
samples <- unique(multi_model$orig.ident)
# Initialize an empty list to collect data
pct_list <- list()
logfc_list <- list()

for (s in samples) {
  message("Processing sample: ", s)
  
  # subset by sample
  sample_obj <- subset(multi_model, subset = orig.ident == s)
  Idents(sample_obj) <- "cell_subclass"
  
  # run FindAllMarkers within this sample (all cells vs. other cells)
  deg <- FindAllMarkers(
    object = sample_obj,
    only.pos = FALSE,
    logfc.threshold = 0,
    min.pct = 0,
    return.thresh = 1
  )
  
  # Merge marker_dict to extract the corresponding markers within your own cell type
  deg_filtered <- merge(deg, gene_to_type,
                        by.x = c("gene", "cluster"),
                        by.y = c("gene", "celltype"))
  
  # Extract logFC and pct.1
  logfc_vec <- setNames(deg_filtered$avg_log2FC, deg_filtered$gene)
  pct1_vec  <- setNames(deg_filtered$pct.1, deg_filtered$gene)
  
  # Align all features and fill in 0 for those that don't have any value
  logfc_list[[s]] <- setNames(rep(0, length(features)), features)
  logfc_list[[s]][names(logfc_vec)] <- logfc_vec
  
  pct_list[[s]] <- setNames(rep(0, length(features)), features)
  pct_list[[s]][names(pct1_vec)] <- pct1_vec
}

# combine into a matrix
logfc_mat <- do.call(cbind, logfc_list)
pct_mat   <- do.call(cbind, pct_list)

# Reorder rows by features
logfc_mat <- logfc_mat[rev(features), ]
pct_mat   <- pct_mat[rev(features), ]

# color gradient function
col_fun_pct <- circlize::colorRamp2(c(0, max(logfc_mat)), c("grey90", "blue"))

# heatmap
pdf(file = "multi_model_mice/plots/heatmap/marker_heatmap_splitby_orig.ident.pdf",width = 5,height = 5)
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

# celltype propotion stackbarplot
multi_model$model_type_simp = "others"
multi_model@meta.data[multi_model$model_type == "Control","model_type_simp"] = "Control"
multi_model@meta.data[multi_model$model_type == "Crush","model_type_simp"] = "Disease model"
multi_model@meta.data[multi_model$model_type == "SOHU","model_type_simp"] = "Disease model"
multi_model@meta.data[multi_model$model_type == "AION","model_type_simp"] = "Disease model"
multi_model@meta.data[multi_model$model_type == "MS","model_type_simp"] = "Disease model"

data <- as.data.frame(table(multi_model$cell_subclass,multi_model$model_type_simp))
data <- data %>% group_by(Var2) %>% mutate(percentage = Freq/sum(Freq)) %>% mutate(label = Freq/sum(Freq)*100)
data$label <- paste(sprintf("%.1f", data$label),"%")
seurat.cols <- c("#479D88","#3C77AF","#D1352B","#BBDD78","#EE934E","#FFBC40","#B383B9","#AECDE1","#8FA4AE","#F5CFE4","#F6948D")

ggplot(data, aes( x = Var2, y=percentage,fill = Var1))+
  geom_col(position = 'stack', width = 0.8)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  theme_bw()+
  scale_fill_manual(values=seurat.cols)+ 
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Percentage",
       fill="celltype",title="")+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,
                             vjust = 0.5),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
  )+ 
  guides(fill=guide_legend(keywidth = 1, keyheight = 1))+ 
  ylim(c(-0.005,1.005))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  geom_text(aes(label=label),size=4,color="black",position = "stack")

ggsave(filename = "multi_model_mice/plots/proportion plot/stack barplot_simp.pdf",width = 4.5,height = 5)

# modelwise barplot grouped by celltype
counts <- table(multi_model@meta.data$cell_subclass, multi_model@meta.data$model_type)
counts <- as.data.frame(counts)
colnames(counts) <- c("cell_type", "model_type", "Freq")
counts$Percent <- counts$Freq / ave(counts$Freq, counts$model_type, FUN = sum)
target_cells <- c("Macrophage", "Microglia", "T cell", "pDC")
counts <- counts %>% filter(cell_type %in% target_cells)
target_models <- c("Crush", "SOHU", "AION", "MS")

plots <- list()

for (model in target_models) {
  
  model_data <- counts %>%
    filter(model_type %in% c("Control", model)) %>%
    mutate(model_type = factor(model_type, levels = c("Control", model)))

  p <- ggplot(model_data, aes(x = cell_type, y = Percent * 100, fill = model_type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Control" = "#1f77b4", model = "#ff7f0e")) +
    labs(
      title = paste("Control vs", model),
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
  
  plots[[model]] <- p
}

combined_plot <- grid.arrange(grobs = plots, ncol = 4)

ggsave(
  filename = "multi_model_mice/plots/proportion plot/modelwise_barplot_grouped_by_celltype.pdf",
  plot = combined_plot,
  width = 8.2, height = 4.8
)

# extract macrophage for scoring
skull_derived_gene <- c("Igkc","Gm36486","Pim1","Gm10076","Mfge8","Ear6","Lars2","Lrp1","mt-Nd3","Fn1",
                        "mt-Nd4l","Slc6a6","Atp13a2","C77080","Man2b1","Stab1","Pld3","Gpc1","Mobp","Ctsd")
blood_derived_gene <- c("Ngp","S100a9","Nrgn","Csf2","Ifng","Tnfsf11","Nkg7","Rgs16","Il17a","Lck","Gimap5",
                        "Gzma","Cd8b1","Cd4","Ptprcap","Cd3e","Lat","Ctsw","Trac","Cd3g")

Idents(multi_model) <- "cell_subclass"
multi_model_macro <- subset(multi_model,ident="Macrophage")
multi_model_macro <- AddModuleScore(multi_model_macro,features = list(skull_derived_gene) ,name = "EAE_Skull_derived_score",assay = "RNA")
multi_model_macro <- AddModuleScore(multi_model_macro,features = list(blood_derived_gene) ,name = "EAE_Blood_derived_score",assay = "RNA")
multi_model_macro@meta.data$EAE_Skull_derived_score <- rescale(multi_model_macro@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
multi_model_macro@meta.data$EAE_Blood_derived_score <- rescale(multi_model_macro@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
pal <- viridis::viridis(n = 10, option = "H", direction = 1)

FeaturePlot(multi_model_macro,features = "EAE_Skull_derived_score",split.by = "model_type",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) & NoAxes() #&  theme(legend.position = "right")
ggsave(filename = "multi_model_mice/plots/featureplot/EAE_Skull_derived_score_macrophage.pdf",width = 18,height = 4)
FeaturePlot(multi_model_macro,features = "EAE_Blood_derived_score",split.by = "model_type",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) & NoAxes() #&  theme(legend.position = "right")
ggsave(filename = "multi_model_mice/plots/featureplot/EAE_Blood_derived_score_macrophage.pdf",width = 18,height = 4)

# vlnplot
multi_model_macro$Skull_derived_score <- multi_model_macro$EAE_Skull_derived_score
multi_model_macro$Blood_derived_score <- multi_model_macro$EAE_Blood_derived_score
meta_df <- multi_model_macro@meta.data
macrophage_df <- meta_df %>%
  filter(cell_subclass == "Macrophage") %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, orig.ident, model_type)

macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")

p <- ggplot(macrophage_long, aes(x = Score_Type, y = Score, fill = Score_Type)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",        # 显示“*”等符号
    comparisons = list(c("Skull_derived_score", "Blood_derived_score")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,              # 控制*符号显示位置
    size=18
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
    axis.title.y = element_text(size = 25, face = "bold"),  # 纵轴标题
    axis.text.x = element_text(size = 25, face = "bold",angle = 25,vjust = 0.85,hjust = 0.75),                  # 横轴刻度
    axis.text.y = element_text(size = 25, face = "bold")                   # 纵轴刻度
  )+ylim(c(0,14))

p

ggsave(filename = "multi_model_mice/plots/featureplot/EAE_Skull_vs_Blood_vlnplot_macrophage.pdf",height = 8.5,width = 4.5)

# Extracting t cells and macrophage for scoring
# skull or blood score
Idents(multi_model) <- "cell_subclass"
multi_model_t <- subset(multi_model,ident=c("T cell","Macrophage"))
DimPlot(multi_model_t,reduction = "umap")
multi_model_t <- AddModuleScore(multi_model_t,features = list(skull_derived_gene) ,name = "EAE_Skull_derived_score" )
multi_model_t <- AddModuleScore(multi_model_t,features = list(blood_derived_gene) ,name = "EAE_Blood_derived_score" )
multi_model_t@meta.data$Skull_derived_score <- rescale(multi_model_t@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
multi_model_t@meta.data$Blood_derived_score <- rescale(multi_model_t@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
FeaturePlot(multi_model_t,features = "Skull_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "~/r4.3_docker_rstudio/multi_models/results/control_sohu_crush_aion_ms_10samples/subtype/mp/plots/featureplot/EAE_Skull_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)
FeaturePlot(multi_model_t,features = "Blood_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "~/r4.3_docker_rstudio/multi_models/results/control_sohu_crush_aion_ms_10samples/subtype/mp/plots/featureplot/EAE_Blood_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)

meta_df <- multi_model_t@meta.data

t_df <- meta_df %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, model_type,cell_subclass)

t_long <- t_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")

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
p
ggsave(filename = "multi_model_mice/plots/featureplot/Skull_vs_blood_score_vlnplot_tcell&macrophage.pdf",height = 8,width = 4)

# mouse brain data
markers <- readRDS("~/r4.3_docker_rstudio/multi_models/results/control_sohu_crush_aion_ms_10samples/brain immune atlas/A single-cell atlas of mouse brain macrophages/brain_data_celltype_Markers.rds")
top_markers_list <- list()
top_n <- 30
all_clusters <- unique(markers$cluster)
sorted_markers <- markers %>%
  arrange(cluster, desc(avg_log2FC))

for (cl in all_clusters) {
  
  cluster_markers <- sorted_markers %>%
    filter(cluster == cl)

  cluster_selected <- c()
  
  for (i in 1:nrow(cluster_markers)) {
    gene_i <- cluster_markers$gene[i]
    logfc_i <- cluster_markers$avg_log2FC[i]

    gene_conflicts <- markers %>% filter(gene == gene_i)
    max_logfc <- max(gene_conflicts$avg_log2FC)
    max_cluster <- gene_conflicts$cluster[which.max(gene_conflicts$avg_log2FC)]

    if (max_cluster == cl && !(gene_i %in% selected_genes)) {
      cluster_selected <- c(cluster_selected, gene_i)
      selected_genes <- c(selected_genes, gene_i)
    }

    if (length(cluster_selected) >= top_n) {
      break
    }
  }

  top_markers_list[[cl]] <- cluster_markers %>%
    filter(gene %in% cluster_selected)
}

# mouse brain data score
bam_gene <- top_markers_list$BAM$gene[1:20]
mono_gene <- top_markers_list$monocyte$gene[1:20]

multi_model_macro <- AddModuleScore(multi_model_macro,features = list(bam_gene) ,name = "BAM_score",assay = "RNA")
multi_model_macro <- AddModuleScore(multi_model_macro,features = list(mono_gene) ,name = "Monocyte_score",assay = "RNA")
multi_model_macro@meta.data$BAM_score <- rescale(multi_model_macro@meta.data$BAM_score1, to = c(0, 10))
multi_model_macro@meta.data$Monocyte_score <- rescale(multi_model_macro@meta.data$Monocyte_score1, to = c(0, 10))
FeaturePlot(multi_model_macro,features = "BAM_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10))
ggsave(filename = "multi_model_mice/plots/featureplot/BAM_score.pdf",width = 3.3,height = 2.8)
FeaturePlot(multi_model_macro,features = "Monocyte_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10))
ggsave(filename = "multi_model_mice/plots/featureplot/Monocyte_score.pdf",width = 3.3,height = 2.8)

# vlnplotvlnplot
meta_df <- multi_model_macro@meta.data

macrophage_df <- meta_df %>%
  filter(cell_subclass == "Macrophage") %>%
  dplyr::select(BAM_score, Monocyte_score, orig.ident, model_type)

macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("BAM_score", "Monocyte_score"),
                      names_to = "Score_Type",
                      values_to = "Score")
macrophage_long$Score_Type <- factor(macrophage_long$Score_Type,levels = c( "Monocyte_score","BAM_score"))

p <- ggplot(macrophage_long, aes(x = Score_Type, y = Score, fill = Score_Type)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",        # 显示“*”等符号
    comparisons = list(c("BAM_score", "Monocyte_score")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,              # 控制*符号显示位置
    size=18
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
    axis.title.y = element_text(size = 25, face = "bold"),  # 纵轴标题
    axis.text.x = element_text(size = 25, face = "bold",angle = 25,vjust = 0.85,hjust = 0.75),                  # 横轴刻度
    axis.text.y = element_text(size = 25, face = "bold")                   # 纵轴刻度
  )+ylim(c(0,14))
p
ggsave(filename = "multi_model_mice/plots/featureplot/BAM_vs_Monocyte_score_vlnplot.pdf",height = 8.5,width = 6)


