#Load packages
library(tidyverse)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(reshape2)
library(dplyr)
library(openxlsx)
library(scDblFinder)
library(BiocParallel)
library(gridExtra)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(future)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
library(CellChat)
source("utils/utilFxns.R")
source("utils/plottingFxns.R")

## Preprocessing
# Load files and generate full object
on_sn_ctrl1 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/goat_on_control_4samples/new/goat-normal-ON1/filtered_feature_bc_matrix/")
on_sn_ctrl2 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/goat_on_control_4samples/new/goat-normal-ON2/filtered_feature_bc_matrix/")
on_sn_7dpi1 <- Read10X("~/r4.3_docker_rstudio/multi_species optic nerve/Exp_Mtx/Gaot_normal_ON snRNA-seq 2samples 180g with 2 7dpi samples/2OS3/filtered_feature_bc_matrix/")
on_sn_7dpi2 <- Read10X("~/r4.3_docker_rstudio/multi_species optic nerve/Exp_Mtx/Gaot_normal_ON snRNA-seq 2samples 180g with 2 7dpi samples/3OS3/filtered_feature_bc_matrix/")
on_sn_1mpi1 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/1mpi/on & chiasm spatial-seq YA2024004SC-2_scRNA_Report/CellRanger/on-1mpi-1/filtered_feature_bc_matrix/")
on_sn_1mpi2 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/1mpi/on & chiasm spatial-seq YA2024004SC-2_scRNA_Report/CellRanger/on-1mpi-2/filtered_feature_bc_matrix/")
on_sn_3mpi1 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/3mpi/on & chiasm snRNA-seq YA2024010SC-3_scRNA_Report/YA2024010SC-3_scRNA_Report/cellranger/goatON-3mpi-A9/filtered_feature_bc_matrix/")
on_sn_3mpi2 <- Read10X("~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/Exp_Mtx/snRNA-seq/3mpi/on & chiasm snRNA-seq YA2024010SC-3_scRNA_Report/YA2024010SC-3_scRNA_Report/cellranger/goatON-3mpi-B3(碎片多，捕获细胞量5000)/filtered_feature_bc_matrix/")

# add animal tags
colnames(on_sn_ctrl1) <- paste("ctrl1", colnames(on_sn_ctrl1), sep = "_")
colnames(on_sn_ctrl2) <- paste("ctrl2", colnames(on_sn_ctrl2), sep = "_")
colnames(on_sn_7dpi1) <- paste("7dpi1", colnames(on_sn_7dpi1), sep = "_")
colnames(on_sn_7dpi2) <- paste("7dpi2", colnames(on_sn_7dpi2), sep = "_")
colnames(on_sn_1mpi1) <- paste("1mpi1", colnames(on_sn_1mpi1), sep = "_")
colnames(on_sn_1mpi2) <- paste("1mpi2", colnames(on_sn_1mpi2), sep = "_")
colnames(on_sn_3mpi1) <- paste("3mpi1", colnames(on_sn_3mpi1), sep = "_")
colnames(on_sn_3mpi2) <- paste("3mpi2", colnames(on_sn_3mpi2), sep = "_")
on_mat <- cbind(on_sn_ctrl1, on_sn_ctrl2, on_sn_7dpi1, on_sn_7dpi2, on_sn_1mpi1, on_sn_1mpi2, on_sn_3mpi1, on_sn_3mpi2)
alldata_merge <- CreateSeuratObject(on_mat, project = "goat_on_muti_injury", min.cells = 3, min.features = 200)

# add animal tag to metadata
alldata_merge@meta.data[alldata_merge$orig.ident=="ctrl1", 'injury_group'] = "control"
alldata_merge@meta.data[alldata_merge$orig.ident=="ctrl2", 'injury_group'] = "control"
alldata_merge@meta.data[alldata_merge$orig.ident=="7dpi1", 'injury_group'] = "7dpi"
alldata_merge@meta.data[alldata_merge$orig.ident=="7dpi2", 'injury_group'] = "7dpi"
alldata_merge@meta.data[alldata_merge$orig.ident=="1mpi1", 'injury_group'] = "1mpi"
alldata_merge@meta.data[alldata_merge$orig.ident=="1mpi2", 'injury_group'] = "1mpi"
alldata_merge@meta.data[alldata_merge$orig.ident=="3mpi1", 'injury_group'] = "3mpi"
alldata_merge@meta.data[alldata_merge$orig.ident=="3mpi2", 'injury_group'] = "3mpi"
rm(on_mat, on_sn_ctrl1, on_sn_ctrl2, on_sn_7dpi1, on_sn_7dpi2, on_sn_1mpi1, on_sn_1mpi2, on_sn_3mpi1, on_sn_3mpi2);gc()
alldata_merge@meta.data$injury_group <- factor(alldata_merge@meta.data$injury_group, levels = c("control", "7dpi", "1mpi", "3mpi"))

# initial cluster
alldata_merge <- ClusterSeurat(alldata_merge,use.SCT=F)

# drop doublet 
sce <- as.SingleCellExperiment(alldata_merge) 
bp <- MulticoreParam(3, RNGseed=1234)
bpstart(bp) 
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3)) 
bpstop(bp)
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE) 
table(sce$scDblFinder.class) 
alldata_merge_sce <- as.Seurat(sce)
DimPlot(alldata_merge_sce,group.by = "scDblFinder.class")
alldata_merge <- alldata_merge_sce
rm(alldata_merge_sce,sce,bp);gc()
alldata_merge <- subset(alldata_merge, scDblFinder.class == "singlet")
Idents(alldata_merge) <- 'injury_group'
VlnPlot(alldata_merge, features = "nCount_RNA", pt.size = 0)
VlnPlot(alldata_merge, features = "nFeature_RNA", pt.size = 0)
DefaultAssay(alldata_merge) <- "RNA"

# filter cells
alldata_merge <- subset(alldata_merge, subset = nCount_RNA < 8000 & nFeature_RNA >200)  

# re-cluster
alldata_merge <- ClusterSeurat(alldata_merge,use.SCT = F)

# find DEG
dir.create("goat_TON_snRNA_seq/Markers/", recursive = T)
Idents(alldata_merge) <- "seurat_clusters"
all_markers <- FindAllMarkers(alldata_merge, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,assay = "RNA")
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "goat_TON_snRNA_seq/Markers/alldata_merge_res0.5_Markers.rds")
openxlsx::write.xlsx(all_markers,file = "goat_TON_snRNA_seq/Markers/alldata_merge_res0.5_Markers.xlsx",rowNames=T)

# annotation for clusters
# delete heterogeneous cell populations
alldata_merge <- subset(alldata_merge,ident=19,invert=T)
alldata_merge$seurat_clusters <- droplevels(alldata_merge$seurat_clusters)

Idents(alldata_merge) <- "seurat_clusters"
alldata_merge@meta.data$cell_class <- "Other"

alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(0,3,9)),]$cell_class = "Astrocyte"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(1,2,6)),]$cell_class = "Oligodendrocyte"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(4,8)),]$cell_class = "Immune cell"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(7,10,12,14)),]$cell_class = "Fibroblast"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(11,15)),]$cell_class = "OPC"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(5,16)),]$cell_class = "VEC"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(13)),]$cell_class = "Mural cell"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(17)),]$cell_class = "Muscle cell"
alldata_merge@meta.data[WhichCells(alldata_merge, idents = c(18)),]$cell_class = "Proliferative cell"

# further analysis of MP subsets
Idents(alldata_merge) <- "cell_class"
alldata_mp <- subset(alldata_merge, ident ="Immune cell")

# run clustering pipeline
DefaultAssay(alldata_mp) <- "RNA"
alldata_mp <- ClusterSeurat(alldata_mp,cluster_resolution = 0.2) # res 0.2

# Find DE markers
alldata_mp_markers <- FindAllMarkers(alldata_mp,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
alldata_mp_markers <- alldata_mp_markers %>% arrange(desc(avg_log2FC))
saveRDS(alldata_mp_markers, file = "goat_TON_snRNA_seq/Markers/alldata_mp_res0.2_Markers.rds")
openxlsx::write.xlsx(alldata_mp_markers, file = "goat_TON_snRNA_seq/Markers/alldata_mp_res0.2_Markers.xlsx",rowNames=F)

# annotation
Idents(alldata_mp) <- "seurat_clusters"
alldata_mp$cell_subclass <- "others"
alldata_mp@meta.data[WhichCells(alldata_mp, idents = c(1,3)),]$cell_subclass = "Microglia"
alldata_mp@meta.data[WhichCells(alldata_mp, idents = c(0,2,4,5,7)),]$cell_subclass = "Macrophage"
alldata_mp@meta.data[WhichCells(alldata_mp, idents = c(6)),]$cell_subclass = "T cell"
alldata_mp@meta.data[WhichCells(alldata_mp, idents = c(8)),]$cell_subclass = "Mast cell"

# Re-import the metadata back into multi_model by exporting the metadata and merging it in excel
meta_mp <- alldata_mp@meta.data
meta_all <- alldata_merge@meta.data
write.xlsx(meta_mp,file = "goat_TON_snRNA_seq/meta_mp.xlsx",rowNames=T)
write.xlsx(meta_all,file = "goat_TON_snRNA_seq/meta_all.xlsx",rowNames=T)
meta_all <- read.xlsx("goat_TON_snRNA_seq/meta_all.xlsx")
alldata_merge$cell_subclass <- meta_all$cell_subclass

# plot for cell_subclass
Idents(alldata_merge) <- "cell_subclass"
seurat_cols <- c("#479D88","#3C77AF","#D1352B","#AECDE1","#BBDD78","#EE934E","#9B5B33","#FFBC40","#B383B9","#8FA4AE","#F5CFE4","#F6948D")
DimPlot(alldata_merge, group.by = "cell_subclass",label = T,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "goat_TON_snRNA_seq/plots/umapplot/umapplot_cellclass.pdf",width = 6.5,height = 5)
DimPlot(alldata_merge, group.by = "cell_subclass",label = T, split.by = "injury_group",ncol = 2,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "goat_TON_snRNA_seq/plots/umapplot/umapplot_cellclass splitby_injury_group.pdf",width = 9,height = 8)
alldata_merge$orig.ident <- factor(alldata_merge$orig.ident, levels = c("ctrl1","ctrl2","7dpi1","7dpi2","1mpi1","1mpi2","3mpi1","3mpi2"))
DimPlot(alldata_merge, group.by = "cell_subclass",label = T, split.by = "orig.ident", ncol = 4,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "goat_TON_snRNA_seq/plots/umapplot/umapplot_cellclass splitby_orig.ident.pdf",width = 14,height = 10)

# Find celltype DE markers
Idents(alldata_merge) <- "cell_subclass"
all_markers <- FindAllMarkers(alldata_merge, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,verbose = T)
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "goat_TON_snRNA_seq/Markers/alldata_merge_Markers_cell_subclass.rds")
openxlsx::write.xlsx(all_markers,file = "goat_TON_snRNA_seq/Markers/alldata_merge_Markers_cell_subclass.xlsx",rowNames=F)

# celltype marker dotplot
olig_markers= c("MBP","MOG")
astro_markers=c("GFAP","AQP4")
fibro_markers=c("FN1","COL1A2")
opc_markers=c("PDGFRA","BRINP3")
vec_markers=c("FLT1","ERG")
mural_markers=c("ACTA2","CALD1")
muscle_cell_markers=c("TNNT3","RYR1")
proliferative_markers=c("TOP2A","CENPF")
micro_markers=c("P2RY12", "SRGAP2")
macro_marker=c("CD36","MSR1")
mast_marker=c("KIT","CPA3")
tcell_marker=c("BCL11B","SKAP1")

dir.create("goat_TON_snRNA_seq/plots/dotplot/",recursive = T)
DotPlot(alldata_merge, features = c(astro_markers,fibro_markers,macro_marker,mast_marker,micro_markers,mural_markers,muscle_cell_markers,olig_markers,opc_markers,proliferative_markers,tcell_marker,vec_markers), group.by = "cell_subclass", assay = "RNA",) + RotatedAxis()
ggsave(filename = "goat_TON_snRNA_seq/plots/dotplot/cell_class_marker_dotplot_cellsubclass_coord_flip.pdf",height = 5,width = 9)

# Cell Ratio Heat Map
# set Idents
Idents(alldata_merge) <- "cell_subclass"
meta_df <- alldata_merge@meta.data

prop_df <- meta_df %>%
  group_by(orig.ident, cell_subclass) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

prop_mat <- dcast(prop_df, cell_subclass ~ orig.ident, value.var = "proportion", fill = 0)
rownames(prop_mat) <- prop_mat$cell_subclass
prop_mat <- as.matrix(prop_mat[, -1])  
prop_mat <- prop_mat[rev(1:nrow(prop_mat)),] 

col_fun_prop <- circlize::colorRamp2(c(0, max(prop_mat)), c("grey90", "blue"))

dir.create("goat_TON_snRNA_seq/plots/heatmap/")
pdf(file = "goat_TON_snRNA_seq/plots/heatmap/cellsubclass_proportion_heatmap_splitby_orig.ident.pdf", 
    width = 3.5, height = 5)

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
  `Mast cell` = c("KIT", "CPA3"),
  Microglia = c("P2RY12", "SRGAP2"),
  `Mural cell` = c("ACTA2", "CALD1"),
  `Muscle cell` = c("TNNT3", "RYR1"),
  Oligodendrocyte = c("MBP", "MOG"),
  OPC = c("PDGFRA", "BRINP3"),
  `Proliferative cell` = c("TOP2A", "CENPF"),
  `T cell` = c("BCL11B", "SKAP1"),
  VEC = c("FLT1", "ERG")
)

gene_to_type <- stack(marker_dict)
colnames(gene_to_type) <- c("gene", "celltype")
features <- gene_to_type$gene
types <- gene_to_type$celltype

samples <- unique(alldata_merge$orig.ident)

pct_list <- list()
logfc_list <- list()

for (s in samples) {
  message("Processing sample: ", s)

  sample_obj <- subset(alldata_merge, subset = orig.ident == s)
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

pdf(file = "goat_TON_snRNA_seq/plots/heatmap/marker_heatmap_splitby_orig.ident.pdf",width = 6,height = 6)
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

# cell marker heatmap
features <- c(astro_markers,fibro_markers,macro_marker,mast_marker,micro_markers,mural_markers,muscle_cell_markers,olig_markers,proliferative_markers,tcell_marker,vec_markers)
options(future.globals.maxSize = 10 * 1024^3)
expr_matrix <- GetAssayData(alldata_merge_scaled, layer = "scale.data",assay = "RNA")[features, ]
Idents(alldata_merge) <- "cell_subclass"
cell_types <- Idents(alldata_merge)  
expr_avg <- aggregate(t(expr_matrix), by = list(cell_types), FUN = mean)
rownames(expr_avg) <- expr_avg$Group.1  
expr_avg <- as.matrix(expr_avg[, -1])  
unique_cell_types <- rownames(expr_avg)
extract.col(DimPlot(alldata_merge,group.by = "cell_subclass"))
seurat.cols <- seurat_cols
if (length(unique_cell_types) > length(seurat.cols)) {
  extra_colors <- grDevices::rainbow(length(unique_cell_types) - length(seurat.cols))
  cell_type_colors <- c(seurat.cols, extra_colors)
} else {
  cell_type_colors <- seurat.cols[1:length(unique_cell_types)]
}

cell_type_colors <- structure(cell_type_colors, names = unique_cell_types)

ha <- rowAnnotation(
  Cell_Type = anno_simple(rownames(expr_avg), col = cell_type_colors)
)

col_fun <- colorRamp2(c(-2, 0, 2), brewer.pal(9, "Purples")[c(2, 5, 9)])

expr_avg <- expr_avg[c(4,5,7,10,1,9,11,3,2,12,8,6),]

pdf(file ="goat_TON_snRNA_seq/plots/heatmap/cell_subclass_marker_heatmap.pdf",width = 11,height = 5 )
Heatmap(
  expr_avg,
  name = "Avg. Scaled Exp.",
  col = col_fun,
  show_row_names = TRUE,  
  show_column_names = TRUE, 
  cluster_rows = FALSE,  
  cluster_columns = FALSE,  
  left_annotation = ha,  
)
dev.off()

# celltype propotion stackbarplot
data <- as.data.frame(table(alldata_merge$cell_subclass,alldata_merge$injury_group))
data <- data %>% group_by(Var2) %>% mutate(percentage = Freq/sum(Freq)) %>% mutate(label = Freq/sum(Freq)*100)
data$label <- paste(sprintf("%.1f", data$label),"%")
extract.col(DimPlot(alldata_merge,group.by = "cell_subclass"))
seurat.cols <- seurat_cols
ggplot(data, aes( x = Var2, y=percentage,fill = Var1))+
  geom_col(position = 'stack', width = 0.8)+
  theme_bw()+
  scale_fill_manual(values=seurat.cols)+ 
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Percentage",
       fill="celltype",title="")+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 0, hjust = 0.5,
                             vjust = 0),
    legend.title=element_text(size=12), 
    legend.text=element_text(size=12))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
  )+ 
  guides(fill=guide_legend(keywidth = 1, keyheight = 1))+ 
  ylim(c(-0.005,1.005))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  geom_text(aes(label=label),size=4,color="black",position = "stack")

ggsave(filename = "goat_TON_snRNA_seq/plots/proportion plot/stack barplot_cell_subclass_short_size_change.pdf",width = 6,height = 4)

# line plot
cell_types_to_plot <- c("Macrophage", "Microglia", "T cell","Mast cell")

counts <- table(alldata_merge@meta.data$cell_subclass, alldata_merge@meta.data$injury_group)
counts <- as.data.frame(counts)
colnames(counts) <- c("cell_class", "injury_group", "Freq")

counts <- counts %>%
  group_by(injury_group) %>%
  mutate(Percent = Freq / sum(Freq) * 100) %>%
  ungroup() %>%
  filter(cell_class %in% cell_types_to_plot)

seurat_cols <- seurat_cols[c(3,5,4,11)]

p <- ggplot(counts, aes(x = injury_group, y = Percent, color = cell_class, group = cell_class)) +
  geom_line(size = 3) +                          
  geom_point(size = 3, shape = 21, fill = "white") +  
  labs(x = "", y = "Percentage (%)", color = "Cell Type") +
  scale_color_manual(values = seurat_cols) + 
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

p

ggsave(plot = p, filename = "goat_TON_snRNA_seq/plots/proportion plot/cellclass_lineplot_selected.pdf", height = 4, width = 6)

# extract macrophage for scoring
# skull or blood score
skull.genes <- c("IGKC","GM36486","PIM1","GM10076","MFGE8","EAR6","LARS2","LRP1","MT-ND3",
                 "FN1","MT-ND4L","SLC6A6","ATP13A2","C77080","MAN2B1","STAB1","PLD3","GPC1","MOBP","CTSD")
blood.genes <- c("NGP","S100A9","NRGN","CSF2","IFNG","TNFSF11","NKG7","RGS16","IL17A","LCK","GIMAP5",
                 "GZMA","CD8B1","CD4","PTPRCAP","CD3E","LAT","CTSW","TRAC","CD3G")
Idents(alldata_merge) <- "cell_subclass"
alldata_merge_macrophage <- subset(alldata_merge,ident="Macrophage")
alldata_merge_macrophage <- AddModuleScore(alldata_merge_macrophage,features = list(skull.genes) ,name = "EAE_Skull_derived_score",assay = "RNA")
alldata_merge_macrophage <- AddModuleScore(alldata_merge_macrophage,features = list(blood.genes) ,name = "EAE_Blood_derived_score",assay = "RNA")
alldata_merge_macrophage@meta.data$EAE_Skull_derived_score <- rescale(alldata_merge_macrophage@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
alldata_merge_macrophage@meta.data$EAE_Blood_derived_score <- rescale(alldata_merge_macrophage@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
pal <- viridis::viridis(n = 10, option = "H", direction = 1)

Idents(alldata_merge_macrophage) <- "injury_group"
alldata_mp_7dpi <- subset(alldata_merge_macrophage,ident=c("control","7dpi"))

FeaturePlot(alldata_mp_7dpi,features = "EAE_Skull_derived_score",split.by = "injury_group",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) & NoAxes() #&  theme(legend.position = "right")
ggsave(filename = "goat_TON_snRNA_seq/plots/geneset score plot/EAE_Skull_derived_score_order_true_7dpi&ctrl.pdf",width = 8,height = 4)
FeaturePlot(alldata_mp_7dpi,features = "EAE_Blood_derived_score",split.by = "injury_group",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10))  & NoAxes() #&  theme(legend.position = "right")
ggsave(filename = "goat_TON_snRNA_seq/plots/geneset score plot/EAE_Blood_derived_score_order_true_7dpi&ctrl.pdf",width = 8,height = 4)

# vlnplot
alldata_mp_7dpi$Skull_derived_score <- alldata_mp_7dpi$EAE_Skull_derived_score
alldata_mp_7dpi$Blood_derived_score <- alldata_mp_7dpi$EAE_Blood_derived_score
meta_df <- alldata_mp_7dpi@meta.data

macrophage_df <- meta_df %>%
  filter(cell_subclass == "Macrophage") %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, orig.ident, injury_group)

macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")

library(ggpubr)
p <- ggplot(macrophage_long, aes(x = Score_Type, y = Score, fill = Score_Type)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",     
    comparisons = list(c("Skull_derived_score", "Blood_derived_score")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,           
    size=18
  ) +
  scale_fill_manual(values = c("#fcae59", "#db498e")) +
  labs(x = "", y = "Score") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")+
  theme_classic2(base_size = 14) + 
  theme(legend.position = "none")+  
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20, face = "bold"),  
    axis.title.y = element_text(size = 25, face = "bold"),  
    axis.text.x = element_text(size = 25, face = "bold",angle = 25,vjust = 0.5,hjust = 0.5),                 
    axis.text.y = element_text(size = 25, face = "bold")                   
  )+ylim(c(0,14))

ggsave(filename = "goat_TON_snRNA_seq/plots/geneset score plot/EAE_Skull_vs_Blood_vlnplot_7dpi&ctrl.pdf",height = 10,width = 7.2)

# Extracting t cells and macrophage for scoring
# skull or blood score
Idents(alldata_merge) <- "cell_subclass"
alldata_t <- subset(alldata_merge,ident=c("T cell","Macrophage"))
alldata_t <- ClusterSeurat(alldata_t,use.SCT = F)
DimPlot(alldata_t,reduction = "umap",group.by = "cell_subclass")
alldata_t <- AddModuleScore(alldata_t,features = list(skull.genes) ,name = "EAE_Skull_derived_score" )
alldata_t <- AddModuleScore(alldata_t,features = list(blood.genes) ,name = "EAE_Blood_derived_score" )
alldata_t@meta.data$Skull_derived_score <- rescale(alldata_t@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
alldata_t@meta.data$Blood_derived_score <- rescale(alldata_t@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
FeaturePlot(alldata_t,features = "Skull_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/results/GZB_reanalysis_20250324/snRNA-seq_drop2ctrl/goat_TON_snRNA_seq/plots/geneset score plot/EAE_Skull_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)
FeaturePlot(alldata_t,features = "Blood_derived_score",order = T,reduction = "umap") & scale_color_gradientn(colors = pal, limits = c(0, 10)) &  theme(legend.position = "right")
ggsave(filename = "~/r4.3_docker_rstudio/multi_injury_time_optic_nerve/results/GZB_reanalysis_20250324/snRNA-seq_drop2ctrl/goat_TON_snRNA_seq/plots/geneset score plot/EAE_Blood_derived_score_order_true_tcell&macrophage.pdf",width = 3.3,height = 2.8)

# vlnplot
meta_df <- alldata_t@meta.data

t_df <- meta_df %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, injury_group,cell_subclass)

t_long <- t_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Skull_derived_score", "Blood_derived_score"),
                      names_to = "Score_Type",
                      values_to = "Score")
ggplot(t_long, aes(x = cell_subclass, y = Score, fill = cell_subclass)) +
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
ggsave(filename = "goat_TON_snRNA_seq/plots/geneset score plot/Skull_vs_blood_score_vlnplot_tcell&macrophage.pdf",height = 8,width = 4)

# volcano plot
## Macrophage
markers <- readRDS("goat_TON_snRNA_seq/Markers/alldata_merge_Markers_cell_subclass.rds")
data<-markers[markers$cluster=="Macrophage",]
data <- data %>%
  mutate(gene_type = case_when(avg_log2FC >= 0.5 & p_val <= 0.05 ~ "up",
                               avg_log2FC <= -0.5 & p_val <= 0.05 ~ "down",
                               TRUE ~ "ns")) 
data$log10_p <- -log10(pmax(data$p_val, 1e-320))
sig_il_genes <- rbind(head(data,n=5),tail(data,n=5))
up_il_genes <- head(data,n=5)
down_il_genes <- tail(data,n=5)
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)       
ggplot(data = data,aes(x = avg_log2FC,y = log10_p)) + 
  geom_point(aes(colour = gene_type), 
             alpha = 1, 
             shape = 16,
             size = 2) + 
  geom_point(data = up_il_genes,
             shape = 21,
             size = 3, 
             fill = "firebrick", 
             colour = "black") + 
  geom_point(data = down_il_genes,
             shape = 21,
             size = 3, 
             fill = "steelblue", 
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed") +
  geom_label_repel(data = sig_il_genes,   
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1) +
  scale_colour_manual(values = cols) + 
  scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                     limits = c(-10, 10)) +
  labs(title = "Macrophage DEGs",
       x = "log2(fold change)",
       y = "-log10(P-value)",
       colour = "Expression \nchange") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
ggsave(filename = "goat_TON_snRNA_seq/plots/volcano/macrophage DEGs volcano_combine_im_7dpi.pdf",height = 5,width = 7)

# enrichment analysis
Idents(alldata_merge) <- "injury_group"
alldata_merge_7dpi <- subset(alldata_merge,ident="7dpi")
Idents(alldata_merge_7dpi) <- "cell_subclass" 
markers <- FindAllMarkers(alldata_merge_7dpi,assay = "RNA", min.pct = 0.25, logfc.threshold = 0,only.pos = F)
ids=bitr(markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
markers=merge(markers,ids,by.x='gene',by.y='SYMBOL')
markers <- markers %>% arrange(desc(avg_log2FC))
markers_fc1 <- markers[markers$avg_log2FC >= 0.5,]
markers_fc1 <- markers_fc1[markers_fc1$p_val <= 0.05,]
table(markers_fc1$cluster)
# split as list
marker_list <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  summarise(
    gene_vec = list(setNames(avg_log2FC, ENTREZID))
  ) %>%
  deframe() 
marker_list_fc1 <- markers_fc1 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  summarise(
    gene_vec = list(setNames(avg_log2FC, ENTREZID))
  ) %>%
  deframe()  

str(marker_list)

# GO富集
ego_macro <- enrichGO(gene    = names(marker_list_fc1$Macrophage),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP", # "BP","MF","ALL"
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

ego_micro <- enrichGO(gene    = names(marker_list_fc1$Microglia),
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP", # "BP","MF","ALL"
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
# check
t(ego_macro@result[1,])
t(ego_micro@result[1,])
# save rds file
saveRDS(ego_macro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GO_BP.rds")
saveRDS(ego_micro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GO_BP.rds")
# save xlsx file
data <- ego_macro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GO_BP.xlsx")
data <- ego_micro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GO_BP.xlsx")

# GO 去冗余
library(GOSemSim)
dim(ego_macro)
# [1] 355   9
dim(ego_micro)
# [1] 54   9

ego_macro_sim <- clusterProfiler::simplify(ego_macro, cutoff=0.7, measure = "Wang",
                                           by="p.adjust", select_fun=min)
ego_micro_sim <- clusterProfiler::simplify(ego_micro, cutoff=0.7, measure = "Wang",
                                           by="p.adjust", select_fun=min)
dim(ego_macro_sim)
# [1] 157   9
dim(ego_micro_sim)
# [1] 39  9


# save rds file
saveRDS(ego_macro_sim,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GO_BP_simplify.rds")
saveRDS(ego_micro_sim,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GO_BP_simplify.rds")
# save xlsx file
data <- ego_macro_sim@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GO_BP_simplify.xlsx")
data <- ego_micro_sim@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GO_BP_simplify.xlsx")

# GSEA
gseago_macro <- gseGO(geneList     = marker_list$Macrophage,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
gseago_macro = setReadable(gseago_macro, OrgDb = org.Hs.eg.db)
t(gseago_macro@result[1,])
gseago_micro <- gseGO(geneList     = marker_list$Microglia,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP",
                      minGSSize    = 100,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
gseago_micro = setReadable(gseago_micro, OrgDb = org.Hs.eg.db)
t(gseago_micro@result[1,])

# save rds file
saveRDS(gseago_macro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GSEA_GO_BP.rds")
saveRDS(gseago_micro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GSEA_GO_BP.rds")
# save xlsx file
data <- gseago_macro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GSEA_GO_BP.xlsx")
data <- gseago_micro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GSEA_GO_BP.xlsx")

# KEGG
ekg_macro <- enrichKEGG(gene         = names(marker_list_fc1$Macrophage),
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
ekg_macro <- setReadable(ekg_macro, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
t(ekg_macro@result[1,])

ekg_micro <- enrichKEGG(gene         = names(marker_list_fc1$Microglia),
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
ekg_micro <- setReadable(ekg_micro, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
t(ekg_micro@result[1,])

# save rds file
saveRDS(ekg_macro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/KEGG.rds")
saveRDS(ekg_micro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/KEGG.rds")
# save xlsx file
data <- ekg_macro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/KEGG.xlsx")
data <- ekg_micro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/KEGG.xlsx")

# GSEA_kegg
gseakg_macro <- gseKEGG(geneList     = marker_list$Macrophage,
                        organism     = 'hsa',
                        minGSSize    = 10,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
gseakg_macro <- setReadable(gseakg_macro, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
t(gseakg_macro@result[1,])

gseakg_micro <- gseKEGG(geneList     = marker_list$Microglia,
                        organism     = 'hsa',
                        minGSSize    = 10,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
gseakg_micro <- setReadable(gseakg_micro, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
t(gseakg_micro@result[1,])

# save rds file
saveRDS(gseakg_macro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GSEA_KEGG.rds")
saveRDS(gseakg_micro,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GSEA_KEGG.rds")
# save xlsx file
data <- gseakg_macro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GSEA_KEGG.xlsx")
data <- gseakg_micro@result
openxlsx::write.xlsx(data,file = "goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GSEA_KEGG.xlsx")

# Visualization of enrichment results
plot_enrichment_rds <- function(root_dir) {
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
  
  cell_folders <- list.dirs(root_dir, full.names = TRUE, recursive = FALSE)
  
  for (folder in cell_folders) {
    message("Processing: ", folder)
    
    rds_files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)
    
    for (rds_path in rds_files) {
      rds_name <- tools::file_path_sans_ext(basename(rds_path))
      message("  -> Reading ", rds_name)
      
      rds_obj <- readRDS(rds_path)
      if (is.null(rds_obj) || length(rds_obj) == 0) next
      
      if (!startsWith(rds_name, "GSEA")) {
        p <- dotplot(rds_obj, showCategory = 20,label_format=100) +
          ggtitle(rds_name)
        ggsave(filename = file.path(folder, paste0(rds_name, "_dotplot.pdf")),
               plot = p, width = 10, height = 8)
      } else {
        
        gsea_subdir <- file.path(folder, rds_name)
        dir.create(gsea_subdir, showWarnings = FALSE)
        
        if (!"Description" %in% colnames(rds_obj@result)) next
        
        for (i in seq_len(nrow(rds_obj@result))) {
          term <- rds_obj@result[i, ]
          term_title <- gsub("[/:*?\"<>|]", "_", term$Description)
          
          p <- tryCatch({
            gseaplot(rds_obj, geneSetID = i, title = term$Description)
          }, error = function(e) NULL)
          
          if (!is.null(p)) {
            ggsave(filename = file.path(gsea_subdir, paste0(term_title, ".pdf")),
                   plot = p, width = 7, height = 5)
          }
        }
      }
    }
  }
}
plot_enrichment_rds("goat_TON_snRNA_seq/enrichment_cell_subclass/")

# GO BP simplify enrichment barplot
## macrophage
data<-read.csv("goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/macrophage_7dpi_GO_BP_selected.csv")
font.size =12
data %>% 
  ggplot(aes(x=forcats::fct_reorder(Description,pvalue,.desc = T),y=-log10(pvalue)))+ 
  geom_bar(stat="identity",width = 0.7)+
  coord_flip()+
  labs(x=NULL) +
  ggtitle("")+
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))+
  geom_col(fill="orange",width = 0.7)
ggsave(filename = "goat_TON_snRNA_seq/plots/enrichment/macrophage GO BP barplot in 7dpi vs control.pdf",width = 5,height = 3.5)

## microglia
data<-read.csv("goat_TON_snRNA_seq/enrichment_cell_subclass/7dpi_logfc0.5/microglia/microglia_7dpi_GO_BP_selected.csv")
font.size =12
data %>% 
  ggplot(aes(x=forcats::fct_reorder(Description,pvalue,.desc = T),y=-log10(pvalue)))+ 
  geom_bar(stat="identity",width = 0.7)+
  coord_flip()+
  labs(x=NULL) +
  ggtitle("")+
  theme_classic() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))+
  geom_col(fill="orange",width = 0.7)
ggsave(filename = "goat_TON_snRNA_seq/plots/enrichment/microglia GO BP barplot in 7dpi vs control.pdf",width = 6,height = 3.5)

# Visualization of GSEA GOBP
gseaRes_macro <- readRDS("goat_TON_snRNA_seq/enrichment_cell_subclass/macrophage/GSEA_GO_BP.rds")
gseaRes_micro <- readRDS("goat_TON_snRNA_seq/enrichment_cell_subclass/microglia/GSEA_GO_BP.rds")
m_gsea_list <- list(gseaRes_macro,gseaRes_micro)

geneSetID <- c("GO:0045321",	"GO:0050778","GO:0001816")
gseaNb(object = gseaRes_macro,
       geneSetID = geneSetID,
       subPlot = 2,
       termWidth = 35,
       addPval = T,
       legend.position = c(0.55,0.8),
       pvalX = 0.98 ,pvalY = 0.9)
ggsave("goat_TON_snRNA_seq/plots/enrichment/macropahge_GSEA_immune.pdf",height = 7,width = 7.5)

GSEAmultiGP(gsea_list = m_gsea_list,
            geneSetID = "GO:0045321",
            exp_name = c("Macrophage","Microglia"),
            curve.col = ggsci::pal_d3()(2)[c(2,1)],
            addPval = T,
            pvalX = 0.99,pvalY = 0.99,
            legend.position = c(0.9,0.7))
ggsave("goat_TON_snRNA_seq/plots/enrichment/macropahge_microglia_GSEA_leukocyte activation.pdf",height = 7,width = 7.5)

# chemokine pairs daotplot
## ligand
Idents(alldata_merge) <- "injury_group"
alldata_merge_7dpi <- subset(alldata_merge,ident="7dpi")
Idents(alldata_merge) <- "cell_subclass"
# 保留关键的几个胶质细胞
Idents(alldata_merge_7dpi) <- "cell_subclass"
alldata_merge_7dpi <- subset(alldata_merge_7dpi,ident='T cell',invert=T)
alldata_merge_7dpi <- subset(alldata_merge_7dpi,ident='Mast cell',invert=T)
alldata_merge_7dpi <- subset(alldata_merge_7dpi,ident='Muscle cell',invert=T)

DotPlot(alldata_merge_7dpi,features = c("CXCL12",
                                        "CX3CL1",
                                        "LOC102182115",
                                        "LOC102174969",
                                        "LOC102180880",
                                        "LOC102175889",
                                        "CCL5",
                                        "LOC102175436"),cols = c("grey", "#ff70b5"))+RotatedAxis()+coord_flip()
ggsave(filename = "goat_TON_snRNA_seq/plots/dotplot/ligand dotplot_only.7dpi.pdf",width = 7,height = 5.5)
## receptor
DotPlot(alldata_merge_7dpi,features = c("CXCR4",
                                        "CX3CR1",
                                        "CXCR2",
                                        "LOC102179993",
                                        "LOC102178953",
                                        "CCR5"),cols = c("grey", "#ff70b5"))+RotatedAxis()+coord_flip()
ggsave(filename = "goat_TON_snRNA_seq/plots/dotplot/receptor dotplot_only.7dpi.pdf",width = 7,height = 5.5)

# cellchat analysis
alldata_mp_olig.opc <- subset(alldata_merge,ident=c("Macrophage","Microglia","OPC","Oligodendrocyte"))
# cellchat analysis (7dpi)
Idents(alldata_mp_olig.opc) <- "injury_group"
alldata_mp_olig.opc_7dpi <- subset(alldata_mp_olig.opc, ident="7dpi")
alldata_mp_olig.opc_7dpi$injury_group <- droplevels(alldata_mp_olig.opc_7dpi$injury_group)
table(alldata_mp_olig.opc_7dpi$injury_group)

Idents(alldata_mp_olig.opc_7dpi) <- "cell_subclass"
data_input <- alldata_mp_olig.opc_7dpi[["RNA"]]@data
meta <- alldata_mp_olig.opc_7dpi@meta.data
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "cell_subclass")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
rm(CellChatDB,CellChatDB.use)
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4) 
options(future.globals.maxSize= 10*1024^3) #写多少数字就是多少GB
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",population.size = T) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.table(df.net,paste0('goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat.df.net_7dpi.xls'),sep = "\t",row.names = F)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file = "goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_7dpi_revised.rds")

## cellchat analysis (control)
Idents(alldata_mp_olig.opc) <- "injury_group"
alldata_mp_olig.opc_ctrl <- subset(alldata_mp_olig.opc, ident="control")
alldata_mp_olig.opc_ctrl$injury_group <- droplevels(alldata_mp_olig.opc_ctrl$injury_group)
table(alldata_mp_olig.opc_ctrl$injury_group)

Idents(alldata_mp_olig.opc_ctrl) <- "cell_subclass"
data_input <- alldata_mp_olig.opc_ctrl[["RNA"]]@data
meta <- alldata_mp_olig.opc_ctrl@meta.data
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "cell_subclass")
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
rm(CellChatDB,CellChatDB.use)
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4) 
options(future.globals.maxSize= 10*1024^3) #写多少数字就是多少GB
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = T,type = "truncatedMean",population.size = T) # revised version
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.table(df.net,paste0('goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat.df.net_ctrl_revised.xls'),sep = "\t",row.names = F)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file = "goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_ctrl.rds")

## cellchat analysis (7dpi vs control)
# Load CellChat object of each dataset and merge them together
cellchat.ctrl <- readRDS("goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_ctrl_revised.rds")
cellchat.ctrl <- netAnalysis_computeCentrality(cellchat.ctrl, slot.name = "netP")
cellchat.7dpi <- readRDS("goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_7dpi_revised.rds")
cellchat.7dpi <- netAnalysis_computeCentrality(cellchat.7dpi, slot.name = "netP")
object.list <- list(injury = cellchat.7dpi,control = cellchat.ctrl)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
# (B) Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

# Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Macrophage")
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Microglia")
# netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Intermediate MP")
gg1

# Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Macrophage", targets.use = c("OPC","Oligodendrocyte"), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Macrophage", targets.use = c("OPC","Oligodendrocyte"), stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave("goat_TON_snRNA_seq/cellchat_mp_olig.opc/information flow of each signaling pathway_macrophage_to_OLs_lineage_7dpi_vs_control.pdf",height = 5,width = 8)

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Microglia", targets.use = c("OPC","Oligodendrocyte"), stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = "Microglia", targets.use = c("OPC","Oligodendrocyte"), stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave("goat_TON_snRNA_seq/cellchat_mp_olig.opc/information flow of each signaling pathway_Microglia_to_OLs_lineage_7dpi_vs_control.pdf",height = 5,width = 7)

# Identify dysfunctional signaling by comparing the communication probabities
gg1 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(3,4),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in 7dpi", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1), targets.use = c(3,4),  comparison = c(1, 2), max.dataset = 2, title.name = "Decreased signaling in 7dpi", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
ggsave("goat_TON_snRNA_seq/cellchat_mp_olig.opc/netVisual_bubble_Macrophage_to_OLs_linage_7dpi_vs_ctrl.pdf",height = 16,width = 13)

# Identify dysfunctional signaling by comparing the communication probabities
gg1 <- netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(3,4),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in 7dpi", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(3,4),  comparison = c(1, 2), max.dataset = 2, title.name = "Decreased signaling in 7dpi", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
ggsave("goat_TON_snRNA_seq/cellchat_mp_olig.opc/netVisual_bubble_Microglia_to_OLs_linage_7dpi_vs_ctrl.pdf",height = 12,width = 13)

# Comparison of outgoing/incoming signals associated with each cell population
cellchat.ctrl <- readRDS("goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_ctrl_revised.rds")
cellchat.ctrl <- netAnalysis_computeCentrality(cellchat.ctrl, slot.name = "netP")
cellchat.7dpi <- readRDS("goat_TON_snRNA_seq/cellchat_mp_olig.opc/sn_cellchat_7dpi_revised.rds")
cellchat.7dpi <- netAnalysis_computeCentrality(cellchat.7dpi, slot.name = "netP")
object.list <- list(injury = cellchat.7dpi,control = cellchat.ctrl)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
i = 1

pathway.union <- c('COLLAGEN','ADGRL','APP','FN1','HGF','LAMININ','NRXN','ADGRG','SEMA3','SerotoninDopamin',
                   'DHEA','CD99','Cholesterol','Glutamate','MPZ','SEMA4','GRN','CADM','PSAP','CD39','PDGF',
                   'SPP1','SEMA6','ADGRB','ApoA',"NCAM",'SLIT','ADGRL','Netrin')
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "incoming", signaling = pathway.union, #pattern = "incoming"
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "incoming", signaling = pathway.union, #pattern = "incoming"
                                        title = names(object.list)[i+1], width = 5, height = 6)
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], 
                                        pattern = "outgoing", signaling = pathway.union, #pattern = "incoming"
                                        title = names(object.list)[i], width = 5, height = 6,color.heatmap = "OrRd")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], 
                                        pattern = "outgoing", signaling = pathway.union, #pattern = "incoming"
                                        title = names(object.list)[i+1], width = 5, height = 6,color.heatmap = "OrRd")
pdf("goat_TON_snRNA_seq/cellchat_mp_olig.opc/singaling heatmap_MP_to_OLs_lineage.pdf",width = 23,height = 9)
draw(ht1 +ht3+ ht2  +ht4, ht_gap = unit(0.5, "cm"))
dev.off()

