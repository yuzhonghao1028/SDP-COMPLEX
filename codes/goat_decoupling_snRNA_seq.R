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
library(scDblFinder)
library(BiocParallel)
source("utils/utilFxns.R")
source("utils/plottingFxns.R")

## Preprocessing
# Load files and generate full object
on_sn_7dpi <- Read10X("~/r4.3_docker_rstudio/data_analysis/SBC/decompression_sn-seq/YA2023250-1_scRNA_Report/A3 injury ON1/filtered_feature_bc_matrix/")
on_sn_dec <- Read10X("~/r4.3_docker_rstudio/data_analysis/SBC/decompression_sn-seq/YA2023250-1_scRNA_Report/A1 decompression ON1/filtered_feature_bc_matrix/")

# add animal tags
colnames(on_sn_7dpi) <- paste("7dpi", colnames(on_sn_7dpi), sep = "_")
colnames(on_sn_dec) <- paste("decoupling", colnames(on_sn_dec), sep = "_")
on_mat <- cbind(on_sn_7dpi,on_sn_dec) #, on_sn_7dpi2,on_sn_7dpi3
alldata_dec <- CreateSeuratObject(on_mat, project = "goat_on_decoupling", min.cells = 3, min.features = 200)

# add animal tag to metadata
alldata_dec@meta.data[alldata_dec$orig.ident=="7dpi", 'injury_group'] = "Injury"
alldata_dec@meta.data[alldata_dec$orig.ident=="decoupling", 'injury_group'] = "Surgery"
rm(on_mat, on_sn_7dpi, on_sn_dec);gc()

# initial cluster
alldata_dec <- ClusterSeurat(alldata_dec,use.SCT=F)

# drop doublet 

sce <- as.SingleCellExperiment(alldata_dec) 
bp <- MulticoreParam(3, RNGseed=1234)
bpstart(bp) 
sce <- scDblFinder(sce, samples="orig.ident", BPPARAM=MulticoreParam(3)) 
bpstop(bp)
# 可以考虑使用BPPARAM参数对其进行多线程处理（假设有足够的RAM）
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE) 
# plotDoubletMap(sce) 
table(sce$scDblFinder.class) 
alldata_dec_sce <- as.Seurat(sce)
DimPlot(alldata_dec_sce,group.by = "scDblFinder.class")
alldata_dec <- alldata_dec_sce
rm(alldata_dec_sce,sce,bp);gc()
alldata_dec <- subset(alldata_dec, scDblFinder.class == "singlet")
Idents(alldata_dec) <- 'injury_group'
VlnPlot(alldata_dec, features = "nCount_RNA", pt.size = 0)
VlnPlot(alldata_dec, features = "nFeature_RNA", pt.size = 0)
DefaultAssay(alldata_dec) <- "RNA"

# filter cells
alldata_dec <- subset(alldata_dec, subset = nCount_RNA < 8000 & nFeature_RNA >200)  
VlnPlot(alldata_dec, features = "nCount_RNA", pt.size = 0, group.by = "injury_group")
VlnPlot(alldata_dec, features = "nFeature_RNA", pt.size = 0, group.by = "injury_group")
# saveRDS(alldata_dec,file = "rds files/all_ON_SeuratObj_drop_doublet_filtered_v2.rds")

# re-cluster
DimPlot(alldata_dec,group.by = "orig.ident",split.by = "orig.ident")
alldata_dec <- ClusterSeurat(alldata_dec,use.SCT = F)

# find DEG
dir.create("goat_decoupling_snRNA_seq/Markers/", recursive = T)
Idents(alldata_dec) <- "seurat_clusters"
all_markers <- FindAllMarkers(alldata_dec, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,assay = "RNA")
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "goat_decoupling_snRNA_seq/Markers/alldata_dec_res0.5_Markers.rds")
openxlsx::write.xlsx(all_markers,file = "goat_decoupling_snRNA_seq/Markers/alldata_dec_res0.5_Markers.xlsx",rowNames=T)

# annotation for clusters
Idents(alldata_dec) <- "seurat_clusters"
alldata_dec@meta.data$cell_class <- "Other"

alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(0,5,9)),]$cell_class = "Astrocyte"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(2,10,14)),]$cell_class = "Oligodendrocyte"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(1,4,8)),]$cell_class = "MP"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(6,12)),]$cell_class = "Fibroblast"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(3,13,17)),]$cell_class = "OPC"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(7,16)),]$cell_class = "VEC"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(15)),]$cell_class = "Mural cell"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(18)),]$cell_class = "T cell"
alldata_dec@meta.data[WhichCells(alldata_dec, idents = c(11)),]$cell_class = "Proliferative cell"

# further analysis of MP subsets
Idents(alldata_dec) <- "cell_class"
dec_mp <- subset(alldata_dec, ident ="MP")

# run clustering pipeline
DefaultAssay(dec_mp) <- "RNA"
dec_mp <- ClusterSeurat(dec_mp,cluster_resolution = 0.2) # res 0.2

# Find DE markers
dec_mp_markers <- FindAllMarkers(dec_mp,assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25,only.pos = T)
dec_mp_markers <- dec_mp_markers %>% arrange(desc(avg_log2FC))
saveRDS(dec_mp_markers, file = "goat_decoupling_snRNA_seq/Markers/dec_mp_res0.2_Markers.rds")
openxlsx::write.xlsx(dec_mp_markers, file = "goat_decoupling_snRNA_seq/Markers/dec_mp_res0.2_Markers.xlsx",rowNames=F)

# annotation
Idents(dec_mp) <- "RNA_snn_res.0.2"
dec_mp$cell_subclass <- "others"
dec_mp@meta.data[WhichCells(dec_mp, idents = c(2)),]$cell_subclass = "Microglia"
dec_mp@meta.data[WhichCells(dec_mp, idents = c(0,1,3)),]$cell_subclass = "Macrophage"

# Re-import the metadata back into multi_model by exporting the metadata and merging it in excel
meta_mp <- dec_mp@meta.data
meta_all <- alldata_dec@meta.data
write.xlsx(meta_mp,file = "goat_decoupling_snRNA_seq/meta_mp.xlsx",rowNames=T)
write.xlsx(meta_all,file = "goat_decoupling_snRNA_seq/meta_all.xlsx",rowNames=T)
meta_all <- read.xlsx("goat_decoupling_snRNA_seq/meta_all.xlsx")
alldata_dec$cell_subclass <- meta_all$cell_subclass

# plot for cell_subclass
# let intermediate MP as macrophage 
Idents(alldata_dec) <- "cell_subclass"
# alldata_dec@meta.data[alldata_dec$cell_subclass=="Intermediate MP","cell_subclass"] = "Macrophage"
seurat_cols <- c("#479D88","#3C77AF","#D1352B","#BBDD78","#EE934E","#FFBC40","#B383B9","#8FA4AE","#F5CFE4","#F6948D")
DimPlot(alldata_dec, group.by = "cell_subclass",label = T,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "goat_decoupling_snRNA_seq/plots/umapplot/umapplot_cellclass.pdf",width = 6.5,height = 5)
DimPlot(alldata_dec, group.by = "cell_subclass",label = T, split.by = "injury_group",ncol = 2,reduction = "umap",cols = seurat_cols)&NoAxes()
ggsave(filename = "goat_decoupling_snRNA_seq/plots/umapplot/umapplot_cellclass splitby_injury_group.pdf",width = 8,height = 4)

# Find celltype DE markers
Idents(alldata_dec) <- "cell_subclass"
all_markers <- FindAllMarkers(alldata_dec, logfc.threshold = 0.25, min.pct = 0.25,only.pos = T,verbose = T)
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "goat_decoupling_snRNA_seq/Markers/alldata_dec_Markers_cell_subclass.rds")
openxlsx::write.xlsx(all_markers,file = "goat_decoupling_snRNA_seq/Markers/alldata_dec_Markers_cell_subclass.xlsx",rowNames=F)

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

dir.create("goat_decoupling_snRNA_seq/plots/dotplot/")
DotPlot(alldata_dec, features = c(astro_markers,fibro_markers,macro_marker,micro_markers,mural_markers,olig_markers,opc_markers,proliferative_markers,tcell_marker,vec_markers), group.by = "cell_subclass", assay = "RNA",) + RotatedAxis()
ggsave(filename = "goat_decoupling_snRNA_seq/plots/dotplot/cell_class_marker_dotplot_coord_flip.pdf",height = 5,width = 9)

# Cell Ratio Heat Map
# set Idents
Idents(alldata_dec) <- "cell_subclass"

meta_df <- alldata_dec@meta.data

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

pdf(file = "goat_decoupling_snRNA_seq/plots/heatmap/cellsubclass_proportion_heatmap_splitby_orig.ident.pdf", 
    width = 1.5, height = 5)

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
  `Proliferative cell` = c("TOP2A", "CENPF"),
  `T cell` = c("BCL11B", "SKAP1"),
  VEC = c("FLT1", "ERG")
)

gene_to_type <- stack(marker_dict)
colnames(gene_to_type) <- c("gene", "celltype")
features <- gene_to_type$gene
types <- gene_to_type$celltype

samples <- unique(alldata_dec$orig.ident)

pct_list <- list()
logfc_list <- list()

for (s in samples) {
  message("Processing sample: ", s)

  sample_obj <- subset(alldata_dec, subset = orig.ident == s)
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

pdf(file = "goat_decoupling_snRNA_seq/plots/heatmap/marker_heatmap_splitby_orig.ident.pdf",width = 1.5,height = 6)
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

# extract macrophage for scoring
# skull or blood score
skull.genes <- c("IGKC","GM36486","PIM1","GM10076","MFGE8","EAR6","LARS2","LRP1","MT-ND3",
                 "FN1","MT-ND4L","SLC6A6","ATP13A2","C77080","MAN2B1","STAB1","PLD3","GPC1","MOBP","CTSD")
blood.genes <- c("NGP","S100A9","NRGN","CSF2","IFNG","TNFSF11","NKG7","RGS16","IL17A","LCK","GIMAP5",
                 "GZMA","CD8B1","CD4","PTPRCAP","CD3E","LAT","CTSW","TRAC","CD3G")

alldata_dec <- AddModuleScore(alldata_dec,features = list(skull.genes) ,name = "EAE_Skull_derived_score",assay = "RNA")
alldata_dec <- AddModuleScore(alldata_dec,features = list(blood.genes) ,name = "EAE_Blood_derived_score",assay = "RNA")
alldata_dec@meta.data$EAE_Skull_derived_score <- rescale(alldata_dec@meta.data$EAE_Skull_derived_score1, to = c(0, 10))
alldata_dec@meta.data$EAE_Blood_derived_score <- rescale(alldata_dec@meta.data$EAE_Blood_derived_score1, to = c(0, 10))
pal <- viridis::viridis(n = 10, option = "H", direction = 1)

# vlnplot
alldata_dec$Skull_derived_score <- alldata_dec$EAE_Skull_derived_score
alldata_dec$Blood_derived_score <- alldata_dec$EAE_Blood_derived_score

meta_df <- alldata_dec@meta.data

macrophage_df <- meta_df %>%
  filter(cell_subclass == "Macrophage") %>%
  dplyr::select(Skull_derived_score, Blood_derived_score, orig.ident, injury_group,cell_subclass)

macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("injury_group"),
                      names_to = "Score_Type",
                      values_to = "group")

ggplot(macrophage_long, aes(x = group, y = Skull_derived_score, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("Injury", "Surgery")),
    label.y = max(macrophage_df$Skull_derived_score, na.rm = TRUE) * 1.2,
    size = 6
  ) +
  facet_wrap(~cell_subclass, nrow = 1) +  # 每个细胞一格
  scale_fill_manual(values = c("#fcae59", "#db498e")) +
  labs(x = "", y = "Skull_derived_score") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold", angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 20, face = "bold")  # facet 标题
  ) + ylim(c(0, 14))

ggsave(filename = "goat_decoupling_snRNA_seq/plots/geneset score plot/Skull_derived_score_injury_vs_surgery_vlnplot_Macrophage.pdf",height = 6,width = 4)

# pathway scoring vlnplot plotting function
plot_geneset_score <- function(data, pathway_name, GOID, output_dir) {

  GOgeneID <- tryCatch({
    
    get(GOID, org.Hs.egGO2ALLEGS) %>%
      mget(org.Hs.egSYMBOL) %>%
      unlist()
  }, error = function(e) {
    message("Error with GOID: ", GOID, " - skipping this GOID.")
    return(NULL) 
  })
  
  if (is.null(GOgeneID)) {
    message("No genes found for ", pathway_name, " (GOID: ", GOID, ").")
    return(NULL)
  }

  data <- AddModuleScore(data, features = list(GOgeneID), name = paste0(pathway_name, "_score"), assay = "RNA")

  data@meta.data[[paste0(pathway_name, "_score")]] <- rescale(data@meta.data[[paste0(pathway_name, "_score1")]], to = c(0, 10))

  meta_df <- data@meta.data

  macrophage_df <- meta_df %>%
    dplyr::select(paste0(pathway_name, "_score"), injury_group, cell_subclass)
  
  macrophage_long <- macrophage_df %>%
    mutate(cell = rownames(.)) %>%
    tidyr::pivot_longer(cols = c("injury_group"),
                        names_to = "Score_Type",
                        values_to = "group")
  
  ggplot(macrophage_long, aes_string(x = "group", y = paste0(pathway_name, "_score"), fill = "group")) +
    geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
    geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      comparisons = list(c("Injury", "Surgery")),
      label.y = max(macrophage_long[[paste0(pathway_name, "_score")]], na.rm = TRUE) * 1.2,
      size = 6
    ) +
    facet_wrap(~cell_subclass, nrow = 1) +  
    scale_fill_manual(values = c("#fcae59", "#db498e")) +
    labs(x = "", y = pathway_name) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 25, face = "bold"),
      axis.text.x = element_text(size = 20, face = "bold", angle = 0, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(size = 20, face = "bold"),
      strip.text = element_text(size = 20, face = "bold")  
    ) + 
    coord_cartesian(ylim = c(0, 14))  
}

pathways <- list(
  list(name = "neuroinflammatory response", GOID = "GO:0150076"),
  list(name = "immune response_activating signaling pathway", GOID = "GO:0002757"),
  list(name = "cytokine production involved in inflammatory response", GOID = "GO:0002534"),
  list(name = "positive regulation of cell adhesion", GOID = "GO:0045785"),
  list(name = "positive regulation of lymphocyte activation", GOID = "GO:0051251")
)

# 输出目录
output_dir <- "goat_decoupling_snRNA_seq/plots/geneset score plot/gobp_pathways/"

# 循环遍历每个通路，生成并保存图像
for (pathway in pathways) {
  pathway_name <- pathway$name
  GOID <- pathway$GOID
  
  # 将空格替换为下划线
  pathway_name <- gsub(" ", "_", pathway_name)
  
  # 调用绘图函数并保存结果
  plot_geneset_score(alldata_dec, pathway_name, GOID, output_dir)
  
  # 保存图像
  ggsave(paste0(output_dir, pathway_name, "_score_injury_vs_surgery_vlnplot.pdf"), height = 6, width = 40)
}
