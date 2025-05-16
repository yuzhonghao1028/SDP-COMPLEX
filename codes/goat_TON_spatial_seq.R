#Load packages
library(Seurat)
library(ClusterGVis)
library(org.Hs.eg.db)
library(ggpubr)

# Load10X_Spatial (for expr matrix) Function
Load10X_Spatial_ver2<-function (data.dir, filename = "filtered_feature_bc_matrix",
                                assay = "Spatial", slice = "slice1", filter.matrix = TRUE,
                                to.upper = FALSE, image = NULL, ...)
{
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'",
            immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  data <- Read10X(data.dir = file.path(data.dir, filename),
                  ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  if (is.null(x = image)) {
    image <- Read10X_Image(image.dir = file.path(data.dir,
                                                 "spatial"), filter.matrix = filter.matrix)
  }
  else {
    if (!inherits(x = image, what = "VisiumV1"))
      stop("Image must be an object of class 'VisiumV1'.")
  }
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}

# Import data
spe.2os3 <- Load10X_Spatial_ver2(
  data.dir = "~/r4.3_docker_rstudio/data_analysis/SBC/spatial-Seq/data_original/2OS3/",
  filename = "filtered_feature_bc_matrix",
  assay = "Spatial",
  slice = "2OS3"
)
spe.3os3 <- Load10X_Spatial_ver2(
  data.dir = "~/r4.3_docker_rstudio/data_analysis/SBC/spatial-Seq/data_original/3OS3/",
  filename = "filtered_feature_bc_matrix",
  assay = "Spatial",
  slice = "3OS3"
)
spe.3od3 <- Load10X_Spatial_ver2(
  data.dir = "~/r4.3_docker_rstudio/data_analysis/SBC/spatial-Seq/data_original/3OD3/",
  filename = "filtered_feature_bc_matrix",
  assay = "Spatial",
  slice = "3OD3"
)

# change orig.ident
spe.2os3<-RenameIdents(spe.2os3,"SeuratProject"="2OS3")
spe.3os3<-RenameIdents(spe.3os3,"SeuratProject"="3OS3")
spe.3od3<-RenameIdents(spe.3od3,"SeuratProject"="3OD3")
spe.2os3$orig.ident<-Idents(spe.2os3)
spe.3os3$orig.ident<-Idents(spe.3os3)
spe.3od3$orig.ident<-Idents(spe.3od3)

# remove region of optic chiasm
optic_chiasm_region <- readRDS("goat_TON_spatial_seq/optic_chiasm_region.rds")
spe.3od3 <- subset(spe.3od3,cells=c(optic_chiasm_region),invert=T)

# merge data
spatial_data <- merge(spe.2os3, y = list(spe.3os3, spe.3od3), add.cell.ids = c("2OS3", "3OS3", "3OD3"))

# NormalizeData
spatial_data <- NormalizeData(spatial_data, normalization.method = "LogNormalize", scale.factor = 10000)
# FindVariableFeatures
spatial_data <- FindVariableFeatures(spatial_data, selection.method = "vst", nfeatures = 2000)
# ScaleData
spatial_data <- ScaleData(spatial_data, features = VariableFeatures(spatial_data))
# basic processing
spatial_data <- RunPCA(spatial_data)
spatial_data <- FindNeighbors(spatial_data, dims = 1:20,reduction = "pca")
spatial_data <- FindClusters(spatial_data,resolution = 0.1, algorithm = 4, method = "igraph")
spatial_data <- RunUMAP(spatial_data, dims = 1:20)

# find DEG
spatial_data <- JoinLayers(spatial_data)
Idents(spatial_data) <- "seurat_clusters"
all_markers<-FindAllMarkers(spatial_data,logfc.threshold = 0.25,only.pos = T,min.pct = 0.25,verbose = T)
all_markers <- all_markers %>% arrange(desc(avg_log2FC))
saveRDS(all_markers, file = "goat_TON_spatial_seq/Markers/spatial_res0.1_Markers.rds")
openxlsx::write.xlsx(all_markers,file = "goat_TON_spatial_seq/Markers/spatial_res0.1_Markers.xlsx",rowNames=T)

# annotation for clusters
Idents(spatial_data) <- "seurat_clusters"
spatial_data@meta.data$region_class <- "Other"

spatial_data@meta.data[WhichCells(spatial_data, idents = c(1)),]$region_class = "Healthy region"
spatial_data@meta.data[WhichCells(spatial_data, idents = c(2)),]$region_class = "Meningeal region"
spatial_data@meta.data[WhichCells(spatial_data, idents = c(3)),]$region_class = "Neuroinflammatory region"

# plot for regions
my_cols <- c("Healthy region" = "#0dbf2f", "Meningeal region" = "#5088da", "Neuroinflammatory region" = "#ed6954")
SpatialDimPlot(spatial_data,group.by = "region_class",label = F, cols = my_cols,pt.size.factor = 2.5)
ggsave(filename = "goat_TON_spatial_seq/plots/umapplot/spatialdimplot_region_class.pdf",height = 5,width = 15)
DimPlot(spatial_data,label = T,cols = my_cols,group.by = "region_class")
ggsave(filename = "goat_TON_spatial_seq/plots/umapplot/dimplot_region_class.pdf",height = 5,width = 6)

# heatmap & enrichment
# find markers for every cluster compared to all remaining cells
# report only the positive ones
Idents(spatial_data) <- "region_class"
spatial.markers.all <- Seurat::FindAllMarkers(spatial_data,
                                              only.pos = TRUE,
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25)

spatial.markers.all <- spatial.markers.all %>% arrange(desc(avg_log2FC))
saveRDS(spatial.markers.all, file = "goat_TON_spatial_seq/Markers/Markers_cell_subclass.rds")
openxlsx::write.xlsx(spatial.markers.all,file = "goat_TON_spatial_seq/Markers/Markers_cell_subclass.xlsx",rowNames=F)

# no average cells
spatial.markers1 <- spatial.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

# retain duplicate diff gene in multiple clusters
st.data <- prepareDataFromscRNA(object = spatial_data,
                                diffData = spatial.markers1,
                                showAverage = FALSE)

# check
str(st.data)

# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

# check
head(enrich)

# add gene name
markGenes = unique(spatial.markers$gene)[sample(1:length(unique(spatial.markers$gene)),15,
                                                replace = F)]

# line plot
visCluster(object = st.data,
           plot.type = "line")

pdf('goat_TON_spatial_seq/plots/visCluster/sc1_all.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           sample.col = my_cols,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:9))
dev.off()

pdf('goat_TON_spatial_seq/plots/visCluster/sc2_all_top10genes.pdf',height = 12,width = 15.2,onefile = F)
visCluster(object = st.data,
           sample.col = my_cols,
           plot.type = "both",
           column_names_rot = 25,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:9),
           go.col = rep(jjAnno::useMyCol("stallion",n = 3),each = 5),
           add.bar = T)
dev.off()

# skull or blood score
skull.genes <- c("IGKC","GM36486","PIM1","GM10076","MFGE8","EAR6","LARS2","LRP1","MT-ND3",
                 "FN1","MT-ND4L","SLC6A6","ATP13A2","C77080","MAN2B1","STAB1","PLD3","GPC1","MOBP","CTSD")
blood.genes <- c("NGP","S100A9","NRGN","CSF2","IFNG","TNFSF11","NKG7","RGS16","IL17A","LCK","GIMAP5",
                 "GZMA","CD8B1","CD4","PTPRCAP","CD3E","LAT","CTSW","TRAC","CD3G")
# skull_derived
exprMatrix <- as.data.frame(spatial_data@assays$Spatial$data)
exprMatrix_seclet <- exprMatrix[rownames(exprMatrix) %in% skull.genes,]
avg_row <- colMeans(exprMatrix_seclet)
exprMatrix_seclet <- rbind(exprMatrix_seclet, avg_row)
exprMatrix_seclet <- t(exprMatrix_seclet[14,])
spatial_data <- AddMetaData(spatial_data, metadata=exprMatrix_seclet, col.name="Avg_skull_derived_gene_Expr")

# blood_derived
exprMatrix <- as.data.frame(spatial_data@assays$Spatial$data)
exprMatrix_seclet <- exprMatrix[rownames(exprMatrix) %in% blood.genes,]
avg_row <- colMeans(exprMatrix_seclet)
exprMatrix_seclet <- rbind(exprMatrix_seclet, avg_row)
exprMatrix_seclet <- t(exprMatrix_seclet[18,])
spatial_data <- AddMetaData(spatial_data, metadata=exprMatrix_seclet, col.name="Avg_blood_derived_gene_Expr")

SpatialFeaturePlot(spatial_data,features = c("Avg_skull_derived_gene_Expr","Avg_blood_derived_gene_Expr"),keep.scale = "all",pt.size.factor = 2.5) # & ggplot2::scale_fill_gradient2(limits = c(0.0, 1.0), breaks = c(0.0, 0.5, 1.0), low = "#FEE8C8", mid = "#FDBB84", high = "#E34A33", midpoint = 0.5)#目前在这里，准备换颜色，突然发现pt.size.factor改了会好看点，前面也可以改改
ggsave(width = 14.7,height = 10.5,filename = "goat_TON_spatial_seq/plots/score plot/spatial_skull & blood_derived_score.pdf")

# vlnplot
spatial_data$region_class_simp <- "others"
spatial_data@meta.data[spatial_data$region_class=="Healthy region","region_class_simp"] = "Parenchyma"
spatial_data@meta.data[spatial_data$region_class=="Neuroinflammatory region","region_class_simp"] = "Parenchyma"
spatial_data@meta.data[spatial_data$region_class=="Meningeal region","region_class_simp"] = "Meninge"

meta_df <- spatial_data@meta.data

macrophage_df <- meta_df %>%
  dplyr::select(Avg_skull_derived_gene_Expr, Avg_blood_derived_gene_Expr, region_class_simp)

macrophage_long <- macrophage_df %>%
  mutate(cell = rownames(.)) %>%
  tidyr::pivot_longer(cols = c("Avg_skull_derived_gene_Expr","Avg_blood_derived_gene_Expr"),
                      names_to = "Score_Type",
                      values_to = "Score")

p <- ggplot(macrophage_long, aes(x = Score_Type, y = Score, fill = Score_Type)) +
  geom_violin(trim = FALSE, scale = "width", adjust = 1.2) +
  geom_boxplot(width = 0.2, fill = "white", outlier.size = 0) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    comparisons = list(c("Avg_skull_derived_gene_Expr","Avg_blood_derived_gene_Expr")),
    label.y = max(macrophage_long$Score, na.rm = TRUE) * 1.2,
    tip.length = 0,
    size = 10
  ) +
  facet_wrap(~region_class_simp, ncol = 2) +  # 按样本分面，ncol 可改
  scale_fill_manual(values = c("#fcae59", "#db498e")) +
  labs(x = "", y = "Score") +
  theme_classic2(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold")  # 分面标题
  )
p
ggsave(filename = "goat_TON_spatial_seq/plots/score plot/EAE_Skull_vs_Blood_vlnplot_region_class_simp.pdf",height =4 ,width = 6)





