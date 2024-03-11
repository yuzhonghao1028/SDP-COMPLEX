library(Seurat)
library(ggplot2)
library(future)
library(dplyr)

# snRNA-seq data input 
control1 <- Read10X("2.1 Exp Matrix/1OD3/filtered_feature_bc_matrix/")
control2 <- Read10X("2.1 Exp Matrix/3OD3/filtered_feature_bc_matrix/")
injury1 <- Read10X("2.1 Exp Matrix/2OS3/filtered_feature_bc_matrix/")
injury2 <- Read10X("2.1 Exp Matrix/3OS3/filtered_feature_bc_matrix/")
seurat.list<-list(control1=control1,control2=control2,injury1=injury1,injury2=injury2)


# crate seurat object
for (i in 1:length(seurat.list)) {
  seurat.list[[i]] <- CreateSeuratObject(seurat.list[[i]],
                                         project = names(seurat.list)[i], 
                                         min.cells = 3,
                                         min.features = 200,
  )
}


# quality control
dir.create("result/QC",recursive = T)
## QC
for (i in 1:length(seurat.list)) {
  VlnPlot(seurat.list[[i]], features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,pt.size = 0)
  ggsave(filename = paste0("result/QC/QC-VlnPlot_",names(seurat.list)[i],".pdf"),width = 6,height = 4.5)
  a <- seurat.list[[i]]@meta.data
  ggplot(data = a,aes(x = a$nCount_RNA))+geom_density(fill="pink")+xlim(c(0,30000))
  ggsave(filename = paste0("result/QC/QC-nCount.density_",names(seurat.list)[i],".pdf"),width = 6,height = 4.5)
  ggplot(data = a,aes(x = a$nFeature_RNA))+geom_density(fill="pink")+xlim(c(0,10000))
  ggsave(filename = paste0("result/QC/QC-nFeature.density_",names(seurat.list)[i],".pdf"),width = 6,height = 4.5)
}


# filter threshold: Feature_RNA > 300 & nFeature_RNA < quantile(seurat_sample, 0.975)
## control1
cat(names(seurat.list)[1],": Before filter ",nrow(seurat.list[[1]]@meta.data),"cells\n") #Before filter : 8940 cells
seurat.list[[1]] <- subset(seurat.list[[1]], 
                           subset = 
                             nFeature_RNA > 300 & 
                             nFeature_RNA < 4234)
cat(names(seurat.list)[1],": After filter ",nrow(seurat.list[[1]]@meta.data),"cells\n") #Before filter : 8714 cells

## control2
cat(names(seurat.list)[2],": Before filter ",nrow(seurat.list[[2]]@meta.data),"cells\n") #Before filter : 6877 cells
seurat.list[[2]] <- subset(seurat.list[[2]], 
                           subset = 
                             nFeature_RNA > 300 & 
                             nFeature_RNA < 5091)
cat(names(seurat.list)[2],": After filter ",nrow(seurat.list[[2]]@meta.data),"cells\n") #Before filter : 6736 cells

## injury1
cat(names(seurat.list)[3],": Before filter ",nrow(seurat.list[[3]]@meta.data),"cells\n") #Before filter : 9427 cells
seurat.list[[3]] <- subset(seurat.list[[3]], 
                           subset = 
                             nFeature_RNA > 300 & 
                             nFeature_RNA < 4778)
cat(names(seurat.list)[3],": After filter ",nrow(seurat.list[[3]]@meta.data),"cells\n") #Before filter : 9106 cells

## injury2
cat(names(seurat.list)[4],": Before filter ",nrow(seurat.list[[4]]@meta.data),"cells\n") #Before filter : 11987 cells
seurat.list[[4]] <- subset(seurat.list[[4]], 
                           subset = 
                             nFeature_RNA > 300 & 
                             nFeature_RNA < 4927)
cat(names(seurat.list)[4],": After filter ",nrow(seurat.list[[4]]@meta.data),"cells\n") #Before filter : 11684 cells


# SCTransform
seurat.list <- lapply(X = seurat.list, FUN = SCTransform )
features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = features)
seurat.obj.anchors <- FindIntegrationAnchors(object.list = seurat.list, normalization.method = "SCT", anchor.features = features)
seurat.obj.combined <- IntegrateData(anchorset = seurat.obj.anchors, normalization.method = "SCT")


# PCA
seurat.obj.combined <- RunPCA(seurat.obj.combined, verbose = FALSE)


# dir.create
dir.create("result/SCT.CCA")


# ElbowPlot
pdf(paste0("result/SCT.CCA/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(seurat.obj.combined,ndims = 50)
dev.off()


# choose dimention
dim.use <- 1:30


# clustering
seurat.obj.combined <- FindNeighbors(seurat.obj.combined, dims = dim.use,reduction = "pca")
seurat.obj.combined <- FindClusters(seurat.obj.combined, resolution = 0.2)


# UMAP reduction
seurat.obj.combined <- RunUMAP(seurat.obj.combined, reduction = "pca", dims = dim.use)


# save RDS file
saveRDS(seurat.obj.combined,file = "result/SCT.CCA/goat.ON.snRNA_SCT.CCA.rds")
