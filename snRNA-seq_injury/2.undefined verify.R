# load packages
library(Seurat)

# load RDS data
seurat.obj.combined <- readRDS(file = "snRNA-seq_injury/results/SCT.CCA/1.data input & merge/goat_ON injury snRNA_SCT.CCA after delete outlier.rds")

# create folder
dir.create("snRNA-seq_injury/results/SCT.CCA/2.undefined verify",recursive = T)

# find undefined subclusters
idents(seurat.obj.combined) <- "celltype_undefined"
seurat.obj.combined <- FindNeighbors(seurat.obj.combined,graph.name = "leiden")
seurat.obj.combined <- FindSubCluster(object = seurat.obj.combined,cluster = "Undefined",subcluster.name = "undefined_subtype",graph.name = "leiden",resolution = 0.1,algorithm = 4)
Dimplot(seurat.obj.combined, group.by= "undefined_subtype")
ggsave("snRNA-seq_injury/results/SCT.CCA/2.undefined verify/undefined_subtype_dimplot_res0.1.pdf",width = 7,height = 5)

# cell cycle scoring
seurat.obj.combined <- CellCycleScoring(seurat.obj.combined,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
VlnPlot(seurat.obj.combined,features = c("S.Score","G2M.Score"),group.by = "undefined_subtype")&xlab("")
ggsave(filename = "snRNA-seq_injury/results/SCT.CCA/2.undefined verify/all.celltype_undefined_split_cc.score.pdf",width = 10,height = 5)

# dotplot
undefined <- subset(seurat.obj.combined,celltype_undefined=="undefined")
DotPlot(undefined,features = c("GFAP","P2RY12","COL3A1","PCDH15"))+RotatedAxis()+xlab("")+ylab("")
ggsave(filename = "snRNA-seq_injury/results/SCT.CCA/2.undefined verify/signature gene dotplot.pdf",width = 5.5,height = 3.5)

# undefined cluster renameidents
Idents(seurat.obj.combined) <- "undefined_subtype"
seurat.obj.combined<-RenameIdents(seurat.obj.combined,
                                    "Oligodendrocyte"="Oligodendrocyte",
                                    "Astrocyte"="Astrocyte",
                                    "Fibroblast"="Fibroblast",
                                    "OPC"="OPC",
                                    "VEC"="VEC",
                                    "MP"="MP",
                                    "Vascular Mural Cell"="Vascular Mural Cell",
                                    "Undefined_0"="Astrocyte",
                                    "Undefined_1"="Astrocyte",
                                    "Undefined_2"="OPC",
                                    "Undefined_3"="MP",
                                    "Undefined_4"="Fibroblast")
seurat.obj.combined$celltype <- Idents(seurat.obj.combined)

# save RDS file
saveRDS(seurat.obj.combined, file = "snRNA-seq_injury/results/SCT.CCA/2.undefined verify/goat_ON injury snRNA_SCT.CCA undefined_verify.rds")
