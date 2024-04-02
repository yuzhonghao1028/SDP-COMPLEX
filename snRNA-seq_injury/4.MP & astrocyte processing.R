# load packages
library(Seurat)
library(AUCell)
library(viridis)
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(edgeR)
library(zinbwave)
library(scran)


# extract seurat color function
extract.col<-function(picture){
  p1 <- picture 
  x<-ggplot_build(p1)
  info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
  info <- unique((arrange(info, group)))
  seurat.cols <<- as.character(info$colour)
  rm(p1,x,info)
}


# load RDS data
seurat.obj.combined <- readRDS(file = "snRNA-seq_injury/results/SCT.CCA/2.undefined verify/goat_ON injury snRNA_SCT.CCA undefined_verify.rds")


# create folder
dir.create("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing",recursive = T)


# Isolation of monocular phagocytes populations
MP <- subset(seurat.obj.combined, celltype=="MP")


# heatmap of top20 blood-derived monocytes and skull derived monocytes signatures
## top 20 upregulated DEGs for blood-derived monocytes and skull derived monocytes from Andrea Cugurra's study
skull.genes <- c("PIM1",
                 "MFGE8",
                 "LARS2",
                 "LRP1",
                 "FN1",
                 "SLC6A6",
                 "LOC102183471",
                 "NHSL3",
                 "MAN2B1",
                 "LOC102172305",
                 "PLD3",
                 "GPC1",
                 "CTSD",
                 "NAGLU",
                 "LOC102184794",
                 "PSAP")
blood.genes <- c("CD3G",
                 "CTSW",
                 "LAT",
                 "CD3E",
                 "PTPRCAP",
                 "CD4",
                 "GZMA",
                 "LOC102184497",
                 "LCK",
                 "IL17A",
                 "RGS16",
                 "NKG7",
                 "TNFSF11",
                 "IFNG",
                 "CSF2",
                 "NRGN",
                 "S100A9")
seurat.obj.combined$celltype_group <- paste0(seurat.obj.combined$celltype,"_",seurat.obj.combined$group)
Idents(seurat.obj.combined)<-"celltype_group"
selected.genes<-c(skull.genes,blood.genes)
hm_average <- as.data.frame(AverageExpression(seurat.obj.combined, assays = "RNA",slot = "scale.data"))
hm_average_selective <- hm_average %>% filter(row.names(hm_average) %in% selected.genes)
hm_average_selective <- hm_average_selective[,c("RNA.MP_injury","RNA.MP_control")]
hm_average_selective$row_names <- row.names(hm_average_selective)
hm_average_selective <- hm_average_selective[order(match(hm_average_selective$row_names, selected.genes)), ]
hm_average_selective$row_names <- NULL
colnames(hm_average_selective)<-c("Ctrl.","Inj.")
pdf(file = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/EAE_Skull&blood gene heatmap.pdf",width = 1.9,height = 3)
Heatmap((hm_average_selective), row_names_gp = gpar(fontsize = 6), col = brewer.pal(9,name = "Purples"),
        cluster_rows = F,cluster_columns = F,name = "Scaled Exp.",column_names_gp = gpar(fontsize = 10))
dev.off()


# AUCell scoring
exprMatrix <- GetAssayData(MP, slot='counts')
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=2, plotStats=TRUE)
## skull derived score
cells_AUC <- AUCell_calcAUC(skull.genes, cells_rankings, aucMaxRank = 500, verbose = TRUE)
data <- t(cells_AUC@assays@data$AUC)
MP <- AddMetaData(MP, metadata=data, col.name="AUCell_skull_derived_score")
## blood derived score
cells_AUC <- AUCell_calcAUC(blood.genes, cells_rankings, aucMaxRank = 500, verbose = TRUE)
data <- t(cells_AUC@assays@data$AUC)
MP <- AddMetaData(MP, metadata=data, col.name="AUCell_blood_derived_score")
## featureplot
FeaturePlot(MP,features = "AUCell_skull_derived_score",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/AUCell_Skull_derived_score_top20.pdf",width = 7.4,height = 2.8)
FeaturePlot(MP,features = "AUCell_blood_derived_score",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/AUCell_blood_derived_score_top20.pdf",width = 7.4,height = 2.8)


# AddModuleScore for MP
## DAM score
pal <- viridis(n = 10, option = "H", direction = 1)
MP <- AddModuleScore(MP,features =list(c("CD9",
                                         "APOE",
                                         "TREM2",
                                         "TYROBP",
                                         "CD63",
                                         "CD68",
                                         "AXL",
                                         "SPP1",
                                         "CTSB",
                                         "CTSD",
                                         "LPL",
                                         "ITGAX",
                                         "CST7",
                                         "CSF1",
                                         "GPNMB",
                                         "IRF8",
                                         "FTH1")),name = "DAM_score" )
FeaturePlot(MP,features = "DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/DAM.score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## Amyloid DAM score
MP <- AddModuleScore(MP,features = list(c("MYO1E",
                                          "ANKH",
                                          "CACNA1A",
                                          "FGF13",
                                          "CTNNA3",
                                          "TMEM163",
                                          "MAMDC2",
                                          "XYLT1",
                                          "APBB2",
                                          "KCNMA1",
                                          "PTCHD1",
                                          "ACACA",
                                          "HIF1A",
                                          "ARHGAP24",
                                          "PRKG1",
                                          "CSF1")),name = "Amyloid_DAM_score")
FeaturePlot(MP,features = "Amyloid_DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/Amyloid_DAM_score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "Amyloid_DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## Myelin_DAM_score
MP <- AddModuleScore(MP,features = list(c("GPNMB",
                                          "APOBEC1",
                                          "ATP6V0D2",
                                          "FLT1",
                                          "FGR",
                                          "MAMDC2",
                                          "DDHD1",
                                          "CERS6",
                                          "BCAR3",
                                          "MDFIC",
                                          "TRPC4",
                                          "MARCH3",
                                          "FMN1",
                                          "MYOF",
                                          "FAM20C",
                                          "NRP2",
                                          "CSF2RA",
                                          "COLEC12")),name = "Myelin_DAM_score")
FeaturePlot(MP,features = "Myelin_DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/Myelin_DAM_score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "Myelin_DAM_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## top10 heatmap for Amyloid and Myelin_DAM_score
Idents(MP) <- "group"
extract.col(DimPlot(MP))
DoHeatmap(MP,features =c("MYO1E",
                         "ANK",
                         "CACNA1A",
                         "FGF13",
                         "ETL4",
                         "CTNNA3",
                         "GM15155",
                         "TMEM163",
                         "MAMDC2",
                         "XYLT1",
                         "APBB2",
                         "KCNMA1",
                         "PTCHD1",
                         
                         "GPNMB",
                         "APOBEC1",
                         "ATP6V0D2",
                         "FLT1",
                         "FGR",
                         "MAMDC2",
                         "DDHD1",
                         "CERS6",
                         "BCAR3",
                         "MDFIC",
                         "CYBB",
                         "TRPC4",
                         "LILRB4A",
                         "MARCH3") ,label = F,slot = "scale.data",group.colors =c(seurat.cols[2],seurat.cols[1]) )+scale_fill_gradientn(colors = c("lightskyblue","white","firebrick3"))
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/amyloid_myelin_dam_score_top10 gene heatmap.png",height = 3.5,width = 3)

## homeostatic score
MP <- AddModuleScore(MP,features = list(c("SELPLG",
                                          "ITGA9",
                                          "LOC102185671",
                                          "RAPGEF5",
                                          "CDH23",
                                          "ZFHX3",
                                          "FGD2",
                                          "FCHSD2",
                                          "PDE3B",
                                          "MYO18B",
                                          "PLXNA4",
                                          "NAV3",
                                          "PHYHD1",
                                          "SRGAP2",
                                          "COL27A1")),name = "Homeostatic_score")
FeaturePlot(MP,features = "Homeostatic_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/homeotatic.score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "Homeostatic_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## MHC II score
GOID <- c("GO:0019886")
GOgeneID <- get(GOID, org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
MP <- AddModuleScore(MP,features =list(GOgeneID),name = "MHC_score")
FeaturePlot(MP,features = "MHC_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MHC_score featureplot(GO:0019886).pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "MHC_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## IFN score
MP <- AddModuleScore(MP,features =list(c("LOC102180655",
                                         "IFIT3",
                                         "CXCL10",
                                         "STAT1",
                                         "IRF7",
                                         "ISG15")),name = "IFN_score" )
FeaturePlot(MP,features = "IFN_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/IFN_score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "IFN_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## Proliferation score
MP <- AddModuleScore(MP,features =list(c("TOP2A",
                                         "MKI67",
                                         "CENPE",
                                         "MCM5",
                                         "BIRC5",
                                         "LOC102183306")),name = "Proliferation_score" )
FeaturePlot(MP,features = "Proliferation_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/Proliferation_score featureplot.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(MP,features = "Proliferation_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)

## Inflammatory response score
GOID <- c("GO:0006954")
GOgeneID <- get(GOID, org.Hs.egGO2ALLEGS) %>% mget(org.Hs.egSYMBOL) %>% unlist()
MP <- AddModuleScore(MP,features =list(GOgeneID),name = "Inflammatory_response_score",assay = "RNA")
FeaturePlot(MP,features = "Inflammatory_response_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
ggsave("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/inflammatory_response_score splitby group featureplot.pdf",height = 3,width = 8)
Idents(MP) <- "group"
VlnPlot(MP,features = "Inflammatory_response_score1",pt.size = 0,cols = c(seurat.cols[2],seurat.cols[1]))+xlab("")+
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
ggsave("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/inflammatory_response_score splitby group vlnplot with box.pdf",height = 3.5,width = 2.5)


# blood derived score vs. skull derived score vlnplot
data <- read.csv("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/EAE skull &amp; blood score split by group.csv")
ggplot(data, aes(x=group, y=score, fill=source)) +
  geom_violin(scale="width",
              adjust = 1, 
              trim = F )+ 
  labs(x = 'Cluster', y = "Expression Level", fill = NULL)+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ 
  theme(panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  theme_classic()+ 
  theme(text = element_text(size = 15))+NoLegend()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/AUCell_EAE_blood&skull derived score splitby_group vlnplot.pdf",width = 4,height = 4)  


# MP subcluster
MP <- FindNeighbors(MP, dims = dim.use,reduction = "pca",graph.name = "subcluster")
MP <- FindClusters(MP, resolution = 0.2, algorithm = 4, method = "igraph", graph.name = "subcluster",random.seed = 123)
DimPlot(MP,group.by = "subcluster_res.0.2",split.by= "group")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MP subcluster res_0.2 splitby_group dimplot.pdf",width = 10.5,height = 5)  
# extract seurat colors
extract.col(Dimplot(MP,group.by="subcluster_res.0.2"))


# MP subcluster propotion stack barplot
data <- as.data.frame(table(MP$subcluster_res.0.2,MP$group))
data <- data %>% group_by(Var2) %>% mutate(percentage = Freq/sum(Freq)) %>% mutate(label = Freq/sum(Freq)*100)
data$label <- paste(sprintf("%.1f", data$label),"%")
extract.col(DimPlot(MP,group.by = "subcluster_res.0.2"))
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
    legend.text=element_text(size=12)
    #legend.position = ' none' 
    #legend.position="bottom" ,
    #legend.background = element_rect(fill="lightblue",size=0.5, linetype="solid", colour ="darkblue")
  )+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
  )+ 
  guides(fill=guide_legend(keywidth = 1, keyheight = 1))+ 
  ylim(c(-0.005,1.005))+
  scale_y_continuous(labels = scales::percent_format(scale = 100))+
  geom_text(aes(label=label),size=4,color="black",position = "stack")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/subcluster propotion stack barplot.pdf",width = 5,height = 6)


# MP subcluster addmodelescore
## microglia score vlnplot
MP <- AddModuleScore(MP,features = list(c("P2RY12","TMEM119","HEXB","SALL1")),name = "microglia_score")
VlnPlot(MP,features = "microglia_score1",pt.size = 0,split.by = "subcluster_res.0.2",group.by = "subcluster_res.0.2",cols = seurat.cols)+
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white")+NoLegend()+labs(y = "Expression Level")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/microglia_score splitby_group vlnplot.pdf",height = 3.45 ,width = 3)
## skull derived score vlnplot
MP <- AddModuleScore(MP,features = list(skull.genes),name = "EAE_Skull_derived_score")
VlnPlot(MP,features = "EAE_Skull_derived_score1",pt.size = 0,split.by = "subcluster_res.0.2",group.by = "subcluster_res.0.2",cols = seurat.cols)+
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white")+NoLegend()+labs(y = "Expression Level")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/EAE_Skull_derived_score splitby_group vlnplot.pdf",height = 3.3,width = 3)

VlnPlot(MP,features = "AUCell_skull_derived_score",pt.size = 0,split.by = "subcluster_res.0.2",group.by = "subcluster_res.0.2",cols = seurat.cols)+
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill="white")+NoLegend()+labs(y = "Expression Level")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/AUCell_skull_derived_score splitby_group vlnplot.pdf",height = 3.3,width = 3)


# MP DEGs analysis (ZINBWAVE method)
## renamed cluster
Idents(MP) <- "subcluster_res.0.2"
MP <- RenameIdents(MP,
                    "1"="1",
                    "2"="others",
                    "3"="others",
                    "4"="others",
                    "5"="others")
MP$cluster1 <- Idents(MP)

Idents(MP) <- "subcluster_res.0.2"
MP <- RenameIdents(MP,
                    "1"="others",
                    "2"="2",
                    "3"="others",
                    "4"="others",
                    "5"="others")
MP$cluster2 <- Idents(MP)

Idents(MP) <- "subcluster_res.0.2"
MP <- RenameIdents(MP,
                    "1"="others",
                    "2"="others",
                    "3"="3",
                    "4"="others",
                    "5"="others")
MP$cluster3 <- Idents(MP)

Idents(MP) <- "subcluster_res.0.2"
MP <- RenameIdents(MP,
                    "1"="others",
                    "2"="others",
                    "3"="others",
                    "4"="4",
                    "5"="others")
MP$cluster4 <- Idents(MP)

Idents(MP) <- "subcluster_res.0.2"
MP <- RenameIdents(MP,
                    "1"="others",
                    "2"="others",
                    "3"="others",
                    "4"="others",
                    "5"="5")
MP$cluster5 <- Idents(MP)

## Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(MP, assay="RNA")
## Filter out lowly expressed genes
## Identify genes that have a count of at least 3 in at least 6 cells
filter <- rowSums(assay(sce)>2)>5
## Generate information on how many genes meet these criteria
table(filter)
## Filter matrix based on genes above
sce.filt <- sce[filter,]
## View matrix stats after filtering
sce.filt
rm(sce,filter)
## Select highly variable genes on which to focus analysis
## Model gene variance
sce.var <- modelGeneVar(sce.filt)
## Using gene variance, identify the top 2000
keep <- getTopHVGs(sce.var, n=2000)
## Filter matrix for only these genes
sce.filt <- sce.filt[keep,]
rm(sce.var,keep)
## Generate observational weights for genes
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
## Isolate weights to pass to DGE calculations
weights <- assay(sce.zinb, "weights")

## subcluster 1 DEGs
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$cluster1) 
design <- model.matrix(~0+condition)
## Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
## Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
## Add observational weights into the edgeR object
dge$weights <- weights
## Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
## Fit model
fit <- glmFit(dge, design)
## Make comparison between cluster1 and others 
colnames(design) <- make.names(colnames(design))
contrast <- makeContrasts(condition1-conditionothers, levels=design) 
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_DGE_result_MP_subcluster1.rds') 
## Summarize results
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster1 0.25lfc.xlsx') 

## subcluster 2 DEGs
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$cluster2) 
design <- model.matrix(~0+condition)
## Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
## Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
## Add observational weights into the edgeR object
dge$weights <- weights
## Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
## Fit model
fit <- glmFit(dge, design)
## Make comparison between cluster2 and others 
colnames(design) <- make.names(colnames(design))
contrast <- makeContrasts(condition2-conditionothers, levels=design) 
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_DGE_result_MP_subcluster2.rds')  
## Summarize results
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster2 0.25lfc.xlsx') 

## subcluster 3 DEGs
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$cluster3) 
design <- model.matrix(~0+condition)
## Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
## Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
## Add observational weights into the edgeR object
dge$weights <- weights
## Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
## Fit model
fit <- glmFit(dge, design)
## Make comparison between cluster3 and others 
colnames(design) <- make.names(colnames(design))
contrast <- makeContrasts(condition3-conditionothers, levels=design) 
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_DGE_result_MP_subcluster3.rds') 
## Summarize results
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster3 0.25lfc.xlsx') 

## subcluster 4 DEGs
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$cluster4) 
design <- model.matrix(~0+condition)
## Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
## Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
## Add observational weights into the edgeR object
dge$weights <- weights
## Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
## Fit model
fit <- glmFit(dge, design)
## Make comparison between cluster4 and others 
colnames(design) <- make.names(colnames(design))
contrast <- makeContrasts(condition4-conditionothers, levels=design) 
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_DGE_result_MP_subcluster4.rds')  
## Summarize results
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster4 0.25lfc.xlsx')

## subcluster 5 DEGs
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$cluster5) 
design <- model.matrix(~0+condition)
## Create edgeR object from model with weights
dge <- DGEList(assay(sce.zinb))
## Caluclate normalization factors with TMM
dge <- calcNormFactors(dge)
## Add observational weights into the edgeR object
dge$weights <- weights
## Estimate dispersion, passing defined design matrix to function
dge <- estimateDisp(dge, design)
## Fit model
fit <- glmFit(dge, design)
## Make comparison between cluster5 and others 
colnames(design) <- make.names(colnames(design))
contrast <- makeContrasts(condition5-conditionothers, levels=design)
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_DGE_result_MP_subcluster5.rds') 
## Summarize results
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file='results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster5 0.25lfc.xlsx') 


# MP subcluster volcano plot
## subcluster1
data<-read.xlsx("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster1 0.25lfc.xlsx")
data <- data %>%
  mutate(gene_type = case_when(logFC >= 0.25 & Pval <= 0.05 ~ "up",
                               logFC <= -0.25 & Pval <= 0.05 ~ "down",
                               TRUE ~ "ns")) 
sig_il_genes <- data %>%
  filter(row.names.genes. %in% c("CD36", "ENPP1", "NYAP2","STAR","GPNMB",
                       "PDZRN4","BICD1","LOC102173852","GJA1","ANK3"))
up_il_genes <- data %>%
  filter(row.names.genes. %in% c("CD36", "ENPP1", "NYAP2","STAR","GPNMB"))
down_il_genes <- data %>%
  filter(row.names.genes. %in% c("PDZRN4","BICD1","LOC102173852","GJA1","ANK3"))
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
        sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
        alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)       
        ggplot(data = data,aes(x = logFC,y = -log10(Pval))) + 
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
          geom_vline(xintercept = c(-0.25, 0.25),
                     linetype = "dashed") +
          geom_label_repel(data = sig_il_genes,   
                           aes(label = row.names.genes.),
                           force = 2,
                           nudge_y = 1) +
          scale_colour_manual(values = cols) + 
          scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                             limits = c(-10, 10)) +
          labs(title = "Cluster1 DEGs in Inj. vs. Ctrl.",
               x = "log2(fold change)",
               y = "-log10(P-value)",
               colour = "Expression \nchange") +
          theme_bw() + # Select theme with a white background  
          theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) 
        ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MP subcluster1 DEGs volcano in injury vs control.pdf",height = 3,width = 4)

## subcluster2
data<-read.xlsx("snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/zinb_sig_genes MP_subcluster2 0.25lfc.xlsx")
data <- data %>%
mutate(gene_type = case_when(logFC >= 0.25 & Pval <= 0.05 ~ "up",
                                       logFC <= -0.25 & Pval <= 0.05 ~ "down",
                                       TRUE ~ "ns")) 
sig_il_genes <- data %>%
filter(row.names.genes. %in% c("SCIN", "CTSL", "PLXNA4","LOC106502447","KMO",
                                "CR2","LOC102169149","SKAP1","SPP1","ADCY8"))
up_il_genes <- data %>%
filter(row.names.genes. %in% c("SCIN", "CTSL", "PLXNA4","LOC106502447","KMO"))
down_il_genes <- data %>%
filter(row.names.genes. %in% c("CR2","LOC102169149","SKAP1","SPP1","ADCY8"))
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
ggplot(data = data,aes(x = logFC,y = -log10(Pval))) + 
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
geom_vline(xintercept = c(-0.25, 0.25),
                     linetype = "dashed") +
geom_label_repel(data = sig_il_genes,   
                           aes(label = row.names.genes.),
                           force = 2,
                           nudge_y = 1) +
scale_colour_manual(values = cols) + 
scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                             limits = c(-10, 10)) +
labs(title = "Cluster2 DEGs in Inj. vs. Ctrl.",
               x = "log2(fold change)",
               y = "-log10(P-value)",
               colour = "Expression \nchange") +
theme_bw() + # Select theme with a white background  
theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) 
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MP subcluster2 DEGs volcano in injury vs control.pdf",height = 3,width = 4)
        
        
# enrichment barplot
## subcluster 1
data<-read.csv("cluster1 gobp.csv")
font.size =12
data %>% 
ggplot(aes(x=forcats::fct_reorder(Term,P.value,.desc = T),y=-log10(P.value)))+ 
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
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MP subcluster1 DEGs GO BP barplot in injury vs control.pdf",width = 6.3,height = 2.4)
        
## subcluster 1
data<-read.csv("cluster2 gobp.csv")
font.size =12
data %>% 
ggplot(aes(x=forcats::fct_reorder(Term,P.value,.desc = T),y=-log10(P.value)))+ 
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
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/MP subcluster1 DEGs GO BP barplot in injury vs control.pdf",width = 8.0,height = 2.4)
        
        
# chemokine pairs daotplot
## ligand
DotPlot(seurat.obj.combined.injury,features = c("CXCL12",
                                                        "CX3CL1",
                                                        "LOC102182115",
                                                        "LOC102174969",
                                                        "LOC102180880",
                                                        "LOC102175889",
                                                        "CCL5",
                                                        "LOC102175436"),cols = c("grey", "#ff70b5"))+RotatedAxis()+coord_flip()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/ligand dotplot_only.injury.pdf",width = 5.5,height = 5.5)
        
## receptor
DotPlot(seurat.obj.combined.injury,features = c("CXCR4",
                                                        "CX3CR1",
                                                        "CXCR2",
                                                        "LOC102179993",
                                                        "LOC102178953",
                                                        "CCR5"),cols = c("grey", "#ff70b5"))+RotatedAxis()+coord_flip()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/receptor dotplot_only.injury.pdf",width = 5.5,height = 5.5)
        
        
# MP cell cycle score vlnplot
VlnPlot(MP,features = "S.Score",pt.size = 0)+ labs(y = "Expression Level")+NoLegend()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/vlnplot s.score for subcluster.pdf",height = 3,width = 3.5)
        
VlnPlot(MP,features = "G2M.Score",pt.size = 0)+ labs(y = "Expression Level")+NoLegend()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/vlnplot G2M.score for subcluster.pdf",height = 3,width = 3.2)
        
VlnPlot(MP,features = "MKI67",pt.size = 0)+ labs(y = "Expression Level")+NoLegend()
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/vlnplot MKI67 for subcluster.pdf",height = 3,width = 3.2)
        
        
# save MP RDS file
saveRDS(MP,file = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/goat_ON injury snRNA_SCT.CCA MP cluster.rds")
                
                
# Astrocyte
## Isolation of monocular phagocytes populations
astro <- subset(seurat.obj.combined, celltype=="Astrocyte")
                
## A1 score
astro<-AddModuleScore(astro,features =list(c("SERPING1",	"LIG1",	"LOC102170144",	"FBLN5",	
                                                "FBLN1",	"LOC102190288",	"FKBP5",	"PSMB8",	
                                                   "SRGN",	"AMIGO2")),name = "A1_score" )
FeaturePlot(astro,features = "A1_score1",split.by = "group",order = T,cols= pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/A1_score.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(astro,features = "A1_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
                
## A2 score
astro<-AddModuleScore(astro,features =list(c("CLCF1",	"TGM1",	"PTX3",	"S100A10",	"SPHK1",	"CD109",
                                               "PTGS2",	"EMP1",	"SLC10A6",	"TM4SF1",	"B3GNT5",	"CD14")),name = "A2_score" )
FeaturePlot(astro,features = "A2_score1",split.by = "group",order = T,cols= pal)& theme(legend.position = "right")
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/A2_score.pdf",width = 6.8,height = 2.5)
## calculate positive rate
p1 <- FeaturePlot(astro,features = "A2_score1",split.by = "group",order = T,cols = pal)& theme(legend.position = "right")
data <- as.data.frame(as.matrix(p1[[1]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
data <- as.data.frame(as.matrix(p1[[2]]$data))
data[,4] <- as.numeric(data[,4])
table(data[,4]>5)
                
# A1 & A2 heatmap
extract.col(DimPlot(astro,group.by = "group"))
DoHeatmap(astro,features =c("SERPING1",	"LIG1",	"LOC102170144",	"FBLN5",	"FBLN1",
                               "LOC102190288",	"FKBP5",	"PSMB8",	"SRGN",	"AMIGO2",
                                "CLCF1",	"TGM1",	"PTX3",	"S100A10",	"SPHK1",	"CD109",
                                 "PTGS2",	"EMP1",	"SLC10A6",	"TM4SF1",	"B3GNT5",	"CD14"),
label = F,group.by = "group" ,assay = "RNA",slot = "scale.data",size = 2,
group.colors = c(seurat.cols[2],seurat.cols[1]))+
scale_fill_gradientn(colors = c("lightskyblue","white","firebrick3"))
ggsave(filename = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/astro.A1&A2 gene heatmap.png",height = 5,width = 3.5)
                
# save Astrocyte RDS file
saveRDS(MP,file = "snRNA_seq_injury/results/SCT.CCA/4.MP & astrocyte processing/goat_ON injury snRNA_SCT.CCA Astrocyte cluster.rds")
