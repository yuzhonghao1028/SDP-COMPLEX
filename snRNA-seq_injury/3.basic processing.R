# load packages
library(Seurat)
library(ggplot2)
library(reshape2)
library(SingleCellExperiment)
library(edgeR)
library(zinbwave)
library(scran)

# load RDS data
seurat.obj.combined <- readRDS(file = "snRNA-seq_injury/results/SCT.CCA/undefined verify/goat_ON injury snRNA_SCT.CCA undefined_verify.rds")

# create folder
dir.create("snRNA-seq_injury/results/SCT.CCA/basic processing",recursive = T)

# dimplot for celltype
seurat.col <- pal_locuszoom("default", alpha = 1)(7)
seurat.obj.combined$celltype <- factor(seurat.obj.combined$celltype,levels = c("Astrocyte","MP","OPC","Fibroblast","Oligodendrocyte","VEC","Vascular Mural Cell"))
DimPlot(seurat.obj.combined,group.by = "celltype",split.by = "group",cols = seurat.col)+NoAxes()+labs(title = "")
ggsave("snRNA-seq_injury/results/SCT.CCA/basic processing/celltype.splitby.group.umap.pdf",width = 12,height = 5)

DimPlot(seurat.obj.combined,group.by = "celltype",cols = seurat.col)+NoAxes()+labs(title = "")
ggsave("result/SCT.CCA/basic processing/celltype.umap.pdf",width = 7.5,height = 5)
DimPlot(seurat.obj.combined,group.by = "celltype",split.by = "group",cols = seurat.col)+NoAxes()+labs(title = "")
ggsave("snRNA-seq_injury/results/SCT.CCA/basic processing/celltype.splitby.group.umap.pdf",width = 12,height = 5)

# findallmarkers
DEGs <- FindAllMarkers(seurat.obj.combined,logfc.threshold = 0.25,min.pct = 0.25,verbose = T)
write.table(deg1,file = "snRNA-seq_injury/results/SCT.CCA/basic processing/celltype_deg.xls",sep = "\t")

# dotplot for marker genes
DotPlot(seurat.obj.combined,features = c("MOBP","MBP","GFAP","ALDH1L1","FN1","COL3A1","LHFPL3","PCDH15",
                                      "ERG","FLT1","P2RY12","CD36","CALD1","ACTA2"), cols = c("grey", "#ff70b5"))+RotatedAxis()+ggplot2:::coord_flip()+ylab("")+xlab("")
ggsave(filename = "snRNA-seq_injury/results/SCT.CCA/basic processing/celltype_marker_dotplot.pdf",width = 6,height = 5.5)

# celltype propotion stackbarplot
data <- as.data.frame(table(seurat.obj.combined$celltype,seurat.obj.combined$group))
data <- data %>% group_by(Var2) %>% mutate(percentage = Freq/sum(Freq)) %>% mutate(label = Freq/sum(Freq)*100)
data$label <- paste(sprintf("%.1f", data$label),"%")

ggplot(data, aes( x = Var2, y=percentage,fill = Var1))+
  geom_col(position = 'stack', width = 0.8)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  theme_bw()+
  scale_fill_manual(values=seurat.col)+ 
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

ggsave(filename = "snRNA-seq_injury/results/SCT.CCA/basic processing/stack barplot.pdf",width = 4.5,height = 5)


# all celltype DEGs analysis in injury vs. control (ZINBWAVE method)
Idents(seurat.obj.combined) <- "celltype"
celltype_name <- as.character(unique(seurat.obj.combined$celltype))
for (i in celltype_name){
cluster_choose <- make.names(paste0(i))
celltype.sub <- subset(seurat.obj.combined, idents=c(i)) 
## Convert to SingleCellExperiment object
sce <- as.SingleCellExperiment(celltype.sub, assay="RNA")

## Filter out lowly expressed genes
## Identify genes that have a count of at least 3 in at least 6 cells
filter <- rowSums(assay(sce)>2)>5
## Generate information on how many genes meet these criteria
table(filter)
## Filter matrix based on genes above
sce.filt <- sce[filter,]
## View matrix stats after filtering
sce.filt
rm(sce)
## Select highly variable genes on which to focus analysis
## Model gene variance
sce.var <- modelGeneVar(sce.filt)
## Using gene variance, identify the top 2000
keep <- getTopHVGs(sce.var, n=2000)
## Filter matrix for only these genes
sce.filt <- sce.filt[keep,]
rm(sce.var)
## Generate observational weights for genes
## Create zero-inflated negative binomial regression model to the data (note: this step is VERY computationally intensive)
## Convert SingleCellExperiment object counts to matrix format
## (oberservational weights can only be computed from SummarizedExperiment or matrix)
sce.filt@assays@data@listData$counts <- as.matrix(sce.filt@assays@data@listData$counts)
sce.zinb <- zinbwave(sce.filt, K=0, epsilon=1e12, observationalWeights=TRUE, verbose=TRUE)
## Isolate weights to pass to DGE calculations
weights <- assay(sce.zinb, "weights")
#Save model with weights
saveRDS(sce.zinb, file=paste0('snRNA-seq_injury/results/SCT.CCA/basic processing/zinb_weights_',cluster_choose,'.rds')) 
## DGE analysis with edgeR
## Create design matrix specifying conditions and sex as variables
condition <- factor(sce.zinb$group)  
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
## Make comparison between AD and NS (by specifying AD first, NS becomes control)
colnames(design)<-make.names(colnames(design))
contrast <- makeContrasts(conditioninjury-conditioncontrol, levels=design) 
## Zero-inflation adjusted F test for assessing DGE
## This object will have all of the DGE results
dge_res <- glmWeightedF(fit, contrast=contrast)
## Save DGE so it doesn't need to be recalculated
saveRDS(dge_res, file=paste0('snRNA-seq_injury/results/SCT.CCA/basic processing/zinb_DGE_result_',cluster_choose,'.rds'))  
## Summarize results
## Basic MD plot to visualize results
## Filter DGEs meeting statistical significance mentioned above
genes <- topTags(dge_res, n=2000, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
## Create dataframe with gene name, log2 fold change, and adjusted p-value
genes.sig <- data.frame(row.names(genes))
genes.sig$Pval <- genes$table$FDR
genes.sig$logFC <- genes$table$logFC
genes.sig.lfc <- subset(genes.sig, genes.sig$logFC >= 0.25 | genes.sig$logFC <= -0.25)
## Save dataframe as csv file
openxlsx::write.xlsx(genes.sig.lfc, file=paste0('snRNA-seq_injury/results/SCT.CCA/basic processing/zinb_sig_genes_',cluster_choose,'-0.25lfc.xlsx')) 
}

# DEGs two-way bar charts
df <- read.csv(file = "snRNA-seq_injury/results/SCT.CCA/basic processing/DEG number_log2fc1.csv")
df <- melt(df) 
ggplot(df, aes(
  x = factor(X,levels = unique(X)),             
  y = ifelse(variable == "Up", value, -value),  
  fill = variable)) +
  geom_bar(stat = 'identity')+                                
  #coord_flip()+                                               
  geom_text(                                                 
    aes(label=value,                                          
        #vjust = ifelse(variable == "Up", -0.5, 1),         
        hjust = ifelse(variable == "Up", -0.4, 1.1)),size=0)+
  scale_y_continuous(                                        
    labels = abs,                                           
    expand = expansion(mult = c(0.1, 0.1)))+                
ylab("# of DEGs")+xlab("")+theme_classic()+
scale_y_continuous(breaks = seq(-900, 900, 100))+
  scale_fill_manual(values=c("tomato",
                             "dodgerblue"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
ggsave(filename = "snRNA-seq_injury/results/SCT.CCA/basic processing/Two-way bar charts_log2fc1.pdf",width = 4,height = 5)
