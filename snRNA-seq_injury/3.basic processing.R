# load packages
library(Seurat)
library(ggplot2)
library(reshape2)

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

# DEGs two-way bar charts
df <- read.csv(file = "DEG number_log2fc1.csv")
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
