#Load packages
library(DESeq2)
library(tibble)
library(openxlsx)
library(data.table)
library(stringr) 
library(ggplot2)
library(ggprism)
options(stringsAsFactors = F)

mydata<- fread("bulk_seq/gene_count_matrix.csv")
any(duplicated(mydata$gene_id))
genename_ref<- fread("bulk_seq/gene_id.csv")
str(mydata)
head(mydata)
mydata<-as.data.frame(mydata)
rownames(mydata)<-mydata[,1]
mydata<-mydata[,-1]
head(mydata)

filter_counts<-mydata[rowSums(mydata)>0,] 
del<-c()
for (i in 1:nrow(filter_counts)) {
  whether <- filter_counts[i,]
  if(length(which(whether==0))>25){
    del<-rbind(del,whether)
  }
} 
mydata<-filter_counts[-which(rownames(filter_counts) %in% rownames(del)),]
rm(del,filter_counts,whether,i);gc()

# Building dds objects
group_list=c(rep(c('BL_3dpi'),5),rep(c('BL_ctrl'),5), rep(c('SB_3dpi'),5),rep(c('SB_ctrl'),5), rep(c('SR_3dpi'),5),rep(c('SR_ctrl'),5))
group_list <- factor(group_list)
str(group_list)
colData <- data.frame(row.names=colnames(mydata), group_list=group_list)
rm(group_list)

dds <- DESeqDataSetFromMatrix(mydata,colData = colData,design = ~ group_list)

dds<-DESeq(dds)

# variance stabilizing transformation (VST)
vsd <- vst(dds, blind = TRUE)
vsd@colData
# PCA
pca_data <- plotPCA(vsd, intgroup = "group_list", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# PCA plot
ggplot(pca_data, aes(PC1, PC2, color = group_list)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

ggplot(pca_data, aes(PC1, PC2, color = group_list, label = rownames(pca_data))) +
  geom_point(size = 4, alpha = 0.8) +  # 适当透明度
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic(base_size = 14) +  # 经典主题 + 适当字体大小
  theme(
    legend.position = "right",  # 图例位置
    legend.title = element_blank(),  # 去掉图例标题
    axis.text = element_text(size = 12),  # 坐标轴字体
    axis.title = element_text(size = 14, face = "bold"),  # 加粗坐标轴标题
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),  # 主要网格线（浅色虚线）
    panel.grid.minor = element_line(color = "grey90", linetype = "dotted")   # 次要网格线
  )+stat_ellipse(level = 0.8)
ggsave(filename = "bulk_seq/plots/pca.pdf",width = 6.2,height = 5)

# FPKM distribution
fpkm <- read.csv("gene_count_matrix_FPKM.csv", row.names = 1)
fpkm <- fpkm[rownames(fpkm) %in% rownames(mydata),]
fpkm_long <- fpkm %>%
  rownames_to_column(var = "GeneID") %>%
  pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "FPKM") %>%
  mutate(
    log2FPKM = log2(FPKM + 1),
    Group = str_extract(Sample, "^(BL_3dpi|BL_ctrl|SB_3dpi|SB_ctrl|SR_3dpi|SR_ctrl)")
  )

group_colors <- c(
  "BL_ctrl" = "#F8766D",
  "BL_3dpi" = "#B79F00",
  "SB_ctrl" = "#00BA38",
  "SB_3dpi" = "#00BFC4",
  "SR_ctrl" = "#619CFF",
  "SR_3dpi" = "#F564E3"
)

ggplot(fpkm_long, aes(x = Sample, y = log2FPKM, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.3)+ # add box whiskers
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", size = 0.2) +  # 去掉离群值
  scale_fill_manual(values = group_colors) +
  ylim(0, 7.5) +
  labs(
    title = "FPKM Expression Distribution Across Samples",
    y = "log2(FPKM + 1)",
    x = NULL,
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(filename = "bulk_seq/plots/boxplot.pdf",height = 5,width = 12)


# Cibersort data input and basic processing
df <- read.csv("bulk_seq/bl_sr_sb_FPKM_data.txt",header = T,sep = "\t")
df <- df[!duplicated(df$Gene), ]
colnames(df) <- make.names(colnames(df), unique = TRUE)
df <- df[complete.cases(df), ]
df <- df[, colnames(df) != ""]
df <- df[, !apply(df, 2, function(x) all(is.na(x)))]
head(df)  # 查看预处理后的数据框
write.table(df,file = "bulk_seq/bl_sr_sb_FPKM_data_processed.txt",sep = "\t",col.names = T,row.names = F)

library(CIBERSORT)
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
mixture_file <- "bulk_seq/bl_sr_sb_FPKM_data_processed_upper.txt"
results <- cibersort(sig_matrix, mixture_file)
saveRDS(results,file = "bulk_seq/results_upper.rds")

# group information
rownames(results)
phenotype = as.data.frame(c(rep("BL_3dpi",5),rep("BL_ctrl",5),rep("SB_3dpi",5),rep("SB_ctrl",5),rep("SR_3dpi",5),rep("SR_ctrl",5)))
colnames(phenotype) <- "group"
phenotype$sample <- rownames(results)

group_list <- phenotype$group %>%
  factor(.)
table(group_list) 

# plot
data <- as.data.frame(results[,1:22])
data$group <- group_list
data$sample <- row.names(data)
write.xlsx(data,file = "bulk_seq/SB_SR_BL_immune_prediction_data.xlsx",rowNames=T)

data_new = melt(data)

colnames(data_new)=c("Group","Sample","Celltype","Composition")  #设置行名
head(data_new)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

problem_celltypes <- data_new %>%
  group_by(Celltype) %>%
  summarise(p_value = tryCatch({
    summary(aov(Composition ~ Group))[[1]][["Pr(>F)"]][1]
  }, error = function(e) NA)) %>%
  filter(is.na(p_value))

print(problem_celltypes)


data_filtered <- data_new %>%
  filter(!(Celltype %in% problem_cells$Celltype))

box_data <- ggplot(data_filtered, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "Immnue Cell Prediction")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#EB7369", "#1CB4B8", "#660099", "#006600", "#CC0066", "#FF9900"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "anova",
                     hide.ns = T)+
  theme(plot.title = element_text(size = 16))  
box_data
ggsave("bulk_seq/SB_SR_BL_immune_prediction_simplify.pdf",height=6,width=12)



