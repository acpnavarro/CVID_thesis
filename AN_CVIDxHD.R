library(dplyr)
library("edgeR")
library("DESeq2")
library("readxl")
library(pheatmap)
library(sm)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(ggrepel)
library(tidyverse) 
library("org.Hs.eg.db")
library("writexl")
###############################################################################
#####count_table#####
countData<- read.xlsx("C:/Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
####sample_table####
sampleTable=read_excel("C:/Users/Anna/Documents/CVID_project/Excel_tables/table2_dds.xlsx")
sampleTable$Condition=as.factor(sampleTable$Condition)
sampleTable$CellType=as.factor(sampleTable$CellType)
sampleTable$Status=as.factor(sampleTable$Status)
sampleTable$Batch=as.factor(sampleTable$Batch)
sampleTable$Group <- as.factor(paste(sampleTable$Condition, sampleTable$Status, sampleTable$CellType, sep = "_"))

###DESeq2###
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleTable,
                              design = ~ Group)
##transformation of count data
normal_dds<- vst(dds)
##save
saveRDS(normal_dds, "C://Users//Anna//Documents//CVID_project//DDS//norm_dds.RDS")
readRDS("C://Users//Anna//Documents//CVID_project//DDS/norm_dds.RDS")
##dds2
dds2 <-DESeq(dds)
resultsNames(dds2)
res <- results(dds2, contrast = c("Group", "CVID_Activated_NaïveBCs", "HD_Activated_NaïveBCs"))

###pheatmap package###
DEGs=data.frame(res[(res$padj)<0.05,])
#save the data frame in excel
write.xlsx(DEGs,"C:/Users/Anna/Documents/CVID_project/DEGs/AN_CVIDxHD/pvalue0.05.xlsx", asTable = TRUE, rowNames=TRUE)
DEG2= DEGs%>% 
  filter(abs(log2FoldChange) > 1.5)%>%
  arrange(desc(abs(log2FoldChange)))
#Filter
top25_down <- DEGs %>%
  filter(log2FoldChange < -1.5) %>%
  arrange(desc(abs(log2FoldChange)))  %>%
  dplyr::slice(1:25)

top25_up<- DEGs %>%
  filter(log2FoldChange > 1.5) %>%
  arrange(desc(abs(log2FoldChange)))  %>%
  dplyr::slice(1:25)

#############
#Filter
top5_down <- DEGs %>%
  filter(log2FoldChange < -1.5) %>%
  arrange(desc(abs(log2FoldChange)))  %>%
  dplyr::slice(1:5)

top5_up<- DEGs %>%
  filter(log2FoldChange > 1.5) %>%
  arrange(desc(abs(log2FoldChange)))  %>%
  dplyr::slice(1:5)
top10 <- rbind (top5_up, top5_down)
write.xlsx(top10,"C:/Users/Anna/Documents/CVID_project/DEGs/AN_CVIDxHD/top10.xlsx", asTable = TRUE, rowNames=TRUE)

########
#data frame with 50 most DEGs according to log2FC
top50 <- rbind (top25_up, top25_down)
write.xlsx(top50,"C:/Users/Anna/Documents/CVID_project/DEGs/AN_CVIDxHD/top50.xlsx", asTable = TRUE, rowNames=TRUE)
##prepare the count table
mat<-assay(normal_dds)[rownames(top50),]
#Z score
s_mat=t(scale(t(mat)))
# select only the activated naive cell samples patients from the count table
s_mat=as.data.frame(s_mat)
s_mat2<-s_mat[, grep("^AN", names(s_mat))]

##prepare annotation sample table
ann_sample<-data.frame(sampleTable$Sample, sampleTable$Group)
colnames(ann_sample)<-c("Sample","Group")
df_ann=ann_sample %>% 
  filter(str_detect(Sample, "^AN"))
df_ann <- data.frame(df_ann, row.names = 1)
df_ann$Group<-gsub("_Activated_NaïveBCs", "", df_ann$Group)

#prepare the symptoms groups for CVID
df_clinical<- read.xlsx("C:/Users/Anna/Documents/CVID_project/Excel_tables/df_clinical2.xlsx", rowNames=TRUE)
df_AN=df_clinical %>% 
  filter(str_detect(row.names(df_clinical), "^AN"))
###df_ann=rownames_to_column(df_AN)
###colnames(df_ann)[colnames(df_ann) == "rowname"] <- "Sample" 
df_ann2=cbind(df_ann,df_AN)
###df_ann2=df_ann2[,-1]

##Colors for the groups
color_ann = list(df_ann2, Group=c("CVID"  = "#048907",
                                  "HD" = "#890486"),
                 Only_Infections=c("no"=	"#4040a3", "yes"="#894a04", " "="grey"))

#heatmap
s_mat2<-as.matrix(s_mat2)
heat_plot <- pheatmap(s_mat2,scale="none", clustering_method="average", name="Z score",
                      cluster_rows = T, cluster_cols = T, 
                      show_rownames=TRUE,
                      row_labels=rownames(top50),
                      clustering_distance_rows = 'euclidean',
                      annotation_col = df_ann2,
                      annotation_colors = color_ann, 
                      fontsize_row = 6,
                      legend_breaks = c(-2, 0, 2),
                      legend_labels = c("Low", "Medium", "High"),
                      angle_col = c("45"), 
                      show_colnames = T,
                      fontsize_col = 6,
                      main = "Activated Naïve B cells") 

#######################################################################
###volcano plot###
df=data.frame(res)
write.xlsx(df,"C:/Users/Anna/Documents/CVID_project/DF/df.AN_CVIDxHD.xlsx", asTable = TRUE, rowNames=TRUE)
# Theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.5), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.5), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"
# if log2Foldchange > 1.5 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$log2FoldChange > 1.5 & df$pvalue < 0.05] <- "UP"
# if log2Foldchange < -1.5 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -1.5 & df$pvalue < 0.05] <- "DOWN"
# Explore a bit
head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])
# create a column for symbol in df 
df$symbol <- rownames(df)
# create a column delabel in df and select the top 20 symbols
df$delabel <- ifelse(df$symbol %in% head(df[order(df$padj), "symbol"], 10), df$symbol, NA)
head(df$delabel)
ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(aes(color=df$diffexpressed), size = 2, shape=20 ) +
  scale_color_manual(values = c("#00AFBB", "grey", "#e2062c"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  coord_cartesian(ylim = c(0, 15), xlim = c(-5, 10)) +
  labs(color = 'Gene expression',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  theme(axis.title = element_text(size = 30,
                                  color = "black",
                                  face = "bold"))+
  theme(legend.title = element_text(color = "black", size = 15), legend.text = element_text(color = "black", size=10))+
  scale_x_continuous(breaks = seq(-10, 10, 5))+
  labs(title="Activated Naïve B cells: CVID x HD")+
  theme(plot.title = element_text(size = 15))+
  geom_text_repel(aes(label = delabel), size = 2.5, max.overlaps = Inf)



ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(aes(color=df$diffexpressed), size = 2, shape=20 ) +
  scale_color_manual(values = c("#00AFBB", "grey", "#e2062c"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  coord_cartesian(ylim = c(0, 15), xlim = c(-5, 10)) +
  labs(color = 'Gene expression',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  theme(axis.title = element_text(size = 30,
                                  color = "black",
                                  face = "bold"))+
  theme(legend.title = element_text(color = "black", size = 15), legend.text = element_text(color = "black", size=10))+
  scale_x_continuous(breaks = seq(-10, 10, 5))+
  geom_text_repel(aes(label = delabel), size = 2.5, max.overlaps = Inf)
  theme(plot.title = element_text(size = 15))


