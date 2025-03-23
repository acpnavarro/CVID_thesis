#####count_table#####
countData<- read.xlsx("C://Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
countData_Naive <- countData[, grep("^AN|^IN", names(countData))]
countData_Bright <- countData[, grep("^AB|^IB", names(countData))]
####sample_table####
sampleTable=read_excel("C:/Users/Anna/Documents/CVID_project/Excel_tables/table2_dds.xlsx")
sampleTable$Condition=as.factor(sampleTable$Condition)
sampleTable$Sample=as.factor(sampleTable$Sample)
sampleTable$CellType=as.factor(sampleTable$CellType)
sampleTable$Status=as.factor(sampleTable$Status)
sampleTable$Batch=as.factor(sampleTable$Batch)
sampleTable$Group <- as.factor(paste(sampleTable$Condition, sampleTable$Status, sampleTable$CellType, sep = "_"))
sampleTable_Naive = sampleTable %>% dplyr::filter(!grepl("^AB|^IB", Sample))
sampleTable_Bright = sampleTable %>% dplyr::filter(!grepl("^AN|^IN", Sample))
###DESeq2###
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleTable,
                              design = ~ Group)
##transformation of count data
normal_dds<- vst(dds)

#ggplot
pcaData <- plotPCA(normal_dds, intgroup = c("Sample", "Condition", "Status", "CellType"), returnData = TRUE)
attributes(pcaData)

percentVar <- round(100 * attr(pcaData, "percentVar"))
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
pca2<-ggplot(pcaData, aes(x = PC1, y = PC2, color=interaction(Condition,CellType), shape=Status)) +
  geom_point(size =3) +
  geom_text(label=dds$Sample,  nudge_x = 1, nudge_y = 1, color = "black", cex=2)+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(shape="Status", colour="Condition/Cell Type")+
  scale_shape_manual(values = c(15,17,22,22,22,22))+
  ggtitle("PCA VST data") 
####
dds <- DESeqDataSetFromMatrix(countData = countData_Bright,
                              colData = sampleTable_Bright,
                              design = ~ Group)
##transformation of count data
normal_dds<- vst(dds)

#ggplot
pcaData <- plotPCA(normal_dds, intgroup = c("Condition", "Status"), returnData = TRUE)
attributes(pcaData)

percentVar <- round(100 * attr(pcaData, "percentVar"))
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
pca2<-ggplot(pcaData, aes(x = PC1, y = PC2, color=Condition, shape=Status)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(shape="Status", colour="Condition")+
  scale_shape_manual(values = c(15,17,22,22,22,22))+
  ggtitle("PCA Bright Memory B cells") 

boxplot(countData$AB.CVID01)

