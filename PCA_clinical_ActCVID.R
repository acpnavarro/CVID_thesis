#####count_table#####
countData<- read.xlsx("C://Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
countData_CVID_Activated <- countData[, grep("^AB.CVID|^AN.CVID", names(countData))]
####sample_table####
sampleTable=read_excel("C:/Users/Anna/Documents/CVID_project/Excel_tables/clinical_table.xlsx")
sampleTable$Diagnosis =as.factor(sampleTable$Diagnosis)
sampleTable$CellType=as.factor(sampleTable$CellType)
sampleTable$ID=as.factor(sampleTable$ID)
sampleTable$Age=as.factor(sampleTable$Age)
sampleTable$Sex =as.factor(sampleTable$Sex)
sampleTable$OnlyInfections <- as.factor(sampleTable$OnlyInfections)
sampleTable$Group <- as.factor(paste(sampleTable$ID,sampleTable$CellType, sampleTable$Age, sampleTable$Sex,sampleTable$OnlyInfections, sep = "_"))
sampleTable_CVID_Activated = sampleTable %>% dplyr::filter(!grepl("^IB|^IN", Sample))
###DESeq2###
dds <- DESeqDataSetFromMatrix(countData = countData_CVID_Activated,
                              colData = sampleTable_CVID_Activated,
                              design = ~ Group)
##transformation of count data
normal_dds<- vst(dds)

#ggplot
pcaData <- plotPCA(normal_dds, intgroup = c( "ID", "Sex", "OnlyInfections", "CellType"), returnData = TRUE)
attributes(pcaData)

percentVar <- round(100 * attr(pcaData, "percentVar"))
theme_set(theme_classic(base_size = 10) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
pca2<-ggplot(pcaData, aes(x = PC1, y = PC2, color=interaction(OnlyInfections,CellType), shape=Sex)) +
  geom_point(size =3) +
  geom_text(label=dds$ID,  nudge_x = 1, nudge_y = 1, color = "black", cex=2)+ 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  labs(shape="Sex", colour="Only Infections - CellType")+
  scale_shape_manual(values = c(15,17,22,22,22,22))+
  ggtitle("PCA Activated Cells - Clinical Data") 

?geom_text  
  ?ggplot
  
