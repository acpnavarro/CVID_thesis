#Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#create a gene list from a df resulted from DESeq2 
df <-read.xlsx("C:/Users/Anna/Documents/CVID_project/DF/df.AB_CVIDxHD.xlsx", rowNames=FALSE)
# log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$row.names
##gene_list_entrez
#Convert gene IDs 
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ids <-as.data.frame(ids)
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$row.names %in% dedup_ids$SYMBOL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZ = dedup_ids$ENTREZID
# Create a vector of the gene unuiverse
entrez_geneList <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(entrez_geneList) <- df2$ENTREZ
# omit any NA values 
entrez_geneList<-na.omit(entrez_geneList)
# sort the list in decreasing order (required for clusterProfiler)
entrez_geneList = sort(entrez_geneList, decreasing = TRUE)
head(entrez_geneList)
###############################################################################

###VENN diagram 1x2###
##read data
venn1 <- read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/INTERS/inters_1x2.xlsx", rowNames=FALSE)

inter_long <- gather(venn1, intersections, symbol, "Group1x2":"Group2", factor_key=T) %>%
  dplyr::filter(symbol != "")

# Annotation
inter_long$Entrez <- mapIds(org.Hs.eg.db,
                            keys=inter_long$symbol,
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")

geneList <- inter_long %>%
  dplyr::select(c(Entrez, intersections)) %>%
  dplyr::filter(!is.na(Entrez))

#Enrichment analyses
set.seed(123)
bp <- compareCluster(geneClusters = Entrez~intersections,
                     data=geneList,
                     fun="enrichGO",
                     OrgDb=organism,
                     universe=names(entrez_geneList),
                     keyType="ENTREZID",
                     ont="BP",
                     pAdjustMethod="BH",
                     minGSSize=10,
                     maxGSSize=500,
                     pvalueCutoff=0.05,
                     readable=TRUE)

p <- as.data.frame(bp) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Group1_Group2")

#kegg
set.seed(123)
kegg_organism = "hsa"
kegg<- compareCluster(geneClusters = Entrez~intersections,
                         data=geneList,
                         organism = kegg_organism,
                         fun="enrichKEGG",
                         universe=names(entrez_geneList),
                         keyType = "kegg",
                         pAdjustMethod="BH",
                         minGSSize=10,
                         maxGSSize=500,
                         pvalueCutoff=0.05)
 
##dotplot
p1 <- as.data.frame(kegg) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG: Group1_Group2")

###############################################################################
###VENN diagram 1x2 UP AND DOWN###
##read data
venn1_ud <- read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/INTERS/1x2_up_down.xlsx", rowNames=FALSE)

inter_long <- gather(venn1_ud, intersections, symbol, "1:2_down":"2_up", factor_key=T) %>%
  dplyr::filter(symbol != "")

# Annotation
inter_long$Entrez <- mapIds(org.Hs.eg.db,
                            keys=inter_long$symbol,
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")

geneList <- inter_long %>%
  dplyr::select(c(Entrez, intersections)) %>%
  dplyr::filter(!is.na(Entrez))

#Enrichment analyses
set.seed(123)
bp <- compareCluster(geneClusters = Entrez~intersections,
                     data=geneList,
                     fun="enrichGO",
                     OrgDb=organism,
                     universe=names(entrez_geneList),
                     keyType="ENTREZID",
                     ont="BP",
                     pAdjustMethod="BH",
                     minGSSize=10,
                     maxGSSize=500,
                     pvalueCutoff=0.05,
                     readable=TRUE)

data_go <- as.data.frame(bp)
write.xlsx(data_go,"C:/Users/Anna/Documents/CVID_project/Venn/GO_1x2_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

p <- as.data.frame(bp) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Group1_Group2")

#kegg
set.seed(123)
kegg_organism = "hsa"
kegg<- compareCluster(geneClusters = Entrez~intersections,
                      data=geneList,
                      organism = kegg_organism,
                      fun="enrichKEGG",
                      universe=names(entrez_geneList),
                      keyType = "kegg",
                      pAdjustMethod="BH",
                      minGSSize=10,
                      maxGSSize=500,
                      pvalueCutoff=0.05)

data_kegg <- as.data.frame(kegg)
write.xlsx(data_kegg,"C:/Users/Anna/Documents/CVID_project/Venn/Kegg_1x2_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

##dotplot
p1 <- as.data.frame(kegg) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG: Group1_Group2")

#DO
set.seed(123)
do<- compareCluster(geneClusters = Entrez~intersections,
                      data=geneList,
                      organism = "hsa",
                      fun="enrichDO",
                      universe=names(entrez_geneList),
                      ont = "DO",
                      pAdjustMethod="BH",
                      minGSSize=10,
                      maxGSSize=500,
                      pvalueCutoff=0.05)
?enrichDO

data_do <- as.data.frame(do)
write.xlsx(data_do,"C:/Users/Anna/Documents/CVID_project/Venn/DO_1x2_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

##dotplot
p1 <- as.data.frame(do) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("DO: Group1_Group2")
#Reactome
set.seed(123)
reactome<- compareCluster(geneClusters = Entrez~intersections,
                      data=geneList,
                      organism = "human",
                      fun="enrichPathway",
                      universe=names(entrez_geneList),
                      pAdjustMethod="BH",
                      minGSSize=10,
                      maxGSSize=500,
                      pvalueCutoff=0.05)
?enrichPathway
data_react<- as.data.frame(reactome)
write.xlsx(data_react,"C:/Users/Anna/Documents/CVID_project/Venn/reactome_1x2_up_down.xlsx", asTable = TRUE, rowNames=FALSE)
#dotplot
p1 <- as.data.frame(reactome) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Reactome: Group1_Group2")

###motif gene sets


?compareCluster
###################################################################################
###VENN diagram 3x4###

##read data
venn2 <- read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/INTERS/inters_3x4.xlsx", rowNames=FALSE)

inter_long <- gather(venn2, intersections, symbol, "Group3x4":"Group4", factor_key=T) %>%
  dplyr::filter(symbol != "")

# Annotation
inter_long$Entrez <- mapIds(org.Hs.eg.db,
                            keys=inter_long$symbol,
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")

geneList2 <- inter_long %>%
  dplyr::select(c(Entrez, intersections)) %>%
  dplyr::filter(!is.na(Entrez))

#Enrichment analyses
set.seed(123)
bp <- compareCluster(geneClusters = Entrez~intersections,
                     data=geneList2,
                     fun="enrichGO",
                     OrgDb=organism,
                     universe=names(entrez_geneList),
                     keyType="ENTREZID",
                     ont="BP",
                     pAdjustMethod="BH",
                     minGSSize=10,
                     maxGSSize=500,
                     pvalueCutoff=0.05,
                     readable=TRUE)
?enrichGO
?gseGO

##dotplot
p2 <- as.data.frame(bp) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Group3_Group4")

#kegg
set.seed(123)
kegg_organism = "hsa"
kegg2 <- compareCluster(geneClusters = Entrez~intersections,
                     data=geneList2,
                     organism = kegg_organism,
                     fun="enrichKEGG",
                     universe=names(entrez_geneList),
                     keyType = "kegg",
                     pAdjustMethod="BH",
                     minGSSize=10,
                     maxGSSize=500,
                     pvalueCutoff=0.05)
##dotplot
p2 <- as.data.frame(kegg2) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG: Group3_Group4")

#emapplot
#emapplot
kegg_1<- pairwise_termsim(kegg) 
kegg_2 <- pairwise_termsim(kegg2)
p2 <- emapplot(kegg_1, cex_category=1.5, layout="kk")
p3 <-emapplot(kegg_2,cex_category=1.5,layout="kk")
cowplot::plot_grid(p2, p3, ncol=2,labels=LETTERS[1:2])
?emapplot

kegg_1<- pairwise_termsim(kegg) 
kegg_2 <- pairwise_termsim(kegg2)
p2 <- emapplot(kegg_1)
p3 <-emapplot(kegg_2)
cowplot::plot_grid(p2, p3, ncol=2)
?emapplot
###VENN diagram 3x4 UP AND DOWN###
##read data
venn2_ud <- read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/INTERS/3x4_up_down.xlsx", rowNames=FALSE)

inter_long <- gather(venn2_ud, intersections, symbol, "3:4_down":"4_up", factor_key=T) %>%
  dplyr::filter(symbol != "")

# Annotation
inter_long$Entrez <- mapIds(org.Hs.eg.db,
                            keys=inter_long$symbol,
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")

geneList <- inter_long %>%
  dplyr::select(c(Entrez, intersections)) %>%
  dplyr::filter(!is.na(Entrez))

#Enrichment analyses
set.seed(123)
bp <- compareCluster(geneClusters = Entrez~intersections,
                     data=geneList,
                     fun="enrichGO",
                     OrgDb=organism,
                     universe=names(entrez_geneList),
                     keyType="ENTREZID",
                     ont="BP",
                     pAdjustMethod="BH",
                     minGSSize=10,
                     maxGSSize=500,
                     pvalueCutoff=0.05,
                     readable=TRUE)

data_go <- as.data.frame(bp)
write.xlsx(data_go,"C:/Users/Anna/Documents/CVID_project/Venn/GO_3x4_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

p <- as.data.frame(bp) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Group3_Group4")

#kegg
set.seed(123)
kegg_organism = "hsa"
kegg<- compareCluster(geneClusters = Entrez~intersections,
                      data=geneList,
                      organism = kegg_organism,
                      fun="enrichKEGG",
                      universe=names(entrez_geneList),
                      keyType = "kegg",
                      pAdjustMethod="BH",
                      minGSSize=10,
                      maxGSSize=500,
                      pvalueCutoff=0.05)

data_kegg <- as.data.frame(kegg)
write.xlsx(data_kegg,"C:/Users/Anna/Documents/CVID_project/Venn/Kegg_3x4_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

##dotplot
p1 <- as.data.frame(kegg) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG: Group3_Group4")

#DO
set.seed(123)
do<- compareCluster(geneClusters = Entrez~intersections,
                    data=geneList,
                    organism = "hsa",
                    fun="enrichDO",
                    universe=names(entrez_geneList),
                    ont = "DO",
                    pAdjustMethod="BH",
                    minGSSize=10,
                    maxGSSize=500,
                    pvalueCutoff=0.05)
?enrichDO

data_do <- as.data.frame(do)
write.xlsx(data_do,"C:/Users/Anna/Documents/CVID_project/Venn/DO_3x4_up_down.xlsx", asTable = TRUE, rowNames=FALSE)

##dotplot
p1 <- as.data.frame(do) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("DO: Group3_Group4")
#Reactome
set.seed(123)
reactome<- compareCluster(geneClusters = Entrez~intersections,
                          data=geneList,
                          organism = "human",
                          fun="enrichPathway",
                          universe=names(entrez_geneList),
                          pAdjustMethod="BH",
                          minGSSize=10,
                          maxGSSize=500,
                          pvalueCutoff=0.05)
?enrichPathway
data_react<- as.data.frame(reactome)
write.xlsx(data_react,"C:/Users/Anna/Documents/CVID_project/Venn/reactome_3x4_up_down.xlsx", asTable = TRUE, rowNames=FALSE)
#dotplot
p1 <- as.data.frame(reactome) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(Cluster) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = Cluster, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Reactome: Group3_Group4")
