#Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#create a gene list from a df resulted from DESeq2 
df <-read.xlsx("C:/Users/Anna/Documents/CVID_project/DF/df.IB_CVIDxHD.xlsx", rowNames=FALSE)
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
##gseGO
set.seed(123)
gse <- gseGO(geneList=entrez_geneList, 
             ont ="BP", 
             keyType = "ENTREZID",
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH",
             eps =0,
             seed=TRUE,
             exponent =1,
             by="fgsea")

data_go <- as.data.frame(gse)
write.xlsx(data_go,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/GO.xlsx", asTable = TRUE, rowNames=TRUE)

#dotplot
p <- as.data.frame(gse) %>%
  mutate(
    num_genes = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = num_genes / setSize) %>%
  dplyr::mutate(type = dplyr::if_else(NES > 0, "upregulated", "downregulated")) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:15) %>%
  ggplot(., aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  facet_grid(~ type) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Nonactivated Bright CVID x HD")

#bar
p <- as.data.frame(gse) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Nonactivated Memory CVID x HD")

#rideplot
ridgeplot(gse) + 
  labs(x= "Enrichment Distribution" ) +
  theme(axis.title.x = element_text(size = 10,
                                    color = "black",
                                    face = "bold")) +
  geom_density_ridges() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") + 
  theme_ridges(grid = FALSE)
#net
pwt<-pairwise_termsim(gse)
emapplot(pwt,showCategory=10) +
  ggtitle("GO: Nonactivated Bright CVID x HD")+
  theme(plot.title= element_text(color="black", size=20, face="bold"))

##############################################################################################
##KEGG
#kegg
set.seed(123)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = entrez_geneList,
               organism     = kegg_organism,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "kegg",
               verbose=TRUE,
               seed=TRUE,
               by="fgsea")

data_kegg <- as.data.frame(kk2)
write.xlsx(data_kegg,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/KEGG.xlsx", asTable = TRUE, rowNames=TRUE)


#dotplot
p <- as.data.frame(kk2) %>%
  mutate(
    num_genes = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = num_genes / setSize) %>%
  dplyr::mutate(type = dplyr::if_else(NES > 0, "upregulated", "downregulated")) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:30) %>%
  ggplot(., aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  facet_grid(~ type) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG: Nonactivated Bright CVID x HD")

#emapplot
pwt<-pairwise_termsim(kk2) 
p1=emapplot(pwt,showCategory=10, node_label="category") +
  ggtitle("KEGG: Nonactivated Bright CVID x HD")+
  theme(plot.title= element_text(color="black", size=10, face="bold"))
?emapplot
# categorySize can be either 'pvalue' or 'geneNum'
p2=cnetplot(kk2, categorySize="pvalue", foldChange=entrez_geneList)
cowplot::plot_grid(p1, p2, ncol=2)

#running score
plot<- gseaplot(kk2, geneSetID = 1, by = "runningScore", title = "KEGG: Nonactivated Bright CVID x HD")
#pathview
library(pathview)
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=entrez_geneList, pathway.id="hsa04660", species = kegg_organism)

d=as.data.frame(kk2)

##REACTOME
set.seed(123)
react <- gsePathway(geneList = entrez_geneList,
                    organism = "human",
                    minGSSize = 10,
                    maxGSSize  = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    verbose=TRUE,
                    seed=TRUE,
                    by="fgsea")
data_react <- as.data.frame(react)
write.xlsx(data_react,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/Reactome.xlsx", asTable = TRUE, rowNames=FALSE)


?gsePathway
p <- as.data.frame(react) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:15) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Reactome: Nonactivated Bright CVID x HD")
#upsetplot
upsetplot(react)
##DOSE
set.seed(123)
do <- gseDO(geneList=entrez_geneList,
            organism = "hsa",
            minGSSize = 10,
            maxGSSize  = 500,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            verbose=TRUE,
            seed=TRUE,
            by="fgsea")
data_do <- as.data.frame(do)
write.xlsx(data_do,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/DO.xlsx", asTable = TRUE, rowNames=FALSE)
p <- as.data.frame(do) %>%
  mutate(
    num_genes = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = num_genes / setSize) %>%
  dplyr::mutate(type = dplyr::if_else(NES > 0, "upregulated", "downregulated")) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  facet_grid(~ type) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("DO: Nonactivated Bright CVID x HD")
#barplot
p <- as.data.frame(do) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("DO: Nonactivated Bright CVID x HD")
##############################################################
###motif gene sets
library(msigdbr)
msigdbr_show_species()
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
#GSEA motifs
C3_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)
em2 <- GSEA(entrez_geneList, TERM2GENE = C3_t2g)
head(em2)

data_em2 <- as.data.frame(em2)
write.xlsx(data_em2,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/Motifs.xlsx", asTable = TRUE, rowNames=FALSE)

p <- as.data.frame(em2) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:15) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Motifs: Nonactivated Bright CVID x HD")

##C7 - immune signatures
C7_t2g <- msigdbr(species = "Homo sapiens", category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C7_t2g)
em3<- GSEA(entrez_geneList, TERM2GENE = C7_t2g)
head(em3)

data_em3 <- as.data.frame(em3)
write.xlsx(data_em3,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/Signatures.xlsx", asTable = TRUE, rowNames=FALSE)

p <- as.data.frame(em3) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Signatures:  Nonactivated Bright CVID x HD")

##hallmarks
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(hall)
em4<- GSEA(entrez_geneList, TERM2GENE = hall)
head(em4)

data_em4 <- as.data.frame(em4)
write.xlsx(data_em4,"C:/Users/Anna/Documents/CVID_project/GSEA/NonactBright_CVIDxHD/Hallmarks.xlsx", asTable = TRUE, rowNames=FALSE)

p <- as.data.frame(em4) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:30) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Hallmark pathways: Nonactivated Bright CVID x HD")


