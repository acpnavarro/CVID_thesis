install.packages('tidyverse')
BiocManager::install("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
install.packages('pheatmap')
install.packages("DOSE")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot", force=TRUE)
install.packages("ggupset")
install.packages("genekitr")
install.packages("ggridges")
BiocManager::install(c("qvalue","GO.db"))
BiocManager::install("pathview")
install.packages("ggplot2")
install.packages("openxlsx")
BiocManager::install("BiocParallel", force=TRUE)
BiocManager::install("GSEA.barplot")
install.packages("wordcloud")
BiocManager::install("fgsea")
###############################################################################
library(wordcloud)
library(dplyr)
library(BiocManager)
library(BiocParallel)
library(ggridges)
library(DOSE)
library(ggridges)
library(genekitr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(openxlsx)
library(dplyr)
library(stringr)
library(fgsea)
library(pathview)
library(ReactomePA)
#Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
#read df
df <-read.xlsx("C:/Users/Anna/Documents/CVID_project/DF/df.IN_CVIDxHD.xlsx", rowNames=TRUE)
# log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- row.names(df)
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)
##gseGO
set.seed(123)
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL",
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


#dotplot
p <- as.data.frame(gse) %>%
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
                  ggtitle("GO pathway enrichment - Nonactivated Naïve CVID x HD")

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
  ggtitle("Nonactivated Naïve CVID x HD")+
  theme(plot.title= element_text(color="black", size=20, face="bold"))

##############################################################################################
##KEGG
#Convert gene IDs for gseKEGG function
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ids <-as.data.frame(ids)
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$SYMBOL %in% dedup_ids$SYMBOL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZ = dedup_ids$ENTREZID
# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZ
# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
#kegg
set.seed(123)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "kegg",
               verbose=TRUE,
               seed=TRUE,
               by="fgsea")
head(kegg_gene_list)
#dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

p <- as.data.frame(kk2) %>%
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
  ggtitle("KEGG - Nonactivated Naïve CVID x HD")


#emapplot
pwt<-pairwise_termsim(kk2) 
emapplot(pwt,showCategory=10, node_label="group") +
  ggtitle("Nonactivated Naïve CVID x HD")+
  theme(plot.title= element_text(color="black", size=10, face="bold.italic"))
?emapplot
# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
#ridgeplot
ridgeplot(kk2) + 
  labs(x= "Enrichment Distribution" ) +
  theme(axis.title.x = element_text(size = 10,
                                    color = "black",
                                    face = "bold")) +
  geom_density_ridges() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") + 
  theme_ridges(grid = FALSE)
#pathview
library(pathview)
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05150", species = kegg_organism)

d=as.data.frame(kk2)

##REACTOME
react <- gsePathway(geneList=kegg_gene_list, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE)
?gsePathway
react2=data.frame(react)
viewPathway("Biosynthesis of specialized proresolving mediators (SPMs)", 
            readable = TRUE, 
            foldChange = geneList)
?viewPathway
##DOSE
do <- gseDO(geneList=kegg_gene_list,
           minGSSize     = 120,
           pvalueCutoff  = 0.05,
           pAdjustMethod = "BH",
           verbose       = FALSE)

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
  ggtitle("DO - Nonactivated Naïve CVID x HD")
#########################
library(enrichplot)

p <- as.data.frame(do) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("KEGG - Nonactivated Naïve CVID x HD")

head(gene_list)
names(kegg_gene_list)
