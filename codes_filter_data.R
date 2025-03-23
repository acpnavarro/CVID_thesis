BiocManager::install("BiomaRT")
library("readxl")
library("writexl")
library(biomaRt)
##########################################################################################
#countdata
data<- read.xlsx("C://Users/Anna/Documents/CVID_project/countData0.xlsx", rowNames=TRUE)
#annotation with Biomart
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")

results <-  getBM(attributes = c("external_gene_name", "gene_biotype"), 
                  filters = "external_gene_name", 
                  values = list(row.names(as.data.frame(data))),
                  mart = human)%>%
  dplyr::distinct(external_gene_name,.keep_all = TRUE)

results$gene_biotype=as.factor(results$gene_biotype)
levels(results$gene_biotype) 
## join the 2 data frames (results + counts)
df <- results %>% dplyr::inner_join(mutate(data, GeneSymbol = rownames(data)), 
                                    by = c("external_gene_name" ="GeneSymbol"))
#filter the pseudogenes
length(grep("pseudogene$", df$gene_biotype))
count = df %>% dplyr::filter(!grepl("pseudogene$", gene_biotype))
#remove rows with the patterns IGK", "IGL", "IGH", "HLA"
remove.list <- paste(c("IGK", "IGL", "IGH", "HLA"), collapse = '|') 
countData = count %>% filter(!grepl(remove.list, external_gene_name)) %>%
  dplyr::select(-c(gene_biotype))%>%
  dplyr::rowwise() %>%
  dplyr::mutate(RowSum = sum(c_across(where(is.numeric)))) %>%
  dplyr::filter(RowSum > 0) %>%
  dplyr::arrange(desc(RowSum)) %>%
  tibble::column_to_rownames("external_gene_name")
countData = countData[,1:65] 
#save the countData in excel
write.xlsx(countData,"C:/Users/Anna/Documents/CVID_project/countData.xlsx", asTable = TRUE, rowNames=TRUE)

length(grep("HLA", countData$external_gene_name))


         