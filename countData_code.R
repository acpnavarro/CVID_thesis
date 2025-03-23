###Prepare the countData
count<- read.table(file="C:\\Users\\Anna\\Documents\\CVID_project\\Counts_FeatureCounts\\counts_featureCounts.txt", header=TRUE, row.names = "Geneid")
#select the columns
x=1:5
new_count<- subset(count, select=-c(x))
#remove part of the names from the columns names
colnames(new_count)<-gsub ("X.media.anna.4TB_2.Anna_Navarro.STAR.","",colnames (new_count))
colnames(new_count)<-gsub ("_Aligned.sortedByCoord.out.bam","",colnames (new_count))
new_count[1:7, 1:7]
#filter low counts
meanLog2CPM <- rowMeans(log2(cpm(new_count) + 1))
hist(meanLog2CPM, main="Histogram of Count Data")
sum(meanLog2CPM <= 1)
new_count <- new_count[meanLog2CPM > 1, ]
dim(new_count)
new_count[1:20, 1:7]
#create a column with gene symbol 
new_count$symbol=mapIds(org.Hs.eg.db, keys=rownames(new_count), keytype="ENSEMBL", column="SYMBOL")
#check for NAs and duplicated gene symbol
sum(duplicated(new_count$symbol)) ### 3825 (out of 17849) duplicated/NAs genes
## Gene Expression data filtering: 
countData <- new_count %>%
  dplyr::filter(!is.na(symbol)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(RowSum = sum(c_across(where(is.numeric)))) %>%
  dplyr::filter(RowSum > 0) %>%
  dplyr::arrange(desc(RowSum)) %>%
  dplyr::distinct(symbol, .keep_all = T) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(-c(RowSum, Gene)) %>%
  tibble::column_to_rownames("symbol")
