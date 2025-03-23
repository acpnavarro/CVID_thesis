library(c)
#####count_table#####
tpm_countTable= read.xlsx("C:/Users/Anna/Documents/CVID_project/tpm_countTable.xlsx", rowNames=TRUE)


df_AB <- tpm_countTable[, grep("^AB.CVID", names(countData))]
write.xlsx(df_AB,"C:/Users/Anna/Documents/CVID_project/df_AB.xlsx", rowNames=TRUE)

df_AB<-read.xlsx("C:/Users/Anna/Documents/CVID_project/df_AB.xlsx", rowNames=TRUE, colNames=TRUE)


tdf_AB=t(df_AB)
group<- row.names(tdf_AB)
data <- cbind(group,tdf_AB)
data <- as.data.frame(data) %>%
  dplyr::mutate(group = if_else(group =="Yes",1,0))


plot(data[,-1],state,xlab="genes",ylab="Probability",xlim=c(0,5))
legend("bottomright",levels(names),col=1:2,pch=1)


#make this example reproducible
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
train <- data[sample, ]
test <- data[!sample, ] 


fit <- glm(group~.,data=train,family=binomial)

#############################################################################
#####count_table#####
tpm_countTable= read.xlsx("/media/anna/4TB_2/Anna_Navarro/tpm_countTable.xlsx", rowNames=TRUE)


df_sample<-read.xlsx("/media/anna/4TB_2/Anna_Navarro/AB_group.xlsx", colNames=TRUE)


tdf_AB=as.data.frame(t(tpm_countTable))%>%
  tibble::rownames_to_column("Sample")%>%
  dplyr::inner_join(.,df_sample)%>%
  tibble::column_to_rownames("Sample")%>%
  dplyr::select(Group,everything())


#make this example reproducible
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(tdf_AB), replace=TRUE, prob=c(0.7,0.3))
train <- tdf_AB[sample, ]
test <- tdf_AB[!sample, ] 


fit <- glm(Group~.,data=train,family=binomial)

summary(fit)

predict(fit, test, type="response")


