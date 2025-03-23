count<- read.xlsx("C:/Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
df_AN<-count[, grep("^AN", names(count))]
tdf_AN<-t(df_AN)
tdf_AN<-as.data.frame(tdf_AN)
Group=rep(c("CVID","HD"),c(11,4))
Group=factor(Group)
#LASSO
library(ncvreg)
res.ncvreg = cv.ncvreg(X=scale(tdf_AN), y=Group, family="binomial")
barplot(coef(res.ncvreg)[-1], las=2)
?cv.ncvreg
plot(res.ncvreg)
summary(res.ncvreg)
#Logistic regression with DNPH1
fit <- glm(Group~KLF4,data=tdf_AN,family=binomial)
?glm
step(fit, direction="backward")
stripchart(tdf_AN$DNPH1~Group,vertical = T, method="jitter",ylab="Read count",main="RAD51D")
#ROC
library(pROC)
roc1=roc(Group, tdf_AN$KLF4) 
plot(roc1, print.auc=TRUE, main="ROC KLF4 gene")
t.test(tdf_AN$KLF4 ~Group)
ci.auc(roc1)  ## 95% CI: 0.7159-1 (DeLong)
data.frame(roc1$thresholds,roc1$sensitivities,roc1$specificities)
coords(roc1,"best", best.method="youden",transpose = FALSE)
coords(roc1, "best", best.method="closest.topleft",transpose = FALSE)
#stripchart with the threshold
stripchart(tdf_AN$KLF4~Group,data=df_AN,vertical = T,xlim=c(0.7,2.3),ylim=c(0,100),method="jitter",ylab="Read Count",col="darkblue", pch=0,cex=1.5, main="KLF4 gene")
abline(h=37.5,lty=2,col="red")
?stripchart
#LOOCV
library(caret)
df3=data.frame(KLF4=tdf_AN$KLF4,Group=Group)
fit4 = train(Group ~ KLF4, data=df3, method="glm",family="binomial",
             trControl = trainControl(method = "LOOCV"))
fit4$results
##########################################################################
#LASSO with glmnet function and then Logistic Regression
count=read.xlsx("/media/anna/4TB_2/Anna_Navarro/countData.xlsx", rowNames=TRUE)
df_AN<-count[, grep("^AN", names(count))]
#data frame with the samples and the condition (0 or 1)
df_clinical <- read.xlsx("C:/Users/Anna/Documents/CVID_project/Excel_tables/table_AN.xlsx", rowNames=TRUE)
#transpose the count table and join the column with the condition
tdf_AN=as.data.frame(t(df_AN))%>%
  tibble::rownames_to_column("Sample")%>%
  dplyr::inner_join(.,df_clinical)%>%
  tibble::column_to_rownames("Sample")%>%
  dplyr::select(Condition,everything())
## Scale and Center the data
tdf_AN.scale <- base::scale(base::data.matrix(tdf_AN[,2:length(tdf_AN)]), center = T, scale = T)
## Condition as factor
y_tdf_AN <- factor(tdf_AN$Condition)
# LASSO
## perform k-fold cross-validation to find optimal lambda value
set.seed(1234)
install.packages("glmnet")
library(glmnet)
?glmnet
cv_model <- glmnet::cv.glmnet(tdf_AN.scale, y_tdf_AN, alpha = 1, family = "binomial", measure = "class")
## find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min

## get coeficients
z <- stats::coef(cv_model, s= best_lambda)
Genes <- c(0, z@Dimnames[[1]][z@i])
z_df <- base::data.frame(coef = z@i, beta =  z@x, Genes = Genes)
z_df$Genes

genes<-as.data.frame(z_df$Genes)
genes<-as.data.frame(genes[-1,1])

#Logistic Regression
set.seed(1234)
fit <- glm(Condition~PIP4K2C + KLHL42 + SFR1,data=tdf_AN,family=binomial)
summary(fit)



count=read.xlsx("/media/anna/4TB_2/Anna_Navarro/countData.xlsx", rowNames=TRUE)
df_AN<-count[, grep("^AN", names(count))]
tdf_AN<-t(df_AN)
tdf_AN<-as.data.frame(tdf_AN)
Group=rep(c("CVID","HD"),c(11,4))
Group=factor(Group)
cv_model <- glmnet(x=scale(tdf_AN), y=Group, alpha = 1, family = "binomial")
?glmnet
best_lambda <- cv_model$lambda.min

library(glmnet)
cvfit <- cv.glmnet(x=scale(tdf_AN), y=Group, nfolds=15, family = "binomial", type.measure = "class")
?cv.glmnet
plot(cvfit)
cvfit$lambda.min
coef(cvfit)
