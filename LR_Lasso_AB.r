count<- read.xlsx("C:/Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
df_AB<-count[, grep("^AB", names(count))]
tdf_AB<-t(df_AB)
tdf_AB<-as.data.frame(tdf_AB)
Group=rep(c("CVID","HD"),c(10,4))
Group=factor(Group)
#LASSO
library(ncvreg)
res.ncvreg = cv.ncvreg(X=scale(tdf_AB), y=Group, family="binomial")
barplot(coef(res.ncvreg)[-1], las=2)
#Logistic regression with CHADL
fit <- glm(Group~NENF+SMIM15+SS18L1,data=tdf_AB,family=binomial)
step(fit, direction="backward")
stripchart(tdf_AB$SMIM15~Group,vertical = T, method="jitter",ylab="Read count",main="TBL2")
#ROC
library(pROC)
roc1=roc(Group, tdf_AB$NENF + tdf_AB$SMIM15+ tdf_AB$SS18L1,plot = TRUE, print.auc = TRUE) ## AUC = 0.83
?roc
plot(roc1)
wilcox.test(tdf_AB$SMIM15 + tdf_AB$NENF~Group) ## p value > 0.076
ci.auc(roc1)
data.frame(roc1$thresholds,roc1$sensitivities,roc1$specificities)
coords(roc1,"best", best.method="youden",transpose = FALSE)
coords(roc1, "best", best.method="closest.topleft",transpose = FALSE)
stripchart(tdf_AB$TBL2~Group,data=df_AB,vertical = T,xlim=c(0.7,2.3),method="jitter",ylab="Read Count",main="TBL2 gene")
abline(h=781.5,lty=2,col="red")
#LOOCV
library(caret)
df3=data.frame(genes= tdf_AB$NENF + tdf_AB$SMIM15,Group=Group)
fit4 = train(Group ~ genes, data=df3, method="glm",family="binomial",
             trControl = trainControl(method = "LOOCV"))
fit4$results




