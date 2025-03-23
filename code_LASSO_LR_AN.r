count=read.xlsx("/media/anna/4TB_2/Anna_Navarro/countData.xlsx", rowNames=TRUE)
df_AN<-count[, grep("^AN", names(count))]
tdf_AN<-t(df_AN)
tdf_AN<-as.data.frame(tdf_AN)
Group=rep(c("CVID","HD"),c(11,4))
Group=factor(Group)
#LASSO
library(ncvreg)
res.ncvreg = cv.ncvreg(X=scale(tdf_AN), y=Group, family="binomial")
barplot(coef(res.ncvreg)[-1], las=2)
#Logistic regression with CHADL
fit <- glm(Group~CHADL,data=tdf_AN,family=binomial)
step(fit, direction="backward")
stripchart(tdf_AN$CHADL~Group,vertical = T, method="jitter",ylab="Read count",main="CHADL")
#ROC
library(pROC)
roc1=roc(Group, tdf_AN$CHADL) ==> AUC = 0.84
plot(roc1)
wilcox.test(tdf_AN$CHADL~Group) ===> p value > 0.05
ci.auc(roc1)
data.frame(roc1$thresholds,roc1$sensitivities,roc1$specificities)
coords(roc1,"best", best.method="youden",transpose = FALSE)
coords(roc1, "best", best.method="closest.topleft",transpose = FALSE)
stripchart(tdf_AN$CHADL~Group,data=df_AN,vertical = T,xlim=c(0.7,2.3),method="jitter",ylab="Read Count",main="CHADL gene")
abline(h=674.5,lty=2,col="red")
##########################################################################

tpm_countTable<- read.xlsx("~/tpm_countTable.xlsx", rowNames=TRUE)
df_clinical <- read.xlsx("~/table_AN.xlsx", rowNames=FALSE)