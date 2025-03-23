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
summary(res.ncvreg)
#Logistic regression with CHADL
fit <- glm(Group~ARHGAP32,data=tdf_AB,family=binomial)
step(fit, direction="backward")
#ROC
library(pROC)
roc1=roc(Group, tdf_AB$ARHGAP32) 
plot(roc1)
plot(roc1, print.auc=TRUE, main="ROC ARHGAP32 gene")
t.test(tdf_AB$ARHGAP32~Group)
ci.auc(roc1)
data.frame(roc1$thresholds,roc1$sensitivities,roc1$specificities)
coords(roc1,"best", best.method="youden",transpose = FALSE)
coords(roc1, "best", best.method="closest.topleft",transpose = FALSE)
stripchart(tdf_AB$ARHGAP32~Group,data=df_AB,vertical = T,xlim=c(0.7,2.3),ylim=c(0,100),method="jitter",ylab="Read Count",col="darkblue", pch=0,cex=1.5, main="ARHGAP32 gene")
abline(h=42,lty=2,col="red")
#LOOC
library(caret)
df3=data.frame(ARHGAP32=tdf_AB$ARHGAP32,Group=Group)
fit4 = train(Group ~ARHGAP32, data=df3, method="glm",family="binomial",
             trControl = trainControl(method = "LOOCV"))
fit4$results
