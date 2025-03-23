count<- read.xlsx("C:/Users/Anna/Documents/CVID_project/countData.xlsx", rowNames=TRUE)
df_AN<-count[, grep("^AN", names(count))]
tdf_AN<-t(df_AN)
tdf_AN<-as.data.frame(tdf_AN)
Group=rep(c("CVID","HD"),c(11,4))
Group=factor(Group)
#LASSO
library(ncvreg)
set.seed(123)
res.ncvreg = cv.ncvreg(X=scale(tdf_AN), y=Group, family="binomial")
?cv.ncvreg
barplot(coef(res.ncvreg)[-1], las=2)
plot(res.ncvreg)
summary(res.ncvreg)
#Logistic regression with DNPH1
fit <- glm(Group~ATP7A,data=tdf_AN,family=binomial)
fit <- glm(Group~ZNF561,data=tdf_AN,family=binomial)
fit <- glm(Group~WDR76,data=tdf_AN,family=binomial)
fit <- glm(Group~KLF4,data=tdf_AN,family=binomial)
fit <- glm(Group~WDR76 + ATP7A,data=tdf_AN,family=binomial)
?glm
step(fit, direction="backward")
stripchart(tdf_AN$ATP7A~Group,vertical = T, method="jitter",ylab="Read count",main="ATP7A")
stripchart(tdf_AN$KLF4~Group,vertical = T, method="jitter",ylab="Read count",main="KLF4")
#ROC
install.packages("pROC")
library(pROC)
roc1=roc(Group, tdf_AN$ATP7A)
roc1=roc(Group, tdf_AN$KLF4)
roc1=roc(Group, tdf_AN$WDR76)
wilcox.test(tdf_AN$ATP7A~Group)
plot(roc1)
ci.auc(roc1) 
data.frame(roc1$thresholds,roc1$sensitivities,roc1$specificities)
coords(roc1,"best", best.method="youden",transpose = FALSE)
coords(roc1, "best", best.method="closest.topleft",transpose = FALSE)
#stripchart with the threshold
stripchart(tdf_AN$KLF4~Group,data=df_AN,vertical = T,xlim=c(0.7,2.3),ylim=c(0,100),method="jitter",ylab="Read Count",main="ATP7A gene")
abline(h=37.5,lty=2,col="red")
#LOOCV
install.packages("caret")
library(caret)
?caretFuncs
df3=data.frame(ATP7A=tdf_AN$ATP7A,Group=Group)
fit4 = train(Group ~ ATP7A, data=df3, method="glm",family="binomial",
             trControl = trainControl(method = "LOOCV"))
fit4$results
df3=data.frame(ZNF561=tdf_AN$ZNF561,Group=Group)
fit4 = train(Group ~ ZNF561, data=df3, method="glm",family="binomial",
             trControl = trainControl(method = "LOOCV"))
fit4$results

test_prob = predict(model_glm, newdata = default_tst, type = "response")
test_roc = roc(Group, tdf_AN$ATP7A, plot = TRUE, print.auc = TRUE)
