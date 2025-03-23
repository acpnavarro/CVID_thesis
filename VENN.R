install.packages("VennDiagram")
install.packages("gplots")
install.packages("ggVennDiagram")
library(gplots)
library(VennDiagram)
library(openxlsx)
library(ggvenn)
library(ggVennDiagram)

# upload the files
Group1 = read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/HD_ANxIN.xlsx",rowNames=FALSE)
Group2 = read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/CVID_ANxIN.xlsx",rowNames=FALSE)
Group3 = read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/HD_ABxIB.xlsx", rowNames=FALSE)
Group4 = read.xlsx("C:/Users/Anna/Documents/CVID_project/Venn/CVID_ABxIB.xlsx",rowNames=FALSE)

##Group1
##Divide among up and down regulated genes
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
Group1$diffexpressed <- "NO"
# if log2Foldchange > 2.5 and pvalue < 0.05, set as "UP"
Group1$diffexpressed[Group1$log2FoldChange > 2.5 & Group1$pvalue < 0.05] <- "UP"
# if log2Foldchange < -2.5 and pvalue < 0.05, set as "DOWN"
Group1$diffexpressed[Group1$log2FoldChange < -2.5 & Group1$pvalue < 0.05] <- "DOWN"
# Explore a bit
head(Group1[order(Group1$padj) & Group1$diffexpressed == 'DOWN', ])
##Group2
##Divide among up and down regulated genes
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
Group2$diffexpressed <- "NO"
# if log2Foldchange > 2.5 and pvalue < 0.05, set as "UP"
Group2$diffexpressed[Group2$log2FoldChange > 2.5 & Group2$pvalue < 0.05] <- "UP"
# if log2Foldchange < -2.5 and pvalue < 0.05, set as "DOWN"
Group2$diffexpressed[Group2$log2FoldChange < -2.5 & Group2$pvalue < 0.05] <- "DOWN"
##Group3
##Divide among up and down regulated genes
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
Group3$diffexpressed <- "NO"
# if log2Foldchange > 2.5 and pvalue < 0.05, set as "UP"
Group3$diffexpressed[Group3$log2FoldChange > 2.5 & Group3$pvalue < 0.05] <- "UP"
# if log2Foldchange < -2.5 and pvalue < 0.05, set as "DOWN"
Group3$diffexpressed[Group3$log2FoldChange < -2.5 & Group3$pvalue < 0.05] <- "DOWN"
##Group4
##Divide among up and down regulated genes
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
Group4$diffexpressed <- "NO"
# if log2Foldchange > 2.5 and pvalue < 0.05, set as "UP"
Group4$diffexpressed[Group4$log2FoldChange > 2.5 & Group4$pvalue < 0.05] <- "UP"
# if log2Foldchange < -2.5 and pvalue < 0.05, set as "DOWN"
Group4$diffexpressed[Group4$log2FoldChange < -2.5 & Group4$pvalue < 0.05] <- "DOWN"

g1_up <- Group1 %>% 
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::select(row.names)

g1_down<- Group1 %>% 
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::select(row.names)

g2_up <- Group2 %>% 
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::select(row.names)

g2_down<- Group2 %>% 
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::select(row.names)

g3_up <- Group3 %>% 
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::select(row.names)

g3_down<- Group3 %>% 
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::select(row.names)

g4_up <- Group4 %>% 
  dplyr::filter(diffexpressed == "UP")%>%
  dplyr::select(row.names)

g4_down<- Group4 %>% 
  dplyr::filter(diffexpressed == "DOWN")%>%
  dplyr::select(row.names)

##
Group1=Group1[,1]
Group2=Group2[,1]
Group3=Group3[,1]
Group4=Group4[,1]
##
g1_up=g1_up[,1]
g1_down=g1_down[,1]
g2_up=g2_up[,1]
g2_down=g2_down[,1]
g3_up=g3_up[,1]
g3_down=g3_down[,1]
g4_up=g4_up[,1]
g4_down=g4_down[,1]


###GROUP1X2 - UP AND DOWN
x <- list(g1_up,g1_down, g2_up, g2_down)
cols <- c("#004599", "#00c2cc", "#e69138", "#e2eb32") 
names(x) <- c("Group1 UP","Group1 DOWN","Group2 UP", "Group2 DOWN")
g<-ggvenn(
  x, 
  fill_color = cols,
  stroke_size = 0.5, set_name_size = 4) 
lg <-legendGrob(labels=c("Group 1 UP = HD ANxIN UP", "Group 1 DOWN = HD ANxIN DOWN",
                         "Group 2 UP = CVID ANxIN UP", "Group 2 DOWN = CVID ANxIN DOWN"),
                pch=rep(19),
                gp=gpar(col=cols, fill="gray", cex=0.6),
                byrow=TRUE)
gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))

###GROUP3X4 - UP AND DOWN
x <- list(g3_up,g3_down, g4_up, g4_down)
cols <- c("#004599", "#00c2cc", "#e69138", "#e2eb32") 
names(x) <- c("Group3 UP","Group3 DOWN","Group4 UP", "Group4 DOWN")
g<-ggvenn(
  x, 
  fill_color = cols,
  stroke_size = 0.5, set_name_size = 4) 
lg <-legendGrob(labels=c("Group 3 UP = HD ABxIB UP", "Group 3 DOWN = HD ABxIB DOWN",
                         "Group 4 UP = CVID ABxIB UP", "Group 4 DOWN = CVID ABxIB DOWN"),
                pch=rep(19),
                gp=gpar(col=cols, fill="gray", cex=0.6),
                byrow=TRUE)
gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))

###SIMPLE COMPARISONS
x <- list(Group1, Group2)
cols <- c("#00c2cc", "#cc0a00") 
names(x) <- c("Group 1","Group 2")
g<-ggvenn(
  x, 
  fill_color = cols,
  stroke_size = 0.5, set_name_size = 4)
lg <-legendGrob(labels=c("Group 1 = HD ANxIN","Group 2 = CVID ANxIN"), pch=rep(19),
                gp=gpar(col=cols, fill="gray"),
                byrow=TRUE)
gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))


x <- list(Group3, Group4)
cols <- c("#540099", "#d84000") 
names(x) <- c("Group 3", "Group 4")
g<-ggvenn(
  x, 
  fill_color = cols,
  stroke_size = 0.5, set_name_size = 4)
lg <-legendGrob(labels=c("Group 3 = HD ABxIB", "Group 4 = CVID ABXIB"), pch=rep(19),
                gp=gpar(col=cols, fill="gray"),
                byrow=TRUE)
gridExtra::grid.arrange(g, lg, ncol = 2, widths = c(4,1))

#extract the genes from the Venn diagram
tmp=venn(x,showSetLogicLabel=TRUE)
isect <- attr(tmp, "intersection")
isect
df<- data.frame(isect[["Group 3:Group 4"]])
df<-write.xlsx(df,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/inters_3x4.xlsx", asTable = TRUE)
df2<- data.frame(isect["Group 3"])
df2<-write.xlsx(df2,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/group3.xlsx", asTable = TRUE)
df3<- data.frame(isect["Group 4"])
df3<-write.xlsx(df3,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/group4.xlsx", asTable = TRUE)


#extrac all the intersections
#1x2
tmp=venn(x,showSetLogicLabel=TRUE)
isect <- attr(tmp, "intersection")
isect
df<- data.frame(isect[["Group1 DOWN:Group2 DOWN"]])
df<-write.xlsx(df,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/1x2_down.xlsx", asTable = TRUE)
df2<- data.frame(isect["Group1 DOWN"])
df2<-write.xlsx(df2,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/1_down.xlsx", asTable = TRUE)
df3<- data.frame(isect["Group2 DOWN"])
df3<-write.xlsx(df3,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/2_down.xlsx", asTable = TRUE)
df4<- data.frame(isect[["Group1 UP:Group2 UP"]])
df4<-write.xlsx(df4,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/1x2_up.xlsx", asTable = TRUE)
df5<- data.frame(isect["Group1 UP"])
df5<-write.xlsx(df5,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/1_up.xlsx", asTable = TRUE)
df6<- data.frame(isect["Group2 UP"])
df6<-write.xlsx(df6,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/2_up.xlsx", asTable = TRUE)

#3x4
tmp=venn(x,showSetLogicLabel=TRUE)
isect <- attr(tmp, "intersection")
isect
df<- data.frame(isect[["Group3 DOWN:Group4 DOWN"]])
df<-write.xlsx(df,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/3x4_down.xlsx", asTable = TRUE)
df2<- data.frame(isect["Group3 DOWN"])
df2<-write.xlsx(df2,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/3_down.xlsx", asTable = TRUE)
df3<- data.frame(isect["Group4 DOWN"])
df3<-write.xlsx(df3,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/4_down.xlsx", asTable = TRUE)
df4<- data.frame(isect[["Group3 UP:Group4 UP"]])
df4<-write.xlsx(df4,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/3x4_up.xlsx", asTable = TRUE)
df5<- data.frame(isect["Group3 UP"])
df5<-write.xlsx(df5,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/3_up.xlsx", asTable = TRUE)
df6<- data.frame(isect["Group4 UP"])
df6<-write.xlsx(df6,"C:/Users/Anna/Documents/CVID_project/Venn/INTERS/4_up.xlsx", asTable = TRUE)
