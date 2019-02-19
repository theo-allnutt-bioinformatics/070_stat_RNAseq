#!/usr/local/bin/Rscript

#EdgeR

stat vs wt
library("locfit")
library("edgeR")
x <- read.delim("counts.tab",row.names="geneID")
group <- factor(c(1,1,1,1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tt<-topTags(qlf,n='all',sort.by="logFC")
write.table(as.data.frame(tt),file="statVwt.edger",sep="\t")

fdr<0.05 = 135 genes


#stat1 vs wt
#cut -f1,2,3,4,8,9,10 counts.tab > stat1.tab

library("locfit")
library("edgeR")
x <- read.delim("stat1.tab",row.names="geneID")
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tt<-topTags(qlf,n='all',sort.by="logFC")
write.table(as.data.frame(tt),file="stat1.edger",sep="\t")

#stat2 vs wt
#cut -f1,5,6,7,8,9,10 counts.tab > stat2.tab

library("locfit")
library("edgeR")
x <- read.delim("stat2.tab",row.names="geneID")
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tt<-topTags(qlf,n='all',sort.by="logFC")
write.table(as.data.frame(tt),file="stat2.edger",sep="\t")

#stat1 vs stat2
#cut -f1,2,3,4,5,6,7 counts.tab > stat1v2.tab

library("locfit")
library("edgeR")
x <- read.delim("stat1v2.tab",row.names="geneID")
group <- factor(c(1,1,1,2,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
tt<-topTags(qlf,n='all',sort.by="logFC")
write.table(as.data.frame(tt),file="stat1v2.edger",sep="\t")


#ALDEX2
#source("https://bioconductor.org/biocLite.R")
#biocLite("ALDEx2")
#devtools::install_github("Bioconductor-mirror/ALDEx2")

library("ALDEx2")
x <- read.delim("counts.tab",row.names="geneID")
group <-c(1,1,1,1,1,1,2,2,2)
tt <- aldex(as.data.frame(x), group,test = "t", mc.samples = 128)
aldex.plot(tt, type = "MW", cutoff = .05 / ncol(x))
write.table(tt,"aldex_stat.txt",sep="\t")

#aldex = 45 genes we.eBH <0.05

