###packages
library(readr)
library(edgeR)
library(RColorBrewer)
library(ggplot2)

getwd()
#setwd

#read in Drosophila read counts "droso_counts.csv"
droso_counts<-read_csv("droso_counts.csv")
colnames(droso_counts)

#add factor: sex (F=female, M= male), stage (L=larva, P=pupa, A=adult)
sex<- c("F","F","F","F","F","F","F","M","M","F","F","M","M","M","M","M","M","F","M","M","F","M","M","F","M")
stage<- c("P","P","P","A","L","A","L","A","A","A","L","P","L","P","P","L","A","L","A","L","A","L","P","P","A")
sex_stage<-paste(sex, stage, sep="_")

#make edgeR object
x<-DGEList(counts = droso_counts[, 2:26],
           genes = droso_counts[, 1], 
           group = sex_stage)
x$samples
table(x$samples$group)
summary(x$samples$lib.size)

##exploring log2 data before normalization 
logcounts <- cpm(x, prior.count = 2, log=TRUE)
hist(logcounts[,4])
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, cex.axis=0.6)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

#MDS plot
fac<- factor(sex_stage)
colours<- brewer.pal(nlevels(fac), "Paired")
names_samples<- c("FP.H94.3","FP.H8.2","FP.H75.1","FA.H184.3","FL.H169.1","FA.H74.3","FL.H74.3","MA.H94.1","MA.H55.3","FA.H8.3","FL.H60.3","MP.H163.2","ML.H8.1","MP.H74.1","MP.H196.3","ML.H60.2","MA.H75.3","FL.H94.2","MA.H75.1","ML.H55.3","FA.H55.2","ML.H75.2","MP.H184.3","FP.H163.1","MA.H74.2")
plotMDS(logcounts, 
        labels= factor(x$samples$group), 
        col= colours[as.numeric(fac)])
plotMDS(logcounts, 
        labels= factor(names_samples), 
        col= colours[as.numeric(fac)])

#subseting data to analyse per stage#
###3rd instar Larva###
counts_larvae<- droso_counts[, c(6,8,12,14,17,19,21,23)]
colnames(counts_larvae)
counts_larvae$Gene_name<- droso_counts$Gene_name
sex_l<-c("F","F","F","M","M","F","M","M") 
#make edgeR object for larval stage
x_l_stage<-DGEList(counts = droso_counts[, c(6,8,12,14,17,19,21,23)],
                   genes = droso_counts[, 1],
                   group = sex_l)
x_l_stage$samples

#normalization for library sizes
cpm(10, mean(x_l_stage$samples$lib.size))
x_l_stage
head(str(x_l_stage))
#filter out genes that are not expressed in at least 3 male or 3 female libraries 
keep_l <- rowSums(cpm(x_l_stage[,c(1,2,3,6)])>0.5) >=3 | 
  rowSums(cpm(x_l_stage[,c(4,5,7,8)])>0.5) >=3  

x_l_stage<- x_l_stage[keep_l, , keep.lib.sizes=FALSE]

#show number of genes kept and filtered out
table(keep_l)
#calculate normalization factor
x_l_stage<- calcNormFactors(x_l_stage)
x_l_stage$samples

fac_l<- factor(sex_l)
colours<- brewer.pal(nlevels(fac_l), "Paired")

plotMDS(x_l_stage, 
        labels= factor(x_l_stage$samples$group), 
        col= colours[as.numeric(fac_l)])

#design a model
design_l <- model.matrix(~ 0 + group, data= x_l_stage$samples)
colnames(design_l)<- levels(x_l_stage$samples$group)
colnames(design_l)

#estimate dispersion
x_l_stage <- estimateDisp(x_l_stage, design_l, robust=TRUE)
summary(x_l_stage$trended.dispersion)
x_l_stage$common.dispersion

#Biological coef. variation plot
plotBCV(x_l_stage)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_l <- glmQLFit(x_l_stage, design_l, robust=TRUE)
head(fit_l$coefficients)
plotQLDisp(fit_l)
summary(fit_l$df.prior)

#make contrast between m and F
FM_l<- makeContrasts(F - M, levels=design_l)

#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
##Conduct genewise statistical tests for a given coefficient or contrast.
res_FM_l <- glmQLFTest(fit_l, contrast=FM_l)

#Identify which genes are significantly differentially expressed from an edgeR fit object
is.de_FM_l <- decideTestsDGE(res_FM_l, p.value=0.05)

#number of Sex-biased genes 
summary(is.de_FM_l)

#get fc values and p-values for (n) genes
larvae_top_fm<- topTags(res_FM_l, n=12000)
larvae_top_fm<-larvae_top_fm$table

plotMD(res_FM_l, 
       status=is.de_FM_l,
       values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15), 
       main = "DGE between females and males at larval stage", 
       hl.cex=0.7 )


###calculate the average logCPM for males and females separtely
##log counts for larval stage
logcounts_l <- cpm(x_l_stage,log=TRUE)
colnames(logcounts_l)<- x_l_stage$samples$group
rownames(logcounts_l)<- x_l_stage$genes$Gene_name

names<-as.character(larvae_top_fm$Gene_name)
logcounts_l<-cbind(rownames(logcounts_l), logcounts_l)
colnames(logcounts_l)[1]<-"Gene_name"
logcounts_l<-as.data.frame(logcounts_l)
logcounts_l_sub<-logcounts_l[logcounts_l$Gene_name %in% names,]
logcounts_l_sub[,-1]<-as.numeric(as.character(unlist(logcounts_l_sub[,-1])))
average_M<-rowMeans(logcounts_l_sub[,c(5,6,8,9)])
average_F<-rowMeans(logcounts_l_sub[,c(2,3,4,7)])
logcounts_l_sub<-cbind(average_F, logcounts_l_sub)
logcounts_l_sub<-cbind(average_M, logcounts_l_sub)

larvae_top_fm<-larvae_top_fm[order(larvae_top_fm$Gene_name),]
logcounts_l_sub<- logcounts_l_sub[order(logcounts_l_sub$Gene_name),]
larvae_top_fm<- cbind(logcounts_l_sub$average_M, larvae_top_fm)
larvae_top_fm<- cbind(logcounts_l_sub$average_F, larvae_top_fm)
colnames(larvae_top_fm)[1]<- "Average_logCPM_F"
colnames(larvae_top_fm)[2]<- "Average_logCPM_M"

## calculating average CPM for males and females and binding it to TOP_fm
counts_l <- cpm(x_l_stage, log = FALSE)
rownames(counts_l)<-x_l_stage$genes$Gene_name
colnames(counts_l)<- x_l_stage$samples$group

counts_l<-cbind(rownames(counts_l), counts_l)
colnames(counts_l)[1]<-"Gene_name"
counts_l<-as.data.frame(counts_l)
counts_l_sub<-counts_l[counts_l$Gene_name %in% names,]
counts_l_sub[,-1]<-as.numeric(as.character(unlist(counts_l_sub[,-1])))
average_cpm_M<-rowMeans(counts_l_sub[,c(5,6,8,9)])
average_cpm_F<-rowMeans(counts_l_sub[,c(2,3,4,7)])
counts_l_sub<-cbind(average_cpm_F, counts_l_sub)
counts_l_sub<-cbind(average_cpm_M, counts_l_sub)

larvae_top_fm<-larvae_top_fm[order(larvae_top_fm$Gene_name),]
counts_l_sub<- counts_l_sub[order(counts_l_sub$Gene_name),]
larvae_top_fm<- cbind(counts_l_sub$average_cpm_M, larvae_top_fm)
larvae_top_fm<- cbind(counts_l_sub$average_cpm_F, larvae_top_fm)
colnames(larvae_top_fm)[1]<- "Average_CPM_F"
colnames(larvae_top_fm)[2]<- "Average_CPM_M"


#gene classification#
#expression_type<- if(hatch_top_fm$Average_logCPM_M <(-3.5) == "F_limited" | hatch_top_fm$Average_logCPM_F <(-3.5)=="M_limited") else(hatch_top_fm$Average_logCPM_M >(-3.5) =="F_biased")
larvae_top_fm$fact<-rep("NA", nrow(larvae_top_fm))

#gene classification#
larvae_top_fm$fact<-rep("NA", nrow(larvae_top_fm))
larvae_top_fm$fact[larvae_top_fm$Average_CPM_M ==(0) & larvae_top_fm$logFC>1] <- "F_limited"
larvae_top_fm$fact[larvae_top_fm$Average_CPM_F ==(0) & larvae_top_fm$logFC< (-1)] <- "M_limited"
larvae_top_fm$fact[larvae_top_fm$logFC>1 & larvae_top_fm$Average_CPM_M > (0)] <- "F_biased"
larvae_top_fm$fact[larvae_top_fm$logFC<(-1) & larvae_top_fm$Average_CPM_F > (0)]<- "M_biased"
larvae_top_fm$fact[larvae_top_fm$logFC>0 & larvae_top_fm$logFC<1]<- "Any_f_biased"
larvae_top_fm$fact[larvae_top_fm$logFC<0 & larvae_top_fm$logFC>-1]<- "Any_m_biased"
larvae_top_fm$fact[larvae_top_fm$FDR >0.05]<- "Not DE"

####### plot
palette(c("palevioletred1","slategray1","firebrick1","lightslateblue", "navy","gainsboro"))

table(larvae_top_fm$fact)
larvae_top_fm$fact<-as.factor(larvae_top_fm$fact)
unbias<-subset(larvae_top_fm, fact == "Not DE")
bias<-subset(larvae_top_fm, fact != "Not DE")
table(bias$fact)
plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", pch=19, 
     ylim = c(-3,17), 
     xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="Gene expression in males and females at the larval stage", 
     cex=0.8)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)

legend("topleft", c("Slightly FB","slightly MB" , "Strongly FB", "Female limited","Male biased", "Male limited", "Not DE"), fill=c(1:7), cex = 0.6,bty = "n")
