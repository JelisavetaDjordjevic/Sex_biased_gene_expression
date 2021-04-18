###load packages
library(readr)
library(edgeR)
library(RColorBrewer)
library(ggplot2)

getwd()

#setwd

#read in Drosophila read counts "droso_counts.csv"
droso_counts<-read_csv("droso_counts.csv")
colnames(droso_counts)

##add factor: sex (F=female, M= male), stage (L=larva, P=pupa, A=adult)
sex<- c("F","F","F","F","F","F","F","M","M","F","F","M","M","M","M","M","M","F","M","M","F","M","M","F","M")
stage<- c("P","P","P","A","L","A","L","A","A","A","L","P","L","P","P","L","A","L","A","L","A","L","P","P","A")
sex_stage<-paste(sex, stage, sep="_")


##make edgeR object
x<-DGEList(counts = droso_counts[, 2:26],
           genes = droso_counts[, 1], 
           group = sex_stage)
x$samples
table(x$samples$group)
summary(x$samples$lib.size)

##exploring log2 data before normalization 
logcounts <- cpm(x, prior.count = 2, log=TRUE)
hist(logcounts[,4])

#boxplot Log2CPM
boxplot(logcounts, xlab="",
        ylab="Log2 counts per million",
        las=2, 
        cex.axis=0.6)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

#MDS plot all samples
fac<- factor(sex_stage)
colours<- brewer.pal(nlevels(fac), "Paired")
name<- c("FP.H8.2", "FP.H94.3"  ,"FP.H75.3"  ,"FL.H169.1" ,"FA.H184.3", "ML.H8.1","FL.H60.3","FA.H8.3", "MP.H163.2", "MA.H94.1" , "MA.H55.3" , "FL.H74.3" , "FA.H74.3"  ,"MA.H74.2" , "MP.H184.3","FP.H163.1", "ML.H75.2",  "ML.H55.3",  "FA.H55.2"  ,"MA.H75.1",  "FL.H94.2",  "MA.H75.3",  "ML.H60.2","MP.H169.3", "MP.H196.3")
plotMDS(logcounts, labels= name, col= colours[as.numeric(fac)])
plotMDS(logcounts, labels= factor(x$samples$group), col= colours[as.numeric(fac)])

#subseting data to analyse per stage#
### Adult ###
counts_adult<- droso_counts[, c(5,7,9,10,11,18,20,22,26)]
colnames(counts_adult)
counts_adult$Gene_name<- droso_counts$Gene_name
sex_a<-c("F","F","M","M","F", "M", "M", "F", "M") 
x_a_stage<-DGEList(counts = droso_counts[, c(5,7,9,10,11,18,20,22,26)],
                   genes = droso_counts[, 1],
                   group = sex_a)
x_a_stage$samples

#normalization
cpm(10, mean(x_a_stage$samples$lib.size))
x_a_stage
head(str(x_a_stage))
keep_a <- rowSums(cpm(x_a_stage[,c(1,2,5,8)])>0.5) >=3 | 
  rowSums(cpm(x_a_stage[,c(3,4,6,7,9)])>0.5) >=3  
x_a_stage<- x_a_stage[keep_a, , keep.lib.sizes=FALSE]
table(keep_a)
x_a_stage<- calcNormFactors(x_a_stage)
x_a_stage$samples

#MDS plot of adult samples
fac_a<- factor(sex_a)
colours<- brewer.pal(nlevels(fac_a), "Paired")
plotMDS(x_a_stage, labels= factor(x_a_stage$samples$group), col= colours[as.numeric(fac_a)])

#design a model
design_a <- model.matrix(~ 0 + group, data= x_a_stage$samples)
colnames(design_a)<- levels(x_a_stage$samples$group)
colnames(design_a)
#estimate dispersion
x_a_stage <- estimateDisp(x_a_stage, design_a, robust=TRUE)
summary(x_a_stage$trended.dispersion)
sqrt(x_a_stage$common.dispersion) #this is BCV value

#plot biological coef. of variation
plotBCV(x_a_stage)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_a <- glmQLFit(x_a_stage, design_a, robust=TRUE)
head(fit_a$coefficients)
plotQLDisp(fit_a)
summary(fit_a$df.prior)

#between m and F
FM_a<- makeContrasts(F - M, levels=design_a)
res_FM_a <- glmQLFTest(fit_a, contrast=FM_a)
is.de_FM_a <- decideTestsDGE(res_FM_a, p.value=0.05)
#number of sex-biased genes
summary(is.de_FM_a)

#save top tags in a table for n genes
adult_top_fm<- topTags(res_FM_a, n=14000)
adult_top_fm<-adult_top_fm$table


plotMD(res_FM_a, 
       status=is.de_FM_a,
       values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15),
       main = "DGE between sexual females and males at pupal stage",
       hl.cex=0.7 )



###calculate the average logCPM 

##To make classify genes into sex-bias categories
logcounts_a <- cpm(x_a_stage,log=TRUE)
colnames(logcounts_a)<- x_a_stage$samples$group
rownames(logcounts_a)<- x_a_stage$genes$Gene_name

names<-as.character(adult_top_fm$Gene_name)
logcounts_a<-cbind(rownames(logcounts_a), logcounts_a)
colnames(logcounts_a)[1]<-"Gene_name"
logcounts_a<-as.data.frame(logcounts_a)

logcounts_a_sub<-logcounts_a[logcounts_a$Gene_name %in% names,]
logcounts_a_sub[,-1]<-as.numeric(as.character(unlist(logcounts_a_sub[,-1])))
average_M<-rowMeans(logcounts_a_sub[,c(4,5,7,8,10)])
average_F<-rowMeans(logcounts_a_sub[,c(2,3,6,9)])
logcounts_a_sub<-cbind(average_F, logcounts_a_sub)
logcounts_a_sub<-cbind(average_M, logcounts_a_sub)

adult_top_fm<-adult_top_fm[order(adult_top_fm$Gene_name),]
logcounts_a_sub<- logcounts_a_sub[order(logcounts_a_sub$Gene_name),]
adult_top_fm<- cbind(logcounts_a_sub$average_M, adult_top_fm)
adult_top_fm<- cbind(logcounts_a_sub$average_F, adult_top_fm)
colnames(adult_top_fm)[1]<- "Average_logCPM_F"
colnames(adult_top_fm)[2]<- "Average_logCPM_M"

## calculating average CPM for males and females and binding it to TOP_fm
counts_a <- cpm(x_a_stage, log = FALSE)
rownames(counts_a)<-x_a_stage$genes$Gene_name
colnames(counts_a)<- x_a_stage$samples$group

counts_a<-cbind(rownames(counts_a), counts_a)
colnames(counts_a)[1]<-"Gene_name"
counts_a<-as.data.frame(counts_a)
counts_a_sub<-counts_a[counts_a$Gene_name %in% names,]
counts_a_sub[,-1]<-as.numeric(as.character(unlist(counts_a_sub[,-1])))
average_cpm_M<-rowMeans(counts_a_sub[,c(4,5,7,8,10)])
average_cpm_F<-rowMeans(counts_a_sub[,c(2,3,6,9)])
counts_a_sub<-cbind(average_cpm_F, counts_a_sub)
counts_a_sub<-cbind(average_cpm_M, counts_a_sub)

adult_top_fm<-adult_top_fm[order(adult_top_fm$Gene_name),]
counts_a_sub<- counts_a_sub[order(counts_a_sub$Gene_name),]
adult_top_fm<- cbind(counts_a_sub$average_cpm_M, adult_top_fm)
adult_top_fm<- cbind(counts_a_sub$average_cpm_F, adult_top_fm)
colnames(adult_top_fm)[1]<- "Average_CPM_F"
colnames(adult_top_fm)[2]<- "Average_CPM_M"

#gene classification#
adult_top_fm$fact<-rep("NA", nrow(adult_top_fm))
adult_top_fm$fact[adult_top_fm$Average_CPM_M ==(0) & adult_top_fm$logFC>2] <- "F_limited"
adult_top_fm$fact[adult_top_fm$Average_CPM_F ==(0) & adult_top_fm$logFC< (-2)] <- "M_limited"
adult_top_fm$fact[adult_top_fm$logFC>2 & adult_top_fm$Average_CPM_M > (0)] <- "F_biased"
adult_top_fm$fact[adult_top_fm$logFC<(-2) & adult_top_fm$Average_CPM_F > (0)]<- "M_biased"
adult_top_fm$fact[adult_top_fm$logFC>0 & adult_top_fm$logFC<2]<- "Any_f_biased"
adult_top_fm$fact[adult_top_fm$logFC<0 & adult_top_fm$logFC>-2]<- "Any_m_biased"
adult_top_fm$fact[adult_top_fm$FDR >0.05]<- "Not DE"

####### plot
palette(c("palevioletred1","slategray1","firebrick1","red4","lightslateblue", "navy","gainsboro"))
table(adult_top_fm$fact)
adult_top_fm$fact<-as.factor(adult_top_fm$fact)
unbias<-subset(adult_top_fm, fact == "Not DE")
bias<-subset(adult_top_fm, fact != "Not DE")
table(bias$fact)
plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", 
     pch=19, 
     ylim = c(-3,17), 
     xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="Gene expression in males and females at the pupal stage", 
     cex=0.8)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), pch=19, cex=0.7)

legend("topleft", c("Slightly FB","slightly MB","Strongly FB", "Female limited","Male biased", "Male limited", "Not DE"), fill=c(1:7), cex = 0.6,bty = "n")

