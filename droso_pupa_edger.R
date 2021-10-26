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
x<-DGEList(counts = droso_counts[, 2:26],genes = droso_counts[, 1], group = sex_stage)
x$samples
table(x$samples$group)
summary(x$samples$lib.size)

##exploring log2 data before normalization 
logcounts <- cpm(x, prior.count = 2, log=TRUE)
hist(logcounts[,4])
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2, cex.axis=0.6)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs")

library(RColorBrewer)
fac<- factor(sex_stage)
colours<- brewer.pal(nlevels(fac), "Paired")
names_samples<- c("FP.H94.3","FP.H8.2","FP.H75.1","FA.H184.3","FL.H169.1","FA.H74.3","FL.H74.3","MA.H94.1","MA.H55.3","FA.H8.3","FL.H60.3","MP.H163.2","ML.H8.1","MP.H74.1","MP.H196.3","ML.H60.2","MA.H75.3","FL.H94.2","MA.H75.1","ML.H55.3","FA.H55.2","ML.H75.2","MP.H184.3","FP.H163.1","MA.H74.2")
plotMDS(logcounts, labels= factor(x$samples$group), col= colours[as.numeric(fac)])
plotMDS(logcounts, labels= factor(names_samples), col= colours[as.numeric(fac)])

#subseting data to analyse per stage#
###pupae###
counts_pupae<- droso_counts[, c(2,3,4,13,15,16,24,25)]
colnames(counts_pupae)
counts_pupae$Gene_name<- droso_counts$Gene_name
sex_p<-c("F","F","F","M","M","M","M","F") 

#make edgeR object for pupal stage
x_p_stage<-DGEList(counts = droso_counts[, c(2,3,4,13,15,16,24,25)],
                   genes = droso_counts[, 1],
                   group = sex_p)
x_p_stage$samples

#normalization
cpm(10, mean(x_p_stage$samples$lib.size))
x_p_stage
head(str(x_p_stage))

#filter out genes that are not expressed in at least 3 male or 3 female libraries 
keep_p <- rowSums(cpm(x_p_stage[,c(1,2,3,8)])>0.5) >=3 | 
        rowSums(cpm(x_p_stage[,c(4,5,6,7)])>0.5) >=3  
x_p_stage<- x_p_stage[keep_p, , keep.lib.sizes=FALSE]

#show number of genes kept and filtered out
table(keep_p)

#calculate normalization factor
x_p_stage<- calcNormFactors(x_p_stage)
x_p_stage$samples

fac_p<- factor(sex_p)
colours<- brewer.pal(nlevels(fac_p), "Paired")
plotMDS(x_p_stage, labels= factor(x_p_stage$samples$group), col= colours[as.numeric(fac_p)])

#design a model
design_p <- model.matrix(~ 0 + group, data= x_p_stage$samples)
colnames(design_p)<- levels(x_p_stage$samples$group)
colnames(design_p)

#estimate dispersion
x_p_stage <- estimateDisp(x_p_stage, design_p, robust=TRUE)
x_p_stage$common.dispersion

#plot Biologicaé coeff. of variation
plotBCV(x_p_stage)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_p <- glmQLFit(x_p_stage, design_p, robust=TRUE)
head(fit_p$coefficients)
plotQLDisp(fit_p)
summary(fit_p$df.prior)

#between m and F
FM_p<- makeContrasts(F - M, levels=design_p)

#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
##Conduct genewise statistical tests for a given coefficient or contrast.
res_FM_p <- glmQLFTest(fit_p, contrast=FM_p)

#Identify which genes are significantly differentially expressed from an edgeR fit object
is.de_FM_p <- decideTestsDGE(res_FM_p, p.value=0.05)

#number of Sex-biased genes
summary(is.de_FM_p)

#get fc values and p-values for (n) genes
pupae_top_fm<- topTags(res_FM_p, n=12000)
pupae_top_fm<-pupae_top_fm$table

plotMD(res_FM_p, 
       status=is.de_FM_p,
       values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15),
       main = "DGE between females and males at pupal stage", 
       hl.cex=0.7 )

#plot different classes of sex-biased genes
###calculate the average logCPM for males and females separtely

##log counts for pupal stage
logcounts_p <- cpm(x_p_stage,log=TRUE)
colnames(logcounts_p)<- x_p_stage$samples$group
rownames(logcounts_p)<- x_p_stage$genes$Gene_name

names<-as.character(pupae_top_fm$Gene_name)
logcounts_p<-cbind(rownames(logcounts_p), logcounts_p)
colnames(logcounts_p)[1]<-"Gene_name"
logcounts_p<-as.data.frame(logcounts_p)
logcounts_p_sub<-logcounts_p[logcounts_p$Gene_name %in% names,]
logcounts_p_sub[,-1]<-as.numeric(as.character(unlist(logcounts_p_sub[,-1])))
average_M<-rowMeans(logcounts_p_sub[,c(5,6,7,8)])
average_F<-rowMeans(logcounts_p_sub[,c(2,3,4,9)])
logcounts_p_sub<-cbind(average_F, logcounts_p_sub)
logcounts_p_sub<-cbind(average_M, logcounts_p_sub)

pupae_top_fm<-pupae_top_fm[order(pupae_top_fm$Gene_name),]
logcounts_p_sub<- logcounts_p_sub[order(logcounts_p_sub$Gene_name),]
pupae_top_fm<- cbind(logcounts_p_sub$average_M, pupae_top_fm)
pupae_top_fm<- cbind(logcounts_p_sub$average_F, pupae_top_fm)
colnames(pupae_top_fm)[1]<- "Average_logCPM_F"
colnames(pupae_top_fm)[2]<- "Average_logCPM_M"

## calculating average CPM for males and females and binding it to TOP_fm
counts_p <- cpm(x_p_stage, log = FALSE)
rownames(counts_p)<-x_p_stage$genes$Gene_name
colnames(counts_p)<- x_p_stage$samples$group

counts_p<-cbind(rownames(counts_p), counts_p)
colnames(counts_p)[1]<-"Gene_name"
counts_p<-as.data.frame(counts_p)
counts_p_sub<-counts_p[counts_p$Gene_name %in% names,]
counts_p_sub[,-1]<-as.numeric(as.character(unlist(counts_p_sub[,-1])))
average_cpm_M<-rowMeans(counts_p_sub[,c(5,6,7,8)])
average_cpm_F<-rowMeans(counts_p_sub[,c(2,3,4,9)])
counts_p_sub<-cbind(average_cpm_F, counts_p_sub)
counts_p_sub<-cbind(average_cpm_M, counts_p_sub)

pupae_top_fm<-pupae_top_fm[order(pupae_top_fm$Gene_name),]
counts_p_sub<- counts_p_sub[order(counts_p_sub$Gene_name),]
pupae_top_fm<- cbind(counts_p_sub$average_cpm_M, pupae_top_fm)
pupae_top_fm<- cbind(counts_p_sub$average_cpm_F, pupae_top_fm)
colnames(pupae_top_fm)[1]<- "Average_CPM_F"
colnames(pupae_top_fm)[2]<- "Average_CPM_M"

#gene classification#
pupae_top_fm$fact<-rep("NA", nrow(pupae_top_fm))
pupae_top_fm$fact[pupae_top_fm$Average_CPM_M ==(0) & pupae_top_fm$logFC>1] <- "F_limited"
pupae_top_fm$fact[pupae_top_fm$Average_CPM_F ==(0) & pupae_top_fm$logFC< (-1)] <- "M_limited"
pupae_top_fm$fact[pupae_top_fm$logFC>1 & pupae_top_fm$Average_CPM_M > (0)] <- "F_biased"
pupae_top_fm$fact[pupae_top_fm$logFC<(-1) & pupae_top_fm$Average_CPM_F > (0)]<- "M_biased"
pupae_top_fm$fact[pupae_top_fm$logFC>0 & pupae_top_fm$logFC<1]<- "Any_f_biased"
pupae_top_fm$fact[pupae_top_fm$logFC<0 & pupae_top_fm$logFC>-1]<- "Any_m_biased"
pupae_top_fm$fact[pupae_top_fm$FDR >0.05]<- "Not DE"

####### plot
palette(c("palevioletred1","slategray1","firebrick1","lightslateblue", "navy","gainsboro"))
table(pupae_top_fm$fact)
pupae_top_fm$fact<-as.factor(pupae_top_fm$fact)
unbias<-subset(pupae_top_fm, fact == "Not DE")
bias<-subset(pupae_top_fm, fact != "Not DE")
table(bias$fact)
plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", 
     pch=19, 
     ylim = c(-3,17), xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="Gene expression in males and females at the pupal stage", cex=0.8)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), pch=19, cex=0.7)


legend("topleft", c("Slightly FB","slightly MB" , "Strongly FB", "Female limited","Male biased", "Male limited", "Not DE"), fill=c(1:7), cex = 0.6,bty = "n")
