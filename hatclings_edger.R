#Hathcling sex biased gene expression
###packages
library(readr)
library(edgeR)
library(RColorBrewer)
getwd()

#import timema_counts
tcm_genome <- read.csv("/Users/jdjordje/Desktop/Jeli/job/genome_analysis/R_scripts/timema_counts.csv")
colnames(tcm_genome)

sex<- c("F","M","M","F","M","F","M","M","M","M","M",
        "F","F","F","F","F","M","M","M","M","M","F","F","F")
stage<- c("H","H","A","H","A","H","H","H","A","A","H","A","H",
          "A","A","A","H","J","J","J","J","J","J","J")
sex_stage<-paste(sex, stage, sep="_")
colnames(tcm_genome)


#subseting data to analyse per stage#
###HATCHLING STAGE###
counts_hatch<- tcm_genome[, c(2,3,5,7,8,9,12,14,18)]
colnames(counts_hatch)
counts_hatch$Gene_name<- tcm_genome$Gene_name
#corresponding sex
sex_h<-c("F","M","F","F","M","M","M","F","M") 

#get egder object DGElist
x_h_stage<-DGEList(counts = tcm_genome[, c(2,3,5,7,8,9,12,14,18)],
                   genes = tcm_genome[, 1],group = sex_h)
x_h_stage$samples

#normalization
#cpm that correspond to the count of 10 for the library size
cpm(10, mean(x_h_stage$samples$lib.size))
head(str(x_h_stage))
#filtering lowly expressed genes
keep_h <- rowSums(cpm(x_h_stage[,c(1,3,4,8)])>0.5) >=3 | 
        rowSums(cpm(x_h_stage[,c(2,5,6,7,9)])>0.5) >=3  

x_h_stage<- x_h_stage[keep_h, , keep.lib.sizes=FALSE]

#show number of genes that are kept for the analysis
table(keep_h)
#calculate normalization factor to scale library sizes
x_h_stage<- calcNormFactors(x_h_stage)
x_h_stage$samples

#mds plot- hatchling
fac_h<- factor(sex_h)
colours<- brewer.pal(nlevels(fac_h), "Paired")
colours<- c("red", "blue")
plotMDS(x_h_stage, 
        labels= factor(x_h_stage$samples$group), 
        col= colours[as.numeric(fac_h)])

#design a model
design_h <- model.matrix(~ 0 + group, data= x_h_stage$samples)
colnames(design_h)<- levels(x_h_stage$samples$group)
colnames(design_h)

#estimate dispersion
x_h_stage <- estimateDisp(x_h_stage, design_h, robust=TRUE)
x_h_stage$trended.dispersion
x_h_stage$common.dispersion
plotBCV(x_h_stage)

####Genewise Negative Binomial Generalized Linear Models with Quasi-likelihood Tests
#Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 
#Conduct genewise statistical tests for a given coefficient or contrast.
fit_h <- glmQLFit(x_h_stage, design_h, robust=TRUE)
head(fit_h$coefficients)
plotQLDisp(fit_h)
summary(fit_h$df.prior)

#between m and F
FM_h<- makeContrasts(F - M, levels=design_h)
res_FM_h <- glmQLFTest(fit_h, contrast=FM_h)
is.de_FM_h <- decideTestsDGE(res_FM_h, p.value=0.05, adjust.method="BH")
#Nukmber of sex-biased genes 
summary(is.de_FM_h)

#get a table with FDR and FC for each gene
hatch_top_fm<- topTags(res_FM_h, n=14000)
topTags(res_FM_h)
hatch_top_fm<-hatch_top_fm$table


plotMD(res_FM_h, 
       status=is.de_FM_h,values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15), 
       main = "DGE between females and males at hatchling stage", 
       hl.cex=0.7 )

##### classification and visualization of sex-biased genes ######
###calculate the average logCPM for males and females separtely, (logCPM to plot)
##logCPM
logcounts_h <- cpm(x_h_stage,log=TRUE)
colnames(logcounts_h)<- c("Tcm_F_H_rep1","Tcm_M_H_rep3",
                          "Tcm_F_H_rep5","Tcm_F_H_rep7", 
                          "Tcm_M_H_rep4","Tcm_M_H_rep6", 
                          "Tcm_M_H_rep2","Tcm_F_H_rep9",
                          "Tcm_M_H_rep8")
rownames(logcounts_h)<- x_h_stage$genes$genes
logcounts_h<- as.table(logcounts_h)
logcounts_h

names<-as.character(hatch_top_fm$genes)
names
logcounts_h<-cbind(rownames(logcounts_h), logcounts_h)
colnames(logcounts_h)[1]<-"Gene_name"
logcounts_h<-as.data.frame(logcounts_h)
logcounts_h_sub<-logcounts_h[logcounts_h$Gene_name %in% names,]
logcounts_h_sub[,-1]<-as.numeric(as.character(unlist(logcounts_h_sub[,-1])))
average_M<-rowMeans(logcounts_h_sub[,c(3,6,7,8,10)])
average_F<-rowMeans(logcounts_h_sub[,c(2,4,5,9)])
logcounts_h_sub<-cbind(average_F, logcounts_h_sub)
logcounts_h_sub<-cbind(average_M, logcounts_h_sub)

hatch_top_fm<-hatch_top_fm[order(hatch_top_fm$genes),]
logcounts_h_sub<- logcounts_h_sub[order(logcounts_h_sub$Gene_name),]
hatch_top_fm<- cbind(logcounts_h_sub$average_M, hatch_top_fm)
hatch_top_fm<- cbind(logcounts_h_sub$average_F, hatch_top_fm)
colnames(hatch_top_fm)[1]<- "Average_logCPM_F"
colnames(hatch_top_fm)[2]<- "Average_logCPM_M"

summary(hatch_top_fm$Average_logCPM_F)
summary(hatch_top_fm$Average_logCPM_M)

## calculating average CPM for males and females and binding it to TOP_fm
counts_h <- cpm(x_h_stage, log = FALSE)
summary(counts_h)
rownames(counts_h)<- x_h_stage$genes$genes
colnames(counts_h)<- sex_h

names_h<-as.character(hatch_top_fm$genes)
counts_h<-cbind(rownames(counts_h), counts_h)
colnames(counts_h)[1]<-"Gene_name"
counts_h<-as.data.frame(counts_h)
counts_h_sub<-counts_h[counts_h$Gene_name %in% names_h,]
counts_h_sub[,-1]<-as.numeric(as.character(unlist(counts_h_sub[,-1])))
average_cpm_M<-rowMeans(counts_h_sub[,c(3,6,7,8,10)])
average_cpm_F<-rowMeans(counts_h_sub[,c(2,4,5,9)])
counts_h_sub<-cbind(average_cpm_F, counts_h_sub)
counts_h_sub<-cbind(average_cpm_M, counts_h_sub)

hatch_top_fm<-hatch_top_fm[order(hatch_top_fm$genes),]
counts_h_sub<- counts_h_sub[order(counts_h_sub$Gene_name),]
hatch_top_fm<- cbind(counts_h_sub$average_cpm_M, hatch_top_fm)
hatch_top_fm<- cbind(counts_h_sub$average_cpm_F, hatch_top_fm)
colnames(hatch_top_fm)[1]<- "Average_CPM_F"
colnames(hatch_top_fm)[2]<- "Average_CPM_M"

#####plot different gene categories####
## => CPM for sex-bias classification
#gene classification#
hatch_top_fm$fact<-rep("NA", nrow(hatch_top_fm))

hatch_top_fm$fact[hatch_top_fm$Average_CPM_M ==(0) & hatch_top_fm$logFC>2] <- "F_limited"
hatch_top_fm$fact[hatch_top_fm$Average_CPM_F ==(0) & hatch_top_fm$logFC< (-2)] <- "M_limited"
hatch_top_fm$fact[hatch_top_fm$logFC>2 & hatch_top_fm$Average_CPM_M > (0)] <- "F_biased"
hatch_top_fm$fact[hatch_top_fm$logFC<(-2) & hatch_top_fm$Average_CPM_F > (0)]<- "M_biased"
hatch_top_fm$fact[hatch_top_fm$logFC>0 & hatch_top_fm$logFC<2]<- "Any_f_biased"
hatch_top_fm$fact[hatch_top_fm$logFC<0 & hatch_top_fm$logFC>-2]<- "Any_m_biased"
hatch_top_fm$fact[hatch_top_fm$FDR >0.05]<- "Not DE"

####### plot
palette(c("palevioletred1","firebrick1","lightslateblue","gainsboro"))

hatch_top_fm$fact<-as.factor(hatch_top_fm$fact)
unbias<-subset(hatch_top_fm, fact == "Not DE")
bias<-subset(hatch_top_fm, fact != "Not DE")
plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", 
     pch=19, 
     ylim = c(-3,17), 
     xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="GE in males and females at the hatchling stage", 
     cex=0.7)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.7)

levels(hatch_top_fm$fact)
table(hatch_top_fm$fact)

legend("topleft", 
       c("Slightly FB", "Female biased","Male biased", "Not DE"), 
       fill=c(1:4), 
       cex = 0.6,
       bty = "n")

