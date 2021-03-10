#Adult sex biased gene expression
###packages
library(readr)
library(edgeR)
library(RColorBrewer)
library(ggplot2)
getwd()

#with tcm genome v8
tcm_genome <- read.csv("https://github.com/JelisavetaDjordjevic/Sex_biased_gene_expression/main/timema_counts.csv")
colnames(tcm_genome)

sex<- c("F","M","M","F","M","F","M","M","M","M","M",
        "F","F","F","F","F","M","M","M","M","M","F","F","F")
stage<- c("H","H","A","H","A","H","H","H","A","A","H","A","H",
          "A","A","A","H","J","J","J","J","J","J","J")
sex_stage<-paste(sex, stage, sep="_")
colnames(tcm_genome)

### ADULTS ###
colnames(tcm_genome)
sex_a_sp<-c("M","M","M","M","F","F","F","F") 
#get edgeR class object
x_a_stage<-DGEList(counts = tcm_genome[, c(4,6,10,11,13,15,16,17)],
                   genes = tcm_genome[, 1], group = sex_a_sp)
x_a_stage$samples

#filtering lowly expressed genes
cpm(10, mean(x_a_stage$samples$lib.size))
keep_a <- rowSums(cpm(x_a_stage[,c(1,2,3,4)])>0.5) >=3 | 
        rowSums(cpm(x_a_stage[,c(5,6,7,8)])>0.5) >=3  
x_a_stage<- x_a_stage[keep_a, , keep.lib.sizes=FALSE]
#filtered and kept genes
table(keep_a)
##normalization for the library sizes
x_a_stage<- calcNormFactors(x_a_stage)
x_a_stage$samples

## MDS plot
fac_a<- factor(sex_a_sp)
colours<- brewer.pal(nlevels(fac_a), "Paired")
plotMDS(x_a_stage, labels= factor(x_a_stage$samples$group),
        main="MDS adult stage", col= colours[as.numeric(fac_a)])

#design a model 
design_adu <- model.matrix(~ 0 + sex_a_sp, data= x_a_stage$samples)
colnames(design_adu)<- levels(x_a_stage$samples$group)
design_adu

#Dispersion
x_a_stage <- estimateDisp(x_a_stage, design_adu, robust=TRUE)
summary(x_a_stage$trended.dispersion)
sqrt(x_a_stage$common.dispersion) #this is BCV value
plotBCV(x_a_stage)

#Fitting a GLM model 
fit_adu <- glmQLFit(x_a_stage, design_adu, robust=TRUE)
head(fit_adu$coefficients)
plotQLDisp(fit_adu)
summary(fit_adu$df.prior)

#between m and F
FM_a<- makeContrasts(F- M, levels=design_adu)
res_FM_a <- glmQLFTest(fit_adu, contrast=FM_a)
is.de_FM_a <- decideTestsDGE(res_FM_a, p.value=0.05)
#number of sex-biased genes
summary(is.de_FM_a)
topTags(res_FM_a)

#Get a table with all the genes their FC and FDR, change n for the number of genes  
adu_top_fm<- topTags(res_FM_a, n= 13000)
adu_top_fm<- adu_top_fm$table
#MD plot
plotMD(res_FM_a, status=is.de_FM_a,values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15), 
       main = "DE genes between females and males at adult stage")
abline(h=c(-2,2), col="black")

###calculate the average logCPM and CPM values for males and females separtely
##to classify genes based on their FC and average expression
##logCPM (to plot)
logcounts_a <- cpm(x_a_stage,log=TRUE)

colnames(logcounts_a)<- sex_a_sp
rownames(logcounts_a)<- x_a_stage$genes$genes

names_adu<-as.character(adu_top_fm$genes)
logcounts_a<-cbind(rownames(logcounts_a), logcounts_a)
colnames(logcounts_a)[1]<-"Gene_name"
logcounts_a<-as.data.frame(logcounts_a)
logcounts_a_sub<-logcounts_a[logcounts_a$Gene_name %in% names_adu,]
logcounts_a_sub[,-1]<-as.numeric(as.character(unlist(logcounts_a_sub[,-1])))
average_M<-rowMeans(logcounts_a_sub[,c(2,3,4,5)])
average_F<-rowMeans(logcounts_a_sub[,c(6,7,8,9)])
logcounts_a_sub<-cbind(average_F, logcounts_a_sub)
logcounts_a_sub<-cbind(average_M, logcounts_a_sub)

adu_top_fm<-adu_top_fm[order(adu_top_fm$genes),]
logcounts_a_sub<- logcounts_a_sub[order(logcounts_a_sub$Gene_name),]
adu_top_fm<- cbind(logcounts_a_sub$average_M, adu_top_fm)
adu_top_fm<- cbind(logcounts_a_sub$average_F, adu_top_fm)
colnames(adu_top_fm)[1]<- "Average_logCPM_F"
colnames(adu_top_fm)[2]<- "Average_logCPM_M"

summary(adu_top_fm$Average_logCPM_F)
summary(adu_top_fm$Average_logCPM_M)

## calculating average CPM for males and females and binding it to TOP_fm 
# CPM => for sex-bias classification
counts_a <- cpm(x_a_stage,log=FALSE)
colnames(counts_a)<- sex_a_sp
rownames(counts_a)<- x_a_stage$genes$genes
names_adu<-as.character(adu_top_fm$genes)
counts_a<-cbind(rownames(counts_a), counts_a)
colnames(counts_a)[1]<-"Gene_name"
counts_a<-as.data.frame(counts_a)
counts_a_sub<-counts_a[counts_a$Gene_name %in% names_adu,]
counts_a_sub[,-1]<-as.numeric(as.character(unlist(counts_a_sub[,-1])))
average_cpm_M<-rowMeans(counts_a_sub[,c(2,3,4,5)])
average_cpm_F<-rowMeans(counts_a_sub[,c(6,7,8,9)])
counts_a_sub<-cbind(average_cpm_F, counts_a_sub)
counts_a_sub<-cbind(average_cpm_M, counts_a_sub)

adu_top_fm<-adu_top_fm[order(adu_top_fm$genes),]
counts_a_sub<- counts_a_sub[order(counts_a_sub$Gene_name),]
adu_top_fm<- cbind(counts_a_sub$average_cpm_M, adu_top_fm)
adu_top_fm<- cbind(counts_a_sub$average_cpm_F, adu_top_fm)
colnames(adu_top_fm)[1]<- "Average_CPM_F"
colnames(adu_top_fm)[2]<- "Average_CPM_M"

#gene classification#
adu_top_fm$fact<-rep("NA", nrow(adu_top_fm))
adu_top_fm$fact[adu_top_fm$Average_CPM_M ==(0) & adu_top_fm$logFC>2] <- "F_limited"
adu_top_fm$fact[adu_top_fm$Average_CPM_F ==(0) & adu_top_fm$logFC< (-2)] <- "M_limited"
adu_top_fm$fact[adu_top_fm$logFC>2 & adu_top_fm$Average_CPM_M > (0)] <- "F_biased"
adu_top_fm$fact[adu_top_fm$logFC<(-2) & adu_top_fm$Average_CPM_F > (0)]<- "M_biased"
adu_top_fm$fact[adu_top_fm$logFC>0 & adu_top_fm$logFC<2]<- "Any_f_biased"
adu_top_fm$fact[adu_top_fm$logFC<0 & adu_top_fm$logFC>-2]<- "Any_m_biased"
adu_top_fm$fact[adu_top_fm$FDR >0.05]<- "Not DE"

####### plot
color<-palette(c("palevioletred1","slategray1","firebrick1",
                 "red4","lightslateblue","navy","gainsboro"))
adu_top_fm$fact<-as.factor(adu_top_fm$fact)
unbias<-subset(adu_top_fm, fact == "Not DE")
bias<-subset(adu_top_fm, fact != "Not DE")
plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", 
     pch=19, 
     ylim = c(-3,17), 
     xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="GE in males and females at the adult stage", 
     cex=0.8)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)

legend("topleft", 
       c("Slightly FB","Slightly MB", "Strongly FB", 
         "Female limited","Strongly MB","Male limited", "Not DE"), 
       fill=c(1:7), 
       cex = 0.6,
       bty = "n")


getwd()
print(sessionInfo())
