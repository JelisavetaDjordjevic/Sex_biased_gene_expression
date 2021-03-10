#Juvenile sex biased gene expression
###packages
library(readr)
library(edgeR)
library(RColorBrewer)
getwd()

#with tcm genome v8
tcm_genome <- read.csv("https://github.com/JelisavetaDjordjevic/Sex_biased_gene_expression/tree/main/Data/timema_counts.csv")
colnames(tcm_genome)

sex<- c("F","M","M","F","M","F","M","M","M","M","M",
        "F","F","F","F","F","M","M","M","M","M","F","F","F")
stage<- c("H","H","A","H","A","H","H","H","A","A","H","A","H",
          "A","A","A","H","J","J","J","J","J","J","J")
sex_stage<-paste(sex, stage, sep="_")
colnames(tcm_genome)

###JUVENILE###
colnames(tcm_genome)
sex_j<-c("M", "M", "M", "M", "F", "F","F")

#sex_j<-c("Tcm_M_J_rep2", "Tcm_F_J_rep2" ,"Tcm_M_J_rep3", "Tcm_M_J_rep4", "Tcm_F_J_rep3", "Tcm_F_J_rep1","Tcm_M_J_rep1") 

#get EdgeR object 
x_j_stage<-DGEList(counts = tcm_genome[, c(19,20,21,22,23,24,25)],
                   genes = tcm_genome[, 1], 
                   group = sex_j)
colnames(x_j_stage)

#filtering genes with low expression
cpm(10, mean(x_j_stage$samples$lib.size))
keep_j <- rowSums(cpm(x_j_stage[,c(1,2,3,4)])>0.5) >=3 | 
  rowSums(cpm(x_j_stage[,c(5,6,7)])>0.5) >=2  

x_j_stage<- x_j_stage[keep_j, , keep.lib.sizes=FALSE]
#number of genes filtered out and tha ones that passed the filtering
table(keep_j)

#normalization for library sizes
x_j_stage<- calcNormFactors(x_j_stage)
x_j_stage$samples

#MDS plot 
fac_j<- factor(sex_j)
colours<- brewer.pal(nlevels(fac_j), "Paired")
plotMDS(x_j_stage, 
        labels= factor(x_j_stage$samples$group), 
        col= colours[as.numeric(fac_j)])

#design a model
design_j <- model.matrix(~ 0 + group, data= x_j_stage$samples)
colnames(design_j)<- levels(x_j_stage$samples$group)
design_j
#disperzion estimation
x_j_stage <- estimateDisp(x_j_stage, design_j, robust=TRUE)
summary(x_j_stage$trended.dispersion)
x_j_stage$common.dispersion
plotBCV(x_j_stage)

#fitting model 
fit_j <- glmQLFit(x_j_stage, design_j, robust=TRUE)
head(fit_j$coefficients)
plotQLDisp(fit_j)
summary(fit_j$df.prior)

#between m and F
FM_j<- makeContrasts(F - M, levels=design_j)
res_FM_j <- glmQLFTest(fit_j, contrast=FM_j)

#number of sex-biased genes 
is.de_FM_j <- decideTestsDGE(res_FM_j, p.value=0.05)
summary(is.de_FM_j)

topTags(res_FM_j, n= 13000)
juv_top_fm<- topTags(res_FM_j, n= 13000)
juv_top_fm<- juv_top_fm$table

plotMD(res_FM_j, 
       status=is.de_FM_j,
       values=c(1,-1),
       col=c("red","blue"),
       legend="topright", 
       ylim=c(-15,15), 
       main = "DE genes between sexual females and males at juvenile stage", 
       cex=0.7 )
abline(h=c(-2,2), col="blue")

##### classification and visualization of sex-biased genes ######
###calculate the average logCPM for males and females separtely, for DE genes => logCPM to plot
x_j_stage$samples
##log counts for juv stage 
logcounts_j <- cpm(x_j_stage,log=TRUE)
colnames(logcounts_j)<-c("Tcm_M_J_rep4","Tcm_M_J_rep2",
                         "Tcm_M_J_rep3","Tcm_M_J_rep1" ,
                         "Tcm_F_J_rep1","Tcm_F_J_rep3",
                         "Tcm_F_J_rep2")
rownames(logcounts_j)<- x_j_stage$genes$genes

names_juv<-as.character(juv_top_fm$genes)
logcounts_j<-cbind(rownames(logcounts_j), logcounts_j)
colnames(logcounts_j)[1]<-"Gene_name"
logcounts_j<-as.data.frame(logcounts_j)
logcounts_j_sub<-logcounts_j[logcounts_j$Gene_name %in% names_juv,]
logcounts_j_sub[,-1]<-as.numeric(as.character(unlist(logcounts_j_sub[,-1])))
average_M<-rowMeans(logcounts_j_sub[,c(2,3,4,5)])
average_F<-rowMeans(logcounts_j_sub[,c(6,7,8)])
logcounts_j_sub<-cbind(average_F, logcounts_j_sub)
logcounts_j_sub<-cbind(average_M, logcounts_j_sub)

juv_top_fm<-juv_top_fm[order(juv_top_fm$genes),]
logcounts_j_sub<- logcounts_j_sub[order(logcounts_j_sub$Gene_name),]
juv_top_fm<- cbind(logcounts_j_sub$average_M, juv_top_fm)
juv_top_fm<- cbind(logcounts_j_sub$average_F, juv_top_fm)
colnames(juv_top_fm)[1]<- "Average_logCPM_F"
colnames(juv_top_fm)[2]<- "Average_logCPM_M"

### calculat average CPM for males and females and bind it to TOP_fm 
## => CPM for sex-bias classification
counts_j <- cpm(x_j_stage, log = FALSE)
summary(counts_j)
rownames(counts_j)<- x_j_stage$genes$genes
colnames(counts_j)<- sex_j

names_juv<-as.character(juv_top_fm$genes)
counts_j<-cbind(rownames(counts_j), counts_j)
colnames(counts_j)[1]<-"Gene_name"
counts_j<-as.data.frame(counts_j)
counts_j_sub<-counts_j[counts_j$Gene_name %in% names_juv,]
counts_j_sub[,-1]<-as.numeric(as.character(unlist(counts_j_sub[,-1])))
average_cpm_M<-rowMeans(counts_j_sub[,c(2,3,4,5)])
average_cpm_F<-rowMeans(counts_j_sub[,c(6,7,8)])
counts_j_sub<-cbind(average_cpm_F, counts_j_sub)
counts_j_sub<-cbind(average_cpm_M, counts_j_sub)

juv_top_fm<-juv_top_fm[order(juv_top_fm$genes),]
counts_j_sub<- counts_j_sub[order(counts_j_sub$Gene_name),]
juv_top_fm<- cbind(counts_j_sub$average_cpm_M, juv_top_fm)
juv_top_fm<- cbind(counts_j_sub$average_cpm_F, juv_top_fm)
colnames(juv_top_fm)[1]<- "Average_CPM_F"
colnames(juv_top_fm)[2]<- "Average_CPM_M"

#gene classification#
juv_top_fm$fact<-rep("NA", nrow(juv_top_fm))
juv_top_fm$fact[juv_top_fm$Average_CPM_M ==(0) & juv_top_fm$logFC>2] <- "F_limited"
juv_top_fm$fact[juv_top_fm$Average_CPM_F ==(0) & juv_top_fm$logFC< (-2)] <- "M_limited"
juv_top_fm$fact[juv_top_fm$logFC>2 & juv_top_fm$Average_CPM_M >(0)] <- "F_biased"
juv_top_fm$fact[juv_top_fm$logFC<(-2) & juv_top_fm$Average_CPM_F >(0)]<- "M_biased"
juv_top_fm$fact[juv_top_fm$logFC>0 & juv_top_fm$logFC<2]<- "Any_f_biased"
juv_top_fm$fact[juv_top_fm$logFC<0 & juv_top_fm$logFC>-2]<- "Any_m_biased"
juv_top_fm$fact[juv_top_fm$FDR >0.05]<- "Not DE"
table(juv_top_fm$fact)

####### plot
palette(c("palevioletred1","slategray1","firebrick1",
          "red4","lightslateblue", "navy","gainsboro"))

juv_top_fm$fact<-as.factor(juv_top_fm$fact)
unbias<-subset(juv_top_fm, fact == "Not DE")
bias<-subset(juv_top_fm, fact != "Not DE")

plot(unbias$Average_logCPM_M~unbias$Average_logCPM_F, 
     col="gainsboro", 
     pch=19, 
     ylim = c(-3,17), 
     xlim = c(-3,17),
     xlab = "Average female expression [logCPM]", 
     ylab = "Average male expression [logCPM]", 
     main="GE in males and females at the juvenile stage", 
     cex=0.7)

points(bias$Average_logCPM_M~bias$Average_logCPM_F, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.7)

legend("topleft", 
       c("Slightly FB","Slightly MB","Strongly FB",
         "Female limited","Strongly MB","Male limited", "Not DE"), 
       fill=c(1:7), 
       cex = 0.6,
       bty = "n")

