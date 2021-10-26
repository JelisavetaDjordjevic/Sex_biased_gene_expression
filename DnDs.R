library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ppcor)
library(RColorBrewer)
library(gplots)


##### ********Import TABLE WITH CPM calculated and dnds ******** #####
### this table contains dnds for Tcm branch,mean(CPM) values for males and females, GC, Log2FC, sex-bias factor, 

all_stages_dnds <- read.delim("https://github.com/JelisavetaDjordjevic/Sex_biased_gene_expression/tree/main/Data/all_stages_dnds.txt")
table(as.factor(all_stages_dnds$fact))

df<- all_stages_dnds
df[is.na(df)] <- "NA"
df<-df[!grepl("NA", df$dnds),]
df$dnds<-as.numeric(as.character(df$dnds))

## change sex-b factor to more general only un/mb/fb
df$fact[df$fact=="Any_m_biased" | df$fact=="M_limited"]<- "M_biased"
df$fact[df$fact=="Any_f_biased" | df$fact=="F_limited"]<- "F_biased"
df$fact[df$fact=="Any_m_biased" | df$fact=="M_limited"]<- "M_biased"
df$fact[df$fact=="Any_f_biased" | df$fact=="F_limited"]<- "F_biased"
df$fact[df$fact=="Any_m_biased" | df$fact=="M_limited"]<- "M_biased"
df$fact[df$fact=="Any_f_biased" | df$fact=="F_limited"]<- "F_biased"

levels(as.factor(df$fact))

#plot  dnds values per sex bias category groupped by developmental stage
#order groups
df$stage <- factor(df$stage,
                   levels = c('hatch','juv','adu'),ordered = TRUE)
#plot
p <- ggplot(df, aes(x=stage, y=dnds, fill=fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("red","deepskyblue1","gray"))
p

#Wilcoxon test to compare dnds values between sex-bias categories
test<-compare_means(dnds ~ fact, group.by="stage", data = df, p.adjust.method = "BH")

####### subset per stage and look for the correlation in strength of the bias (FC) and dnds
##### split male and female biased genes

##adult male biased 

sub_adu_m<-subset(all_stages_dnds, dnds<10 & stage=="adu")

sub_adu_m<- sub_adu_m[sub_adu_m$fact=="Not DE" | sub_adu_m$fact=="Any_m_biased"|sub_adu_m$fact=="M_biased"| sub_adu_m$fact=="M_limited",]
sub_adu_m$abslogFC<- abs(sub_adu_m$logFC)

###### without the Not DE (for the plot)
sub_adu_m_sb<- sub_adu_m[sub_adu_m$fact=="Any_m_biased"|sub_adu_m$fact=="M_biased"| sub_adu_m$fact=="M_limited",]
sub_adu_m_sb$abslogFC<- abs(sub_adu_m_sb$logFC)
levels(as.factor(sub_adu_m_sb$fact))

#is there a correlation between dnds and logfc
hist(sub_adu_m$dnds, breaks = 100)
hist(sub_adu_m$abslogFC, breaks = 100)

#plot male biased genes
color<-palette(c("#b8e0ff","#3877ff","#001d5c","#DEDEDE"))

sub_adu_m$fact<-as.factor(sub_adu_m$fact)
unbias<-subset(sub_adu_m, fact == "Not DE")
bias<-subset(sub_adu_m, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gray95", 
     pch=19, 
     ylim = c(0,1), 
     xlim = c(0,15),
     xlab ="|FC|", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)
levels(bias$fact)


#partial correlation test between fc and dnds, when taking into account GC and average expression CPM
#subset data columns for the test (dnds,absfc, average exp,GC bias)
#Only use subset with sex-biased genes

sub_ma<- sub_adu_m_sb[,c(6,11,7,10)]
pcor(sub_ma, method = c("spearman"))
p=c(0.000000e+00, 1.603566e-04, 9.745023e-01, 4.865218e-20)
p.adjust(p, method ="fdr", n = length(p))


##adult female biased 
sub_adu_f<-subset(all_stages_dnds,dnds<10 & stage=="adu")

sub_adu_f<- sub_adu_f[sub_adu_f$fact=="Not DE" | sub_adu_f$fact=="Any_f_biased"|sub_adu_f$fact=="F_biased"| sub_adu_f$fact=="F_limited",]
sub_adu_f$abslogFC<- abs(sub_adu_f$logFC)

#without NOT DE genes - for the statistical test
sub_adu_f_sb<- sub_adu_f[sub_adu_f$fact=="Any_f_biased"|sub_adu_f$fact=="F_biased"| sub_adu_f$fact=="F_limited",]
sub_adu_f_sb$abslogFC<- abs(sub_adu_f_sb$logFC)
levels(as.factor(sub_adu_f_sb$fact))

#is there a correlation between dnds and logfc
hist(sub_adu_f$dnds, breaks = 100)
hist(sub_adu_f$logFC, breaks = 100)

#plot
#color<-palette(c("palevioletred1","firebrick1","gainsboro"))
color<-palette(c("#ffd1d1","#fa0000","#DEDEDE"))

sub_adu_f$fact<-as.factor(sub_adu_f$fact)
unbias<-subset(sub_adu_f, fact == "Not DE")
bias<-subset(sub_adu_f, fact != "Not DE")
table(bias$fact)
plot(unbias$dnds~unbias$abslogFC, 
     col="gray95", 
     pch=19, 
     ylim = c(0,1), 
     xlim = c(0,15),
     xlab ="FC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)


#subset data for columns for test (absfc, average exp,gc, dnds)
sub_fa<- sub_adu_f_sb[,c(6,7,10,11)]
pcor(sub_fa, method = c("spearman"))
p<-c(0.00E+00,	5.38E-01,	2.17E-08,	6.70E-12)
p.adjust(p, method ="fdr", n = length(p))


##juvenile male and female biased
##female biased 
sub_juv_f<-subset(all_stages_dnds, dnds<10 & stage=="juv")
sub_juv_f<- sub_juv_f[sub_juv_f$fact=="Not DE" | sub_juv_f$fact=="Any_f_biased"|sub_juv_f$fact=="F_biased"| sub_juv_f$fact=="F_limited",]
sub_juv_f$abslogFC<- abs(sub_juv_f$logFC)

#without NotDE
sub_juv_f_sb<- sub_juv_f[sub_juv_f$fact=="Any_f_biased"|sub_juv_f$fact=="F_biased"| sub_juv_f$fact=="F_limited",]
sub_juv_f_sb$abslogFC<- abs(sub_juv_f_sb$logFC)

##male biased
sub_juv_m<-subset(all_stages_dnds, dnds<10 & stage=="juv")
sub_juv_m<- sub_juv_m[sub_juv_m$fact=="Not DE" | sub_juv_m$fact=="Any_m_biased"|sub_juv_m$fact=="M_biased"| sub_juv_m$fact=="M_limited",]
sub_juv_m$abslogFC<- abs(sub_juv_m$logFC)

#without Not DE
sub_juv_m_sb<- sub_juv_m[sub_juv_m$fact=="Any_m_biased"|sub_juv_m$fact=="M_biased"| sub_juv_m$fact=="M_limited",]
sub_juv_m_sb$abslogFC<- abs(sub_juv_m_sb$logFC)

levels(as.factor(sub_juv_f$fact))
levels(as.factor(sub_juv_m$fact))

#is there a correlation between dnds and logfc
hist(sub_juv_f$dnds, breaks = 100)
hist(sub_juv_f$logFC, breaks = 100)
hist(sub_juv_m$dnds, breaks = 100)
hist(sub_juv_m$abslogFC, breaks = 100)

#plot
#color<-palette(c("palevioletred1","firebrick1","gainsboro"))
color<-palette(c("#ffd1d1","#fa0000","#DEDEDE"))
sub_juv_f$fact<-as.factor(sub_juv_f$fact)
unbias<-subset(sub_juv_f, fact == "Not DE")
bias<-subset(sub_juv_f, fact != "Not DE")
table(bias$fact)
plot(unbias$dnds~unbias$abslogFC, 
     col="gray95", 
     pch=19, 
     ylim = c(0,1), 
     xlim = c(0,15),
     xlab ="FC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)


#plot male-b
#color<-palette(c("slategray1","lightslateblue","navy","gainsboro"))
color<-palette(c("#3877ff","#001d5c","#DEDEDE"))
sub_juv_m$fact<-as.factor(sub_juv_m$fact)
unbias<-subset(sub_juv_m, fact == "Not DE")
bias<-subset(sub_juv_m, fact != "Not DE")
table(bias$fact)
plot(unbias$dnds~unbias$abslogFC, 
     col="gray95", 
     pch=19, 
     ylim = c(0,1), 
     xlim = c(0,15),
     xlab ="|FC|", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)


#subset data for columns for test (absfc, average exp,gc, dnds)
#juvenile male
sub_mj<- sub_juv_m_sb[,c(6,7,10,11)]
pcor(sub_mj, method = c("spearman"))
p<-c(0,0.324760949,0.000963596,0.451646396)
p.adjust(p, method ="fdr", n = length(p))

#subset data for columns for test (absfc, average exp,gc, dnds)
##juvenile famale
sub_fj<- sub_juv_f_sb[,c(6,7,10,11)]
pcor(sub_fj, method = c("spearman"))
p<-c(0.0000000, 0.9316842, 0.8495992, 0.2863937)
p.adjust(p, method ="fdr", n = length(p))

##hatchling female biased 
sub_hatch_f<-subset(all_stages_dnds, dnds<10 & stage=="hatch")
sub_hatch_f<- sub_hatch_f[sub_hatch_f$fact=="Not DE" | sub_hatch_f$fact=="Any_f_biased"|sub_hatch_f$fact=="F_biased"| sub_hatch_f$fact=="F_limited",]
sub_hatch_f$abslogFC<- abs(sub_hatch_f$logFC)

#without notDE
sub_hatch_f_sb<- sub_hatch_f[sub_hatch_f$fact=="Any_f_biased"|sub_hatch_f$fact=="F_biased"| sub_hatch_f$fact=="F_limited",]
sub_hatch_f_sb$abslogFC<- abs(sub_hatch_f_sb$logFC)
levels(as.factor(sub_hatch_f$fact))

#plot
color<-palette(c("#fa0000","#DEDEDE"))
sub_hatch_f$fact<-as.factor(sub_hatch_f$fact)
unbias<-subset(sub_hatch_f, fact == "Not DE")
bias<-subset(sub_hatch_f, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gray95", 
     pch=19, 
     ylim = c(0,1), 
     xlim = c(0,15),
     xlab ="FC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)

sub_fh<- sub_hatch_f_sb[,c(6,7,10,11)]

pcor(sub_fh, method = c("spearman"))