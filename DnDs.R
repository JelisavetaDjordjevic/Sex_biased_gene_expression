###packages
library(readr)
library(edgeR)
library(ggplot2)
library(ppcor)
library(RColorBrewer)
library(ggpubr)

getwd()

##### ********Import TABLE WITH CPM calculated and dnds ******** #####
### this table contains dnds and mean(CPM) values for males and females and overal
all_stages<-read.delim("all_stages_dnds.txt", 
                       header = TRUE, 
                       sep = "\t", 
                       quote = "\"", 
                       fill = TRUE, 
                       comment.char = "", 
                       stringsAsFactors = FALSE)

all_stages$stage <- factor(all_stages$stage,
                           levels = c('hatch','juv','adu'),ordered = TRUE)

### Filtering genes with very high dNdS, setting up the treshhold to 10
sub<-subset(all_stages, dnds<10)
summary(sub$dnds)
hist(sub$dnds, breaks = 300)

#boxplot- dnds at 3 stages 
#First classify genes in broader categories, female-biased/male-biased/un-biased)

sub$fact[sub$fact=="F_limited"]<- "F_biased"
sub$fact[sub$fact=="M_limited"]<- "M_biased"
sub$fact[sub$fact=="M_biased"]<- "M_biased"
sub$fact[sub$fact=="F_biased"]<- "F_biased"
sub$fact[sub$fact=="Any_f_biased"]<- "F_biased"
sub$fact[sub$fact=="Any_m_biased"]<- "M_biased"
levels(as.factor(sub$fact))

p <- ggplot(sub, aes(x=stage, y=dnds, fill=fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("red","deepskyblue1","gray"))
p

##test if DnDs differ between sex-bias categories
#group by "stage"
#wilcoxon test

test<-compare_means(dnds ~ fact, 
                    group.by="stage", 
                    data = sub, 
                    p.adjust.method = "BH")

###### Test correlation beetwen dnds and logfc
#### subset per stage and sex-bias 

### *adult male biased*
###!!!! Split the genes <0 logFC (like this we subset as well un-biased genes!)
sub_adu_m<-subset(all_stages, dnds<10 & stage=="adu")

###### Exclude Un-biased genes for the statistical test
sub_adu_m<- sub_adu_m[sub_adu_m$fact=="Any_m_biased"|
                        sub_adu_m$fact=="M_biased"| 
                        sub_adu_m$fact=="M_limited",]
sub_adu_m$abslogFC<- abs(sub_adu_m$logFC)
levels(as.factor(sub_adu_m$fact))

#plot male biased genes- adult stage
color<-palette(c("slategray1","lightslateblue","navy","gainsboro"))
sub_adu_m$fact<-as.factor(sub_adu_m$fact)
unbias<-subset(sub_adu_m, fact == "Not DE")
bias<-subset(sub_adu_m, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gainsboro", 
     pch=19, 
     ylim = c(0,5), 
     xlim = c(0,15),
     xlab ="|LogFC|", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)
levels(bias$fact)

#pairwise partial correlation test for each pair of variables fc and dnds and cpm
#subset data for columns for test (absfc, average exp, average expr in males, dnds)
sub_ma<- sub_adu_m[,c(6,7,9,10)]
pcor(sub_ma, method = c("spearman"))

### *adult female biased*
sub_adu_f<-subset(all_stages, dnds<10 & stage=="adu")
#without un-biased 
sub_adu_f<- sub_adu_f[sub_adu_f$fact=="Any_f_biased"|
                        sub_adu_f$fact=="F_biased"| 
                        sub_adu_f$fact=="F_limited",]
sub_adu_f$abslogFC<- abs(sub_adu_f$logFC)

levels(as.factor(sub_adu_f$fact))

#plot
color<-palette(c("palevioletred1","firebrick1","gainsboro"))
sub_adu_f$fact<-as.factor(sub_adu_f$fact)
unbias<-subset(sub_adu_f, fact == "Not DE")
bias<-subset(sub_adu_f, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gainsboro", 
     pch=19, 
     ylim = c(0,5), 
     xlim = c(0,15),
     xlab ="LogFC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)

#subset data for columns for test (absfc, average exp, average expr in males, dnds)
sub_fa<- sub_adu_f[,c(6,7,8,10)]
pcor(sub_fa, method = c("spearman"))

## *Juvenile male and female biased*
###female biased 
sub_juv_f<-subset(all_stages, dnds<10 & stage=="juv")
#without NotDE
sub_juv_f<- sub_juv_f[sub_juv_f$fact=="Any_f_biased"|
                        sub_juv_f$fact=="F_biased"| 
                        sub_juv_f$fact=="F_limited",]
sub_juv_f$abslogFC<- abs(sub_juv_f$logFC)

##male biased
sub_juv_m<-subset(all_stages, dnds<10 & stage=="juv")

#without Not DE 
sub_juv_m<- sub_juv_m[sub_juv_m$fact=="Any_m_biased"|
                        sub_juv_m$fact=="M_biased"| 
                        sub_juv_m$fact=="M_limited",]
sub_juv_m$abslogFC<- abs(sub_juv_m$logFC)

levels(as.factor(sub_juv_f$fact))
levels(as.factor(sub_juv_m$fact))

#plot - female biased- juvenile stage
color<-palette(c("palevioletred1","firebrick1","gainsboro"))
sub_juv_f$fact<-as.factor(sub_juv_f$fact)
unbias<-subset(sub_juv_f, fact == "Not DE")
bias<-subset(sub_juv_f, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gainsboro", 
     pch=19, 
     ylim = c(0,5), 
     xlim = c(0,15),
     xlab ="LogFC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)


#plot male-biased - Juvenile stage
color<-palette(c("slategray1","lightslateblue","navy","gainsboro"))
sub_juv_m$fact<-as.factor(sub_juv_m$fact)
unbias<-subset(sub_juv_m, fact == "Not DE")
bias<-subset(sub_juv_m, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gainsboro", 
     pch=19, 
     ylim = c(0,5), 
     xlim = c(0,15),
     xlab ="LogFC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)


##subset data for columns for test (absfc, average exp, average expr in males, dnds)
#male juvenile
sub_mj<- sub_juv_m[,c(6,7,9,10)]
pcor(sub_mj, method = c("spearman"))
#female juvenile
sub_fj<- sub_juv_f[,c(6,7,8,10)]
pcor(sub_fj, method = c("spearman"))


#### *Hatchling female biased* 
sub_hatch_f<-subset(all_stages, dnds<10 & stage=="hatch")
#without notDE
sub_hatch_f<- sub_hatch_f[sub_hatch_f$fact=="Any_f_biased"|sub_hatch_f$fact=="F_biased"| sub_hatch_f$fact=="F_limited",]
sub_hatch_f$abslogFC<- abs(sub_hatch_f$logFC)
levels(as.factor(sub_hatch_f$fact))

#plot
color<-palette(c("palevioletred1","firebrick1","gainsboro"))
sub_hatch_f$fact<-as.factor(sub_hatch_f$fact)
unbias<-subset(sub_hatch_f, fact == "Not DE")
bias<-subset(sub_hatch_f, fact != "Not DE")
plot(unbias$dnds~unbias$abslogFC, 
     col="gainsboro", 
     pch=19, 
     ylim = c(0,5), 
     xlim = c(0,15),
     xlab ="LogFC", 
     ylab = "Dn/DS", 
     cex=0.8)
points(bias$dnds~bias$abslogFC, 
       col=as.factor(bias$fact), 
       pch=19, 
       cex=0.9)
#statistical test
sub_fh<- sub_hatch_f[,c(6,7,8,10)]
pcor(sub_fh, method = c("spearman"))
