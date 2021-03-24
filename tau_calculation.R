###packages
library(readr)
library(edgeR)
library(ggplot2)
library(ggpubr)
getwd()

#### getting table to calculate Tau

#import table with stage, bias, FDR...
#import table with counts
all_stages<-read.delim("all_stages", 
                       header = TRUE, 
                       sep = "\t", 
                       quote = "\"",
                       dec = ".", 
                       fill = TRUE, 
                       comment.char = "", 
                       stringsAsFactors = FALSE)
#counts
tcm_genome <- read.csv("/Users/jdjordje/Desktop/Jeli/job/genome_analysis/R_scripts/timema_counts.csv")
levels(as.factor(all_stages$fact))

#### get uniq list of genes for all 3 stages
uniq_list_genes<- unique(all_stages$genes)
uniq_list_genes<-factor(uniq_list_genes)
length(uniq_list_genes)
length(all_stages$genes)

#Timema counts
colnames(tcm_genome)

sex<- c("F","M","M","F","M","F","M","M","M","M","M",
        "F","F","F","F","F","M","M","M","M","M","F","F","F")
stage<- c("H","H","A","H","A","H","H","H","A","A","H","A","H",
          "A","A","A","H","J","J","J","J","J","J","J")
sex_stage<-paste(sex, stage, sep="_")

#get Digital Gene Expression data -class 
x<-DGEList(counts = tcm_genome[, 2:25],
           genes = tcm_genome[, 1], 
           group = sex_stage)
x$samples
x$samples$group

################ *calculating Tau* ##################

##for calculating tau dont log CPM (because of negative values)
##get CPM values
counts_all <- cpm(x,log=FALSE)
colnames(counts_all)<- sex_stage
counts_all<-as.data.frame(counts_all)
fac<-factor(as.character(x$genes$genes))
counts_all$gene<-fac

##get counts for genes that are expressed in all three stages
## by matching the list of unique genes and the counts
sub_counts<-counts_all[counts_all$gene %in% uniq_list_genes,]

#Calculate median expression in males and females per stage
#Males adults
median_M_adu<-sub_counts[,c(3,5,9,10)]
median_M_adu<-apply(median_M_adu, 1, median)
#Females adult
median_F_adu<-sub_counts[,c(12,14,15,16)]
median_F_adu<-apply(median_F_adu, 1, median)
#Males juveniles
median_M_juv<-sub_counts[,c(18,19,20,21)]
median_M_juv<-apply(median_M_juv, 1, median)
#Females juveniles
median_F_juv<-sub_counts[,c(22,23,24)]
median_F_juv<-apply(median_F_juv, 1, median)
#Males hatchling
median_M_hat<-sub_counts[,c(2,7,8,11,17)]
median_M_hat<-apply(median_M_hat, 1, median)
#females hatchling
median_F_hat<-sub_counts[,c(1,4,6,13)]
median_F_hat<-apply(median_F_hat, 1, median)

### Bind all the median expressions together
sub_counts<-cbind(median_F_adu, sub_counts)
sub_counts<-cbind(median_M_adu, sub_counts)
sub_counts<-cbind(median_F_juv, sub_counts)
sub_counts<-cbind(median_M_juv, sub_counts)
sub_counts<-cbind(median_F_hat, sub_counts)
sub_counts<-cbind(median_M_hat, sub_counts)

#Keep only median expressions and gene names
exp<-sub_counts[, c(1:6,31)]

#Expression females, matrix
exp_f<-sub_counts[, c(2,4,6,31)]
rownames(exp_f)<- exp_f$gene
exp_f_mat<- exp_f[, c(1:3)]
exp_f_mat<- as.matrix(exp_f_mat)

#Expression males, matrix
exp_m<-sub_counts[, c(1,3,5,31)]
rownames(exp_m)<- exp_m$gene
exp_m_mat<- exp_m[, c(1:3)]
exp_m_mat<- as.matrix(exp_m_mat)

#Run Tau function 
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    } 
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  } 
  return(res)
}

#Stage specificity of gene expression in females
Tau_values_f<-apply(exp_f_mat, 1, fTau)
summary(Tau_values_f)

#stage specificity of gene expression in males
Tau_values_m<-apply(exp_m_mat, 1, fTau)
summary(Tau_values_m)

#Get tau values for males and females 
Tau_values<- cbind(Tau_values_f,Tau_values_m)
Tau_values<-as.data.frame(Tau_values)
Tau_values<- cbind(sub_counts$gene,Tau_values)

#create a empty column for males and females in all_stages data frame
all_stages$tau_f <- rep("NA", nrow(all_stages))
all_stages$tau_m <- rep("NA", nrow(all_stages))

#match the two data frames tau_values and all_stages
colnames(Tau_values)[1]<-"gene"
all_stages$tau_f <- Tau_values$Tau_values_f[match(all_stages$genes,Tau_values$gene)]
all_stages$tau_m <- Tau_values$Tau_values_m[match(all_stages$genes,Tau_values$gene)]

########## Tau for general sex bias categories: Female, male and UB

levels(as.factor(all_stages$fact))
all_stages$fact[all_stages$fact=="Any_f_biased"]<- "F_biased"
all_stages$fact[all_stages$fact=="Any_m_biased"]<- "M_biased"
all_stages$fact[all_stages$fact=="F_limited"]<- "F_biased"
all_stages$fact[all_stages$fact=="M_limited"]<- "M_biased"

levels(as.factor(all_stages$fact))
all_stages$stage <- factor(all_stages$stage,
                           levels = c('hatch','juv','adu'),ordered = TRUE)

p <- ggplot(all_stages, 
            aes(x=stage, 
                y=tau_m, 
                fill=fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("red","deepskyblue1","gray"))
p

p <- ggplot(all_stages, 
            aes(x=stage, 
                y=tau_f, 
                fill=fact)) + 
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("red","deepskyblue1","gray"))
p

### ***  statistics  ***  ####

###wilcoxon test###
hist(all_stages$tau_m, breaks = 100)
hist(all_stages$tau_f, breaks = 100)
shapiro.test(all_stages$tau_m[0:5000])

test_tauf<-compare_means(tau_f ~ fact, 
                         group.by="stage", 
                         data = all_stages, 
                         p.adjust.method = "BH")

test_taum<-compare_means(tau_m ~ fact, 
                         group.by="stage", 
                         data = all_stages, 
                         p.adjust.method = "BH")
test<- rbind(test_taum, test_tauf)
