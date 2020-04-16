library(ggfortify)
library(survival)
library(survminer)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(matrixStats)
library(gplots)
######################################## preprocessing ##########################################
options(stringsAsFactors = FALSE)

WES<-read.csv('~/Box/CRC-iMM-nonMet/features.summary.w.clin.SK.csv')[,1:25] #Table S2
TSO500<-read.csv('~/Box/CRC-iMM-nonMet/TSO500.csv') #Table S4

WES<-WES[-which(WES$Sample=='CR698'),] # low quality sample

WES$sidedness[!(WES$sidedness=='right' | WES$sidedness=='left')]<-'other'
WES$Adj<-as.factor(WES$Adj)
WES$Age<-as.numeric(WES$Age)

TIL<-read.csv('Box/CRC-iMM-nonMet/TILs/deconv.txt',sep = '\t') #Table S?
TIL<-TIL[TIL$sample%in%WES$Sample,]
TIL[TIL<0]<-0
HERV<-read.delim('~/Box/CRC-iMM-nonMet/HERV/normalized.counts.tsv') #Table S1
hERV <- HERV[,which(colnames(HERV)%in%WES$Sample)]
hERV[hERV < 0.5] <- 0
hERV<-hERV[which(rowMedians(as.matrix(hERV))>0.1),]

TPM<-read.csv('Box/CRC-iMM-nonMet/expression/TPM.tsv',stringsAsFactors = FALSE, sep = '\t') # Table S3
rownames(TPM)<-TPM$ENSMBL
TPM<-TPM[,colnames(TPM)%in%WES$Sample]

combined<-cbind(WES,TIL[,2:4],colMedians(as.matrix(hERV)))
combined$Stage[grep("III",combined$Stage)]<-'III'
combined$Stage[-grep("III",combined$Stage)]<-'II'
colnames(combined)[29]<-'median.hERV'
combined$total.TIL<-as.numeric(rowSums(combined[,26:28]))

GZMB<-as.numeric(TPM["ENSG00000100453",])
PRF<-as.numeric(TPM["ENSG00000180644",])
PD1<-as.numeric(TPM["ENSG00000188389",])
FOXP3<-as.numeric(TPM["ENSG00000049768",])
CTL<-sqrt(PRF*GZMB)
EOMES<-as.numeric(TPM["ENSG00000163508",])
IFNG<-as.numeric(TPM["ENSG00000111537",])


combined$CTL<-CTL
combined$PD1<-PD1
combined$Treg<-FOXP3
combined$EOMES<-EOMES
combined$IFNG<-IFNG

######### feature correlations ###########
# 1. TIL vs. MSI/TMB
boxplot(combined$CD8~combined$MSI)
boxplot(combined$CD4~combined$MSI)
boxplot(combined$CD19~combined$MSI)

plot(combined$CD8,combined$all.TMB)
plot(combined$CD4,combined$all.TMB)
plot(combined$CD19,combined$all.TMB)

a<-combined$CD8[combined$MSI=='MSS']
b<-combined$CD8[combined$MSI=='MSH']

a<-combined$EOMES[combined$MSI=='MSS']
b<-combined$EOMES[combined$MSI=='MSH']


a<-sqrt(combined$EOMES*combined$IFNG)[combined$MSI=='MSS']
b<-sqrt(combined$EOMES*combined$IFNG)[combined$MSI=='MSH']

a<-sqrt(combined$EOMES*combined$IFNG*combined$CD8)[combined$MSI=='MSS']
b<-sqrt(combined$EOMES*combined$IFNG*combined$CD8)[combined$MSI=='MSH']

a<-sqrt(combined$CTL)[combined$MSI=='MSS']
b<-sqrt(combined$CTL)[combined$MSI=='MSH']


boxplot(cbind(a,b), names= c('MSS','MSH'), col = 'red', outline = FALSE, ylim = c(0,2.1))
stripchart(cbind(a,b)~cbind(rep(1,length(a)),rep(2,length(b))),
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

wilcox.test(a,b)


a<-combined$CD8[combined$Matastasis==1]
b<-combined$CD8[combined$Matastasis==0]


# 2. TIL vs hERV
plot(combined$CD8,combined$median.hERV)
cor.test(combined$CD8,combined$median.hERV, method = 'spearman')
plot(combined$CD4,combined$median.hERV)
cor.test(combined$CD4,combined$median.hERV, method = 'spearman')
plot(combined$CD19,combined$median.hERV)
cor.test(combined$CD19,combined$median.hERV, method = 'spearman')

my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)

#M<-cor(combined[,c(11,12,13,26:33)],method = 'spearman')
M<-cor(combined[,c(2,3,5:10,26:29,32:35)],method = 'spearman')

heatmap.2(M, col = my_palette, cexRow = 0.7, cexCol = 0.7)

cor.test(combined$total.TIL,combined$purity,method = 'spearman')

HIF1a.id<-which(rownames(TPM)=='ENSG00000100644')
plot(as.numeric(TPM[HIF1a.id,]),combined$median.hERV)
cor.test(as.numeric(TPM[HIF1a.id,]),combined$median.hERV, method = 'spearman')
cor.test(as.numeric(TPM[HIF1a.id,]),combined$CD8)


M<-t(cor(combined[,c(26:29)], t(hERV),method = 'spearman'))

boxplot(M, names= c('CD8','CD4','CD19','median.hERV'), col = 'red', outline = FALSE, ylim = c(-0.3,1) , ylab = 'Spearman correlation')
stripchart(cbind(M)~cbind(rep(1,length(M[,1])),rep(2,length(M[,1])),rep(3,length(M[,1])),rep(4,length(M[,1]))),
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.1)


##### cutoff #####
clin.fac<-combined
n<-30
num<-c(2:3,5:13,25:28,30:35)
for(i in 1:length(num)){
  idx<-num[i]
  clin.fac[,idx]<-as.factor(ifelse(clin.fac[,idx] > quantile(clin.fac[,idx],prob=1-n/100),'high','low'))
}

n<-30
num<-c(29)
for(i in 1:length(num)){
  idx<-num[i]
  clin.fac[,idx]<-as.factor(ifelse(clin.fac[,idx] > quantile(clin.fac[,idx],prob=1-n/100),'high','low'))
}

##########
fit <- surv_fit(Surv(OS, Dead) ~ clin.fac$median.hERV, data = clin.fac, plotTable=TRUE)
autoplot(fit, risk.table=TRUE)
surv_pvalue(fit)

fit <- surv_fit(Surv(RFS, relapse) ~ clin.fac$median.hERV, data = clin.fac, plotTable=TRUE)
autoplot(fit, risk.table=TRUE)
surv_pvalue(fit)

#############
n<-30
hERV.P<-rep(0,dim(hERV)[1])
for (i in 1:dim(hERV)[1]) {
  combined$median.hERV<-as.numeric(hERV[i,])
  clin.fac[,29]<-combined$median.hERV
  clin.fac[,29]<-as.factor(ifelse(clin.fac[,29] > quantile(clin.fac[,29],prob=1-n/100),'high','low'))
  fit <- surv_fit(Surv(OS, Dead) ~ clin.fac$median.hERV, data = clin.fac, plotTable=TRUE)
  hERV.P[i]<-surv_pvalue(fit)$pval
}

combined$median.hERV<-as.numeric(hERV[701,])
clin.fac[,29]<-combined$median.hERV
clin.fac[,29]<-as.factor(ifelse(clin.fac[,29] > quantile(clin.fac[,29],prob=1-n/100),'high','low'))

### hERVs HR ###
n<-30
hERV.HR<-rep(0,dim(hERV)[1])
for (i in 1:dim(hERV)[1]) {
  combined$median.hERV<-as.numeric(hERV[i,])
  clin.fac[,29]<-combined$median.hERV
  clin.fac[,29]<-as.factor(ifelse(clin.fac[,29] > quantile(clin.fac[,29],prob=1-n/100),'high','low'))
  clin.fac %>% analyse_survival(vars(RFS, Dead), by=median.hERV) -> tmp
  hERV.HR[i]<-exp(tmp$coxph$coefficients)
}

barplot(sort(hERV.HR), xlab = 'HERVs', ylab = 'HR')
abline(h=1, col = 'red', lty = 2)
#############


clin.fac %>% analyse_survival(vars(OS, Dead), median.hERV) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac %>%
       analyse_survival(vars(RFS, relapse), 
                        by=median.hERV)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("median.HERV", "median.HERV"),
                                     legend=c('right','right')
                    )) %>%
print

###################### multi level survival #######
#clin.fac.1<-clin.fac[clin.fac$median.hERV=='high',]
clin.fac.1<-clin.fac
clin.fac.1$median.hERV='CD8+/hERV+'
clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$CD8 =='high']='CD8+/hERV-'
clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$CD8 =='low']='CD8-/hERV-'
clin.fac.1$median.hERV[clin.fac$median.hERV=='high' & clin.fac$CD8 =='low']='CD8-/hERV+'
clin.fac.1$median.hERV=as.factor(clin.fac.1$median.hERV)


#clin.fac.1<-clin.fac
#clin.fac.1$median.hERV='CD8+/MSH'
#clin.fac.1$median.hERV[clin.fac$MSI=='MSS' & clin.fac$CD8 =='high']='CD8+/MSS'
#clin.fac.1$median.hERV[clin.fac$MSI=='MSS' & clin.fac$CD8 =='low']='CD8-/MSS'
#clin.fac.1$median.hERV[clin.fac$MSI=='MSH' & clin.fac$CD8 =='low']='CD8-/MSH'
#clin.fac.1$median.hERV=as.factor(clin.fac.1$median.hERV)

#clin.fac.1<-clin.fac
#clin.fac.1$median.hERV='MSH/hERV+'
#clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$MSI=='MSS']='MSS/hERV-'
#clin.fac.1$median.hERV[clin.fac$median.hERV=='high' & clin.fac$MSI=='MSS']='MSS/hERV+'
#clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$MSI=='MSH']='MSH/hERV-'
#clin.fac.1$median.hERV=as.factor(clin.fac.1$median.hERV)




clin.fac.1 %>% analyse_survival(vars(OS, Dead), median.hERV, cox_reference_level='CD8-/hERV+') -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     risk.table.height=0.15,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.1 %>%
       analyse_survival(vars(RFS, relapse), 
                        by=median.hERV, cox_reference_level='CD8-/hERV+')
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("CD8/hERV status"),
                                          legend=c('right','right')
                    )) %>%
  print
######################




clin.fac %>% analyse_multivariate(vars(OS, Dead),
                                  covariates = vars(Stage,HLA.A.normal,HLA.B.normal,HLA.C.normal,median.hERV,PD1,Age,Adj,IFNG,Sex,sidedness,MSI, CD8,CD4,CD19,Matastasis,EOMES,Treg)) -> result
#print(result)

forest_plot(result)

clin.fac %>% analyse_multivariate(vars(RFS, relapse),
                                  covariates = vars(Stage,HLA.A.normal,HLA.B.normal,HLA.C.normal,median.hERV,PD1,Age,Adj,IFNG,Sex,sidedness,MSI, CD8,CD4,CD19,Matastasis,EOMES,Treg)) -> result
#print(result)

forest_plot(result)

##########################

df<-clin.fac
map(vars(Stage,HLA.A.normal,HLA.B.normal,HLA.C.normal,median.hERV,PD1,Age,Adj,IFNG,Sex,sidedness,MSI, CD8,CD4,CD19,Matastasis,EOMES,Treg), function(by)
{
  analyse_multivariate(df,
                       vars(OS, Dead),
                       covariates = list(by))
}) -> result

forest_plot(result)

map(vars(Stage,HLA.A.normal,HLA.B.normal,HLA.C.normal,median.hERV,PD1,Age,Adj,IFNG,Sex,sidedness,MSI, CD8,CD4,CD19,Matastasis,EOMES,Treg), function(by)
{
  analyse_multivariate(df,
                       vars(RFS, Dead),
                       covariates = list(by))
}) -> result

forest_plot(result)
########################### clinopathological vs WES/WTS

clin.fac.1<-clin.fac


clin.fac.1<-clin.fac
clin.fac.1$median.hERV='clinicopathological+'
clin.fac.1$median.hERV[clin.fac$Stage=="III" | clin.fac$Age=='high' | clin.fac$sidedness=='right']='clinicopathological-'

clin.fac.1<-clin.fac
clin.fac.1$median.hERV='clinicopathological+'
clin.fac.1$median.hERV[(clin.fac$Stage=="III" | clin.fac$Age=='high' | clin.fac$sidedness=='right') & (clin.fac$median.hERV=='high' & clin.fac$CD8 =='low')]='clinicopathological-/WTS-'
clin.fac.1$median.hERV[(clin.fac$Stage=="III" | clin.fac$Age=='high' | clin.fac$sidedness=='right') & !(clin.fac$median.hERV=='high' & clin.fac$CD8 =='low')]='clinicopathological-/WTS+'

clin.fac.1<-clin.fac
clin.fac.1$median.hERV='WTS-'
clin.fac.1$median.hERV[(clin.fac$CD4=='low') & !(clin.fac$median.hERV=='high' & clin.fac$CD8 =='low')]='low/WTS+'
clin.fac.1$median.hERV[(clin.fac$CD4=="high") & !(clin.fac$median.hERV=='high' & clin.fac$CD8 =='low')]='high/WTS+'

clin.fac.1=clin.fac.1[!(clin.fac.1$median.hERV)=='WTS-',]

clin.fac.1$median.hERV=as.factor(clin.fac.1$median.hERV)

clin.fac.1 %>% analyse_survival(vars(OS, Dead), median.hERV) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, relapse), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.1 %>%
       analyse_survival(vars(RFS, relapse), 
                        by=median.hERV)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("clinicopathological/WTS status"),
                                          legend=c('right','right')
                    )) %>%
  print


############# 


clin.fac$median.hERV <- as.factor(ifelse(MEs$MEgrey>0,'high','low'))

clin.fac %>% analyse_survival(vars(OS, Dead), median.hERV) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac %>%
       analyse_survival(vars(RFS, Dead), 
                        by=median.hERV)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c(rep("MEbrown",2)),
                                          legend=c('right','right')
                    )) %>%
  print




clin.fac.1<-clin.fac
clin.fac.1$median.hERV='CD8+/hERV+'
clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$CD8 =='high']='CD8+/hERV-'
clin.fac.1$median.hERV[clin.fac$median.hERV=='low' & clin.fac$CD8 =='low']='CD8-/hERV-'
clin.fac.1$median.hERV[clin.fac$median.hERV=='high' & clin.fac$CD8 =='low']='CD8-/hERV+'
clin.fac.1$median.hERV=as.factor(clin.fac.1$median.hERV)



clin.fac.1 %>% analyse_survival(vars(OS, Dead), median.hERV) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac.1 %>%
       analyse_survival(vars(RFS, Dead), 
                        by=median.hERV)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("CD8/hERV status"),
                                          legend=c('right','right')
                    )) %>%
  print


############## Chemo resistance #######

clin.fac[clin.fac$median.hERV=='high',] %>% analyse_survival(vars(OS, Dead), Adj) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac[clin.fac$median.hERV=='low',] %>%
       analyse_survival(vars(OS, Dead), 
                        by=Adj)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("median.HERV: high\nAdjuvant Treatment", "median.HERV: low\nAdjuvant Treatment"),
                                          legend=c('right','right')
                    )) %>%
  print


##### CNV #####
armName<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.arm.amp.tsv",stringsAsFactors = FALSE)[,1]
arm.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.arm.amp.tsv",stringsAsFactors = FALSE)
idx<-which(colnames(arm.amp)%in%colnames(hERV))
arm.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.arm.amp.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(arm.amp)), names = arm.amp$armName, las = 2, ylim = c(0,1))

arm.del<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.arm.del.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(arm.del)), names = arm.del$armName, las = 2, ylim = c(0,1))


bandName<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.band.amp.tsv",stringsAsFactors = FALSE)[,1]
band.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.band.amp.tsv",stringsAsFactors = FALSE)
idx<-which(colnames(band.amp)%in%colnames(hERV))
band.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.band.amp.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(band.amp)), names = band.amp$bandName, las = 2, ylim = c(0,1), cex.names = 0.1)

band.del<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.band.del.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(band.del)), names = band.del$bandName, las = 2, ylim = c(0,1), cex.names = 0.1)

geneName<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.gene.amp.tsv",stringsAsFactors = FALSE)[,1]
gene.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.gene.amp.tsv",stringsAsFactors = FALSE)
idx<-which(colnames(gene.amp)%in%colnames(hERV))
gene.amp<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.gene.amp.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(gene.amp)), names = gene.amp$geneName, las = 2, ylim = c(0,1), cex.names = 0.1)

gene.del<-read.delim("~/Box/CRC-iMM-nonMet/CNV/FC/summary.gene.del.tsv",stringsAsFactors = FALSE)[,idx]
#barplot(rowMedians(as.matrix(gene.del)), names = gene.del$geneName, las = 2, ylim = c(0,1), cex.names = 0.1)


#barplot(rowMeans(as.matrix(arm.del>0.2)), names = arm.del$armName, las = 2, cex.names  = 0.6, ylim = c(0,0.8), ylab = 'frequency of arm deletion')
#barplot(rowMeans(as.matrix(arm.amp>0.2)), names = arm.amp$armName, las = 2, cex.names  = 0.6, ylim = c(0,0.8), ylab = 'frequency of arm amplification')

barplot(rowMeans(as.matrix(arm.del>0.2)), names = arm.del$armName, las = 2, cex.names  = 0.6, ylim = c(0,0.8), ylab = 'frequency of arm deletion')
barplot(rowMeans(as.matrix(arm.amp>0.2)), names = arm.amp$armName, las = 2, cex.names  = 0.6, ylim = c(0,0.8), ylab = 'frequency of arm amplification')

idx<-which(colnames(combined)%in%c('CD8','CD4','CD19','median.hERV'))

### arm ###
CNV.cor.amp<-cor(t(data.matrix(arm.amp)),(combined[,idx]), method = 'spearman')
CNV.cor.del<-cor(t(data.matrix(arm.del)),(combined[,idx]), method = 'spearman')


my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
CNV.cor.amp[is.na(CNV.cor.amp)]<-0
CNV.cor.del[is.na(CNV.cor.del)]<-0
rownames(CNV.cor.amp)<-armName
rownames(CNV.cor.del)<-armName

heatmap.2(CNV.cor.amp, col = my_palette, cexRow = 0.7, cexCol = 0.7, na.rm = TRUE)
heatmap.2(CNV.cor.del, col = my_palette, cexRow = 0.7, cexCol = 0.7, na.rm = TRUE)

### band ###
CNV.cor.amp<-cor(t(data.matrix(band.amp)),(combined[,idx]), method = 'spearman')
CNV.cor.del<-cor(t(data.matrix(band.del)),(combined[,idx]), method = 'spearman')


my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
CNV.cor.amp[is.na(CNV.cor.amp)]<-0
CNV.cor.del[is.na(CNV.cor.del)]<-0
rownames(CNV.cor.amp)<-bandName
rownames(CNV.cor.del)<-bandName

heatmap.2(CNV.cor.amp, col = my_palette, cexRow = 0.07, cexCol = 0.7, na.rm = TRUE)
heatmap.2(CNV.cor.del, col = my_palette, cexRow = 0.07, cexCol = 0.7, na.rm = TRUE)

### gene ###
CNV.cor.amp<-cor(t(data.matrix(gene.amp)),(combined[,idx]), method = 'spearman')
CNV.cor.del<-cor(t(data.matrix(gene.del)),(combined[,idx]), method = 'spearman')


my_palette <- colorRampPalette(c("yellow","purple"))(n = 1000)
CNV.cor.amp[is.na(CNV.cor.amp)]<-0
CNV.cor.del[is.na(CNV.cor.del)]<-0
rownames(CNV.cor.amp)<-geneName
rownames(CNV.cor.del)<-geneName

heatmap.2(CNV.cor.amp, col = my_palette, cexRow = 0.7, cexCol = 0.7, na.rm = TRUE)
heatmap.2(CNV.cor.del, col = my_palette, cexRow = 0.7, cexCol = 0.7, na.rm = TRUE)

##### hERV clustering
M<-cor(data.matrix(hERV), method = 'spearman')
ht<-heatmap.2(M , col = my_palette, cexCol = 0.78, cexRow = 1,, trace = "none")
clin.fac$hERV.cluster<-as.factor(ifelse(seq(113)%in%ht$rowInd[1:41],'cluster 1', 'cluster 2'))
clin.fac$hERV.cluster[ht$rowInd[1:41]]<-'cluster 1'
clin.fac$hERV.cluster[-ht$rowInd[1:41]]<-'cluster 2'


clin.fac %>% analyse_survival(vars(OS, Dead), hERV.cluster) -> OS.result
#clin.fac %>% analyse_survival(vars(RFS, Dead), median.hERV) -> RFS.result


default_args <- list(break.time.by="breakByYear",
                     xlab=c(".OS.months"),
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     ggtheme=ggplot2::theme_bw(10))

list(OS.result,
     clin.fac %>%
       analyse_survival(vars(RFS, relapse), 
                        by=hERV.cluster)
) %>%
  kaplan_meier_grid(nrow = 2, default_args,
                    break.time.by="breakByYear",
                    mapped_plot_args=list(legend.title=c("cluster ID", "cluster ID"),
                                          legend=c('right','right')
                    )) %>%
  print


a=combined$median.hERV[ht$rowInd[1:41]]
b=combined$median.hERV[-ht$rowInd[1:41]]

boxplot(cbind(a,b), names= c('cluster 1','cluster 2'), col = 'red', outline = FALSE, ylim = c(0,2.5))
stripchart(cbind(a,b)~cbind(rep(1,length(a)),rep(2,length(b))),
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)


