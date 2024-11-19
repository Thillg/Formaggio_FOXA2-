#!/usr/bin/Rscript
library(ggplot2)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggrepel)
library(org.Hs.eg.db)
library(msigdbr)
library(gridExtra)
library(ggpubr)

txdb = makeTxDbFromGFF("./Utils/GRCh38_genes_igenome.gtf", 
                       format="gtf", 
                       dataSource = "UCSC", 
                       organism = "Homo sapiens")


mdata = read.csv("./meta/peaks.csv")
mdata$name = paste(mdata$cellLine, mdata$treatment, mdata$target, "RECHIP_R", mdata$replicate, sep="_")
length(unique(mdata$name)) == nrow(mdata)
### Import and clean the peaks
peakList = lapply(mdata$file, readPeakFile)

peakList = lapply(peakList, function (x)x[seqnames(x) %in% c(paste0("chr", 1:22), "chrX", "chrY")])
for ( i in 1:length(peakList)){
  seqlevels(peakList[[i]])= c(paste0("chr", 1:22), "chrX", "chrY")
}
names(peakList) = mdata$name

### remove outlier peaks in length 
dir.create("./Peaks/", showWarnings = F)
pdf("./Peaks/cleaning_by_length.pdf")
for ( name in names(peakList)){
  thr = ifelse(grepl("H3K27", name), 5000, 1000)
  p<-ggplot(data.frame(x=width(peakList[[name]]))) + geom_histogram(aes(x=x, fill=x<= thr), bins=100) + geom_vline(xintercept = thr) + ggtitle(name)
  print(p)
  peakList[[name]] =peakList[[name]][width(peakList[[name]]) <= thr,] 
}
dev.off()

for ( name in names(peakList) ){
  peakList[[name]]$sample = name
  peakList[[name]]$peak_name = paste0(name, 1:length(peakList[[name]]))
}
### Set a threshold of 50 for the q-value
## GT: This is too stringent so excluded for compound analysis T - LNCaP, PC3 and 
#VCAP, LuCaP, MSKPCa1
pdf("./Peaks/cleaning_by_qvalue.pdf")
 for ( name in names(peakList)){
   thr = ifelse(grepl("H3K27", name), 25, 0)
   thr = ifelse(grepl("_AR", name), 20, thr)
   p<-ggplot(data.frame(x=peakList[[name]]$V9 )) + geom_histogram(aes(x=x, fill=x > thr), bins=100) + geom_vline(xintercept = thr) + ggtitle(name)
   print(p)
   peakList[[name]] =peakList[[name]][peakList[[name]]$V9 >= thr,] 
 }
 dev.off()

 for ( name in names(peakList) ){
   peakList[[name]]$sample = name
   peakList[[name]]$peak_name = paste0(name, 1:length(peakList[[name]]))
 }

mdata$sample = paste(mdata$cellLine , mdata$treatment, mdata$target, sep="_")
  
pdf("./Peaks/venn_replicates.pdf")
for (sam in  unique(mdata$sample) ){
  if ( sum(mdata$sample == sam) > 1 ){
    vennplot(peakList[mdata$name[mdata$sample == sam]])
    if ( sum(mdata$sample == sam) == 2 ){
      tmp = rbind( as.data.frame(peakList[[mdata$name[mdata$sample == sam][1]]]), 
                   as.data.frame(peakList[[mdata$name[mdata$sample == sam][2]]]))
      nrow(tmp)
      tmp$overlap = tmp$peak_name %in% subsetByOverlaps(peakList[[mdata$name[mdata$sample == sam][1]]], peakList[[mdata$name[mdata$sample == sam][2]]] )$peak_name |  tmp$peak_name %in% subsetByOverlaps(peakList[[mdata$name[mdata$sample == sam][2]]], peakList[[mdata$name[mdata$sample == sam][1]]] )$peak_name
      p<- ggplot(tmp) + geom_histogram(aes(x=V9, fill = overlap), bins=50)+xlim(c(0,500)) +facet_grid(rows=vars(sample))
      print(p)
    }
  }
}
dev.off()

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply( peakList, getTagMatrix, windows=promoter)

unique(mdata$name)
summary(peakList)
# 
# 
 for ( sam in unique(mdata$sample)){
   if ( sum(mdata$sample == sam) == 2 ){
     peakList[[sam]]=subsetByOverlaps(peakList[[mdata$name[mdata$sample == sam][1]]],peakList[[mdata$name[mdata$sample == sam][2]]] )
   }else {
     peakList[[sam]]=peakList[[mdata$name[mdata$sample == sam][1]]]
   }
 }
# 
 summary(peakList)
peakList = peakList[unique(mdata$sample)]
colnames(mdata)
mdata=unique(mdata[,c( "sample", "cellLine","target", "treatment") ])


saveRDS(peakList, "./objects/peakList.RDS")
saveRDS(mdata, "./objects/mdata_peaks.RDS")


peakAnnoList = lapply( peakList , annotatePeak, tssRegion=c(-3000, 3000),
                                           TxDb=txdb, annoDb="org.Hs.eg.db")
names(peakList)

order_peaks = names(peakList)

#order_peaks = c("MSKPCa_WT_FOXA1_R_1","MSKPCa_WT_FOXA2_R_1")

#OR

#order_peaks = c("PDX_145_WT_FOXA1_R_1","PDX_145_WT_FOXA2_R_1","PDX_145_WT_IGG_R_1")
#OR
# order_peaks = c("VCAP_WT_FOXA1_R_1","VCAP_WT_H3K27AC_R_1","VCAP_CTR_FOXA1_R_2",   
#                 "VCAP_CTR_H3K27AC_R_2","VCAP_F2OE_FOXA1_R_1","VCAP_F2OE_FOXA2_R_1",  
#                 "VCAP_F2OE_H3K27AC_R_1","VCAP_F2OE_FOXA1_R_2","VCAP_F2OE_FOXA2_R_2",
#                 "VCAP_F2OE_H3K27AC_R_2")


# order_peaks = c( "LnCaP_WT_FOXA1_R_1","LnCaP_WT_FOXA1_R_2","LnCaP_T_FOXA1_R_2","LnCaP_T_FOXA1_R_1", 
#                  "PC3_WT_FOXA1_R_2","PC3_WT_FOXA1_R_1","PC3_WT_FOXA2_R_2","PC3_WT_FOXA2_R_1",  
#                  "PC3_T_FOXA1_R_1","PC3_T_FOXA2_R_1","PC3_T_FOXA1_R_2","PC3_T_FOXA2_R_2") 

#OR
# 
# order_peaks = c( "LnCaP_WT_FOXA1", "LnCaP_T_WT_FOXA1",  
#                  "PC3_WT_FOXA1", "PC3_WT_FOXA2",    
#                  "PC3_T_WT_FOXA1","PC3_T_WT_FOXA2")  


all(names(peakList)%in% order_peaks)
names(peakList)[!names(peakList) %in% order_peaks]

pdf("./Peaks/plotAnnoBar.pdf")
plotAnnoBar(peakAnnoList[order_peaks])
dev.off()

### save the bed files of the peaks - not required to create a new bed folder if it already exists

#dir.create("./bed/", showWarnings = F)
for ( name in names(peakList) ){
  tmp = as.data.frame(peakList[[name]])
  write.table(tmp[,c("seqnames","start", "end","peak_name", "V9", "strand")], file=paste0("./bed/", name, ".bed"), row.names = F, col.names = F,quote=F, sep="\t" )
} 
### extract the TSS of the individual genes
tmp = as.data.frame(genes(txdb))
TSS_flank = 3000
tmp$end[tmp$strand == "+"] = tmp$start[tmp$strand=="+"] + TSS_flank
tmp$start[tmp$strand == "+"] = tmp$start[tmp$strand=="+"] - TSS_flank
tmp$start[tmp$strand == "-"] = tmp$end[tmp$strand=="-"] - TSS_flank
tmp$end[tmp$strand == "-"] = tmp$end[tmp$strand=="-"] + TSS_flank
### In case of large flanks ( not with 3k )
tmp$start[tmp$start < 0 ]= 0 
tmp$end[tmp$end < 0 ]= 0 
tmp$start = as.integer(tmp$start)
tmp$end = as.integer(tmp$end)

write.table(tmp[,c("seqnames","start", "end", "gene_id" )], file=paste0("./bed/H660_T_genes_tss_",TSS_flank,".bed"), row.names = F, col.names = F,quote=F, sep="\t" )


##GT: Plots FOXA1 Vs FOXA2 peaks from Claudio PC3, H660 and LnCAP ====
bigWig_dat = readRDS("./PeaksAnalysis_Claudio/objects_Claudio/BigWig_dat.RDS")


PC3_Corr_FOXA1_H3K = round(cor(bigWig_dat$PC3$PC3_WT_H3K27AC, bigWig_dat$PC3$PC3_WT_FOXA1), 2)
PC3_Corr_FOXA2_H3K = round(cor(bigWig_dat$PC3$PC3_WT_FOXA2, bigWig_dat$PC3$PC3_WT_H3K27AC), 2)
PC3_Corr_FOXA1_AT = round(cor(bigWig_dat$PC3$AT, bigWig_dat$PC3$PC3_WT_FOXA1), 2)
PC3_Corr_FOXA2_AT = round(cor(bigWig_dat$PC3$PC3_WT_FOXA2, bigWig_dat$PC3$AT), 2)
  
pdf("./png/PC3_WT_FOXA1_FOXA2_peaks.pdf")
p = ggplot(bigWig_dat$PC3)+geom_point(aes(x=PC3_WT_FOXA2, y=PC3_WT_FOXA1))+
  #scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,3)) + xlim(c(0,3)) + theme_classic() +
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$PC3$PC3_WT_FOXA2, bigWig_dat$PC3$PC3_WT_FOXA1), 2)))+
  ggtitle("PC3") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$PC3$PC3_WT_FOXA2, y = bigWig_dat$PC3$PC3_WT_FOXA1),
             method = "pearson", label.x = 0)
dev.off()

pdf("./png/PC3_WT_FOXA1_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$PC3)+geom_point(aes(x=PC3_WT_FOXA2, y=PC3_WT_FOXA1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,3)) + xlim(c(0,3)) + theme_classic() +
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$PC3$PC3_WT_FOXA2, bigWig_dat$PC3$PC3_WT_FOXA1), 2)))+
  ggtitle("PC3") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$PC3$PC3_WT_FOXA2, y = bigWig_dat$PC3$PC3_WT_FOXA1),
             method = "pearson", label.x = 0)


dev.off()

pdf("./png/PC3_WT_FOXA1_FOXA2_peaks_H3K.pdf")
ggplot(bigWig_dat$PC3)+geom_point(aes(x=PC3_WT_FOXA2, y=PC3_WT_FOXA1, color = PC3_WT_H3K27AC))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$PC3$PC3_WT_FOXA1))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$PC3$PC3_WT_FOXA2, bigWig_dat$PC3$PC3_WT_FOXA1), 2)))+
  xlab("PC3 FOXA2") +
  ylab("PC3 FOXA1")
dev.off()

H660_Corr_FOXA1_H3K = round(cor(bigWig_dat$H660$H660_WT_H3K27AC, bigWig_dat$H660$H660_WT_FOXA1), 2)
H660_Corr_FOXA2_H3K = round(cor(bigWig_dat$H660$H660_WT_FOXA2, bigWig_dat$H660$H660_WT_H3K27AC), 2)
H660_Corr_FOXA1_AT = round(cor(bigWig_dat$H660$AT, bigWig_dat$H660$H660_WT_FOXA1), 2)
H660_Corr_FOXA2_AT = round(cor(bigWig_dat$H660$H660_WT_FOXA2, bigWig_dat$H660$AT), 2)
  
pdf("./png/H660_WT_FOXA1_FOXA2_peaks.pdf")
ggplot(bigWig_dat$H660)+geom_point(aes(x=H660_WT_FOXA2, y=H660_WT_FOXA1))+
  #scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$H660$H660_WT_FOXA1))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$H660$H660_WT_FOXA2, bigWig_dat$H660$H660_WT_FOXA1), 2)))+
  xlab("H660 FOXA2") +
  ylab("H660 FOXA1")
dev.off()

pdf("./png/H660_WT_FOXA1_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$H660)+geom_point(aes(x=H660_WT_FOXA2, y=H660_WT_FOXA1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$H660$H660_WT_FOXA1))) + theme_classic() +
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$H660$H660_WT_FOXA2, bigWig_dat$H660$H660_WT_FOXA1), 2)))+
  ggtitle("H660") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$H660$H660_WT_FOXA2, y = bigWig_dat$H660$H660_WT_FOXA1),
             method = "pearson", label.x = 0)

dev.off()

pdf("./png/H660_WT_FOXA1_FOXA2_peaks_H3K.pdf")
ggplot(bigWig_dat$H660)+geom_point(aes(x=H660_WT_FOXA2, y=H660_WT_FOXA1, color = H660_WT_H3K27AC))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$H660$H660_WT_FOXA1))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$H660$H660_WT_FOXA2, bigWig_dat$H660$H660_WT_FOXA1), 2)))+
  xlab("H660 FOXA2") +
  ylab("H660 FOXA1")
dev.off()

LnCaP_Corr_FOXA1_H3K = round(cor(bigWig_dat$LnCaP$LnCaP_WT_H3K27AC, bigWig_dat$LnCaP$LnCaP_WT_FOXA1), 2)
LnCaP_Corr_FOXA2_H3K = round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2,bigWig_dat$LnCaP$LnCaP_WT_H3K27AC), 2)
LnCaP_Corr_FOXA1_AT = round(cor(bigWig_dat$LnCaP$AT, bigWig_dat$LnCaP$LnCaP_WT_FOXA1), 2)
LnCaP_Corr_FOXA2_AT = round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$AT), 2)


pdf("./png/LnCAP_OE_FOXA1_FOXA2_peaks.pdf")
ggplot(bigWig_dat$LnCaP)+geom_point(aes(x=LnCaP_OE_FOXA2, y=LnCaP_OE_FOXA1))+
  ylim(c(0,max(bigWig_dat$LnCaP$LnCaP_OE_FOXA1))) + theme_classic()+
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$LnCaP_OE_FOXA1), 2)))+
  xlab("LNCaP - FOXA2: FOXA2") +
  ylab("LNCaP - FOXA2: FOXA1")
dev.off()

pdf("./png/LnCAP_OE_FOXA1_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$LnCaP)+geom_point(aes(x=LnCaP_OE_FOXA2, y=LnCaP_OE_FOXA1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$LnCaP$LnCaP_OE_FOXA1))) + theme_classic()+
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$LnCaP_OE_FOXA1), 2)))+
  xlab("LNCaP-FOXA2 FOXA2") +
  ylab("LNCaP-FOXA2 FOXA1")
dev.off()

pdf("./png/LnCAP_WT_FOXA1_OE_FOXA2_peaks.pdf")
ggplot(bigWig_dat$LnCaP)+geom_point(aes(x=LnCaP_OE_FOXA2, y=LnCaP_WT_FOXA1))+
  ylim(c(0,max(bigWig_dat$LnCaP$LnCaP_WT_FOXA1))) + theme_classic()+
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$LnCaP_WT_FOXA1), 2)))+
  xlab("LNCaP-FOXA2 FOXA2") +
  ylab("LNCaP-Control FOXA1")
dev.off()

pdf("./png/LnCAP_WT_FOXA1_OE_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$LnCaP)+geom_point(aes(x=LnCaP_OE_FOXA2, y=LnCaP_WT_FOXA1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$LnCaP$LnCaP_WT_FOXA1))) + theme_classic()+
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$LnCaP_WT_FOXA1), 2)))+
  ggtitle("LNCaP") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$LnCaP$LnCaP_OE_FOXA2, y = bigWig_dat$LnCaP$LnCaP_WT_FOXA1),
             method = "pearson", label.x = 0)


dev.off()

pdf("./png/LnCAP_WT_FOXA1_OE_FOXA2_peaks_H3K.pdf")
ggplot(bigWig_dat$LnCaP)+geom_point(aes(x=LnCaP_OE_FOXA2, y=LnCaP_WT_FOXA1, color = LnCaP_WT_H3K27AC))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$LnCaP$LnCaP_WT_FOXA1))) + theme_classic()+
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$LnCaP$LnCaP_OE_FOXA2, bigWig_dat$LnCaP$LnCaP_WT_FOXA1), 2)))+
  xlab("LNCaP-FOXA2 FOXA2") +
  ylab("LNCaP-Control FOXA1")
dev.off()


##GT: Plots FOXA1 Vs FOXA2 peaks LuCaP-145, MSKCa1 and VCap ====
bigWig_dat = readRDS("./PeaksAnalysis_test/objects/BigWig_dat_LuCap.RDS")

LuCaP_Corr_FOXA1_AT = round(cor(bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1, bigWig_dat$WT$AT), 2)
LuCaP_Corr_FOXA2_AT = round(cor(bigWig_dat$WT$PDX_145_WT_FOXA2_WT_1, bigWig_dat$WT$AT), 2)


pdf("./png/LuCap_WT_FOXA1_FOXA2_peaks.pdf")
ggplot(bigWig_dat$WT)+geom_point(aes(x=PDX_145_WT_FOXA2_WT_1, y=PDX_145_WT_FOXA1_WT_1))+
  #scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$WT$PDX_145_WT_FOXA2_WT_1, bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1), 2)))+
  xlab("LuCap-145 FOXA2") +
  ylab("LuCap-145 FOXA1")
dev.off()

pdf("./png/LuCap_WT_FOXA1_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$WT)+geom_point(aes(x=PDX_145_WT_FOXA2_WT_1, y=PDX_145_WT_FOXA1_WT_1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1))) + theme_classic() +
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$WT$PDX_145_WT_FOXA2_WT_1, bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1), 2)))+
  ggtitle("LuCaP-145") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$WT$PDX_145_WT_FOXA2_WT_1, y = bigWig_dat$WT$PDX_145_WT_FOXA1_WT_1),
             method = "pearson", label.x = 0)


dev.off()


bigWig_dat = readRDS("./PeaksAnalysis_test/objects/BigWig_dat_MSK.RDS")

MSK_Corr_FOXA1_AT = round(cor(bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1, bigWig_dat$WT$AT), 2)
MSK_Corr_FOXA2_AT = round(cor(bigWig_dat$WT$MSKPCa1_WT_FOXA2_WT_1, bigWig_dat$WT$AT), 2)


pdf("./png/MSKPCa1_WT_FOXA1_FOXA2_peaks.pdf")
ggplot(bigWig_dat$WT)+geom_point(aes(x=MSKPCa1_WT_FOXA2_WT_1, y=MSKPCa1_WT_FOXA1_WT_1))+
  ylim(c(0,max(bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$WT$MSKPCa1_WT_FOXA2_WT_1, bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1), 2)))+
  xlab("MSKPCa1 FOXA2") +
  ylab("MSKPCa1 FOXA1")
dev.off()

pdf("./png/MSKPCa1_WT_FOXA1_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$WT)+geom_point(aes(x=MSKPCa1_WT_FOXA2_WT_1, y=MSKPCa1_WT_FOXA1_WT_1, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1))) + theme_classic() +
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$WT$MSKPCa1_WT_FOXA2_WT_1, bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1), 2)))+
  ggtitle("MSKPCa1") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$WT$MSKPCa1_WT_FOXA2_WT_1, y = bigWig_dat$WT$MSKPCa1_WT_FOXA1_WT_1),
             method = "pearson", label.x = 0)

dev.off()


bigWig_dat = readRDS("./PeaksAnalysis_test/objects/BigWig_dat_VCAP.RDS")

VCAP_Corr_FOXA1_H3K = round(cor(bigWig_dat$OE_WT$VCAP_CTR_H3K27AC_OE_WT_2, bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2), 2)
VCAP_Corr_FOXA2_H3K = round(cor(bigWig_dat$OE_WT$VCAP_CTR_H3K27AC_OE_WT_2, bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2), 2)
VCAP_Corr_FOXA1_AT = round(cor(bigWig_dat$OE_WT$AT, bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2), 2)
VCAP_Corr_FOXA2_AT = round(cor(bigWig_dat$OE_WT$AT, bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2), 2)


pdf("./png/VCAP_WT_FOXA1_OE_FOXA2_peaks.pdf")
ggplot(bigWig_dat$OE_WT)+geom_point(aes(x=VCAP_F2OE_FOXA2_OE_WT_2, y=VCAP_CTR_FOXA1_OE_WT_2))+
  ylim(c(0,max(bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2))) + theme_classic() +
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2, bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2), 2)))+
  xlab("VCAP-FOXA2 FOXA2") +
  ylab("VCAP-Control FOXA1")
dev.off()

pdf("./png/VCAP_WT_FOXA1_OE_FOXA2_peaks_AT.pdf")
p = ggplot(bigWig_dat$OE_WT)+geom_point(aes(x=VCAP_F2OE_FOXA2_OE_WT_2, y=VCAP_CTR_FOXA1_OE_WT_2, color = AT))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2))) + theme_classic()+
  #ggtitle(paste0("CORR: ", round(cor(bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2, bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2), 2)))+
  ggtitle("VCAP") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  xlab("Peak score FOXA2") +
  ylab("Peak score FOXA1")

p + stat_cor(aes(x = bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2, y = bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2),
             method = "pearson", label.x = 0)
dev.off()

pdf("./png/VCAP_WT_FOXA1_OE_FOXA2_peaks_H3K.pdf")
ggplot(bigWig_dat$OE_WT)+geom_point(aes(x=VCAP_F2OE_FOXA2_OE_WT_2, y=VCAP_CTR_FOXA1_OE_WT_2, color = VCAP_CTR_H3K27AC_OE_WT_2 ))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0.5) +
  ylim(c(0,max(bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2))) + theme_classic()+
  ggtitle(paste0("CORR: ", round(cor(bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2, bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2), 2)))+
  xlab("VCAP-FOXA2 FOXA2") +
  ylab("VCAP-Control FOXA1")
dev.off()

corr_tab = as.data.frame(rbind(PC3_Corr_FOXA1_H3K,PC3_Corr_FOXA2_H3K, PC3_Corr_FOXA1_AT,PC3_Corr_FOXA2_AT,
                 H660_Corr_FOXA1_H3K,H660_Corr_FOXA2_H3K,H660_Corr_FOXA1_AT,H660_Corr_FOXA2_AT,
                 LnCaP_Corr_FOXA1_H3K,LnCaP_Corr_FOXA2_H3K,LnCaP_Corr_FOXA1_AT,LnCaP_Corr_FOXA2_AT,
                 LuCaP_Corr_FOXA1_AT,LuCaP_Corr_FOXA2_AT,MSK_Corr_FOXA1_AT,MSK_Corr_FOXA2_AT,
                 VCAP_Corr_FOXA1_H3K,VCAP_Corr_FOXA2_H3K,VCAP_Corr_FOXA1_AT,VCAP_Corr_FOXA2_AT))


names(corr_tab) = "Corr_Coeff"

pdf("./png/CORR_FOXA1_FOXA2_H3K_AT.pdf")
grid.table(corr_tab)
dev.off()
