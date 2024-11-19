#!/usr/bin/Rscript
library(ggplot2)
library(clusterProfiler)
library(msigdbr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggrepel)
library(org.Hs.eg.db)
library(stringr)
library(dplyr)

txdb = makeTxDbFromGFF("./Utils/GRCh38_genes_igenome.gtf", 
                       format="gtf", 
                       dataSource = "UCSC", 
                       organism = "Homo sapiens")

#setwd("~/work/FOXA2_final/hg38/")

peakList = readRDS("./objects/peakList.RDS")
mdata = readRDS("./objects/mdata_peaks.RDS")

summary(peakList)
summary(mdata)

### Load the BigWigAverage files
bwmeta = read.csv("./meta/bigWigAverage.csv")
for ( f in bwmeta$file) {
  if ( ! file.exists(f)){
    print(f)
  }
}
bwmeta$name = paste(bwmeta$cellLine, bwmeta$treatment, bwmeta$target, bwmeta$comparison,  bwmeta$replicate, sep="_")
bwmeta$sample = paste(bwmeta$cellLine, bwmeta$treatment, bwmeta$target, bwmeta$comparison,bwmeta$replicate, sep="_")
files = bwmeta$file
names(files)= bwmeta$name
bigWigAvg = lapply(files, function(x){ read.table(x, col.names = c("name", "size", "covered", "sum", "mean0", "mean")) })
score_col = "mean0"
bigWig_dat = list()

for ( comparison in unique(bwmeta$comparison)){
  samples= unique(bwmeta$sample[bwmeta$comparison == comparison ])
  for ( sample in samples){
    tmp = NULL
    mask = bwmeta$sample == sample 
    if ( sum(mask) > 1 ){
      for ( name in bwmeta$name[mask]){
        if ( is.null(tmp)){
          tmp = data.frame(name = bigWigAvg[[name]]$name)
        }
        tmp[,name]= bigWigAvg[[name]][,score_col]
      }
      tmp[,sample] = apply(tmp[,bwmeta$name[bwmeta$sample == sample]],1, mean)
      tmp = tmp[,c("name", sample)]
    } else {
      name = bwmeta$name[mask]
      tmp = data.frame(name = bigWigAvg[[name]]$name)
      tmp[,sample] =  bigWigAvg[[name]][,score_col]
    }
    # possible rescaling
    # tmp[,sample] = log2(tmp[,sample]+1)
    # tmp[,sample]  = tmp[,sample] / mean(tmp[,sample])
    if ( is.null(bigWig_dat[[comparison]])){
      bigWig_dat[[comparison]]=tmp
    } else {
      bigWig_dat[[comparison]] = merge(bigWig_dat[[comparison]],tmp, by="name")   
    }
  }
}



seq_files = list.files("./bigWigAverage/", pattern = "sequences*")
seq_obj = list()
for ( fname in seq_files ){
  seq_dat = read.table(paste0("./bigWigAverage/",fname), col.names =c("name","seq") )
  seq_dat$C_freq = str_count(seq_dat$seq, "C") / nchar(seq_dat$seq)
  seq_dat$G_freq = str_count(seq_dat$seq, "G") / nchar(seq_dat$seq)
  seq_dat$A_freq = str_count(seq_dat$seq, "A") / nchar(seq_dat$seq)
  seq_dat$T_freq = str_count(seq_dat$seq, "T") / nchar(seq_dat$seq)
  seq_dat$CG = seq_dat$C_freq + seq_dat$G_freq
  seq_dat$AT = seq_dat$T_freq + seq_dat$A_freq
  seq_obj[[fname]] = seq_dat
}


### OE_WT ====

basedir= "./PeakAnalysis/"
dir.create(basedir, showWarnings = F, recursive = T)

#Rep2 VCAP
bigWig_dat$OE_WT$FOXA1_change = log2((bigWig_dat$OE_WT$VCAP_F2OE_FOXA1_OE_WT_2+1)/ (bigWig_dat$OE_WT$VCAP_CTR_FOXA1_OE_WT_2+1))
bigWig_dat$OE_WT$H3K27AC_change = log2((bigWig_dat$OE_WT$VCAP_F2OE_H3K27AC_OE_WT_2+1)/ (bigWig_dat$OE_WT$VCAP_CTR_H3K27AC_OE_WT_2+1))
bigWig_dat$OE_WT$FOXA2_change = log2((bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2+1)/ (bigWig_dat$OE_WT$VCAP_WT_INPUT_OE_WT_2 +1))
bigWig_dat$OE_WT$FOX_ratio = bigWig_dat$OE_WT$VCAP_F2OE_FOXA1_OE_WT_2 / (bigWig_dat$OE_WT$VCAP_F2OE_FOXA1_OE_WT_1  + bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_2 )
bigWig_dat$OE_WT$FOX_ratio[is.na(bigWig_dat$OE_WT$FOX_ratio)]=0
bigWig_dat$OE_WT = merge(bigWig_dat$OE_WT, seq_obj$sequences_VCAP_FOXA1_FOXA2_H3K_R_2.tsv, by="name")

saveRDS(bigWig_dat, "./objects/bigWig_dat_VCAP.RDS")

#LuCap
bigWig_dat$WT = merge(bigWig_dat$WT, seq_obj$sequences_PDX_145_FOXA1_FOXA2_R_1.tsv, by="name")
saveRDS(bigWig_dat, "./objects/bigWig_dat_LuCap.RDS")

#MSKPCa1
bigWig_dat$WT = merge(bigWig_dat$WT, seq_obj$sequences_MSKPCa_FOXA1_FOXA2_R_1.tsv, by="name")
saveRDS(bigWig_dat, "./objects/bigWig_dat_MSK.RDS")

### identify Synergetic and Anatagonist peaks in VCAP WT

bigWig_dat$OE_WT$type = "-"
synergy_thr = 0
bigWig_dat$OE_WT$type[ bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_1 > 0 & bigWig_dat$OE_WT$FOXA1_change > synergy_thr &  bigWig_dat$OE_WT$H3K27AC_change > synergy_thr ] = "Synergy"
bigWig_dat$OE_WT$type[ bigWig_dat$OE_WT$VCAP_F2OE_FOXA2_OE_WT_1 > 0 & bigWig_dat$OE_WT$FOXA1_change < -synergy_thr &  bigWig_dat$OE_WT$H3K27AC_change < -synergy_thr ] = "Antagonism"

pdf("./PeakAnalysis/AntagonismSynergy_VCAP_R1.pdf")
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA1_change, y=H3K27AC_change, color=type), alpha=0.8)
# ggplot(bigWig_dat$F1_F2)+ geom_boxplot(aes(x=type, y=AT, color=type), alpha=0.8)
# ggplot(bigWig_dat$F1_F2)+ geom_boxplot(aes(x=type, y=AT, fill=AR_overlap), alpha=0.8)
# ggplot(bigWig_dat$F1_F2)+ geom_bar(aes(x=type, fill=AR_overlap))
dev.off()

pdf("./PeakAnalysis/FOXA1_OE_Vs_WT_VCAP_R1.pdf")
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=VCAP_WT_FOXA1_OE_WT_1, y=VCAP_F2OE_FOXA1_OE_WT_1, color=VCAP_F2OE_H3K27AC_OE_WT_1), 
                                     alpha=0.8) + 
  scale_color_gradient2(low="blue", high="red", midpoint = 0)
# ggplot(bigWig_dat$F1_F2)+ geom_boxplot(aes(x=type, y=AT, color=type), alpha=0.8)
# ggplot(bigWig_dat$F1_F2)+ geom_boxplot(aes(x=type, y=AT, fill=AR_overlap), alpha=0.8)
# ggplot(bigWig_dat$F1_F2)+ geom_bar(aes(x=type, fill=AR_overlap))
dev.off()

peaks = read.table("./merged_bed/VCAP_FOXA1_FOXA2_H3K_R_1.bed", col.names = c("chr","start","end","name"))
synergy = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Synergy"],]
antagonism = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Antagonism"],]


peaks_WT_vs_OE = peaks
write.table(synergy, "./PeakAnalysis/VCAP_synergy_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(antagonism, "./PeakAnalysis/VCAP_antagonism_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)

#### identify increased and decreased peaks in resistant cell line
synergy_thr=0
bigWig_dat$OE_WT$type = "-"
bigWig_dat$OE_WT$type[bigWig_dat$OE_WT$FOXA1_change > synergy_thr & bigWig_dat$OE_WT$FOXA2_change > synergy_thr ] = "Increased"
#& bigWig_dat$OE_WT$H3K27AC_change > synergy_thr
bigWig_dat$OE_WT$type[bigWig_dat$OE_WT$FOXA1_change < -synergy_thr & bigWig_dat$OE_WT$FOXA2_change < -synergy_thr ] = "Decreased"
#& bigWig_dat$OE_WT$H3K27AC_change > -synergy_thr
bigWig_dat$OE_WT$type[bigWig_dat$OE_WT$FOXA1_change < -synergy_thr & bigWig_dat$OE_WT$FOXA2_change > synergy_thr ] = "Competition_F2"
#& bigWig_dat$OE_WT$H3K27AC_change > -synergy_thr
bigWig_dat$OE_WT$type[bigWig_dat$OE_WT$FOXA1_change > synergy_thr & bigWig_dat$OE_WT$FOXA2_change < -synergy_thr  ] = "Competition_F1"
#&  bigWig_dat$OE_WT$H3K27AC_change > synergy_thr

pdf("./PeakAnalysis/IncreaseDecreased_VCAP_R1.pdf")
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA1_change, y=FOXA2_change, color=type))
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA1_change, y=FOXA2_change, color=H3K27AC_change))+ 
  scale_color_gradient2(low="blue", high="red", midpoint = 0)
#ggplot(bigWig_dat$ROE_vs_OE)+ geom_boxplot(aes(x=type, y=AT, color=type))
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA2_change, y=H3K27AC_change, color=FOXA1_change))+ 
  scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA1_change, y=FOXA2_change, color=H3K27AC_change))+ 
  facet_grid(rows = vars(type))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$OE_WT)+ geom_point(aes(x=FOXA1_change, y=H3K27AC_change, color=FOXA2_change))+ 
  facet_grid(rows = vars(type))+
  scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$OE_WT)+ geom_bar(aes(x=type, fill=H3K27AC_change > 0 ))
ggplot(bigWig_dat$OE_WT)+ geom_bar(aes(x=type, fill=FOXA1_change > 0 ))
dev.off()


peaks = read.table("./merged_bed/VCAP_FOXA1_FOXA2_H3K_R_1.bed",col.names = c("chr","start","end","name"))
increased = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Decreased"],]
Competition_F2 = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Competition_F2"],]
Competition_F1 = peaks[peaks$name %in% bigWig_dat$OE_WT$name[bigWig_dat$OE_WT$type == "Competition_F1"],]

#peaks_ROE_vs_OE = read.table("./merged_bed/LnCaP_R_OE.bed",col.names = c("chr","start","end","name"))

write.table(increased, "./PeakAnalysis/VCAP_increased_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/VCAP_decreased_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(Competition_F2, "./PeakAnalysis/VCAP_Competition_F2_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(Competition_F1, "./PeakAnalysis/VCAP_Competition_F1_R_1.bed", sep="\t", quote=F, row.names = F, col.names = F)

# peakListProcessed = list(
#   synergy = GRanges(synergy),
#   antagonism = GRanges(antagonism),
#   increased = GRanges(increased),
#   decreased = GRanges(decreased),
#   Competition_F2 = GRanges(Competition_F2),
#   Competition_F1 = GRanges(Competition_F1),
#   OE_vs_WT = GRanges(peaks_WT_vs_OE),
#   ROE_vs_OE = GRanges(peaks_ROE_vs_OE)
# )


### FoldChange T Vs DMSO
#PC3_FOXA1_R1 ====
basedir= "./PeakAnalysis/PC3_FOXA1_R1/"
dir.create(basedir, showWarnings = F, recursive = T)

bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_PC3_FOXA1_R_1.tsv , by="name")


peaks = read.table('./merged_bed/PC3_FOXA1_R_1.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change = log2((bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1+1)/ (bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1+1))


#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1 = "-"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 >= 0 ] = "Increased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1 == "Increased" & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0 &
                                         bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5] = "Decreased_H3K"

#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"


pdf(paste0(basedir, "/PC3_FOXA1_IncreasedDecreased_and_AT.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_PC3_FOXA1_R1 != "-") %>%
  group_by(type_PC3_FOXA1_R1) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA1_R1 != "-"), aes(x=AT, y=..count.., fill=type_PC3_FOXA1_R1)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_PC3_FOXA1_R1))

dev.off()

increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R1 == "Decreased"],]

write.table(increased, "./PeakAnalysis/PC3_FOXA1_R1_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/PC3_FOXA1_R1_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL

#PC3_FOXA1_R2 ====
basedir= "./PeakAnalysis/PC3_FOXA1_R2/"
dir.create(basedir, showWarnings = F, recursive = T)
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_PC3_FOXA1_R2.tsv , by="name")

peaks = read.table('./merged_bed/PC3_FOXA1_R_2.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA1_change = log2((bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_2+1)/ (bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_2+1))

#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2 = "-"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Increased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5  ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5  ] = "Decreased_H3K"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"

pdf(paste0(basedir, "/PC3_FOXA1_IncreasedDecreased_and_AT_R2.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_PC3_FOXA1_R2 != "-") %>%
  group_by(type_PC3_FOXA1_R2) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA1_R2 != "-"), aes(x=AT, y=..count.., fill=type_PC3_FOXA1_R2)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_PC3_FOXA1_R2))
dev.off()


increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA1_R2 == "Decreased"],]

write.table(increased, "./PeakAnalysis/PC3_FOXA1_R2_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/PC3_FOXA1_R2_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL


#PC3_FOXA2_R1 ====
basedir= "./PeakAnalysis/PC3_FOXA2_R1/"
dir.create(basedir, showWarnings = F, recursive = T)
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_PC3_FOXA2_R_1.tsv , by="name")

peaks = read.table('./merged_bed/PC3_FOXA2_R_1.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA2_change = log2((bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_1+1)/ (bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_1+1))

#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1 = "-"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA2_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5] = "Increased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA2_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA2_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5  ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA2_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5  ] = "Decreased_H3K"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"

pdf(paste0(basedir, "/PC3_FOXA2_IncreasedDecreased_and_AT_R1.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_PC3_FOXA2_R1 != "-") %>%
  group_by(type_PC3_FOXA2_R1) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA2_R1 != "-"), aes(x=AT, y=..count.., fill=type_PC3_FOXA2_R1)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_PC3_FOXA2_R1))
dev.off()


increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R1 == "Decreased"],]

write.table(increased, "./PeakAnalysis/PC3_FOXA2_R1_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/PC3_FOXA2_R1_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL

#PC3_FOXA2_R2 ====
basedir= "./PeakAnalysis/PC3_FOXA2_R2/"
dir.create(basedir, showWarnings = F, recursive = T)

bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_PC3_FOXA2_R_2.tsv , by="name")

peaks = read.table('./merged_bed/PC3_FOXA2_R_2.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA2_change = log2((bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_2+1)/ (bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_2+1))

#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2 = "-"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA2_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5] = "Increased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA2_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA2_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2[bigWig_dat$T_Vs_DMSO$PC3_R2_FOXA2_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA2_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA2_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Decreased_H3K"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"

pdf(paste0(basedir, "/PC3_FOXA2_IncreasedDecreased_and_AT_R2.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_PC3_FOXA2_R2 != "-") %>%
  group_by(type_PC3_FOXA2_R2) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA2_R2 != "-"), aes(x=AT, y=..count.., fill=type_PC3_FOXA2_R2)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_PC3_FOXA2_R2))
dev.off()


increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_PC3_FOXA2_R2 == "Decreased"],]

write.table(increased, "./PeakAnalysis/PC3_FOXA2_R2_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/PC3_FOXA2_R2_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL


#LnCap_FOXA1_R1 ====
basedir= "./PeakAnalysis/LnCap_FOXA1_R1/"
dir.create(basedir, showWarnings = F, recursive = T)

bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_LnCap_FOXA1_R_1.tsv , by="name")

peaks = read.table('./merged_bed/LnCap_FOXA1_R_1.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$LnCap_R1_FOXA1_change = log2((bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_1+1)/ (bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_1+1))

#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1 = "-"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1[bigWig_dat$T_Vs_DMSO$LnCap_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Increased"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1[bigWig_dat$T_Vs_DMSO$LnCap_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1[bigWig_dat$T_Vs_DMSO$LnCap_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1[bigWig_dat$T_Vs_DMSO$LnCap_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_1 >= 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Decreased_H3K"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"

pdf(paste0(basedir, "/LnCap_FOXA1_IncreasedDecreased_and_AT_R1.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_LnCap_FOXA1_R1 != "-") %>%
  group_by(type_LnCap_FOXA1_R1) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_LnCap_FOXA1_R1 != "-"), aes(x=AT, y=..count.., fill=type_LnCap_FOXA1_R1)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_LnCap_FOXA1_R1))
dev.off()

increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R1 == "Decreased"],]

write.table(increased, "./PeakAnalysis/LnCap_FOXA1_R1_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/LnCap_FOXA1_R1_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL

#LnCap_FOXA1_R2
basedir= "./PeakAnalysis/LnCap_FOXA1_R2/"
dir.create(basedir, showWarnings = F, recursive = T)

bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, seq_obj$sequences_LnCap_FOXA1_R_2.tsv , by="name")

peaks = read.table('./merged_bed/LnCap_FOXA1_R_2.bed', sep="\t", col.names = c("chr", "start", "end", "name"))
bigWig_dat$T_Vs_DMSO = merge(bigWig_dat$T_Vs_DMSO, peaks, by="name")
bigWig_dat$T_Vs_DMSO$LnCap_R2_FOXA1_change = log2((bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_2+1)/ (bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_2+1))

#### identify increased and decreased peaks with treatement
synergy_thr=0.5
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2 = "-"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2[bigWig_dat$T_Vs_DMSO$LnCap_R2_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Increased"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2[bigWig_dat$T_Vs_DMSO$LnCap_R2_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Increased_H3K"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2[bigWig_dat$T_Vs_DMSO$LnCap_R2_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 < 0.5 ] = "Decreased"
bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2[bigWig_dat$T_Vs_DMSO$LnCap_R2_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$LnCaP_T_FOXA1_T_Vs_DMSO_2 >= 0 & bigWig_dat$T_Vs_DMSO$LnCaP_DMSO_FOXA1_T_Vs_DMSO_2 > 0 & bigWig_dat$T_Vs_DMSO$PC3_WT_H3K27AC_T_Vs_DMSO_1 >= 0.5 ] = "Decreased_H3K"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change > synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 > 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 <= 0 ] = "New"
#bigWig_dat$T_Vs_DMSO$type[bigWig_dat$T_Vs_DMSO$PC3_R1_FOXA1_change < -synergy_thr & bigWig_dat$T_Vs_DMSO$PC3_T_FOXA1_T_Vs_DMSO_1 <= 0 & bigWig_dat$T_Vs_DMSO$PC3_DMSO_FOXA1_T_Vs_DMSO_1 > 0  ] = "Disappeared"

pdf(paste0(basedir, "/LnCap_FOXA1_IncreasedDecreased_and_AT_R2.pdf"))
dummy <- bigWig_dat$T_Vs_DMSO %>%
  filter(type_LnCap_FOXA1_R2 != "-") %>%
  group_by(type_LnCap_FOXA1_R2) %>%
  summarize(mean = mean(AT))

ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_LnCap_FOXA1_R2 != "-"), aes(x=AT, y=..count.., fill=type_LnCap_FOXA1_R2)) +
  geom_density(alpha=0.4) + geom_vline(data = dummy, aes(xintercept = mean, color = type_LnCap_FOXA1_R2))
dev.off()

increased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2 == "Increased"],]
decreased = peaks[peaks$name %in% bigWig_dat$T_Vs_DMSO$name[bigWig_dat$T_Vs_DMSO$type_LnCap_FOXA1_R2 == "Decreased"],]

write.table(increased, "./PeakAnalysis/LnCap_FOXA1_R2_increased.bed", sep="\t", quote=F, row.names = F, col.names = F)
write.table(decreased, "./PeakAnalysis/LnCap_FOXA1_R2_decreased.bed", sep="\t", quote=F, row.names = F, col.names = F)

bigWig_dat$T_Vs_DMSO$chr = bigWig_dat$T_Vs_DMSO$start = bigWig_dat$T_Vs_DMSO$end = NULL
bigWig_dat$T_Vs_DMSO$seq = bigWig_dat$T_Vs_DMSO$C_freq = bigWig_dat$T_Vs_DMSO$G_freq =
  bigWig_dat$T_Vs_DMSO$A_freq = bigWig_dat$T_Vs_DMSO$T_freq = bigWig_dat$T_Vs_DMSO$CG = bigWig_dat$T_Vs_DMSO$AT = NULL


saveRDS(bigWig_dat, "./objects/bigWigAverage_thresh0.5.RDS")

bigWig_dat = readRDS("./objects/bigWigAverage_thresh0.5.RDS")

#ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA1_R1 != "-"))+ geom_bar(aes(x=type_PC3_FOXA1_R1, fill = type_PC3_FOXA1_R1))

#Plots ====
pdf("./PeakAnalysis/IncreaseDecreased_PC3_FOXA1_H3K.pdf")
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=PC3_DMSO_FOXA1_T_Vs_DMSO_1, y=PC3_T_FOXA1_T_Vs_DMSO_1, color=type_PC3_FOXA1_R1),
                                         alpha = 1/5, size = 1) +
  xlab("PC3_DMSO_FOXA1_R1") + ylab("PC3_T_FOXA1_R1") + theme_classic()
  #facet_grid(rows = vars(type))
  #scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA1_R1 != "-"))+ geom_bar(aes(x=type_PC3_FOXA1_R1, fill = type_PC3_FOXA1_R1),alpha = 1/5) +
                                                                              scale_fill_manual(values=c("blue", "green","Orange","Pink"))
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=PC3_DMSO_FOXA1_T_Vs_DMSO_2, y=PC3_T_FOXA1_T_Vs_DMSO_2, color=type_PC3_FOXA1_R2),
                                         alpha = 1/5, size = 1) +
  xlab("PC3_DMSO_FOXA1_R2") + ylab("PC3_T_FOXA1_R2") + theme_classic()
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA1_R2 != "-"))+ geom_bar(aes(x=type_PC3_FOXA1_R2, fill = type_PC3_FOXA1_R2), alpha = 1/5)+
                                                                              scale_fill_manual(values=c("blue", "green","Orange","Pink"))
dev.off()

pdf("./PeakAnalysis/IncreaseDecreased_PC3_FOXA2_H3K.pdf")
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=PC3_DMSO_FOXA2_T_Vs_DMSO_1, y=PC3_T_FOXA2_T_Vs_DMSO_1, color=type_PC3_FOXA2_R1),
                                         alpha = 1/5, size = 1) +
  xlab("PC3_DMSO_FOXA2_R1") + ylab("PC3_T_FOXA2_R1") + theme_classic()
#facet_grid(rows = vars(type))
#scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA2_R1 != "-"))+ geom_bar(aes(x=type_PC3_FOXA2_R1, fill = type_PC3_FOXA2_R1),alpha = 1/5) +
  scale_fill_manual(values=c("blue", "green","Orange","Pink"))
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=PC3_DMSO_FOXA2_T_Vs_DMSO_2, y=PC3_T_FOXA2_T_Vs_DMSO_2, color=type_PC3_FOXA2_R2),
                                         alpha = 1/5, size = 1) +
  xlab("PC3_DMSO_FOXA2_R2") + ylab("PC3_T_FOXA2_R2") + theme_classic()
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_PC3_FOXA2_R2 != "-"))+ geom_bar(aes(x=type_PC3_FOXA2_R2, fill = type_PC3_FOXA2_R2), alpha = 1/5)+
  scale_fill_manual(values=c("blue", "green","Orange","Pink"))
dev.off()

pdf("./PeakAnalysis/IncreaseDecreased_LnCap_FOXA1_H3K.pdf")
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=LnCaP_DMSO_FOXA1_T_Vs_DMSO_1, y=LnCaP_T_FOXA1_T_Vs_DMSO_1, color=type_LnCap_FOXA1_R1),
                                         alpha = 1/5, size = 1) +
  xlab("LnCap_DMSO_FOXA1_R1") + ylab("LnCap_T_FOXA1_R1") + theme_classic()
#facet_grid(rows = vars(type))
#scale_color_gradient2(low="blue", high="red", midpoint = 0)
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_LnCap_FOXA1_R1 != "-"))+ geom_bar(aes(x=type_LnCap_FOXA1_R1, fill = type_LnCap_FOXA1_R1),alpha = 1/5) +
  scale_fill_manual(values=c("blue", "green","Orange","Pink"))
ggplot(bigWig_dat$T_Vs_DMSO)+ geom_point(aes(x=LnCaP_DMSO_FOXA1_T_Vs_DMSO_2, y=LnCaP_T_FOXA1_T_Vs_DMSO_2, color=type_LnCap_FOXA1_R2),
                                         alpha = 1/5, size = 1) +
  xlab("LnCap_DMSO_FOXA1_R2") + ylab("LnCap_T_FOXA1_R2") + theme_classic()
ggplot(bigWig_dat$T_Vs_DMSO %>% filter(type_LnCap_FOXA1_R2 != "-"))+ geom_bar(aes(x=type_LnCap_FOXA1_R2, fill = type_LnCap_FOXA1_R2), alpha = 1/5)+
  scale_fill_manual(values=c("blue", "green","Orange","Pink"))
dev.off()
