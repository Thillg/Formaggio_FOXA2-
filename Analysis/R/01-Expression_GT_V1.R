#!/usr/bin/Rscript
library("ggplot2")
#BiocManager::install("tximport")
library(tximport)
library(DESeq2)
library(msigdbr)
library(clusterProfiler)
library("RColorBrewer")
library(pheatmap)
library(dplyr)


### Functions 
performDE = function(md, lev_names ){
  fnames = md$file
  names(fnames) = md$name
  txi = tximport(fnames, type="salmon", tx2gene = tx2gene)
  md$condition = factor(md$condition, levels= lev_names)
  if ( length(unique(md$batch)) > 1 ) {
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = md,
                                    design = ~batch + condition)   
  } else {
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = md,
                                    design = ~ condition)
  }
  
  dds = dds[rowSums(counts(dds) > 10) >= 3 ,]
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}

performDEtx = function(md, lev_names ){
  fnames = md$file
  names(fnames) = md$name
  txi = tximport(fnames, type="salmon",txOut = T)
  md$condition = factor(md$condition, levels= lev_names)
  if ( length(unique(md$batch)) > 1 ) {
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = md,
                                    design = ~batch + condition)   
  } else {
    dds <- DESeqDataSetFromTximport(txi,
                                    colData = md,
                                    design = ~ condition)
  }
  
  dds = dds[rowSums(counts(dds) > 10) >= 3 ,]
  dds <- DESeq(dds)
  res <- results(dds)
  return(res)
}
pcaPlot <- function(indata=NULL,  colorRef=NULL,  ngenes=500, pca=NULL){
  if ( is.null(pca) ){
    rv_all= rowVars(indata)
    pcagenes_all = rownames(indata)[order(rv_all, decreasing = TRUE)]
    selected_all = pcagenes_all[1:ngenes]
    pca = prcomp( t(indata[selected_all,]))  
  }
  percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
  print(ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], color = colorRef ))+geom_point(aes(x=PC1, y=PC2, color=color), size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%")) +
          theme_classic()+ theme (legend.title = element_blank()) + coord_fixed())
  return(pca)
}

### 

mdata = read.csv("./meta/expression.csv")
tx2gene = read.table("./meta/salmon_tx2gene.tsv")
mdata$name = paste(mdata$cellLine, mdata$condition, mdata$replicate, sep="_")

fnames = mdata$file
names(fnames) = mdata$name
for ( f in fnames ){
  if ( ! file.exists(f)){
    print(f)
  }
}

txi = tximport(fnames, type = "salmon", tx2gene = tx2gene)
mdata$group = paste(mdata$cellLine, mdata$condition, sep="_")
mdata$batch = factor(mdata$batch)

dds = DESeqDataSetFromTximport(txi, colData = mdata, design=~ group)
dds = DESeq(dds)

rawcount = data.frame(assay(dds))

ncount = assay(vst(dds, blind=F))
base_dir= "./expression_analysis/"
dir.create(base_dir, showWarnings = F)
pdf(paste0(base_dir, "PCA_DKD_T.pdf"))
pca = pcaPlot(ncount, paste(mdata$cellLine, mdata$condition, sep="_"), 1000 ) 
dev.off()

pdf(paste0(base_dir, "PCA_batch_effect.pdf"))
pca = pcaPlot(ncount, paste(mdata$cellLine, mdata$condition, sep="_"), 1000 ) 
percentVar <- setNames(object = round(100* pca$sdev^2/sum(pca$sdev^2)), nm = colnames(pca$x))
print(ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name = mdata$cellLine, 
                        batch = factor(mdata$batch) ))+geom_point(aes(x=PC1, y=PC2, color=name, shape=batch), size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%")))
print(ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name = mdata$cellLine, condition = mdata$condition,
                        batch = factor(mdata$batch) ))+geom_point(aes(x=PC1, y=PC2, color=condition, shape=batch), size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%")))

vsd = vst(dds, blind=F)
mat <- assay(vsd)
mm <- model.matrix(~group, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
assay(vsd) <- mat
ncount = mat
pca = pcaPlot(assay(vsd), paste(mdata$cellLine, mdata$condition, sep="_"), 1000 )
print(ggplot(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], name =  mdata$cellLine, 
                        batch = factor(mdata$batch) ))+geom_point(aes(x=PC1, y=PC2, color=name, shape=batch), size=3)+ylab(paste("PC2",percentVar[2],"%"))+xlab(paste("PC1", percentVar[1],"%")))

dev.off()

dir.create("./objects", showWarnings = F)
saveRDS(ncount, "./objects/ncount_corrected.RDS")
saveRDS(txi$counts, "./objects/raw_counts.RDS")

saveRDS(mdata, "./objects/expression.meta.RDS")

DEtx_results =list()
DEtx_results$PC3_T = performDEtx(mdata[mdata$cellLine == "PC3" & mdata$condition %in% c("DMSO","T"),],c("DMSO","T"))
DEtx_results$PC3_DKD = performDEtx(mdata[mdata$cellLine == "PC3" & mdata$condition %in% c("PLKO", "DKD"),],c("PLKO", "DKD"))
#DEtx_results$H660_SH1 = performDEtx(mdata[mdata$cellLine == "H660" & mdata$condition %in% c("WT", "SH1"),],c("WT", "SH1"))
#DEtx_results$H660_SH3 = performDEtx(mdata[mdata$cellLine == "H660" & mdata$condition %in% c("WT", "SH3"),],c("WT", "SH3"))
#DEtx_results$LnCaP = performDEtx(mdata[mdata$cellLine == "LnCaP" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
#DEtx_results$LnCaP_E = performDEtx(mdata[mdata$cellLine == "LnCaP_E" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
#DEtx_results$LnCaP_R = performDEtx(mdata[mdata$cellLine == "LnCaP_RES" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
#tmp = mdata[(mdata$cellLine == "LnCaP_RES" | mdata$cellLine == "LnCaP" ) & mdata$condition == "OE",]
#tmp$condition = tmp$cellLine
#DEtx_results$LnCaP_R_OE = performDEtx(tmp, c("LnCaP", "LnCaP_RES"))

saveRDS(DEtx_results, "./objects/DEtx_results.RDS")


DE_results =list()
DE_results$PC3_T = performDE(mdata[mdata$cellLine == "PC3" & mdata$condition %in% c("DMSO","T"),],c("DMSO","T"))
DE_results$PC3_DKD = performDE(mdata[mdata$cellLine == "PC3" & mdata$condition %in% c("PLKO", "DKD"),],c("PLKO", "DKD"))
# DE_results$H660_SH1 = performDE(mdata[mdata$cellLine == "H660" & mdata$condition %in% c("WT", "SH1"),],c("WT", "SH1"))
# DE_results$H660_SH3 = performDE(mdata[mdata$cellLine == "H660" & mdata$condition %in% c("WT", "SH3"),],c("WT", "SH3"))
# DE_results$LnCaP = performDE(mdata[mdata$cellLine == "LnCaP" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
# DE_results$LnCaP_E = performDE(mdata[mdata$cellLine == "LnCaP_E" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
# DE_results$LnCaP_R = performDE(mdata[mdata$cellLine == "LnCaP_RES" & mdata$condition %in% c("OE","WT"),], c("WT", "OE"))
# tmp = mdata[(mdata$cellLine == "LnCaP_RES" | mdata$cellLine == "LnCaP" ) & mdata$condition == "OE",]
# tmp$condition = tmp$cellLine
# DE_results$LnCaP_R_OE = performDE(tmp, c("LnCaP", "LnCaP_RES"))
# tmp = mdata[(mdata$cellLine == "LnCaP_RES"  & mdata$condition == "OE") | (mdata$cellLine == "LnCaP" & mdata$condition == "WT"),]
# tmp
# DE_results$LnCaP_ROE_P = performDE(tmp, c("OE", "WT"))
# tmp = mdata[(mdata$cellLine == "LnCaP_RES" | mdata$cellLine == "LnCaP" ) & mdata$condition == "WT",]
# tmp$condition = tmp$cellLine
# DE_results$LnCaP_PR = performDE(tmp, c("LnCaP", "LnCaP_RES"))

saveRDS(DE_results, "./objects/DE_results.RDS")

### Save the rawcounts, ncounts and DE results as CSV files
write.csv(rawcount, paste0(base_dir, "rawcounts_PC3_DKD_T.csv"), row.names = T, quote = F)

ncount_df = as.data.frame(ncount)
ncount_df$gene_id = rownames(ncount_df)
write.csv(ncount_df, paste0(base_dir, "ncounts_", name, ".csv"), row.names = F, quote = F)

#Run Expression.R from Caludio_Common/hg38 and save LNCaP count data
base_dir= "./expression_analysis/"
write.csv(ncount_df, paste0(base_dir, "ncounts_LNCaP_F2_Enza_ARCC4.csv"), row.names = F, quote = F)

for ( name in names(DE_results)){
  tmp = as.data.frame(DE_results[[name]])
  tmp= tmp[!is.na(tmp$padj),]
  tmp= tmp[order(tmp$padj),]
  tmp$gene = rownames(tmp)
  write.csv(tmp, paste0(base_dir, "DE_results_", name, ".csv"), row.names = F, quote = F)
}


#GT
### Add Gene names to results table
gene_info = read.csv("./meta/gene_info.csv")

for ( name in names(DE_results)){
  tmp = as.data.frame(DE_results[[name]])
  tmp= tmp[!is.na(tmp$padj),]
  tmp= tmp[order(tmp$padj),]
  tmp$gene_id = rownames(tmp)
  tmp = merge(tmp,gene_info, by = "gene_id")
  write.csv(tmp, paste0(base_dir, "DE_results_genesymbol_", name, ".csv"), row.names = F, quote = F)
}


#Plot FOXA1 and FOXA2 expression, Heatmap and scatter plot
#FOXA1 FOXA2 Expression ====
gene_info = read.csv("./meta/gene_info.csv")
ncount_df = read.csv("./expression_analysis/ncounts_PC3_DKD.csv")
ncount_df = merge(ncount_df,gene_info, by = "gene_id")

pdf(paste0(base_dir, "FOXA1_FOXA2_DKD_T.pdf"))
ncounts_FOX = ncounts %>%
  select(-c(gene_id, gene_type)) %>%
  filter(gene_name == "FOXA1" | gene_name == "FOXA2") %>%
  tidyr::gather(Sample,N_Counts, PC3_PLKO_1:PC3_T_3)

ncounts_FOX$Treatment = rep(c("PC3_PLKO","PC3_DKD","PC3_DMSO","PC3_T"), each = 6)

p = ggplot(data = ncounts_FOX, aes(x=factor(Treatment, level=c('PC3_PLKO', 'PC3_DKD', 
                                                               'PC3_DMSO', 'PC3_T')), y= N_Counts)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) + facet_wrap(~gene_name) +
  ggtitle("FOXA1 and FOXA2 expression") +
  xlab("") + ylab("Normalised Counts") + theme_classic()
p

dev.off()


#Plots
PC3_T =  read.csv("./expression_analysis/DE_results_genesymbol_PC3_T.csv")
PC3_KDK = read.csv("./expression_analysis/DE_results_genesymbol_PC3_DKD.csv")

library(dplyr)
library(plotly)
library(htmlwidgets)
library(hrbrthemes)
library(ggpubr)
library(survival)

#Density Plot ====
PC3_T_Sel = PC3_T %>%
  filter (gene_type == "protein_coding" | gene_type == "lncRNA" | gene_type == "miRNA") %>%
  select(gene_name,log2FoldChange,padj) %>% 
  mutate(Sample = "PC3_T", Sig = case_when(
    padj < 0.05 ~ "Sig",
    padj >= 0.05 ~ "Non-Sig"))
#names(PC3_T_Sel) = c("gene_name","log2FoldChange_T","padj_T","Sample")

PC3_KDK_Sel = PC3_KDK %>%
  filter (gene_type == "protein_coding" | gene_type == "lncRNA" | gene_type == "miRNA") %>%
  select(gene_name,log2FoldChange,padj) %>% 
  mutate(Sample = "PC3_KDK", Sig = case_when(
    padj < 0.05 ~ "Sig",
    padj >= 0.05 ~ "Non-Sig"))
#names(PC3_KDK_Sel) = c("gene_name","log2FoldChange_KDK","padj_KDK","Sample")

PC3_T_KDK = rbind(PC3_T_Sel,PC3_KDK_Sel)

base_dir= "./expression_analysis/"
pdf(paste0(base_dir, "Log2FCDensity_DKD_T.pdf"))
p1 <- ggplot(data=PC3_T_KDK, aes(x=log2FoldChange, group=Sample, fill=Sample)) +
  geom_density(adjust=1.5, alpha=.4) +theme_classic() + theme(legend.title = element_blank())
p1

p2 <- ggplot(data = filter(PC3_T_KDK, Sig == "Sig"), aes(x=log2FoldChange, group=Sample, fill=Sample)) +
  geom_density(adjust=1.5, alpha=.4) +
  ggtitle("padj < 0.05") +theme_classic() + theme(legend.title = element_blank())
p2
dev.off()

#Scatter Plot ====
Intersect_T_KDK = PC3_T_Sel %>%
  inner_join(PC3_KDK_Sel, by = "gene_name", keep = TRUE, suffix = c("_T", "_DKD"),)
#All_T_KDK = PC3_T_Sel %>%
#  full_join(PC3_KDK_Sel, by = "gene_name", keep = TRUE, suffix = c("_T", "_DKD"),)


# All_T_KDK[is.na(All_T_KDK$gene_name_T),]$log2FoldChange_T = 0
# All_T_KDK[is.na(All_T_KDK$gene_name_T),]$padj_T = 1
# All_T_KDK[is.na(All_T_KDK$gene_name_T),]$Sample_T = "PC3_T"
# All_T_KDK[is.na(All_T_KDK$gene_name_T),]$Sig_T = "Non-Sig"
# All_T_KDK[is.na(All_T_KDK$gene_name_T),]$gene_name_T = All_T_KDK[is.na(All_T_KDK$gene_name_T),]$gene_name_DKD
# 
# All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$log2FoldChange_DKD = 0
# All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$padj_DKD = 1
# All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$Sample_DKD = "PC3_DKD"
# All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$Sig_DKD = "Non-Sig"
# All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$gene_name_DKD = All_T_KDK[is.na(All_T_KDK$gene_name_DKD),]$gene_name_T

All_T_KDK = Intersect_T_KDK
pdf(paste0(base_dir, "Log2FC_DKD_T_Contourback_New.pdf"))
p = All_T_KDK %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD )) + 
  #geom_point(alpha = 1/4) +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()

p + geom_point(alpha = 1/15) + stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
                                        method = "pearson", label.x = -8)

p + geom_point(alpha = 1/15) + stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
                                        method = "pearson", label.x = -8) + 
  geom_density_2d(bins = 30, aes(color = ..level..), show.legend = FALSE) 

p + geom_point(alpha = 1/3) + stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
                                       method = "pearson", label.x = -6) + 
  geom_density_2d_filled(bins = 9,show.legend = FALSE, alpha = 0.7)  +
  scale_fill_brewer(palette = "Greys") + geom_density_2d(colour = "black")+
  lims(x= c(-6, 6), y = c(-6, 6))


p + geom_density_2d(bins = 30, aes(color = ..level..), show.legend = FALSE)  
p + geom_density_2d_filled(bins = 9,show.legend = FALSE) + theme(legend.position = "none") +
  scale_fill_brewer(palette = "Greys") + lims(x= c(-3, 3), y = c(-3, 3))


p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig", "red4", "grey")) %>%
ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD, color = Color )) + 
  # geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + geom_point(alpha = 1/3) +
  scale_color_identity() +  
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
             method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)") 

p + geom_density_2d(bins = 50) + geom_point(alpha = 1/3) +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)") 


p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig", "red4", "grey")) %>%
  filter(Color == "red4") %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD )) + 
  # geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + geom_density_2d(bins = 50, aes(color = ..level..)) + theme(legend.position = "none") 
p + geom_density_2d_filled(bins = 9) + theme(legend.position = "none") +
  scale_fill_brewer(palette = "Blues") + lims(x= c(-3, 3), y = c(-3, 3))
p + geom_density_2d(bins = 30, colour = "black") + geom_point(alpha = 1/3, color = "grey40") +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)")  


p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig" &
                          abs(log2FoldChange_T) >= 1 & abs(log2FoldChange_DKD) >= 1, "red4", "grey")) %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD, color = Color )) + 
  #geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + geom_point(alpha = 1/4) + scale_color_identity() +
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05& abs[log2FC] > 1 (T and DKD)")
p + geom_density_2d(bins = 50) + geom_point(alpha = 1/3) +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)") + geom_density_2d(bins = 50) 

p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig" &
                          abs(log2FoldChange_T) >= 1 & abs(log2FoldChange_DKD) >= 1, "red4", "grey")) %>%
  filter(Color == "red4") %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD)) + 
  #geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + geom_density_2d(aes(color = ..level..), bins = 50) + theme(legend.position = "none") 
p + geom_density_2d_filled(bins = 9) + theme(legend.position = "none") +
  scale_fill_brewer(palette = "Blues") + lims(x= c(-6, 5), y = c(-6, 5))
p  + geom_density_2d(bins = 30, colour = "black") + geom_point(alpha = 1/3, color = "grey40") +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)")


p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig" &
                          abs(log2FoldChange_T) >= 2 & abs(log2FoldChange_DKD) >= 2, "red4", "grey")) %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD, color = Color)) + 
  geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
             method = "pearson", label.x = -8)  +
  ggtitle("padj < 0.05 & abs[log2FC] > 2 (T and DKD)")
p + geom_density_2d(bins = 50) + geom_point(alpha = 1/3) +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)") 


p = All_T_KDK %>%
  mutate(Color = ifelse(Sig_T == "Sig" & Sig_DKD == "Sig" &
                          abs(log2FoldChange_T) >= 2 & abs(log2FoldChange_DKD) >= 2, "red4", "grey")) %>%
  filter(Color == "red4") %>%
  ggplot(aes(x = log2FoldChange_T, y = log2FoldChange_DKD )) + 
  #geom_point(alpha = 1/4) + scale_color_identity() +
  xlab("T - log2FC") + ylab("DKD - log2FC") + theme_classic()
p + geom_density_2d(aes(color = ..level..), bins = 50) + 
  theme(legend.position = "none")  + lims(x= c(-10, 10), y = c(-10, 10))
p + geom_density_2d_filled(bins = 9) + theme(legend.position = "none") +
  scale_fill_brewer(palette = "Blues") + lims(x= c(-10, 10), y = c(-6, 6))
p + geom_density_2d(bins = 30, colour = "black") + geom_point(alpha = 1/3, color = "grey40") +
  scale_color_identity() + 
  stat_cor(aes(x = log2FoldChange_T, y = log2FoldChange_DKD),
           method = "pearson", label.x = -8) +
  ggtitle("padj < 0.05 (T and DKD)") 

dev.off()


  
#, label=gene_name_T + geom_text(aes(label=ifelse(log2FoldChange_T < -2.5 & log2FoldChange_DKD < -2.5 ,
#as.character(gene_name_T),'')),hjust=1,vjust=1,size = 7.5/.pt)+

#ggsave("./expression_analysis/log2Foldchange_T_KDK.pdf", p, width =12, height =12)


##Just tried concordance - Not that useful
dd = All_T_KDK %>%
  filter(Sig_T == "Sig" & Sig_DKD == "Sig"  ) %>%
  mutate(Aver = rowMeans(select(.,log2FoldChange_T,log2FoldChange_DKD ), na.rm = TRUE),
         Diff = log2FoldChange_T - log2FoldChange_DKD )

# x1 and y2 are both continuous variables         
c = concordance(log2FoldChange_T ~ log2FoldChange_DKD, data= dd)
c

ggplot(dd, aes(x=Aver, y=Diff)) + 
  geom_point()

