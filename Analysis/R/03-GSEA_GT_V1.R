#First time installation ====
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("AnnotationDbi")
BiocManager::install("qusage")
BiocManager::install("EGSEA")
BiocManager::install("GeneOverlap")
BiocManager::install("clusterProfiler")

#Load libraries ====
neededlibraries <- c("edgeR", "org.Hs.eg.db", "ape", "ggplot2", "reshape","RColorBrewer","plotly","pheatmap","DESeq2",
                     "AnnotationDbi","org.Hs.eg.db","ggrepel","qusage","EGSEA","GeneOverlap","reshape2", "clipr", "svglite", 
                     "readtext","ggpubr", "clusterProfiler","tidyr")
lapply(neededlibraries, require, character.only = TRUE)

#Load data ====
#From Claudio setwd - I:\IOR-Groups\GL_Theurillat\Common\Common_Claudio\FOXA2_final\hg38
ncount = readRDS("./objects/ncount_corrected.RDS")
mdata = readRDS("./objects/expression.meta.RDS")
DE_results = readRDS("./objects/DE_results.RDS")
rcount = readRDS("./objects/raw_counts.RDS")
g_info = read.csv("./meta/gene_info.csv")

ncount_df = as.data.frame(ncount)
names(ncount_df)
salmon_merged_gene_counts = ncount_df %>%
  select("LnCaP_WT_1","LnCaP_WT_2","LnCaP_WT_3","LnCaP_E_WT_1",
  "LnCaP_E_WT_2","LnCaP_E_WT_3","LnCaP_RES_WT_1","LnCaP_RES_WT_2","LnCaP_RES_WT_3")

salmon_merged_gene_counts = ncount_df

salmon_merged_gene_counts$gene_id = rownames(salmon_merged_gene_counts)
salmon_merged_gene_counts = merge(salmon_merged_gene_counts,g_info, by = "gene_id")

counts = salmon_merged_gene_counts[,2:10] %>%
  mutate(across(is.numeric, round))

reduced_counts = salmon_merged_gene_counts[ rowSums(counts) > 0, ] 
SYMBOL = reduced_counts$gene_name

entrez <- mapIds(org.Hs.eg.db,
                 keys = SYMBOL,
                 column = "ENTREZID",
                 keytype = "SYMBOL",
                 multiVals = "first")
gene_annot = data.frame(SYMBOL=SYMBOL, ENTREZ=entrez)

categorys <- c("HALLMARK_ANDROGEN_RESPONSE","HALLMARK_E2F_TARGETS","HALLMARK_ESTROGEN_RESPONSE_EARLY",
               "HALLMARK_G2M_CHECKPOINT","HALLMARK_HYPOXIA","HALLMARK_MITOTIC_SPINDLE","HALLMARK_MTORC1_SIGNALING",
               "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2","HALLMARK_MYOGENESIS")

cat <- c("ANDROGEN_RESPONSE","E2F_TARGETS","ESTROGEN_RESPONSE_EARLY",
               "G2M_CHECKPOINT","HYPOXIA","MITOTIC_SPINDLE","MTORC1_SIGNALING",
               "MYC_TARGETS_V1","MYC_TARGETS_V2","MYOGENESIS")


# GSEA LNCaP ====
LnCaP = as.data.frame(DE_results$LnCaP)
LnCaP$gene_id = rownames(LnCaP)

LnCaP = merge(LnCaP,g_info,by = "gene_id")
LnCaP = merge(LnCaP,gene_annot,by.x = "gene_name", by.y = "SYMBOL")
LnCaP = LnCaP %>% drop_na(ENTREZ)


geneList<- LnCaP$stat
names(geneList)<-as.character(LnCaP$ENTREZ)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
head(geneList)

library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
CRVSR <- GSEA(geneList,
              TERM2GENE = m_t2g,
              pvalueCutoff = 0.05,
              eps=1e-50,
              verbose = FALSE, seed = T)

res_LNCaP = CRVSR@result
write.csv(res_LNCaP,'I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/expression_analysis/GSEA/LNCaP_CTR_Hallmark.csv', 
          row.names = FALSE)

dotplot(CRVSR, showCategory=categorys ,  
        orderBy="p.adjust", x = "enrichmentScore") 


# GSEA LNCaP_E ====
LnCaP_E = as.data.frame(DE_results$LnCaP_E)
LnCaP_E$gene_id = rownames(LnCaP_E)

LnCaP_E = merge(LnCaP_E,g_info,by = "gene_id")
LnCaP_E = merge(LnCaP_E,gene_annot,by.x = "gene_name", by.y = "SYMBOL")
LnCaP_E = LnCaP_E %>% drop_na(ENTREZ)


geneList<- LnCaP_E$stat
names(geneList)<-as.character(LnCaP_E$ENTREZ)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
head(geneList)

library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
CR1VSR <- GSEA(geneList,
              TERM2GENE = m_t2g,
              pvalueCutoff = 0.05,
              eps=1e-50,
              verbose = FALSE, seed = T)

res_LNCaP_E = CR1VSR@result
write.csv(res_LNCaP_E,'I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/expression_analysis/GSEA/LNCaP_E_Hallmark.csv', 
          row.names = FALSE)

dotplot(CRVSR, showCategory=categorys ,  
        orderBy="p.adjust", x = "enrichmentScore") 


# GSEA LNCaP_R ====
LnCaP_R = as.data.frame(DE_results$LnCaP_R)
LnCaP_R$gene_id = rownames(LnCaP_R)

LnCaP_R = merge(LnCaP_R,g_info,by = "gene_id")
LnCaP_R = merge(LnCaP_R,gene_annot,by.x = "gene_name", by.y = "SYMBOL")
LnCaP_R = LnCaP_R %>% drop_na(ENTREZ)


geneList<- LnCaP_R$stat
names(geneList)<-as.character(LnCaP_R$ENTREZ)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
head(geneList)

library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
CR1VSCR <- GSEA(geneList,
              TERM2GENE = m_t2g,
              pvalueCutoff = 0.05,
              eps=1e-50,
              verbose = FALSE, seed = T)

res_LNCaP_R = CR1VSCR@result
write.csv(res_LNCaP_R,'I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/expression_analysis/GSEA/LNCaP_R_Hallmark.csv', 
          row.names = FALSE)

dotplot(CRVSR, showCategory=categorys ,  
        orderBy="p.adjust", x = "enrichmentScore") 


### VERY IMPORTANT - Merge plots ====
library(DOSE)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)


x.1 = CRVSR
x.2 = CR1VSR
x.3 = CR1VSCR

## count the gene number for both results
gene_count.x1 <- x.1@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x2 <- x.2@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
gene_count.x3 <- x.3@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)


## merge with the original dataframes
dot_df.x1<- left_join(x.1@result, gene_count.x1, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x2<- left_join(x.2@result, gene_count.x2, by = "ID") %>% mutate(GeneRatio = count/setSize)
dot_df.x3<- left_join(x.3@result, gene_count.x3, by = "ID") %>% mutate(GeneRatio = count/setSize)

## merge the two results
library(clusterProfiler)
merged.res <- as.data.frame(merge_result(list(DMSO =dot_df.x1, ENZA=dot_df.x2, ARCC4 = dot_df.x3)))

## merged.res <- rbind(dot_df.x1, dot_df.x2) #This merging works but it does **not** include source of results (i.e. 'fgsea' or 'dose')

## Set up/downregulation
merged.res$type = "upregulated"
merged.res$type[merged.res$NES < 0] = "downregulated"

merged.res$SIG = -log10(merged.res$p.adj)
p <- ggplot(merged.res , aes(x = enrichmentScore, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = enrichmentScore)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2(limits=c(-0.8, 0.8), low = "darkblue",high = "#b50404", mid = "white") +
  ylab(NULL) + facet_grid(.~Cluster) +
  scale_size(range = c(1,10)) 
p
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_LNCaP_E_R_ALL.pdf', p, width = 12, height =8)
#readtext("plot.svg")-> temp
#write_clip(temp$text, object_type = "character")

p_sel <- ggplot(subset(merged.res,Description %in%categorys) , aes(x = enrichmentScore, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = enrichmentScore)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2(limits=c(-0.8, 0.8), low = "darkblue",high = "#b50404", mid = "white") +
  ylab(NULL) + facet_grid(.~Cluster) +
  scale_size(range = c(1,10)) 
p_sel
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_LNCaP_E_R_Sel.pdf', p_sel, width = 12, height =8)


p <- ggplot(merged.res, aes(x = NES, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = NES)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2(limits=c(-3.5, 3.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + facet_grid(.~Cluster) +
  scale_size(range = c(3,12)) 
p
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_NES_LNCaP_E_R_ALL.pdf', p, width =18, height = 10)
#readtext("plot.svg")-> temp
#write_clip(temp$text, object_type = "character")

p_sel <- ggplot(subset(merged.res,Description %in%categorys), aes(x = NES, y = fct_reorder(Description, SIG))) + 
  geom_point(aes(size = SIG, color = NES)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank()) +
  scale_colour_gradient2(limits=c(-3.5, 4.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + facet_grid(.~Cluster,scales='free_x', space='free_x') +
  scale_size(range = c(3,12)) 
  
p_sel
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_NES_LNCaP_E_R_Sel_New.pdf', p_sel, width =12, height = 4)


merged.res = merged.res %>% 
  separate_wider_delim(Description, "_", names = c("Hallmark", "Desc"),
                       too_many = "merge", too_few = "align_start")


p_sel <- ggplot(subset(merged.res,Desc %in%cat), aes(x = Cluster, y = fct_reorder(Desc, SIG))) + 
  geom_point(aes(size = SIG, color = NES)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 12),axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5)) +
  scale_colour_gradient2(limits=c(-3.5, 4.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + 
  scale_size(range = c(3,12)) + labs(size ='-log10(adj.p)', title = 'LNCaP-FOXA2 vs LNCaP-CTR') 

p_sel
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_NES_LNCaP_E_R_Sel_Merged.pdf', p_sel, width =5, height = 5)

#Plot PCa Profiler Result ====

F2_high_low = read.table('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/resultsF2HIGHVSLOWARnegativePCa.tsv',
                         sep = '\t', stringsAsFactors = FALSE,header = TRUE)

F2_high_low_ENTREZ = merge(F2_high_low,gene_annot,by.x= "Gene.Name",by.y= "SYMBOL")
geneList<- F2_high_low_ENTREZ$statistic

names(geneList)<-as.character(F2_high_low_ENTREZ$ENTREZ)
geneList<-geneList[!duplicated(names(geneList))]
geneList<-sort(geneList,decreasing=TRUE)
head(geneList)

library(msigdbr)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)
PCa_Prof <- GSEA(geneList,
               TERM2GENE = m_t2g,
               pvalueCutoff = 1, 
               eps=1e-50,
               verbose = FALSE, seed = T)

res_PCa_Prof = PCa_Prof@result
write.csv(res_PCa_Prof,'I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/expression_analysis/GSEA/PCaProf_Hallmark.csv', 
          row.names = FALSE)


library(DOSE)
library(ggplot2)
library(dplyr)
library(stringr)
library(forcats)


x.1 = PCa_Prof
# x.2 = CR1VSR
# x.3 = CR1VSCR

## count the gene number for both results
gene_count.x1 <- x.1@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
# gene_count.x2 <- x.2@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
# gene_count.x3 <- x.3@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)


## merge with the original dataframes
dot_df.x1<- left_join(x.1@result, gene_count.x1, by = "ID") %>% mutate(GeneRatio = count/setSize)
# dot_df.x2<- left_join(x.2@result, gene_count.x2, by = "ID") %>% mutate(GeneRatio = count/setSize)
# dot_df.x3<- left_join(x.3@result, gene_count.x3, by = "ID") %>% mutate(GeneRatio = count/setSize)

## merge the two results
# library(clusterProfiler)
# merged.res <- as.data.frame(merge_result(list(DMSO =dot_df.x1, Enza=dot_df.x2, ARCC4 = dot_df.x3)))

## merged.res <- rbind(dot_df.x1, dot_df.x2) #This merging works but it does **not** include source of results (i.e. 'fgsea' or 'dose')

## Set up/downregulation
dot_df.x1$type = "upregulated"
dot_df.x1$type[dot_df.x1$NES < 0] = "downregulated"

dot_df.x1$SIG = -log10(dot_df.x1$p.adj)

#Selection1
categorys <- c("HALLMARK_E2F_TARGETS",
               "HALLMARK_G2M_CHECKPOINT",
               "HALLMARK_MITOTIC_SPINDLE",
               "HALLMARK_MYC_TARGETS_V1",
               "HALLMARK_MYC_TARGETS_V2")
#subset(dot_df.x1,Description %in%categorys)
categorys = c("HALLMARK_E2F_TARGETS","HALLMARK_G2M_CHECKPOINT",
              "HALLMARK_MYC_TARGETS_V1","HALLMARK_MYC_TARGETS_V2")
dot_df.x1_sel = subset(dot_df.x1,Description %in% categorys)

#Selection2
top10 = dot_df.x1 %>% slice_max(NES, n = 10)
bottom10 = dot_df.x1 %>% slice_min(NES, n = 10)
dot_df.x1_sel = rbind(top10,bottom10)

dot_df.x1_sel = dot_df.x1_sel %>% 
            separate_wider_delim(Description, "_", names = c("Hallmark", "Desc"),
                                 too_many = "merge", too_few = "align_start")




p_sel <- ggplot(dot_df.x1_sel, aes(x = NES, y = fct_reorder(Desc, SIG))) + 
  geom_point(aes(size = SIG, color = NES)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 10),axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), axis.text.y = element_text(size = 12), plot.title = element_text(hjust = 0.5)) +
  scale_colour_gradient2(limits=c(0, 4.5), low = "#3953A4", high = "#b50404") +
  ylab(NULL) + scale_x_continuous(minor_breaks = seq(-4, 6, 2),limits = c(1,4)) +
  scale_size(range = c(3,12)) + labs(size ='-log10(adj.p)', title = expression(paste('FOXA2'^high,' vs FOXA2'^low, '(AR Negative CRPCs)')))

p_sel
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_NES_F2High_low_Sel.pdf', p_sel, width =4, height = 3)
ggsave('I:/IOR-Groups/GL_Theurillat/Gayathri Thillaiyampalam/FOXA2/RNASeq_DE/png/GSEA_NES_F2High_low_Sel_legend.pdf', p_sel, width =4, height = 4)

