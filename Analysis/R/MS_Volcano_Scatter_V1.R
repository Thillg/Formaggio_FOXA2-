#!/usr/bin/Rscript
library("ggplot2")
#BiocManager::install("tximport")
#BiocManager::install("SparseArray")
library(tximport)
library(clusterProfiler)
library("RColorBrewer")
library(dplyr)
library(ggpubr)
library(matrixStats)


PC3_DKD =  read.csv("./PC3_MS_DKD_PLKO.csv")
PC3_F1 = read.csv("./PC3_MS_F1_PLKO.csv")
PC3_F2 = read.csv("./PC3_MS_F2_PLKO.csv")

PC3_DKD_F1 = read.csv("./PC3_MS_DKD_F1.csv")
PC3_DKD_F2 = read.csv("./PC3_MS_DKD_F2.csv")

TF_info = read.csv("./TF.csv")

PC3_MS_T = read.csv("./PC3_MS.csv")
LNCaP_MS_T = read.csv("./LNCaP_MS.csv")
H660_MS_T = read.csv("./H660_MS.csv")

PC3_BioID = read.csv("./F1andF2_BioID_PC3.csv")
LNCaP_BioID = read.csv("./F1andF2_BioID_LNCaP_Filtered.csv")
H660_BioID = read.csv("./F1andF2_BioID_H660.csv")

PC3_MS = list(PC3_DKD,PC3_F1,PC3_F2,PC3_DKD_F1,PC3_DKD_F2,PC3_MS_T)
names(PC3_MS) = c("PC3_DKD_PLKO","PC3_F1_PLKO","PC3_F2_PLKO","PC3_DKD_F1","PC3_DKD_F2",
                  "PC3_MS_T")

MS_T = list(PC3_MS_T,LNCaP_MS_T,H660_MS_T)
names(MS_T) = c("PC3_MS_T","LNCaP_MS_T","H660_MS_T")

#VolcanoPlots ====
library(EnhancedVolcano)

for (n in names(PC3_MS))
{
  PC3_MS[[n]] = PC3_MS[[n]] %>%
    mutate(TF = case_when(
      gene_name %in% TF_info$TranscriptionFactors ~ "Yes" ,
      TRUE ~ "No"),
      absDiff = abs(Difference)) %>%
    group_by(gene_name) %>%
    filter(absDiff == max(absDiff), na.rm=TRUE)
}

for (n in names(MS_T))
{
  MS_T[[n]] = MS_T[[n]] %>%
    mutate(TF = case_when(
      gene_name %in% TF_info$TranscriptionFactors ~ "Yes" ,
      TRUE ~ "No"),
      absDiff = abs(Difference)) %>%
    group_by(gene_name) %>%
    filter(absDiff == max(absDiff), na.rm=TRUE)
}

TF_Sel = c("FOXA1","FOXA2","JUNB","ETV5","ETV1","HOXA13",
           "HOXB13","FOSL2","ID1")

Overlap = c("PRR12", "AHDC1", "CUX1", "ZNF384", "MYBL2", 
            "TCF20", "JMJD1C", "ARID1A", "ARID1B","PIAS1", "PIAS2", "PIAS3",
            "BCOR", "NCOR1", "NCOR2", "PHF21A", "TLE3", "ELMSAN1")


#New Volcanoplot with labels on one side ====
#https://biostatsquid.com/volcano-plots-r-tutorial/
#https://ggrepel.slowkow.com/articles/examples.html
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

plotVolcano_BioID_label <- function(BioID,MS,n,title,subtitle)
{
  BioID_F1 = BioID %>%
    slice_max(F1, n = n)
  
  BioID_F2 = BioID %>%
    slice_max(F2, n = n)
  
  BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names, "FOXA1", "FOXA2", "NKX2-1")
  BioID_F1_F2 = unique(BioID_F1_F2)
  
  #uncomment for list
  # final_plot = list()
  #  for (n in names(MS))
  #  {
  #  res = MS[[n]]
  
  res = MS
  
  TF_Sig = res %>%
    filter(TF == "Yes" & abs(Difference) >=1 & q_Value < 0.05) 
  
  BioID_F1_F2_Sig = res %>%
    filter(gene_name %in% BioID_F1_F2 & abs(Difference) >=0.5 & q_Value < 0.05) 
  
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
  res$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  res$diffexpressed[res$Difference > 0.5 & res$q_Value < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  res$diffexpressed[res$Difference < -0.5 & res$q_Value < 0.05] <- "DOWN"
  
  p_min = min(res[which(res$q_Value > 0),'q_Value'])
  res[which(res[,'q_Value'] == 0), 'q_Value'] <- p_min
  
  f_down = res %>%
    filter(gene_name %in% BioID_F1_F2 & Difference <=-0.5 & q_Value < 0.05) 
  FOXA1_FOXA2 = res %>%
    filter(gene_name == "FOXA1" | gene_name == "FOXA2" )
  FOXA1 = res %>%
    filter(gene_name == "FOXA1")
  f_up = res %>%
    filter(gene_name %in% BioID_F1_F2 & Difference >=0.5 & q_Value < 0.05)
  
  f_down_F1_F2 <- rbind(f_down, FOXA1_FOXA2)  
  f_down_F1 <- rbind(f_down, FOXA1)  
  
  # Add threshold lines
  final_plot = ggplot(data = res, aes(x = Difference, y = -log10(q_Value))) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(alpha = 4/5, size = 2.5, shape=19,
               aes(col = diffexpressed)) + # Biostatsquid theme
    theme_set(theme_classic(base_size = 14) +
                theme(
                  axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                  axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                  plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=14),
                  axis.text.y = element_text(size=14),
                )) + 
    scale_color_manual(values = c("royalblue", "grey40", "red"),
                       labels = c("Downregulated", "Not significant", "Upregulated")) +
    
    labs(color = 'diffexpressed', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    ggtitle(paste0(title," ",subtitle)) + # Plot title 
    # geom_label_repel(data=res %>%
    #                    filter(gene_name %in% BioID_F1_F2 & abs(Difference) >=0.5 & q_Value < 0.05), # Filter data first
    #                  aes(label=gene_name),color = "grey20", #max.overlaps = 5000, 
    #                  size= 5,
    #                  box.padding = 0.8, max.overlaps = Inf, direction = "x")

    geom_label_repel(data=f_up, # Filter data first
                      aes(label= gene_name),color = "grey20", #max.overlaps = 5000, 
                      size= 5, box.padding = 0.5, max.overlaps = Inf, force = 0.1,
                      nudge_x  = 4 - f_up$Difference , 
                      direction  = "y", hjust = 0) +  
    geom_label_repel(data=f_down, # Filter data first
                     aes(label=gene_name),color = "grey20", #max.overlaps = 5000, 
                     size= 5, box.padding = 0.5, max.overlaps = Inf, force = 0.1,
                     #nudge_x      = -1,
                     nudge_x  = -5 - f_down$Difference , 
                     direction    = "y" ,
                     hjust = 1, segment.size = 0.2) +
    ylim (0,5)  + xlim(-6,6) + guides(color="none") 

  
  return(final_plot)
  
}

PC3_T_BioID_50_lab = plotVolcano_BioID_label(PC3_BioID,MS_T$PC3_MS_T,50,"PC3_MS_T_BioID", "Top50") #xlim = c(-6,6), ylim = c(0,5), data=f_down
LNCaP_T_BioID_50_lab = plotVolcano_BioID_label(LNCaP_BioID,MS_T$LNCaP_MS_T,50,"LNCaP_MS_T_BioID_NoMito", "Top50") #xlim = c(-6,6), ylim = c(0,5), data=f_down
H660_T_BioID_50_lab = plotVolcano_BioID_label(H660_BioID,MS_T$H660_MS_T,50,"H660_MS_T_BioID", "Top50") #xlim = c(-6,6), ylim = c(0,4), data=f_down_F1_F2

PC3_MS_BioID_50_lab_grey = plotVolcano_BioID_label(PC3_BioID,PC3_MS$PC3_DKD_PLKO,50,"PC3_MS_DKD_BioID","Top50") #ylim (0,4)  + xlim(-6,6),  data=f_down


pdf("./PC3_LNCaP_H660_MS_T_VolcanoPlot_p05_FC05_BioID50_label_New4.1.pdf", 
    width = 10, height = 10 )
PC3_T_BioID_50_lab
LNCaP_T_BioID_50_lab
H660_T_BioID_50_lab
dev.off()

pdf("./H660_MS_T_EnhancedVolcanoPlot_p05_FC05_BioID50_label_New2.1.pdf", 
    width = 10, height = 10 )
H660_T_BioID_50_lab
dev.off()

pdf("./PC3_MS_DKD_VolcanoPlot_p05_FC05_BioID50_label_New4.1.pdf", 
    width = 10, height = 10 )
PC3_MS_BioID_50_lab_grey
dev.off()

plotVolcano_BioID <- function(BioID,MS,n,title,subtitle)
{
  BioID_F1 = BioID %>%
    slice_max(F1, n = n)
  
  BioID_F2 = BioID %>%
    slice_max(F2, n = n)
  
  BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names)
  BioID_F1_F2 = unique(BioID_F1_F2)
  
  res = MS
  
  TF_Sig = res %>%
    filter(TF == "Yes" & abs(Difference) >=1 & q_Value < 0.05) 
  
  BioID_F1_F2_Sig = res %>%
    filter(gene_name %in% BioID_F1_F2 & abs(Difference) >=1 & q_Value < 0.05) 
  res = res %>%
    mutate(Regulation = case_when(
      q_Value < 0.05 & Difference > 0  ~ "Up" ,
      q_Value < 0.05 & Difference < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    BioID = case_when(
      gene_name %in% BioID_F1_F2  ~ "Yes" ,
      TRUE ~ "No"
    ) )
  
  t = data.frame(Y = c(sum(res$BioID == "Yes" & res$Regulation == "Up"),
                       sum(res$BioID == "Yes" & res$Regulation == "Down")),
                 N = c(sum(res$BioID == "No" & res$Regulation == "Up"),
                       sum(res$BioID == "No" & res$Regulation == "Down")))
  test <- chisq.test(t)
  p = test$p.value
  
  BioID_F1_F2_Filt = setdiff(BioID_F1_F2,Overlap)
  
  # # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
  # res$diffexpressed <- "NO"
  # # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  # res$diffexpressed[res$Difference > 0.5 & res$q_Value < 0.05] <- "UP"
  # # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  # res$diffexpressed[res$Difference < -0.5 & res$q_Value < 0.05] <- "DOWN"
  
  p_min = min(res[which(res$q_Value > 0),'q_Value'])
  res[which(res[,'q_Value'] == 0), 'q_Value'] <- p_min
  
  f_down = res %>%
    filter(gene_name %in% BioID_F1_F2 & Difference <=-0.5 & q_Value < 0.05) 
  FOXA1_FOXA2 = res %>%
    filter(gene_name == "FOXA1" | gene_name == "FOXA2" )
  FOXA1 = res %>%
    filter(gene_name == "FOXA1")
  f_up = res %>%
    filter(gene_name %in% BioID_F1_F2 & Difference >=0.5 & q_Value < 0.05)
  
  f_down_F1_F2 <- rbind(f_down, FOXA1_FOXA2)  
  f_down_F1 <- rbind(f_down, FOXA1)  
  
  res <- res %>%
    mutate(diffexpressed  = case_when(
      gene_name %in% BioID_F1_F2_Filt  ~ "BioID_LSpecific" ,
      gene_name %in% Overlap ~ "BioID_Core",
      TRUE ~ "Others"
    )) %>%
    arrange(desc(diffexpressed))
  # keyvals <- ifelse(
  #   res$gene_name %in% BioID_F1_F2_Filt, 'royalblue',
  #   ifelse(res$gene_name %in% Overlap, 'red4',
  #          'grey70'))
  # keyvals[is.na(keyvals)] <- 'black'
  # names(keyvals)[keyvals == 'grey70'] <- 'Others'
  # names(keyvals)[keyvals == 'royalblue'] <- 'BioID Lineage Specific'
  # names(keyvals)[keyvals == 'red4'] <- 'BioID Core proteins'
  # 
  # red <- which(keyvals == 'red4')
  # blue <- which(keyvals == 'royalblue')
  # r_b = c(blue,red)
  # notblue <- setdiff(1:length(keyvals), r_b)
  # res <- res[c(notblue, blue,red),]
  # keyvals <- keyvals[c(notblue, blue,red)]
  
  # Add threshold lines
  final_plot = ggplot(data = res, aes(x = Difference, y = -log10(q_Value))) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(alpha = 4/5, shape=19,
               aes(col = diffexpressed, size = factor(diffexpressed))) + # Biostatsquid theme
    theme_set(theme_classic(base_size = 14) +
                theme(
                  axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
                  axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
                  plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size=14),
                  axis.text.y = element_text(size=14)
                )) + 
    scale_color_manual(values = c("red4","royalblue", "grey70" ),
                       labels = c("BioID Core Proteins","BioID Lineage Specific", "Others")) +
    scale_size_manual(values = c(4, 4, 2)) + 
    
    labs(color = 'diffexpressed', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
    ggtitle(paste0(title," ",subtitle)) + # Plot title 
    ylim (0,3.5)  + xlim(-2,2) +
    annotate("text", x=-1.5, y=3.5, label=paste0("Chi-square\np-value: ", signif(p,3))) +
    guides(color="none", size="none")
    
  
  return(final_plot)
  
}

PC3_T_BioID_100 = plotVolcano_BioID(PC3_BioID,MS_T$PC3_MS_T,100,"PC3_MS_T_BioID", "Top100") #xlim = c(-6,6), ylim = (0,4.5), annotate("text", x=-4, y=4.5,..)
LNCaP_T_BioID_100 = plotVolcano_BioID(LNCaP_BioID,MS_T$LNCaP_MS_T,100,"LNCaP_MS_T_BioID_NoMito", "Top100") #xlim = c(-6,6), ylim = (0,4.5), annotate("text", x=-4, y=4.5,..)
H660_T_BioID_100 = plotVolcano_BioID(H660_BioID,MS_T$H660_MS_T,100,"H660_MS_T_BioID", "Top100") #xlim = c(-2,2),ylim = (0,3.5) annotate("text", x=-1.5, y=3.5,..)

PC3_MS_BioID_100_DKD_PLKO = plotVolcano_BioID(PC3_BioID,PC3_MS$PC3_DKD_PLKO,100,"PC3_MS_DKD_BioID","Top100") #ylim (0,3.5)  + xlim(-4,4), annotate("text", x=-2, y=3.5,..)
PC3_MS_BioID_100_F1 = plotVolcano_BioID(PC3_BioID,PC3_MS$PC3_F1_PLKO,100,"PC3_MS_F1_BioID","Top100") #ylim (0,3)  + xlim(-4,4), annotate("text", x=-2, y=3,..)
PC3_MS_BioID_100_F2 = plotVolcano_BioID(PC3_BioID,PC3_MS$PC3_F2_PLKO,100,"PC3_MS_F2_BioID","Top100") #ylim (0,3.5)  + xlim(-4,4), annotate("text", x=-2, y=3.5,..)

write.csv(PC3_T_BioID_100$data, "./PC3_MS_Plotdata.csv", quote=F, row.names = F)
write.csv(LNCaP_T_BioID_100$data, "./LNCaP_MS_Plotdata.csv", quote=F, row.names = F)
write.csv(H660_T_BioID_100$data, "./H660_MS_Plotdata.csv", quote=F, row.names = F)

pdf("./PC3_LNCaP_H660_MS_T_VolcanoPlot_p05_FC05_BioID50_Coloured_Core_Enrichment_New4.5.pdf",
    height = 8, width = 10)
PC3_T_BioID_100
LNCaP_T_BioID_100
H660_T_BioID_100
dev.off()

pdf("./PC3_MS_DKD_VolcanoPlot_p05_FC05_BioID50_Coloured_Core_Enrichment_New4.1.pdf",
    height = 8, width = 8)
PC3_MS_BioID_100_DKD_PLKO
PC3_MS_BioID_100_F1
PC3_MS_BioID_100_F2
dev.off()

#2D Scatter ====
#DKD_F1 vs DKD_F2
PC3_DKD_F = merge(PC3_MS$PC3_DKD_F1,PC3_MS$PC3_DKD_F2,by="gene_name", suffixes = c("_F1", "_F2"))

BioID_F1 = PC3_BioID %>%
  slice_max(F1, n = 100) %>%
  filter(Gene.names != "FOXA1" & Gene.names != "FOXA2" )

BioID_F2 = PC3_BioID %>%
  slice_max(F2, n = 100) %>%
  filter(Gene.names != "FOXA1" & Gene.names != "FOXA2")

BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names)
BioID_F1_F2 = unique(BioID_F1_F2)

library(ggrepel)
library(ggpp)

p = PC3_DKD_F %>%
  mutate(Color = ifelse(q_Value_F1 < 0.05 & q_Value_F2 < 0.05   , "red", "grey"),
         my_label = ifelse(gene_name %in% BioID_F1_F2 & q_Value_F1 < 0.05 & q_Value_F2 < 0.05,
                           gene_name, "")) %>%
  ggplot(aes(x=Difference_F1, y=Difference_F2, color = Color))+ geom_point(alpha = 4/5, size = 1.5) +
  xlab("DKD vs shFOXA1 - log2FoldChange") + ylab("DKD vs shFOXA2 - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_F1, y = Difference_F2),
           method = "pearson", label.x = -4)+ theme_classic() + xlim (-4,3) + ylim(-3,3) +
  # geom_label_repel(aes(label = my_label), box.padding = 0.1, 
  #                  max.overlaps = 5000, point.size = NA, size = 5, colour = "grey20", 
  #                  position = 
  #                    position_nudge_center(x = 0.5,
  #                                          y = -0.5,
  #                                          center_x = 0,
  #                                          center_y = 2,
  #                                          direction = "radial")) +
  geom_label_repel(aes(label = my_label),
                   size= 5,colour = "grey20", 
                   box.padding = 0.7, max.overlaps = 5000,
                   position = 
                     position_nudge_center(x = 0.5,
                                           y = -0.5,
                                           center_x = 0,
                                           center_y = 1,
                                           direction = "radial")) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
p


pdf("./DKD_F1_F2_ScatterPlot_BioID100_label_New_2.3.pdf", width = 12, height = 12)
p
dev.off()

p_trimmed = PC3_DKD_F %>%
  mutate(Color = ifelse(q_Value_F1 < 0.05 & q_Value_F2 < 0.05   , "red4", "grey"),
         my_label = ifelse(gene_name %in% BioID_F1_F2 & q_Value_F1 < 0.05 & q_Value_F2 < 0.05,
                           gene_name, "")) %>%
  ggplot(aes(x=Difference_F1, y=Difference_F2, color = Color))+ geom_point(alpha = 4/5, size = 2) +
  xlab("DKD vs shFOXA1 - log2FoldChange") + ylab("DKD vs shFOXA2 - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_F1, y = Difference_F2),
           method = "pearson", label.x = -4)+ theme_classic() + xlim (-4,2) + ylim(-4,4) +
  # geom_label_repel(aes(label = my_label), box.padding = 0.1, 
  #                  max.overlaps = 6000, point.size = NA, size = 5, colour = "black") +
  geom_label_repel(aes(label = my_label),
                   size= 6,
                   box.padding = 0.8, max.overlaps = Inf) +
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
p_trimmed


pdf("./DKD_F1_F2_ScatterPlot_BioID100_label_New_2.pdf")
#width = 10, height = 10)
p
#p_trimmed
dev.off()

ggsave(paste0("./DKD_F1_F2_ScatterPlot_BioID100_label_New_2.pdf"),p)

#, height = 10, width = 10)

#Highlight BioID
p = PC3_DKD_F %>%
  mutate(Color = ifelse(gene_name %in% BioID_F1_F2,
                        "red4", "grey70")) %>% arrange(Color) %>%
  ggplot(aes(x=Difference_F1, y=Difference_F2, color = Color))+ geom_point(alpha = 1/4, size = 1) +
  xlab("DKD vs shFOXA1 - log2FoldChange") + ylab("DKD vs shFOXA2 - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_F1, y = Difference_F2),
           method = "pearson", label.x = -6)+ theme_classic()
p

ggsave(paste0("./DKD_F1_F2_ScatterPlot_BioID100_Colour.pdf"),p, height = 8, width = 8)

#F1_PLKO vs F2_PLKO
PC3_F1_F2_PLKO = merge(PC3_F1,PC3_F2,by="gene_name", suffixes = c("_F1", "_F2"))

p = PC3_F1_F2_PLKO %>%
  mutate(Color = ifelse(q_Value_F1 < 0.05 & q_Value_F2 < 0.05   , "red4", "grey")) %>%
  ggplot(aes(x=Difference_F1, y=Difference_F2, color = Color))+ geom_point(alpha = 1/5, size = 1) +
  xlab("shFOXA1 vs PLKO - log2FoldChange") + ylab("shFOXA2 vs PLKO - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_F1, y = Difference_F2),
           method = "pearson", label.x = -8)+ theme_classic()
p

ggsave(paste0("./F1_F2_PLKO_ScatterPlot.pdf"),p, height = 8, width = 8)

#DKD vs T
PC3_DKD_T = merge(PC3_MS$PC3_DKD_PLKO,MS_T$PC3_MS_T,by.x ="gene_name",by.y = "gene_name", suffixes = c("_DKD", "_T"))

BioID_F1 = PC3_BioID %>%
  slice_max(F1, n = 100)

BioID_F2 = PC3_BioID %>%
  slice_max(F2, n = 100)

BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names)
BioID_F1_F2 = unique(BioID_F1_F2)

#Corr sig
p = PC3_DKD_T %>%
  mutate(Color = ifelse(q_Value_DKD < 0.05 & q_Value_T < 0.05,
                        "red4", "grey")) %>%
  ggplot(aes(x=Difference_DKD, y=Difference_T, color = Color))+ geom_point(alpha = 1/5, size = 1) +
  xlab("DKD - log2FoldChange") + ylab("T - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_DKD, y = Difference_T),
           method = "pearson", label.x = -8)+ theme_classic()
p

ggsave(paste0("./PC3_MS_DKD_T_ScatterPlot.pdf"),p, height = 8, width = 8)

p = PC3_DKD_T %>%
  mutate(Color = ifelse(q_Value_DKD < 0.05 & q_Value_T < 0.05,
                        "red4", "grey"),
         my_label = ifelse(gene_name %in% BioID_F1_F2 & q_Value_DKD < 0.05 & q_Value_T < 0.05,
                           gene_name, "")) %>%
  ggplot(aes(x=Difference_DKD, y=Difference_T, color = Color))+ geom_point(alpha = 1/5, size = 1) +
  xlab("DKD - log2FoldChange") + ylab("T - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_DKD, y = Difference_T),
           method = "pearson", label.x = -8)+ theme_classic() +
  geom_label_repel(aes(label = my_label), box.padding = 0.1, 
                   max.overlaps = 6000, point.size = NA, size = 2, colour = "black")
p

ggsave(paste0("./PC3_MS_DKD_T_ScatterPlot_BioID100_label.pdf"),p, height = 8, width = 8)

#Corr bioID
p = PC3_DKD_T %>%
  mutate(Color = ifelse(gene_name %in% BioID_F1_F2,
                        "red4", "grey70")) %>% arrange(Color) %>%
  ggplot(aes(x=Difference_DKD, y=Difference_T, color = Color))+ geom_point(alpha = 1/4, size = 1.5) +
  xlab("DKD - log2FoldChange") + ylab("T - log2FoldChnage") + scale_color_identity() +
  stat_cor(aes(x = Difference_DKD, y = Difference_T),
           method = "pearson", label.x = -8)+ theme_classic()
p

ggsave(paste0("./PC3_MS_DKD_T_ScatterPlot_BioID100_colour.pdf"),p, height = 8, width = 8)

#Enrichment analysis ====

plotEnrich <- function(BioID,MS,n)
{
  BioID_F1 = BioID %>%
    slice_max(F1, n = n)
  
  BioID_F2 = BioID %>%
    slice_max(F2, n = n)
  
  BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names)
  BioID_F1_F2 = unique(BioID_F1_F2)
  
  final_plot = list()
  res = MS
  res = res %>%
    mutate(Regulation = case_when(
      q_Value < 0.05 & Difference > 0  ~ "Up" ,
      q_Value < 0.05 & Difference < 0 ~ "Down",
      TRUE ~ "NS"
    ),
    BioID = case_when(
      gene_name %in% BioID_F1_F2  ~ "Yes" ,
      TRUE ~ "No"
    ) )
  
  t = data.frame(Y = c(sum(res$BioID == "Yes" & res$Regulation == "Up"),
                       sum(res$BioID == "Yes" & res$Regulation == "Down")),
                 N = c(sum(res$BioID == "No" & res$Regulation == "Up"),
                       sum(res$BioID == "No" & res$Regulation == "Down")))
  test <- chisq.test(t)
  p = test$p.value
  
  final_plot$bar = res %>%
    mutate(Regulation_ordered = factor(Regulation, levels=c("Down", "Up", "NS"))) %>%
    ggplot( aes(x = Regulation_ordered, fill = BioID )) + 
    geom_bar(position="fill") +
    scale_fill_brewer() + labs(title =  n) +
    annotate("text", x=1, y=1.1, label=paste0("p-value: ", signif(p,3)))
  
  Res_Sig_Down =res %>% 
    filter(gene_name %in% BioID_F1_F2) %>% 
    mutate(Regul_BioID = case_when(q_Value < 0.05 & Difference > 0  ~ "Up" ,
                                   q_Value < 0.05 & Difference < 0 ~ "Down",
                                   TRUE ~ "NS"
    ))
  
  # Data transformation
  df <- Res_Sig_Down %>% 
    group_by(Regul_BioID) %>% # Variable to be transformed
    count() %>% 
    ungroup() %>% 
    mutate(perc = `n` / sum(`n`)) %>% 
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))
  
  final_plot$pie = ggplot(df, 
                          aes(x = "", y = perc, fill = Regul_BioID)) +
    geom_col(color = "black") +
    geom_label(aes(label = labels), color = c(1, "white", "white"),
               position = position_stack(vjust = 0.5),
               show.legend = FALSE) +
    guides(fill = guide_legend(title = "MS")) +
    scale_fill_viridis_d() +
    coord_polar(theta = "y") + 
    theme_void()
  
  return(final_plot)
}


r_T_PC3 = plotEnrich(PC3_BioID,MS_T$PC3_MS_T,100)
r_T_LNCaP = plotEnrich(LNCaP_BioID,MS_T$LNCaP_MS_T,100)
r_T_H660 = plotEnrich(H660_BioID,MS_T$H660_MS_T,100)

pdf("./PC3_LNCaP_H660_MS_T_EnrichPlot_BarPie_p05_FC05_BioID100.pdf")
r_T_PC3$bar + r_T_PC3$pie
r_T_LNCaP$bar + r_T_LNCaP$pie
r_T_H660$bar + r_T_H660$pie
dev.off()

plotEnrich_list <- function(BioID,MS,n)
{
  BioID_F1 = BioID %>%
    slice_max(F1, n = n)
  
  BioID_F2 = BioID %>%
    slice_max(F2, n = n)
  
  BioID_F1_F2 = c(BioID_F1$Gene.names, BioID_F2$Gene.names)
  BioID_F1_F2 = unique(BioID_F1_F2)
  
  final_plot = list()
  for (n in names(MS))
  {
    res = MS[[n]]
    res = res %>%
      mutate(Regulation = case_when(
        q_Value < 0.05 & Difference > 0.5  ~ "Up" ,
        q_Value < 0.05 & Difference < 0.5 ~ "Down",
        TRUE ~ "NS"
      ),
      BioID = case_when(
        gene_name %in% BioID_F1_F2  ~ "Yes" ,
        TRUE ~ "No"
      ) )
    
    t = data.frame(Y = c(sum(res$BioID == "Yes" & res$Regulation == "Up"),
                         sum(res$BioID == "Yes" & res$Regulation == "Down")),
                   N = c(sum(res$BioID == "No" & res$Regulation == "Up"),
                         sum(res$BioID == "No" & res$Regulation == "Down")))
    test <- chisq.test(t)
    p = test$p.value
    
    final_plot[[n]]$bar = res %>%
      mutate(Regulation_ordered = factor(Regulation, levels=c("Down", "Up", "NS"))) %>%
      ggplot( aes(x = Regulation_ordered, fill = BioID )) + 
      geom_bar(position="fill") +
      scale_fill_brewer() + labs(title =  n) +
      annotate("text", x=1, y=1.1, label=paste0("p-value: ", signif(p,3)))
    
    Res_Sig_Down =res %>% 
      filter(gene_name %in% BioID_F1_F2) %>% 
      mutate(Regul_BioID = case_when(q_Value < 0.05 & Difference > 0.5  ~ "Up" ,
                                     q_Value < 0.05 & Difference < 0.5 ~ "Down",
                                     TRUE ~ "NS"
      ))
    
    # Data transformation
    df <- Res_Sig_Down %>% 
      group_by(Regul_BioID) %>% # Variable to be transformed
      count() %>% 
      ungroup() %>% 
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    
    final_plot[[n]]$pie = ggplot(df, 
                                 aes(x = "", y = perc, fill = Regul_BioID)) +
      geom_col(color = "black") +
      geom_label(aes(label = labels), color = c(1, "white", "white"),
                 position = position_stack(vjust = 0.5),
                 show.legend = FALSE) +
      guides(fill = guide_legend(title = "MS")) +
      scale_fill_viridis_d() +
      coord_polar(theta = "y") + 
      theme_void()
  }
  
  return(final_plot)
}


r_DKD = plotEnrich_list(PC3_BioID,PC3_MS,100)

pdf("./PC3_LNCaP_H660_MS_DKD_EnrichPlot_BarPie_p05_FC05_BioID100.pdf")
r_DKD$PC3_DKD_PLKO$bar + r_DKD$PC3_DKD_PLKO$pie
r_DKD$PC3_F1_PLKO$bar + r_DKD$PC3_F1_PLKO$pie
r_DKD$PC3_F2_PLKO$bar + r_DKD$PC3_F2_PLKO$pie
r_DKD$PC3_DKD_F1$bar + r_DKD$PC3_DKD_F1$pie
r_DKD$PC3_DKD_F2$bar + r_DKD$PC3_DKD_F2$pie
dev.off()

