#!/usr/bin/Rscript
library("ggplot2")
#BiocManager::install("tximport")
library(tximport)
library(clusterProfiler)
library("RColorBrewer")
library(dplyr)
library(ggpubr)


PC3_BioID = read.csv("./F1andF2_BioID_PC3.csv")
LNCaP_BioID = read.csv("./F1andF2_BioID_LNCaP_Filtered.csv")
H660_BioID = read.csv("./F1andF2_BioID_H660.csv")

PC3_LNCaP_H660_BioID = list(PC3_BioID,LNCaP_BioID,H660_BioID)
names(PC3_LNCaP_H660_BioID) = c("PC3_BioID","LNCaP_BioID","H660_BioID")


#ScatterPlots Correlation ====
library(ggrepel)
library(ggpp)
plotScatter_BioID <- function(BioID,title)
{

p = BioID %>%
  #mutate(Color = ifelse(q_Value_F1 < 0.05 & q_Value_F2 < 0.05   , "red4", "grey"),
  #       my_label = ifelse(gene_name %in% BioID_F1_F2 & q_Value_F1 < 0.05 & q_Value_F2 < 0.05,
  #                         gene_name, "")) %>%
  ggplot(aes(x=F1, y=F2))+ geom_point(alpha = 4/5, size = 1, color = "grey30") +
  xlab("FOXA1 - log2FoldChange") + ylab("FOXA2 - log2FoldChnage") + labs(title = title) +
  #scale_color_identity() +
  stat_cor(aes(x = F1, y = F2),
           method = "pearson", label.x = 1)+ theme_classic() + labs(colour="") +

   geom_label_repel(data=BioID %>% filter(Category != "NA"), # Filter data first
                   aes(label=Gene.names, colour = factor(Category)), #max.overlaps = 5000, 
                   size= 5,
                  box.padding = 0.8, max.overlaps = Inf) +
  
  scale_colour_manual(values = c("grey30", "royalblue")) + xlim(0,9) + ylim(0,9) +
  theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))
  #geom_label_repel(aes(label = my_label), box.padding = 0.1, 
  #                 max.overlaps = 5000, point.size = NA, size = 3, colour = "black")


return(p)
}

p_PC3 = plotScatter_BioID(PC3_LNCaP_H660_BioID$PC3_BioID,"PC3 BioID") #xylim = (0,10)
p_LNCaP = plotScatter_BioID(PC3_LNCaP_H660_BioID$LNCaP_BioID,"LNCaP BioID") #xylim = (0,8)
p_H660 = plotScatter_BioID(PC3_LNCaP_H660_BioID$H660_BioID,"H660 BioID") #xylim = (0,9)
pdf("./PC3_LNCaP_H660_F1_F2_LabelSelected_Trimmed_New_2.2.pdf", width = 10, height = 10)
p_PC3
p_LNCaP
p_H660
dev.off()

#VennDiagram ====
# install.packages("ggVennDiagram")
library(ggVennDiagram)
#install.packages("UpSetR")
library(UpSetR)
library(VennDiagram)
library(ggvenn)

plotVenn_Upset_BioID <- function(BioID,n,t)
{
  BioID_F1_list = BioID_F2_list = list()
  for (na in names(BioID))
  {
    BioID_F1 = BioID[[na]] %>%
      slice_max(F1, n = n)
    
    BioID_F1_list[[paste0(na,"_F1")]] = BioID_F1$Gene.names
    
    BioID_F2 = BioID[[na]] %>%
      slice_max(F2, n = n)
    
    BioID_F2_list[[paste0(na,"_F2")]] = BioID_F2$Gene.names
  }
  names(BioID_F1_list) = names(BioID_F2_list) = c("PC3", "LNCaP","H660")
  
  p = list()
  
  # Venn diagram with white
  p$F1$VennBlank = ggvenn(BioID_F1_list, 
                     fill_color = c("white", "white", "white"),
                     #fill_alpha = 0.3,
                     stroke_size = 0.5, set_name_size = 4
  ) + labs(title =  t, subtitle = "FOXA1")
  p$F2$VennBlank = ggvenn(BioID_F2_list, 
                     fill_color = c("white", "white", "white"),
                     #fill_alpha = 0.3,
                     stroke_size = 0.5, set_name_size = 4
  ) + labs(title =  t, subtitle = "FOXA2")
  
  # Venn diagram with custom colors
  p$F1$Venn = ggvenn(BioID_F1_list, 
    fill_color = c("orange", "lightpink1", "mediumpurple1"),
    fill_alpha = 0.3,
    stroke_size = 0.5, set_name_size = 4
  ) + labs(title =  t, subtitle = "FOXA1")
  p$F2$Venn = ggvenn(BioID_F2_list, 
                     fill_color = c("orange", "lightpink1", "mediumpurple1"),
                     fill_alpha = 0.3,
                     stroke_size = 0.5, set_name_size = 4
  ) + labs(title =  t, subtitle = "FOXA2")
  
  # Upsetplot 
  p$F1$Upset = upset(fromList(BioID_F1_list), order.by = "freq",
                      set_size.show = FALSE)
    #scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
    #labs(title =  t, subtitle = "FOXA1")
  
  p$F2$Upset = upset(fromList(BioID_F2_list), order.by = "freq",
                     set_size.show = FALSE)
    #scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
    #labs(title =  t, subtitle = "FOXA2")
  
  return(p)
}

F_F_50 = plotVenn_Upset_BioID(PC3_LNCaP_H660_BioID,50,"BioID Top50")
F_F_100 = plotVenn_Upset_BioID(PC3_LNCaP_H660_BioID,100,"BioID Top100")

pdf("./PC3_LNCaP_H660_Venn_newcolor_BioID100.pdf")
F_F_100$F1$Venn 
F_F_100$F1$VennBlank
F_F_100$F2$Venn 
F_F_100$F2$VennBlank
dev.off()

pdf("./PC3_LNCaP_H660_Venn_Upset_BioID50_BioID100.pdf")
F_F_50$F1$Venn  
F_F_50$F1$Upset
F_F_50$F2$Venn
F_F_50$F2$Upset
F_F_100$F1$Venn 
F_F_100$F1$Upset
F_F_100$F2$Venn
F_F_100$F2$Upset
dev.off()

#GSEA GO for subgroups ====
library(BioVenn)
bioVenn_BioID <- function(BioID,n)
{
  BioID_F1_list = BioID_F2_list = list()
  for (na in names(BioID))
  {
    BioID_F1 = BioID[[na]] %>%
      slice_max(F1, n = n)
    
    BioID_F1_list[[paste0(na,"_F1")]] = BioID_F1$Gene.names
    
    BioID_F2 = BioID[[na]] %>%
      slice_max(F2, n = n)
    
    BioID_F2_list[[paste0(na,"_F2")]] = BioID_F2$Gene.names
  }
  
  BioID_F1_F2 = list(F1 = BioID_F1_list,F2 = BioID_F2_list)
  
  return(BioID_F1_F2)
}

BioID_F1_F2 = bioVenn_BioID(PC3_LNCaP_H660_BioID,100)
  
biovenn_BioID_F1 <- draw.venn(BioID_F1_F2$F1$PC3_BioID_F1, BioID_F1_F2$F1$LNCaP_BioID_F1, 
                              BioID_F1_F2$F1$H660_BioID_F1, 
                             title =  "bioID FOXA1", 
                             subtitle = "Top 100", 
                             xtitle = "", 
                             ytitle = "", 
                             ztitle = "",
                             nrtype="abs", nr_s = 2, nr_fb = 1,
                             x_c = "grey",y_c = "skyblue",z_c = "royalblue",bg_c = "white",
                             output = "creen", filename = NULL)



write.csv(biovenn_BioID_F1$x_only,paste0("./FOXA1_PC3_Only.csv"), row.names = FALSE )
write.csv(biovenn_BioID_F1$y_only,paste0("./FOXA1_LNCaP_Only.csv"), row.names = FALSE)
write.csv(biovenn_BioID_F1$z_only,paste0("./FOXA1_H660_Only.csv"), row.names = FALSE)

write.csv(c(biovenn_BioID_F1$xy_only,biovenn_BioID_F1$xz_only,biovenn_BioID_F1$yz_only, biovenn_BioID_F1$xyz),
          paste0("./FOXA1_PC3_LNCaP_H660_Union.csv"), row.names = FALSE )

write.csv(biovenn_BioID_F1$xyz,
          "./FOXA1_PC3_LNCaP_H660_Overlap.csv", row.names = FALSE)


biovenn_BioID_F2 <- draw.venn(BioID_F1_F2$F2$PC3_BioID_F2, BioID_F1_F2$F2$LNCaP_BioID_F2, 
                              BioID_F1_F2$F2$H660_BioID_F2, 
                              title =  "bioID FOXA2", 
                              subtitle = "Top 100", 
                              xtitle = "", 
                              ytitle = "", 
                              ztitle = "",
                              nrtype="abs", nr_s = 2, nr_fb = 1,
                              x_c = "grey",y_c = "skyblue",z_c = "royalblue",bg_c = "white",
                              output = "creen", filename = NULL)

write.csv(biovenn_BioID_F2$x_only,paste0("./FOXA2_PC3_Only.csv"), row.names = FALSE )
write.csv(biovenn_BioID_F2$y_only,paste0("./FOXA2_LNCaP_Only.csv"), row.names = FALSE)
write.csv(biovenn_BioID_F2$z_only,paste0("./FOXA2_H660_Only.csv"), row.names = FALSE)

write.csv(c(biovenn_BioID_F2$xy_only,biovenn_BioID_F2$xz_only,biovenn_BioID_F2$yz_only, biovenn_BioID_F2$xyz),
          paste0("./FOXA2_PC3_LNCaP_H660_Union.csv"), row.names = FALSE )

write.csv(biovenn_BioID_F2$xyz,
          "./FOXA2_PC3_LNCaP_H660_Overlap.csv", row.names = FALSE)


FOXA1_PC3_LNCaP_H660 = list(read.csv("./FOXA1_PC3_Only.csv",  col.names = c("FOXA1_PC3_Only"))$FOXA1_PC3_Only,
                           read.csv("./FOXA1_LNCaP_Only.csv",  col.names = c("FOXA1_LNCaP_Only"))$FOXA1_LNCaP_Only,
                           read.csv("./FOXA1_H660_Only.csv",  col.names = c("FOXA1_H660_Only"))$FOXA1_H660_Only,
                           read.csv("./FOXA1_PC3_LNCaP_H660_Union.csv",  col.names = c("FOXA1_PC3_LNCaP_H660_Union"))$FOXA1_PC3_LNCaP_H660_Union,
                           read.csv("./FOXA1_PC3_LNCaP_H660_Overlap.csv",  col.names = c("FOXA1_PC3_LNCaP_H660_Overlap"))$FOXA1_PC3_LNCaP_H660_Overlap)
names(FOXA1_PC3_LNCaP_H660) = c("PC3_Only","LNCaP_Only", "H660_Only", "Union",
                               "Overlap")

FOXA2_PC3_LNCaP_H660 = list(read.csv("./FOXA2_PC3_Only.csv",  col.names = c("FOXA2_PC3_Only"))$FOXA2_PC3_Only,
                         read.csv("./FOXA2_LNCaP_Only.csv",  col.names = c("FOXA2_LNCaP_Only"))$FOXA2_LNCaP_Only,
                         read.csv("./FOXA2_H660_Only.csv",  col.names = c("FOXA2_H660_Only"))$FOXA2_H660_Only,
                         read.csv("./FOXA2_PC3_LNCaP_H660_Union.csv",  col.names = c("FOXA2_PC3_LNCaP_H660_Union"))$FOXA2_PC3_LNCaP_H660_Union,
                         read.csv("./FOXA2_PC3_LNCaP_H660_Overlap.csv",  col.names = c("FOXA2_PC3_LNCaP_H660_Overlap"))$FOXA2_PC3_LNCaP_H660_Overlap)
names(FOXA2_PC3_LNCaP_H660) = c("PC3_Only","LNCaP_Only", "H660_Only", "Union","Overlap")

Overlap = list(FOXA1 = FOXA1_PC3_LNCaP_H660$Overlap, FOXA2 = FOXA2_PC3_LNCaP_H660$Overlap )
Overlap = list(FOXA1_FOXA2_Excl = intersect(FOXA1_PC3_LNCaP_H660$Overlap,FOXA2_PC3_LNCaP_H660$Overlap))
               #FOXA1_FOXA2_Incl = c(intersect(FOXA1_PC3_LNCaP_H660$Overlap,FOXA2_PC3_LNCaP_H660$Overlap),"FOXA1","FOXA2"))


library(tidyverse)  
library(msigdbr)

gsets = list( 
  HALLMARKS = msigdbr(species="Homo sapiens", category = "H"),
  KEGG = msigdbr(species="Homo sapiens", category = "C2", subcategory = "CP:KEGG"),
  REACTOME= msigdbr(species="Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
  PID= msigdbr(species="Homo sapiens", category = "C2", subcategory = "CP:PID"),
  BP= msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:BP"),
  CC= msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:CC"),
  MF= msigdbr(species="Homo sapiens", category = "C5", subcategory = "GO:MF")
)


plotEnricher <- function(PC3_LNCaP_H660,Bait,ggset,gsub_ext)
{
  PC3_LNCaP_H660_GO = list()
  TERM2GENE_Sel = gsets[[ggset]]
  for ( n in names(PC3_LNCaP_H660))
  {
    tmp = enricher(PC3_LNCaP_H660[[n]], pvalueCutoff  = 0.05,
                   pAdjustMethod = "BH", TERM2GENE = TERM2GENE_Sel[,c("gs_name", "gene_symbol")])
    
    PC3_LNCaP_H660_GO[[n]]= tmp@result
  }
  
  merged.res = as.data.frame(merge_result(PC3_LNCaP_H660_GO))
  merged.res$SIG = -log10(merged.res$p.adjust)
  merged.res$GOBP = gsub(gsub_ext, "", merged.res$Description)
  
  plots = list()
  plots$p <- ggplot(merged.res, aes(x = Cluster, y = fct_reorder(GOBP, SIG))) + 
    geom_point(aes(size = Count, color = SIG)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 10),axis.ticks.x=element_blank(),
          axis.title.x=element_blank(), axis.text.y = element_text(size = 8), plot.title = element_text(hjust = 0.5)) +
    scale_colour_gradient2(limits=c(-3.5, 30), low = "#3953A4", mid = "grey", high = "#b50404") +
    ylab(NULL) + 
    scale_size(range = c(3,8)) + labs(color ='-log10(adj.p)',size = "Count",  title = paste0(Bait,'_PC3_LNCaP_H660 BioID100 GO',ggset)) 
  #p
  categories = c("ANDROGEN_RECEPTOR_SIGNALING_PATHWAY",
                 "HORMONE_MEDIATED_SIGNALING_PATHWAY",
                 "INTRACELLULAR_RECEPTOR_SIGNALING_PATHWAY",
                 "RESPONSE_TO_STEROID_HORMONE",
                 "PROTEIN_SUMOYLATION"
                 )
  
  plots$p_sig <- ggplot(subset(merged.res, p.adjust < 0.05 & GOBP %in% categories), aes(x = Cluster, y = fct_reorder(GOBP, SIG))) + 
    geom_point(aes(size = Count, color = SIG)) +
    #theme_bw(base_size = 10) +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin(), axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),axis.text.x = element_text(size = 12), axis.title.y = element_blank(),  axis.text.y=element_blank(),
          plot.title = element_text(hjust = 0.5),
          ) +
    #scale_colour_gradient2(limits=c(-3.5, 30), low = "#3953A4", mid = "grey", high = "#b50404") +
    scale_colour_gradient2( low = "#3953A4", high = "#b50404") +
    ylab(NULL)  + guides(size = guide_legend(order=1)) + 
    coord_flip() + theme_classic() +
    scale_size(range = c(3,10)) + labs(color ='-log10(adj.p)',size = "Count") #,  title = paste0(Bait,'_PC3_LNCaP_H660 BioID100 GO',ggset)) 
  #p_sig
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, size = 10),
  
  #return (plots)
  return(PC3_LNCaP_H660_GO)
}

p_FOXA1 = plotEnricher(FOXA1_PC3_LNCaP_H660,"FOXA1","BP","GOBP_")
p_FOXA1_CC = plotEnricher(FOXA1_PC3_LNCaP_H660,"FOXA1","CC","GOCC_")
p_FOXA1_MF = plotEnricher(FOXA1_PC3_LNCaP_H660,"FOXA1","MF","GOMF_")
pdf("./FOXA1_PC3_LNCaP_H660_GO_All.pdf", height = 100 , width = 20)
p_FOXA1$p
p_FOXA1_CC$p
p_FOXA1_MF$p
dev.off()

pdf("./FOXA1_PC3_LNCaP_H660_GO_Sig.pdf", height = 20 , width = 20)
p_FOXA1$p_sig
p_FOXA1_CC$p_sig
p_FOXA1_MF$p_sig
dev.off()

p_FOXA2 = plotEnricher(FOXA2_PC3_LNCaP_H660,"FOXA2", "BP","GOBP_")
p_FOXA2_CC = plotEnricher(FOXA2_PC3_LNCaP_H660,"FOXA2", "CC","GOCC_")
p_FOXA2_MF = plotEnricher(FOXA2_PC3_LNCaP_H660,"FOXA2", "MF","GOMF_")
pdf("./FOXA2_PC3_LNCaP_H660_GO_All.pdf", height = 100 , width = 20)
p_FOXA2$p
p_FOXA2_CC$p
p_FOXA2_MF$p
dev.off()

pdf("./FOXA2_PC3_LNCaP_H660_GO_Sig.pdf", height = 20 , width = 20)
p_FOXA2$p_sig
p_FOXA2_CC$p_sig
p_FOXA2_MF$p_sig
dev.off()

res = plotEnricher(Overlap,"F1_F2","BP","GOBP_")
res_PC3_LNcaP_H660_GOBP = res$FOXA1_FOXA2_Excl
write.csv(res_PC3_LNcaP_H660_GOBP, "./FOXA1_FOXA2_Overlap_PC3_LNCaP_H660_Overlap_GOBP.csv", row.names = FALSE)
p_FOXA1_FOXA2_BB = plotEnricher(Overlap,"F1_F2","BP","GOBP_")
pdf("./FOXA1_FOXA2_Overlap_PC3_LNCaP_H660_Overlap_GOBP_Sel.pdf", height = 3.5 , width = 5)
p_FOXA1_FOXA2_BB$p_sig
dev.off()
pdf("./FOXA1_FOXA2_Overlap_PC3_LNCaP_H660_Overlap_GOBP_Sel_legend.pdf", height = 3.5 , width = 8)
p_FOXA1_FOXA2_BB$p_sig
dev.off()


p_FOXA1_H = plotEnricher(FOXA1_PC3_LNCaP_H660,"FOXA1","HALLMARKS","HALLMARK_")
pdf("./FOXA1_PC3_LNCaP_H660_HALLMARK.pdf")
p_FOXA1_H$p
p_FOXA1_H$p_sig
dev.off()

p_FOXA2_H = plotEnricher(FOXA2_PC3_LNCaP_H660,"FOXA2","HALLMARKS","HALLMARK_")
pdf("./FOXA2_PC3_LNCaP_H660_HALLMARK.pdf")
p_FOXA2_H$p
p_FOXA2_H$p_sig
dev.off()

#Heatmap / table for overlapping proteins ====
FOXA1_Overlap = read.csv("./FOXA1_PC3_LNCaP_H660_Overlap.csv", 
                         col.names = c("FOXA1_PC3_LNCaP_H660_Overlap"))$FOXA1_PC3_LNCaP_H660_Overlap
FOXA2_Overlap = read.csv("./FOXA2_PC3_LNCaP_H660_Overlap.csv", 
                         col.names = c("FOXA2_PC3_LNCaP_H660_Overlap"))$FOXA2_PC3_LNCaP_H660_Overlap

#Overlap = c(intersect(FOXA1_Overlap,FOXA2_Overlap),"FOXA1","FOXA2")
Overlap = intersect(FOXA1_Overlap,FOXA2_Overlap)
Overlap_BioID = list(PC3 = PC3_BioID %>%
                       filter(Gene.names %in% Overlap),
                     H660 = H660_BioID %>%
                       filter(Gene.names %in% Overlap),
                     LNCaP = LNCaP_BioID %>%
                       filter(Gene.names %in% Overlap))

library(tidyverse)
Overlap_BioID_tab = Overlap_BioID %>% reduce(full_join, by='Gene.names',
                                             suffix = c("_PC3","_H660"))

#Manual selection - 2 entries for NCOR1
Overlap_BioID_tab = #Overlap_BioID_tab[1:20,] %>%
  Overlap_BioID_tab[1:18,] %>%
  select(Gene.names,F1_PC3,F2_PC3,F1_H660,F2_H660,F1,F2)
names(Overlap_BioID_tab) = c("Gene.names","F1_PC3","F2_PC3","F1_H660",
                             "F2_H660","F1_LNCaP","F2_LNCaP")

annotation_col = data.frame(
  Bait = rep(c("FOXA1","FOXA2"),3),
  Cellline = c(rep("PC3",2),rep("H660",2),rep("LNCaP",2)))
rownames(annotation_col) = names(Overlap_BioID_tab)[-1]


annotation_row = data.frame(
  Type = c(rep("TFs",5),rep("TF activators",2),rep("SWI/SNF",2),
           rep("SUMO ligases",3),rep("Co-repressors",6)))

rownames(annotation_row) =c("PRR12", "AHDC1", "CUX1", "ZNF384", "MYBL2", 
                            "TCF20", "JMJD1C", "ARID1A", "ARID1B","PIAS1", "PIAS2", "PIAS3",
                            "BCOR", "NCOR1", "NCOR2", "PHF21A", "TLE3", "ELMSAN1")
#ELMSAN1 / MIDEAS

ann_colors = list(
  Bait = c(FOXA1="red4", FOXA2="blue"),
  Cellline = c(PC3 = "khaki2",H660 = "mediumpurple1", LNCaP = "lightpink1"),
  Type = c(TFs = "salmon2", `TF activators` = "orange4",`SWI/SNF` = "yellow4",
           `SUMO ligases` = "plum4", `Co-repressors` = "wheat2" ))

Overlap_BioID_tab = column_to_rownames(Overlap_BioID_tab, var = "Gene.names")
Overlap_BioID_tab = Overlap_BioID_tab[rownames(annotation_row), ]
dat_mat = Overlap_BioID_tab
dat_mat_scaled = scale(dat_mat)

library(pheatmap) #will not work if tidyverse is called after
p = pheatmap(dat_mat, cellheight = 12, color = colorRampPalette(c("white","firebrick3"))(50),
             annotation_col = annotation_col,annotation_colors = ann_colors, annotation_row = annotation_row,
         cluster_rows = F, cluster_cols = F,annotation_names_row = F, cellwidth = 10,
         show_colnames = F, border_color = "NA", gaps_col = c(2,4))
#brewer.pal(n = 9, name = "Blues")

library(gridExtra)
library(grid)
BioID_type_tab = data.frame( Type = c(rep("TFs",7),rep("TF activators",2),rep("SWI/SNF",2),
                                          rep("SUMO ligases",3),rep("Co-repressors",6)),
                                 Proteins = c("PRR12", "AHDC1", "CUX1", "ZNF384", "MYBL2", "FOXA1", "FOXA2",
                                              "TCF20", "JMJD1C", "ARID1A", "ARID1B","PIAS1", "PIAS2", "PIAS3",
                                              "BCOR", "NCOR1", "NCOR2", "PHF21A", "TLE3", "ELMSAN1"))

BioID_type_tab = data.frame( Type = c("TFs","TF activators","SWI/SNF","SUMO ligases","Co-repressors"),
                             Proteins = c("PRR12, AHDC1, CUX1, ZNF384, MYBL2, FOXA1, FOXA2",
                             "TCF20, JMJD1C", "ARID1A, ARID1B" , "PIAS1, PIAS2, PIAS3",
                                          "BCOR, NCOR1, NCOR2, PHF21A, TLE3, ELMSAN1"))


pdf("./FOXA1_FOXA2_PC3_LNCaP_H660_Heatmap.pdf", height = 8 , width = 6)
p
dev.off()

pdf("./FOXA1_FOXA2_PC3_LNCaP_H660_Table.pdf", height = 4 , width = 6)
grid.table(BioID_type_tab, rows = NULL)
dev.off()
