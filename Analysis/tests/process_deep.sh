#!/bin/bash

THREADS=16
source /data/GT/miniconda3/etc/profile.d/conda.sh
simg=/data/GT/singularity/nfcore-chipseq-1.2.2.img

FOXA1_1=/data/LORENZI/ChIP_seq_22_11_15/results_hg19/bwa/mergedLibrary/bigwig/LnCaP1_FOXA1_R1.bigWig
FOXA1_2=/data/LORENZI/ChIP_seq_22_11_15/results_hg19/bwa/mergedLibrary/bigwig/LnCaP1_FOXA1_R2.bigWig

computeMatrix reference-point --averageTypeBins max --missingDataAsZero --referencePoint center -a 1000 -b 1000 -p ${THREADS} --outFileName "./FOXA1_R1.matrix.gz" -S $FOXA1_1 $FOXA1_F1 --samplesLabel "FOXA1_WT_R1" "FOXA1_OE_R1" -R FOXA2_with_FOXA1_H3_decreasing_sorted.bed FOXA2_with_FOXA1_H3_increasing_sorted.bed

plotHeatmap -m ./FOXA1_R1.matrix.gz -o ./FOXA1_R1_sorted.pdf --sortRegions descend --legendLocation none
plotProfile -m ./FOXA1_R1.matrix.gz -out ./FOXA1_R1_sorted_ymax4_mergedplotProf.png --perGroup --colors red blue --yMax 5 5
plotProfile -m ./FOXA1_R1.matrix.gz -out ./FOXA1_R1_sorted_mergedplotProf.png --perGroup --colors red blue

#Updated colors
plotHeatmap -m ./png/PC3_F1_Merged/PC3.matrix.gz -o ./png/PC3_F1_Merged/PC3_Reds_short.pdf --sortRegions descend --legendLocation none --colorMap Reds --heatmapHeight 15
plotHeatmap -m ./png/PC3_F2_Merged/PC3.matrix.gz -o ./png/PC3_F2_Merged/PC3_F2_Blues_short.pdf --sortRegions descend --legendLocation none --colorMap Blues --heatmapHeight 15
plotHeatmap -m ./png/PC3_H3K_Merged/PC3.matrix.gz -o ./png/PC3_H3K_Merged/PC3_H3K_Greens_short.pdf --sortRegions descend --legendLocation none --colorMap Greens  --heatmapHeight 15

plotProfile -m ./png/PC3_F1_Merged/PC3.matrix.gz -out ./png/PC3_F1_Merged/PC3_F1_prof_newcolor.pdf --perGroup --colors black orange darkcyan
plotProfile -m ./png/PC3_DMSO_T_F1_R1.matrix.gz -out ./png/PC3_DMSO_T_F1_R1_prof_newcolor.pdf --perGroup --colors black darkred


