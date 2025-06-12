# Formaggio_FOXA2
Code and Data for FOXA2 project ("Targeting FOXA1 and FOXA2 Disrupts the Lineage-Specific Oncogenic Output Program in Prostate Cancer ")

## Directory structure
```bash
└─ Analysis/:
|   └─ R/: helper scripts to process the results from ChIPseq and RNAseq pipelines and generate plots
|   └─ tests/: commandline scripts to run nf-core pipelines (ChIPseq and RNAseq) and create heatmaps (deepTools) 
└─ data/:
|   └─ meta: input files/tables to run the scripts
|   └─ Utils: Other utilities required to execute the scripts
````
## Analysis
R (version 4.3.1) scripts executed in RStudio IDE are provided to reproduce the analyses in the paper. Additionally, bash scripts containing a sequence of commands to execute nextflow nf-core pipelines and deepTools.




