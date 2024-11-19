# Formaggio_FOXA2
Code and Data for FOXA2 project

## Directory structure
└─ notebooks/:
|   └─ R/: helper scripts to process the results from ChIPseq and RNAseq pipelines and create plots
|   └─ tests/: commandline scripts to run nf-core pipelines (ChIPseq and RNAseq) and create heatmaps (deepTools) 
|   └─ BioID_Score.R : helper script to calculate ssGSEA of BioID core proteins
└─ data/:
|   └─ meta: input files/tables to run the scripts
|   └─ Utils: Other utilities required to execute the scripts
