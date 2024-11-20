#!/bin/bash
### exit when error
set -e

pipeline_version="1.2.2"

genome="GRCh38"
threads="8"
max_mem="50.GB"
igenome="/data/LORENZI/Utils/references"

input_csv="./FOXA2_project/chip_seq_20-10-2023/input_test_ReChip.csv"
out_dir="./FOXA2_project/chip_seq_20-10-2023/results_hg38_narrow_ReChip_ChipF1/"

##For broadpeak
#[ -d ${out_dir} ] || NXF_VER=22.10.1 ./nextflow run nf-core/chipseq  --fragment_size 100 -r ${pipeline_version} --outdir ${out_dir} \
#   --input ${input_csv}  -profile singularity   --single_end --genome ${genome} --igenomes_base ${igenome} --skip_consensus_peaks --max_cpus ${threads} --max_memory ${max_mem}


##For narrowpeak
[ -d ${out_dir} ] || NXF_VER=22.10.1 ./nextflow run nf-core/chipseq  --fragment_size 100 -r ${pipeline_version} --outdir ${out_dir} \
   --input ${input_csv}  -profile singularity --deseq2_vst  --single_end --genome ${genome} --igenomes_base ${igenome} --max_cpus ${threads} --max_memory ${max_mem} --narrow_peak

