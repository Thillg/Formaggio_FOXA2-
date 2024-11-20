#!/bin/bash

BASE_UTILS="/data/LORENZI/Utils"
BASE_DIR=$(realpath ./RNA_seq_19_07_2024/)
OUT_DIR="${BASE_DIR}/output"

pipeline_version="3.6"
threads="32"
max_mem="100.GB"
GENCODE_VERSION="39"

reference_gtf="${BASE_UTILS}/gencode_v${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf"
reference_genome="${BASE_UTILS}/gencode_v${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa"

[ -d ${OUT_DIR} ] || NXF_VER=22.10.1 ./nextflow run nf-core/rnaseq -profile singularity -r ${pipeline_version} --aligner star_salmon --max_cpus ${threads} --max_memory ${max_mem} --input input.csv \
                   --outdir ${OUT_DIR}/  --fasta ${reference_genome}  --gtf ${reference_gtf} --gencode

rm -fr ./work
