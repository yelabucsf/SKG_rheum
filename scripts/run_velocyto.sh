#!/bin/bash
#$ -cwd    ## use current working directory
#$ -j yes  ## merge stdout and stderr
#$ -o temp/log/
#$ -e temp/log/
#$ -t 1-8


tasks=(0 \
Lane1_1a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane2_1b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane3_2a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane4_2b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane5_3a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane6_3b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane7_4a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 \
Lane8_4b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 )
  
input="${tasks[$SGE_TASK_ID]}"


mkdir scvelo_anlysis
cd scvelo_analysis

rm -r $input
mkdir $input
cd $input

bamfile=/path/to/cellranger/outputs/"$input"/outs/possorted_genome_bam.bam
samtools sort -m M -t [tagname] -O BAM -@ 4 -o cellsorted_possorted_genome_bam.bam $bamfile

velocyto run \
/path/to/cellranger/outputs/"$input"/outs/possorted_genome_bam.bam \
/path/to/reference/mm10_withGFP/genes/genes.gtf \
-m /path/to/reference/mm10_rmsk.gtf \ 
-b /path/to/cellranger/outputs/"$input"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
-o "$input"_velocyto \
--samtools-threads 4 \
--samtools-memory 64000


## End-of-job summary
qstat -j $JOB_ID
