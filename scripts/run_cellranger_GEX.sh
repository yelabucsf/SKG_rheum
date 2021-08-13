#!/bin/bash
#$ -cwd    ## use current working directory
#$ -j yes  ## merge stdout and stderr
#$ -o temp/log/
#$ -e temp/log/
#$ -t 1-8
tasks=(0 Lane1_1a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 Lane2_1b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 Lane3_2a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 Lane4_2b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 Lane5_3a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05 Lane6_3b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 Lane7_4a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5GEXL_2018_09_05  Lane8_4b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5GEXL_2018_09_05 ) 
input="${tasks[$SGE_TASK_ID]}"



cellranger count \
--id=$input \
--fastqs=/path/to/fastq/files/$input/ \
--transcriptome=/path/to/reference/mm10_withGFP/ \
--sample=$input \
--expect-cells=25000 \
--localcores=4 \
--localmem=256 \
--nosecondary


## End-of-job summary
qstat -j $JOB_ID
