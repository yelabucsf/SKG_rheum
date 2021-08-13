#!/bin/bash
#$ -cwd    ## use current working directory
#$ -j yes  ## merge stdout and stderr
#$ -o temp/log/
#$ -e temp/log/
#$ -t 1-8
tasks=(0 Lane1_1a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5TEL_2018_11_28 Lane2_1b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5TEL_2018_11_28 Lane3_2a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5TEL_2018_11_28 Lane4_2b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5TEL_2018_11_28 Lane5_3a_SKG_CD4Naive_GFPhi_Aug17_10xVDJ_5TEL_2018_11_28 Lane6_3b_SKG_CD4Naive_GFPlo_Aug17_10xVDJ_5TEL_2018_11_28 Lane7_4a_WT_CD4Naive_GFPhi_Aug17_10xVDJ_5TEL_2018_11_28 Lane8_4b_WT_CD4Naive_GFPlo_Aug17_10xVDJ_5TEL_2018_11_28 ) 
input="${tasks[$SGE_TASK_ID]}"



cellranger vdj \
--id=$input \
--fastqs=/path/to/fastq/files/$input/ \
--reference=refdata-cellranger-vdj-GRCm38-alts-ensembl-3.1.0/ \
--sample=$input \
--localcores=4 \
--localmem=256


## End-of-job summary
qstat -j $JOB_ID
