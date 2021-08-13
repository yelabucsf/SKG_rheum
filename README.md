---
title: "Arthritogenic T cells have a transcriptional program of activation and a repertoire pruned by superantigen."
date: "August 13, 2021"
output:
  html_document:
    toc: true
    keep_md: true
---

This repo contains the code for the analyses and figures from "Arthritogenic T cells have a transcriptional program of activation and a repertoire pruned by superantigen". 

This document is divided into six sections. The input and output files and jupyter notebooks are listed and described first (1. Directory). The next sections describe the experiment and analysis for the bulk RNA sequencing data (2. Bulk RNA Sequencing Analysis) and for the single cell RNA sequencing data in three sections (3. Single Cell RNA Seq - Cell sub-type and T.4N_Nr4a1 Analysis, 4. Trajectory Analysis, and 5 TCR analysis). For each section, the jupyter notebooks that go along with each analysis are listed along with the section headers within the notebook to facilitate easily finding code for a particular figure/analysis. The last section (6. Other Software Versions) details software versions not provided in the previous sections.


### 1. Directory

    Jupyter Notebooks

**1_SKG_RA_bulk_RNA_seq_analysis.ipynb**  
**2_SKG_RA_single_cell_preprocessing_mouse.ipynb**  
**3_SKG_RA_single_cell_clustering.ipynb**  
**4_SKG_RA_single_cell_sub_type_profiling.ipynb**  
**5_SKG_RA_single_cell_T_4_Nr4a1_analysis.ipynb**  
**6_SKG_RA_single_cell_trajectory_analysis.ipynb**  
**7_SKG_RA_single_cell_TRA_clonotype.ipynb**  
**8_SKG_RA_single_cell_TRBV.ipynb**  
**9_SKG_RA_single_cell_MAST.ipynb**  



    /adata_object
**adata_only_T_cells.h5ad**: *scanpy anndata object with processed data*  
**adata_only_T_cells_scvelo_Nr4a1_high_cluster_dynamical_050321.h5ad**: *scvelo anndata object with trajectory analysis*  


    /custom_reference_input_files
    
**GFP_gene_info.txt**: *eGfp transcript info*  
**GFP_sequence.txt**: *eGfp transcript sequence*
    
    /data


**bulk_rna_seq_meta_data.csv**: *meta data for bulk RNA seq count matrix*  
**bulk_RNAseq_data.csv**: *bulk RNA seq count matrix*  

 
**sc_RNA_sample_sheet.csv**: *meta data for single cell RNA seq samples*  
**non_T_barcodes.npy**: *numpy file with contaminating B cell barcodes to remove for scRNA-seq data*   
**G1_S.csv**: *Mus musculus G1/S gene list from REACTOME*  
**G2_M.csv**: *Mus musculus G2/M gene list from REACTOME*  
**LN_protein_data.csv**: *Vb frequency for each subgroup*  
**Vb11_combo_joint.csv**   
**Vb14_combo_joint.csv**  
**Vb3_combo_joint.csv**  
**Vb5_combo_joint.csv**  
**Vb6_combo_joint.csv**  
**Vb8_combo_joint.csv**  


    /results 
- /bulk_RNA_seq

Note: data_S1_diff_exp_Group1_Group2.csv is for Group 1 v Group 2 i.e. positive log2FoldChange is for genes UP in Group 1 versus Group 2  
Column Annotations:  
X - Gene name  
baseMean—The average of the normalized count values, dividing by size factors, taken over all samples  
log2FoldChange–The effect size estimate  
lfcSE–The standard error estimate for the log2 fold change estimate  
stat–The value of the test statistic for the gene or transcript  
pvalue–P-value of the test for the gene or transcript   
padj–Adjusted P-value for multiple testing for the gene or transcript

**data_S1_diff_exp_WT_low_SKG_low.csv**
**data_S1_diff_exp_WT_low_SKG_low.csv**
**data_S1_diff_exp_WT_high_SKG_high.csv**
**data_S1_diff_exp_SKG_low_SKG_high.csv**

**data_S2_heatmap_gene_list_with_modules.csv**: *Ordered list of genes in heat map with module annotations*  
**Gene_list_for_heatmap_2021_02_v1_orig.csv**: *List of genes to annotate in the heatmap*  
**Gene_list_for_heatmap_2021_02_v1.csv**: *Same file as above with added column "module" to indicate module assignment*  

**dotplot_FEA_pathways_for_gene_modules.csv**: *Curated list of enriched GO:BP or KEGG pathways for each gene module*

**fig_1_diff_exp_SKG_high_v_WT_High_ranked.rnk**: - *Ranked list for differential expression of SKG High v WT High*     

|
|


- /single_cell_RNA_seq
  + /correlations
  
**gene_correlations_for_hvg_mouse_all_cells_spearman.csv**      
**gene_correlations_for_hvg_mouse_Nr4a1_hi_cluster_spearman.csv**  

  + /differential_expression
  
**data_S4_scRNAseq_diff_genes_by_cluster.csv**

**data_S5_diff_genes_Nr4a1_high_cluster_versus_other_cells.csv**

**data_S5_diff_genes_SKG_High_v_WT_High_Nr4a1_cluster.csv**

**data_S6_diff_exp_Tnfrsf9_pos_Egr2_pos_Nr4a1_hi_cluster.csv**



  + /ranked_lists

**SKG_High_v_WT_High_Nr4a1_high_cluster.rnk**

**Nr4a1_high_cluster_Egr_high_v_Tnf_high_genes.rnk**

**Nr4a1_cluster_stage_1_v_stage_4.rnk**


  + /trajectory

**data_S7_top_300_heatmap_gene_list.txt**  

  + /TCR

**data_S8_gini_coefficients.csv**  
**TCR_data.pickle**  
**data_S9_TRBV_paired_test_nr4a1_cluster.csv**  
**select_TRBV_sample_frequencies_Nr4a1_high_cluster.csv**  
**select_TRBV_sample_frequencies_all_cells.csv**  
**volcano_plot_gene_list_to_label.csv**  
**TRBV3_MAST_SKG_high_v_SKG_low.csv**  
**TRBV19_MAST_SKG_high_v_SKG_low.csv**  


    /scripts

**run_cellranger_GEX.sh**   
**run_cellranger_TCR.sh**  
**run_velocyto.sh**  

    

### 2. Bulk RNA Sequencing Analysis

    Sequencing (3 batches run on 3 lanes of HiSeq 2500)

Batch 1 (H94843)  
2a_SKGNur_CD4Naive_GFPlo  
2b_SKGNur_CD4Naive_GFPhi  
3a_SKGNur_CD4Naive_GFPlo  
4a_WTNur_CD4Naive_GFPlo  
4b_WTNur_CD4Naive_GFPhi


Batch 2 (H95020)

5a_WTNur_CD4Naive_GFPlo  
5b_WTNur_CD4Naive_GFPhi  
6a_WTNur_CD4Naive_GFPlo  
6b_WTNur_CD4Naive_GFPhi  


Batch 3 (H96272)  
1a_SKGNur_CD4Naive_GFPlo  
1b_SKGNur_CD4Naive_GFPhi  
3b_SKGNur_CD4Naive_GFPhi  

    Results  

[] Add sequencing reads per sample  


    Analysis

**1_SKG_RA_bulk_RNA_seq_analysis.ipynb**

Software versions:

  
python: 
pandas v.1.1.3
numpy v.1.19.2
rpy2 v.2.9.4
matplotlib v.3.3.1

R: 
VennDiagram v.1.6.20
pheatmap v.1.0.12
hash v.2.2.6.1
ggplot2 v.3.3.2
DESeq2 v.1.22.2

Sections:  
  + Data Processing  
  + PCA Analysis  
  + Differential Expression  
  + Heatmap  
  + GO Plot  
  + Module Distribution  
  + Volcano plot SKG High v WT High 
  + Ranked list for GSEA  
  + Other

### 3. Single Cell RNA Seq - Cell sub-type and T.4N_Nr4a1 Analysis

    Sequencing (8 wells of 5'10x-VDJ sequenced on NovaSeq 6000)

**GEX**

Sample | PE Sequencing Reads
-------|---------------------
1       |745,890,946
2       |727,971,779
3       |784,588,458
4       |747,029,488
5       |840,784,64
6       |831,011,967
7       |892,260,945
8       |715,376,870



**TCR**


Sample|#PE Sequencing Reads
------|--------------
1|209,571,333
2|438,975,429
3|403,195,137
4|366,362,663
5|419,581,007
6|278,298,845
7|510,594,226
8|404,954,306


    Alignment

**GEX**

/scripts/run_cellranger_GEX.sh  
Used cellranger v.3.0.1 and mm10_withGFP transcriptome.

mm10_withGFP transcriptome creation: 

  + Input files: /custom_reference_input_files  


```bash
#Concatenate GFP sequence and GFP gene description to files from refdata-cellranger-mm10-3.0.0
cat GFP_sequence.txt>>genome.fa
cat GFP_gene_info.txt>>genes.gtf
#Create mm10_withGFP reference
cellranger mkref --genome=mm10_withGFP --fasta=genome.fa --genes=genes.gtf
```
**TCR**

/scripts/run_cellranger_TCR.sh

    Analysis

**2_SKG_RA_single_cell_preprocessing_mouse.ipynb**

Software versions:
scanpy v.1.4.3

Preprocessing to create scanpy objects for each 10x lane cellranger output

**3_SKG_RA_single_cell_clustering.ipynb**


Software versions:
scanpy v.1.4.3

Data normalization and clustering

**4_SKG_RA_single_cell_sub_type_profiling.ipynb**

Software versions:
scanpy v.1.5.1

Sections:    

  + Create cell sub types and run diff exp  
  + UMAPs by subgroup  
  + Density UMAPs by subgroup  
  + Stacked Violin Plots and Matrix Plots  
  + Dot Plot for Cell Type Markers  
  + Distribution over cell sub-types by subgroup  
  + Scoring for bulk RNA Seq modules  
  + Calculate Cell Cycle  
  + Save adata object  


**5_SKG_RA_single_cell_T_4_Nr4a1_analysis.ipynb**

Software versions:
scanpy v.1.7.1

Sections:

  + UMAP colored by T.4 Nr4a1 cluster  
  + Differential Expression for T.4 Nr4a1 cluster
  + Volcano Plot for T.4N Nr4a1 cluster v other cells
  + UMAP for SKG High v WT High in T.4N Nr4a1 cluster
  + Differential expression for SKG High v WT High in T.4N Nr4a1 high cluster
  + Volcano Plot for SKG High v WT High in T.4N Nr4a1 cluster
  + UMAPs for gene markers
  + Correlation Heatmaps
  + Volcano Plot Egr2 High v Tnfrsf9 High
  + Cell Cycle Analysis
  

### 4. Trajectory Analysis

    Analysis Pipeline

Create loom files from cellranger BAM files using velocyto (/scripts/run_velocyto.sh). Use loom files as input for scvelo.
  
    Analysis  

**6_SKG_RA_single_cell_trajectory_analysis.ipynb**


Software versions:
velocyto v.0.17.17 loompy v.2.0.17 scvelo v.0.2.1 scanpy v.1.5.1

Sections:

  + UMAPs for Egr2 and Tnfrsf9 expression  
  + Merge data and run scvelo in dynamical mode  
  + UMAP overlays  
  + Visualize latent time distributions  
  + Separate stages of latent time distribution with Gaussian Mixture Model  
  + Heatmap of expression top genes for modelling latent time  
  + Differential expression between Stage 1 and Stage 4  
  + Visualize smoothed gene expression over latent time
  + Run PAGA
  + Save adata object

### 5. TCR Analysis

    Analysis
    
    
**7_SKG_RA_single_cell_TRA_clonotype.ipynb**

Software versions:
scanpy v.1.5.1

Sections:

  + Load TCR Data 
  + Add TRAV data to adata
  + Add TRBV data to adata
  + Add Clonotype data to adata
  + Gini Coefficient Analysis
  + Barplot paired TCR coverage
  + Filter cells based on TRAV
  + Calculate double TRA frequency
  + Filter cells based on TRAV (cont.)
  + Plot TRAV abundance by subgroup
  + TRAV diff between SKG High and SKG Low
  + Analysis for T.4 Nr4a1 hi cluster
  

**8_SKG_RA_single_cell_TRBV.ipynb**

Software versions:
scanpy v.1.5.1

Sections:

  + Filter based on TRBVs and dual TRAs 
  + Plot TRBV abundance by subgroup
  + TRBV diff between SKG High and SKG Low
  + Barplots for TRBV frequencies by subgroup
  + Analysis for Protein Data
  + Analysis for T.4N Nr4a1 cluster
  

**9_SKG_RA_single_cell_MAST.ipynb**

Software versions:
python: 
rpy2 v.3.3.2
anndata2ri v.1.0.4

R: 
MAST v.1.12.0


Sections:

  + Remove Dual TRA cells
  + Subset adata to SKG TRBV3 cells
  + Set up adata object for MAST
  + Run MAST for TRBV3
  + Save results
  + Volcano plot for TRBV3
  + Setup TRBV19 adata object for MAST
  + Run MAST for TRBV19
  + Volcano plot for TRBV19
  
  
### 6. Other Software Versions

Python3

R v.3.5.1

Jupiter Notebook version info:  
jupyter-notebook : 6.0.3  
ipykernel        : 5.1.4  
jupyter lab      : 1.2.6  

