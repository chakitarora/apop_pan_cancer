# apop_pan_cancer
research project

WorkFlow
1. TCGA data i.e <cancer>_final.csv to <cancer>_HRmedian_raw.csv : PRR_total_matrix.R
  note- used raw tcga data for this
  
2. <cancer>_HRmedian_raw.csv to <cancer>_matrix : PRR_cancerwise.ipynb
  note - provided a genelist for apoptotic_genes
  
3.<cancer>_matrix to Apop_pancan.csv : done manually.
  
4. Apop_pancan to fav and unfav csvs : pancancer_survival.ipynb
 
5. from fav unfav to 01 type csvs : done manually

6. used pan_can_gpm_bpm_common.ipynb : for pairwise bpm/gpm geneset intersection analysis.

7. apop_pan_can_sim_intersect.ipynb : finds intersecting genes between two cancers, also lists HR etc for those genes
   cancan_cutoffwise.ipynb : lists no of gpm/bpm genes in 33 cancers at different HR cutoffs.
   dendro.ipynb : for creating dendrograms

8. similarity_matrix_pancan.ipynb : code to evaluate 8 diff similarity coefficients (33x33 matrices, density plots and single column pairwise files)