# Kanoski-Lab-Vagus-HPCd-ACh-manuscript-code
This is MATLAB_R2022a code pertaining to the analysis of fiber photometry data conducted for the Vagus HPCd ACh manuscript

Installation time for MATLAB_R2022a with all statistical, linear algebra, and plotting packages on a standard i7 intel CPU chip with 8GB of RAM is approximately ~30-45 minutes. 

Neurophotometrics_LTL.m is the MATLAB script that performs the primary anlysis of raw photometry files, it performs the isosbestic signal subtraction and correction. 

Refeeding_meal_bout_analysis.m is the script that performs the quantification of changes in ACh release within active eating bouts interbout intervals and extracts data to perform AUC quantification on pre consumption and post consumption periods. 

refeeeding_revisions_analyses.m performs all analyses of HPCd ACh release timing around behavioral switchpoints included in Supplementary figure 2.

Photometry_SFS.m performs all analyses for analyzing HPCd ACh release during the food location working memory task, including the mixed-effects model for ACh-Speed coupling. 
