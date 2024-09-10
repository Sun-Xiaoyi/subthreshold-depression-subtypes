# subthreshold-depression-subtypes
This repository provides core code and relevant toolboxes for data analysis in the article "Connectome-based growth models reveal individual heterogeneity and neurophysiological subtypes of subthreshold depression".

## Overview
Contents include standalone software, source code, and demo data. All data necessary to reproduce our results have been publicly available, including age, sex, functional connectivity strengths (FCS), intermediate results during the analysis, and data for visualizing the figures.

## Installation Guide
All the scripts and the toolboxes can be executed by adding the appropriate environment paths according to the corresponding code language. Use the "add path" function for Matlab, the "pip install" function for Python, and the "install.packages()" function for R. These procedures are not time-consuming.

## Code and Data
1. Data preprocessing<br>
- SPM 12 (www.fil.ion.ucl.ac.uk/spm/)<br>
- SeeCAT (www.nitrc.org/projects/seecat)<br>
- Matlab 2020a (https://www.mathworks.com/products/matlab.html)<br>
2. Estimate the normative models of FCS on the healthy participants and quantify individual deviations for participants with subthreshold depression (StD)<br>
- step01_normative_deviations.py<br>
- PCNtoolkit 0.28 (https://github.com/amarquand/PCNtoolkit)<br>
- Python 3.8.18 (https://www.python.org)
- data/data_forNorm/age_sex_HC_Final.txt<br>
   Age and sex of all healthy participants<br>
- data/data_forNorm/age_sex_StD_Final.txt<br>
   Age and sex of all StD participants<br>
- data/data_forNorm/zFCS_HC_forNorm.txt<br>
   FCS values of all healthy participants<br>
- data/data_forNorm/zFCS_StD_forNorm.txt<br>
   FCS values of all StD participants<br>
3. Heterogeneous patterns of individual deviations in the functional Connectome in StD participants<br>
- step02_overlap_map.m<br>
- BrainNet Viewer 1.62 (https://www.nitrc.org/projects/bnv) <br>
- data/data_results/Z_estimate.txt<br>
 Deviation values of StD participants<br>
- data/data_results/z_base.txt<br>
 Baseline deviation values of StD participants<br>
- data/data_results/z_treatment.mat<br>
 Baseline and follow-up deviation values of StD participants who underwent 8 weeks of bright light therapy (BLT)<br>
- data/data_results/pro_dev_pos.nii<br>
Spatial overlap map of extreme positive deviations<br>
- data/data_results/pro_dev_pos.nii<br>
Spatial overlap map of extreme negative deviations<br>
4. Identify StD subtypes based on the individual deviation patterns<br>
- step03_subtype_NbClust.R<br>
- NbClust 3.0 (https://www.rdocumentation.org/packages/NbClust)<br>
- R 3.6.1 (https://www.r-project.org)<br>
- step04_subtype_results.m<br>
- data/data_results/clus_base.mat<br>
Subtyping results of StD participants<br>
- data/data_results/clus_treatment.mat<br>
Subtyping results of StD participants who underwent 8 weeks of BLT<br>
5. Subtype differences in deviation patterns<br>
- step05_subtype_brain_base.m<br>
- data/data_results/z_clus1.nii<br>
Mean deviation map of subtype 1<br>
- data/data_results/z_clus2.nii<br>
Mean deviation map of subtype 2<br>
6. Subtype differences in clinical manifestations<br>
- step06_subtype_clinic_base.m<br>
- IBM SPSS Statistics 20 (https://www.ibm.com/spss)
7. Subtype differences in gene expression profiles of FCS deviations<br>
- step07_1_gene_expression.py<br>
- abagen 0.1.3 (https://github.com/rmarkello/abagen) <br>
- step07_2_z_surrogate.py<br>
- step07_3_corr_gene_deviation.m<br>
- Metascape (https://metascape.org/gp/index.html#/main/step1) <br>
- Due to size limitations, the script-relevant files and data can be found at https://pan.bnu.edu.cn/l/w1fD5X <br>
8. Subtype differences in treatment response to BLT<br>
- step08_1_subtype_clinic_follow.m<br>
- step08_2_subtype_clinic_follow.R<br>
- step09_1_subtype_brain_follow.m<br>
- step09_2_subtype_brain_follow.R<br>
- step10_subtype_brain_clinic_follow.m<br>
- LIBSVM 3.3 (https://www.csie.ntu.edu.tw/~cjlin/libsvm/)<br>
- Pattern_Regression (https://github.com/ZaixuCui/Pattern_Regression_Clean) <br>
- data/data_results/brain_F_allFDR.nii<br>
The brain regions with significant subtype-by-treatment interaction effect<br>

## Questions, suggestions and improvements
Email: xy.xiaoyisun@gmail.com
