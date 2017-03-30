# Bafstu
Bafstu Internship at LUMC, Human Genetics, BioSemantics Group

# "Signature in blood and muscle of mdx mice on a transcriptomics level"
The used scripts are mostly written in R which is a statistical programming language which is available to download at https://cran.r-project.org/bin/windows/base/. An open-source program named Rstudio(https://www.rstudio.com/) has been used to create and run the created scripts. 

For performing the normalization on the data, the following script was used:
  Normalization_MDX_Data_Final.R
  
Thereafter, WGCNA was performed on the normalized data. This process took 3 steps which are performed in 3 different scripts:
  WGCNA_Preprocessing_1_Final.R : for the preprocessing of the used dataset
  WGCNA_Modules_Eigengenes_2_Edited.R: creation of network and modules of the used dataset
  WGCNA_Retrieval_Genelist_Modules_3_Final.R: retrieval of the genelist of the significant modules and correlation between modules and traits
  
  
 For performing the global test, the following script has been used:
  Global_Test_Final.R
  
