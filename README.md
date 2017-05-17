# "Signature in blood and muscle of mdx mice on a transcriptomics level"
Bafstu Internship at LUMC, Human Genetics, BioSemantics Group

The used scripts are mostly written in R which is a statistical programming language which is available to download at https://cran.r-project.org/bin/windows/base/. 
An open-source program named Rstudio(https://www.rstudio.com/) has been used to create and run the created scripts. 

For performing the normalization on the data, the following script has been used:
  
    Normalization_MDX_Data.R
    
For performing the global test, the following script has been used:
 
      Global_Test.R
  
  
Thereafter, WGCNA was performed on the normalized data. This process took 3 steps which are performed in 3 different scripts:
  
    WGCNA_Preprocessing_1.R : for the preprocessing of the used dataset
    WGCNA_Modules_Eigengenes_2.R: creation of network and modules of the used dataset
    WGCNA_Retrieval_Genelist_Modules_3.R: retrieval of the genelist of the significant modules and correlation between modules and traits
  

