# "Signature in blood and muscle of mdx mice on a transcriptomics level"
Bafstu Internship at LUMC, Human Genetics, BioSemantics Group

The used scripts are mostly written in R which is a statistical programming language which is available to download at https://cran.r-project.org/bin/windows/base/. 
An open-source program named Rstudio(https://www.rstudio.com/) has been used to create and run the created scripts. 

The used scripts for the analysis are located in the directory' Scripts'.
The used transcriptomics mdx dataset is available at:

    https://www.dropbox.com/s/1tyjw2844s6eanj/batch1_batch2_all.counts.tsv?dl=0 

For performing the normalization on the data, the following script has been used (for performing only normalization on dystrophic genoypes, 'WT' should be removed in line 64):

    Normalization_MDX_Data_Final.R
    
The normalization script will produce multiple files which are needed for performing the global test and WGCNA. Those files are; 

    Count_<Tissue name>_Normalized.csv
    Count_<Tissue name>_Log_Transformed.csv
    genotypes_<Tissue name>.RData
    
For performing the global test, the following script has been used:

    Global_Test_Final.R
    
After running the script, the following files are produced:

    counts_data_<Tissue name>.txt
    Pathway_Gene_Symbol.csv # All pathways with their genes as symbols
    Ensemble_Map.csv # Mapping file to convert Ensembl ID's to Symbols
    <genotype1>_vs_<genotype2>_<Tissue name>.txt # Global test results
    <genotype1>_vs_<genotype2>_<Tissue name>.csv # Table format pathway enrichment results
    <genotype1>_vs_<genotype2>_<Tissue name>.pdf # Figure format pathway enrichment results
    
Thereafter, WGCNA was performed on the normalized data. The preprocessing - step 1 of WGCNA - is done with the following script.

    WGCNA_Preprocessing_1_Final.R : for the preprocessing of the used dataset
    
The first script will produce the following files:

    sampleClustering_<Tissue name>.pdf
    sampleDendogram_<Tissue name>.pdf
    Mice-01-dataInput-<Tissue name>.RData # File needed for step 2 and 3 of WGCNA
    WGCNA_Modules_Eigengenes_2_Edited.R: creation of network and modules of the used dataset
    WGCNA_Retrieval_Genelist_Modules_3_Final.R: retrieval of the genelist of the significant modules and correlation between modules and traits

After preprocessing the data, the network is created with the second script for WGCNA:

    WGCNA_Modules_Eigengenes_2b.R
The second script will produce the following files:

    Topology_Connectivity_<Tissue name>.pdf
    Dissimilarity_Matrix_<Tissue name>.pdf
    Gene_Dendogram_Module_Colors_<Tissue name>_<Chosen minModuleSize>.pdf
    Clustering_Module_Eigengenes_<Tissue name>_<Chosen minModuleSize>.pdf
    Dynamic_Cut_Merged_Dendogram_<Tissue name>_<Chosen minModuleSize>_<Chosen MEDissThreshold>.pdf
    Mice-02-networkConstruction_<Tissue name>.RData # Needed for the third step in WGCNA

In the final and third step of WGCNA, the genes and ME's get retrieved with the script:

    WGCNA_Retrieval_Genelist_Modules_3.R
    
The third script will produce the following files:

    Module_Trait_Relationship_<Tissue name>.pdf
    Module_Trait_Relationship_Genotypes_<Tissue name>.pdf
    Module_Trait_Relationship_MEs.pdf
    Module_Membership_Gene_Sinificance_<Tissue name>.pdf
    Entrez_Map.csv # Mappingfile to convert symbols to Entrez ID's
    Module_<Module name>_<Tissue name>_Genelist.txt
    Module_<Module name>_<Tissue name>_HubGenelist.txt
    Normalized_CountData_HubGenes_<module name>_<Tissue name>.csv
    
The Linear Regression is performed with the script:

    Linear_Regression_MDX.R
    
The results of the linear regression are printed out. 
