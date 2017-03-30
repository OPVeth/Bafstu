#!/usr/bin/Rscript
# This script has been created for the normalization of the MDX mouse dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# Normalization and Log Transformation - 31-10-2016

source("http://bioconductor.org/biocLite.R")
library(edgeR)
library(matrixStats)
setwd("~/Dropbox/Bafstu/Normalization")
options(stringsAsFactors = FALSE)

preProcessingData <- function(data){
    # The following parameters are loaded in:
    #     data: The CSV file containing the data
    # First, the genes are retrieved from the first column of the dataset.
    # Afterwards, the columns with the header containing 'feature' are removed, 
    # because they only contain gene symbols and no actual count data. 
    # The table is transposed and the rows get the names of the Gene Symbols. 
    # The preprocessed data gets returned. 
    genes <- data[, 1]
    data <- data[, which(colnames(data) != c("feature", "feature.1"))]
    data <- as.data.frame(t(data))
    names(data) <- genes
    return(data)
}

performNormalization <- function(miceData, traitDataTissue, tissue){
    # The following parameters are loaded in:
    #     miceData: preprocessed expression data
    #     traitDataTissue: The CSV file containing the trait data
    #     tissue: a String object containing the name of the tissue
    # The function getTrait is called and this results with a specific selection
    # of the original trait data file. 
    # Thus, containing only the relevant information about the tissue in question
    # and selected genotypes and age. 
    # The samples of the trait data and dataset itself are compared, and the
    # resulting dataset is given to the function getNormalizedData() with tissue. 
    # The resulting matrix, containing the normalized data gets returned. 
    traitData <- getTrait(traitDataTissue, tissue, miceData)
    exprData <- miceData[which(rownames(miceData) %in% traitData$deCODEbarcode),]
    exprDataNormal <- getNormalizedData(exprData, tissue)
    return(exprDataNormal)
}

getTrait <- function(traitDataTissue, tissue, miceData){
    # The following parameters are loaded in:
    #     traitDataTissue: The CSV file containing the trait data
    #     tissue: a String object containing the name of the tissue
    #     miceData: preprocessed expression data
    # Two variables are created, age and genotypes, and these two are freely
    # changeable to will. Age is an integer which indicates which week you
    # want to study and genotypes indicates the genotypes which you're
    # interested in. 
    # Thereafter, the genotypes in the strain column of the traitData get rid
    # of the '/utrn', thus changing the format. 
    # This results in 'mdx+-' rather then mdx/utrn+/-. 
    # The rows which age and genotypes matched with the chosen
    # variable value are retrieved and afterwards compared to each other.
    # This results in rows which match both qualifications of age and genotype.
    # The rows of the resulting trait dataset which are also present in 
    # miceData get retrieved and are saved and returned afterwards
    age <- 30
    genotypes <- list("WT", "mdx", "mdx++", "mdx+-")
    traitDataTissue$strain1 <- gsub('utrn', '', gsub('/', '', traitDataTissue$strain1))
    agePos <- grep(age, traitDataTissue$age)
    genPos <- which(traitDataTissue$strain1 %in% genotypes)
    traitData <- traitDataTissue[intersect(agePos, genPos), ]
    genotypes <- traitData[which(traitData$deCODEbarcode %in% rownames(miceData)),]$strain1
    save(genotypes, file=paste("genotypes_", tissue, "WT.RData", sep=""))
    return(traitData)
}

getNormalizedData <- function(dataTissue, tissue){
    # The following parameters are loaded in:
    #     dataTissue: The CSV file containing the trait data
    #     tissue: a String object containing the name of the tissue
    # The RData file containing the genotypes of a specific tissue is 
    # loaded in.
    # Afterwards, a DGE object is created with the transposed dataset and
    # genotypes as 'group'. 
    # A filtering is afterwards used: if the count is lower then
    # 10 for 40% of the samples, then the gene is removed from the dataset. 
    # The resulting dataset gets normalized with calcNormFactors() and 
    # the parameter 'TMM'.
    # The TMM normalized dataset is saved and the final normalized is 
    # retrieved with cpm().
    # The final normalised dataset is returned. 
    load(paste("genotypes_", tissue, "WT.RData", sep=""))
    dge <- DGEList(t(data.matrix(dataTissue)), group=genotypes)
    keep <- rowSums(dge$counts > 10) >= nrow(dataTissue) * 0.4
    dge <- dge[keep, keep.lib.sizes=FALSE]
    tmm <- calcNormFactors(dge, method= "TMM")
    save(tmm, file=paste("Gene_names_", tissue, "WT.RData", sep=""))
    normalizedData <- cpm(tmm)
    return(normalizedData)
}

createBoxPlots <- function(tissueDataNormal, dataTissue, tissue){
    # The following parameters are loaded in:
    #     tissueDataNormal: A data frame of the normalized count data
    #     dataTissue: The CSV file containing the trait data
    #     tissue: a String object containing the name of the tissue
    # The RData file containing the genotypes of the tissue is loaded in. 
    # Afterwards, two boxplots are created, consisting of the normalized
    # and not normalized dataset. 
    # The result gets saved in a file. 
    load(paste("genotypes_", tissue, ".RData", sep=""))
    pdf(paste("Boxplot_Normalization", tissue, ".pdf", sep=""))
    par(mfrow = c(1, 2))
    boxplot(log2(tissueDataNormal+1), main=paste("TMM - Normalized",tissue," Data"), las=2)
    boxplot(log2(t(data.matrix(dataTissue))+1), main=paste("Not-normalized ",tissue," Data"), las=2)
    dev.off()
}

saveNormalData <- function(tissueDataNormal, tissue){
  # The following parameters are loaded in:
  #     tissueDataNormal: A data frame of the normalized count data
  #     tissue: a String object containing the name of the tissue
  # Two files are created, the first one containing the normalized
  # dataset and the second the log transformed and normalized dataset.
  write.csv2(tissueDataNormal, paste("Count_",tissue,"_Normalized.csv", sep="")) #Normalized data
  write.csv2(log2(tissueDataNormal + 1), paste("Count_",tissue,"_Normalized_Log_Transformed.csv", sep=""))
  logtissueDataNormal <- log2(tissueDataNormal + 1)
  save(logtissueDataNormal, file=paste("Count_",tissue,"_Normalized_Log_Transformed.RData", sep=""))
}

main <- function(){
    # The count data file is read in and afterwards, a vector is created
    # containing the two types of tissue as strings.
    # For every tissue in the list, the matching trait data is firstly
    # read in. 
    # Thereafter, the function performNormalization() is given
    # the original data, the data containing the traits and the tissue 
    # as a string. This results with a data frame containing the 
    # normalized data. 
    # Afterwards, createBoxplots() gets the same parameters as the 
    # function before. 
    # Finally the normalized data is saved with the function saveNormalData.
    miceData <- preProcessingData(read.csv2("~/Dropbox/Bafstu/batch1_batch2_all.counts.tsv", sep=" "))
    tissues <- c("Muscle", "Blood")
    for (x in 1:length(tissues)){
      if (tissues[x] == "Blood"){
        traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[1:96,]
      } else{
        traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[97:147,]
      }
      print("Normalizing...")
      normalizedData <- performNormalization(miceData, traitData, tissues[x])
      print("Creating Box plots...")
      createBoxPlots(normalizedData, traitData, tissues[x])
      print("Saving normalized Data")
      saveNormalData(normalizedData, tissues[x]) #Blood
    }
}
main()
# Session Information
# R version 3.3.2 (2016-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 18
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# other attached packages:
#  [1] matrixStats_0.51.0   edgeR_3.12.1         limma_3.26.9         BiocInstaller_1.20.3 globaltest_5.24.0    survival_2.40-1      tidyr_0.6.0          pander_0.6.0         annotables_0.1.1     dplyr_0.5.0         
#  [11] zoo_1.7-13           devtools_1.12.0      gtools_3.5.0         rWikiPathways_0.0.1  urltools_1.6.0       RCurl_1.95-4.8       bitops_1.0-6         jsonlite_1.1         plyr_1.8.4           curl_2.2            
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.7          tools_3.3.2          digest_0.6.10        annotate_1.48.0      memoise_1.0.0        RSQLite_1.0.0        tibble_1.2           lattice_0.20-34      Matrix_1.2-7.1       DBI_0.5-1           
#  [11] parallel_3.3.2       withr_1.0.2          S4Vectors_0.8.11     IRanges_2.4.8        stats4_3.3.2         triebeard_0.3.0      grid_3.3.2           Biobase_2.30.0       R6_2.2.0             AnnotationDbi_1.32.3
#  [21] XML_3.98-1.5         magrittr_1.5         splines_3.3.2        BiocGenerics_0.16.1  assertthat_0.1       xtable_1.8-2        