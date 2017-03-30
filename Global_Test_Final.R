#!/usr/bin/Rscript
# This script has been created for performing the global test on the MDX mouse dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# Global Testing with Pathway Subsets - 30-03-2016

library(curl)
library(plyr)
library(jsonlite)
library(rWikiPathways)
library(gtools)
library(jsonlite)
library(devtools)
library(zoo)
library(dplyr)
library(annotables)
library(pander)
library(tidyr)
library(globaltest)
setwd("~/Dropbox/Bafstu/Global Test")

preProcessingData <- function(dataTissue, tissue){
  # The following parameters loaded in:
  #     dataTissue: a data frame containing the count data
  #     tissue: a String object
  # The rows of the first column is used as the rownames of
  # dataTissue.
  # Afterwards, the first column is removed and transposed. 
  # Thereafter, a Rdata is loaded in which contains the genotypes
  # corresponding to the dataset.
  # Afterwards, the function 'createMap()' is called. 
  # Thereafter, a CSV file is read in which contains gene symbols
  # and their corrresponding Ensembl ID's. 
  # This file is given to the function 'getSymbols' with 'dataTissue'. 
  # The resulting data frame gets an extra column named 'Group'
  # containing the genotypes. 
  # The rownames of the data frame get removed and afterwards, the frame
  # gets returned.
  genes <- dataTissue$X
  dataTissue <- dataTissue[, 2:ncol(dataTissue)]
  rownames(dataTissue) <- genes
  dataTissue <- as.data.frame(t(dataTissue))
  load(paste("~/Dropbox/Bafstu/Normalization/genotypes_", tissue,".RData", sep=""))
  createMap()
  mapData <- read.csv2("Ensemble_Map.csv", sep=";")
  dataTissue <- dataTissue[,which(names(dataTissue) %in% mapData$SYMBOL)]
  geneSymbols <- colnames(dataTissue)
  dataTissue$Group <- genotypes
  rownames(dataTissue) <- NULL
  return(dataTissue)
}


createMap <- function(){
  # A data frame is created wth two columns containing all gene Symbols
  # and their ENSEMBL ID's. 
  # The column names are changed into 'ENSEMBL' and'SYMBOL' and the frame
  # is written to a CSV file 
  mapData <- data.frame(grcm38$ensgene, grcm38$symbol) ##Mouse
  colnames(mapData) <- c("ENSEMBLE", "SYMBOL")
  write.table(mapData, file="Ensemble_Map.csv", sep=";", row.names=F)
}

createPathwayMappingFile <- function(mousePathways){
  # The following parameter will be used:
  #     mousePathways: data frame containing genes and pathways
  # Two variables are created and both set on NULL. 
  #   With a for-loop a double list is created with the
  #   pathway as the key and the genes with Esembl ID's as its values. 
  #     The second for-loop checks if the gene is a NA value, if not, the gene ID is matched
  #     with a symbol. If there are multiple matches, the first one is used. The resulting
  #     Symbol is afterwards used as the gene ID.
  #   If the loop runs for the first time, zooMatrixGenes is equal to a
  #   column containing genes of the pathway. 
  #   If that's not the case, zooGenes is equalized to this column
  #   and is added to zooMatrixGenes.
  # After the loop, zooMatrixGenes is converted to a data frame and its
  # columns get the name of the corresponding pathway. 
  # Thereafter, the data frame gets returned. 
  listGenes <- NULL
  zooMatrixGenes <- NULL
  mapData <- read.csv2("Ensemble_Map.csv", sep=";")
  for (x in 1:length(mousePathways)) {
    listGenes[[mousePathways[x]]] <- getXrefList(pathway = mousePathways[x], systemCode="En") #Changeable
    genes <- listGenes[[mousePathways[x]]]
    for (y in 1:length(genes)){
      if (!(is.na(genes[y]))){
         matching <- which(mapData$ENSEMBLE %in% genes[y])
         if (length(matching) > 1){
              genes[y] <- as.character(mapData$SYMBOL[matching[1]])
          } else if (length(matching) == 1){
              genes[y] <- as.character(mapData$SYMBOL[matching])
          } 
      }
    }
    if (x == 1){
      zooMatrixGenes <- zoo(do.call(cbind, list(genes)))
    } else{
      zooGenes <- zoo(do.call(cbind, list(genes)))
      zooMatrixGenes <- merge(zooMatrixGenes, zooGenes)
    }
  }
  dataPathwayGenes <- data.frame(zooMatrixGenes)
  colnames(dataPathwayGenes) <- paste(mousePathways, listPathways(organism="Mus musculus")$name)
  write.csv2(dataPathwayGenes, "Pathway_Gene_Symbol.csv", row.names=F)
  return(dataPathwayGenes)
}

## From: Ekrem Sabir's code: GT_Mouse_meta_WestfallYoung_final.R
getPathways <- function(dataTissue, dataPathwayGenes){
  # The following parameters are loaded in:
  #     dataTissue: a data frame containing count data
  #     dataPathwayGenes: a data frame containing the pathway and mathing symbols
  # Two empty vectors are defined and the levels of dataPathwayGenes are removed.
  #   Thereafter, a for-loop checks if any of the rows in a column contains NA values. If that's not the case, the genes in dataPathwayGenes
  #   are matched with dataTissue column names. 
  #   If there's a match, the pathway gets the name of the pathway as key
  #   and the matching genes as values.
  # All NaN values get removed from the data frame. Thereafter, a filtering
  # is used which removes pathways with
  # less then 2 matching genes. 
  # The resulting data frame is returned. 
  pathways <- c()
  pathwaysList <- c()
  dataPathwayGenes <- droplevels(dataPathwayGenes)
  for(y in 1:ncol(dataPathwayGenes)){
    for(x in 1:length(dataPathwayGenes[y])){
      if (!(is.na(dataPathwayGenes[[y]][x]))){
        matching <- which(colnames(dataTissue) %in% dataPathwayGenes[[y]])
        if (length(matching) > 0){
          pathwaysList[[colnames(dataPathwayGenes[y])]] <- colnames(dataTissue)[matching]
        }
        else{
          pathwaysList[[colnames(dataPathwayGenes[y])]] <- NA
        }
      }
    }
  }
  pathwaysList <- lapply(pathwaysList, function(x) x[!is.na(x)])
  pathwaysList <- Filter(function(x) length(x) > 2, pathwaysList)
  return(pathwaysList)
}

doGT <- function(genotype1, genotype2, dataTissue, tissue){
  # The following parameters are loaded in:
  #     genotype1: a String object
  #     genotype2: a String object
  #     dataTissue: a data frame containing count data
  #     tissue: a String object
  # A list is created containing genotype1 and genotype2. 
  # Thereafter, it is checked which rows of the 'Group' column of
  # dataTissue matches with the genotypes in the created list. 
  # Only the matching rows are retrieved and given to gt() with all
  # columns of dataTraits except the last one. 
  # The resulting GT statistics are returned.
  listGenotypes <- list(genotype1, genotype2)
  dataTissue <-  as.vector(dataTissue[which(dataTissue$Group %in% listGenotypes),])
  gtResults <- gt(as.factor(dataTissue$Group), dataTissue[, 1:ncol(dataTissue)-1], permutations=10000)
  write.table(result(gtResults), file = paste(listGenotypes[1], "_vs_", listGenotypes[2],"_",tissue,".txt", sep=""), sep = "\t")
  return(gtResults)
}

performGT <- function(dataTissue, tissue){
  # The following parameter is loaded in:
  #     dataTissue: a data frame containing count data
  #     tissue: a String object
  # The function doGT() is called with dataTissue and different 
  # genotypes and prints every result. 
  print(doGT("WT", "mdx", dataTissue, tissue))
  print(doGT("WT", "mdx+-", dataTissue, tissue))
  print(doGT("WT", "mdx++", dataTissue, tissue))
  print(doGT("mdx", "mdx+-", dataTissue, tissue))
  print(doGT("mdx", "mdx++", dataTissue, tissue))
  print(doGT("mdx+-", "mdx++", dataTissue, tissue))
  #getStatistics(doGT("WT", "mdx", dataTissue))
}

performGTPathway <- function(pathwaysTissue, dataTissue, tissue){
  # The following parameters are loaded in:
  #     pathwaysTissue: a data frame containing the pathways
  #     dataTissue: a data frame containing count data
  #     tissue: a String object
  # The function doGTPathway() is called with dataTissue, 
  # pathwaysTissue, different genotypes and prints every result.)
  #subsets <- read.table("Pathway_subsets_multinomial.txt", sep = "\t")
  #genotypes <- c("WT", "mdx", "mdx+-", "mdx++")
  #for (x in length(pathwaysTissue)){
    doGTPathway("WT", "mdx", pathwaysTissue, dataTissue, tissue)
    doGTPathway("WT", "mdx+-", pathwaysTissue, dataTissue, tissue)
    doGTPathway("WT", "mdx++", pathwaysTissue, dataTissue, tissue)
    result(doGTPathway("mdx", "mdx+-", pathwaysTissue, dataTissue, tissue))
    result(doGTPathway("mdx", "mdx++", pathwaysTissue, dataTissue, tissue))
    result(doGTPathway("mdx+-", "mdx++", pathwaysTissue, dataTissue, tissue))
    #getStatistics(doGT_muscle("WT", "mdx"))
   # }
}

getStatistics <- function(gtValue){
  # The following parameter is loaded in:
  #     gtValue: a Global Test object
  # The features are retrieved of the object and returned.
  return(features(gtValue), plot=F)
}

doGTPathway <- function(genotype1, genotype2, pathwaysTissue, dataTissue, tissue){
  # The following parameters are loaded in:
  #     genotype1: a String object
  #     genotype2: a String object
  #     pathwaysTissue: a data frame containing the pathways
  #     dataTissue: a data frame containing count data
  #     tissue: a String object
  # A list is created containing genotype1 and genotype2. 
  # Thereafter, it is checked which rows of the 'Group' column of
  # dataTissue matches with the genotypes in the created 
  # list and given to gt() with all columns of dataTraits except the last
  # one with each genotype comparison, for every pathway in pathwaysTissue and with a
  # 'multinomial' based method. 
  # If the p-value < 0.05 after a Holm correction, each genotype comparison,
  # pathwaysTissue, the posiion and dataTissue is given to saveResultsGT().
  # Thereafter a regular GT() is performed on every pathway. A Holm correction is performed on the global test data and saved 
  # in a .CSV file. The results are returned afterwards.
  listGenotypes <- list(genotype1, genotype2)
  positions <- which(dataTissue$Group %in% listGenotypes)
  genotypes <- c("WT", "mdx", "mdx++", "mdx+-")
  exactRows <- dataTissue[positions,]$Group
  for (x in 1:length(pathwaysTissue)){
    gtResults <- gt(as.factor(exactRows), dataTissue[,1:ncol(dataTissue)-1], permutations=10000, subsets = pathwaysTissue[[x]], model = "multinomial")
    if (p.adjust(p.value(gtResults)) < 0.05){
      saveResultsGT(genotypes[1], genotypes[2], pathwaysTissue, x, dataTissue)
      saveResultsGT(genotypes[1], genotypes[3], pathwaysTissue, x, dataTissue)
      saveResultsGT(genotypes[1], genotypes[4], pathwaysTissue, x, dataTissue)
      saveResultsGT(genotypes[2], genotypes[3], pathwaysTissue, x, dataTissue)
      saveResultsGT(genotypes[2], genotypes[4], pathwaysTissue, x, dataTissue)
      saveResultsGT(genotypes[3], genotypes[4], pathwaysTissue, x, dataTissue)
    }
  }
  gtResults <- gt(as.factor(exactRows), dataTissue[positions, 1:ncol(dataTissue)-1], permutations=10000, subsets = pathwaysTissue[[pathways]])
  gtResults_Holm <- p.adjust(gtResults, method= "holm")
  write.table(result(gtResults_Holm), file=paste(genotype1, "vs", genotype2, "_", tissue, ".csv"))
  covariates(gtResults_Holm, legend= c(paste("associated with ", genotype1), paste("associated with ", genotype2)), pdf=paste(listGenotypes[1],"_vs_",listGenotypes[2],tissue, ".pdf"))
  return(gtResults)
}

saveResultsGT <- function(genotype1, genotype2, pathway, pos, dataTissue){
  # The following parameters are read in:
  #     genotype1: a String object
  #     genotype2: a String object
  #     pathway: a list
  #     pos: an Integer
  #     dataTissue: a data frame containing count data
  # A list is created of both Strings and the indexes which matches with the list
  # within dataTissue are retrieved.
  # The rows with the right genotypes of dataTissue are retrieved.
  # Thereafter gt() is performed with exactRowsSecond, the right rows and every column except the last of dataTissue and
  # one pathway with the method 'logistic'. 
  # If the resulting p-value <0.05 the result is saved in a .TXT file.
  listGenotypes <- list(genotype1, genotype2)
  positions <- which(dataTissue$Group %in% listGenotypes)
  exactRowsSecond <- dataTissue[positions,]$Group
  gtObject <- gt(as.factor(exactRowsSecond), dataTissue[positions,1:ncol(dataTissue)-1], permutations=10000, subsets = pathway[[pos]], model="logistic")
  if (p.value(gtObject) < 0.05){
   write.table(paste(genotype1, genotype2, p.value(gtObject), names(pathway)[pos], sep=";"), file="Results_Multinomial_Full.txt", append = T, col.names = F, row.names = F)
  }
}


main <- function(){
  # A vector containing two String object is created and all the
  # ID's of the pathways of the Mus musculus are retrieved.
  # The list of genes is then given to the function createPathwayMappingFile.
  # Thereafter, a CSV file containing the genes and pathway data is read in.
  #   For each tissue, a specific CSV file is read in containing
  #   the trait data. 
  #   Thereafter, the normalized expression data of the tissue is read in
  #   and used for the function 'preProcessingData' with the tissue as a 
  #   String object. 
  #   The resulting preprocessed data is saved in a variable. 
  #   The function 'performGT' is called with dataTissue, resulting in print
  #   statements of the global test results.
  #   Thereafter, the matching pathways are retrieved with the function 'getPathways'
  #   with dataTissue and dataPathwayGenes.
  #   Afterwards, the function 'performGTPathway' is called with pathwaysTissue
  #   and dataTissue. Thus, performing Global Test with the pathways as subsets. 
  tissues <- c("Blood")
  mousePathways <- listPathways(organism="Mus musculus")$id 
  createPathwayMappingFile(mousePathways)
  dataPathwayGenes <- read.csv2("Pathway_Gene_Symbol.csv")
  for(x in 1:length(tissues)){
    if (tissues[x] == "Blood"){
      traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[1:96,]
    } else{
      traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[97:147,]
    }
    dataTissue <- data.frame(read.csv2(paste("~/Dropbox/Bafstu/Normalization/Count_",tissues[x], "_Normalized_Log_Transformed.csv", sep="")))
    print("Data preprocessing...")
    dataTissue <- preProcessingData(dataTissue, tissues[x])
    write.table(dataTissue, file=paste("counts_data_", tissues[x],".txt", sep=""))
    dataTissue <- read.table(paste("counts_data_", tissues[x],".txt", sep=""), header=TRUE, check.names=FALSE)
    attach(dataTissue)
    print("Performing global test...")
    performGT(dataTissue, tissues[x])
    print("Performing global test with pathways...")
    pathwaysTissue <- getPathways(dataTissue, dataPathwayGenes)
    performGTPathway(pathwaysTissue, dataTissue, tissues[x])
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
