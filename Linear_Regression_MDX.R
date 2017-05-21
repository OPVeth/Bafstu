#!/usr/bin/Rscript
# This script has been created for the linear regression analysis of the dystrophic modules of the MDX mouse dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# Normalization and Log Transformation - 2-05-2017

setwd("")
library(edgeR)
library(limma)

readData <- function(module, tissue){
  # The following parameters are loaded in:
  #     module: a String object
  #     tissue: a String object
  # The hubgenes symbols are retrieved of a specific 
  # tissue and module and thereafter returned.
  data <- read.table(paste("Module_", module, "_", tissue, "_HubGenelist.txt", sep= ""))
  return(data)
}

getTraitData <- function(datTraits){
  # The following parameters are loaded in:
  #     traitData: a data frame containing sample information
  # The mdx column is edited from ntegers to strings indicating
  # the corresponding genotype.
  # Thereafter the rows with the genotype 'WT' and' mdx' are
  # retrieved and returned.
  for (x in 1:length(datTraits$mdx)){
    if (datTraits$mdx[x] == 0){
      datTraits$mdx[x] <- "WT"
    } else if (datTraits$mdx[x] == 1){
      datTraits$mdx[x] <- "mdx"
    } else if (datTraits$mdx[x] == 2){
      datTraits$mdx[x] <- "mdx++"
    } else if (datTraits$mdx[x] == 3){
      datTraits$mdx[x] <- "mdx+-"
    } else {
      datTraits$mdx[x] <- ""
    }
  }
  return(datTraits[which(datTraits$mdx %in% c("WT", "mdx")),])
}

main <- function(){
  # A vector is created containing the two types of tissue as strings.
  # Depending on the tissue, a vector is created with specific Strings
  # indicating the dystrophic module names.
  #   The RData file of the WGCNA preprocessing script with all genotypes is read in
  #   and thereafter, count data file is read in. 
  #     Preprocessing of the count data is performed before the retrieval of the
  #     count data of each hubgene from the file of the specific tissue and module.
  #     Thereafter, a design is created based on the count data and mdx column.
  #     After fitting a model and performing Bayes, the p-value regarding the
  #     linear regression is printed.
  tissues <- c("Blood", "Muscle")
  for (x in 1:length(tissues)){
    if (tissues[x] == "Blood"){
      modules <- c("green", "magenta", "yellow", "lightcyan", "brown", "blue")
    } else {
      modules <- c("darkgreen", "magenta", "tan", "lightcyan", "salmon", "yellow", "cyan", "midnightblue")
    }
    load(paste("Mice-01-dataInput-", tissues[x], ".RData", sep=""))
    datExpr <- read.csv2(paste("Count_", tissues[x], "_Normalized.csv", sep=""))
    rownames(datExpr) <- datExpr$X
    datExpr <- datExpr[, 2:ncol(datExpr)]
    datExpr <- as.data.frame(t(datExpr))
    datTraits <- getTraitData(datTraits)
    for (y in 1:length(modules)){
      data <- readData(modules[y], tissues[x])
      countData <- datExpr[which(rownames(datExpr) %in% rownames(datTraits)), which(colnames(datExpr) %in% data[,1])]
      design <- model.matrix(~ datTraits$mdx, data=countData)
      lm_test <- lmFit(t(data.matrix(countData)), design)
      lm_bayes <- eBayes(lm_test, trend=F)
      print(topTable(lm_bayes))
    }
  }
}
main()
# session Information
# R version 3.3.2 (2016-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 18
#
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=nl_NL.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
#
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
#  [1] annotables_0.1.1      magrittr_1.5          zoo_1.7-14            WGCNA_1.51            RSQLite_1.1-2         fastcluster_1.1.22    dynamicTreeCut_1.63-1 matrixStats_0.52.0   
#  [9] BiocInstaller_1.20.3  edgeR_3.12.1          limma_3.26.9         
#
# loaded via a namespace (and not attached):
# [1] splines_3.3.2         lattice_0.20-35       colorspace_1.3-2      htmltools_0.3.5       stats4_3.3.2          base64enc_0.1-3       survival_2.41-2       foreign_0.8-67       
# [9] DBI_0.6-1             BiocGenerics_0.16.1   RColorBrewer_1.1-2    foreach_1.4.3         plyr_1.8.4            stringr_1.2.0         munsell_0.4.3         gtable_0.2.0         
# [17] htmlwidgets_0.8       codetools_0.2-15      memoise_1.0.0         latticeExtra_0.6-28   Biobase_2.30.0        knitr_1.15.1          IRanges_2.4.8         doParallel_1.0.10    
# [25] parallel_3.3.2        AnnotationDbi_1.32.3  htmlTable_1.9         preprocessCore_1.32.0 Rcpp_0.12.10          acepack_1.4.1         scales_0.4.1          backports_1.0.5      
# [33] checkmate_1.8.2       S4Vectors_0.8.11      Hmisc_4.0-2           gridExtra_2.2.1       impute_1.44.0         ggplot2_2.2.1         digest_0.6.12         stringi_1.1.3        
# [41] grid_3.3.2            tools_3.3.2           lazyeval_0.2.0        tibble_1.3.0          Formula_1.2-1         cluster_2.0.6         GO.db_3.2.2           Matrix_1.2-8         
# [49] data.table_1.10.4     iterators_1.0.8       rpart_4.1-10          nnet_7.3-12       