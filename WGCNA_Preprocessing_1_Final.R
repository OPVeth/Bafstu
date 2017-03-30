# This script has been created for performing the WGCNA on the MDX dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# WGCNA - Preprocessing - Script 1 - 14-02-2016

setwd("~/Dropbox/Bafstu/WGCNA")
library(WGCNA)
library(matrixStats)
library(zoo)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

preProcessingData <- function(miceData, tissue){
    # The following parameters are loaded in:
    #     miceData: The CSV file containing the data
    #     tissue: a String object
    # The rows of the first column is used as the rownames of miceData. 
    # Afterwards, the first column is removed and transposed. 
    # Thereafter, the data frame gets returned.
    rownames(miceData) <- miceData$X
    miceData <- miceData[, 2:ncol(miceData)]
    miceDataFlipped <- as.data.frame(t(miceData))
    return(miceDataFlipped)
}

doClustering <- function(dataTissue, tissue){
    # The following parameters are loaded in:
    #     dataTissue: a data frame containing the count data
    #     tissue: a String object
    # A file is loaded in, containing a list of strings which are
    # the matching genotypes.
    # Afterwards, each of the genotype gets a number before it
    # ( ie. "1 WT", "2 WT", "3 mdx") and is used as the rownames
    # of the dataTissue. Clustering is then performed, and the result
    # is saved in a PDF file.
    # The cluster object variable gets returned. 
    load(paste("~/Dropbox/Bafstu/Normalization/genotypes_",tissue,".RData",
               sep=""))
    rownames(dataTissue) <- paste(seq(1, length(genotypes)),  genotypes)
    samplesTreeMice <- hclust(dist(dataTissue), method = "average")
    pdf(paste("sampleClustering_",tissue,".pdf", sep=""), width = 12, height = 9)
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(samplesTreeMice, main = paste("Sample clustering to detect outliers in ",tissue, sep=""),
         sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    dev.off()
    return(samplesTreeMice)
}

changeClusters <- function(dataTissue, treeTissue, cutHeight){
    # The following parameters are loaded in:
    #     dataTissue: a data frame containing the count data
    #     treeTissue: a Hierarchical cluster object
    #     cutHeight: an integer
    # The original cluster object is cut depended on the given
    # cutHeight.
    # Afterwards, a specific cluster can be chosen with 'clust'. 
    # Only the remaining samples will be saved. 
    # The resulting dataset gets returned. 
    clust <- cutreeStatic(treeTissue, cutHeight = cutHeight, minSize = 10)
    table(clust)
    keepSamples <- (clust==0)
    datExprTissue <- dataTissue[keepSamples, ]
    return(datExprTissue)
}

createGenotypeColumn <- function(traitData){
    # The following parameters are loaded in:
    #     traitData: a data frame containing sample information
    # A column is created and added to traitData with the genotypes
    # in the format 'mdx' instead of 'mdx/utrn+/+'.
    # The edited traitData gets returned. 
    traitData[, 'mdx'] <- gsub('utrn', '', gsub('/', '', traitData[, 'strain1']))
    for (x in 1:nrow(traitData)){
      if (traitData[x, 'mdx'] == 'WT'){
        traitData[x, 'mdx'] <- '0'
      }
      if (traitData[x, 'mdx'] == 'mdx'){
        traitData[x, 'mdx'] <- '1'
      } else if (traitData[x, 'mdx'] == 'mdx++'){
        traitData[x, 'mdx'] <- '2'
      } else if (traitData[x, 'mdx'] == 'mdx+-'){
        traitData[x, 'mdx'] <- '3'
      }
      else {
        traitData[x, 'mdx'] <- '4'
      }
    }
    return(traitData)
}

getTraitData <- function(traitData, datExprTissue){
    # The following parameters are loaded in:
    #     traitData: a data frame containing sample information
    #     datExprTissue: a data frame containing the count data
    # Only the 'mdx' column is selected and saved in a new data frame.
    # Afterwards, only the samples are saved which are (aanwezig)
    # in datExprTissue.
    # The rownames of the resulting data frame gets replaced by
    # the rownames of datExprTissue.
    # The resulting data frame gets returned. 
    allTraits <- traitData[, c(11:12, 21)]
    datTraits <- allTraits[which(traitData$deCODEbarcode %in% rownames(datExprTissue)),]
    rownames(datTraits) <- rownames(datExprTissue)
    return(datTraits)
}

doClusteringAfterFiltering <- function(datExpr, tissue, datTraits, traitData){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the count data
    #     tissue: a String object
    #     datTraits: a data frame containing the trait data
    #     traitData: a data frame containing sample information
    # The genotypes from the samples which are in datTraits are retrieved. 
    # Afterwards, each of the genotype gets a number before it
    # ( ie. "1 WT", "2 WT", "3 mdx") and is used as the rownames
    # of the dataTissue. 
    # The dataTissue gets clustered and a dendogram is made.
    # The resulting figure is saved.
    collectGarbage()
    genotypes <- datTraits[,3]
    rownames(datExpr) <- paste(seq(1, length(genotypes)),  genotypes)
    sampleTree <- hclust(dist(datExpr), method = "average")
    traitColors <- numbers2colors(data.matrix(datTraits[1:(nrow(datTraits)),]), signed = T) 
    pdf(paste("sampleDendogram_", tissue, ".pdf", sep=""), width = 12, height = 9)
    plotDendroAndColors(sampleTree, traitColors,
                        groupLabels = names(datTraits),
                        main = paste("Sample dendrogram and trait heatmap of ", tissue))
    dev.off()
}

saveResults <- function(datExpr, datTraits, tissue){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the count data
    #     datTraits: a data frame containing the trait data
    #     tissue: a String object
    # The three parameters are saved.
    save(datExpr, datTraits, file=paste("Mice-01-dataInput-",tissue,".RData", sep=""))
}

main <- function(){
    # Two list are created containing two strings and integers and is
    # followed by a for-loop.
    # The part of a trait file is chosen which contains the data of
    # the tissue, dependent on the given tissue.
    # Afterwards, a specific cutHeight is chosen for each tissue.
    # Thereafter, a file containing the count data of the tissue
    # is read in and given with the tissue as string object to the
    # function 'preProcessingData'.
    # This results with a preprocessed data frame. 
    # The clustering is performed on the resulted data frame with
    # the function 'doClustering' with also the tissue as a String
    # object as a parameter. A hierarchical clustering object returns.
    # The function 'changeClusters' is called with the preprocessed
    # dataset, clustering object and cutHeight. 
    # The traitData is then given to the function
    # 'createGenotypeColumn'. The function 'getTraitData' is called
    # afterwards with traitData and datExprTissue. 
    # A data frame containing the trait data of the samples gets
    # returned and is given to the function 'doClusteringAfterFiltering'
    # with the String object of the tissue, datTraits and traitData.
    # After the clustering of the resulting data, the final dataset gets
    # saved with 'saveResults()' with datExprTissue, datTraits, and the
    # tissue as a String object.
    tissues <- c("Blood", "Muscle")
    cutHeight <- c(80000, 70000)
    for (x in 1:length(tissues)){
      if (tissues[x] == "Blood"){
        traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[1:96,]
      } else{
        traitData <- read.csv2("~/Dropbox/Bafstu/Spitali_Neuromics.csv")[97:147,]
      }
      miceData  <- read.csv2(paste("~/Dropbox/Bafstu/Normalization/Count_",tissues[x],"_Normalized.csv", sep=""))
      print("Preprocessing...")
      miceData <- preProcessingData(miceData, tissues[x])
      print("Clustering for removing outliers...")
      samplesTreeMice <- doClustering(miceData, tissues[x])
      print("Cutting clusters...")
      datExprTissue <- changeClusters(miceData, samplesTreeMice, cutHeight[x])
      traitData <- createGenotypeColumn(traitData)
      print("Retrieving corresponding trait data...")
      datTraits <- getTraitData(traitData, datExprTissue)
      print("Clustering resulting data...")
      doClusteringAfterFiltering(datExprTissue, tissues[x], datTraits, traitData)
      print("Saving results...")
      saveResults(datExprTissue, datTraits, tissues[x])
    }
}
main()
# R version 3.3.2 (2016-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 18
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=nl_NL.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=nl_NL.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=nl_NL.UTF-8 LC_IDENTIFICATION=C       
# attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
# other attached packages:
#  [1] zoo_1.7-13            matrixStats_0.51.0    WGCNA_1.51            RSQLite_1.0.0         DBI_0.5-1             fastcluster_1.1.21    dynamicTreeCut_1.63-1
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.7           RColorBrewer_1.1-2    plyr_1.8.4            iterators_1.0.8       tools_3.3.2           digest_0.6.10         rpart_4.1-10          preprocessCore_1.32.0 htmlTable_1.7        
#  [10] tibble_1.2            gtable_0.2.0          lattice_0.20-34       Matrix_1.2-7.1        foreach_1.4.3         parallel_3.3.2        gridExtra_2.2.1       stringr_1.1.0         knitr_1.15           
#  [19] cluster_2.0.5         S4Vectors_0.8.11      IRanges_2.4.8         stats4_3.3.2          grid_3.3.2            nnet_7.3-12           impute_1.44.0         data.table_1.9.6      Biobase_2.30.0       
#  [28] AnnotationDbi_1.32.3  survival_2.40-1       foreign_0.8-67        latticeExtra_0.6-28   Formula_1.2-1         magrittr_1.5          GO.db_3.2.2           ggplot2_2.2.0         htmltools_0.3.5      
#  [37] Hmisc_4.0-0           scales_0.4.1          codetools_0.2-15      splines_3.3.2         BiocGenerics_0.16.1   assertthat_0.1        colorspace_1.3-0      stringi_1.1.2         acepack_1.4.1        
#  [46] lazyeval_0.2.0        doParallel_1.0.10     munsell_0.4.3         chron_2.3-47         
