#!/usr/bin/Rscript
# This script has been created for performing the WGCNA on the MDX dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# WGCNA - Modules and Eigengenes - Script 2 - 31-01-2017

library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
setwd("")#Directory of interest

getSoftThreshold <- function(datExpr){
    #  The following parameter is loaded in:
    #     datExpr: a data frame containing the count data
    # A vector of integers ranging from 1-10 and 12-20 with steps of two is created.
    # Afterwards, the soft thresholding is performed based on the given powers. The
    # result is returned. 
    memory.limit(size = 60000)
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, 
                             networkType = 'signed') 
    return(sft)
}

createTopologyThreshold <- function(sft, tissue){
    # The following parameters are loaded in:
    #     sft: a Soft Threshold object
    #     tissue: a String object
    # A vector of integers ranging from 1-10 and 12-20 with steps of two are created.
    # Thereafter, two plots are created which represent the topology connectivity and the 
    # mean connectivity of each power. The plots are saved as a PDF file. 
    powers <- c(c(1:10), seq(from = 12, to=20, by=2))
    sizeGrWindow(9, 5)
    par(mfrow = c(1,2))
    cex1 = 0.9
    pdf(paste("Topology_Connectivity_",tissue,".pdf", sep=""))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
         main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
    dev.off()
}

doTOM <- function(sft, datExpr, power){
    # The following parameters are loaded in:
    #     sft: a Soft Threshold object
    #     datExpr: a data frame containing the count data
    #     power: an integer
    # First the adjacency is calculated of datExpr with the assigned power and afterwards, 
    # the TOM is calculated. One is distracted of the TOM, resulting in the dissimilarity
    # matrix of the TOM. Because of the big sizes of the TOM and adjacency, they get removed.
    adjacency <- adjacency(datExpr, power = power, type='signed')
    TOM <- TOMsimilarity(adjacency, TOMType = 'signed')
    dissTOM <- 1-TOM
    rm(adjacency)
    rm(TOM)
    return(dissTOM)
}

doClustering <- function(dissTOM, tissue){
    # The following parameters are loaded in:
    #     dissTOM: matrix
    #     tissue: a String object
    # The dissTOM gets clustered and its result is placed in a plot. 
    # The dissimilarity plot is afterwards saved in a PDF file. 
    # The hierarchical cluster object 'geneTree' is then returned. 
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    sizeGrWindow(12,9)
    pdf(paste("Dissimilarity_Matrix_",tissue,".pdf", sep=""))
    plot(geneTree, xlab="", sub="", main = 
           paste("Gene clustering on TOM-based dissimilarity on ",tissue),
         labels = FALSE, hang = 0.04)
    dev.off()
    return(geneTree)
}

createDendogram <- function(dissTOM, geneTree, deepSplit, cutHeight, minModuleSize, tissue){
    # The following parameters are loaded in:
    #     dissTOM: matrix
    #     geneTree: a hierarchical clustering object
    #     deepSplit: an integer
    #     cutHeight: an integer
    #     minModuleSize: an integer
    #     tissue: a String object
    # First, the hierarchical cluster object gets cut with
    # the given cutHeight, deepSplit, minModuleSize and dissTOM.
    # The resulting modules get their corresponding colors with
    # the function 'labels2colors'. 
    # Afterwards,  the modules and hierarchical cluster object are
    # plotted in a dendogram. The figure gets saved and the list of
    #the color labels of the modules gets returned.
    dynamicMods <- cutreeDynamic(dendro = geneTree, cutHeight = cutHeight, distM = dissTOM,
                                       deepSplit = deepSplit, pamRespectsDendro = FALSE,
                                       minClusterSize = minModuleSize, method="hybrid")
    dynamicColors <- labels2colors(dynamicMods)
    sizeGrWindow(8,6)
    pdf(paste("Gene_Dendogram_Module_Colors_",tissue,"_", minModuleSize,".pdf", sep=""))
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = paste("Gene dendrogram and module colors in",tissue))
    dev.off()
    return(dynamicColors)
}

getModuleEigengenes <- function(datExpr, dynamicColors){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the count data
    #     dynamicColors: a list of the module colors as strings
    # ModuleEigengenes are calculated of the datExpr with dynamicColors
    # as the corresponding colors.
    # Thereafter, the eigengenes are retrieved and returned. 
    MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
    MEs <- MEList$eigengenes
    return(MEs)
}

doClusteringMEs <- function(MEs, minModuleSize, tissue){
    # The following parameters are loaded in:
    #     MEs: a list of eigengenes
    #     minModuleSize: an integer
    #     tissue: a String object
    # Correlation of the MEs is calculated and with one destracting of this,
    # the dissimilarity value results.
    # This value is used for custering and gets clustered. 
    # The figure is saved in a PDF file. 
    MEDiss <- 1-cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    sizeGrWindow(7, 6)
    pdf(paste("Clustering_Module_Eigengenes_",tissue,"_",minModuleSize,".pdf", sep=""))
    plot(METree, main = "Clustering of module eigengenes",
         xlab = "", sub = "")
    dev.off()
}

mergeModules <- function(datExpr, dynamicColors, MEDissThreshold, minModuleSize, tissue, geneTree){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the count data
    #     dynamicColors: a list of the module colors as strings 
    #     MEDissThreshold: an integer
    #     minModuleSize: an integer
    #     tissue: a String object
    #     geneTree: a hierarchical cluster object
    # First, the modules of datExpr gets clustered with the colors within
    # dynamicColors and with the merging threshold of MEDissThreshold. 
    # The colors of the merged modules are retrieved afterwards, and get returned. 
    merge <- mergeCloseModules(data.matrix(datExpr), dynamicColors, 
                               cutHeight = MEDissThreshold, verbose = 3)
    mergedColors <- merge$colors
    return(mergedColors)
}

createDendogramMerged <- function(geneTree, tissue, minModuleSize, MEDissThreshold, dynamicColors, mergedColors){
    # The following parameters are loaded in:
    #     geneTree: a hierarchical cluster object
    #     tissue: a String object
    #     minModuleSize: an integer
    #     MEDissThreshold: an integer
    #     dynamicColors: a list of the module colors as strings 
    #     mergedColors: a list of the module colors as strings after merging
    # A dendogram is produced of the geneTree with the non-merged and merged modules. 
    # The figure is saved as a PDF file afterwards.
    sizeGrWindow(12, 9)
    pdf(file = paste("Dynamic_Cut_Merged_Dendogram_",tissue,"_",minModuleSize,"_",
                     MEDissThreshold,".pdf", sep=""), wi = 9, he = 6)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    dev.off()
}

getModuleLabels <- function(mergedColors){
    # The following parameter is loaded in:
    #     mergedColors: a list of module colors as strings after merging
    # A vector containing different colors, with 'grey' as the first index
    # is created and afterwards matched with the mergedColors. 
    # Each matching position is reduced by one and the resulting list matching
    # colors get returned. 
    colorOrder <- c("grey", standardColors(50))
    moduleLabels <- match(mergedColors, colorOrder)-1
    return(moduleLabels)
}

saveResults <- function(MEs, moduleLabels, moduleColors, geneTree, tissue){
    # The following parameters are loaded in:
    #     MEs: a list of eigengenes
    #     moduleLabels: a list of module colors as strings after matching
    #     moduleColors: a list of module colors as string
    #     geneTree: a hierarchical cluster object
    #     tissue: a String object
    # The parameters are saved in a RData file.
    save(MEs, moduleLabels, moduleColors, geneTree, 
         file = paste("Mice-02-networkConstruction_",tissue,".RData", sep=""))
}

main <- function(){
    # Multiple vectors are created, containing integers which are tissue specific. 
    # These variables play a role in the functions afterwards and thus, the creation
    # of the modules and their eigengenes. 
    #   For each tissue, the corresponding RData file from the first script is loaded
    #   in and afterwards getSoftThresholding gets called with datExpr which contains
    #   the expression data. 
    #   The resulting softThresholding object is given with the tissueas a String object
    #   to 'createoftThreshold'. 
    #   The function 'doTOM()' is called with the softThresholding object, datExpr and 
    #   the corresponding power integer. 
    #   The resulting dissimilarity matrix is given to the function 'doClustering'. 
    #   A hierarchical cluster object is thus created. 
    #   Afterwards, createDendogram() gets called with geneTree, the corresponding deepSplit,
    #   cutHeight, minModuleSize and tissue. 
    #   A list of colors as strings is retrieved from the function and given to the function
    #   'getModuleEigengenes' with datExpr.
    #   A clustering of the moduleEigengenes is afterwards performed with the parameters MEs,
    #   minModulSie and the tissue as a String object.
    #   The function 'mergeModules is afterwards called with datExpr, dynamicColors, 
    #   MEDissThreshold, minModuleSize, tissues[x] and geneTree.
    #   This resulted in a list of the colors of the merged modules and is given to the function
    #   'createDendogramMerged', together with the tissue as a String object, the corresponding
    #   minModuleSize, MEDissThreshold, dynamicColors and mergedColors.
    #   The new module eigengenes are retrieved from the resulting list of merged modules and
    #   given to the function 'getModuleLabels'.
    #   A list of the colors of the modules is produced and the MEs, moduleLabels, moduleColors
    #   and geneTree get saved after calling the fuction saveResults(). 
    tissues <- c("Blood", "Muscle")
    powers <- c(20, 16)
    minModuleSize <- c(30, 30)
    deepSplit <- c(0, 2)
    MEDissThreshold <- c(0.15, 0.0001)
    cutHeight <- c(0.995, 0.999)
    for (x in 1:length(tissues)){
      load(file = paste("Mice-01-dataInput-",tissues[x],".RData", sep=""))
      print("Softhresholding...")
      sft <- getSoftThreshold(datExpr)
      createTopologyThreshold(sft, tissues[x])
      print("Performing TOM...")
      dissTOM <- doTOM(sft, datExpr, powers[x])
      print("Clustering...")
      geneTree <- doClustering(dissTOM, tissues[x])
      print("Creating dendogram...")
      dynamicColors <- createDendogram(dissTOM, geneTree, deepSplit[x], cutHeight[x], 
                                      minModuleSize, tissues[x])
      print("Get Module Eigengenes...")
      MEs <- getModuleEigengenes(datExpr, dynamicColors)
      print("Clustering Module Eigengenes...")
      doClusteringMEs(MEs, minModuleSize, tissues[x])
      print("Merging modules...")
      merge <- mergeCloseModules(data.matrix(datExpr), dynamicColors, 
                                 cutHeight = MEDissThreshold, verbose = 3)
      mergedColors <- merge$colors
      print("Creation merged dendogram...")
      createDendogramMerged(geneTree, tissues[x], minModuleSize, MEDissThreshold, 
                           dynamicColors, mergedColors)
      mergedMEs <- merge$newMEs
      moduleLabels <- getModuleLabels(mergedColors)
      print("Saving results...")
      saveResults(MEs, moduleLabels, mergedColors, geneTree, tissues[x])
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
#  [1] Rcpp_0.12.7           compiler_3.3.2        RColorBrewer_1.1-2    plyr_1.8.4            iterators_1.0.8       tools_3.3.2           digest_0.6.10         rpart_4.1-10          preprocessCore_1.32.0
#  [10] htmlTable_1.7         tibble_1.2            gtable_0.2.0          lattice_0.20-34       Matrix_1.2-7.1        foreach_1.4.3         parallel_3.3.2        gridExtra_2.2.1       stringr_1.1.0        
#  [19] knitr_1.15            cluster_2.0.5         S4Vectors_0.8.11      IRanges_2.4.8         stats4_3.3.2          grid_3.3.2            nnet_7.3-12           impute_1.44.0         data.table_1.9.6     
#  [28] Biobase_2.30.0        AnnotationDbi_1.32.3  survival_2.40-1       foreign_0.8-67        latticeExtra_0.6-28   Formula_1.2-1         magrittr_1.5          GO.db_3.2.2           ggplot2_2.2.0        
#  [37] htmltools_0.3.5       Hmisc_4.0-0           scales_0.4.1          codetools_0.2-15      splines_3.3.2         BiocGenerics_0.16.1   assertthat_0.1        colorspace_1.3-0      stringi_1.1.2        
#  [46] acepack_1.4.1         lazyeval_0.2.0        doParallel_1.0.10     munsell_0.4.3         chron_2.3-47
