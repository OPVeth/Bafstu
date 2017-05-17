#!/usr/bin/Rscript
# This script has been created for performing the WGCNA on the MDX dataset
#    s1084813, Olga Veth of Hogeschool Leiden
#
# WGCNA - Retrieval Genes of Modules - Script 3 - 30-03-2017

library(WGCNA)
library(magrittr)
library(annotables)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
setwd("~/Dropbox/Bafstu/WGCNA")#Directory of interest

orderModuleEigengenes <- function(datExpr, moduleColors){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the count data
    #     moduleColors: a list of module colors as string
    # The moduleEigengenes are retrieved of datExpr and moduleColors
    # and get ordered. 
    # The ordered module eigengenes are returned.
    MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs <- orderMEs(MEs0)
    return(MEs)
}

createHeatmap <- function(moduleTraitCor, nSamples, tissue, MEs, datTraits){
    # The following parameters are loaded in:
    #     moduleTraitCor: correlation matrix of modules and traits
    #     nSamples: an integer
    #     tissue: a String object
    #     MEs: a list of eigengenes
    #     datTraits:  a data frame containing the trait data
    # A T-test is performed on moduleTraitCor with nSamples. 
    # Thereafter, the resulting values are changed to this format: 
    #     '<correlation value>'
    #     (<p-value>)
    # The resulting matrix gets the same dimension as moduleTraitCor
    # and afterwards a heatmap is created based on the matrix, datExpr,
    # moduleTraitCor and MEs. 
    # This figure is saved in a PDF file afterwards.
    moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples) #P-values
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n",
                        signif(moduleTraitPvalue, 1), "", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3)) 
    pdf(file = paste("Module_Trait_Relationship_",tissue,".pdf"), wi = 9, he = 6)
    labeledHeatmap(Matrix = moduleTraitPvalue,
                   xLabels = names(datTraits),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships of",tissue,"tissue"))
    dev.off()
}

createHeatmapGenotype <- function(datTraits, tissue, MEs, datExpr){
  # The following parameters are loaded in:
  #     tissue: a String object
  #     MEs: a list of eigengenes
  #     datTraits:  a data frame containing the trait data
  # The 'mdx' column is changed back into the corresponding genotype.
  # Afterwards, the 'mdx' column is changed into a factor.
  #   For every module eigengene, their significance gets 
  #   calculated as a p.value with the Kruskal-Wallis test.
  #   Afterwards, the resulting values are corrected with 
  #   fdr (False Discovery Rate). If these values are lower then 0.05
  #   they are saved into a variable. The variable gets leveled and
  #   combined.
  #     Afterwards a function is created which retrieves the module
  #     eigengenes, and checked if it is equal to the first and
  #     second position of combinations which is given as a parameter.
  #     The Wilcoxon test is performed on the eigen genes of the
  #     combinations and the p-value is eventually returned. 
  #   The resulting p-values are corrected with FDR.
  # The matrix of p-values is checked if there are p-values below 0.05
  # and thus sigificant. Thereafter a matrix containing one row with all the adjusted p-values
  # gets created with the column names of textMatrix.
  # Afterwards the part with the significant p-values gets binded as a matrix.
  # A heatmap is finally created of the matrix.
  for (x in 1:length(datTraits$mdx)){
    if (datTraits$mdx[x] == 0){
      datTraits$mdx[x] <- "WT"
    } 
    else if (datTraits$mdx[x] == 1){
      datTraits$mdx[x] <- "mdx"
    } else if (datTraits$mdx[x] == 2){
      datTraits$mdx[x] <- "mdx++"
    } else if (datTraits$mdx[x] == 3){
      datTraits$mdx[x] <- "mdx+-"
    } else {
      datTraits$mdx[x] <- ""
    }
  }
  group <- as.factor(datTraits$mdx)
  p.values <- sapply(colnames(MEs), function(col) kruskal.test(MEs[, col], group)$p.value)
  p.values.adjusted <- p.adjust(p.values, 'fdr')
  
  combinations <- combn(levels(group), 2)
  textMatrix <- sapply(colnames(MEs), function(col) {
    eigengene <- MEs[, col]
    p.values <- apply(combinations, 2, function(combination) {
      group1_eigengene <- eigengene[which(group == combination[1])]
      group2_eigengene <- eigengene[which(group == combination[2])]
      return(wilcox.test(group1_eigengene, group2_eigengene)$p.value)
    })
    return(p.adjust(p.values, 'fdr'))
  })
  mat <- (textMatrix < 0.05) * 1
  sigMat <- matrix(data=p.values.adjusted, nrow = 1)
  colnames(sigMat) <- colnames(textMatrix)
  textMatrix <- rbind(sigMat, textMatrix)
  
  mat <- rbind((sigMat < 0.05) * 1, mat)
  cols <- c('Significance', apply(combinations, 2, function(c) paste(c[1], 'vs', c[2])))
  par(mar = c(6, 8.5, 3, 3))
  pdf(file = paste("Module_Trait_Relationship_Genotype_",tissue,".pdf", sep=""), wi = 9, he = 6)
  labeledHeatmap(Matrix = t(mat),
                 xLabels = cols,
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = t(signif(textMatrix, 2)),
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 plotLegend = FALSE,
                 main = paste("Module-genotype relationships"))
  dev.off()

}

createHeatmapMEs <- function(datTraits, tissue, MEs, datExpr){
  # The following parameters are loaded in:
  #     datTraits:  a data frame containing the trait data
  #     tissue: a String object
  #     datExpr:  a data frame containing the count data
  # A CSV file with the ME's of the muscle is read in, and the first column is skipped.
  # The datTraits and datExpr of the blood is saved in another variable.
  # Thereafter the traitData is read in and a RData file of muscle is loaded in.
  # Only the rows of blood and muscle are retrieved and then their mouseID in traitData.
  # Only the matching mouseID's rows are thereafter retrieved.
  # After calculating the correlation between the ME's of both tissue, the p-value
  # gets also calculated. 
  # The moduleTraitPvalue gets corrected for the amount of numbers and is saved in a variable.
  # A heatmap is created between the ME's of the muscle and modules of blood and afterwards
  # saved in a PDF file.
  MEs_Muscle <- data.frame(read.csv2("MEs_Muscle_Week30.csv")[, c(-1)])
  datTraits_Blood <- datTraits
  datExpr_Blood <- datExpr
  load(file = paste("Mice-01-dataInput-Muscle.RData", sep=""))
  traitData <- read.csv2("Spitali_Neuromics.csv")
  traitBlood <- traitData[which(traitData$deCODEbarcode %in% rownames(datTraits_Blood)),]$mouseID
  traitMuscle <- traitData[which(traitData$deCODEbarcode %in% rownames(datTraits)),]$mouseID
  MEs_Muscle <- MEs_Muscle[which(traitMuscle %in% traitBlood),]
  moduleTraitCor <- cor(MEs, MEs_Muscle, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_Blood))
  textMatrix <- paste(signif(moduleTraitPvalue, 1), sep = "")
  par(mar = c(6, 8.5, 3, 3)) 
  pdf(file = "Module_Trait_Relationship_MEs.pdf", wi = 9, he = 6)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(MEs_Muscle),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(60),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships of blood tissue"))
  dev.off()
}

getCorTrait <- function(datTraits, MEs, datExpr, nSamples, weight){
    # The following parameters are loaded in:
    #     datTraits: a data frame containing the count data
    #     MEs: the modules
    #     datExpr: a data frame containing the trait data
    #     nSamples: an integer
    #     weight: a data frame containing the trait data
    # The correlation between the expression data and the trait gets
    # calculated. 
    # The names of the resulting matrix consists of a 'GS. ', followed
    # by the weight names. 
    # The T-test is performed on the resulting correlation matrix. 
    # This results in another data frame which names are changed into
    # 'p.GS' with the names of the weight afterwards. 
    # The correlation matrix gets returned.
    geneTraitSignificance <- as.data.frame(cor(datExpr, weight, use = "p"))
    names(geneTraitSignificance) <- paste("GS.", names(weight), sep="")
    GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) 
    names(GSPvalue) <- paste("p.GS.", names(weight), sep="")
    return(geneTraitSignificance)
}

getCorME <- function(datExpr, MEs, nSamples, modNames){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the trait data
    #     MEs: the modules
    #     nSamples: an integer
    #     modNames:  a String object
    # The correlation between the expression data and the module eigengenes
    # gets calculated and the t-test is performed on the resulting data.
    # The names of the resulting matrix are changed into 'p.MM' combined 
    # with the module names.
    # Afterwards, the matrix is returned. 
    geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    names(MMPvalue) <- paste("p.MM", modNames, sep="")
    return(MMPvalue)
}

makeScatterPlot <- function(tissue, geneModuleMembership, moduleGenes, column, geneTraitsSignificance, module){
    # The following parameters are loaded in:
    #     tissue: a String object
    #     geneModuleMembership: correlation matrix between expression data and
    #                           module eigengenes
    #     moduleGenes: a module
    #     column: a String object
    #     geneTraitsSignificance: correlation matrix between trait data and
    #                           module eigengenes
    #     module: a String object
    # A scatterplot is created based on the given module and saved as a PDF file.
    par(mfrow = c(1,1))
    print(dim(geneModuleMembership))
    pdf(file = paste("Module_Membership_Gene_Significance_",tissue,".pdf"), wi = 9, he = 6)
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitsSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for body weight", 
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,
                       abline=TRUE, abline.color = "red", abline.lty = 1)
    dev.off()
}

calculateMean <- function(moduleColors, datExpr){
    # The following parameters are loaded in:
    #     moduleColors: a list of all colors of the modules
    #     datExpr: a data frame containing the count data
    # A variable, 'calMean' is set on NULL and with a for-loop, 
    # the amount of all the genes of a specific module are calculated
    # and saved in a vector.
    # Thereafter, the mean is calculated of the vector and shown.
    calMean <- NULL
    for(i in 1:length(moduleColors)){
      calMean$genes[i] <- length(names(datExpr)[moduleColors==moduleColors[i]])
    }
    print(unique(moduleColors))
    print(unique(calMean$genes))
    print(mean(calMean$genes))
}

createMap <- function(){
  # A data frame is created wth two columns containing all gene Symbols
  # and their ENTREZ ID's. 
  # The column names are changed into 'ENTREZ' and'SYMBOL' and the data
  # frame is written to a CSV file 
  mapData <- data.frame(grcm38$entrez, grcm38$symbol) ##Mouse
  colnames(mapData) <- c("ENTREZ", "SYMBOL")
  write.table(mapData, file="Entrez_Map.csv", sep=";")
}

getSymbols <- function(dataTissue, mapData){
  # The following parameters are loaded in:
  #     dataTissue: a data frame containing the count data
  #     mapData: a data frame
  # Only the columns of dataTissue which match with the symbols within mapData
  # are selected. 
  # Within a for-loop it is checked which symbol matches with the column name.
  # If there's more than one matche, then the first match will be used. 
  # Otherwise, the matching one will be used.
  # After the for-loop, dataTissue gets returned. 
  for (x in 1:length(colnames(dataTissue))){
    location <- which(mapData$SYMBOL %in% colnames(dataTissue)[x])
    if (length(location) > 1){
      colnames(dataTissue)[x] <- as.character(mapData$ENSEMBLE[location[1]])
    } else if (length(location) == 1) {
      colnames(dataTissue)[x] <- as.character(mapData$ENSEMBLE[location]) 
    }
  }
  return(dataTissue)
}

saveGeneList <- function(datExpr, moduleColors, module, tissue, genes){
    # The following parameters are loaded in:
    #     datExpr: a data frame containing the trait data
    #     moduleColors: a list of all colors of the modules
    #     module: a String object
    #     tissue: a String object
    #     genes: a String object, can be empty depending if hubgenes
    #            are given or not
    # The function 'createMap()' is called and afterwards, a mapping file is read in.
    # A string is created and thereafter, the 'genes' is checked if it contains a String.
    #   If that's not the case, all genes of the given module are retrieved and the created string
    #   is equalized to an empty string. 
    #   Afterwards all the genes are saved into a TXT file
    #   after their gene symbols have been converted to an Entrez ID.
    #createMap()
    #mapData <- read.csv2("Entrez_Map.csv", sep=";")
    mapData <- read.csv2("Ensemble_Map.csv", sep=";")
    mapData2 <- read.csv2("Ensemble_Human_Map.csv", sep=";")
    mapData3 <- read.csv2("Entrez_Map.csv", sep=";")
    mapData4 <- read.csv2("Entrez_Human_Map.csv", sep=";")
    hub <- "Hub"
    modules <- c("pink", "blue", "brown", "green", "magenta", "tan", "grey60", "lightyellow", "darkgreen", "yellow", "greenyellow", "red", "turquoise", "royalblue", "cyan", "midnightblue", "lightgreen", "black", "salmon", "darkred", "purple")  
    for (x in 1:length(modules)){
      module <- modules[x]
      #if (genes == "x"){
        genes <- names(datExpr)[moduleColors==module]
        genes2 <- names(datExpr)[moduleColors==module]
        genes3 <- names(datExpr)[moduleColors==module]
        hub <- ""
      #}
      for (x in 1:length(genes)){
        location <- which(mapData$SYMBOL %in% genes[x])
        if (length(location) > 1){
          genes[x] <- as.character(mapData$ENSEMBLE[location[1]])
        } else if (length(location) == 1) {
          genes[x] <- as.character(mapData$ENSEMBLE[location]) 
        }
      }
      genes <- toupper(genes)
      for (x in 1:length(genes)){
        location <- which(mapData2$SYMBOL %in% genes[x])
        if (length(location) > 1){
          genes[x] <- as.character(mapData2$ENSEMBLE[location[1]])
        } else if (length(location) == 1) {
          genes[x] <- as.character(mapData2$ENSEMBLE[location]) 
        }
        else if (!grepl("EN", genes[x])) {
          genes[x] <- "NA"
        }
      }
      for (x in 1:length(genes3)){
        location <- which(mapData3$SYMBOL %in% genes3[x])
        if (length(location) > 1){
          genes3[x] <- as.character(mapData3$ENTREZ[location[1]])
        } else if (length(location) == 1) {
          genes3[x] <- as.character(mapData3$ENTREZ[location]) 
        }
      }
      for (x in 1:length(genes3)){
        location <- which(mapData4$SYMBOL %in% genes3[x])
        if (length(location) > 1){
          genes3[x] <- as.character(mapData4$ENTREZ[location[1]])
        } else if (length(location) == 1) {
          genes3[x] <- as.character(mapData4$ENTREZ[location]) 
        }
        else if (grepl('^[A-Za-z]+$', genes3[x])) {
          genes3[x] <- "NA"
        }
      }
      genes2 <- genes2[!is.na(genes2)] # Symbol
      moduleNames <- rep(module, length(genes2))
      genesData <- t(rbind(genes2, genes, genes3, moduleNames)) #Symbol; Ensembl; Entrez
      write.csv(genesData, file=paste("Dystrophic_Genotypes/Blood/ID/Module_", module, "_", tissue, ".csv"))
      }
    #write(genes, file=paste("Blood/Ensembl/Module_",module,"_",tissue,"_", hub, "Genelist.txt", sep=""))
}

writeCountData <- function(datExpr, genes, module, tissue, datTraits){
  # The following parameters are loaded in:
  #     datExpr:  a data frame containing the trait data
  #     genes: a vector of String objects
  #     module: a String object
  #     datTraits: a data frame containing the count data
  for (x in 1:length(datTraits$mdx)){
    #if (datTraits$mdx[x] == 0){
    #  datTraits$mdx[x] <- "WT"
    #}
    if (datTraits$mdx[x] == 1){
      datTraits$mdx[x] <- "mdx"
    } else if (datTraits$mdx[x] == 2){
      datTraits$mdx[x] <- "mdx++"
    } else if (datTraits$mdx[x] == 3){
      datTraits$mdx[x] <- "mdx+-"
    } else {
      datTraits$mdx[x] <- ""
    }
  }
  countData <- datExpr[, which(colnames(datExpr) %in% genes)]
  countData[,6] <- datTraits$mdx
  write.csv2(countData, file = paste("Normalized_CountData_HubGenes_", module, "_", tissue, ".csv", sep=""))
}

main <- function(){
    tissues <- c("Muscle")
    modules <- c("salmon") #Color of module
    for(x in 1:length(tissues)){
      load(file = paste("Mice-01-dataInput-",tissues[x],".RData", sep=""))
      load(file = paste("Mice-02-networkConstruction_",tissues[x],".RData", sep=""))
      nSamples <- nrow(datExpr)
      print("Ordering module eigengenes...")
      print(moduleColors)
      MEs <- orderModuleEigengenes(datExpr, moduleColors)
      if (tissues[x] == "Muscle"){
        write.csv2(MEs, "MEs_Muscle_Week30.csv")
      }
      print("Creating heatmap...")
      weight <- as.data.frame(datTraits$mdx) 
      names(weight) = "weight"
      modNames <- substring(names(MEs), 3)
      print("Creating heatmap...")
      moduleTraitCor <- cor(MEs, datTraits, use = "p")
      print("Calculating correlation between expression data and trait...")
      saveGeneList(datExpr, moduleColors, modules[x], tissues[x], "x")
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
