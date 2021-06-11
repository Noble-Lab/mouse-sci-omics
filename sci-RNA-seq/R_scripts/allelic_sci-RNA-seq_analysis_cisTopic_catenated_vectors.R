##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200508 allelic sci-RNA-seq analysis using cisTopic
#          with concatened chromosome-specific vectors
#          for F121 cells only
#
##################################################
##################################################
##################################################



##################################################
# Set up
#
# R/3.6.1
#
# https://github.com/aertslab/cisTopic
# devtools::install_github("aertslab/cisTopic")
# install.packages("umap")
#
##################################################



rm(list=ls())
gc()
sessionInfo()
    









##################################################
##################################################
##################################################
# Packages and functions

# require(monocle)
require(data.table)
require(plotly)
require(scatterplot3d)
require(ggplot2)



##################################################
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk
t_col <- function(color, percent = 50, name = NULL) {
    #	  color = color name
    #	percent = % transparency
    #	   name = an optional name for the color
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    ## Save the color
    invisible(t.col)
}
    









##################################################
##################################################
##################################################
# Input

# 20200720
# min_counts_per_chromosome=0
min_counts_per_chromosome=10

projectName = "allelic_sci-RNA-seq"

inputDir = paste("../data/allelic_sci-RNA-seq_catenated_count_matrices/")
    
mainOutputDir = paste("../results/allelic_sci-RNA-seq_cisTopic_catenated_allelic_vectors/")
if (min_counts_per_chromosome>0) mainOutputDir = paste0( mainOutputDir, ".min_counts_per_chromosome", min_counts_per_chromosome )
if ( ! file.exists( mainOutputDir ) ) {
  if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")

tpCols2=c("red", "blue", "orange", "green", "black")
names(tpCols2) = c("d0", "d11", "d3", "d7", "NPC")

pchVal=20
pchcexVal=1.25

binCountThresh=1
minCells=1
minPeaks=1

whichChrom="genome"
# whichChrom="autosomes"
# whichChrom="chr1_and_chrX"
# whichChrom="chrX"
# whichChrom="chr1"










##################################################
##################################################
##################################################
# Load counts data

if (whichChrom == "genome") {
    
    count.table.OG = c()
    for (theChrom in c( paste0( "chr", 1:19 ), "chrX" ) ) { 
        inputFileName = paste0("allelic_sci-RNA-seq.UMImatrix.filtered.minExpr0.1.minCells5.F121.",
                               theChrom, 
                               ".stringentlyFiltered",
                               ifelse(min_counts_per_chromosome>0, paste0(".min_counts_per_chromosome", min_counts_per_chromosome), ""), 
                               ".tsv.gz" )
        countsFile = file.path(inputDir, inputFileName)
        count.table.OG = rbind( count.table.OG,
                                as.data.frame(fread(countsFile, header=TRUE, sep="\t") ) )
    }
    
} else if (whichChrom == "autosomes") {
    
    count.table.OG = c()
    for (theChrom in c( paste0( "chr", 1:19 ) ) ) { 
        inputFileName = paste0("allelic_sci-RNA-seq.UMImatrix.filtered.minExpr0.1.minCells5.F121.",
                               theChrom, 
                               ".stringentlyFiltered",
                               ifelse(min_counts_per_chromosome>0, paste0(".min_counts_per_chromosome", min_counts_per_chromosome), ""), 
                               ".tsv.gz" )
        countsFile = file.path(inputDir, inputFileName)
        count.table.OG = rbind( count.table.OG,
                                as.data.frame(fread(countsFile, header=TRUE, sep="\t") ) )
    }
    
} else if (whichChrom == "chr1_and_chrX") {
    
    count.table.OG = c()
    for (theChrom in c( "chr1", "chrX" ) ) { 
        inputFileName = paste0("allelic_sci-RNA-seq.UMImatrix.filtered.minExpr0.1.minCells5.F121.",
                               theChrom, 
                               ".stringentlyFiltered",
                               ifelse(min_counts_per_chromosome>0, paste0(".min_counts_per_chromosome", min_counts_per_chromosome), ""), 
                               ".tsv.gz" )
        countsFile = file.path(inputDir, inputFileName)
        count.table.OG = rbind( count.table.OG,
                                as.data.frame(fread(countsFile, header=TRUE, sep="\t") ) )
    }
    
} 
else { 
    inputFileName = paste0("allelic_sci-ATAC-seq.UMImatrix.filtered.minExpr0.1.minCells5.F121.",
                           whichChrom, ".tsv.gz" )
    countsFile = file.path(inputDir, inputFileName)
    
    count.table.OG = as.data.frame(fread(countsFile, header=TRUE, sep="\t"))
}

rownames(count.table.OG) = count.table.OG[,1]
count.table.OG = count.table.OG[,-1]

dim(count.table.OG)
# 20966  1048 = genome
# 20308  1048 = autosomes
# 2126 1048 = chr1_and_chrX










##################################################
##################################################
##################################################
# Run cisTopic

library(cisTopic)
library(Rtsne)
library(fastcluster)
library(NMF)
library(data.table)
library(ComplexHeatmap)

# generateModels = FALSE
generateModels = TRUE

numModels = NULL
# numModels = 15

analyzeTopicsOnly = FALSE
# analyzeTopicsOnly = TRUE

if (analyzeTopicsOnly) {
    generateModels = FALSE 
}

# Filter out top/bottom 5%? 
# No, already pre-filtered.
# filter90 = TRUE
filter90 = FALSE

# Filter lowly covered and biased peaks?
# No, already pre-filtered.
# filterBiasedPeaks = TRUE
filterBiasedPeaks = FALSE










##################################################
##################################################
# Filter data

dim(count.table.OG)

# New loci
locIDs.OG <- rownames(count.table.OG) 

count.matrix = count.table.OG
# rownames(count.matrix)

if (filterBiasedPeaks) {
    count.matrix = count.matrix[ !(refBiasedPeakBoolIdx | altBiasedPeakBoolIdx), ] 
}
dim(count.matrix)
locIDs <- rownames(count.matrix) 

# only consider middle 90% of cells
if (filter90) {
    acct <- colSums(count.matrix != 0)
    quant = 0.05
    low_lim <- quantile(acct, quant)
    upp_lim <- quantile(acct, 1-quant)
    sum_f = acct[(acct > low_lim)]
    sum_f = as.matrix(sum_f[(sum_f < upp_lim)])
    length(rownames(sum_f))
    count.matrix <- count.matrix[,rownames(sum_f)]
    cell_list <- as.matrix(count.matrix[rownames(sum_f)])
    # rownames(count.matrix) = locIDs.OG
    dim(count.matrix)
}










##################################################
##################################################
# Ouput folder and cisTopic object

# Results folder
nameText = paste( projectName, 
                  paste(".", whichChrom, sep=""),  
                  ifelse(filterBiasedPeaks, ".filterBiasedPeaks", ""),
                  ifelse(filter90, ".filter90", ""),
                  ifelse(binCountThresh>1, paste("_binCountThresh", binCountThresh, sep=""), ""),
                  ifelse(minCells>1, paste("_minCells", minCells, sep=""), ""),
                  ifelse(minPeaks>1, paste("_minPeaks", minPeaks, sep=""), ""),
                  sep="" ) 
currOutputDir = file.path( mainOutputDir, 
                           paste("cisTopicResults.", 
                                 nameText, 
                                 sep="" ))
if ( ! file.exists( currOutputDir ) ) {
    if ( ! dir.create( currOutputDir, showWarnings=FALSE, recursive=TRUE )  ) currOutputDir = mainOutputDir
}

if (!generateModels) {
    cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
} else {
    cisTopicObject <- createcisTopicObject(count.matrix,
                                           project.name=paste(projectName, sep=""), 
                                           keepCountsMatrix = TRUE,
                                           is.acc=binCountThresh,
                                           min.regions=minPeaks,
                                           min.cells=minCells)
    
    
    
    ##################################################
    ##################################################
    # runModel
    
    cisTopicObject <- runCGSModels(cisTopicObject,  topic=c(seq(5,40,5)), seed=999, nCores=1, burnin = 250, iterations = 500)
    # cisTopicObject <- runCGSModels(cisTopicObject,  topic=c(seq(45,70,5)), seed=999, nCores=1, burnin = 250, iterations = 500)
    
    
    
    ##################################################
    ##################################################
    # modelSelection and results
    
    PDFfilename = file.path(currOutputDir, paste("cisTopic.modelSelection.", nameText, 
                                                 ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                 ".pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    logLikelihoodByIter(cisTopicObject)
    cisTopicObject <- selectModel(cisTopicObject, select = NULL, type = "maximum", 
                                  keepBinaryMatrix = TRUE, keepModels = TRUE)
    # select
    # Number of topics of the selected model. 
    # If NULL, the model with the best log likelihood is picked.
    # 
    # type 	
    # Method for automatic selection of the best number of topics. 
    # By default, we use 'derivative' which calculates the second derivative in each point 
    # and takes the maximum (recommended with runWarpLDAModels, where curves are less noisy). 
    # Alternatively, if set to 'maximum', the model with the maximum log-likelihood is selected 
    # (recommended with runCGSModels). <- *************************************************
    # For WarpLDA models, the minimum perplexity can be also used to select the best model.
    
    dev.off()
    
    
    
    ##################################################
    ##################################################
    # Embeddings
    
    cisTopicObject <- runPCA(cisTopicObject, target='cell')  # specify target in v0.2.0
    # cisTopicObject <- runtSNE(cisTopicObject, seed=999, target='cell')
    cisTopicObject <- runUmap(cisTopicObject, seed=999, target='cell')  # Available from 0.2.0
    
    
    
    ##################################################
    ##################################################
    # save RDS        
    saveRDS(cisTopicObject, file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")))
    
}



if (!analyzeTopicsOnly) {
    
    ##################################################
    ##################################################
    # Setup metadata
    
    ##########        
    # time-points
    timepoints = as.matrix(sub("F121_", "", sub("\\..*$", "", colnames(count.matrix)))) 
    rownames(timepoints) <- colnames(count.matrix)
    colnames(timepoints) <- "timepoints"
    
    table(timepoints)
    # timepoints
    # d0 d11  d3  d7 NPC 
    # 489 145 127  78 209 

    ##########        
    # Colors
    colsToUse = list()
    colsToUse[['timepoints']] = tpCols2
    
    # Transparent colors
    tcolsToUse = list()
    tcolsToUse[['timepoints']] = apply(as.matrix(tpCols2), 1, t_col, 65)
    
    
    ##########        
    # Add metadata
    cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = timepoints)

    # rm(count.matrix)
    
    # names(attributes(cisTopicObject))
    # # [1] "count.matrix"             "binary.count.matrix"     
    # # [3] "is.acc"                   "models"                  
    # # [5] "selected.model"           "log.lik"                 
    # # [7] "dr"                       "calc.params"             
    # # [9] "cell.names"               "cell.data"               
    # # [11] "region.names"             "region.ranges"           
    # # [13] "region.data"              "binarized.cisTopics"     
    # # [15] "signatures"               "binarized.regions.to.Rct"
    # # [17] "binarized.RcisTarget"     "binarized.rGREAT"        
    # # [19] "cistromes.ctx"            "cistromes.regions"       
    # # [21] "cistromes.genes"          "project.name"            
    # # [23] "other"                    "version"                 
    # # [25] "class"  
    # 
    # names(attributes(cisTopicObject)[["cell.data"]])
    # [1] "nCounts"    "nAcc"       "timepoints" 
    
    
    ##########       
    # counts to plot
    countValues2Plot = c('nCounts', "nAcc")

    
    
    ##################################################
    ##################################################
    # save data
    
    saveRDS(cisTopicObject, file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")))
    # cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
    
    
    
    ##################################################
    ##################################################
    # Heatmap 
    # In v0.2.0
    # Error in cellTopicHeatmap(cisTopicObject, method = "probability", distfun = "correlation") : 
    #     Please, install ComplexHeatmap: source("https://bioconductor.org/biocLite.R") 
    # biocLite("ComplexHeatmap")
    
    PDFfilename = file.path(currOutputDir, paste("cisTopic.heatmap.correlation.probability.", nameText, 
                                                 ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                 ".pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    cellTopicHeatmap(cisTopicObject, method='Probability', 
                     colorBy=c('timepoints'),
                     colVars=colsToUse,
                     clustering_distance_rows="spearman") 
    dev.off()
    
    
    
    ##################################################
    ##################################################
    # Plot topics and embeddings
    
    # for (embedMethod in c("PCA", "tSNE", "Umap")) {
    for (embedMethod in c("PCA", "Umap")) {
        
        # embedMethod="PCA"
        # embedMethod="tSNE"
        print(embedMethod)
        
        PDFfilename = file.path(currOutputDir, paste("cisTopic.plotFeatures.cell.", embedMethod, ".counts.metadata.probability.", nameText, 
                                                     ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                     ".pdf", sep=""))
        pdf(PDFfilename, height=10, width=10)
        # Counts
        par(mfrow=c(2,2))
        #plotCellStates(cisTopicObject,  # replaced in v0.2.0
        plotFeatures(cisTopicObject, 
                     method=embedMethod,
                     target='cell', 
                     colorBy=countValues2Plot,
                     topic_contr=NULL)
        # timepoints and genders; females and males
        par(mfrow=c(2,2))
        plotFeatures(cisTopicObject, 
                     method=embedMethod, 
                     target='cell', 
                     colorBy=c('timepoints'),
                     colVars=tcolsToUse,
                     topic_contr=NULL)
       # topics
        par(mfrow=c(2,2))
        plotFeatures(cisTopicObject, 
                     method=embedMethod, 
                     target='cell', 
                     # topics='all',
                     topic_contr='Probability') # topic_contr='Zscore')
        dev.off()
        
    }
    
    
    
    ##################################################
    ##################################################
    # 20190712 
    # Save topic by cell matrix
    
    # Code from cellTopicHeatmap() and modelMatSelection()
    modelMatSelection.gb <- function (object, target, method, all.regions = FALSE) {
        # object = cisTopicObject; target="cell"; method=theMethod
        
        if (length(object@selected.model) < 1) {
            stop("Please, run selectModel() first.")
        }
        if (target == "cell") {
            if (method == "Z-score") {
                modelMat <- scale(object@selected.model$document_expects, 
                                  center = TRUE, scale = TRUE)
            }
            else if (method == "Probability") {
                alpha <- object@calc.params[["runCGSModels"]]$alpha/length(object@selected.model$topic_sums)
                modelMat <- apply(object@selected.model$document_sums, 
                                  2, function(x) {
                                      (x + alpha)/sum(x + alpha)
                                  })
            }
            else {
                stop("Incorrect method selected. Chose method between \"Z-score\" and \"Probability\".")
            }
            colnames(modelMat) <- object@cell.names
            rownames(modelMat) <- paste0("Topic", 1:nrow(modelMat))
        }
        else if (target == "region") {
            if (!all.regions) {
                if (length(object@binarized.cisTopics) < 1) {
                    stop("Please, use binarizecisTopics() first for defining the high confidence regions for dimensionality reduction!")
                }
                else {
                    regions <- unique(unlist(lapply(object@binarized.cisTopics, 
                                                    rownames)))
                }
            }
            topic.mat <- object@selected.model$topics
            if (method == "NormTop") {
                normalizedTopics <- topic.mat/(rowSums(topic.mat) + 
                                                   1e-05)
                modelMat <- apply(normalizedTopics, 2, function(x) x * 
                                      (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
            }
            else if (method == "Z-score") {
                modelMat <- scale(object@selected.model$topics, center = TRUE, 
                                  scale = TRUE)
            }
            else if (method == "Probability") {
                beta <- object@calc.params[["runModels"]]$beta
                topic.mat <- object@selected.model$topics
                modelMat <- (topic.mat + beta)/rowSums(topic.mat + 
                                                           beta)
            }
            else {
                stop("Incorrect method selected. Chose \"NormTop\", \"Z-score\" and \"Probability\".")
            }
            colnames(modelMat) <- object@region.names
            rownames(modelMat) <- paste0("Topic", 1:nrow(modelMat))
            if (!all.regions) {
                modelMat <- modelMat[, regions]
            }
        }
        else {
            stop("Please, provide target=\"cell\" or \"region\".")
        }
        return(modelMat)
    }



    # Code from cellTopicHeatmap() and modelMatSelection()
    for (theMethod in c('Z-score', 'Probability')) {
        # theMethod='Z-score'
        # theMethod='Probability'
        
        # topic.mat <- modelMatSelection(cisTopicObject, "cell", method="")
        #
        # For some reason this doesn't run anymore on server version of cisTopic (although same version (0.2.0) as on my laptop)
        # On server, .modelMatSelection is a hidden method:
        #
        # cisTopic::.modelMatSelection
        # Error: '.modelMatSelection' is not an exported object from 'namespace:cisTopic'

        # names(cisTopicObject@calc.params)
        # 
        # cisTopicObject@selected.model$topics
        # dim(cisTopicObject@selected.model$document_expects)
        #
        # attributes(cisTopicObject)[["cell.data"]]
        #
        # cisTopicObject@region.names
        
        topic.mat <- modelMatSelection.gb(cisTopicObject, "cell", method=theMethod)
        # dim(topic.mat)
   
        rownames(topic.mat) <- paste("Topic", seq(1, nrow(topic.mat)))
        
        # all(colnames(count.matrix) == cisTopicObject@cell.names)
        # NOPE, not always
        # chrX sciRNA dropped some cells for some reason.

        # Not necessary
        # cell_list.annotated <- paste( cisTopicObject@cell.names, timepointList, genderList, sep=".")
        # colnames(topic.mat) = cell_list.annotated
        
        matrixfilename = file.path(currOutputDir, paste(nameText, ".topicsXcells.", theMethod, ".tsv", sep=""))
        write.table(topic.mat, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
        try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
    }

}
