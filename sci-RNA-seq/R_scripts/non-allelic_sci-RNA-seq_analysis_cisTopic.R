##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200325 non-allelic sci-RNA-seq analysis using cisTopic
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

projectName = "non-allelic_sci-RNA-seq""

inputDir = "../data/non-allelic_sci-RNA-seq_count_matrices/"

mainOutputDir = "../results/non-allelic_sci-RNA-seq_cisTopic/"
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
minGenes=1

whichChroms="genome"
# whichChroms="autosomes"
# whichChroms="chrX"
# whichChroms="chr1"










##################################################
##################################################
##################################################
# Load filtered CDS from file

# count.matrix.OG = read.table(gzfile(countsFile), header=TRUE, sep="\t")
count.matrix.OG = as.data.frame(data.table::fread(paste("gzip -dc ", inputCountMatrix, sep=""), header=TRUE, sep="\t"))
dim(count.matrix.OG)
# # 52636  7770
# 52636 10186

head(count.matrix.OG[,1:3])
tail(count.matrix.OG[,1:3])

rownames(count.matrix.OG) = count.matrix.OG[,1]
count.matrix.OG = count.matrix.OG[,-1]

fullGeneDescs = rownames(count.matrix.OG) # Save gene names for later
geneLocations = sub("\\|.*$", "", fullGeneDescs)
dupGeneLocations = duplicated(geneLocations)
sum(dupGeneLocations)
# # 124
# 7

count.matrix.OG = count.matrix.OG[!dupGeneLocations, ]
fullGeneDescs = fullGeneDescs[!dupGeneLocations]
geneLocations = geneLocations[!dupGeneLocations]

rownames(count.matrix.OG) = geneLocations

dim(count.matrix.OG)
# # 52512  7769
# 52629 10185










##################################################
##################################################
##################################################
# Add metadata

# head(pData(cds.OG))
colnames(count.matrix.OG) 

plate2IDs = sub("_.*$", "", colnames(count.matrix.OG) ) # pData(cds.OG)$sample
unique(plate2IDs)
# "d0"           "d3-11andNPCs"

table(plate2IDs)
# d0 d3-11andNPCs 
# 3751         6434 

cellIDs = colnames(count.matrix.OG)  # pData(cds.OG)$Cell
cellCount = length(cellIDs)
# 10185

table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))
# d0 d3-11andNPCs 
# 3751         6434 

sampleSex = sub("^.*\\.", "", colnames(count.matrix.OG))
table(sampleSex)
# female   male 
# 4700   5485 

sampleTimePoints = strsplit(colnames(count.matrix.OG), split="\\.")
sampleTimePoints = sapply(sampleTimePoints, "[", 2 )
table(sampleTimePoints)
# d0  d11   d3   d7  NPC 
# 3751 1672 1721 1797 1244 

table(sampleTimePoints, sampleSex)
# sampleSex
# sampleTimePoints female male
# d0    1943 1808
# d3     960  761
# d7     575 1222
# d11    651 1021
# NPC    571  673

sexList.OG <- factor( sampleSex, levels=c("female", "male") )
timepointList.OG <- factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )

names(sexList.OG) <- colnames(count.matrix.OG) 
names(timepointList.OG) <- colnames(count.matrix.OG) 

# Colors
colsToUse = list()
colsToUse[['timepoints']] = tpCols2

colsToUse[['sex']] = c('red', 'blue')
names(colsToUse[['sex']]) = c('female', 'male')

# Transparent colors
tcolsToUse = list()
tcolsToUse[['timepoints']] = apply(as.matrix(tpCols2), 1, t_col, 75)

tcolsToUse[['sex']] = apply(as.matrix(c('red', 'blue')), 1, t_col, 75)
names(tcolsToUse[['sex']]) = c('female', 'male')

tcolsToUse[['females']] = c(t_col("red", 75), "#00000000")
names(tcolsToUse[['females']]) = c('female', 'male')

tcolsToUse[['males']] = c("#00000000", t_col("blue", 75))
names(tcolsToUse[['males']]) = c('female', 'male')

tcolsToUse[['d0females']] = c(t_col(tpCols2['d0'], 50), "#00000000")
names(tcolsToUse[['d0females']]) = c('d0females', 'other')

tcolsToUse[['d3females']] = c(t_col(tpCols2['d3'], 50), "#00000000")
names(tcolsToUse[['d3females']]) = c('d3females', 'other')

tcolsToUse[['d7females']] = c(t_col(tpCols2['d7'], 50), "#00000000")
names(tcolsToUse[['d7females']]) = c('d7females', 'other')

tcolsToUse[['d11females']] = c(t_col(tpCols2['d11'], 50), "#00000000")
names(tcolsToUse[['d11females']]) = c('d11females', 'other')

tcolsToUse[['d0males']] = c(t_col(tpCols2['d0'], 50), "#00000000")
names(tcolsToUse[['d0males']]) = c('d0males', 'other')

tcolsToUse[['d3males']] = c(t_col(tpCols2['d3'], 50), "#00000000")
names(tcolsToUse[['d3males']]) = c('d3males', 'other')

tcolsToUse[['d7males']] = c(t_col(tpCols2['d7'], 50), "#00000000")
names(tcolsToUse[['d7males']]) = c('d7males', 'other')

tcolsToUse[['d11males']] = c(t_col(tpCols2['d11'], 50), "#00000000")
names(tcolsToUse[['d11males']]) = c('d11males', 'other')










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

generateModels = TRUE
# generateModels = FALSE

numModels = NULL
# numModels = 15

analyzeTopicsOnly = FALSE
# analyzeTopicsOnly = TRUE

if (analyzeTopicsOnly) {
    generateModels = FALSE 
}

theChroms2use = c("genome")
# theChroms2use = c("genome", "autosomes", "chrX", "chr1")
# theChroms2use = c("chrX")

whichChroms="genome"
# whichChroms="autosomes"
# whichChroms="chrX"
# whichChroms="chr1"

whichCells = "all"
# 20210316
# whichCells = "F121only"
# 20210413
# whichCells = "F123only"

for (whichChroms in theChroms2use) {
    print(whichChroms)
    
    ##################################################
    ##################################################
    ##################################################
    # Ouput folder and cisTopic object
    
    # Results folder
    # whichChroms="genome"
    nameText = paste( projectName, 
                      ifelse( whichCells!="all", paste(".", whichCells, sep=""), "" ),
                      paste(".", whichChroms, sep=""),                
                      ifelse(binCountThresh>1, paste("_binCountThresh", binCountThresh, sep=""), ""),
                      ifelse(minCells>1, paste("_minCells", minCells, sep=""), ""),
                      ifelse(minGenes>1, paste("_minGenes", minGenes, sep=""), ""),
                      sep="" ) 
    currOutputDir = file.path( mainOutputDir, paste("cisTopicResults.", nameText, sep="" ))
    # cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
    if ( ! file.exists( currOutputDir ) ) {
        if ( ! dir.create( currOutputDir, showWarnings=FALSE, recursive=TRUE )  ) currOutputDir = mainOutputDir
    }
    
    setwd(currOutputDir)
    getwd()
    
    # save.image() # To .RData
    # load(file = ".RData")
    
    
    
    
    
    
   
   
   
    
    ##################################################
    ##################################################
    ##################################################
    # Preprocess and filter count matrix
    
    locIDs.OG = rownames(count.matrix.OG)
    # These will be the same for count.matrix
    
    # only consider middle 90% of coverage/gene
    acct <- colSums(count.matrix.OG != 0)
    quant = 0.05
    low_lim <- quantile(acct, quant)
    upp_lim <- quantile(acct, 1-quant)
    sum_f = acct[(acct > low_lim)]
    sum_f = as.matrix(sum_f[(sum_f < upp_lim)])
    count.matrix <- count.matrix.OG[,rownames(sum_f)]
    # cell_list <- as.matrix(cell_list.OG[rownames(sum_f)])
    # rownames(count.matrix) = locIDs.OG
    dim(count.matrix)
    # 52629  9160 - sciRNA
    
    # count.matrix[1:5, 1:5]
    
    # batchNumList = batchNumList.OG[rownames(sum_f)]
    sexList = sexList.OG[rownames(sum_f)]
    timepointList = timepointList.OG[rownames(sum_f)]
    table( sexList, timepointList)
    # genome
    # sexList   d0   d3   d7  d11  NPC
    # female 1628  901  523  588  518
    # male   1576  704 1159  959  604










    ##################################################
    ##################################################
    ##################################################    
    # filter for chromosomes
    if (whichChroms != "genome" & !testing) {
        if (whichChroms == "autosomes") {
            locBoolIdx = ! grepl("chrX", locIDs.OG)
        } else {
            locBoolIdx = grepl(paste(whichChroms, ":", sep=""), locIDs.OG)
        }
        count.matrix = count.matrix[locBoolIdx,]
        print(dim(count.matrix))
    }
    # chrX 2613 9160
    # autosomes 50016  9160










    ##################################################
    ##################################################
    ##################################################
    # 20210316
    # filter for sex
    
    cellIDs = colnames(count.matrix)
    
    if (whichCells == "F121only") {
        cellBoolIdx = grepl(paste("\\.female", sep=""), cellIDs)
        count.matrix = count.matrix[,cellBoolIdx]
        print(dim(count.matrix))
        # [1] 52629  4158
        
        sexList = sexList[cellBoolIdx]
        timepointList = timepointList[cellBoolIdx]
        table( sexList, timepointList)
        # sexList    d0   d3   d7  d11  NPC
        # female 1628  901  523  588  518
        # male      0    0    0    0    0

    } else if (whichCells == "F123only") { # 20210413
        cellBoolIdx = grepl(paste("\\.male", sep=""), cellIDs)
        count.matrix = count.matrix[,cellBoolIdx]
        print(dim(count.matrix))
        # [1] 52629  5002
        
        sexList = sexList[cellBoolIdx]
        timepointList = timepointList[cellBoolIdx]
        table( sexList, timepointList)
    #     sexList    d0   d3   d7  d11  NPC
    #     female    0    0    0    0    0
    #     male   1576  704 1159  959  604
    }










    ##################################################
    ##################################################
    ##################################################
    # cisTopic

    if (!generateModels) {
        cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
    } else {
        cisTopicObject <- createcisTopicObject(count.matrix,
                                               project.name=paste(projectName, sep=""), 
                                               keepCountsMatrix = TRUE,
                                               is.acc=binCountThresh,
                                               min.regions=minGenes,
                                               min.cells=minCells)
        
        
        
        ##################################################
        ##################################################
        # runModel
        
        # V3
        cisTopicObject <- runCGSModels(cisTopicObject, topic=c(seq(5,40,5)), seed=999, nCores=4, burnin = 250, iterations = 500)
        # cisTopicObject <- runModels(cisTopicObject, topic=c(seq(50,70,5)), seed=999, nCores=4, burnin = 250, iterations = 500)
        
        
        
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
        dev.off()
        
        
        
        ##################################################
        ##################################################
        # Embeddings
        
        cisTopicObject <- runPCA(cisTopicObject, target='cell')  # specify target in v0.2.0
        # cisTopicObject <- runtSNE(cisTopicObject, seed=999, target='cell')
        cisTopicObject <- runUmap(cisTopicObject, seed=999, target='cell')  # Available from 0.2.0
        
        
        
        ##################################################
        ##################################################
        # save data
        
        saveRDS(cisTopicObject, file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")))
        # cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
        
    }
    
    
    
    # Number of topics
    # nTopics = attributes(cisTopicObject)[["log.lik"]][['topics']][which.max(attributes(cisTopicObject)[["log.lik"]][['LL']])]
    nTopics = 45
    
    cisTopicObject <- selectModel(cisTopicObject, select = nTopics, type = "maximum", 
                                  keepBinaryMatrix = TRUE, keepModels = TRUE)
    
    
    
    if (!analyzeTopicsOnly) {
        
        ##################################################
        ##################################################
        # Setup metadata

        # time-points
        timepoints = as.matrix(timepointList)
        rownames(timepoints) <- colnames(count.matrix)
        colnames(timepoints) <- "timepoints"
        table(timepoints)
        # d0  d11   d3   d7  NPC 
        # 3204 1547 1605 1682 1122 
        
        # 20210316
        # F121 only 
        # d0  d11   d3   d7  NPC 
        # 1628  588  901  523  518 
        
        # 20210413
        # F123 only 
        # d0  d11   d3   d7  NPC 
        # 1576  959  704 1159  604 
        
        # sex
        sex = as.matrix(sexList)
        rownames(sex) <- colnames(count.matrix)
        colnames(sex) <- "sex"
        table(sex)
        # female   male 
        # 4158   5002 
        
        # 20210316
        # F121 only 
        # female 
        # 4158 
        
        # 20210413
        # F123 only 
        # male 
        # 5002 

        table(timepoints, sex)
        # timepoints female male
        # d0    1628 1576
        # d11    588  959
        # d3     901  704
        # d7     523 1159
        # NPC    518  604
        
        # 20210316
        # F121 only 
        # timepoints female
        # d0    1628
        # d11    588
        # d3     901
        # d7     523
        # NPC    518
        
        # 20210413
        # F123 only 
        # timepoints male
        # d0  1576
        # d11  959
        # d3   704
        # d7  1159
        # NPC  604
        
        # Hack to facet females and males
        # females and males
        females = sex
        colnames(females) <- "females"
        males = sex
        colnames(males) <- "males"
        
        # Hack to facet females and males by timepoint
        # females and males
        d0females =  as.matrix(rep("other", length(sex)))
        d0females[timepoints=="d0" & sex=="female"] = "d0females"
        colnames(d0females) <- "d0females"
        rownames(d0females) <- colnames(count.matrix)
        
        d3females = as.matrix(rep("other", length(sex)))
        d3females[timepoints=="d3" & sex=="female"] = "d3females"
        colnames(d3females) <- "d3females"
        rownames(d3females) <- colnames(count.matrix)
        
        d7females = as.matrix(rep("other", length(sex)))
        d7females[timepoints=="d7" & sex=="female"] = "d7females"
        colnames(d7females) <- "d7females"
        rownames(d7females) <- colnames(count.matrix)
        
        d11females = as.matrix(rep("other", length(sex)))
        d11females[timepoints=="d11" & sex=="female"] = "d11females"
        colnames(d11females) <- "d11females"
        rownames(d11females) <- colnames(count.matrix)
        
        d0males = as.matrix(rep("other", length(sex)))
        d0males[timepoints=="d0" & sex=="male"] = "d0males"
        colnames(d0males) <- "d0males"
        rownames(d0males) <- colnames(count.matrix)
        
        d3males = as.matrix(rep("other", length(sex)))
        d3males[timepoints=="d3" & sex=="male"] = "d3males"
        colnames(d3males) <- "d3males"
        rownames(d3males) <- colnames(count.matrix)
        
        d7males = as.matrix(rep("other", length(sex)))
        d7males[timepoints=="d7" & sex=="male"] = "d7males"
        colnames(d7males) <- "d7males"
        rownames(d7males) <- colnames(count.matrix)
        
        d11males = as.matrix(rep("other", length(sex)))
        d11males[timepoints=="d11" & sex=="male"] = "d11males"
        colnames(d11males) <- "d11males"
        rownames(d11males) <- colnames(count.matrix)
        
        # Colors
        colsToUse = list()
        colsToUse[['timepoints']] = tpCols2
        
        colsToUse[['sex']] = c('red', 'blue')
        names(colsToUse[['sex']]) = c('female', 'male')
        
        # Transparent colors
        tcolsToUse = list()
        tcolsToUse[['timepoints']] = apply(as.matrix(tpCols2), 1, t_col, 75)
        
        tcolsToUse[['sex']] = apply(as.matrix(c('red', 'blue')), 1, t_col, 75)
        names(tcolsToUse[['sex']]) = c('female', 'male')
        
        tcolsToUse[['females']] = c(t_col("red", 75), "#00000000")
        names(tcolsToUse[['females']]) = c('female', 'male')
        
        tcolsToUse[['males']] = c("#00000000", t_col("blue", 75))
        names(tcolsToUse[['males']]) = c('female', 'male')
        
        tcolsToUse[['d0females']] = c(t_col(tpCols2['d0'], 50), "#00000000")
        names(tcolsToUse[['d0females']]) = c('d0females', 'other')
        
        tcolsToUse[['d3females']] = c(t_col(tpCols2['d3'], 50), "#00000000")
        names(tcolsToUse[['d3females']]) = c('d3females', 'other')
        
        tcolsToUse[['d7females']] = c(t_col(tpCols2['d7'], 50), "#00000000")
        names(tcolsToUse[['d7females']]) = c('d7females', 'other')
        
        tcolsToUse[['d11females']] = c(t_col(tpCols2['d11'], 50), "#00000000")
        names(tcolsToUse[['d11females']]) = c('d11females', 'other')
        
        tcolsToUse[['d0males']] = c(t_col(tpCols2['d0'], 50), "#00000000")
        names(tcolsToUse[['d0males']]) = c('d0males', 'other')
        
        tcolsToUse[['d3males']] = c(t_col(tpCols2['d3'], 50), "#00000000")
        names(tcolsToUse[['d3males']]) = c('d3males', 'other')
        
        tcolsToUse[['d7males']] = c(t_col(tpCols2['d7'], 50), "#00000000")
        names(tcolsToUse[['d7males']]) = c('d7males', 'other')
        
        tcolsToUse[['d11males']] = c(t_col(tpCols2['d11'], 50), "#00000000")
        names(tcolsToUse[['d11males']]) = c('d11males', 'other')
        
        # counts to plot
        if (whichChroms == "genome") {
            countValues2Plot = c('nCounts', "nAcc")
        } else {
            countValues2Plot = c('nCounts', "nAcc")
        }
        
        # Add metadata
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = timepoints)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = sex)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = females)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = males)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d0females)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d3females)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d7females)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d11females)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d0males)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d3males)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d7males)
        cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = d11males)
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
        # # [1] "nCounts"     "nAcc"           "timepoints" "sex"  
        
        
        
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
                         colorBy=c('sex', 'timepoints'),
                         colVars=colsToUse,
                         clustering_distance_rows="spearman")
        dev.off()
        
        
        
        ##################################################
        ##################################################
        # Plot topics and embeddings
        
        for (embedMethod in c("PCA", "Umap")) {
        # for (embedMethod in c("PCA", "tSNE", "Umap")) {
            
            # if (!testing && embedMethod=="Umap") {
            #     next 
            # }
            
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
            # timepoints and sex; females and males
            par(mfrow=c(2,2))
            plotFeatures(cisTopicObject, 
                         method=embedMethod, 
                         target='cell', 
                         colorBy=c('timepoints', 'sex', 'females', 'males'),
                         colVars=tcolsToUse,
                         topic_contr=NULL)
            # female by timepoint
            par(mfrow=c(2,2))
            plotFeatures(cisTopicObject, 
                         method=embedMethod, 
                         target='cell', 
                         colorBy=c('d0females', 'd3females', 'd7females', 'd11females'),
                         colVars=tcolsToUse,
                         topic_contr=NULL)
            # female by timepoint
            par(mfrow=c(2,2))
            plotFeatures(cisTopicObject, 
                         method=embedMethod, 
                         target='cell', 
                         colorBy=c('d0males', 'd3males', 'd7males', 'd11males'),
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
        
    }
    
    
    
    ##################################################
    ##################################################
    # 20210319
    # Plot topics and embeddings V2 - select topics
    
    selectTopics = list()
    
    if (whichCells == "all" ) {
        
        # F121 + F123
        selectTopics[["genome"]] = c( 52, 9,  # ESCs
                                      18, 47, # female ESCs/d3
                                      60, # male ESCs
                                      40, # female ESC/d3/d7,
                                      4, 3, # female d3
                                      35, # male d3
                                      24, # d7
                                      29, # d11
                                      14, 58, 64 # NPCs - also 5 and 6
        )
        selectTopics[["autosomes"]] = c(27,9,10,30,16,12,6)
        selectTopics[["chrX"]] = c(19,20,26,28,8,18,16,15)
        selectTopics[["chr1"]] = c(15,7,9,14,24,11,20,23)
        
        whichTopics =  c( whichChroms )
        
    } else if (whichCells == "F121only") {
        
        # F121 only 
        selectTopics[["genome"]] = c(
            3,32,         # ESCs
            40,           # ESCs / d3
            4,30,         # d3
            43,           # d3 / d7
            42,45,        # d7/d11 1
            9,25,         # d7/d11 2
            6             # NPCs
        )
        selectTopics[["d7d11NPCs"]] = c(
            23, 25,       # d7
            9, 42, # 45,    # d7/d11 2
            6             # NPCs
        )
        selectTopics[["d3d7d11NPCs"]] = c(
            4,21,30,43, # d3/d7
            42,45,      # d7/d11 1
            9,23,25,    # d7/d11 2
            6           # NPCs
        )
        
        whichTopics =  c( whichChroms, "d7d11NPCs", "select" )
        whichTopics = whichChroms
        # whichTopics = "d7d11NPCs"
        # whichTopics = "select"
        
    } else if (whichCells == "F123only") {
        
        # F123 only 
        selectTopics[["genome"]] = c(
            16, 37,       # ESCs
            30,           # ESCs, d3
            14,           # d3 1
            39,           # d3 2
            26, 35,       # d7
            1,            # d7, d11
            12,           # d11
            3, 22         # NPCs
        )
        
        whichTopics =  c( whichChroms )
        whichTopics = whichChroms
        
    }
    
    
    for ( theTopics in whichTopics ) {
        # theTopics = "d7d11NPCs"
        # theTopics = "d3d7d11NPCs"
    
        for (embedMethod in c("PCA", "Umap")) {
            # for (embedMethod in c("PCA", "tSNE", "Umap")) {
            
            # if (!testing && embedMethod=="Umap") {
            #     next 
            # }
            
            print(embedMethod)
            
            PDFfilename = file.path(currOutputDir, paste("cisTopic.plotFeatures.cell.", embedMethod, ".counts.metadata.probability.", nameText, 
                                                         ifelse(theTopics != whichChroms, paste0(".", theTopics), ""),
                                                         ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                         ".V2.pdf", sep=""))
            pdf(PDFfilename, height=10, width=10)
            # Counts
            par(mfrow=c(2,2))
            #plotCellStates(cisTopicObject,  # replaced in v0.2.0
            # timepoints 
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
                         topics=selectTopics[[theTopics]],
                         topic_contr='Probability') # topic_contr='Zscore')
            dev.off()
            
        }
    }










    ##################################################
    ##################################################
    # 20210317
    # Analysis of the regulatory topics:
    # Defining topics

    # To analyze the regions included in the cisTopics,
    # the first step is always to derive a score that
    # evaluates how likely is for a region to belong to a topic.
    # getRegionsScores() calculates these scores based on
    # the proportion of region specific assignments to a topic.
    # These scores can be rescaled into the range [0,1],
    # which will be useful for the binarization step
    # (as it will force data to follow a gamma distribution shape).
    # This information is stored in the region.data slot.
    cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)


    # binarizecisTopics() allows to select the top regions based on two methods:
    # method = "Predefined": to select a predefined number of regions (determined by the cutoffs argument)
    # method = "GammaFit" (default): to automatically select a threshold based on a fit of the scores to a gamma distribution.
    # This is recommended when using method="NormTop" and scale=TRUE in getRegionScores().
    # Note that the probability threshold must be provided by the user and it must be taken after the density
    # (based on the fitted gamma distribution) is stabilised (i.e. in the tail of the distribution).
    PDFfilename = file.path(currOutputDir, paste("cisTopic.binarizecisTopics.", nameText,
                                                 ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                 ".pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    par(mfrow=c(2,2))
    cisTopicObject <- binarizecisTopics(cisTopicObject, method = "GammaFit", thrP=0.975, plot=TRUE)
    dev.off()

    names(attributes(cisTopicObject)[["binarized.cisTopics"]])

    # The regions sets selected and distributions for each cisTopic can then be analized in different ways (examples in next sections).
    # They can also be exported to bed files to analyze with external tools:
    bedDir = file.path( currOutputDir, paste("cisTopicBeds.", nameText, sep="" ))
    if ( ! file.exists( bedDir ) ) {
        if ( ! dir.create( bedDir, showWarnings=FALSE, recursive=TRUE )  ) bedDir = currOutputDir
    }
    getBedFiles(cisTopicObject, path=bedDir)



    ##################################################
    ##################################################
    # Analysis of the regulatory topics
    # Annotation to genes and GO terms
    
    # BiocManager::install("rGREAT")
    require('rGREAT')
    
    rerun_great = FALSE
    # rerun_great = TRUE

    geneHits=list()
    geneHits[["genome"]] = 5
    # geneHits[["autosomes"]] = 5
    # geneHits[["chrX"]] = 1
    # geneHits[["chr1"]] = 1
    
    sigLevel=list()
    sigLevel[["genome"]] = 0.05
    # sigLevel[["autosomes"]] = 0.01
    # sigLevel[["chrX"]] = 1
    # sigLevel[["chr1"]] = 1
    
    if ( rerun_great ) {
        ##########
        # Run great
        cisTopicObject <- GREAT(cisTopicObject, genome='mm10', fold_enrichment=2, geneHits=geneHits[[whichChroms]], 
                                sign=sigLevel[[whichChroms]], request_interval=30)
        names(attributes(cisTopicObject)[["binarized.rGREAT"]])
    }
    
    
    ##########
    # Plot results
    PDFfilename = file.path(currOutputDir, paste("cisTopic.ontologyDotPlot.", nameText,
                                                 ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                 ".pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    par(mfrow=c(2,2))
    for (theTopic in c(1:nTopics)) {
        ontologyDotPlot(cisTopicObject, top=geneHits[[whichChroms]], topics=theTopic, var.y='name', order.by='Binom_Adjp_BH')
    }
    dev.off()
    
    
    ##########
    # Plot results
    PDFfilename = file.path(currOutputDir, paste("cisTopic.ontologyDotPlot.", nameText,
                                                 ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                 ".orderByBinom_Fold_Enrichment.pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    par(mfrow=c(2,2))
    for (theTopic in c(1:nTopics)) {
        ontologyDotPlot(cisTopicObject, top=geneHits[[whichChroms]], topics=theTopic, var.y='name', order.by='Binom_Fold_Enrichment')
    }
    dev.off()
    
    
    ##########
    # Topics together
    
    for ( theTopics in whichTopics ) {
        # theTopics = "d7d11NPCs"
        # theTopics = "d3d7d11NPCs"
        
        for ( topN in c( 5, 3 ) ) {
            # topN = 3
            
            for ( orderBy in c( 'Binom_Adjp_BH', 'Binom_Fold_Enrichment' ) ) {
                # orderBy='Binom_Fold_Enrichment'
                
                PDFfilename = file.path(currOutputDir, paste("cisTopic.ontologyDotPlot.selectTopics.", nameText, 
                                                             ifelse(theTopics != whichChroms, paste0(".", theTopics), ""),
                                                             ifelse(topN != 5, paste0(".top",  topN), ""),
                                                             ifelse(orderBy != 'Binom_Adjp_BH', paste0(".orderBy",  orderBy), ""),
                                                             ifelse(!is.null(numModels), paste(".numModels=", numModels, sep=""), ""),
                                                             ".pdf", sep=""))
                pdf(PDFfilename, height=8, width=8)
                ontologyDotPlot(cisTopicObject, top=topN, topics=selectTopics[[theTopics]], var.y='name', order.by=orderBy)
                dev.off()
            }
        }
    }
    
    
    
    ##################################################
    ##################################################
    # 20190712 
    # Save topic by cell matrix
    
    # Code from cellTopicHeatmap() and modelMatSelection() on my laptop version of cisTopic
    modelMatSelection.gb <- function (object, target, method, all.regions = FALSE) 
    {
        if (length(object@selected.model) < 1) {
            stop("Please, run selectModel() first.")
        }
        if (target == "cell") {
            if (method == "Z-score") {
                modelMat <- scale(object@selected.model$document_expects, 
                                  center = TRUE, scale = TRUE)
            }
            else if (method == "Probability") {
                alpha <- object@calc.params[["runModels"]]$alpha/length(object@selected.model$topic_sums)
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
        
        topic.mat <- modelMatSelection.gb(cisTopicObject, "cell", method=theMethod)
        dim(topic.mat)
        # 60 6983 -- genome
        
        rownames(topic.mat) <- paste("Topic", seq(1, nrow(topic.mat)))
        all(colnames(count.matrix) == cisTopicObject@cell.names)
        # TRUE
        # batchNumList = batchNumList.OG[rownames(sum_f)]
        # sexList = sexList.OG[rownames(sum_f)]
        # timepointList = timepointList.OG[rownames(sum_f)]
        cell_list.annotated <- paste( cisTopicObject@cell.names, timepointList, sexList, sep=".")
  
        colnames(topic.mat) = cell_list.annotated
        
        matrixfilename = file.path(currOutputDir, paste(nameText, ".topicsXcells.filteredTopBottom", quant, ".", theMethod, ".tsv", sep=""))
        write.table(topic.mat, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
        try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
    }
    
    
    
    ##################################################
    ##################################################
    # save data
    
    saveRDS(cisTopicObject, file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")))
    # cisTopicObject = readRDS( file=file.path(currOutputDir, paste("cisTopicObject.", nameText, ".rds", sep="")) )
    
}
