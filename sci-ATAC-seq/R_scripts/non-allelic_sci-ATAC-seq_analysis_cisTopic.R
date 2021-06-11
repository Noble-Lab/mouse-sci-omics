##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200326 non-allelic sci-ATAC-seq analysis using cisTopic
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

projectName = "non-allelic_sci-ATAC-seq""

inputDir = "../data/non-allelic_sci-ATAC-seq_count_matrices/"

# 20200516
# For consistency, better to use same matrix as used for monocle analysis
usePregeneratedMatrix=TRUE
# usePregeneratedMatrix=FALSE
pregenerated_count_matrix = "../data/non-allelic_sci-ATAC-seq_mouseDiff_Nmasked.countsMatrix.filtered.500UMIs.final.tsv.gz"

mainOutputDir = "../results/non-allelic_sci-ATAC-seq_cisTopic/"
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










##################################################
##################################################
##################################################
# Load barcode plate xref

barcodeInfoFile = file.path(inputDir, "cellIDs.unique.txt")

barcodeInfo = as.data.frame(fread(barcodeInfoFile, header=FALSE, sep="\t"))
dim(barcodeInfo)
# 5586    3

unique(barcodeInfo$cellType)
# "Fd10"  "F_NPC" "Fd4"   "Md10"  "Fd7"   "Md0"   "Md7"  "Fd0"   "Md4" 

colnames(barcodeInfo) = c("barcode", "wellID",  "cellType")

barcodeInfo$timepoint = sub("^_", "", sub("^F", "", sub("^M", "", barcodeInfo$cellType)))
barcodeInfo$sex = "female" 
barcodeInfo$sex[grepl("M", barcodeInfo$cellType)] = "male"

# Fix time points
barcodeInfo$timepoint[barcodeInfo$timepoint=='d4'] = 'd3' 
barcodeInfo$timepoint[barcodeInfo$timepoint=='d10'] = 'd11' 

head(barcodeInfo)
# barcode wellID cellType timepoint sex
# 1 CTGAAGCTCGCTTATCGTGGCGGAGGTCCAGGACGT     7G     Fd10       d11 female
# 2 ATTCAGAACATTATTATTTTAATTAGCCAGGCGAAG     5E    F_NPC       NPC female
# 3 CTGAAGCTAGCCGTAGTTTACGATCATCTAATCTTA     7F      Fd4        d3 female
# 4 TCTCGCGCGATCTTCGCATAGTTCTAGATAATCTTA    11F     Md10       d11   male
# 5 GAGATTCCTATTGCTGGATAGTTCTAGAGTACTGAC     4H      Fd7        d7 female
# 6 CGCTCATTAAGTAGCTCAGTAATGATCGTAATCTTA     3F      Fd4        d3 female










if( !usePregeneratedMatrix ) {
	
	##################################################
	##################################################
	##################################################
	# Load counts data and transform into CDS

	countsFile = file.path(inputDir, "sciATAC_mouseDiff.counts.txt")

	count.table.OG = as.data.frame(fread(countsFile, header=FALSE, sep="\t"))
	dim(count.table.OG)
	# 12,266,315        3

	colnames(count.table.OG) = c("locus", "barcode", "counts")
	head(count.table.OG)










	##################################################
	##################################################
	##################################################
	# Transform table into count matrix - unfiltered

	uniqBarcodes.OG = unique(count.table.OG$barcode)
	length(uniqBarcodes.OG)
	# 5586

	uniqLoci.OG = unique(count.table.OG$locus)
	length(uniqLoci.OG)
	# 157060

	count.matrix.OG = data.table::dcast(count.table.OG, locus ~ barcode, value.var="counts", 
										drop=FALSE, fill=0)
	rownames(count.matrix.OG) = count.matrix.OG[,1]
	count.matrix.OG = count.matrix.OG[,-1]

	dim(count.matrix.OG)
	# 157,060   5,586

	# count.matrix.OG[1:5, 1:5]

	# Rename loci
	locIDs.OG <- rownames(count.matrix.OG) 
	locIDs.OG = sub("_", "-", sub("_", ":", locIDs.OG))
	rownames(count.matrix.OG) = locIDs.OG

	# cell IDs
	cell_list.OG <- colnames(count.matrix.OG)



} else {

    countsFile = file.path(pregenerated_count_matrix)

    count.matrix.OG = as.data.frame(fread(countsFile, header=FALSE, sep="\t"))

    rownames(count.matrix.OG) = count.matrix.OG[,1]
    count.matrix.OG = count.matrix.OG[,-1]
    
    colnames(count.matrix.OG) = count.matrix.OG[1,]
    count.matrix.OG = count.matrix.OG[-1,]
    
    dim(count.matrix.OG)
    # 157060   2727
    
    locIDs.OG = rownames(count.matrix.OG)
    
    cell_list.OG.annotated <- colnames(count.matrix.OG)
    cell_list.OG = sub("\\..*", "", cell_list.OG.annotated)
    
    # Setup cell metadata
    bcIdx = match( cell_list.OG, barcodeInfo$barcode )
    all( cell_list.OG == barcodeInfo$barcode[bcIdx] )
    # TRUE
    
    # time-points
    timepoints.OG = as.matrix(barcodeInfo$timepoint[bcIdx])
    rownames(timepoints.OG) <- colnames(count.matrix.OG)
    colnames(timepoints.OG) <- "timepoints"
    table(timepoints.OG)
    # d0  d11   d3   d7  NPC 
    # 367  464 1098  564  234 
    
    # sex
    sex.OG = as.matrix(barcodeInfo$sex[bcIdx])
    rownames(sex.OG) <- colnames(count.matrix.OG)
    colnames(sex.OG) <- "sex"
    table(sex.OG)
    # female   male 
    # 1810    917 
    
    table(timepoints.OG, sex.OG)
    # timepoints.OG female male
    # d0     227  140
    # d3     945  153
    # d7     107  457
    # d11    297  167
	# NPC    234    0

}










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

theChroms2use = c("genome", "autosomes", "chrX", "chr1")

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
    nameText = paste( projectName, 
                      ifelse( whichCells!="all", paste(".", whichCells, sep=""), "" ),
                      paste(".", whichChroms, sep=""),                
                      ifelse(binCountThresh>1, paste("_binCountThresh", binCountThresh, sep=""), ""),
                      ifelse(minCells>1, paste("_minCells", minCells, sep=""), ""),
                      ifelse(minPeaks>1, paste("_minPeaks", minPeaks, sep=""), ""),
                      sep="" ) 
    currOutputDir = file.path( mainOutputDir, paste("cisTopicResults.", nameText, 
                                                    ifelse(rerunCisTopic.20210048, ".20210048", ""),
                                                    sep="" ))
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
    
    if( !usePregeneratedMatrix ) {
        ##################################################
        ##################################################
        # Filter data
        
        # only consider middle 90% of coverage/peak
        acct <- colSums(count.matrix.OG != 0)
        quant = 0.05
        low_lim <- quantile(acct, quant)
        upp_lim <- quantile(acct, 1-quant)
        sum_f = acct[(acct > low_lim)]
        sum_f = as.matrix(sum_f[(sum_f < upp_lim)])
        count.matrix <- count.matrix.OG[,rownames(sum_f)]
        cell_list <- as.matrix(cell_list.OG[rownames(sum_f)])
        cell_list.annotated <- as.matrix(cell_list.OG.annotated[rownames(sum_f)])
        rownames(count.matrix) = locIDs.OG
        dim(count.matrix)
        # # # 51880  1213 - sciHiC
        # # 7021 2747 - sciATAC: chr12
        # # 156,369   2797 - sciATAC
        # 
        # 157060   4975 - sciATAC N-masked
        # 157060   4975 - sciATAC N-masked - chrX
        
        # count.matrix[1:5, 1:5]
        
        # gc()
        # rm(count.matrix.OG, count.table.OG)
        # gc()
    
    } else { 
    
        count.matrix <- data.matrix( count.matrix.OG )
        cell_list <- as.matrix(cell_list.OG)
        cell_list.annotated <- as.matrix(cell_list.OG.annotated)
        #  rownames(count.matrix) = locIDs.OG
        
    }
    
    
    
    
    
    
    
    
    
        
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
    # chrX 4364 2727 - filtered
    # autosomes 152,696   2727 
    
    
    
    
    
    
    
    
    
    
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
        # 157060   1810
        
        cell_list <- cell_list[cellBoolIdx]
        cell_list.annotated <- cell_list.annotated[cellBoolIdx]
    
    } else if (whichCells == "F123only") { # 20210413
        cellBoolIdx = grepl(paste("\\.male", sep=""), cellIDs)
        count.matrix = count.matrix[,cellBoolIdx]
        print(dim(count.matrix))
        # 157060    917
        
        cell_list <- cell_list[cellBoolIdx]
        cell_list.annotated <- cell_list.annotated[cellBoolIdx]
    } 
    
    
    
    
    
    
    
    
    
    
    ##################################################
    ##################################################
    ##################################################
    # Setup metadata
    
    # # Could also do this:
    # timepoints <- as.matrix(timepoints.OG[rownames(sum_f), 1]) 
    # sex <- as.matrix(sex.OG[rownames(sum_f), 1])
    
    # Must reset the bcIndx for filtered matrix
    bcIdx = match( cell_list, barcodeInfo$barcode )
    all( cell_list == barcodeInfo$barcode[bcIdx] )
    # TRUE
    
    # time-points
    timepoints = as.matrix(barcodeInfo$timepoint[bcIdx])
    rownames(timepoints) <- cell_list.annotated
    colnames(timepoints) <- "timepoints"
    table(timepoints)
    # 200 UMIs
    # d0  d11   d3   d7  NPC
    # 367  464 1098  564  234
    
    # 20210316
    # F121only
    # d0 d11  d3  d7 NPC 
    # 227 297 945 107 234
    
    # 20210413
    # F123only
    # d0 d11  d3  d7 
    # 140 167 153 457 
    
    
    # sex
    sex = as.matrix(barcodeInfo$sex[bcIdx])
    rownames(sex) <- cell_list.annotated
    colnames(sex) <- "sex"
    table(sex)
    # female   male 
    # 1810    917 
    
    # 20210316
    # F121only
    # female 
    # 1810 
    
    # 20210413
    # F123only
    # male 
    # 917 
    
    table(timepoints, sex)
    # timepoints female male
    # d0     227  140
    # d11    297  167
    # d3     945  153
    # d7     107  457
    # NPC    234    0
    
    # 20210316
    # F121only
    # timepoints female
    # d0     227
    # d11    297
    # d3     945
    # d7     107
    # NPC    234
    
    # 20210413
    # F123only
    # timepoints male
    # d0   140
    # d11  167
    # d3   153
    # d7   457
    
    
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
    rownames(d0females) <- cell_list.annotated
    
    d3females = as.matrix(rep("other", length(sex)))
    d3females[timepoints=="d3" & sex=="female"] = "d3females"
    colnames(d3females) <- "d3females"
    rownames(d3females) <- cell_list.annotated
    
    d7females = as.matrix(rep("other", length(sex)))
    d7females[timepoints=="d7" & sex=="female"] = "d7females"
    colnames(d7females) <- "d7females"
    rownames(d7females) <- cell_list.annotated
    
    d11females = as.matrix(rep("other", length(sex)))
    d11females[timepoints=="d11" & sex=="female"] = "d11females"
    colnames(d11females) <- "d11females"
    rownames(d11females) <- cell_list.annotated
    
    d0males = as.matrix(rep("other", length(sex)))
    d0males[timepoints=="d0" & sex=="male"] = "d0males"
    colnames(d0males) <- "d0males"
    rownames(d0males) <- cell_list.annotated
    
    d3males = as.matrix(rep("other", length(sex)))
    d3males[timepoints=="d3" & sex=="male"] = "d3males"
    colnames(d3males) <- "d3males"
    rownames(d3males) <- cell_list.annotated
    
    d7males = as.matrix(rep("other", length(sex)))
    d7males[timepoints=="d7" & sex=="male"] = "d7males"
    colnames(d7males) <- "d7males"
    rownames(d7males) <- cell_list.annotated
    
    d11males = as.matrix(rep("other", length(sex)))
    d11males[timepoints=="d11" & sex=="male"] = "d11males"
    colnames(d11males) <- "d11males"
    rownames(d11males) <- cell_list.annotated
    
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

    # save.image() # To .RData
    # load(file = ".RData")
    
    
    
    
    
    
    
    
    
    
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
                                               min.regions=minPeaks,
                                               min.cells=minCells)
        
        
        
        ##################################################
        ##################################################
        # runModel
        
        # 20210408
        # Rerun just for optimal 25 topics b/c getting following error when running rGreat:
        # Error: subscript is a NSBS object that is incompatible with the currentsubsetting operation
        if(rerunCisTopic.20210048) {
            cisTopicObject.20210408 <- runCGSModels(cisTopicObject,  topic=c(25), seed=999, nCores=4, burnin = 250, iterations = 500)
            cisTopicObject.OG = cisTopicObject; cisTopicObject = cisTopicObject.20210408
        } else {
            # V3
            cisTopicObject <- runCGSModels(cisTopicObject,  topic=c(seq(5,40,5)), seed=999, nCores=4, burnin = 250, iterations = 500)
            # cisTopicObject <- runModels(cisTopicObject,  topic=c(seq(15,45,5)), seed=999, nCores=4, burnin = 250, iterations = 500)
        }
        
        
        
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
        
        # save.image() # To .RData
        # load(file = ".RData")
    }
    
    
    
    # Number of topics
    nTopics = attributes(cisTopicObject)[["log.lik"]][['topics']][which.max(attributes(cisTopicObject)[["log.lik"]][['LL']])]
    # nTopics = 20
    
    cisTopicObject <- selectModel(cisTopicObject, select = nTopics, type = "maximum", 
                                  keepBinaryMatrix = TRUE, keepModels = TRUE)
     
    
    
    if (!analyzeTopicsOnly) {
        
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
                         clustering_distance_rows="spearman") + 
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
            # male by timepoint
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
    if(rerunCisTopic.20210048) {
        
        if (whichCells == "all" ) {
            # F121 + F123
            selectTopics[["genome"]] = c(6, 3,              # ESCs
                                         15, 17, 10,        # ESCs/d3
                                         21,                # ESCs/d3/d7
                                         1, 12, 24,         # d3
                                         14, 20,            # d3/d7
                                         2, 9,              # female d7/d11
                                         5, 19,             # male d7/d11
                                         4,                 # NPCs
                                         25, 18             # All
            )
            
            whichTopics =  c( whichChroms )
            
        }
        
    } else { 
        
        if (whichCells == "all" ) {
            
            # F121 + F123
            selectTopics[["genome"]] = c(13, 23, 18,       # ESCs
                                         4, 11,            # d3
                                         6, 10, 25, 7, 19, #d7/d11
                                         22,               # NPCs
                                         5, 7, 20          # All
            )
            selectTopics[["autosomes"]] = c(27,9,10,30,16,12,6)
            selectTopics[["chrX"]] = c(15, 16, 17, 21, 20, 23, # ESCs and d3
                                       9, 10,                  # d7, d11
                                       24                      # NPCS
            )
            selectTopics[["chr1"]] = c(15,7,9,14,24,11,20,23)
            
            whichTopics =  c( whichChroms )
            
        } else if (whichCells == "F121only" ) {
            
            # F121 only
            selectTopics[["genome"]] = c(
                7, 20,      # ESCs
                5, 6,       # ESCs/d3
                11, 14, 16, # d3
                18, 3, 19,  # d7/d11
                8           # NPCs
            )
            selectTopics[["d7d11NPCs"]] = c(
                3, 12,      # d7
                18, 19,     # d11
                8           # NPCs
            )
            selectTopics[["d3d7d11NPCs"]] = c(
                14,        # d3
                18, 3, 19, # d11/d7
                8          # NPCs
            )
            selectTopics[["d3d7d11NPCs.V2"]] = c(
                4,         # ESC/d3
                14, 11,    # d3
                18, 3, 19, # d11/d7
                8          # NPCs
            )
            
            whichTopics =  c( whichChroms, "d7d11NPCs", "d3d7d11NPCs", "d3d7d11NPCs.V2" )
            
        } else if (whichCells == "F123only" ) {
            
            # F123 only
            selectTopics[["genome"]] = c(
                9, 12,     # ESCs
                8, 13,     # ESCs/d3
                1, 10, 11, # d3
                15,        # d7 1
                6,         # d7 2
                5          # d11
            )
            
            whichTopics =  c( whichChroms )
        }
        
    }
    
    theTopics = whichChroms
    
    
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
    # geneHits[["chrX"]] = 2
    # geneHits[["chr1"]] = 2
 
    sigLevel=list()
    sigLevel[["genome"]] = 0.05
    # sigLevel[["autosomes"]] = 0.01
    # sigLevel[["chrX"]] = 1
    # sigLevel[["chr1"]] = 1
    
   
    
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
    
    for ( theTopics in whichTopics) {
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
                pdf(PDFfilename, height=8, width=16)
                ontologyDotPlot(cisTopicObject, top=topN, topics=selectTopics[[theTopics]], var.y='name', order.by=orderBy)
                dev.off()
            }
        }
    }
    
    
    
    ##################################################
    ##################################################
    # 20190712 
    # save topic by cell matrix
    
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
        # 25 4975 -- genome
        
        rownames(topic.mat) <- paste("Topic", seq(1, nrow(topic.mat)))
        colnames(topic.mat) <- cisTopicObject@cell.names

        # Setup cell metadata
        bcIdx = match( colnames(topic.mat), barcodeInfo$barcode )
        all( colnames(topic.mat) == barcodeInfo$barcode[bcIdx] )
        # TRUE
        
        # time-points
        timepoints.OG = as.matrix(barcodeInfo$timepoint[bcIdx])
        rownames(timepoints.OG) <- colnames(topic.mat)
        colnames(timepoints.OG) <- "timepoints"
        table(timepoints.OG)
        # d0  d11   d3   d7  NPC 
        # 574 1059 1811 1174  357
        
        # sex
        sex.OG = as.matrix(barcodeInfo$sex[bcIdx])
        rownames(sex.OG) <- colnames(topic.mat)
        colnames(sex.OG) <- "sex"
        table(sex.OG)
        # female   male 
        # 3147   1828 
        
        table(timepoints.OG, sex.OG)
        # timepoints.OG female male
        # d0     344  230
        # d11    637  422
        # d3    1447  364
        # d7     362  812
        # NPC    357    0
        
        # cell_list.OG = sub("\\..*", "", cell_list.OG)
        cell_list.OG.annotated = paste(colnames(topic.mat), timepoints.OG, sex.OG, sep=".")
        
        colnames(topic.mat) = cell_list.OG.annotated
        
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


