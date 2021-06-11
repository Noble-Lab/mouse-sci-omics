##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200203 allelic sciHic contact decay analysis
#
##################################################
##################################################
##################################################



##################################################
# Set up
# 
# R/3.5.0
#
##################################################



rm(list=ls())
gc()
sessionInfo()
    









##################################################
##################################################
##################################################
# Packages and functions

require(viridis)

require(corrplot)

# install.packages("factoextra")
require(factoextra)

#  Barnes-Hut-SNE 
require(Rtsne)

require(smoother)
require(binhf)

require(gplots)

# require(dynamicTreeCut)

require(scatterplot3d)
# http://www.sthda.com/english/wiki/scatterplot3d-3d-graphics-r-software-and-data-visualization
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')



##################################################
# Functions

# https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# Load latest version of heatmap.3 function
source("heatmap.3.GB.R")

prettyBPs <- function(x) {
  cx = paste(as.character(round(x,1)), "bp", sep="")
  cx[x>=10^3] = paste(as.character(round(x[x>=10^3]/10^3,2)), "kb", sep="")
  cx[x>=10^6] = paste(as.character(round(x[x>=10^6]/10^6,1)), "Mb", sep="")
  return(cx)
}

setAlpha <- function(col, alpha=1) {
  if(missing(col))  stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))
}

getProportions <- function(x) {
  totx = sum(x, na.rm=TRUE)
  return(x/totx)
}

addNZval <- function(x, d=0.001) {
  return(x + (d * min(x[x!=0], na.rm=TRUE)))
}

# https://gist.github.com/Jfortin1/72ef064469d1703c6b30
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

# library(pdist)
euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

binarize_vector <- function(x){
  meanx = mean(x, na.rm=TRUE)
  nux = x
  nux[ x > meanx ] = 1
  nux[ x <= meanx ] = 0
  return(nux)
}










##################################################
##################################################
##################################################
# Inputs and options

projectName = "allelic_sci-Hi-C"

inputDir = "../data/allelic_contact_decay_profiles/"

# 20200803
patski100=FALSE
# patski100=TRUE

# 20200921 - Needed to generate genome-wide allelic CDPs for TsixStop libs
F121sOnly=FALSE
# F121sOnly=TRUE

# 20200921 - Needed to generate genome-wide allelic CDPs for TsixStop libs
TsixStopOnly=FALSE
# TsixStopOnly=TRUE

# 20200921 - Needed to generate genome-wide allelic CDPs for F123 libs
# *** NOTE: Only chr1 used to check for nThresh cells, not both chr1 and chrX
F123sOnly=FALSE
# F123sOnly=TRUE

# 20200922 - Needed to generate genome-wide allelic CDPs for F121 amd F123 libs
# *** NOTE: Only chr1 used to check for nThresh cells, not both chr1 and chrX
# ***       So saved in a seperate folder.
F121andF123=FALSE
# F121andF123=TRUE

# 20210323
oldF121andTsixStop=FALSE
# oldF121andTsixStop=TRUE

# 20210417
F121andTsixStop=FALSE
# F121andTsixStop=TRUE

# 20210304
plotDecayCurvesByCluster=FALSE
# plotDecayCurvesByCluster=TRUE
#
if (plotDecayCurvesByCluster) {
  # Only for differentiated cells.
  NPCsByCluster=TRUE; ESCsByCluster=FALSE
  # And pluripotent cells.
  # ESCsByCluster=TRUE; NPCsByCluster=FALSE
} else {
  NPCsByCluster=FALSE; ESCsByCluster=FALSE
}

# 20210311
plotF121DecayCurvesByClusterJointly=FALSE
# plotF121DecayCurvesByClusterJointly=TRUE
if (plotF121DecayCurvesByClusterJointly) {
  F121sOnly=TRUE
}

# 20210316
plotDecayCurvesAndMapsByXCIstatus=FALSE
# plotDecayCurvesAndMapsByXCIstatus=TRUE

# if ( sum(patski100, F121sOnly, TsixStopOnly, F123sOnly) > 1 ) exit()

supermainOutputDir = "../results/allelic_sci-Hi-C/"
#
# 20200803
# All other analysis was performed using an nThresh=50 so this needs it's own folder:
if(patski100) supermainOutputDir = paste0( supermainOutputDir, ".patski100" )
#
# 20200922 - Needed to generate genome-wide allelic CDPs for F121 amd F123 libs
# ***       So saved in a seperate folder.
if(F121andF123)  supermainOutputDir = paste0( supermainOutputDir, ".F121andF123" )
#
if(plotDecayCurvesByCluster || plotF121DecayCurvesByClusterJointly)  supermainOutputDir = paste0( supermainOutputDir, ".DecayCurvesByCluster" )

if ( ! file.exists( supermainOutputDir ) ) {
  if ( ! dir.create( supermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) supermainOutputDir = getwd()
}
setwd(supermainOutputDir)
getwd()

if ( F121sOnly ) {
    
    libs=c(
      "F121_d0",
      "F121_d3",
      "F121_d7",
      "F121_d11",
      "F121_NPC"
    )
    
    libCols=c(
      # F121
      "red",
      "orange",
      "green",
      "blue",
      "black"
    )
    
    names(libCols) = libs
    
    nThreshes = c(50)  
    
    # chroms2analyze = c("chrX", "chr1", "chr4", "chr19", "all")
    # chroms2analyze = c( paste( "chr", c( 1, "X", 4, 19), sep="" ) )
    # chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
    chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )
    
    # 20210311
    if (plotF121DecayCurvesByClusterJointly) {
      chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
    }

} else if (TsixStopOnly) {
  # 20200921
  
  libs = c("TsixStop_d0",
           "TsixStop_d3",
           "TsixStop_d7",
           "TsixStop_d11",
           "TsixStop_NPC")
  
  libCols = c(
    # Tsix
    "red",
    "orange",
    "green",
    "blue",
    "black")
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
  # chroms2analyze = c("chrX", "chr1", "chr4", "chr19", "all")
  # chroms2analyze = c( paste( "chr", c( 1, "X", 4, 19), sep="" ) )
  # chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )

} else if ( F121andTsixStop ) {
  # 20210414
  
  libs=c(
    "F121_d0",
    "F121_d3",
    "F121_d7",
    "F121_d11",
    "F121_NPC",
    "TsixStop_d0",
    "TsixStop_d3",
    "TsixStop_d7",
    "TsixStop_d11",
    "TsixStop_NPC"
  )
  
  libCols=c(
    # F121
    "red",
    "orange",
    "green",
    "blue",
    "black",
    # Tsix
    "red",
    "orange",
    "green",
    "blue",
    "black"
  )
  
  names(libCols) = libs
  
  nThreshes = c(50)  
  
  # chroms2analyze = c("chrX", "chr1", "chr4", "chr19", "all")
  # chroms2analyze = c( paste( "chr", c( 1, "X", 4, 19), sep="" ) )
  chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  # chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )

} else if (F123sOnly) { 
  # 20200921
  # *** NOT: Only chr1 used to check for nThresh cells, not both chr1 and chrX

  libs=c( "F123_d0", "F123_NPC" )
  libCols=c( "red", "black" )
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
  # chroms2analyze = c("chrX", "chr1", "chr4", "chr19", "all")
  # chroms2analyze = c( paste( "chr", c( 1, "X", 4, 19), sep="" ) )
  # chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )
  
} else if (F121andF123) { 
  # 20200922
  # *** NOT: Only chr1 used to check for nThresh cells, not both chr1 and chrX
  
  libs=c(
    "F121_d0",
    "F121_d3",
    "F121_d7",
    "F121_d11",
    "F121_NPC",
    "F123_d0", 
    "F123_NPC" )
  
  libCols=c(
    # F121
    "red",
    "orange",
    "green",
    "blue",
    "black",
    "red", 
    "black" )
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
} else if (patski100)  {
  
  # 20200803 - Patski min 100 contacts LMD and Spearman analysis 
  libs=c("Patski")
  libCols=c("gray")
  names(libCols) = libs

  nThreshes = c(100)
  
  chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  
} else if (plotDecayCurvesByCluster)  {
  
  # 20210304
  if (NPCsByCluster) {
    libs=c(
      # "Patski",
      # 
      # "EBdiff_F121_d20",
      # "EBdiff_TsixStop_d20", 
      
      "TsixStop_NPC",
      "F121_NPC"
    )
    
    libCols=c(
      # "gray",
      # 
      # # oldTsix,
      # "darkblue",
      # "darkblue",
      
      # Tsix
      "black",
      
      # F121
      "black"
      
    )
  }
  
  # 20210305
  if (ESCsByCluster) {
    libs=c(
      # "EBdiff_F121_d0",
      # "EBdiff_TsixStop_d0", 
      
      "TsixStop_d0",
      "F121_d0"
    )
    
    libCols=c(
      # # oldTsix,
      # "red",
      # "red",
      
      # Tsix
      "red",
      
      # F121
      "red"
      
    )
  } 
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
  chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  
} else if ( oldF121andTsixStop ) {
  # 20210323
  
  libs=c(
    "EBdiff_F121_d0",
    "EBdiff_F121_d7",
    "EBdiff_F121_d20",
    # "EBdiff_F123_d0",
    # "EBdiff_F123_d7",
    # "EBdiff_F123_d20",
    "EBdiff_TsixStop_d0",
    "EBdiff_TsixStop_d7",
    "EBdiff_TsixStop_d20"
  )
  
  # Use this to color Patski libraries instead.
  libCols=c(
    # oldTsix,
    "red",
    "orange",
    "darkblue",
    # "red",
    # "orange",
    # "darkblue",
    "red",
    "orange",
    "darkblue"
  )
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
  # chroms2analyze = c("chrX", "chr1", "chr4", "chr19", "all")
  # chroms2analyze = c( paste( "chr", c( 1, "X", 4, 19), sep="" ) )
  # chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )
  
}  else {
  libs=c(
    "Patski",
    
    # # "mES_EBdiff_D0",
    # # "mES_EBdiff_D7",
    # # "mES_EBdiff_D20",
    
    "EBdiff_F121_d0",
    "EBdiff_F121_d7",
    "EBdiff_F121_d20",
    # "EBdiff_F123_d0",
    # "EBdiff_F123_d7",
    # "EBdiff_F123_d20",
    "EBdiff_TsixStop_d0",
    "EBdiff_TsixStop_d7",
    "EBdiff_TsixStop_d20",
    
    "TsixStop_d0",
    "TsixStop_d3",
    "TsixStop_d7",
    "TsixStop_d11",
    "TsixStop_NPC",
    
    "F121_d0_rep1",
    "F121_d0_rep2",
    
    "F121_d0",
    "F121_d3",
    "F121_d7",
    "F121_d11",
    "F121_NPC",
    
    "F123_d0",
    "F123_NPC"
  )
  
  # Use this to color Patski libraries instead.
  libCols=c(
    "gray",
    
    # oldTsix,
    "red",
    "orange",
    "darkblue",
    # "red",
    # "orange",
    # "darkblue",
    "red",
    "orange",
    "darkblue",
    
    # Tsix
    "red",
    "orange",
    "green",
    "blue",
    "black",
    
    # F121_reps
    "orange",
    "darkred",
    
    # F121
    "red",
    "orange",
    "green",
    "blue",
    "black",
    
    # F123
    "red",
    "black"
  )
  
  names(libCols) = libs
  
  nThreshes = c(50)
  
  # chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
  chroms2analyze = c( paste( "chr", c( c(1:19), "X" ), sep="" ) )
  
}



# 20210316
if( plotDecayCurvesAndMapsByXCIstatus ) {
  chroms2analyze = c( paste( "chr", c( 1, "X"), sep="" ) )
}
  
# alleleShapes=c(16, 17)
alleleShapes=c(21, 24)
names(alleleShapes) = c("ref", "alt")

alleleCols=c("darkred", "darkblue")
names(alleleCols) = c("ref", "alt")

# In-silico cell phasing over the cell cycle
# We counted the number of cis (intra-chromosomal) contacts per cell, binning contacts by distance into logarithmic bins 
# (143 bins, first one for contacts distanced < 1Kb, then each bin covers an exponent step of 0.125, using base 2). 
# log2(1000) ~ 9.96 ~10 
# length(2^seq(10, 27.625, 0.125)/(10^6)) + 1 (trans) = 143
# 2^27.625 = 206,992,033
# 2^9.875 = 939.0121
#
# Actually saved data to 2^27.750 = 225,726,413 Mb, or 143 bins + 1 trans = 144
# binDist = c( 2^seq(10, 27.750, 0.125), Inf)
# length(binDist)
# # 144
binDist = c( 2^seq(10, 27.750, 0.125))
length(binDist)
# 143

heatCols1 = heat.colors(256)
heatCols2 <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                    "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                    "#4393C3", "#2166AC", "#053061"))(256))
heatColsM = viridis(100, option="magma")
heatColsV = viridis(100, option="D")










# ---------------------
# 20210304 

if (plotDecayCurvesByCluster) {
  
  #######################
  # Load kmeans group IDs
  # For plotting allelic CDPs by cell cycle state (based on autosomal CDP clusters).
  # And removing mitotic cells (based on autosomal CDP clusters).
  
  otherProjectName="non-allelic_sci-Hi-C"
  
  clusteringOutputDir = file.path("../data/autosomal_non-allelic_cell_cycle_cluster_ids/")
  
  nThresh=1000 # min total contact counts per cell

  bins2Use="50kb-8Mb"
  # bins2Use="1kb-80Mb"
  
  useScaled=TRUE
  # useMethod="euclidean"
  useMethod="spearman"
  
  kCenters = 4
  # kCenters = 8
  
  autosomalCDPs=TRUE
  
  kClusters = list()
  
  kclustFiledir = file.path(clusteringOutputDir, 
  	paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), 
  	bins2Use, paste( "kclust", ".k=", kCenters, sep="" ), "clusterIDs")
  if ( ! file.exists( kclustFiledir ) ) {
    if ( ! dir.create( kclustFiledir, showWarnings=FALSE, recursive=TRUE )  ) kclustFiledir = kclustMaindir
  }

  kClustDesc = paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""),
                      ".", bins2Use,
                      ".", ifelse(useScaled, "scaled", "unscaled"),
                      ".", paste("kclust.dist", useMethod, sep=""),
                      ".", paste("kCenters", kCenters, sep=""),
                      sep="")
  
  for (alib in libs) {
    
    nameText = paste(alib, kClustDesc, sep=".")
    kclustFileName = file.path(kclustFiledir, paste(otherProjectName, ".", nameText, ".tsv", sep=""))
    kClusters[[alib]] = read.table(file=kclustFileName, sep="\t", header=FALSE)
    colnames(kClusters[[alib]]) = c("cellID", "clusterID")
    print( alib )
    print( table(kClusters[[alib]]$clusterID) )
    
    # # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
    # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    # # kClustColors =  rainbow(kCenters)
    # set.seed(1234)
    # kClustColors =  sample(color, kCenters)
    kClustColors = colorRampPalette(c("#00ffff", "#001e1e"))(kCenters)
    names(kClustColors) = paste(1:kCenters)
    
    # To add transparency to colors
    alphaLvl = 0.25
    txparentColsTable = cbind(t(col2rgb(kClustColors, alpha=FALSE)), alphaLvl*255)/255
    colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
    txparentCols4Dists = c()
    for (colRowIdx in 1:nrow(txparentColsTable)) {
      colRowAsList = as.list(txparentColsTable[colRowIdx,])
      txparentCols4Dists = c(txparentCols4Dists, do.call(rgb, colRowAsList))
    }
    alphaLvl = 0.75
    txparentColsTable = cbind(t(col2rgb(kClustColors, alpha=FALSE)), alphaLvl*255)/255
    colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
    txparentCols4Strip = c()
    for (colRowIdx in 1:nrow(txparentColsTable)) {
      colRowAsList = as.list(txparentColsTable[colRowIdx,])
      txparentCols4Strip = c(txparentCols4Strip, do.call(rgb, colRowAsList))
    }
    
  } # alib
  
}
  
# [1] "F121_d0"
# 1    2    3    4 
# 93   58  411 1641 

# [1] "F121_d3"
# 1   2   3   4 
# 85 143 311 355 

# [1] "F121_d7"
# 1   2   3   4 
# 80 143 260 249 

# [1] "F121_d11"
# 1   2   3   4 
# 33  66 170 339 

# [1] "F121_NPC"
# 1   2   3   4 
# 28  48 432 546 

# [1] "TsixStop_d0"
# 1   2   3   4 
# 14   8 108 212 

# [1] "TsixStop_d3"
# 1   2   3   4 
# 66 100 183 490 

# [1] "TsixStop_d7"
# 1   2   3   4 
# 22  86 214 500 

# [1] "TsixStop_d11"
# 1   2   3   4 
# 18  60 150 243 

# [1] "TsixStop_NPC"
# 1    2    3    4 
# 70   68  622 1775 

# 20210304 
# ---------------------










# ---------------------
# 20210311

if (plotF121DecayCurvesByClusterJointly) {
  
  for (combined_label in c("F121_combined", "TsixStop_combined") ) {
    
    #######################
    # Load kmeans group IDs
    # For plotting allelic CDPs by cell cycle state (based on autosomal CDP clusters).
    # And removing mitotic cells (based on autosomal CDP clusters).
    
    otherProjectName="non-allelic_sci-Hi-C"
    
    clusteringOutputDir = file.path("../data/autosomal_non-allelic_cell_cycle_cluster_ids/")
    
    nThresh=1000 # min total contact counts per cell

    bins2Use="50kb-8Mb"
    # bins2Use="1kb-80Mb"
  
    useScaled=TRUE
    # useMethod="euclidean"
    useMethod="spearman"
  
    kCenters = 4
    # kCenters = 8
    
    autosomalCDPs=TRUE
    
    kClusters = list()
    
    kclustFiledir = file.path(clusteringOutputDir, 
    	paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), 
    	bins2Use, paste( "kclust", ".k=", kCenters, sep="" ), "clusterIDs")
    if ( ! file.exists( kclustFiledir ) ) {
      if ( ! dir.create( kclustFiledir, showWarnings=FALSE, recursive=TRUE )  ) kclustFiledir = kclustMaindir
    }
    
    kClustDesc = paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""),
                        ".", bins2Use,
                        ".", ifelse(useScaled, "scaled", "unscaled"),
                        ".", paste("kclust.dist", useMethod, sep=""),
                        ".", paste("kCenters", kCenters, sep=""),
                        sep="")
    
    for (alib in combined_label) {
      
      nameText = paste(alib, kClustDesc, sep=".")
      kclustFileName = file.path(kclustFiledir, paste(otherProjectName, ".", nameText, ".tsv", sep=""))
      kClusters[[alib]] = read.table(file=kclustFileName, sep="\t", header=FALSE)
      colnames(kClusters[[alib]]) = c("cellID", "clusterID")
      print( alib )
      print( table(kClusters[[alib]]$clusterID) )
      
      # # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
      # color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      # # kClustColors =  rainbow(kCenters)
      # set.seed(1234)
      # kClustColors =  sample(color, kCenters)
      kClustColors = colorRampPalette(c("#00ffff", "#001e1e"))(kCenters)
      names(kClustColors) = paste(1:kCenters)
      
      # To add transparency to colors
      alphaLvl = 0.25
      txparentColsTable = cbind(t(col2rgb(kClustColors, alpha=FALSE)), alphaLvl*255)/255
      colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
      txparentCols4Dists = c()
      for (colRowIdx in 1:nrow(txparentColsTable)) {
        colRowAsList = as.list(txparentColsTable[colRowIdx,])
        txparentCols4Dists = c(txparentCols4Dists, do.call(rgb, colRowAsList))
      }
      alphaLvl = 0.75
      txparentColsTable = cbind(t(col2rgb(kClustColors, alpha=FALSE)), alphaLvl*255)/255
      colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
      txparentCols4Strip = c()
      for (colRowIdx in 1:nrow(txparentColsTable)) {
        colRowAsList = as.list(txparentColsTable[colRowIdx,])
        txparentCols4Strip = c(txparentCols4Strip, do.call(rgb, colRowAsList))
      }
      
    } # alib
    
    # 1    2    3    4 
    # 281  357 1060 3793 
    
  }
  
}
  
# 20210311
# ---------------------










##################################################
##################################################
##################################################
# 20200205 
# Multichromosome analysis: load all the data

matrixFileName = list()
matrixData = list()

for (whichChr in chroms2analyze) {
  # 20200203
  # whichChr=""
  # whichChr="chr1"
  # whichChr="chrX"
  
  cat("\n\n\n----------\n----------\n", whichChr, "\n", sep="")
  
  ##################################################
  # Load decay data
  
  matrixData[[whichChr]] = list()
  
  for(lib in libs) {
    for(allele in c("ref", "alt")) {
      alleliclib=paste(lib, allele, sep=".")
      # Counts file
      matrixFileName[[alleliclib]] <- paste(inputDir,
                                         "allelicContactDecayData", ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), "/",
                                         alleliclib, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), ".logbinContactDecay.mm10.txt.gz", sep="")
      matrixData[[whichChr]][[alleliclib]] <- read.table(gzfile(matrixFileName[[alleliclib]]), header=TRUE)
      rownames(matrixData[[whichChr]][[alleliclib]]) = matrixData[[whichChr]][[alleliclib]][,1]
      matrixData[[whichChr]][[alleliclib]] = matrixData[[whichChr]][[alleliclib]][,-1] # Remove bin ID
      colnames(matrixData[[whichChr]][[alleliclib]]) = sub("\\.", "-", colnames(matrixData[[whichChr]][[alleliclib]]))
      cat(alleliclib, "\t", paste(dim(matrixData[[whichChr]][[alleliclib]]), collapse="\t"), "\n", sep="")
      print(summary(colSums(matrixData[[whichChr]][[alleliclib]], na.rm=TRUE)))
      print(summary(colSums(matrixData[[whichChr]][[alleliclib]][-nrow(matrixData[[whichChr]][[alleliclib]]),], na.rm=TRUE)))
      
      # 20210311
      if (plotF121DecayCurvesByClusterJointly) {
        joint_alleliclib = paste(combined_label, allele, sep=".")
        if (alleliclib %in% c("F121_d0.ref", "F121_d0.alt")) {
          matrixData[[whichChr]][[joint_alleliclib]] = matrixData[[whichChr]][[alleliclib]]
        } else {
        matrixData[[whichChr]][[joint_alleliclib]] = cbind( matrixData[[whichChr]][[joint_alleliclib]],
                                                      matrixData[[whichChr]][[alleliclib]] )
        }
        # Just keep combined data
        matrixData[[whichChr]][[alleliclib]] = NULL
      }

    }
  }
  
  # matrixData[[whichChr]][["F121_R11_pooled.alt"]][1:5,1:5]
  # colSums(matrixData[[whichChr]][[alleliclib]], na.rm=TRUE)[1:5]
  
}  

names(matrixData[[whichChr]])
for(sample in names(matrixData[[whichChr]])) {
  cat(sample, ":\t", length(colnames(matrixData[[whichChr]][[sample]])), "\n", sep="")
}

nrow(matrixData[[whichChr]][[sample]])
# [1] 145 - Include trans and total contacts (rows 144 and 145).

# chrX:
# # mES_EBdiff_D0.ref:	2678
# # mES_EBdiff_D0.alt:	2678
# # mES_EBdiff_D7.ref:	1098
# # mES_EBdiff_D7.alt:	1098
# # mES_EBdiff_D20.ref:	1119
# # mES_EBdiff_D20.alt:	1119
# 
# Patski.ref:	4493
# Patski.alt:	4493
# EBdiff_F121_d0.ref:	1479
# EBdiff_F121_d0.alt:	1479
# EBdiff_F121_d7.ref:	89
# EBdiff_F121_d7.alt:	89
# EBdiff_F121_d20.ref:	301
# EBdiff_F121_d20.alt:	301
# EBdiff_F123_d0.ref:	1199
# EBdiff_F123_d0.alt:	1199
# EBdiff_F123_d7.ref:	62
# EBdiff_F123_d7.alt:	62
# EBdiff_F123_d20.ref:	414
# EBdiff_F123_d20.alt:	414
# EBdiff_TsixStop_d0.ref:	795
# EBdiff_TsixStop_d0.alt:	795
# EBdiff_TsixStop_d7.ref:	152
# EBdiff_TsixStop_d7.alt:	152
# EBdiff_TsixStop_d20.ref:	404
# EBdiff_TsixStop_d20.alt:	404
# TsixStop_d0.ref:	342
# TsixStop_d0.alt:	342
# TsixStop_d3.ref:	839
# TsixStop_d3.alt:	839
# TsixStop_d7.ref:	822
# TsixStop_d7.alt:	822
# TsixStop_d11.ref:	471
# TsixStop_d11.alt:	471
# TsixStop_NPC.ref:	2535
# TsixStop_NPC.alt:	2535
# F121_d0_rep1.ref:	147
# F121_d0_rep1.alt:	147
# F121_d0_rep2.ref:	2056
# F121_d0_rep2.alt:	2056
# F121_d0.ref:	2203
# F121_d0.alt:	2203
# F121_d3.ref:	894
# F121_d3.alt:	894
# F121_d7.ref:	732
# F121_d7.alt:	732
# F121_d11.ref:	608
# F121_d11.alt:	608
# F121_NPC.ref:	1054
# F121_NPC.alt:	1054

# 20200721
# F121_d0.ref:	2203
# F121_d0.alt:	2203
# F121_d3.ref:	894
# F121_d3.alt:	894
# F121_d7.ref:	732
# F121_d7.alt:	732
# F121_d11.ref:	608
# F121_d11.alt:	608
# F121_NPC.ref:	1054
# F121_NPC.alt:	1054

# 20200921
# TsixStop_d0.ref:	342
# TsixStop_d0.alt:	342
# TsixStop_d3.ref:	839
# TsixStop_d3.alt:	839
# TsixStop_d7.ref:	822
# TsixStop_d7.alt:	822
# TsixStop_d11.ref:	471
# TsixStop_d11.alt:	471
# TsixStop_NPC.ref:	2535
# TsixStop_NPC.alt:	2535

# 20200921
# F123_d0.ref:	1670
# F123_d0.alt:	1670
# F123_NPC.ref:	1030
# F123_NPC.alt:	1030

# 20200922
# F121_d0.ref:	2203
# F121_d0.alt:	2203
# F121_d3.ref:	894
# F121_d3.alt:	894
# F121_d7.ref:	732
# F121_d7.alt:	732
# F121_d11.ref:	608
# F121_d11.alt:	608
# F121_NPC.ref:	1054
# F121_NPC.alt:	1054
# F123_d0.ref:	1670
# F123_d0.alt:	1670
# F123_NPC.ref:	1030
# F123_NPC.alt:	1030

# 20210304
# TsixStop_NPC.ref:	2535
# TsixStop_NPC.alt:	2535
# F121_NPC.ref:	1054
# F121_NPC.alt:	1054

# 20210305
# TsixStop_d0.ref:	342
# TsixStop_d0.alt:	342
# F121_d0.ref:	2203
# F121_d0.alt:	2203

# 20210311
# combined.ref:	5491
# combined.alt:	5491

# 20210323
# EBdiff_F121_d0.ref:	1479
# EBdiff_F121_d0.alt:	1479
# EBdiff_F121_d7.ref:	89
# EBdiff_F121_d7.alt:	89
# EBdiff_F121_d20.ref:	301
# EBdiff_F121_d20.alt:	301
# EBdiff_TsixStop_d0.ref:	795
# EBdiff_TsixStop_d0.alt:	795
# EBdiff_TsixStop_d7.ref:	152
# EBdiff_TsixStop_d7.alt:	152
# EBdiff_TsixStop_d20.ref:	404
# EBdiff_TsixStop_d20.alt:	404

# # ------------------------------------
# # 20200805
# # Aside: test for Newcomb-Benford law:
# fc_list = list()
# for (fc in 1:9) {
#   fc_list[[as.character(fc)]] = 0 
# }
# 
# for(sample in names(matrixData[[whichChr]])) {
#   fc = substring( length(colnames(matrixData[[whichChr]][[sample]])), 1, 1)
#   fc_list[[fc]] = fc_list[[fc]] + 1
# }
# 
# barplot(unlist(fc_list))
# barplot((unlist(fc_list) / sum(unlist(fc_list)) * 100))
# 
# # ------------------------------------










##################################################
##################################################
##################################################
# Generate alleleic CDP data
# Also plot coverage plots and distance distributions

##################################################
# Select cells w/ sufficient coverage across chromosomes of interest.
# 1000 contacts / 40 alleles = 25 contacts per homologue.

allelicCDPcompList = list()
# allelicCDPcompList = allelicCDPcompList.all # If running additional nThreshes.
dataMatrices = list()

# 20210311
if (plotF121DecayCurvesByClusterJointly) {
  libs = c(combined_label)
}

for( nThresh in nThreshes ) {
  # nThresh = 50
  # nThresh=10000 # Vijay genome-wide 
  # nThresh=0 # Ritambhara
  # nThresh=100 # 20200803 - Patski min 100 contacts LMD and Spearman analysis 

  ##################################################
  # First select cells w/ sufficient coverage:
  # Filter for cells in both libraries that pass filter and are common to both alleles.
  
  dataMatrix.OG = list()
  matrixCsums.OG = list()
  plotShapes.OG = list()
  passFilterIdx = list()
  passFilterIdx.bothAlleles = list()
  dataMatrix.PF = list()
  matrixCsums.PF = list()
  plotShapes.PF = list()
  
  cellIDs.combo.biallelic = list()
  cell.libOfOrigin.biallelic = list()
  cell.combo.cols = list()
  
  # Initialize boolean vectors for inclusion of cells that pass filter across chromosomes.
  passFilterIdx.bothAlleles.allChroms = list()
  for(lib in libs) {
    allele="ref"
    alleliclib=paste(lib, allele, sep=".")
    passFilterIdx.bothAlleles.allChroms[[lib]] = rep(TRUE, ncol(data.matrix(matrixData[[1]][[alleliclib]])))
  }
  
  # Populate boolean vectors for inclusion of cells that pass filter across chromosomes.
  # As well as data matrices to subselect from.
  #
  # 20200721
  # Only use thresholds on chrX and chr1 to select cells used for all chromosomes
  # Otherwise too few cells.
  #
  for (whichChr in chroms2analyze) {
    
    cat("\n----------\n", whichChr, "\t", nThresh, "\n", sep="")
    
    dataMatrix.OG[[whichChr]] = list()
    matrixCsums.OG[[whichChr]] = list()
    plotShapes.OG[[whichChr]] = list()
    passFilterIdx[[whichChr]] = list()
    passFilterIdx.bothAlleles[[whichChr]] = list()
    dataMatrix.PF[[whichChr]] = list()
    matrixCsums.PF[[whichChr]] = list()
    plotShapes.PF[[whichChr]] = list()
    
    cellIDs.combo.biallelic[[whichChr]] = list()
    cell.libOfOrigin.biallelic[[whichChr]] = list()
    cell.combo.cols[[whichChr]] = list()
    
    for(lib in libs) {
      # for(lib in libs) {
      for(allele in c("ref", "alt")) {
        alleliclib=paste(lib, allele, sep=".")
        dataMatrix.OG[[whichChr]][[alleliclib]] = data.matrix(matrixData[[whichChr]][[alleliclib]][-c(144,145),])
        # 20200212 Note: Removing trans and total contacts ----------------------------------------->^^^^^^^^^^
        matrixCsums.OG[[whichChr]][[alleliclib]] = colSums(dataMatrix.OG[[whichChr]][[alleliclib]], na.rm=TRUE)
        names(matrixCsums.OG[[whichChr]][[alleliclib]]) = colnames(dataMatrix.OG[[whichChr]][[alleliclib]])
        passFilterIdx[[whichChr]][[alleliclib]] = matrixCsums.OG[[whichChr]][[alleliclib]] >= nThresh
        plotShapes.OG[[whichChr]][[alleliclib]] = rep(alleleShapes[allele], length(matrixCsums.OG[[whichChr]][[alleliclib]]))
      }
      passFilterIdx.bothAlleles[[whichChr]][[lib]] = passFilterIdx[[whichChr]][[paste(lib, "ref", sep=".")]] & passFilterIdx[[whichChr]][[paste(lib, "alt", sep=".")]]
      
      # 20200721
      # Only use thresholds on chrX and chr1 to select cells used for all chromosomes
      # Otherwise too few cells.
      if (F123sOnly || F121andF123) { 
        
        if ( whichChr %in% c("chr1") ) {  
          passFilterIdx.bothAlleles.allChroms[[lib]] = passFilterIdx.bothAlleles.allChroms[[lib]] & passFilterIdx.bothAlleles[[whichChr]][[lib]] 
          print(c(lib, sum( passFilterIdx.bothAlleles.allChroms[[lib]] )))
        }
        
      } else {
        
        if ( whichChr %in% c("chr1", "chrX") ) {  
          passFilterIdx.bothAlleles.allChroms[[lib]] = passFilterIdx.bothAlleles.allChroms[[lib]] & passFilterIdx.bothAlleles[[whichChr]][[lib]] 
          print(c(lib, sum( passFilterIdx.bothAlleles.allChroms[[lib]] )))
        }
        
      }
      
    }
  }
  
  # ----------
  #   chr1	50
  # [1] "Patski" "2800"  
  # [1] "EBdiff_F121_d0" "947"           
  # [1] "EBdiff_F121_d7" "77"            
  # [1] "EBdiff_F121_d20" "207"            
  # [1] "EBdiff_TsixStop_d0" "625"               
  # [1] "EBdiff_TsixStop_d7" "125"               
  # [1] "EBdiff_TsixStop_d20" "237"                
  # [1] "TsixStop_d0" "306"        
  # [1] "TsixStop_d3" "440"        
  # [1] "TsixStop_d7" "574"        
  # [1] "TsixStop_d11" "282"         
  # [1] "TsixStop_NPC" "2079"        
  # [1] "F121_d0_rep1" "109"         
  # [1] "F121_d0_rep2" "1823"        
  # [1] "F121_d0" "1932"   
  # [1] "F121_d3" "735"    
  # [1] "F121_d7" "594"    
  # [1] "F121_d11" "496"     
  # [1] "F121_NPC" "841"     
  # 
  # ----------
  #   chrX	50
  # [1] "Patski" "2055"  
  # [1] "EBdiff_F121_d0" "393"           
  # [1] "EBdiff_F121_d7" "47"            
  # [1] "EBdiff_F121_d20" "89"             
  # [1] "EBdiff_TsixStop_d0" "291"               
  # [1] "EBdiff_TsixStop_d7" "51"                
  # [1] "EBdiff_TsixStop_d20" "80"                 
  # [1] "TsixStop_d0" "89"         
  # [1] "TsixStop_d3" "86"         
  # [1] "TsixStop_d7" "159"        
  # [1] "TsixStop_d11" "83"          
  # [1] "TsixStop_NPC" "975"         
  # [1] "F121_d0_rep1" "65"          
  # [1] "F121_d0_rep2" "921"         
  # [1] "F121_d0" "986"    
  # [1] "F121_d3" "435"    
  # [1] "F121_d7" "243"    
  # [1] "F121_d11" "164"     
  # [1] "F121_NPC" "325"  
  
  # 20200803 Patski 100 contacts
  # ----------
  #   chr1	100
  # [1] "Patski" "2029"  
  # 
  # ----------
  #   chrX	100
  # [1] "Patski" "1161"  
  
  # 20200921 F123only
  # ----------
  #   chr1	50
  # [1] "F123_d0" "1494"   
  # [1] "F123_NPC" "876"    
  
  # 20200922 F121andF123
  # ----------
  #   chr1	50
  # [1] "F121_d0" "1932"   
  # [1] "F121_d3" "735"    
  # [1] "F121_d7" "594"    
  # [1] "F121_d11" "496"     
  # [1] "F121_NPC" "841"     
  # [1] "F123_d0" "1494"   
  # [1] "F123_NPC" "876"    
  
  # 20210304
  # ----------
  #   chr1	50
  # [1] "TsixStop_NPC" "2079"        
  # [1] "F121_NPC" "841"     
  # 
  # ----------
  #   chrX	50
  # [1] "TsixStop_NPC" "975"         
  # [1] "F121_NPC" "325"  
  
  # 20210305
  # ----------
  #   chr1	50
  # [1] "TsixStop_d0" "306"        
  # [1] "F121_d0" "1932"   
  # 
  # ----------
  #   chrX	50
  # [1] "TsixStop_d0" "89"         
  # [1] "F121_d0" "986"  
  
  # 20210311
  # ----------
  #   chr1	50
  # [1] "combined" "4598"    
  # 
  # ----------
  #   chrX	50
  # [1] "combined" "2153"  

  
  # 20210311
  if (!plotF121DecayCurvesByClusterJointly) {
    # For combining data
    whichLibs = c( libs, "combined" ) 
  } else {
    whichLibs = c( libs ) 
  }
  
  whichAllelicLibs = c(apply(expand.grid(whichLibs, c("ref", "alt")), 1, paste, collapse="."))
  
  # Use vectors for inclusion of cells that pass filter across chromosomes:
  # Only include cells where both alleles pass filter across chromosomes.
  for (whichChr in chroms2analyze) {
    for(lib in libs) {
      print(c(whichChr, lib, sum(passFilterIdx.bothAlleles.allChroms[[lib]])))
      for(allele in c("ref", "alt")) {
        alleliclib=paste(lib, allele, sep=".")
        dataMatrix.PF[[whichChr]][[alleliclib]] = dataMatrix.OG[[whichChr]][[alleliclib]][ ,which(passFilterIdx.bothAlleles.allChroms[[lib]]) ]
        matrixCsums.PF[[whichChr]][[alleliclib]] = matrixCsums.OG[[whichChr]][[alleliclib]][  passFilterIdx.bothAlleles.allChroms[[lib]] ]
        plotShapes.PF[[whichChr]][[alleliclib]] = plotShapes.OG[[whichChr]][[alleliclib]][ passFilterIdx.bothAlleles.allChroms[[lib]] ]
        
        # 20200229
        # all(colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) == names(matrixCsums.PF[[whichChr]][[alleliclib]]))
        cellIDs.combo.biallelic[[whichChr]][[alleliclib]] = paste(names(matrixCsums.PF[[whichChr]][[alleliclib]]), "_" , alleliclib)
        cell.libOfOrigin.biallelic[[whichChr]][[alleliclib]] =rep(alleliclib, length(matrixCsums.PF[[whichChr]][[alleliclib]]))
      }
      # 20200229
      cell.combo.cols[[whichChr]][[alleliclib]] = rep(alleleCols[[allele]], length(cell.libOfOrigin.biallelic[[whichChr]][[alleliclib]]))
    }
    
    # 20210311
    if (!plotF121DecayCurvesByClusterJointly) {
      
      # 20190529
      # Combine cells from all libraries for each allele
      # Only including cells where both alleles pass filter.
      for(allele in c("ref", "alt")) {
        combolabel = paste("combined", allele, sep=".")
        matrixCsums.PF[[whichChr]][[combolabel]] = c()
        dataMatrix.PF[[whichChr]][[combolabel]] = c()
        plotShapes.PF[[whichChr]][[combolabel]] = c()
        cellIDs.combo.biallelic[[whichChr]][[combolabel]] = c()
        cell.libOfOrigin.biallelic[[whichChr]][[combolabel]] = c()
        for(alleliclib in apply(expand.grid(libs, allele), 1, paste, collapse=".")) {
          names(matrixCsums.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
          matrixCsums.PF[[whichChr]][[combolabel]] = c( matrixCsums.PF[[whichChr]][[combolabel]], matrixCsums.PF[[whichChr]][[alleliclib]] )
          colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, colnames(dataMatrix.PF[[whichChr]][[alleliclib]]), sep=".")
          dataMatrix.PF[[whichChr]][[combolabel]] = cbind( dataMatrix.PF[[whichChr]][[combolabel]], dataMatrix.PF[[whichChr]][[alleliclib]] )
          names(plotShapes.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
          plotShapes.PF[[whichChr]][[combolabel]] = c( plotShapes.PF[[whichChr]][[combolabel]], plotShapes.PF[[whichChr]][[alleliclib]] )
          
          # all(colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) == names(matrixCsums.PF[[whichChr]][[alleliclib]]))
          cellIDs.combo.biallelic[[whichChr]][[combolabel]] = c(cellIDs.combo.biallelic[[whichChr]][[combolabel]], paste(names(matrixCsums.PF[[whichChr]][[alleliclib]]), "_" , alleliclib))
          cell.libOfOrigin.biallelic[[whichChr]][[combolabel]] = c(cell.libOfOrigin.biallelic[[whichChr]][[combolabel]], rep(alleliclib, length(matrixCsums.PF[[whichChr]][[alleliclib]])))
        }
        cell.combo.cols[[whichChr]][[combolabel]] = rep(alleleCols[[allele]], length(cell.libOfOrigin.biallelic[[whichChr]][[combolabel]]))
      }

    } # !plotF121DecayCurvesByClusterJointly  
      
    matrixCsums.biallelic = list()
    # combinedAlleleLabels = list()
    for(lib in whichLibs ) {
      matrixCsums.biallelic[[lib]] = rbind(
        matrixCsums.PF[[whichChr]][[paste(lib, "ref", sep=".")]], # [passFilterIdx.bothAlleles[[whichChr]][[lib]]],
        matrixCsums.PF[[whichChr]][[paste(lib, "alt", sep=".")]]) # [passFilterIdx.bothAlleles[[whichChr]][[lib]]])
      rownames(matrixCsums.biallelic[[lib]]) = c("ref", "alt")
    }
      
    # 20210311
    if (!plotF121DecayCurvesByClusterJointly) {
      
      # Combine cells from all libraries and alleles
      # Only including cells where both alleles pass filter.
      matrixCsums.PF[[whichChr]][["combined"]] = c()
      dataMatrix.PF[[whichChr]][["combined"]] = c()
      plotShapes.PF[[whichChr]][["combined"]] = c()
      cellIDs.combo = c()
      cell.libOfOrigin = c()
      for(alleliclib in apply(expand.grid(libs, c("ref", "alt")), 1, paste, collapse=".")) {
        # Renaming done above
        # names(matrixCsums.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
        matrixCsums.PF[[whichChr]][["combined"]] = c( matrixCsums.PF[[whichChr]][["combined"]], matrixCsums.PF[[whichChr]][[alleliclib]] )
        # colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, colnames(dataMatrix.PF[[whichChr]][[alleliclib]]), sep=".")
        dataMatrix.PF[[whichChr]][["combined"]] = cbind( dataMatrix.PF[[whichChr]][["combined"]], dataMatrix.PF[[whichChr]][[alleliclib]] )
        # names(plotShapes.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
        plotShapes.PF[[whichChr]][["combined"]] = c( plotShapes.PF[[whichChr]][["combined"]], plotShapes.PF[[whichChr]][[alleliclib]] )
        
        # all(colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) == names(matrixCsums.PF[[whichChr]][[alleliclib]]))
        cellIDs.combo = c(cellIDs.combo, paste(names(matrixCsums.PF[[whichChr]][[alleliclib]]), "_" , alleliclib))
        cell.libOfOrigin = c(cell.libOfOrigin, rep(alleliclib, length(matrixCsums.PF[[whichChr]][[alleliclib]])))
      }
      cell.combo.cols[[whichChr]][["combined"]] = rep(alleleCols[["ref"]], length(cell.libOfOrigin))
      cell.combo.cols[[whichChr]][["combined"]][grepl(".alt$", cell.libOfOrigin)] = alleleCols[["alt"]]
      
      # Don't rename here, b/c complicates date assignment
      # names(matrixCsums.PF[[whichChr]][["combined"]]) = cellIDs.combo
      # colnames(dataMatrix.PF[[whichChr]][["combined"]]) = cellIDs.combo
      
    } # !plotF121DecayCurvesByClusterJointly
     
    # 20200229
    # Combine cells from each library and alleles
    # Only including cells where both alleles pass filter.
    for ( lib in libs ) {
      matrixCsums.PF[[whichChr]][[lib]] = c()
      dataMatrix.PF[[whichChr]][[lib]] = c()
      plotShapes.PF[[whichChr]][[lib]] = c()
      cellIDs.combo = c()
      cell.libOfOrigin = c()
      for(alleliclib in apply(expand.grid(alleliclib, c("ref", "alt")), 1, paste, collapse=".")) {
        # Renaming done above
        # names(matrixCsums.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
        matrixCsums.PF[[whichChr]][[lib]] = c( matrixCsums.PF[[whichChr]][[lib]], matrixCsums.PF[[whichChr]][[alleliclib]] )
        # colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, colnames(dataMatrix.PF[[whichChr]][[alleliclib]]), sep=".")
        dataMatrix.PF[[whichChr]][[lib]] = cbind( dataMatrix.PF[[whichChr]][[lib]], dataMatrix.PF[[whichChr]][[alleliclib]] )
        # names(plotShapes.PF[[whichChr]][[alleliclib]]) = paste(alleliclib, names(matrixCsums.PF[[whichChr]][[alleliclib]]), sep=".")
        plotShapes.PF[[whichChr]][[lib]] = c( plotShapes.PF[[whichChr]][[lib]], plotShapes.PF[[whichChr]][[alleliclib]] )
        
        # all(colnames(dataMatrix.PF[[whichChr]][[alleliclib]]) == names(matrixCsums.PF[[whichChr]][[alleliclib]]))
        cellIDs.combo = c(cellIDs.combo, paste(names(matrixCsums.PF[[whichChr]][[alleliclib]]), "_" , alleliclib))
        cell.libOfOrigin = c(cell.libOfOrigin, rep(alleliclib, length(matrixCsums.PF[[whichChr]][[alleliclib]])))
      }
      cell.combo.cols[[whichChr]][[lib]] = rep(alleleCols[["ref"]], length(cell.libOfOrigin))
      cell.combo.cols[[whichChr]][[lib]][grepl(".alt$", cell.libOfOrigin)] = alleleCols[["alt"]]
    }
      
  } # whichChr
  
  tpCols = c(libCols, "black")
  names(tpCols) = c( names(libCols), "combined")
  
  # 2020921
  # Don't analyze combined data - no point
  whichLibs = c( libs )
  whichAllelicLibs = c(apply(expand.grid(whichLibs, c("ref", "alt")), 1, paste, collapse="."))
  

  
  
  
  
  
  
  
  
  ##################################################
  # Perform actual analysis
  
  nThreshLabel = as.character(nThresh)
  
  allelicCDPcompList[[nThreshLabel]] = list()
  dataMatrices[[nThreshLabel]] = list()
  
  for (whichChr in chroms2analyze) {
    # whichChr = "chr1"
    
    cat("\n----------\n", whichChr, "\t", nThresh, "\n", sep="")
    
    allelicCDPcompList[[nThreshLabel]][[whichChr]] = list()
    dataMatrices[[nThreshLabel]][[whichChr]] = list()
    
    ##################################################
    # Ouput folder
    subsupermainOutputDir = file.path(supermainOutputDir, 
                                      paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
    if ( ! file.exists( subsupermainOutputDir ) ) {
      if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
    }
    
    ##################################################
    # Contact decay rate from 50kbTo8Mb
    
    # # Proportion of cis reads in 1-20kb to 60kb-1Mb
    # 
    # # propCis20kb = list()
    # # for(lib in libs) {
    # #   propCis20kb[[lib]] = colSums(matrixCsums.PF[[whichChr]][[lib]][binDist > 1024 & binDist < 20000,], na.rm=TRUE) / 
    # #     colSums(matrixCsums.PF[[whichChr]][[lib]][binDist > 60000 & binDist < 1000000,], na.rm=TRUE)
    # # }
    # 
    # propDecayRate50kbTo8Mb = list()
    # propDecayRate50kbTo8Mb.allelic = list()
    # for(lib in whichAllelicLibs) {
    #   propDecayRate50kbTo8Mb[[lib]] = (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)),]) 
    #                                    - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)),])) /
    #     (log10(8000000)-log10(500000))  
    #   
    #   # Allelic
    #   propDecayRate50kbTo8Mb.allelic[[lib]] = rbind(
    #     (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)), 1:(ncol(dataMatrix.PF[[whichChr]][[lib]])/2)]) 
    #      - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)), 1:(ncol(dataMatrix.PF[[whichChr]][[lib]])/2)])) /
    #     (log10(8000000)-log10(500000)),
    #     (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)), (ncol(dataMatrix.PF[[whichChr]][[lib]])/2+1):ncol(dataMatrix.PF[[whichChr]][[lib]])]) 
    #      - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)), (ncol(dataMatrix.PF[[whichChr]][[lib]])/2+1):ncol(dataMatrix.PF[[whichChr]][[lib]])])) /
    #       (log10(8000000)-log10(500000)) 
    #   )
    # 
    # }
    
    decayRate50kbTo8Mb = list()
    decayRate50kbTo8Mb.allelic = list()
    for(lib in whichAllelicLibs) {
      decayRate50kbTo8Mb[[lib]] = (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)),]) 
                                   - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)),])) /
        (log10(8000000)-log10(50000))  
      
      # # Biallelic
      # decayRate50kbTo8Mb.allelic[[lib]] = rbind(
      #   (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)), 1:(ncol(dataMatrix.PF[[whichChr]][[lib]])/2)]) 
      #    - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)), 1:(ncol(dataMatrix.PF[[whichChr]][[lib]])/2)])) /
      #     (log10(8000000)-log10(50000)),
      #   (log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-50000)), (ncol(dataMatrix.PF[[whichChr]][[lib]])/2+1):ncol(dataMatrix.PF[[whichChr]][[lib]])]) 
      #    - log10(dataMatrix.PF[[whichChr]][[lib]][which.min(abs(binDist-8000000)), (ncol(dataMatrix.PF[[whichChr]][[lib]])/2+1):ncol(dataMatrix.PF[[whichChr]][[lib]])])) /
      #     (log10(8000000)-log10(50000)) 
      # )
    }
    
    
    
    
    
    
    
    
    
    
    ##################################################
    # Coverage level ordering and coloring
    # Colors ordered by cell coverage
    
    covCols = list()
    covOrd = list()
    
    covCols.biallelic = list()
    covOrd.biallelic  = list()
    
    for(lib in whichAllelicLibs) {
      covOrd[[lib]] = order(matrixCsums.PF[[whichChr]][[lib]])
      names(covOrd[[lib]]) = names(matrixCsums.PF[[whichChr]][[lib]])[ covOrd[[lib]] ]
      # covCols[[lib]] = viridis(length(matrixCsums.PF[[whichChr]][[lib]]), alpha=0.5)
      covCols[[lib]][covOrd[[lib]]] = viridis(length(matrixCsums.PF[[whichChr]][[lib]]), alpha=0.5)
      names(covCols[[lib]]) = names(covOrd[[lib]]) 
      
      # biallelic
      # 20190529
      if (lib!="all.combined") {
        covOrd.biallelic[[lib]] = order(colSums(matrixCsums.biallelic[[sub("\\..*", "", lib)]], na.rm=TRUE))
        names(covOrd.biallelic[[lib]]) = colnames(matrixCsums.biallelic[[sub("\\..*", "", lib)]])[ covOrd.biallelic[[lib]] ]
        # covCols.biallelic[[lib]] = viridis(ncol(matrixCsums.biallelic[[sub("\\..*", "", lib)]]), alpha=0.5)
        covCols.biallelic[[lib]][covOrd.biallelic[[lib]]] = viridis(ncol(matrixCsums.biallelic[[sub("\\..*", "", lib)]]), alpha=0.5)
        names(covCols.biallelic[[lib]]) = names(covOrd.biallelic[[lib]]) 
        
      }
    }
    
    # # check
    # head(names(covOrd[[lib]]))
    # head(matrixCsums.PF[[whichChr]][[lib]][ covOrd[[lib]] ])
    # 
    # barplot(rep(1,length( covCols[[lib]] )), col=covCols[[lib]] )
    # barplot(rep(1,length( covCols[[lib]] )), col=covCols[[lib]][covOrd[[lib]]] )
    # 
    # barplot(matrixCsums.PF[[whichChr]][[lib]], col= covCols[[lib]], las=2, cex.names=0.5 )
    # barplot(matrixCsums.PF[[whichChr]][[lib]][ covOrd[[lib]] ], col= covCols[[lib]][covOrd[[lib]]], las=2, cex.names=0.5 )
    
    
    
    
    
    
    
    
    
    
    ##################################################
    # Library labeling
    
    # matrixColIBCs = list()
    matrixColTPs = list()
    matrixColTPcols = list()
    timepointOrd = list()
    timepointThenCovOrd = list()
    # timepointThenDecayOrd = list()
    
    # NOT NECESSARY
    # matrixColIBCs.biallelic = list()
    # matrixColTPs.biallelic = list()
    # matrixColTPcols.biallelic = list()
    # timepointOrd.biallelic = list()
    # timepointThenCovOrd.biallelic = list()
    # timepointThenDecayOrd.biallelic = list()
    
    for(lib in whichAllelicLibs) {
      matrixColTPs[[lib]] = sub("\\..*$", "", cell.libOfOrigin.biallelic[[whichChr]][[lib]])
      matrixColTPcols[[lib]] = rep(tpCols[sub("\\..*", "", lib)], length(matrixColTPs[[lib]]))
      # matrixColTPcols[[lib]] = rep("black", length(matrixColTPs[[lib]]))
      # matrixColTPcols[[lib]][matrixColTPs[[lib]]=="ML2"] = tpCols["ML2"]
      # matrixColTPcols[[lib]][matrixColTPs[[lib]]=="PL1"] = tpCols["PL1"]
      # matrixColTPcols[[lib]][matrixColTPs[[lib]]=="PL2"] = tpCols["PL2"]
      names(matrixColTPcols[[lib]]) = matrixColTPs[[lib]]
      
      timepointOrd[[lib]] = order(cell.libOfOrigin.biallelic[[whichChr]][[lib]])
      
      # Order by timepoint then by coverage
      timepointThenCovOrd[[lib]] = order(cell.libOfOrigin.biallelic[[whichChr]][[lib]])
      names(timepointThenCovOrd[[lib]]) = matrixColTPs[[lib]][ covOrd[[lib]] ][timepointThenCovOrd[[lib]]]
      
      # # Order by timepoint then by decay
      # timepointThenDecayOrd[[lib]] = order(cell.libOfOrigin)
      # names(timepointThenDecayOrd[[lib]]) = matrixColTPs[[lib]][ order(decayRate50kbTo8Mb[[lib]]) ][timepointThenDecayOrd[[lib]]]
    }   
    
    
    
    
    
    
    
    
    
    
    ##################################################
    # Coverage plots
    
    covplotsMaindir = file.path(subsupermainOutputDir, "coveragePlots")
    if ( ! file.exists( covplotsMaindir ) ) {
      if ( ! dir.create( covplotsMaindir, showWarnings=FALSE, recursive=TRUE )  ) covplotsMaindir = subsupermainOutputDir
    }
    
    for(lib in whichLibs) {
      
      orderLib = paste(lib, "ref", sep=".")
      
      # Ranked by timepoint, then coverage, colored by timepoint
      PDFfilename = file.path(covplotsMaindir, paste(projectName, ".coverage.barplots.", lib, ".pdf", sep=""))
      pdf(PDFfilename, height=12, width=12)
      par(mfrow=c(3,1))
      
      for(allele in c("ref", "alt")) {
        
        allelicLib = paste(lib, allele, sep=".")
        barplot(matrixCsums.PF[[whichChr]][[allelicLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], 
                col=matrixColTPcols[[orderLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], border=NA, names.arg=NA, 
                xlab=paste(length(matrixCsums.PF[[whichChr]][[allelicLib]]), " cells >= ", nThresh, " for both alleles", sep=""), 
                ylab="total cis coverage", 
                main=paste(allele, "coverage ranked by ", orderLib, " coverage", sep=""))
      }
      
      barplot(matrixCsums.PF[[whichChr]][[paste(lib, "ref", sep=".")]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]] - matrixCsums.PF[[whichChr]][[paste(lib, "alt", sep=".")]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], 
              col=matrixColTPcols[[orderLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], border=NA, names.arg=NA, 
              xlab=paste(length(matrixCsums.PF[[whichChr]][[orderLib]]), " cells >= ", nThresh, " for both alleles", sep=""), 
              ylab="total cis coverage", 
              main=paste("Differential coverage ranked by ", orderLib, "coverage", sep=""),
              sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "))
      
      dev.off()
      
    }
    
    for(lib in whichLibs) {
      
      orderLib = paste(lib, "ref", sep=".")
      
      # ylim at mean + 3sd
      allLog10counts =  log10( c( matrixCsums.PF[[whichChr]][[paste(lib, "ref", sep=".")]] + 1, matrixCsums.PF[[whichChr]][[paste(lib, "alt", sep=".")]] + 1 ) )
      ylims = c(0, 10 ^ (mean(allLog10counts) + 3 *sd(allLog10counts)) ) 
      
      # Ranked by timepoint, then coverage, colored by timepoint
      PDFfilename = file.path(covplotsMaindir, paste(projectName, ".coverage.barplots.", lib, ".3sd.pdf", sep=""))
      pdf(PDFfilename, height=12, width=12)
      par(mfrow=c(3,1))
      
      for(allele in c("ref", "alt")) {
        
        allelicLib = paste(lib, allele, sep=".")
        barplot(matrixCsums.PF[[whichChr]][[allelicLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], 
                col=matrixColTPcols[[orderLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], border=NA, names.arg=NA, 
                xlab=paste(length(matrixCsums.PF[[whichChr]][[allelicLib]]), " cells >= ", nThresh, " for both alleles", sep=""), 
                ylab="total cis coverage", ylim=ylims,
                main=paste(allele, "coverage ranked by ", orderLib, " coverage", sep=""))
      }
      
      # ylim at mean + 3sd
      allLog10counts =  log10( abs(c(matrixCsums.PF[[whichChr]][[paste(lib, "ref", sep=".")]] - matrixCsums.PF[[whichChr]][[paste(lib, "alt", sep=".")]])) + 1 )
      ylims = c( -10 ^ (mean(allLog10counts) + 3 *sd(allLog10counts)), 10 ^ (mean(allLog10counts) + 3 *sd(allLog10counts)) ) 
      
      barplot(c(matrixCsums.PF[[whichChr]][[paste(lib, "ref", sep=".")]] - matrixCsums.PF[[whichChr]][[paste(lib, "alt", sep=".")]])[covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], 
              col=matrixColTPcols[[orderLib]][covOrd[[orderLib]]][timepointThenCovOrd[[orderLib]]], border=NA, names.arg=NA, 
              xlab=paste(length(matrixCsums.PF[[whichChr]][[orderLib]]), " cells >= ", nThresh, " for both alleles", sep=""), 
              ylab="total cis coverage", ylim=ylims,
              main=paste("Differential coverage ranked by ", orderLib, "coverage", sep=""),
              sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "))
      
      dev.off()
      
    }
    
    
    # # Ranked by timepoint, then coverage, colored by timepoint -- stacked by allele
    # PDFfilename = file.path(covplotsMaindir, paste(projectName, ".coverage.barplots.rankAndColByTimepoint.", nameText, ".allelic.pdf", sep=""))
    # pdf(PDFfilename, height=6, width=12)
    # barplot(matrixCsums.biallelic[[lib]][ ,covOrd.biallelic[[lib]] ][ ,timepointThenCovOrd.biallelic[[lib]] ], 
    #         #col=matrixColTPcols.biallelic[[lib]][covOrd.biallelic[[lib]]][timepointThenCovOrd.biallelic[[lib]]], border=NA, names.arg=NA, 
    #         col=c("red","blue"), border=NA, axisnames = FALSE, 
    #         xlab=paste(ncol(matrixCsums.biallelic[[lib]]), " cell >= ", nThresh, " for both alleles", sep=""), 
    #         ylab="total cis coverage", 
    #         main=paste("ranked by timepoint then coverage: ", titleText, sep=""))
    # legend("topleft", legend=rownames(matrixCsums.biallelic[[lib]]), fill=c("red","blue"), bty="n")
    # # barplot(matrixCsums.biallelic[[lib]][,covOrd.biallelic[[lib]]][timepointThenCovOrd.biallelic[[lib]]], 
    # #         col=matrixColTPcols.biallelic[[lib]][covOrd.biallelic[[lib]]][timepointThenCovOrd.biallelic[[lib]]], border=NA, names.arg=NA, add=TRUE)
    # dev.off()
    
    
    
    
    
    
    
    
    
    
    ##################################################
    # Preprocess data
    
    # for (bins2Use in c("50kb-8Mb", "1kb-Cis", "1kb-Cis.merge10")) {
    for (bins2Use in c("1kb-100Mb.merge10")) {
      # bins2Use="1kb-100Mb.merge10"
      # bins2Use="1kb-Cis.merge10"
      # bins2Use"1kb-8Mb.merge10"
      # bins2Use="50kb-8Mb" # Vijay
      # bins2Use="1kb-80Mb"
      
      # Also work:
      # bins2Use="20kb-8Mb"
      # bins2Use="20kb-20Mb"
      # bins2Use="20kb-Cis"
      # bins2Use="20kb-5Mb"
      #
      # Dont' work
      # bins2Use="50kb-8Mb+Trans"
      # bins2Use="0-10Mb"
      # bins2Use="allCis"
      # bins2Use="all"
      
      cat("\n", whichChr, "\t", nThresh, "\t", bins2Use, "\n", sep="")
      
      allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]] = list()
      dataMatrices[[nThreshLabel]][[whichChr]][[bins2Use]] = list()
      
      mainOutputDir = file.path(subsupermainOutputDir, 
                                bins2Use)
      if ( ! file.exists( mainOutputDir ) ) {
        if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
      }
      
      matrixRvar = list()
      dataMatrix = list()
      dataMatrix.scaled = list()
      dataMatrix.prop = list()
      
      for(lib in whichAllelicLibs) {
        
        if (bins2Use=="50kb-8Mb") {
          # Only bins from 50kb to 8Mb
          binBoolIdx = binDist >= 50000 & binDist <= 8000000
        } else if (bins2Use %in% c("1kb-8Mb", "1kb-8Mb.merge10") ) {
          # Only bins from 1kb to 8Mb
          binBoolIdx = binDist > 1024 & binDist <= 8000000
        } else if (bins2Use=="50kb-8Mb+Trans") {
          # Only bins from 50kb to 8Mb
          binBoolIdx = ( binDist >= 50000 & binDist <= 8000000 ) | c(rep(FALSE, length(binBoolIdx)-1), TRUE)
        } else if (bins2Use=="1kb-50Mb") {
          binBoolIdx = binDist > 1024 & binDist <= 50000000
        } else if (bins2Use %in% c("1kb-100Mb", "1kb-100Mb.merge10") ) {
          binBoolIdx = binDist > 1024 & binDist <= 100000000
        } else if (bins2Use %in% c("1kb-Cis", "1kb-Cis.merge10") ) {
          # Only bins from 20kb to 20Mb
          binBoolIdx = ( binDist > 1024) & c(rep(TRUE, length(binBoolIdx)-1), FALSE)
        } else if (bins2Use=="20kb-Cis") {
          # Only bins from 20kb to 20Mb
          binBoolIdx = ( binDist >= 20000 ) & c(rep(TRUE, length(binBoolIdx)-1), FALSE)
        } else if (bins2Use=="allCis") {
          binBoolIdx = c(rep(TRUE, length(binDist)-1), FALSE)
        }
        print(sum(binBoolIdx))

        dataMatrix[[lib]] = dataMatrix.PF[[whichChr]][[lib]][binBoolIdx, ]
        # dataMatrix[[lib]][is.na(dataMatrix[[lib]])] = 0
        
        if (grepl( "merge10", bins2Use )) {
          altSeqIdx = seq(1, sum(binBoolIdx), 10)
          # altSeqIdx
          dataMatrix[[lib]] = dataMatrix.PF[[whichChr]][[lib]][binBoolIdx, ][altSeqIdx, ]
          dataMatrix[[lib]][is.na(dataMatrix[[lib]])] = 0
          for ( i in 1:10) {
            newValues = dataMatrix.PF[[whichChr]][[lib]][binBoolIdx, ][altSeqIdx[1:(length(altSeqIdx)-1)] + i, ]
            newValues[is.na(newValues)] = 0
            dataMatrix[[lib]][1:(length(altSeqIdx)-1), ] = dataMatrix[[lib]][1:(length(altSeqIdx)-1), ] + newValues
          }
          # Add extras 
          i = 1
          while(altSeqIdx[length(altSeqIdx)]+i <= sum(binBoolIdx)) {
            print(i)
            newValues = dataMatrix.PF[[whichChr]][[lib]][binBoolIdx, ][altSeqIdx[length(altSeqIdx)] + i, ]
            newValues[is.na(newValues)] = 0
            dataMatrix[[lib]][length(altSeqIdx), ] = dataMatrix[[lib]][length(altSeqIdx), ] + newValues
            i = i + 1
          }
        }
        
        print(any(is.na(dataMatrix[[lib]])))
        print(any(dataMatrix[[lib]]==0))
        
      } # lib
      
      # Don't do this here... messes up coverage and sampple indexing.
      # Besides, shouldn't have ot now.
      # 
      # # Remove haplocells with zero counts from both allleles
      # zeroIdx = c()
      # for(lib in whichAllelicLibs) {
      #   # Remove zero cells
      #   zeroIdx = c(zeroIdx, which(colSums(dataMatrix[[lib]], na.rm=TRUE)==0))
      # }
      # zeroIdx = unique(zeroIdx)
      # print(zeroIdx)
      # if(length(zeroIdx)>0) {
      #   for(lib in whichAllelicLibs) {
      #     dataMatrix[[lib]] = dataMatrix[[lib]][, -zeroIdx]
      #   }
      # }
      # dim(dataMatrix[[lib]])
      
      for(lib in whichAllelicLibs) {
        
        # *** CAN'T DO THIS WITHOUT SCREWING UP PLOT DIMENSION CONSISTENCY AND COMPATIBILITY ***  
        # # PCA doesn't work on zero variance rows
        # # Unlikely to have zero variance columns
        # matrixRvar[[lib]] = apply(dataMatrix[[lib]], 1, var, na.rm=TRUE)
        # zeroVarRowIdx = matrixRvar[[lib]] == 0 | is.na(matrixRvar[[lib]])
        # dataMatrix[[lib]] = dataMatrix[[lib]][!zeroVarRowIdx,]
        
        # columnwise scaling
        dataMatrix.scaled[[lib]] = scale(dataMatrix[[lib]]) 
        # dataMatrix.scaled[[lib]][ dataMatrix.scaled[[lib]]==0 ] = NaN
        print(any(is.na(dataMatrix.scaled[[lib]])))
        print(which(apply(dataMatrix.scaled[[lib]], 2, function(x){ any(is.na(x)) })))
        
        # columnwise proportions
        dataMatrix.prop[[lib]] = apply(dataMatrix[[lib]], MARGIN=2, FUN=getProportions)
        # to avoid log10(0) *** BUT BETTER TO JUST NaN the bin ***
        # dataMatrix.prop[[lib]] = apply(dataMatrix.prop[[lib]], MARGIN=2, FUN=addNZval)
        # dataMatrix.prop[[lib]][ dataMatrix.prop[[lib]]==0 ] = NaN
        dataMatrix.prop[[lib]][ dataMatrix.prop[[lib]]==0 ] = min(dataMatrix.prop[[lib]][ dataMatrix.prop[[lib]]!=0 ], na.rm=TRUE)/1000
        dataMatrix.prop[[lib]] = log10(dataMatrix.prop[[lib]])
        print(any(is.na(dataMatrix.prop[[lib]])))
        print(which(apply(dataMatrix.prop[[lib]], 2, function(x){ any(is.na(x)) })))
        
        # print(lib)
        # print(dataMatrix[[lib]][1:3, 1:2])
      }

      # head(dataMatrix[[lib]][,1:3] )
      # tail(dataMatrix[[lib]][,1:3] )
      
      binDistDesc = prettyBPs(binDist[binBoolIdx])
      if (length(binDistDesc) > 100) {
        binDistDesc[seq(2,length(binDistDesc),2)] = ""
      }
      binDistDesc
      
      prettyBPs(binDist)
      # [1] "1.02kb"   "1.12kb"   "1.22kb"   "1.33kb"   "1.45kb"   "1.58kb"   "1.72kb"   "1.88kb"   "2.05kb"   "2.23kb"  
      # [11] "2.44kb"   "2.66kb"   "2.9kb"    "3.16kb"   "3.44kb"   "3.76kb"   "4.1kb"    "4.47kb"   "4.87kb"   "5.31kb"  
      # [21] "5.79kb"   "6.32kb"   "6.89kb"   "7.51kb"   "8.19kb"   "8.93kb"   "9.74kb"   "10.62kb"  "11.59kb"  "12.63kb" 
      # [31] "13.78kb"  "15.02kb"  "16.38kb"  "17.87kb"  "19.48kb"  "21.25kb"  "23.17kb"  "25.27kb"  "27.55kb"  "30.05kb" 
      # [41] "32.77kb"  "35.73kb"  "38.97kb"  "42.49kb"  "46.34kb"  "50.54kb"  "55.11kb"  "60.1kb"   "65.54kb"  "71.47kb" 
      # [51] "77.94kb"  "84.99kb"  "92.68kb"  "101.07kb" "110.22kb" "120.19kb" "131.07kb" "142.94kb" "155.87kb" "169.98kb"
      # [61] "185.36kb" "202.14kb" "220.44kb" "240.39kb" "262.14kb" "285.87kb" "311.74kb" "339.96kb" "370.73kb" "404.28kb"
      # [71] "440.87kb" "480.77kb" "524.29kb" "571.74kb" "623.49kb" "679.92kb" "741.46kb" "808.56kb" "881.74kb" "961.55kb"
      # [81] "1Mb"      "1.1Mb"    "1.2Mb"    "1.4Mb"    "1.5Mb"    "1.6Mb"    "1.8Mb"    "1.9Mb"    "2.1Mb"    "2.3Mb"   
      # [91] "2.5Mb"    "2.7Mb"    "3Mb"      "3.2Mb"    "3.5Mb"    "3.8Mb"    "4.2Mb"    "4.6Mb"    "5Mb"      "5.4Mb"   
      # [101] "5.9Mb"    "6.5Mb"    "7.1Mb"    "7.7Mb"    "8.4Mb"    "9.1Mb"    "10Mb"     "10.9Mb"   "11.9Mb"   "12.9Mb"  
      # [111] "14.1Mb"   "15.4Mb"   "16.8Mb"   "18.3Mb"   "20Mb"     "21.8Mb"   "23.7Mb"   "25.9Mb"   "28.2Mb"   "30.8Mb"  
      # [121] "33.6Mb"   "36.6Mb"   "39.9Mb"   "43.5Mb"   "47.5Mb"   "51.7Mb"   "56.4Mb"   "61.5Mb"   "67.1Mb"   "73.2Mb"  
      # [131] "79.8Mb"   "87Mb"     "94.9Mb"   "103.5Mb"  "112.9Mb"  "123.1Mb"  "134.2Mb"  "146.4Mb"  "159.6Mb"  "174.1Mb" 
      # [141] "189.8Mb"  "207Mb"    "225.7Mb" 
      
      sum(binBoolIdx)
      # 132 - when bin to use is "1kb-100Mb.merge10"
      
      prettyBPs(binDist[binBoolIdx])
      # [1] "1.12kb"   "1.22kb"   "1.33kb"   "1.45kb"   "1.58kb"   "1.72kb"   "1.88kb"  
      # [8] "2.05kb"   "2.23kb"   "2.44kb"   "2.66kb"   "2.9kb"    "3.16kb"   "3.44kb"  
      # [15] "3.76kb"   "4.1kb"    "4.47kb"   "4.87kb"   "5.31kb"   "5.79kb"   "6.32kb"  
      # [22] "6.89kb"   "7.51kb"   "8.19kb"   "8.93kb"   "9.74kb"   "10.62kb"  "11.59kb" 
      # [29] "12.63kb"  "13.78kb"  "15.02kb"  "16.38kb"  "17.87kb"  "19.48kb"  "21.25kb" 
      # [36] "23.17kb"  "25.27kb"  "27.55kb"  "30.05kb"  "32.77kb"  "35.73kb"  "38.97kb" 
      # [43] "42.49kb"  "46.34kb"  "50.54kb"  "55.11kb"  "60.1kb"   "65.54kb"  "71.47kb" 
      # [50] "77.94kb"  "84.99kb"  "92.68kb"  "101.07kb" "110.22kb" "120.19kb" "131.07kb"
      # [57] "142.94kb" "155.87kb" "169.98kb" "185.36kb" "202.14kb" "220.44kb" "240.39kb"
      # [64] "262.14kb" "285.87kb" "311.74kb" "339.96kb" "370.73kb" "404.28kb" "440.87kb"
      # [71] "480.77kb" "524.29kb" "571.74kb" "623.49kb" "679.92kb" "741.46kb" "808.56kb"
      # [78] "881.74kb" "961.55kb" "1Mb"      "1.1Mb"    "1.2Mb"    "1.4Mb"    "1.5Mb"   
      # [85] "1.6Mb"    "1.8Mb"    "1.9Mb"    "2.1Mb"    "2.3Mb"    "2.5Mb"    "2.7Mb"   
      # [92] "3Mb"      "3.2Mb"    "3.5Mb"    "3.8Mb"    "4.2Mb"    "4.6Mb"    "5Mb"     
      # [99] "5.4Mb"    "5.9Mb"    "6.5Mb"    "7.1Mb"    "7.7Mb"    "8.4Mb"    "9.1Mb"   
      # [106] "10Mb"     "10.9Mb"   "11.9Mb"   "12.9Mb"   "14.1Mb"   "15.4Mb"   "16.8Mb"  
      # [113] "18.3Mb"   "20Mb"     "21.8Mb"   "23.7Mb"   "25.9Mb"   "28.2Mb"   "30.8Mb"  
      # [120] "33.6Mb"   "36.6Mb"   "39.9Mb"   "43.5Mb"   "47.5Mb"   "51.7Mb"   "56.4Mb"  
      # [127] "61.5Mb"   "67.1Mb"   "73.2Mb"   "79.8Mb"   "87Mb"     "94.9Mb"  
      
      # Merge bins
      if (grepl( "merge10", bins2Use)) {
        altSeqIdx = seq(1, sum(binBoolIdx), 10)
        binDistDesc = binDistDesc[altSeqIdx]
      }
      # [1] "1.12kb"   "2.66kb"   "6.32kb"   "15.02kb"  "35.73kb"  "84.99kb"  "202.14kb" "480.77kb" "1.1Mb"    "2.7Mb"   
      # [11] "6.5Mb"    "15.4Mb"   "36.6Mb"   "87Mb"   
      # This is the same as the first column above. Yay!
      
      dataMatrices[[nThreshLabel]][[whichChr]][[bins2Use]][["scaled"]] = dataMatrix.scaled
      dataMatrices[[nThreshLabel]][[whichChr]][[bins2Use]][["prop"]] = dataMatrix.prop
      
      
      
      
      
      
      
      
      
      
      ##################################################
      # Correlate ref and alt alleles - Spearman
      # https://stackoverflow.com/questions/6713973/how-do-i-calculate-correlation-between-corresponding-columns-of-two-matrices-and
      
      cat("\n\nSpearman correlation of allelic CDPs\n")
      
      allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]] = list()

      cor.spearman.pw_matrix <- function(mat1, mat2){
        resVector = c()
        for(i in 1:ncol(mat1)) {
          colDist = cor(mat1[,i], mat2[,i], method="spearman", use="pairwise.complete")
          resVector = c(resVector, colDist)
        }
        names(resVector) = colnames(mat1)
        return(resVector)
      }
      
      corMaindir = file.path(mainOutputDir, "distanceDistributions.refVsAlt", "corSpearman")
      if ( ! file.exists( corMaindir ) ) {
        if ( ! dir.create( corMaindir, showWarnings=FALSE, recursive=TRUE )  ) corMaindir = mainOutputDir
      }
      
      # 20190529
      for (lib in whichLibs) {
        allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]][[lib]] = list()
      }
      
      for (useScaled in c(TRUE, FALSE)) {
        
        allelicCDPcomp = list()
        
        # 20190529
        for (lib in whichLibs) {
          print(lib)
          
          if (useScaled) {
            # allelicCDPcomp[[lib]] =  diag(cor(dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 
            #                               dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 
            #                               method="spearman", use="pairwise.complete"))
            allelicCDPcomp[[lib]] = cor.spearman.pw_matrix(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
                                                           dataMatrix.scaled[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]][[lib]][["scaled"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]][[lib]][["scaled"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          } else {
            # allelicCDPcomp[[lib]] =  diag(cor(dataMatrix.prop[[paste(lib, "ref", sep=".")]], 
            #                               dataMatrix.prop[[paste(lib, "alt", sep=".")]], 
            #                               method="spearman", use="pairwise.complete"))
            allelicCDPcomp[[lib]] = cor.spearman.pw_matrix(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
                                                           dataMatrix.scaled[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]][[lib]][["prop"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Spearman"]][[lib]][["prop"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          }
          
          if (! all(is.na(allelicCDPcomp[[lib]]))) {
            
            dp <- density(allelicCDPcomp[[lib]], na.rm=TRUE)
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            maxDensity = min(15, max(c(dp$y), na.rm=TRUE))
            
            titleText = paste("Histogram of ref vs alt Spearman correlations\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltCorDensity.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid", xlim=c(-1,1),
                 main=titleText, ylim=c(0, maxDensity))
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
          }
          
        }
      }
      
      
      
      
      
      
      
      
      
      
      ##################################################
      # Correlate ref and alt alleles - Pearson
      # https://stackoverflow.com/questions/6713973/how-do-i-calculate-correlation-between-corresponding-columns-of-two-matrices-and
      
      cat("\n\nPearson correlation of allelic CDPs\n")
      
      allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]] = list()
      
      cor.pearson.pw_matrix <- function(mat1, mat2){
        resVector = c()
        for(i in 1:ncol(mat1)) {
          colDist = cor(mat1[,i], mat2[,i], method="pearson", use="pairwise.complete")
          resVector = c(resVector, colDist)
        }
        names(resVector) = colnames(mat1)
        return(resVector)
      }
      
      corMaindir = file.path(mainOutputDir, "distanceDistributions.refVsAlt", "corPearson")
      if ( ! file.exists( corMaindir ) ) {
        if ( ! dir.create( corMaindir, showWarnings=FALSE, recursive=TRUE )  ) corMaindir = mainOutputDir
      }
      
      # 20190529
      for (lib in whichLibs) {
        allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]][[lib]] = list()
      }
      
      for (useScaled in c(TRUE, FALSE)) {
        
        allelicCDPcomp = list()
        
        # 20190529
        for (lib in whichLibs) {
          print(lib)
          
          if (useScaled) {
            # allelicCDPcomp[[lib]] =  diag(cor(dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 
            #                               dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 
            #                               method="pearson", use="pairwise.complete"))
            allelicCDPcomp[[lib]] = cor.pearson.pw_matrix(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
                                                          dataMatrix.scaled[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]][[lib]][["scaled"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]][[lib]][["scaled"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          } else {
            # allelicCDPcomp[[lib]] =  diag(cor(dataMatrix.prop[[paste(lib, "ref", sep=".")]], 
            #                               dataMatrix.prop[[paste(lib, "alt", sep=".")]], 
            #                               method="pearson", use="pairwise.complete"))
            allelicCDPcomp[[lib]] = cor.pearson.pw_matrix(dataMatrix.prop[[paste(lib, "ref", sep=".")]],
                                                          dataMatrix.prop[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]][[lib]][["prop"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["Pearson"]][[lib]][["prop"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          }
          
          if (!all(is.na(allelicCDPcomp[[lib]]))) {
            
            dp <- density(allelicCDPcomp[[lib]], na.rm=TRUE)
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            maxDensity = min(15, max(dp$y, na.rm=TRUE))
            
            titleText = paste("Histogram of ref vs alt Pearson correlations\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltCorDensity.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid", xlim=c(-1,1),
                 main=titleText, ylim=c(0, maxDensity))
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
          }
          
        }
      }
      
      
      
      
      
      
      
      
      
      
      ##################################################
      # Distance between ref and alt alleles
      # https://stackoverflow.com/questions/35106567/how-to-calculate-euclidean-distance-between-two-matrices-in-r#35107198
      
      cat("\n\nEuclidean distance between allelic CDPs\n")
      
      allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]] = list()
      
      dist.pw_matrix <- function(mat1, mat2){
        resVector = c()
        for(i in 1:ncol(mat1)) {
          colDist = euclidean_distance(mat1[,i], mat2[,i])
          resVector = c(resVector, colDist)
        }
        names(resVector) = colnames(mat1)
        return(resVector)
      }
      
      corMaindir = file.path(mainOutputDir, "distanceDistributions.refVsAlt", "eucDist")
      if ( ! file.exists( corMaindir ) ) {
        if ( ! dir.create( corMaindir, showWarnings=FALSE, recursive=TRUE )  ) corMaindir = mainOutputDir
      }
      
      # 20190529
      for (lib in whichLibs) {
        allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]][[lib]] = list()
      }
      
      for (useScaled in c(TRUE, FALSE)) {
        
        allelicCDPcomp = list()
        
        # 20190529
        for (lib in whichLibs) {
          print(lib)
          
          if (useScaled) {
            # allelicCDPcomp[[lib]] =  as.vector(pdist(dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 
            #                                      dataMatrix.scaled[[paste(lib, "alt", sep=".")]])@dist)
            # allelicCDPcomp[[lib]] = diag(outer(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
            #                                dataMatrix.scaled[[paste(lib, "alt", sep=".")]],
            #                                Vectorize(euclidean_distance)))
            allelicCDPcomp[[lib]] = dist.pw_matrix(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
                                                   dataMatrix.scaled[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]][[lib]][["scaled"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]][[lib]][["scaled"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          } else {
            # allelicCDPcomp[[lib]] =  as.vector(pdist(dataMatrix.prop[[paste(lib, "ref", sep=".")]], 
            #                                      dataMatrix.prop[[paste(lib, "alt", sep=".")]])@dist)
            # allelicCDPcomp[[lib]] = diag(outer(dataMatrix.prop[[paste(lib, "ref", sep=".")]],
            #                                dataMatrix.prop[[paste(lib, "alt", sep=".")]],
            #                                Vectorize(euclidean_distance)))
            allelicCDPcomp[[lib]] = dist.pw_matrix(dataMatrix.prop[[paste(lib, "ref", sep=".")]],
                                                   dataMatrix.prop[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]][[lib]][["prop"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["eucdist"]][[lib]][["prop"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          }
          
          if (!all(is.na(allelicCDPcomp[[lib]]))) {
            
            dp <- density(allelicCDPcomp[[lib]], na.rm=TRUE)
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            maxDensity = min(15, max(dp$y, na.rm=TRUE))
            
            titleText = paste("Histogram of ref vs alt euclidean distances\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltDistanceDensity.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid",
                 main=titleText, ylim=c(0, maxDensity))
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
            
            
            ################
            # True histogram
            maxDistance = ceiling(max(allelicCDPcomp[[lib]], na.rm=TRUE))
            
            # maxDistance = 1
            dp <- hist(allelicCDPcomp[[lib]], plot=FALSE, breaks=seq(0, maxDistance, 
                                                                     ifelse( useScaled, 0.1, 0.1 )))
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            # maxDensity = length(na.omit(allelicCDPcomp[[lib]])) * 0.5
            maxDensity = max(dp$counts, na.rm=TRUE)
            
            titleText = paste("Histogram of ref vs alt binarizdd euclidean distances\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltDistanceHistogram.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid", xlim=c(0,maxDistance),
                 main=titleText, ylim=c(0, maxDensity), 
                 ylab="counts", xlab="Euclidean distance")
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
            
          }
        }
      }  
      
      
      
      
      
      
      
      
      
      
      ##################################################
      # Distance between ref and alt alleles - binarized
      # https://stackoverflow.com/questions/35106567/how-to-calculate-euclidean-distance-between-two-matrices-in-r#35107198
      
      cat("\n\nBinarized distance between allelic CDPs\n")
      
      allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]] = list()
      
      bindist.pw_matrix <- function(mat1, mat2){
        resVector = c()
        for(i in 1:ncol(mat1)) {
          binMat1 = mat1[,i]
          colDist = dist(rbind(binarize_vector(mat1[,i]), binarize_vector(mat2[,i])), "binary")
          resVector = c(resVector, colDist)
        }
        names(resVector) = colnames(mat1)
        return(resVector)
      }
      
      corMaindir = file.path(mainOutputDir, "distanceDistributions.refVsAlt", "binEucDist")
      if ( ! file.exists( corMaindir ) ) {
        if ( ! dir.create( corMaindir, showWarnings=FALSE, recursive=TRUE )  ) corMaindir = mainOutputDir
      }
      
      # 20190529
      for (lib in whichLibs) {
        allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]][[lib]] = list()
      }
      
      for (useScaled in c(TRUE, FALSE)) {
        
        allelicCDPcomp = list()
        
        # 20190529
        for (lib in whichLibs) {
          print(lib)
          
          if (useScaled) {
            # allelicCDPcomp[[lib]] =  as.vector(pdist(dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 
            #                                      dataMatrix.scaled[[paste(lib, "alt", sep=".")]])@dist)
            # allelicCDPcomp[[lib]] = diag(outer(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
            #                                dataMatrix.scaled[[paste(lib, "alt", sep=".")]],
            #                                Vectorize(euclidean_distance)))
            allelicCDPcomp[[lib]] = bindist.pw_matrix(dataMatrix.scaled[[paste(lib, "ref", sep=".")]],
                                                      dataMatrix.scaled[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]][[lib]][["scaled"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]][[lib]][["scaled"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          } else {
            # allelicCDPcomp[[lib]] =  as.vector(pdist(dataMatrix.prop[[paste(lib, "ref", sep=".")]], 
            #                                      dataMatrix.prop[[paste(lib, "alt", sep=".")]])@dist)
            # allelicCDPcomp[[lib]] = diag(outer(dataMatrix.prop[[paste(lib, "ref", sep=".")]],
            #                                dataMatrix.prop[[paste(lib, "alt", sep=".")]],
            #                                Vectorize(euclidean_distance)))
            allelicCDPcomp[[lib]] = bindist.pw_matrix(dataMatrix.prop[[paste(lib, "ref", sep=".")]],
                                                      dataMatrix.prop[[paste(lib, "alt", sep=".")]])
            allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]][[lib]][["prop"]] = allelicCDPcomp[[lib]]
            names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["bindist"]][[lib]][["prop"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
          }
          
          if (! all(is.na(allelicCDPcomp[[lib]]))) {
            
            dp <- density(allelicCDPcomp[[lib]], na.rm=TRUE)
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            maxDensity = min(15, max(c(dp$y), na.rm=TRUE))
            maxDistance = 1
            
            titleText = paste("Histogram of ref vs alt binary distances\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltDistanceDensity.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid", xlim=c(0,maxDistance),
                 main=titleText, ylim=c(0, maxDensity))
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
            
            
            ################
            # True histogram
            
            # maxDistance = max(5, ceiling(max(allelicCDPcomp[[lib]], na.rm=TRUE)))
            maxDistance = 1
            dp <- hist(allelicCDPcomp[[lib]], plot=FALSE, breaks=seq(0, maxDistance, 0.05))
            # alleliclib = paste(lib, "ref", sep=".")
            # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
            # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
            # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
            # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
            # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
            # maxDensity = length(na.omit(allelicCDPcomp[[lib]])) * 0.5
            maxDensity = max(dp$counts, na.rm=TRUE)
            
            titleText = paste("Histogram of ref vs alt binarizdd euclidean distances\n", 
                              lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                              " - ", bins2Use, "\n", whichChr, sep="")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            PDFfilename = file.path(corMaindir, paste(projectName, ".refVsAltDistanceHistogram.", nameText, ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dp, lwd=1, lty="solid", xlim=c(0,maxDistance),
                 main=titleText, ylim=c(0, maxDensity), 
                 ylab="counts", xlab="binary distance")
            # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
            # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
            # lines(dp0, col=tpCols["d0"], lwd=2)
            # lines(dp3, col=tpCols["d3"], lwd=2)
            # lines(dp7, col=tpCols["d7"], lwd=2)
            # lines(dp11, col=tpCols["d11"], lwd=2)
            legend("topleft", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                       # paste( "d0 (n = ", nd0, ")", sep=""),
                                       # paste( "d3 (n = ", nd3, ")", sep=""),
                                       # paste( "d7 (n = ", nd7, ")", sep=""),
                                       # paste( "d11 (n = ", nd11, ")", 
                                       sep=""),
                   # fill=c("black", tpCols), 
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
          }
        }  
      }
      
      
      
      
      
      
      
      
      
      
      ##################################################
      # 20204021 - 20200511
      # Long-range to Mid-range Difference (LMD)
      
      cat("\n\nlong-to-mid difference per allelic CDPs\n")
      
      if (bins2Use == "1kb-100Mb.merge10") {
        
        ##############
        ##############
        ##############
        #
        # Testing LMD strategies
        # And regions
        # F121_d0 vs F121_NPC and Patski
        
        testing = FALSE
        if (testing) {
          
          LRbins = c( "102", "112", "122" )
          MRbins = c(  "52", "62", "72" )
          
          if (useScaled) {
            
            #####
            # corr
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              # mmCDP.ref = apply( dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              # mmCDP.alt = apply( dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              mmCDP.ref = dataMatrix.scaled[[paste(lib, "ref", sep=".")]]
              mmCDP.alt = dataMatrix.scaled[[paste(lib, "alt", sep=".")]]
              
              hist( diag(cor(mmCDP.ref[c(LRbins, MRbins),], 
                             mmCDP.alt[c(LRbins, MRbins),])) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile(  diag(cor(mmCDP.ref[c(LRbins, MRbins),], mmCDP.alt[c(LRbins, MRbins),])), 
                                            p=0.1 )
              } else {
                print( c(lib, 
                         round( sum( diag(cor(mmCDP.ref[c(LRbins, MRbins),], mmCDP.alt[c(LRbins, MRbins),])) <= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }      
              
            } 
            # [1] "F121_NPC" "73.23"   
            # [1] "Patski" "61.8"  
            
            
            #####
            # dist
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              mmCDP.ref = apply( dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              # mmCDP.ref = dataMatrix.scaled[[paste(lib, "ref", sep=".")]]
              # mmCDP.alt = dataMatrix.scaled[[paste(lib, "alt", sep=".")]]
              
              hist( dist.pw_matrix(mmCDP.ref[c(LRbins, MRbins),], mmCDP.alt[c(LRbins, MRbins),]) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile(  dist.pw_matrix(mmCDP.ref[c(LRbins, MRbins),], mmCDP.alt[c(LRbins, MRbins),]), 
                                            p=0.9 )
              } else {
                print( c(lib, 
                         round( sum( dist.pw_matrix(mmCDP.ref[c(LRbins, MRbins),], mmCDP.alt[c(LRbins, MRbins),]) >= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }
              
            }
            # [1] "F121_NPC" "79.08"   
            # [1] "Patski" "67.54" 
            
            
            #####
            # Diff of diffs
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              mmCDP.ref = apply( dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              par(mfrow=c(1,2))
              hist( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) )
              hist( ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) )
              
              par(mfrow=c(1,1))
              hist( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                           ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile( abs( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                                               ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) ), 
                                           p=0.9 )
              } else {
                print( c(lib, 
                         round( sum( abs( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                                         ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) ) >= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }
              
            }
            # [1] "F121_NPC" "82.46"   
            # [1] "Patski" "78.64" 
            
            
            #####
            # Ratio of ratios
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              mmCDP.ref = apply( dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              par(mfrow=c(1,2))
              hist( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) )
              hist( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) )
              
              par(mfrow=c(1,1))
              hist( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                      ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile( abs( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                                             ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) ), 
                                           p=0.9 )
              } else {
                print( c(lib, 
                         round( sum( abs( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                                         ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) ) >= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }
              
            }
            # [1] "F121_NPC" "74.15"   
            # [1] "Patski" "64.72" 
            
            
          } else {
            
            
            #####
            # Diff of diffs
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              mmCDP.ref = apply( 10^dataMatrix.prop[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( 10^dataMatrix.prop[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              # mmCDP.ref = 10^dataMatrix.prop[[paste(lib, "ref", sep=".")]]
              # mmCDP.alt = 10^dataMatrix.prop[[paste(lib, "alt", sep=".")]]
              
              par(mfrow=c(1,2))
              hist( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) )
              hist( ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) )
              
              par(mfrow=c(1,1))
              hist( ( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                       ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) ) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile( abs( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                                               ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) ), 
                                           p=0.9 )
              } else {
                print( c(lib, 
                         round( sum( ( ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                                         ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) ) >= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }
              
            }
            # [1] "F121_NPC" "82.46"   
            # [1] "Patski" "78.54"
            
            # w/o min-max normalization:
            # [1] "F121_NPC" "80.92"   
            # [1] "Patski" "76.55" 
            
            
            #####
            # Ratio of ratios
            
            for ( lib in c( "F121_d0",  "F121_NPC", "Patski") ) {
              
              mmCDP.ref = apply( 10^dataMatrix.prop[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( 10^dataMatrix.prop[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              par(mfrow=c(1,2))
              hist( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) )
              hist( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) )
              
              par(mfrow=c(1,1))
              hist( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                      ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) )
              
              if (lib == "F121_d0") {
                F121d010pctFPR = quantile( abs( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                                             ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) ), 
                                           p=0.9 )
              } else {
                print( c(lib, 
                         round( sum( abs( log2( ( colSums(mmCDP.ref[LRbins,]) / colSums(mmCDP.ref[MRbins,]) ) /
                                         ( colSums(mmCDP.alt[LRbins,]) / colSums(mmCDP.alt[MRbins,]) ) ) ) >= F121d010pctFPR ) / 
                                  ncol(mmCDP.ref) * 100, 2) ) )
              }
              
            }
            # [1] "F121_NPC" "74.15"   
            # [1] "Patski" "64.67"
            
          }
          
        } # testing
        
        ##############
        ##############
        ##############
        
        
        LMDmaindir = file.path(mainOutputDir, "distanceDistributions.refVsAlt", "LMD")
        if ( ! file.exists( LMDmaindir ) ) {
          if ( ! dir.create( LMDmaindir, showWarnings=FALSE, recursive=TRUE )  ) LMDmaindir = mainOutputDir
        }
        
        min_max_norm <- function(x) {
          return((x- min(x)) /(max(x)-min(x)))
        }
        
        LRbins = c( "102", "112", "122" )
        MRbins = c(  "52", "62", "72" )
        
        # which( as.character( which(binBoolIdx) ) %in% LRbins )
        LRdescIdx = which( as.character( which(binBoolIdx)[altSeqIdx] ) %in% LRbins )
        MRdescIdx = which( as.character( which(binBoolIdx)[altSeqIdx] ) %in% MRbins )
        
        binDistDesc[ c(LRdescIdx, tail(LRdescIdx, 1) + 1) ] 
        # "6.5Mb"  "15.4Mb" "36.6Mb" "87Mb"  
        
        binDistDesc[ c(MRdescIdx, tail(MRdescIdx, 1) + 1) ]
        # "84.99kb"  "202.14kb" "480.77kb" "1.1Mb"  
        
        allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]] = list()
        
        # 20190529
        for (lib in whichLibs) {
          allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]][[lib]] = list()
        }
        
        for (useScaled in c(TRUE, FALSE)) {
          
          allelicCDPcomp = list()
          
          # 20190529
          for (lib in whichLibs) {
            print(lib)
            
            # lib = "F121_d0"
            # lib = "F121_NPC"
            # lib = "Patski"
            
            if (useScaled) {
              
              #####
              # Diff of diffs
              
              # min max normalization
              mmCDP.ref = apply( dataMatrix.scaled[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( dataMatrix.scaled[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              allelicCDPcomp[[lib]] = ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) 
              
              allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]][[lib]][["scaled"]] = allelicCDPcomp[[lib]]
              names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]][[lib]][["scaled"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
              
            } else {
              
              # min max normalization
              mmCDP.ref = apply( 10^dataMatrix.prop[[paste(lib, "ref", sep=".")]], 2, min_max_norm ) 
              mmCDP.alt = apply( 10^dataMatrix.prop[[paste(lib, "alt", sep=".")]], 2, min_max_norm ) 
              
              #####
              # Diff of diffs
              
              colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) 
              
              allelicCDPcomp[[lib]] = ( colSums(mmCDP.ref[LRbins,]) - colSums(mmCDP.ref[MRbins,]) ) -
                ( colSums(mmCDP.alt[LRbins,]) - colSums(mmCDP.alt[MRbins,]) ) 
              
              allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]][[lib]][["prop"]] = allelicCDPcomp[[lib]]
              names(allelicCDPcompList[[nThreshLabel]][[whichChr]][[bins2Use]][["LMD"]][[lib]][["prop"]]) = sub("ref." , "", colnames(dataMatrix.scaled[[paste(lib, "ref", sep=".")]]))
            }
            
            if (! all(is.na(allelicCDPcomp[[lib]]))) {
              
              dp <- density(allelicCDPcomp[[lib]], na.rm=TRUE)
              # alleliclib = paste(lib, "ref", sep=".")
              # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
              # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
              # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
              # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
              # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
              maxDensity = min(15, max(c(dp$y), na.rm=TRUE))
              maxDistance = max(abs(dp$x), na.rm=TRUE)
              
              titleText = paste("Histogram of ref vs alt LMD differences\n", 
                                lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                                " - ", bins2Use, "\n", whichChr, sep="")
              nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
              PDFfilename = file.path(LMDmaindir, paste(projectName, ".refVsAltDistanceDensity.", nameText, ".pdf", sep=""))
              pdf(PDFfilename, height=7, width=7)
              defparmar = par()$mar
              par(mar=par()$mar+c(0,0,2,0))
              plot(dp, lwd=1, lty="solid", xlim=c(-maxDistance, maxDistance),
                   main=titleText, ylim=c(0, maxDensity))
              # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
              # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
              #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
              # lines(dp0, col=tpCols["d0"], lwd=2)
              # lines(dp3, col=tpCols["d3"], lwd=2)
              # lines(dp7, col=tpCols["d7"], lwd=2)
              # lines(dp11, col=tpCols["d11"], lwd=2)
              legend("topright", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                         # paste( "d0 (n = ", nd0, ")", sep=""),
                                         # paste( "d3 (n = ", nd3, ")", sep=""),
                                         # paste( "d7 (n = ", nd7, ")", sep=""),
                                         # paste( "d11 (n = ", nd11, ")", 
                                         sep=""),
                     # fill=c("black", tpCols), 
                     border=NA, bty='n')
              par(mar=defparmar)
              dev.off()
              
              
              ################
              # True histogram
              
              maxDistance = max(5, ceiling(max(abs(allelicCDPcomp[[lib]]), na.rm=TRUE)))
              # maxDistance = 1
              
              dp <- hist(allelicCDPcomp[[lib]], plot=FALSE, breaks=seq(-maxDistance, maxDistance, 0.05))
              # alleliclib = paste(lib, "ref", sep=".")
              # dp0 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d0"]); nd0 = sum(matrixColTPs[[alleliclib]]=="d0")
              # dp3 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d3"]); nd3 = sum(matrixColTPs[[alleliclib]]=="d3")
              # dp7 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d7"]); nd7 = sum(matrixColTPs[[alleliclib]]=="d7")
              # dp11 <- density(allelicCDPcomp[[lib]][matrixColTPs[[alleliclib]]=="d11"]); nd11 = sum(matrixColTPs[[alleliclib]]=="d11")
              # maxDensity = min(12, max(c(dp$y, dp0$y, dp3$y, dp7$y, dp11$y), na.rm=TRUE))
              # maxDensity = length(na.omit(allelicCDPcomp[[lib]])) * 0.5
              maxDensity = max(dp$counts, na.rm=TRUE)
              
              titleText = paste("Histogram of ref vs alt LMD differences\n", 
                                lib, " contact decay ", ifelse(useScaled, "scaled", "prop"), 
                                " - ", bins2Use, "\n", whichChr, sep="")
              nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
              PDFfilename = file.path(LMDmaindir, paste(projectName, ".refVsAltDistanceHistogram.", nameText, ".pdf", sep=""))
              pdf(PDFfilename, height=7, width=7)
              defparmar = par()$mar
              par(mar=par()$mar+c(0,0,2,0))
              plot(dp, lwd=1, lty="solid", xlim=c(-maxDistance, maxDistance),
                   main=titleText, 
                   ylab="counts", xlab="allelic LMD difference")
              # abline(v = mean(dp$x, na.rm=TRUE), lty="dashed")
              # text(mean(dp$x, na.rm=TRUE), maxDensity*0.05, pos=2,
              #      paste("mean = ",  round(mean(dp$x, na.rm=TRUE), 2), sep=""))
              # lines(dp0, col=tpCols["d0"], lwd=2)
              # lines(dp3, col=tpCols["d3"], lwd=2)
              # lines(dp7, col=tpCols["d7"], lwd=2)
              # lines(dp11, col=tpCols["d11"], lwd=2)
              legend("topright", legend=c(paste("all", " (n = ", dim(dataMatrix.scaled[[paste(lib, "ref", sep=".")]])[2], ")", sep=""),
                                         # paste( "d0 (n = ", nd0, ")", sep=""),
                                         # paste( "d3 (n = ", nd3, ")", sep=""),
                                         # paste( "d7 (n = ", nd7, ")", sep=""),
                                         # paste( "d11 (n = ", nd11, ")", 
                                         sep=""),
                     # fill=c("black", tpCols), 
                     border=NA, bty='n')
              par(mar=defparmar)
              dev.off()
            }
          }  
        }
        
      } # if "1kb-100Mb.merge10"
 
    } # bins2Use
    
  } # whichChr
  
} # nThresh

allelicCDPcompList.all = allelicCDPcompList










##################################################
##################################################
##################################################
# Analyze and plot alleleic CDP data


##################################################
# Heatmaps of contact decay

for ( nThresh in names(allelicCDPcompList) ) {
  
  for ( whichChr in names(allelicCDPcompList[[nThresh]]) ) {
    
    subsupermainOutputDir = file.path(supermainOutputDir, 
                                      paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
    if ( ! file.exists( subsupermainOutputDir ) ) {
      if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
    }
    
    for( bins2Use in names(allelicCDPcompList[[nThresh]][[whichChr]]) ) {
      
      mainOutputDir = file.path(subsupermainOutputDir, 
                                bins2Use)
      if ( ! file.exists( mainOutputDir ) ) {
        if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
      }
      
      cat("\n\nContact decay heat maps\n")
      
      heatMaindir = file.path(mainOutputDir, "decayMaps")
      if ( ! file.exists( heatMaindir ) ) {
        if ( ! dir.create( heatMaindir, showWarnings=FALSE, recursive=TRUE )  ) heatMaindir = mainOutputDir
      }
      
      matrix2Plot.combined = list()
      matrix2Plot.bin.combined = list()
      
      for (lib in whichLibs) {
        print(lib)
        
        matrix2Plot.combined[[lib]] = list()
        matrix2Plot.bin.combined[[lib]] = list()
        
        for (useScaled in c(TRUE, FALSE)) {
          matrix2Plot.combined[[lib]][[ifelse(useScaled, "scaled", "prop")]] = c()
          matrix2Plot.bin.combined[[lib]][[ifelse(useScaled, "scaled", "prop")]] = c()
        }
      }
      
      for (alleliclib in whichAllelicLibs) {
        print(alleliclib)
        
        for (useScaled in c(TRUE, FALSE)) {
          
          # if (!useScaled) {
          #   matrix2Use = dataMatrix.prop[[alleliclib]]
          # } else {
          #   matrix2Use = dataMatrix.scaled[[alleliclib]]
          # }
          if (useScaled) {
            matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
          } else {
            matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
          }
          
          ##################################################
          # Sort cells by timepoint and coverage
          
          # head(matrixCsums.PF[[whichChr]][[alleliclib]][ covOrd[[alleliclib]] ])
          # tail(matrixCsums.PF[[whichChr]][[alleliclib]][ covOrd[[alleliclib]] ])
          
          matrix2Plot = matrix2Use[ ,covOrd[[alleliclib]][timepointThenCovOrd[[alleliclib]]] ]
          cellLabels = rep("", ncol(matrix2Plot))
          cellLabels[1] = colnames(matrix2Plot)[1]
          cellLabels[length(cellLabels)] = colnames(matrix2Plot)[length(cellLabels)]
          
          rownames(matrix2Plot) = paste(alleliclib, rownames(matrix2Plot), sep=".")
          
          TPcolbar = matrixColTPcols[[alleliclib]][covOrd[[alleliclib]]][timepointThenCovOrd[[alleliclib]]]
          covcolbar = covCols[[alleliclib]][covOrd[[alleliclib]]][timepointThenCovOrd[[alleliclib]]]
          rColBars = t(cbind(TPcolbar, covcolbar))
          rownames(rColBars)=c("sample", "coverage")
          
          titleText = paste(whichChr, " ", alleliclib, " ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameText = paste(whichChr, alleliclib, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          print(dim(matrix2Plot))
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          cat(ifelse(useScaled, "scaled", "prop"), ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(useScaled) {
            cols2Use = heatCols2
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = heatColsV
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          if ( ! file.exists( heatOutdir ) ) {
            if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
          }
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=binDistDesc, labRow=cellLabels, margins=c(5,0),
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       symm=ifelse(useScaled,TRUE, FALSE),
                       symkey=ifelse(useScaled,TRUE, FALSE),
                       symbreaks=ifelse(useScaled, TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleText,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       lmat=rbind(c(0, 0, 5, 4), c(0,1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
          
          # Combine CDPs from both alleles
          matrix2Plot.combined[[sub("\\..*$", "", alleliclib)]][[ifelse(useScaled, "scaled", "prop")]] = rbind(
            matrix2Plot.combined[[sub("\\..*$", "", alleliclib)]][[ifelse(useScaled, "scaled", "prop")]], matrix2Plot
          )
          
          
          ##################################################
          # Binarized CDPs
          
          matrix2Plot = apply(matrix2Plot, 2, binarize_vector)
          
          titleText = paste(whichChr, " ", alleliclib, " binarized ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameText = paste(whichChr, alleliclib, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          print(dim(matrix2Plot))
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          cat(ifelse(useScaled, "scaled", "prop"), ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(useScaled) {
            cols2Use = c("lightblue", "darkred") # heatCols2[c(1,length(heatColsV))]
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = c("yellow", "darkblue")# heatColsV[c(1,length(heatColsV))] 
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          if ( ! file.exists( heatOutdir ) ) {
            if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
          }
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.binerized.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=binDistDesc, labRow=cellLabels, margins=c(5,0),
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       symm=ifelse(useScaled,TRUE, FALSE),
                       symkey=ifelse(useScaled,TRUE, FALSE),
                       symbreaks=ifelse(useScaled, TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleText,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       lmat=rbind(c(0, 0, 5, 4), c(0,1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
          
          # Combine CDPs from both alleles
          matrix2Plot.bin.combined[[sub("\\..*$", "", alleliclib)]][[ifelse(useScaled, "scaled", "prop")]] = rbind(
            matrix2Plot.bin.combined[[sub("\\..*$", "", alleliclib)]][[ifelse(useScaled, "scaled", "prop")]], matrix2Plot
          )
          
        }
      }
      
      
      ##################################################
      # Plot combined heatmaps
      
      for (lib in whichLibs) {
        print(lib)
        
        for (useScaled in c(TRUE, FALSE)) {
          
          matrix2Plot = matrix2Plot.combined[[lib]][[ifelse(useScaled, "scaled", "prop")]]
          print(dim(matrix2Plot))
          
          TPcolbar = matrixColTPcols[[paste(lib, "ref", sep=".")]][covOrd[[paste(lib, "ref", sep=".")]]][timepointThenCovOrd[[paste(lib, "ref", sep=".")]]]
          # covcolbar = covCols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
          # rColBars = t(cbind(TPcolbar, covcolbar))
          # rownames(rColBars)=c("sample", "coverage")
          rColBars = t(cbind(TPcolbar))
          rownames(rColBars)=c("sample")
          
          alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
          rownames(alleleBars)=c("allele")
          colnames(alleleBars)=rownames(matrix2Plot)
          
          titleText = paste(whichChr, " ", lib, " bothAlleles ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameText = paste(whichChr, lib, "bothAlleles", ifelse(useScaled, "scaled", "prop"), sep=".")
          
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          cat(ifelse(useScaled, "scaled", "prop"), ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(useScaled) {
            cols2Use = heatCols2
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = heatColsV
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          if ( ! file.exists( heatOutdir ) ) {
            if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
          }
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), labRow=cellLabels, margins=c(5,0),
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       ColSideColors=t(alleleBars), ColSideColorsSize=4,
                       symm=ifelse(useScaled,TRUE, FALSE),
                       symkey=ifelse(useScaled,TRUE, FALSE),
                       symbreaks=ifelse(useScaled, TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleText,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                       lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(1, 0.2, 5), lwid=c(0.2, 0.2, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
          
          ##################################################
          # Combined binarized CDPs
          
          matrix2Plot = matrix2Plot.bin.combined[[lib]][[ifelse(useScaled, "scaled", "prop")]]
          print(dim(matrix2Plot))
          
          TPcolbar = matrixColTPcols[[paste(lib, "ref", sep=".")]][covOrd[[paste(lib, "ref", sep=".")]]][timepointThenCovOrd[[paste(lib, "ref", sep=".")]]]
          # covcolbar = covCols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
          # rColBars = t(cbind(TPcolbar, covcolbar))
          # rownames(rColBars)=c("sample", "coverage")
          rColBars = t(cbind(TPcolbar))
          rownames(rColBars)=c("sample")
          
          alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
          rownames(alleleBars)=c("allele")
          colnames(alleleBars)=rownames(matrix2Plot)
          
          titleText = paste(whichChr, " ", lib, " bothAlleles binarized ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameText = paste(whichChr, lib, "bothAlleles", ifelse(useScaled, "scaled", "prop"), sep=".")
          
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          cat(ifelse(useScaled, "scaled", "prop"), ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(useScaled) {
            cols2Use = c("lightblue", "darkred") # heatCols2[c(1,length(heatColsV))]
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = c("yellow", "darkblue")# heatColsV[c(1,length(heatColsV))] 
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          if ( ! file.exists( heatOutdir ) ) {
            if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
          }
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.binarized.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), labRow=cellLabels, margins=c(5,0),
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       ColSideColors=t(alleleBars), ColSideColorsSize=4,
                       symm=ifelse(useScaled,TRUE, FALSE),
                       symkey=ifelse(useScaled,TRUE, FALSE),
                       symbreaks=ifelse(useScaled, TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleText,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                       lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(1, 0.2, 5), lwid=c(0.2, 0.2, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
        } # scale or prop
      } # lib
      
    } # bins2Use
    
  } # whichChr
  
} # nThresh










##################################################
# Plot Decay curves

for ( nThresh in names(allelicCDPcompList) ) {
  
  for ( whichChr in names(allelicCDPcompList[[nThresh]]) ) {
    
    subsupermainOutputDir = file.path(supermainOutputDir, 
                                      paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
    if ( ! file.exists( subsupermainOutputDir ) ) {
      if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
    }
    
    for( bins2Use in names(allelicCDPcompList[[nThresh]][[whichChr]]) ) {
      
      mainOutputDir = file.path(subsupermainOutputDir, 
                                bins2Use)
      if ( ! file.exists( mainOutputDir ) ) {
        if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
      }  
      
      decayCurveMaindir = file.path(mainOutputDir, "decayCurves")
      if ( ! file.exists( decayCurveMaindir ) ) {
        if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
      }
      
      
      ######################
      # Plot all ref and alt curves side-by-side
      
      lwd2Use = 0.5
      alpha2Use=0.1
      
      for (useScaled in c(FALSE, TRUE)) {
        
        decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "paired")
        if ( ! file.exists( decayOutdir ) ) {
          if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
        }
        
        if (!useScaled) {
          ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
        }
        
        for (lib in whichLibs) {
          print(lib)
          
          meanCDP = list()
          
          titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.", nameText, ".refValt.pdf", sep=""))
          pdf(PDFfilename, height=4, width=12)
          par(mfrow=c(1,2))
          
          for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
            
            # if (!useScaled) {
            #   matrix2Use = dataMatrix.prop[[alleliclib]]
            # } else {
            #   matrix2Use = dataMatrix.scaled[[alleliclib]]
            # }
            if (useScaled) {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
            } else {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
            }
            
            nCells = ncol(matrix2Use)
            if (is.null(nCells) || nCells==0) next
            
            # alaph2Use =  1/nCells
            col2Use = setAlpha("grey", alpha2Use)
            
            for (iCell in 1:nCells) {
              y = matrix2Use[,iCell]
              if (iCell==1){
                if (useScaled) {
                  ylabText = "Scaled Cis counts"
                  ylims = c(-3, +3)
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = yat
                } else {
                  ylabText = "Proportion of Cis contact (log10)"
                  # ylims = c(-3, -1.5)
                  ylims = ylims.prop
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = signif(10^yat, 2)
                }
                plot( y, type="l", col=col2Use, lwd=lwd2Use,
                      main=alleliclib, xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                      ylim=ylims, xaxt="n", yaxt="n")
                axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                axis(2, at=yat, labels=ylabels)
              } else {
                lines(y, col=col2Use, lwd=lwd2Use)
              }
            }
            
            nCells = ncol(matrix2Use)
            if (is.null(nCells) || nCells==0) next
            
            y = rowMeans(matrix2Use, na.rm = TRUE)
            lines(y, col="black", lwd=1.5)
            meanCDP[[alleliclib]] = y
          } # alleliclib
          
          title(main=paste(titleText, " (n=", nCells, ")", sep=""), outer=TRUE, cex.main=1.2, line=-1)
          interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
          legend("topright", legend=paste("spearman = ", round(interAllelicCDPcorr, 2), sep=""), bty='n')
          
          dev.off()
          
        } # libs
        
      } # useScaled
      
      
      ######################
      # Plot all ref and alt curves overlaid
      
      lwd2Use = 0.5
      alpha2Use=0.01
      
      for (useScaled in c(FALSE, TRUE)) {
        
        decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "overlaid")
        if ( ! file.exists( decayOutdir ) ) {
          if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
        }
        
        if (!useScaled) {
          ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
        }
        
        for (lib in whichLibs) {
          print(lib)
          
          meanCDP = list()
          
          titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.overlaid.", nameText, ".refValt.pdf", sep=""))
          pdf(PDFfilename, height=4, width=6)
          
          nuPlotflag=TRUE
          for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
            print(alleliclib)
            
            # if (!useScaled) {
            #   matrix2Use = dataMatrix.prop[[alleliclib]]
            # } else {
            #   matrix2Use = dataMatrix.scaled[[alleliclib]]
            # }
            if (useScaled) {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
            } else {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
            }
            
            nCells = ncol(matrix2Use)
            if (is.null(nCells) || nCells==0) next
            
            # alaph2Use =  1/nCells
            col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
            
            for (iCell in 1:nCells) {
              y = matrix2Use[,iCell]
              if (iCell==1 && nuPlotflag){
                if (useScaled) {
                  ylabText = "Scaled Cis counts"
                  ylims = c(-3, +3)
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = yat
                } else {
                  ylabText = "Proportion of Cis contact (log10)"
                  # ylims = c(-3, -1.5)
                  ylims = ylims.prop
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = signif(10^yat, 2)
                }
                plot( y, type="l", col=col2Use, lwd=lwd2Use,
                      main=paste(titleText, " (n=", nCells, ")", sep=""),
                      xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                      ylim=ylims, xaxt="n", yaxt="n")
                axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                axis(2, at=yat, labels=ylabels)
              } else {
                lines(y, col=col2Use, lwd=lwd2Use)
              }
            }
            
            nCells = ncol(matrix2Use)
            if (is.null(nCells) || nCells==0) next
            
            y = rowMeans(matrix2Use, na.rm = TRUE)
            lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
            meanCDP[[alleliclib]] = y
            nuPlotflag = FALSE
            
          } # alleliclib
          
          alleliclib = sub("\\.alt$", ".ref", alleliclib)
          # if (!useScaled) {
          #   matrix2Use = dataMatrix.prop[[alleliclib]]
          # } else {
          #   matrix2Use = dataMatrix.scaled[[alleliclib]]
          # }
          if (useScaled) {
            matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
          } else {
            matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
          }
          y = rowMeans(matrix2Use, na.rm = TRUE)
          lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
          legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
          
          interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
          legend("topright", legend=paste("spearman = ", round(interAllelicCDPcorr, 2), sep=""), bty='n')
          
          dev.off() 
          
        } # libs
        
      } # useScaled
      
      
      ######################
      # 20200315
      # Plot all ref and alt curves overlaid 
      # Average of spearman correlations of CDPs rather than correlation of average CDPs.
      
      lwd2Use = 0.5
      alpha2Use=0.01
      
      for (useScaled in c(FALSE, TRUE)) {
        
        decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "overlaid.averageSpearman")
        if ( ! file.exists( decayOutdir ) ) {
          if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
        }
        
        if (!useScaled) {
          ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
        }
        
        for (lib in whichLibs) {
          print(lib)
          
          meanCDP = list()
          
          titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.overlaid.v2.", nameText, ".refValt.pdf", sep=""))
          pdf(PDFfilename, height=4, width=6)
          
          nuPlotflag=TRUE
          for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
            print(alleliclib)
            
            # if (!useScaled) {
            #   matrix2Use = dataMatrix.prop[[alleliclib]]
            # } else {
            #   matrix2Use = dataMatrix.scaled[[alleliclib]]
            # }
            if (useScaled) {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
            } else {
              matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
            }
            
            nCells = ncol(matrix2Use)
            if (is.null(nCells) || nCells==0) next
            
            # alaph2Use =  1/nCells
            col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
            
            for (iCell in 1:nCells) {
              y = matrix2Use[,iCell]
              if (iCell==1 && nuPlotflag){
                if (useScaled) {
                  ylabText = "Scaled Cis counts"
                  ylims = c(-3, +3)
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = yat
                } else {
                  ylabText = "Proportion of Cis contact (log10)"
                  # ylims = c(-3, -1.5)
                  ylims = ylims.prop
                  yat = seq(ylims[1], ylims[2], 0.5)
                  ylabels = signif(10^yat, 2)
                }
                plot( y, type="l", col=col2Use, lwd=lwd2Use,
                      main=paste(titleText, " (n=", nCells, ")", sep=""),
                      # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                      xlab="", ylab="",
                      ylim=ylims, xaxt="n", yaxt="n")
                axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                axis(2, at=yat, labels=ylabels)
              } else {
                lines(y, col=col2Use, lwd=lwd2Use)
              }
            }
            y = rowMeans(matrix2Use, na.rm = TRUE)
            lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
            meanCDP[[alleliclib]] = y
            nuPlotflag = FALSE
          } # alleliclib
          
          alleliclib = sub("\\.alt$", ".ref", alleliclib)
          # if (!useScaled) {
          #   matrix2Use = dataMatrix.prop[[alleliclib]]
          # } else {
          #   matrix2Use = dataMatrix.scaled[[alleliclib]]
          # }
          # if (useScaled) {
          #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
          # } else {
          #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
          # }
          # y = rowMeans(matrix2Use, na.rm = TRUE)
          # lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
          # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
          # interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
          
          # if (!useScaled) {
          #   matrix2Use.ref = dataMatrix.prop[[paste(lib, "ref", sep=".")]]
          #   matrix2Use.alt = dataMatrix.prop[[paste(lib, "alt", sep=".")]]
          # } else {
          #   matrix2Use.ref = dataMatrix.scaled[[paste(lib, "ref", sep=".")]]
          #   matrix2Use.alt = dataMatrix.scaled[[paste(lib, "alt", sep=".")]]
          # }
          if (useScaled) {
            matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "ref", sep=".")]]
            matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "alt", sep=".")]]
          } else {
            matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "ref", sep=".")]]
            matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "alt", sep=".")]]
          }

          CDPspearmans = c()
          nCells = ncol(matrix2Use)
          for(iCell in 1:nCells) {
            CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                              method='spearman', use='pairwise.complete.obs')
            CDPspearmans = c(CDPspearmans, CDPspearman)
          }
          
          meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
          medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
          
          legend("bottom", legend=c(
            paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
            paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
            bty='n')
          
          dev.off() 
          
        } # libs
        
      } # useScaled
 
    } # bins2Use
    
  } # whichChr
  
} # nThresh










# ---------------------
# 20210304 

##################################################
# Plot Decay curves by cluster

if (plotDecayCurvesByCluster || plotF121DecayCurvesByClusterJointly) {  
  
  for ( nThresh in names(allelicCDPcompList) ) {
    
    for ( whichChr in names(allelicCDPcompList[[nThresh]]) ) {
      
      subsupermainOutputDir = file.path(supermainOutputDir, 
                                        paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
      if ( ! file.exists( subsupermainOutputDir ) ) {
        if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
      }
      
      for( bins2Use in names(allelicCDPcompList[[nThresh]][[whichChr]]) ) {
        
        mainOutputDir = file.path(subsupermainOutputDir, 
                                  bins2Use)
        if ( ! file.exists( mainOutputDir ) ) {
          if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
        }  
        
        decayCurveMaindir = file.path(mainOutputDir, paste( "decayCurves_perClusterGroup.kCenters", kCenters, sep="") )
        if ( ! file.exists( decayCurveMaindir ) ) {
          if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
        }
        
        
        ######################
        # Plot all ref and alt curves side-by-side
        
        lwd2Use = 0.5
        alpha2Use=0.1
        
        for (useScaled in c(FALSE, TRUE)) {
          
          decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "paired")
          if ( ! file.exists( decayOutdir ) ) {
            if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
          }
          
          if (!useScaled) {
            ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
          }
          
          for (lib in whichLibs) {
            print(lib)
            
            # meanCDP = list()
            
            titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            
            PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.", nameText, ".refValt.pdf", sep=""))
            # pdf(PDFfilename, height=4, width=12)
            # par(mfrow=c(1,2))
            pdf(PDFfilename, height=8, width=4*kCenters)
            par(mfcol=c(2,1*kCenters))
            
            for ( kid in 1:kCenters) {
              
              meanCDP = list()
              
              for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
                
                # if (!useScaled) {
                #   matrix2Use = dataMatrix.prop[[alleliclib]]
                # } else {
                #   matrix2Use = dataMatrix.scaled[[alleliclib]]
                # }
                if (useScaled) {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
                } else {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
                }
                
                # head(colnames(matrix2Use))
                # [1] "TsixStop_NPC.alt.R4_lane1-AAGCCGGT.ACTGCTAC"
                # [2] "TsixStop_NPC.alt.R4_lane1-AAGCCGGT.AGGTCTCA"
                # [3] "TsixStop_NPC.alt.R4_lane1-AAGCCGGT.ATGCTTAA"
                # [4] "TsixStop_NPC.alt.R4_lane1-AAGCCGGT.TGGCAGTT"
                # [5] "TsixStop_NPC.alt.R4_lane1-ACCACCAC.AGGTCTCA"
                # [6] "TsixStop_NPC.alt.R4_lane1-ACCACCAC.GATTCTTA"
                # 
                # head( kClusters[[lib]]$cellID )
                # [1] R4_lane1-AAGCCGGT.AACGCGTA R4_lane1-AAGCCGGT.ACAGACTG
                # [3] R4_lane1-AAGCCGGT.ACTGCTAC R4_lane1-AAGCCGGT.ACTTGATA
                # [5] R4_lane1-AAGCCGGT.AGCAACCT R4_lane1-AAGCCGGT.AGCATATG
                
                # Get cells per cluster group:
                # Match revision:
                # x = c("a", "b", "c", "d")
                # y = c("a", "x", "d")
                #
                # XinY = match(x, y) # Where is x in y? What are the y indices for each x value?
                # # 1 NA NA  3
                #
                # x[!is.na(XinY)]
                # "a" "d"
                #
                # y[na.omit(XinY)]
                # "a" "d"
                # 20210310
                if (plotF121DecayCurvesByClusterJointly) {
                  CDPsInClustersIdx = match( colnames(matrix2Use), 
                                             sub( "(.*)\\.(.*)\\.(.*)", "\\2.\\3", kClusters[[lib]]$cellID ) )
                } else {
                  CDPsInClustersIdx = match( sub(paste0(alleliclib, "."), "", colnames(matrix2Use)), 
                                             kClusters[[lib]]$cellID )
                }
                CDPscoresInClusters = kClusters[[lib]]$clusterID[na.omit(CDPsInClustersIdx)]
                matrix2Use = matrix2Use[ , CDPscoresInClusters == kid]
                
                nCells = ncol(matrix2Use)
                if (is.null(nCells) && length(matrix2Use) == 14) nCells = 1
                if (is.null(nCells) || nCells==0) next
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha("grey", alpha2Use)
                
                for (iCell in 1:nCells) {
                  if (nCells == 1) {
                    y = matrix2Use
                  } else {
                    y = matrix2Use[,iCell]
                  }
                  if (iCell==1){
                    if (useScaled) {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste0( alleliclib, " (n=", nCells, ")" ),
                          xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                
                nCells = ncol(matrix2Use)
                if (is.null(nCells) && length(matrix2Use) == 14) nCells = 1
                if (is.null(nCells) || nCells==0) next
                
                if (nCells == 1) {
                  y = matrix2Use
                } else {
                  y = rowMeans(matrix2Use, na.rm = TRUE)
                }
                lines(y, col="black", lwd=1.5)
                meanCDP[[alleliclib]] = y
              } # alleliclib
              
              interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
              legend("topright", legend=paste("spearman = ", round(interAllelicCDPcorr, 2), "; n = ", nCells, sep=""), bty='n')
              
            } # kCenters
            
            title(main=paste(titleText), outer=TRUE, cex.main=1.2, line=-1)
            dev.off()
            
          } # libs
          
        } # useScaled
        
        
        ######################
        # Plot all ref and alt curves overlaid
        
        lwd2Use = 0.5
        alpha2Use=0.01
        
        for (useScaled in c(FALSE, TRUE)) {
          
          decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "overlaid")
          if ( ! file.exists( decayOutdir ) ) {
            if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
          }
          
          if (!useScaled) {
            ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
          }
          
          for (lib in whichLibs) {
            print(lib)
            
            # meanCDP = list()
            
            titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            
            PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.overlaid.", nameText, ".refValt.pdf", sep=""))
            # pdf(PDFfilename, height=4, width=6)
            pdf(PDFfilename, height=4, width=4*kCenters)
            par(mfcol=c(1,1*kCenters))
                
            for ( kid in 1:kCenters) {
              
              meanCDP = list()
              nuPlotflag=TRUE
              
              for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
                print(alleliclib)
                
                # if (!useScaled) {
                #   matrix2Use = dataMatrix.prop[[alleliclib]]
                # } else {
                #   matrix2Use = dataMatrix.scaled[[alleliclib]]
                # }
                if (useScaled) {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
                } else {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
                }
                
                # Get cells per cluster group:
                # Match revision:
                # x = c("a", "b", "c", "d")
                # y = c("a", "x", "d")
                #
                # XinY = match(x, y) # Where is x in y? What are the y indices for each x value?
                # # 1 NA NA  3
                #
                # x[!is.na(XinY)]
                # "a" "d"
                #
                # y[na.omit(XinY)]
                # "a" "d"
                # 20210310
                if (plotF121DecayCurvesByClusterJointly) {
                  CDPsInClustersIdx = match( colnames(matrix2Use), 
                                             sub( "(.*)\\.(.*)\\.(.*)", "\\2.\\3", kClusters[[lib]]$cellID ) )
                } else {
                  CDPsInClustersIdx = match( sub(paste0(alleliclib, "."), "", colnames(matrix2Use)), 
                                             kClusters[[lib]]$cellID )
                }
                CDPscoresInClusters = kClusters[[lib]]$clusterID[na.omit(CDPsInClustersIdx)]
                matrix2Use = matrix2Use[ , CDPscoresInClusters == kid]
                
                nCells = ncol(matrix2Use)
                if (is.null(nCells) && length(matrix2Use) == 14) nCells = 1
                if (is.null(nCells) || nCells==0) next
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                
                for (iCell in 1:nCells) {
                  if (nCells == 1) {
                    y = matrix2Use
                  } else {
                    y = matrix2Use[,iCell]
                  }
                  if (iCell==1 && nuPlotflag){
                    if (useScaled) {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " (n=", nCells, ")", sep=""),
                          xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                
                nCells = ncol(matrix2Use)
                if (is.null(nCells) && length(matrix2Use) == 14) nCells = 1
                if (is.null(nCells) || nCells==0) next
                
                if (nCells == 1) {
                  y = matrix2Use
                } else {
                  y = rowMeans(matrix2Use, na.rm = TRUE)
                }
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                meanCDP[[alleliclib]] = y
                nuPlotflag = FALSE
                
              } # alleliclib
              
              # alleliclib = sub("\\.alt$", ".ref", alleliclib)
              # # if (!useScaled) {
              # #   matrix2Use = dataMatrix.prop[[alleliclib]]
              # # } else {
              # #   matrix2Use = dataMatrix.scaled[[alleliclib]]
              # # }
              # if (useScaled) {
              #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
              # } else {
              #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
              # }
              # y = rowMeans(matrix2Use, na.rm = TRUE)
              # lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
              # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
              
              interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
              legend("topright", legend=paste("spearman = ", round(interAllelicCDPcorr, 2), "; n = ", nCells, sep=""), bty='n')
 
            } # kCenters

            dev.off() 
            
          } # libs

        } # useScaled
        
        
        ######################
        # 20200315
        # Plot all ref and alt curves overlaid 
        # Average of spearman correlations of CDPs rather than correlation of average CDPs.
        
        lwd2Use = 0.5
        alpha2Use=0.01
        
        for (useScaled in c(FALSE, TRUE)) {
          
          decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"), "overlaid.averageSpearman")
          if ( ! file.exists( decayOutdir ) ) {
            if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
          }
          
          if (!useScaled) {
            ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
          }
          
          for (lib in whichLibs) {
            print(lib)
            
            # meanCDP = list()
            
            titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
            nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=".")
            
            PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.overlaid.v2.", nameText, ".refValt.pdf", sep=""))
            # pdf(PDFfilename, height=4, width=6)
            pdf(PDFfilename, height=4, width=4*kCenters)
            par(mfcol=c(1,1*kCenters))
            
            for ( kid in 1:kCenters) {
              
              meanCDP = list()
              nuPlotflag=TRUE
              
              for (alleliclib in paste(lib, c("ref", "alt"), sep=".")) {
                print(alleliclib)
                
                # if (!useScaled) {
                #   matrix2Use = dataMatrix.prop[[alleliclib]]
                # } else {
                #   matrix2Use = dataMatrix.scaled[[alleliclib]]
                # }
                if (useScaled) {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
                } else {
                  matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
                }
                
                # Get cells per cluster group:
                # Match revision:
                # x = c("a", "b", "c", "d")
                # y = c("a", "x", "d")
                #
                # XinY = match(x, y) # Where is x in y? What are the y indices for each x value?
                # # 1 NA NA  3
                #
                # x[!is.na(XinY)]
                # "a" "d"
                #
                # y[na.omit(XinY)]
                # "a" "d"
                # 20210310
                if (plotF121DecayCurvesByClusterJointly) {
                  CDPsInClustersIdx = match( colnames(matrix2Use), 
                                             sub( "(.*)\\.(.*)\\.(.*)", "\\2.\\3", kClusters[[lib]]$cellID ) )
                } else {
                  CDPsInClustersIdx = match( sub(paste0(alleliclib, "."), "", colnames(matrix2Use)), 
                                             kClusters[[lib]]$cellID )
                }
                CDPscoresInClusters = kClusters[[lib]]$clusterID[na.omit(CDPsInClustersIdx)]
                matrix2Use = matrix2Use[ , CDPscoresInClusters == kid]
                
                nCells = ncol(matrix2Use)
                if (is.null(nCells) && length(matrix2Use) == 14) nCells = 1
                if (is.null(nCells) || nCells==0) next
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                
                for (iCell in 1:nCells) {
                  if (nCells == 1) {
                    y = matrix2Use
                  } else {
                    y = matrix2Use[,iCell]
                  }
                  if (iCell==1 && nuPlotflag){
                    if (useScaled) {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " (n=", nCells, ")", sep=""),
                          # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          xlab="", ylab="",
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                if (nCells == 1) {
                  y = matrix2Use
                } else {
                  y = rowMeans(matrix2Use, na.rm = TRUE)
                }
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                meanCDP[[alleliclib]] = y
                nuPlotflag = FALSE
              } # alleliclib
              
              alleliclib = sub("\\.alt$", ".ref", alleliclib)
              # if (!useScaled) {
              #   matrix2Use = dataMatrix.prop[[alleliclib]]
              # } else {
              #   matrix2Use = dataMatrix.scaled[[alleliclib]]
              # }
              # if (useScaled) {
              #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
              # } else {
              #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
              # }
              # y = rowMeans(matrix2Use, na.rm = TRUE)
              # lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
              # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
              # interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
              
              # if (!useScaled) {
              #   matrix2Use.ref = dataMatrix.prop[[paste(lib, "ref", sep=".")]]
              #   matrix2Use.alt = dataMatrix.prop[[paste(lib, "alt", sep=".")]]
              # } else {
              #   matrix2Use.ref = dataMatrix.scaled[[paste(lib, "ref", sep=".")]]
              #   matrix2Use.alt = dataMatrix.scaled[[paste(lib, "alt", sep=".")]]
              # }
              if (useScaled) {
                matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "ref", sep=".")]]
                matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "alt", sep=".")]]
              } else {
                matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "ref", sep=".")]]
                matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "alt", sep=".")]]
              }
              
              # Get cells per cluster group:
              # Match revision:
              # x = c("a", "b", "c", "d")
              # y = c("a", "x", "d")
              #
              # XinY = match(x, y) # Where is x in y? What are the y indices for each x value?
              # # 1 NA NA  3
              #
              # x[!is.na(XinY)]
              # "a" "d"
              #
              # y[na.omit(XinY)]
              # "a" "d"
              # 20210310
              if (plotF121DecayCurvesByClusterJointly) {
                CDPsInClustersIdx = match( colnames(matrix2Use), 
                                           sub( "(.*)\\.(.*)\\.(.*)", "\\2.\\3", kClusters[[lib]]$cellID ) )
              } else {
                CDPsInClustersIdx = match( sub(paste0(alleliclib, "."), "", colnames(matrix2Use)), 
                                           kClusters[[lib]]$cellID )
              }
              CDPscoresInClusters = kClusters[[lib]]$clusterID[na.omit(CDPsInClustersIdx)]
              matrix2Use.ref = matrix2Use.ref[ , CDPscoresInClusters == kid]
              matrix2Use.alt = matrix2Use.alt[ , CDPscoresInClusters == kid]
              
              CDPspearmans = c()
              nCells = ncol(matrix2Use.ref)
              if (is.null(nCells) && length(matrix2Use.ref) == 14) nCells = 1
              for(iCell in 1:nCells) {
                if (nCells == 1) {
                  CDPspearman = cor(matrix2Use.ref, matrix2Use.alt,
                                    method='spearman', use='pairwise.complete.obs')
                } else {
                  CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                                    method='spearman', use='pairwise.complete.obs')
                }
                CDPspearmans = c(CDPspearmans, CDPspearman)
              }
              
              meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
              medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
              
              legend("bottom", legend=c(
                paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
                paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
                bty='n')
              
              interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
              legend("topright", legend=paste("cluster ", kid, " spearman = ", round(interAllelicCDPcorr, 2), "; n = ", nCells, sep=""), bty='n')
              
            } # kCenters              
            
            dev.off() 
            
            } # libs

        } # useScaled
        
        
      } # bins2Use
      
    } # whichChr
    
  } # nThresh

  exit() # No further analysis 
  
} # Decay curves by cluster 

# 20210304 
# ---------------------










##################################################
# Projections

for ( nThresh in names(allelicCDPcompList) ) {
  
  for ( whichChr in names(allelicCDPcompList[[nThresh]]) ) {
    
    subsupermainOutputDir = file.path(supermainOutputDir, 
                                      paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
    if ( ! file.exists( subsupermainOutputDir ) ) {
      if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
    }
    
    for( bins2Use in names(allelicCDPcompList[[nThresh]][[whichChr]]) ) {
      
      mainOutputDir = file.path(subsupermainOutputDir, 
                                bins2Use)
      if ( ! file.exists( mainOutputDir ) ) {
        if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
      }  
      
      projectionMaindir = file.path(mainOutputDir, "projections")
      if ( ! file.exists( projectionMaindir ) ) {
        if ( ! dir.create( projectionMaindir, showWarnings=FALSE, recursive=TRUE )  ) projectionMaindir = mainOutputDir
      }
      
      
      ######################
      # PCA
      
      for (useScaled in c(TRUE)) {
        
        projectionOutdir = file.path(projectionMaindir, ifelse(useScaled, "scaled", "prop"))
        if ( ! file.exists( projectionOutdir ) ) {
          if ( ! dir.create( projectionOutdir, showWarnings=FALSE, recursive=TRUE )  ) projectionOutdir = projectionMaindir
        }
        
        for (lib in libs) {
          # for ( lib in grep("_rep", grep("^F121_", libs, value=TRUE), invert=TRUE, value=TRUE ) ) {
          print(lib)
          
          matrix2Use = c()
          cols2Use = c()
          shapes2Use = c()
          for ( allele in names(alleleCols) ) {
            
            lib.allele = paste( lib, allele, sep="." )
            
            # if (!useScaled) {
            #   matrix2Use = cbind(matrix2Use, dataMatrix.prop[[lib.allele]])
            #   cols2Use = c(cols2Use, rep(alleleCols[sub("^.*\\.", "", lib.allele)],
            #                              ncol(dataMatrix.prop[[lib.allele]])))
            #   shapes2Use = c(shapes2Use,  rep(alleleShapes[sub("^.*\\.", "", lib.allele)],
            #                                   ncol(dataMatrix.prop[[lib.allele]])))
            # } else {
            #   matrix2Use = cbind(matrix2Use, dataMatrix.scaled[[lib.allele]])
            #   cols2Use = c(cols2Use, rep(alleleCols[sub("^.*\\.", "", lib.allele)],
            #                              ncol(dataMatrix.scaled[[lib.allele]])))
            #   shapes2Use = c(shapes2Use,  rep(alleleShapes[sub("^.*\\.", "", lib.allele)],
            #                                   ncol(dataMatrix.scaled[[lib.allele]])))
            # }
            if (useScaled) {
              matrix2Use = cbind(matrix2Use, dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[lib.allele]])
              cols2Use = c(cols2Use, rep(alleleCols[sub("^.*\\.", "", lib.allele)],
                                         ncol(dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[lib.allele]])))
              shapes2Use = c(shapes2Use,  rep(alleleShapes[sub("^.*\\.", "", lib.allele)],
                                              ncol(dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[lib.allele]])))
            } else {
              matrix2Use = cbind(matrix2Use, dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[lib.allele]])
              cols2Use = c(cols2Use, rep(alleleCols[sub("^.*\\.", "", lib.allele)],
                                         ncol(dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[lib.allele]])))
              shapes2Use = c(shapes2Use,  rep(alleleShapes[sub("^.*\\.", "", lib.allele)],
                                              ncol(dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[lib.allele]])))
            }              
            
          }
          
          
          ######################
          # PCA
          
          thisProjectionOutdir = file.path(projectionOutdir, "PCA")
          if ( ! file.exists( thisProjectionOutdir ) ) {
            if ( ! dir.create( thisProjectionOutdir, showWarnings=FALSE, recursive=TRUE )  ) thisProjectionOutdir = projectionOutdir
          }
          
          prcompRes <- prcomp(t(matrix2Use), scale = TRUE)
          
          titleText = paste(lib, whichChr, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(lib, whichChr, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          PDFfilename = file.path(thisProjectionOutdir, paste(projectName, ".PCA.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=6)
          
          plot(prcompRes$x[, "PC1"],
               prcompRes$x[, "PC2"],
               col=cols2Use,
               pch=20, # shapes2Use,
               xlab="PC1",
               ylab="PC2",
               main=paste0(titleText, " (n=", ncol(matrix2Use), ")" ) )
          legend("topleft", legend=paste0( lib, " ", names(alleleCols) ), fill=alleleCols, bty='n')
          
          dev.off()
          
          
          ######################
          # UMAP
          
          thisProjectionOutdir = file.path(projectionOutdir, "UMAP")
          if ( ! file.exists( thisProjectionOutdir ) ) {
            if ( ! dir.create( thisProjectionOutdir, showWarnings=FALSE, recursive=TRUE )  ) thisProjectionOutdir = projectionOutdir
          }
          
          require(umap)
          set.seed(123)
          
          # 20200421 
          # Use first 3 PCs
          umapRes <- umap(t(matrix2Use), method="naive")
          # umapRes <- umap(prcompRes$x[,1:10], method="naive")
          
          titleText = paste(lib, whichChr, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(lib, whichChr, ifelse(useScaled, "scaled", "prop"), sep=".")
          
          PDFfilename = file.path(thisProjectionOutdir, paste(projectName, ".UMAP.", nameText, ".pdf", sep=""))
          pdf(PDFfilename, height=6, width=6)
          
          plot(umapRes$layout[, 1],
               umapRes$layout[, 2],
               col=cols2Use,
               pch=20, # shapes2Use,
               xlab="Dim1",
               ylab="Dim2",
               main=paste0(titleText, " (n=", ncol(matrix2Use), ")" ) )
          legend("topright", legend=paste0( lib, " ", names(alleleCols) ), fill=alleleCols, bty='n')
          
          dev.off()
          
        } # libs
        
      } # useScaled










##################################################
# Plot allelic comparisons between chromosomes

useInterphaseOnly = FALSE
# useInterphaseOnly = TRUE

if(useInterphaseOnly) {
  allelicCDPcompList = allelicCDPcompList.interphase
} else {
  allelicCDPcompList = allelicCDPcompList.all
}

for ( aThresh in names(allelicCDPcompList) ) {
  
  compMaindir = file.path(supermainOutputDir, paste( "countThresh", aThresh, ".chrXvsChr1", 
                                                     ifelse(useInterphaseOnly, ".interphaseOnly", ""), 
                                                     sep="" ))
  if ( ! file.exists( compMaindir ) ) {
    if ( ! dir.create( compMaindir, showWarnings=FALSE, recursive=TRUE )  ) compMaindir = supermainOutputDir
  }
  
  for( aBins2Use in names(allelicCDPcompList[[aThresh]][["chrX"]]) ) {
    # for( aCompType in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]]) ) {
    for( aCompType in c( "Spearman", "LMD")) {
      
      compSubdir = file.path(compMaindir, aCompType)
      if ( ! file.exists( compSubdir ) ) {
        if ( ! dir.create( compSubdir, showWarnings=FALSE, recursive=TRUE )  ) compSubdir = compMaindir
      }
      
      for( aLib in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]) ) {
        # for( aCDPnormType in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[aLib]]) ) {
        for( aCDPnormType in c( "scaled" ) ) {  
          
          ################
          # Joint density plots
          
          # titleText = paste("Histogram of ref vs alt ", aCompType, "\n", 
          #                   aLib, " contact decay ", aCDPnormType, "\n",
          #                   "chr1 vs chrX ", aBins2Use, " ", aThresh, " contacts", sep="")
          # print(titleText) # }}}}}
          
          NAflag = TRUE
          dpList = list()
          for ( aChr in c("chr1", "chrX") ) {
            anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
            if (all(is.na(anAllelicCDPcomp)) || length(na.omit(anAllelicCDPcomp))<2) {
              NAflag = FALSE
              next 
            }
            dpList[[aChr]] = density(anAllelicCDPcomp, na.rm=TRUE)
          }
          
          if (NAflag) {
            
            maxDensity = max(c(dpList[[1]]$y, dpList[[2]]$y), na.rm=TRUE)
            if ( aCompType == "LMD" ) {
              distRange = c(-max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE), 
                            max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE))
            } else if ( grepl("dist", aCompType) ) {
              distRange = c(0, max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE))
            } else {
              distRange = c(min(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE),
                            max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE))
            }
            
            titleText = paste("Histogram of ref vs alt ", aCompType, "\n", 
                              aLib, " contact decay ", aCDPnormType, "\n",
                              "chr1 vs chrX ", aBins2Use, " ", aThresh, " contacts", sep="")
            nameText = paste(aLib, aBins2Use, aCDPnormType, aCompType, aThresh, sep=".")
            PDFfilename = file.path(compSubdir, paste(projectName, ".refVsAlt.", nameText, ".densityPlot.pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dpList[["chrX"]], lwd=1, lty="solid", xlim=distRange, ylim=c(0, maxDensity),
                 main=titleText, xlab="Spearman correlation between allelic CDPs", col="red")
            # abline(v = median(dpList[["chrX"]]$x, na.rm=TRUE), lty="dashed")
            # text(median(dpList[["chrX"]]$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("med: ",  round(median(dpList[["chrX"]]$x, na.rm=TRUE), 2), sep=""))
            lines(dpList[["chr1"]]$x, dpList[["chr1"]]$y, lwd=1, lty="solid", col="blue")
            # abline(v = median(dpList[["chr1"]]$x, na.rm=TRUE), lty="dashed", col="blue")
            # text(median(dpList[["chr1"]]$x, na.rm=TRUE), maxDensity*0.15, pos=2,
            #      paste("med: ",  round(median(dpList[["chr1"]]$x, na.rm=TRUE), 2), sep=""),
            #      col="blue")
            legend("topleft", legend=c(paste( "chrX n = ", dpList[["chrX"]]$n, sep=""),
                                        paste( "chr1 n = ", dpList[["chr1"]]$n, sep="")),
                   fill=c("red", "blue"),
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
            
            
            ################
            # Joint histograms
            
            if ( grepl("dist", aCompType) ) {
              
              maxDistance = ceiling(max(c(allelicCDPcompList[[aThresh]][["chr1"]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]], 
                                          allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]), 
                                        na.rm=TRUE))
              
              dpList = list()
              for ( aChr in c("chr1", "chrX") ) {
                anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
                dpList[[aChr]] = hist(anAllelicCDPcomp, plot=FALSE, 
                                      breaks=seq(0, maxDistance, ifelse(aCompType=="bindist", 0.05, 0.1)))
              }
              
              maxDensity = max(c(dpList[[1]]$counts, dpList[[2]]$counts), na.rm=TRUE)
              
              titleText = paste("Histogram of ref vs alt ", aCompType, "\n", 
                                aLib, " contact decay ", aCDPnormType, "\n",
                                "chr1 vs chrX ", aBins2Use, " ", aThresh, " contacts", sep="")
              nameText = paste(aLib, aBins2Use, aCDPnormType, aCompType, aThresh, sep=".")
              PDFfilename = file.path(compSubdir, paste(projectName, ".refVsAlt.", nameText, ".histogram.pdf", sep=""))
              pdf(PDFfilename, height=7, width=7)
              defparmar = par()$mar
              par(mar=par()$mar+c(0,0,2,0))
              plot(dpList[["chrX"]],  xlim=c(0, maxDistance), ylim=c(0, maxDensity),
                   main=titleText, ylab="counts", xlab="Distance between allelic CDPs",
                   col=rgb(t(col2rgb("red"))/255, alpha=0.5, maxColorValue = 1))
              # abline(v = median(dpList[[1]]$x, na.rm=TRUE), lty="dashed")
              # text(median(dpList[["chrX"]]$x, na.rm=TRUE), maxDensity*0.05, pos=2,
              #      paste("med: ",  round(median(dpList[["chrX"]]$x, na.rm=TRUE), 2), sep=""))
              plot(dpList[["chr1"]], xlim=c(0, maxDistance), col=rgb(t(col2rgb("blue"))/255, alpha=0.5, maxColorValue = 1), add=TRUE)
              # abline(v = median(dpList[["chr1"]]$x, na.rm=TRUE), lty="dashed", col="blue")
              # text(median(dpList[["chr1"]]$x, na.rm=TRUE), maxDensity*0.15, pos=2,
              #      paste("med: ",  round(median(dpList[["chr1"]]$x, na.rm=TRUE), 2), sep=""),
              #      col="blue")
              legend("topleft", legend=c(paste( "chrX n = ", sum(dpList[["chrX"]]$counts), sep=""),
                                          paste( "chr1 n = ", sum(dpList[["chr1"]]$counts), sep="")),
                     fill=c(rgb(t(col2rgb("red"))/255, alpha=0.5, maxColorValue = 1), 
                            rgb(t(col2rgb("blue"))/255, alpha=0.5, maxColorValue = 1)),
                     border=NA, bty='n')
              par(mar=defparmar)
              dev.off()
              
            } # if dist
            
          }
        }
      }
    }
  }
}










##################################################
# Calculate AUC scores: chrX vs chr1
# And XCI threshold calculations using individual sample thresholds.

useInterphaseOnly = FALSE
# useInterphaseOnly = TRUE

if(useInterphaseOnly) {
  allelicCDPcompList = allelicCDPcompList.interphase
  # bipartiteResults = bipartiteResults.interphase
} else {
  allelicCDPcompList = allelicCDPcompList.all
  # bipartiteResults = bipartiteResults.all
}

tenFPRrateThresholds = list()

nSteps = 100
allelicCDPaucList.all = list()

for ( aThresh in names(allelicCDPcompList) ) {
  
  compMaindir = file.path(supermainOutputDir, paste( "countThresh", aThresh, ".chrXvsChr1", 
                                                     ifelse(useInterphaseOnly, ".interphaseOnly", ""), 
                                                     sep="" ))
  if ( ! file.exists( compMaindir ) ) {
    if ( ! dir.create( compMaindir, showWarnings=FALSE, recursive=TRUE )  ) compMaindir = supermainOutputDir
  }
 
  # bipartiteCompMaindir = file.path(supermainOutputDir, paste( "countThresh", aThresh, ".compareCDPscores2bipartiteScores", 
  #                                                             ifelse(useInterphaseOnly, ".interphaseOnly", ""), 
  #                                                             sep="" ))
  # if ( ! file.exists( bipartiteCompMaindir ) ) {
  #   if ( ! dir.create( bipartiteCompMaindir, showWarnings=FALSE, recursive=TRUE )  ) bipartiteCompMaindir = supermainOutputDir
  # }
  
  tenFPRrateThresholds[[aThresh]] = list()
  
  for( aBins2Use in names(allelicCDPcompList[[aThresh]][["chrX"]]) ) {
    # for( aCompType in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]]) ) {
    
    tenFPRrateThresholds[[aThresh]][[aBins2Use]] = list()
    
    for( aCompType in c( "Spearman", "LMD" )) {
    
      compSubdir = file.path(compMaindir, aCompType)
      if ( ! file.exists( compSubdir ) ) {
        if ( ! dir.create( compSubdir, showWarnings=FALSE, recursive=TRUE )  ) compSubdir = compMaindir
      }
      
      tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]] = list()
      
      for( aLib in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]) ) {
        
        tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]] = list()
        
        # for( aCDPnormType in rev(names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        for( aCDPnormType in c( "scaled" ) ) {  
          
          NAflag = TRUE
          anAllelicCDPcomp = list()
          for ( aChr in c("chr1", "chrX") ) {
            anAllelicCDPcomp[[aChr]] = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
            if (all(is.na(anAllelicCDPcomp[[aChr]])) || length(na.omit(anAllelicCDPcomp[[aChr]]))<2) {
              NAflag = FALSE
              next 
            }
          }
          
          if (NAflag) {
            
            AUCtable = c()
            if ( aCompType == "LMD") {
              distRange = c(min(abs(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]])), na.rm=TRUE),
                            max(abs(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]])), na.rm=TRUE))
              print(distRange)
              aStep = abs(distRange[2]-distRange[1])/nSteps
              for(i in seq(nSteps, 0)) {
                aStepThresh = distRange[1] + (i*aStep)
                TPR = sum(abs(anAllelicCDPcomp[["chrX"]]) >= aStepThresh)/length(anAllelicCDPcomp[["chrX"]])
                FPR = sum(abs(anAllelicCDPcomp[["chr1"]]) >= aStepThresh)/length(anAllelicCDPcomp[["chr1"]])
                NumCells = paste( "(", sum(anAllelicCDPcomp[["chrX"]] >= aStepThresh),
                                  " / ", length(anAllelicCDPcomp[["chrX"]]), ")",  sep="")
                print(c(i, aStepThresh, TPR, FPR, NumCells))
                AUCtable = rbind(AUCtable, c(aStepThresh, FPR, TPR, NumCells))
              }
            } else if ( grepl("dist", aCompType) ) {
              distRange = c(min(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]]), na.rm=TRUE),
                            max(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]]), na.rm=TRUE))
              print(distRange)
              aStep = abs(distRange[2]-distRange[1])/nSteps
              for(i in seq(nSteps, 0)) {
                aStepThresh = distRange[1] + (i*aStep)
                TPR = sum(anAllelicCDPcomp[["chrX"]] >= aStepThresh)/length(anAllelicCDPcomp[["chrX"]])
                FPR = sum(anAllelicCDPcomp[["chr1"]] >= aStepThresh)/length(anAllelicCDPcomp[["chr1"]])
                NumCells = paste( "(", sum(anAllelicCDPcomp[["chrX"]] >= aStepThresh),
                                  " / ", length(anAllelicCDPcomp[["chrX"]]), ")",  sep="")
                print(c(i, aStepThresh, TPR, FPR, NumCells))
                AUCtable = rbind(AUCtable, c(aStepThresh, FPR, TPR, NumCells))
              }
            } else {
              distRange = c(min(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]]), na.rm=TRUE),
                            max(c(anAllelicCDPcomp[["chrX"]], anAllelicCDPcomp[["chr1"]]), na.rm=TRUE))
              print(distRange)
              aStep = abs(distRange[2]-distRange[1])/nSteps
              for(i in seq(0, nSteps)) {
                aStepThresh = distRange[1] + (i*aStep)
                TPR = sum(anAllelicCDPcomp[["chrX"]] <= aStepThresh)/length(anAllelicCDPcomp[["chrX"]])
                FPR = sum(anAllelicCDPcomp[["chr1"]] <= aStepThresh)/length(anAllelicCDPcomp[["chr1"]])
                NumCells = paste( "(", sum(anAllelicCDPcomp[["chrX"]] <= aStepThresh),
                                  " / ", length(anAllelicCDPcomp[["chrX"]]), ")",  sep="")
                print(c(i, aStepThresh, TPR, FPR, NumCells))
                AUCtable = rbind(AUCtable, c(aStepThresh, FPR, TPR, NumCells))
              }
              
            }
            colnames(AUCtable) = c("thresh", "FPR", "TPR", "NumCells")
            if ( aCompType == "LMD") {
              rownames(AUCtable) = round( rev( seq(distRange[1], distRange[2], aStep) ), 2)
            } else if ( grepl("dist", aCompType) ) {
              rownames(AUCtable) = round( rev( seq(distRange[1], distRange[2], aStep) ), 2)
            } else {
              rownames(AUCtable) = round( seq(distRange[1], distRange[2], aStep), 2) 
            }
            aucDesc = paste(aThresh, aBins2Use, aCDPnormType, aCompType, sep=".") 
            allelicCDPaucList.all[[paste(aLib, aucDesc, sep=".") ]] = AUCtable
            
            tenFPRrateIdx = which.min(abs(as.numeric(AUCtable[,"FPR"]) - 0.1))
            tenFPRrate = as.numeric(AUCtable[,"TPR"])[tenFPRrateIdx]
            tenFPRcellNum = AUCtable[,"NumCells"][tenFPRrateIdx]
            # 20201118
            # tenFPRrateThresh = as.numeric(rownames(AUCtable)[tenFPRrateIdx])
            tenFPRrateThresh = as.numeric(AUCtable[,"thresh"])[tenFPRrateIdx]
            
            tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]] = tenFPRrateThresh
            
            ################
            # Overlap analysis with biparite-called XCI cells
            
            # Match revision:
            # x = c("a", "b", "c", "d")
            # y = c("a", "x", "d")
            #
            # XinY = match(x, y) # Where is x in y? What are the y indices for each x value?
            # # 1 NA NA  3
            #
            # x[!is.na(XinY)]
            # "a" "d"
            #
            # y[na.omit(XinY)]
            # "a" "d"
            
            aChr = "chrX"
            anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
            if ( aCompType == "LMD") {
              XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
            } else if ( grepl("dist", aCompType) ) {
              XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
            } else {
              XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
            }
            # head(anAllelicCDPcomp)
            cat( aCompType, ":", print(sum(XCIboolIdx)), "\n", sep="")
            
            # CDPinRBImatchIdx = match( names(anAllelicCDPcomp), bipartiteResults$cellID )
            # print( sum( is.na( CDPinRBImatchIdx ) ) )
            # 
            # RBIinDPmatchIdx = match( bipartiteResults$cellID, names(anAllelicCDPcomp) )
            # print( sum( is.na( RBIinDPmatchIdx ) ) )
            # 
            # XCIbyCDPinRBImatchIdx = match( names(anAllelicCDPcomp)[XCIboolIdx],
            #                                bipartiteResults$cellID )
            # table(! is.na(XCIbyCDPinRBImatchIdx) ) 
            # 
            # XCIbyRBIinCPDmatchIdx = match( bipartiteResults$cellID[bipartiteResults$XCI],
            #                                names(anAllelicCDPcomp) )
            # table(! is.na(XCIbyRBIinCPDmatchIdx))
            # 
            # XCIbyCDPinXCIbyRBImatchIdx = match( names(anAllelicCDPcomp)[XCIboolIdx],
            #                                     bipartiteResults$cellID[bipartiteResults$XCI] )
            # table(! is.na(XCIbyCDPinXCIbyRBImatchIdx) ) 
            # 
            # XCIbyRBIinXCIbyCPDmatchIdx = match( bipartiteResults$cellID[bipartiteResults$XCI],
            #                                     names(anAllelicCDPcomp)[XCIboolIdx] )
            # table(! is.na(XCIbyRBIinXCIbyCPDmatchIdx))
            # 
            # titleText = paste(aLib, " - correlation of CDP comparisons with composite relative bipartite scores\n", aucDesc, sep="")
            # nameText = paste(aLib, aucDesc, "comparisonWithCompositeRelativeBipartiteScores", sep=".")
            # PDFfilename = file.path(bipartiteCompMaindir, paste(projectName, ".", nameText, 
            #                                                     ifelse(useInterphaseOnly, ".interphaseOnly", ""),
            #                                                     ".pdf", sep=""))
            # pdf(PDFfilename, height=8, width=8)
            # plot(anAllelicCDPcomp[!is.na(CDPinRBImatchIdx)], 
            #      bipartiteResults$XCIscore[na.omit(CDPinRBImatchIdx)],
            #      xlab="CDP comparison score", ylab="Composite relative bipartite score",
            #      main=titleText, pch=20, col="grey")
            # points(anAllelicCDPcomp[XCIboolIdx][!is.na(XCIbyCDPinRBImatchIdx)], 
            #        bipartiteResults$XCIscore[na.omit(XCIbyCDPinRBImatchIdx)],
            #        pch=20, col='blue')
            # points(anAllelicCDPcomp[na.omit(XCIbyRBIinCPDmatchIdx)], 
            #        bipartiteResults$XCIscore[bipartiteResults$XCI][!is.na(XCIbyRBIinCPDmatchIdx)],
            #        pch=20, col='red')
            # points(anAllelicCDPcomp[XCIboolIdx][!is.na(XCIbyCDPinXCIbyRBImatchIdx)], 
            #        bipartiteResults$XCIscore[bipartiteResults$XCI][na.omit(XCIbyCDPinXCIbyRBImatchIdx)],
            #        pch=20, col='black')
            # legend(ifelse( grepl("dist", aCompType), "topleft", "topright" ),
            #        legend = c( paste( "Spearman all (n=", sum(!is.na(CDPinRBImatchIdx)), "): ", 
            #                           round( cor(anAllelicCDPcomp[!is.na(CDPinRBImatchIdx)], 
            #                                      bipartiteResults$XCIscore[na.omit(CDPinRBImatchIdx)],
            #                                      use="pairwise.complete.obs", method="spearman"), 2),
            #                           sep=""),
            #                    paste( "XCI by CDP (n=", sum(!is.na(XCIbyCDPinRBImatchIdx)), ")", sep=""),
            #                    paste( "XCI by bipartite (n=", sum(!is.na(XCIbyRBIinCPDmatchIdx)), ")", sep=""),
            #                    paste( "XCI by both (n=", sum(!is.na(XCIbyCDPinXCIbyRBImatchIdx)), ")", sep="") ),
            #        fill=c("grey", "blue", "red", "black"),
            #        bty='n', cex=0.8)
            # dev.off()
            
            
            ################
            # Joint density plots
            
            NAflag = TRUE
            dpList = list()
            for ( aChr in c("chr1", "chrX") ) {
              anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
              if (all(is.na(anAllelicCDPcomp)) || length(na.omit(anAllelicCDPcomp))<2) {
                NAflag = FALSE
                next 
              }
              dpList[[aChr]] = density(anAllelicCDPcomp, na.rm=TRUE)
            }
            
            if (NAflag) {

              maxDensity = max(c(dpList[[1]]$y, dpList[[2]]$y), na.rm=TRUE)
              distRange = c(min(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE),
                            max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE))
              
              titleText = paste(aLib, " - CDP comparison distributions\n", 
                                aucDesc, " - chr1 vs chrX", sep="")
              nameText = paste(aLib, aucDesc, "densityPlot", sep=".")
              PDFfilename = file.path(compSubdir, paste(projectName, ".refVsAlt.withThresholds.", nameText, 
                                                         ifelse(useInterphaseOnly, ".interphaseOnly", ""),
                                                         ".pdf", sep=""))
              pdf(PDFfilename, height=7, width=7)
              defparmar = par()$mar
              par(mar=par()$mar+c(0,0,2,0))
              plot(dpList[["chrX"]], lwd=1, lty="solid", xlim=distRange, ylim=c(0, maxDensity),
                   main=titleText, col="red")
              # abline(v = median(dpList[["chrX"]]$x, na.rm=TRUE), lty="dashed")
              # text(median(dpList[["chrX"]]$x, na.rm=TRUE), maxDensity*0.05, pos=2,
              #      paste("med: ",  round(median(dpList[["chrX"]]$x, na.rm=TRUE), 2), sep=""))
              lines(dpList[["chr1"]]$x, dpList[["chr1"]]$y, lwd=1, lty="solid", col="blue")
              # abline(v = median(dpList[["chr1"]]$x, na.rm=TRUE), lty="dashed", col="blue")
              # text(median(dpList[["chr1"]]$x, na.rm=TRUE), maxDensity*0.15, pos=2,
              #      paste("med: ",  round(median(dpList[["chr1"]]$x, na.rm=TRUE), 2), sep=""),
              #      col="blue")
              abline(v=tenFPRrateThresh, lty="dashed")
              if ( aCompType == "LMD") abline(v=-tenFPRrateThresh, lty="dashed")
              text(tenFPRrateThresh, maxDensity*0.95, pos=2,
                   paste(round(tenFPRrate*100, 2), "% ", tenFPRcellNum, "\n",
                         "@10%FPR threshold (",  round(tenFPRrateThresh, 2), ")", sep=""),
                   col="black", cex=0.8)
              legend("topleft", legend=c(paste( "chrX n = ", dpList[["chrX"]]$n, sep=""),
                                         paste( "chr1 n = ", dpList[["chr1"]]$n, sep="")),
                     fill=c("red", "blue"),
                     border=NA, bty='n')
              par(mar=defparmar)
              dev.off()
              
            } # NAflag
            
          } # NAflag
        } # aCDPnormType
      } # aLib
    } # aCompType
  } # aBins2Use
} # aThresh










##################################################
# AUC plots

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}


plotAUCs <- function(AUClist, cols2Use=rainbow(length(AUClist)), titleText="", legcex=0.9, shortLeg=TRUE, noDetails=FALSE) {
  # AUClist = allelicCDPaucList
  # AUClist = allelicCDPaucList.subset
  auc = c()
  tenFPRrate = c()
  tenFPRcellNum = c()
  tenFPRrateThresh = c()
  for(i in seq(1, length(AUClist)) ) {
    auc[i] = simple_auc(as.numeric(AUClist[[i]][,"TPR"]), 
                        as.numeric(AUClist[[i]][,"FPR"]))
    tenFPRrateIdx = which.min(abs(as.numeric(AUClist[[i]][,"FPR"]) - 0.1))
    tenFPRrate[i] = as.numeric(AUClist[[i]][,"TPR"])[tenFPRrateIdx]
    tenFPRcellNum[i] = AUClist[[i]][,"NumCells"][tenFPRrateIdx]
    # 20201118
    # tenFPRrateThresh[i] = as.numeric(rownames(AUClist[[i]])[tenFPRrateIdx])
    tenFPRrateThresh[i] = as.numeric(AUClist[[i]][,"thresh"])[tenFPRrateIdx]
  }
  aucOrderIdx = order(auc, decreasing=TRUE)
  plot(0, 0, xlim = c(0,1), ylim = c(0,1), type="n",
       ylab="chrX (TPR)", 
       xlab="chr1 (FPR)",
       main=paste(titleText, "\nallelic CDP-based proportions of cells called inactive", sep=""))
  for(i in seq(1, length(AUClist))[aucOrderIdx] ) {
    lines(as.numeric(AUClist[[i]][,"FPR"]),  as.numeric(AUClist[[i]][,"TPR"]), 
          lwd=2, col = cols2Use[i], type='l', pch=20)
  }
  if (noDetails) {
    legNames = sub("\\..*$", "", names(AUClist))
  } else {
    legNames = names(AUClist)
  }
  if (shortLeg) {
    legend("bottomright", legend=paste( legNames, ": ", "auc=", round(auc, 2), sep="")[aucOrderIdx], 
           fill=cols2Use[aucOrderIdx], cex=legcex, bty='n')
  } else {
    if (noDetails) {
      legtext = legend=paste( legNames, ": ",
                              "auc=", round(auc, 2), " ",
                              "@10%FPR=", round(tenFPRrate, 2), " ", tenFPRcellNum,  sep="")[aucOrderIdx]
    } else {
      legtext = paste( legNames, ": ",
                       "auc=", round(auc, 2), " ",
                       "@10%FPR=", round(tenFPRrate, 2), " ", tenFPRcellNum, " ", 
                       "10%FPRthresh=", round(tenFPRrateThresh, 2),  sep="")[aucOrderIdx]
    }
    legend("bottomright", legend=legtext, 
           fill=cols2Use[aucOrderIdx], cex=legcex, bty='n')
  }
  abline(v=0.1, lty="dashed", col="red")
}


plotCalledXiCounts <- function(AUClist, cols2Use=rainbow(length(AUClist)), titleText="", legcex=0.9, shortLeg=TRUE) {
  # AUClist = allelicCDPaucList
  XiCellCount = list()
  # maxCellCount = c()
  thresholdCellCount = c()
  tenFPRrate = c()
  tenFPRcellNum = c()
  tenFPRrateThresh = c()
  for(i in seq(1, length(AUClist)) ) {
    XiCellCount[[i]] = as.numeric(sub(" / .*$", "", sub("^\\(", "", AUClist[[i]][,"NumCells"])))
    # maxCellCount[i] = max(XiCellCount[[i]], na.rm=TRUE)
    tenFPRrateIdx = which.min(abs(as.numeric(AUClist[[i]][,"FPR"]) - 0.1))
    tenFPRrate[i] = as.numeric(AUClist[[i]][,"TPR"])[tenFPRrateIdx]
    tenFPRcellNum[i] = AUClist[[i]][,"NumCells"][tenFPRrateIdx]
    # 20201118
    # tenFPRrateThresh[i] = as.numeric(rownames(AUClist[[i]])[tenFPRrateIdx])
    tenFPRrateThresh[i] = as.numeric(AUClist[[i]][,"thresh"])[tenFPRrateIdx]
    thresholdCellCount[i] = XiCellCount[[i]][tenFPRrateIdx]
  }
  orderIdx = order(thresholdCellCount, decreasing=TRUE)
  maxMaxCellCount = max(thresholdCellCount, na.rm=TRUE) * 1.2
  plot(0, 0, xlim = c(0,0.25), ylim = c(0,maxMaxCellCount), type="n",
       ylab="# Xi", 
       xlab="chr1 (FPR)",
       main=paste(titleText, "\nallelic CDP-based numbers of cells called inactive", sep=""))
  for(i in seq(1, length(AUClist))[orderIdx] ) {
    lines(as.numeric(AUClist[[i]][,"FPR"]), XiCellCount[[i]], col = cols2Use[i], type='l', pch=20)
  }
  abline(v=0.1, lty="dashed")
  if (shortLeg) {
    legend("topleft", legend=paste( names(AUClist), ": ", "#Xi=", round(thresholdCellCount, 2), sep="")[orderIdx], 
           fill=cols2Use[orderIdx], cex=legcex, bty='n')
  } else {
    legend("topleft", legend=paste( names(AUClist), ": ",
                                    "#Xi=", round(thresholdCellCount, 2), " ",
                                    "@10%FPR=", round(tenFPRrate, 2), " ", tenFPRcellNum, " ", 
                                    "10%FPRthresh=", round(tenFPRrateThresh, 2),  sep="")[orderIdx], 
           fill=cols2Use[orderIdx], cex=legcex, bty='n')
  }
}


plotCalledXiCountsTopN <- function(AUClist, cols2Use=rainbow(length(AUClist)), titleText="", legcex=0.9, shortLeg=TRUE, topN=25) {
  # AUClist = allelicCDPaucList
  XiCellCount = list()
  # maxCellCount = c()
  thresholdCellCount = c()
  tenFPRrate = c()
  tenFPRcellNum = c()
  tenFPRrateThresh = c()
  for(i in seq(1, length(AUClist)) ) {
    XiCellCount[[i]] = as.numeric(sub(" / .*$", "", sub("^\\(", "", AUClist[[i]][,"NumCells"])))
    # maxCellCount[i] = max(XiCellCount[[i]], na.rm=TRUE)
    tenFPRrateIdx = which.min(abs(as.numeric(AUClist[[i]][,"FPR"]) - 0.1))
    tenFPRrate[i] = as.numeric(AUClist[[i]][,"TPR"])[tenFPRrateIdx]
    tenFPRcellNum[i] = AUClist[[i]][,"NumCells"][tenFPRrateIdx]
    # 20201118
    # tenFPRrateThresh[i] = as.numeric(rownames(AUClist[[i]])[tenFPRrateIdx])
    tenFPRrateThresh[i] = as.numeric(AUClist[[i]][,"thresh"])[tenFPRrateIdx]    
    thresholdCellCount[i] = XiCellCount[[i]][tenFPRrateIdx]
  }
  orderIdx = order(thresholdCellCount, decreasing=TRUE)
  maxMaxCellCount = max(thresholdCellCount, na.rm=TRUE) * 1.2
  plot(0, 0, xlim = c(0,0.25), ylim = c(0,maxMaxCellCount), type="n",
       ylab="# Xi", 
       xlab="chr1 (FPR)",
       main=paste(titleText, "\nallelic CDP-based numbers of cells called inactive", sep=""))
  for(i in seq(1, length(AUClist))[orderIdx][1:topN] ) {
    lines(as.numeric(AUClist[[i]][,"FPR"]), XiCellCount[[i]], col = cols2Use[i], type='l', pch=20)
  }
  abline(v=0.1, lty="dashed")
  if (shortLeg) {
    legend("bottomright", legend=paste( names(AUClist), ": ", "#Xi=", round(thresholdCellCount, 2), sep="")[orderIdx][1:topN], 
           fill=cols2Use[orderIdx], cex=legcex, bty='n')
  } else {
    legend("bottomright", legend=paste( names(AUClist), ": ",
                                        "#Xi=", round(thresholdCellCount, 2), " ",
                                        "@10%FPR=", round(tenFPRrate, 2), " ", tenFPRcellNum, " ", 
                                        "10%FPRthresh=", round(tenFPRrateThresh, 2),  sep="")[orderIdx][1:topN], 
           fill=cols2Use[orderIdx], cex=legcex, bty='n')
  }
}


aucMaindir = file.path(supermainOutputDir, paste( "AUCanalysis", 
                                                  ".countThresh", aThresh, ".chrXvsChr1", 
                                                  ifelse(useInterphaseOnly, ".interphaseOnly", ""), 
                                                  sep="" ))
if ( ! file.exists( aucMaindir ) ) {
  if ( ! dir.create( aucMaindir, showWarnings=FALSE, recursive=TRUE )  ) aucMaindir = supermainOutputDir
}


##### AUCs #####

for( aCompType in c( "Spearman", "LMD" )) {
  
  aucSubdir = file.path(aucMaindir, aCompType)
  if ( ! file.exists( aucSubdir ) ) {
    if ( ! dir.create( aucSubdir, showWarnings=FALSE, recursive=TRUE )  ) aucSubdir = aucMaindir
  }
  
  allelicCDPaucList = allelicCDPaucList.all[ grepl(paste0(".", aCompType), names(allelicCDPaucList.all)) ]
  
  ### All curves
  PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.AUCs.all.pdf", sep=""))
  pdf(PDFfilename, height=10, width=10)
  plotAUCs(allelicCDPaucList, titleText="All", legcex=1, shortLeg=FALSE, noDetails=TRUE)
  dev.off()
  
  PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.counts.all.pdf", sep=""))
  pdf(PDFfilename, height=10, width=10)
  plotCalledXiCounts(allelicCDPaucList, titleText="All", legcex=1)
  dev.off()
  
  # 20200803
  if (! patski100 ) {
  
    ### Subsets of curves
    PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.AUCs.subsets.pdf", sep=""))
    pdf(PDFfilename, height=5*5, width=8*2)
    par(mfrow=c(5,2))
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^TsixStop", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="TsixStop cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^F121", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="F121 cells: countThresh100 scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^EBdiff_TsixStop", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="EBdiff TsixStop cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^EBdiff_F121", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="EBdiff F121 cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^TsixStop", names(allelicCDPaucList)) |
                                                   grepl("^EBdiff_TsixStop", names(allelicCDPaucList))   ]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="Both TsixStop cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^F121", names(allelicCDPaucList)) |
                                                   grepl("^EBdiff_F121", names(allelicCDPaucList))   ]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="Both F121 cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("[Dd]0", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="D0 cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("[Dd]0", names(allelicCDPaucList)) | 
                                                   grepl("d3", names(allelicCDPaucList)) ]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="D0 and D3 cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("[Dd]7", names(allelicCDPaucList)) | 
                                                   grepl("d11", names(allelicCDPaucList)) ]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="D7 and D11 cells: scaled Spearman", shortLeg=FALSE)
    
    allelicCDPaucList.subset = allelicCDPaucList[grepl("NPC", names(allelicCDPaucList)) | 
                                                   grepl("Patski", names(allelicCDPaucList)) |
                                                   grepl("d11", names(allelicCDPaucList))  |
                                                   grepl("d20", names(allelicCDPaucList))]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, titleText="Differentiated cells: scaled Spearman", shortLeg=FALSE)
    
    dev.off()
    
  }
  
  
  ##### For figures #####
  
  # 20200803
  if (patski100 ) {
  
  ### 100 Patski
  
  if (aThresh == "100") {
    PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.AUCs.Patski.4fig.pdf", sep=""))
    pdf(PDFfilename, height=8, width=8)
    allelicCDPaucList.subset = allelicCDPaucList[grepl("^Patski", names(allelicCDPaucList)) & grepl(".100.", names(allelicCDPaucList) ) ]
    # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
    allelicCDPaucCols.subset <- "black" # rainbow(length(allelicCDPaucList.subset))
    print(length(allelicCDPaucCols.subset))
    plotAUCs(allelicCDPaucList.subset, cols2Use=allelicCDPaucCols.subset, titleText="Patski cells: countThresh100 scaled Spearman", shortLeg=TRUE)
    dev.off()
  }
    
  } else {
    
    ### TsixStop
    if (aThresh == "50") {
      PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.AUCs.TsixStop.4fig.pdf", sep=""))
      pdf(PDFfilename, height=8, width=8)
      
      allelicCDPaucList.subset = allelicCDPaucList[grepl("^Tsix", names(allelicCDPaucList))]
      # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
      allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
      print(length(allelicCDPaucCols.subset))
      plotAUCs(allelicCDPaucList.subset, titleText="TsixStop cells: countThresh50 scaled Spearman", shortLeg=FALSE)
      dev.off()
      
      ### F121
      PDFfilename = file.path(aucSubdir, paste(projectName, ".", aCompType, ".Xi-calling.AUCs.F121.4fig.pdf", sep=""))
      pdf(PDFfilename, height=8, width=8)
      
      allelicCDPaucList.subset = allelicCDPaucList[grepl("^F121", names(allelicCDPaucList)) & 
                                                     ! grepl("rep", names(allelicCDPaucList))]
      # allelicCDPaucCols.subset = allelicCDPaucCols[ names(allelicCDPaucList.subset) ]
      allelicCDPaucCols.subset <- rainbow(length(allelicCDPaucList.subset))
      print(length(allelicCDPaucCols.subset))
      plotAUCs(allelicCDPaucList.subset, titleText="F121 cells: countThresh50 scaled Spearman", shortLeg=FALSE)
      dev.off()
    }  
    
  }
  
} # aCompType










##################################################
# 20200406
# Joint density plots based on NPC threshold

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(allelicCDPcompList) ) {
  
  compMaindir = file.path(supermainOutputDir, paste( "countThresh", aThresh, ".chrXvsChr1.differentiatedThreshold", 
                                                      ifelse(useInterphaseOnly, ".interphaseOnly", ""), 
                                                      sep="" ))
  if ( ! file.exists( compMaindir ) ) {
    if ( ! dir.create( compMaindir, showWarnings=FALSE, recursive=TRUE )  ) compMaindir = supermainOutputDir
  }
  
  for( aBins2Use in names(allelicCDPcompList[[aThresh]][["chrX"]]) ) {
    # for( aCompType in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]]) ) {
    
    for( aCompType in c( "Spearman", "LMD" )) {
    
      compSubdir = file.path(compMaindir, aCompType)
      if ( ! file.exists( compSubdir ) ) {
        if ( ! dir.create( compSubdir, showWarnings=FALSE, recursive=TRUE )  ) compSubdir = compMaindir
      }
      
      for( aLib in names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]) ) {
        
        # for( aCDPnormType in rev(names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        for( aCDPnormType in c( "scaled" ) ) {  
          
          print(aLib)
          whatLib = sub("_$", "", 
                        ifelse(!grepl("_", aLib), aLib, 
                               sub("[^_]+$", "", 
                                   sub("_rep.$", "", aLib) )))
          tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
          
          if ( aCompType == "LMD" ) {
            maxDiff = 1
            for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
              thisMax = max(abs(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]]), na.rm=TRUE)
              maxDiff = ifelse(thisMax > maxDiff, thisMax, maxDiff)
            }
            Xlims = c(-maxDiff, maxDiff)
          } else  if ( aCompType == "Spearman" ) {
            maxDiff = 1
            for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
              thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
              maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
            }
            Xlims = c(maxDiff, 1)
          }
            
          tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
          aChr = "chrX"
          anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
          if ( aCompType == "LMD" ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else if ( grepl("dist", aCompType) ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else {
            XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
          }
          # head(anAllelicCDPcomp)
          
          tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
          tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
          print(c(tenFPRcellNum, tenFPRrate))
          
          
          ################
          # Joint density plots
          
          NAflag = TRUE
          dpList = list()
          for ( aChr in c("chr1", "chrX") ) {
            anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
            if (all(is.na(anAllelicCDPcomp)) || length(na.omit(anAllelicCDPcomp))<2) {
              NAflag = FALSE
              next 
            }
            dpList[[aChr]] = density(anAllelicCDPcomp, na.rm=TRUE)
          }
          
          if (NAflag) {
            
            maxDensity = max(c(dpList[[1]]$y, dpList[[2]]$y), na.rm=TRUE)
            distRange = c(min(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE),
                          max(c(dpList[[1]]$x, dpList[[2]]$x), na.rm=TRUE))
            
            titleText = paste(aLib, " - CDP comparison distributions\n", 
                              aCompType, " - chr1 vs chrX", sep="")
            nameText = paste(aLib, aCompType, "densityPlot", sep=".")
            PDFfilename = file.path(compSubdir, paste(projectName, ".refVsAlt.withThresholds.", nameText, 
                                                      ifelse(useInterphaseOnly, ".interphaseOnly", ""),
                                                      ".pdf", sep=""))
            pdf(PDFfilename, height=7, width=7)
            defparmar = par()$mar
            par(mar=par()$mar+c(0,0,2,0))
            plot(dpList[["chrX"]], lwd=1, lty="solid", xlim=Xlims, ylim=c(0, maxDensity),
                 main=titleText, col="red")
            # abline(v = median(dpList[["chrX"]]$x, na.rm=TRUE), lty="dashed")
            # text(median(dpList[["chrX"]]$x, na.rm=TRUE), maxDensity*0.05, pos=2,
            #      paste("med: ",  round(median(dpList[["chrX"]]$x, na.rm=TRUE), 2), sep=""))
            lines(dpList[["chr1"]]$x, dpList[["chr1"]]$y, lwd=1, lty="solid", col="blue")
            # abline(v = median(dpList[["chr1"]]$x, na.rm=TRUE), lty="dashed", col="blue")
            # text(median(dpList[["chr1"]]$x, na.rm=TRUE), maxDensity*0.15, pos=2,
            #      paste("med: ",  round(median(dpList[["chr1"]]$x, na.rm=TRUE), 2), sep=""),
            #      col="blue")
            if ( aCompType == "LMD") abline(v=-tenFPRrateThresh, lty="dashed")
            abline(v=tenFPRrateThresh, lty="dashed")
            if ( aCompType == "LMD") {
              legend("topright", 
                     legend = paste(tenFPRrate, "% ", tenFPRcellNum, "\n",
                                    "@10%FPR threshold (",  round(tenFPRrateThresh, 2), ")", sep=""),
                     bty='n' )
            } else {
              text(tenFPRrateThresh, maxDensity*0.95, pos=2,
                   paste(tenFPRrate, "% ", tenFPRcellNum, "\n",
                         "@10%FPR threshold (",  round(tenFPRrateThresh, 2), ")", sep=""),
                   col="black", cex=0.8)
            }
            legend("topleft", legend=c(paste( "chrX n = ", dpList[["chrX"]]$n, sep=""),
                                       paste( "chr1 n = ", dpList[["chr1"]]$n, sep="")),
                   fill=c("red", "blue"),
                   border=NA, bty='n')
            par(mar=defparmar)
            dev.off()
            
          } # NAflag
        } # aCDPnormType
      } # aLib
    } # aCompType
  } # aBins2Use
} # aThresh










##################################################
# 20200318 
# Decay curves for active vs inactive cells

decayCurveMaindir = file.path(mainOutputDir, "decayCurves")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

######################
# 20200315
# Plot all ref and alt curves overlaid 
# Average of spearman correlations of CDPs rather than correlation of average CDPs.

lwd2Use = 0.5
alpha2Use=0.01

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    for( aCompType in c( "Spearman", "LMD")) {
      
      for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
        
        # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        for( aCDPnormType in c( "scaled" ) ) {  
          
          decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI", sep="."), aCompType)
          if ( ! file.exists( decayOutdir ) ) {
            if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
          }
          
          tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
          aChr = "chrX"
          anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
          if ( aCompType == "LMD" ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else if ( grepl("dist", aCompType) ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else {
            XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
          }
          # head(anAllelicCDPcomp)
          print(sum(XCIboolIdx))
          
          if (aCDPnormType == "prop") {
            ylims.prop = range(c( unlist(dataMatrices[[aThresh]][["chr1"]][[bins2Use]][[aCDPnormType]]),
                                  unlist(dataMatrices[[aThresh]][["chrX"]][[bins2Use]][[aCDPnormType]]) ),na.rm=TRUE)
          }
          
          for (whichChr in c("chrX", "chr1")) {
            
            for  (plotXCI in c(TRUE, FALSE)) {
              
              meanCDP = list()
              
              titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
              nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
              
              PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayCurve.overlaid.v2.", 
                                                         nameText, ".refValt.pdf", sep=""))
              pdf(PDFfilename, height=4, width=6)
              
              nuPlotflag=TRUE
              for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
                print(alleliclib)
                
                if(plotXCI) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,XCIboolIdx]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                }
                
                nCells = ncol(matrix2Use)
                print(nCells)
                if (is.null(nCells) || nCells==0) {
                  if(grepl("alt", alleliclib) && nuPlotflag) {
                    plot(1, type="n", axes=F, xlab="", ylab="",
                         main=paste(titleText, " (n=", nCells, ")", sep=""))
                  }
                  meanCDP[[alleliclib]] = NaN
                  next
                }
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                
                for (iCell in 1:nCells) {
                  y = matrix2Use[,iCell]
                  if (iCell==1 && nuPlotflag){
                    if (aCDPnormType == "scaled") {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " (n=", nCells, ")", sep=""),
                          # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          xlab="", ylab="",
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                meanCDP[[alleliclib]] = y
                nuPlotflag = FALSE
              } # alleliclib
              
              alleliclib = sub("\\.alt$", ".ref", alleliclib)
              if(plotXCI) {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,XCIboolIdx]
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
              }
              
              nCells = ncol(matrix2Use)
              if (! is.null(nCells) && nCells>0) {
                
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
                
                interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
                
                if(plotXCI) {
                  matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,XCIboolIdx]
                  matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,XCIboolIdx]
                } else {
                  matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,!XCIboolIdx]
                  matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,!XCIboolIdx]
                }
                
                CDPspearmans = c()
                for(iCell in 1:nCells) {
                  CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                                    method='spearman', use='pairwise.complete.obs')
                  CDPspearmans = c(CDPspearmans, CDPspearman)
                }
                
                meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
                medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
                
                legend("bottom", legend=c(
                  paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
                  paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
                  bty='n')
                
              }
                
              dev.off() 

            } # plotXCI
          } # whichChr
        } # aCDPnormType
      } # aLib
    } # aCompType
  } # aBins2Use
} # aThresh










##################################################
# 20200406 
# Decay curves for active vs inactive cells using differentiated threshold

decayCurveMaindir = file.path(mainOutputDir, "decayCurves")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

######################
# 20200315
# Plot all ref and alt curves overlaid 
# Average of spearman correlations of CDPs rather than correlation of average CDPs.

lwd2Use = 0.5
alpha2Use=0.01

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    for( aCompType in c( "Spearman", "LMD")) {
      
      for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
        
        # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        for( aCDPnormType in c( "scaled" ) ) { 

          decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold", sep="."), aCompType)
          if ( ! file.exists( decayOutdir ) ) {
            if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
          }

          print(aLib)
          whatLib = sub("_$", "", 
                        ifelse(!grepl("_", aLib), aLib, 
                               sub("[^_]+$", "", 
                                   sub("_rep.$", "", aLib) )))
          tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
          
          # maxDiff = 1
          # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
          #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
          #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
          # }
          
          tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
          aChr = "chrX"
          anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
          if ( aCompType == "LMD" ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else if ( grepl("dist", aCompType) ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else {
            XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
          }
          # head(anAllelicCDPcomp)
          
          # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
          # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
          # print(c(tenFPRcellNum, tenFPRrate))
          
          if (aCDPnormType == "prop") {
            ylims.prop = range(c( unlist(dataMatrices[[aThresh]][["chr1"]][[bins2Use]][[aCDPnormType]]),
                                  unlist(dataMatrices[[aThresh]][["chrX"]][[bins2Use]][[aCDPnormType]]) ),na.rm=TRUE)
          }
          
          for (whichChr in c("chrX", "chr1")) {
            
            for  (plotXCI in c(TRUE, FALSE)) {
              
              meanCDP = list()
              
              titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
              nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
              
              PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayCurve.overlaid.v2.", 
                                                         nameText, ".refValt.pdf", sep=""))
              pdf(PDFfilename, height=4, width=6)
              
              nuPlotflag=TRUE
              for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
                print(alleliclib)
                
                if(plotXCI) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,XCIboolIdx]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                }
                
                nCells = ncol(matrix2Use)
                print(nCells)
                if (is.null(nCells) || nCells==0) {
                  if(grepl("alt", alleliclib) && nuPlotflag) {
                    plot(1, type="n", axes=F, xlab="", ylab="",
                         main=paste(titleText, " (n=", nCells, ")", sep=""))
                  }
                  meanCDP[[alleliclib]] = NaN
                  next
                }
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                
                for (iCell in 1:nCells) {
                  y = matrix2Use[,iCell]
                  if (iCell==1 && nuPlotflag){
                    if (aCDPnormType == "scaled") {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " (n=", nCells, ")", sep=""),
                          # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          xlab="", ylab="",
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                meanCDP[[alleliclib]] = y
                nuPlotflag = FALSE
                
              } # alleliclib
              
              alleliclib = sub("\\.alt$", ".ref", alleliclib)
              if(plotXCI) {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,XCIboolIdx]
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
              }
              
              nCells = ncol(matrix2Use)
              if (!is.null(nCells) && nCells>0) {
 
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
                
                interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
                
                if(plotXCI) {
                  matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,XCIboolIdx]
                  matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,XCIboolIdx]
                } else {
                  matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,!XCIboolIdx]
                  matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,!XCIboolIdx]
                }
                
                CDPspearmans = c()
                for(iCell in 1:nCells) {
                  CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                                    method='spearman', use='pairwise.complete.obs')
                  CDPspearmans = c(CDPspearmans, CDPspearman)
                }
                
                meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
                medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
                
                legend("bottom", legend=c(
                  paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
                  paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
                  bty='n')
              }
              
             dev.off()
              
            } # plotXCI
          } # whichChr
        } # aCDPnormType
      } # aLib
    } # aCompType
  } # aBins2Use
} # aThresh










##################################################
# 20200422 
# Decay curves for active vs inactive cells using differentiated threshold
# Using LMD to indicate allelic direction of XCI

# *** NOTE: CDP curves for both chr1 and chrX plotted in same chrX directory ***
aChr = "chrX"

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", aChr, sep="")), sep=""))
# if ( ! file.exists( subsupermainOutputDir ) ) {
#   if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
# }

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
# if ( ! file.exists( mainOutputDir ) ) {
#   if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
# }

decayCurveMaindir = file.path(mainOutputDir, "decayCurves")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

decayCurveTableMaindir = file.path(mainOutputDir, "decayCurvesTables")
if ( ! file.exists( decayCurveTableMaindir ) ) {
  if ( ! dir.create( decayCurveTableMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveTableMaindir = mainOutputDir
}


######################
# 20200315
# Plot all ref and alt curves overlaid 
# Average of spearman correlations of CDPs rather than correlation of average CDPs.
# 
# 2020425
# Tabulate

lwd2Use = 0.5
alpha2Use=0.01

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    for( aCompType in c( "Spearman", "LMD")) {
      
      # 2020425
      propXCItable = list()
      propXCItable[["scaled"]] = c()
      propXCItable[["prop"]] = c()
      
      for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
        
        # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        # for( aCDPnormType in c( "scaled" ) ) {  
        aCDPnormType = "scaled"
        
        decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele", sep="."), aCompType)
        if ( ! file.exists( decayOutdir ) ) {
          if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
        }

        # aLib = "F121_NPC"
        print(aLib)
        whatLib = sub("_$", "", 
                      ifelse(!grepl("_", aLib), aLib, 
                             sub("[^_]+$", "", 
                                 sub("_rep.$", "", aLib) )))
        tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
        
        # maxDiff = 1
        # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
        #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
        #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
        # }
        
        tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
        anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
        if ( aCompType == "LMD" ) {
          XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
        } else if ( grepl("dist", aCompType) ) {
          XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
        } else {
          XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
        }
        # head(anAllelicCDPcomp)
        
        LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
        refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
        altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
        
        # 2020425
        propXCIrow = c( aLib, aCompType, aCDPnormType, aBins2Use, aThresh, 
                        length(refXCIboolIdx), 
                        sum(!XCIboolIdx),
                        sum(XCIboolIdx),
                        sum(refXCIboolIdx),
                        sum(altXCIboolIdx),
                        sum(XCIboolIdx) - sum(refXCIboolIdx | altXCIboolIdx),
                        sum(!XCIboolIdx) / length(XCIboolIdx),
                        sum(XCIboolIdx) / length(refXCIboolIdx),
                        sum(refXCIboolIdx) / length(refXCIboolIdx),
                        sum(altXCIboolIdx) / length(altXCIboolIdx),
                        ( sum(XCIboolIdx) - sum(refXCIboolIdx | altXCIboolIdx) ) / length(XCIboolIdx) )
        
        propXCItable[[aCDPnormType]] = rbind( propXCItable[[aCDPnormType]],
                                              propXCIrow)
        
        # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
        # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
        # print(c(tenFPRcellNum, tenFPRrate))
        
        if (aCDPnormType == "prop") {
          ylims.prop = range(c( unlist(dataMatrices[[aThresh]][["chr1"]][[bins2Use]][[aCDPnormType]]),
                                unlist(dataMatrices[[aThresh]][["chrX"]][[bins2Use]][[aCDPnormType]]) ),na.rm=TRUE)
        }
        
        # *** NOTE: CDP curves for both chr1 and chrX plotted in same chrX directory ***
        for (whichChr in c("chrX", "chr1")) {
          
          for (plotXCI in c(TRUE, FALSE)) {
            
            for  (plotRefAlt in c(TRUE, FALSE)) {
              
              meanCDP = list()
              
              titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
              nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
              
              if (plotXCI) {
                titleText = paste(titleText, ifelse(plotRefAlt, "refXCI", "altXCI"), sep=" ")
                nameText = paste(nameText, ifelse(plotRefAlt, "refXCI", "altXCI"), sep=".")
              } else {
                if(!plotRefAlt) next
              }
              
              PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayCurve.overlaid.v2.", 
                                                         nameText, ".refValt.pdf", sep=""))
              pdf(PDFfilename, height=4, width=6)
              
              nuPlotflag=TRUE
              for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
                print(alleliclib)
                
                if(plotXCI) {
                  if(plotRefAlt) {
                    matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                  } else {
                    matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                  }
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                }
                
                nCells = ncol(matrix2Use)
                print(nCells)
                if (is.null(nCells) || nCells==0) {
                  if(grepl("alt", alleliclib) && nuPlotflag) {
                    plot(1, type="n", axes=F, xlab="", ylab="",
                         main=paste(titleText, " (n=", nCells, ")", sep=""))
                  }
                  meanCDP[[alleliclib]] = NaN
                  next
                }
                
                # alaph2Use =  1/nCells
                col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                
                for (iCell in 1:nCells) {
                  y = matrix2Use[,iCell]
                  if (iCell==1 && nuPlotflag){
                    if (aCDPnormType == "scaled") {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( y, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " (n=", nCells, ")", sep=""),
                          # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          xlab="", ylab="",
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                  } else {
                    lines(y, col=col2Use, lwd=lwd2Use)
                  }
                }
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                meanCDP[[alleliclib]] = y
                nuPlotflag = FALSE
              } # alleliclib
              
              alleliclib = sub("\\.alt$", ".ref", alleliclib)
              
              if(plotXCI) {
                if(plotRefAlt) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                }
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
              }
              
              nCells = ncol(matrix2Use)
              if (!is.null(nCells) && nCells>0) {
                y = rowMeans(matrix2Use, na.rm = TRUE)
                lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
                interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
                
                if(plotXCI) {
                  if(plotRefAlt) {
                    matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,refXCIboolIdx]
                    matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,refXCIboolIdx]
                  } else {
                    matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,altXCIboolIdx]
                    matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,altXCIboolIdx]
                  }
                } else {
                  matrix2Use.ref = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "ref", sep=".")]][,!XCIboolIdx]
                  matrix2Use.alt = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[paste(aLib, "alt", sep=".")]][,!XCIboolIdx]
                  
                }
                
                CDPspearmans = c()
                for(iCell in 1:nCells) {
                  CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                                    method='spearman', use='pairwise.complete.obs')
                  CDPspearmans = c(CDPspearmans, CDPspearman)
                }
                meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
                medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
                
                legend("bottom", legend=c(
                  paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
                  paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
                  bty='n')
              }
              
              dev.off()
              
            } # plotRefAlt
          } # plotXCI
        } # whichChr
        
        # } # aCDPnormType
        
      } # aLib
      
      # 2020425
      colnames( propXCItable[[aCDPnormType]] ) = c("dataset",
                                                   "CDP_comparison_method",
                                                   "normalization_method",
                                                   "bin_size",
                                                   "contact_count_threshold",
                                                   "total_cells",
                                                   "non_XCI_cells",
                                                   "XCI_cells",
                                                   "ref_XCI_cells",
                                                   "alt_XCI_cells",
                                                   "unassigned_XCI_cells",
                                                   "non_XCI_cell_proportion",
                                                   "XCI_cell_proportion",
                                                   "ref_XCI_cell_proportion",
                                                   "alt_XCI_cell_proportion",
                                                   "unassigned_XCI_cell_proportion" )
      
      decayTablesOutdir = file.path(decayCurveTableMaindir, paste( aCompType, aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele", sep="."))
      if ( ! file.exists( decayTablesOutdir ) ) {
        if ( ! dir.create( decayTablesOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayTablesOutdir = decayCurveTableMaindir
      }
 
      for( aCDPnormType in names(propXCItable) ) {  
        nameText = paste("XCIcallingResults", aCompType, aCDPnormType, sep=".")
        
        tableFileName = file.path(decayTablesOutdir, paste(projectName, ".", nameText, ".tsv", sep=""))
        write.table(x=propXCItable[[aCDPnormType]] ,
                    file=tableFileName, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
      }

    } # aCompType
  } # aBins2Use
} # aThresh










##################################################
# 20210316
# Decay MAPS for active vs inactive cells using differentiated threshold
# Using LMD to indicate allelic direction of XCI
#
# BUT NOT BINARIZED VERSIONS.
# AND ONLY FOR BOTH ALLELES ALONGSIDE ONE ANOTHER.
# AND ONLY FOR LMD (aka LMD) XCI status.
# AND ONLY FOR SCALED CDPs.

# *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
aChr = "chrX"

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", aChr, sep="")), sep=""))
# if ( ! file.exists( subsupermainOutputDir ) ) {
#   if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
# }

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
# if ( ! file.exists( mainOutputDir ) ) {
#   if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
# }

decayCurveMaindir = file.path(mainOutputDir, "decayMaps.V2")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    # for( aCompType in c( "LMD")) {
    aCompType = "LMD"
    
    for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
      # aLib = "TsixStop_NPC"
      
      # # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
      # for( aCDPnormType in c( "scaled" ) ) {  
      aCDPnormType = "scaled"
      
      decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele", sep="."), aCompType)
      if ( ! file.exists( decayOutdir ) ) {
        if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
      }
      
      print(aLib)
      whatLib = sub("_$", "", 
                    ifelse(!grepl("_", aLib), aLib, 
                           sub("[^_]+$", "", 
                               sub("_rep.$", "", aLib) )))
      tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
      
      # maxDiff = 1
      # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
      #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
      #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
      # }
      
      tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
      anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
      if ( aCompType == "LMD" ) {
        XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
      } else if ( grepl("dist", aCompType) ) {
        XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
      } else {
        XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
      }
      # head(anAllelicCDPcomp)
      
      LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
      refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
      altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
      
      # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
      # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
      # print(c(tenFPRcellNum, tenFPRrate))
    
      # *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
      for (whichChr in c("chrX", "chr1")) {
        
        for (plotXCI in c(FALSE, TRUE)) {
          
          for  (plotRefAlt in c(TRUE, FALSE)) {
            # whichChr="chrX"; plotXCI=TRUE; plotRefAlt=TRUE;
            
            if (!plotXCI && plotRefAlt==FALSE) next;
            
            titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
            nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
            
            if (plotXCI) {
              titleText = paste(titleText, ifelse(plotRefAlt, "refXCI", "altXCI"), sep=" ")
              nameText = paste(nameText, ifelse(plotRefAlt, "refXCI", "altXCI"), sep=".")
            } else {
              if(!plotRefAlt) next
            }
            
            matrix2Plot.combined = c()
            
            for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
              # alleliclib = "TsixStop_NPC.ref"
              print(alleliclib)
              
              if(plotXCI) {
                if(plotRefAlt) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][refXCIboolIdx][ timepointThenCovOrd[[alleliclib]][refXCIboolIdx] ]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][altXCIboolIdx][ timepointThenCovOrd[[alleliclib]][altXCIboolIdx] ]
                }
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                # covOrd2Use = covOrd[[alleliclib]][!XCIboolIdx][ timepointThenCovOrd[[alleliclib]][!XCIboolIdx] ]
              }
              
              if ( is.null(matrix2Use) || length(matrix2Use) == 0 || length(matrix2Use) == 14 ) next
              # length(matrix2Use) == 14 ==> matrix has one row
              # ncol(matrix2Use) < 25  heatmaps look bad
              
              # Combine CDPs from both alleles
              matrix2Plot.combined = rbind( matrix2Plot.combined, matrix2Use )
              # if ( grepl( ".ref", alleliclib)) {
              #   covOrd2Use.combined = covOrd2Use
              # }
              
              # Only plot combined plots 
              # matrix2Plot = matrix2Use
              # # ##################################################
              # # # Sort cells by timepoint and coverage
              # # 
              # # # head(matrixCsums.PF[[whichChr]][[alleliclib]][ covOrd[[alleliclib]] ])
              # # # tail(matrixCsums.PF[[whichChr]][[alleliclib]][ covOrd[[alleliclib]] ])
              # # 
              # # # matrix2Plot = matrix2Use[ ,covOrd[[alleliclib]][timepointThenCovOrd[[alleliclib]]] ]
              # # matrix2Plot = matrix2Use[ ,covOrd2Use ]
              # # cellLabels = rep("", ncol(matrix2Plot))
              # # cellLabels[1] = colnames(matrix2Plot)[1]
              # # cellLabels[length(cellLabels)] = colnames(matrix2Plot)[length(cellLabels)]
              # 
              # rownames(matrix2Plot) = paste(alleliclib, rownames(matrix2Plot), sep=".")
              # 
              # # # TPcolbar = matrixColTPcols[[alleliclib]][covOrd[[alleliclib]]][timepointThenCovOrd[[alleliclib]]]
              # # # covcolbar = covCols[[alleliclib]][covOrd[[alleliclib]]][timepointThenCovOrd[[alleliclib]]]
              # # TPcolbar = matrixColTPcols[[alleliclib]][covOrd2Use]
              # # covcolbar = covCols[[alleliclib]][covOrd2Use]
              # # rColBars = t(cbind(TPcolbar, covcolbar))
              # # rownames(rColBars)=c("sample", "coverage")
              # TPcolbar = matrixColTPcols[[alleliclib]]
              # rColBars = t(cbind(TPcolbar))
              # rownames(rColBars)=c("sample")
              # 
              # titleTextb = paste(titleText, " ", alleliclib, " (n = ", dim(matrix2Plot)[2], ")", sep="")
              # nameTextb = paste(nameText, alleliclib, sep=".")
              # 
              # print(dim(matrix2Plot))
              # minVal = min(matrix2Plot, na.rm=TRUE)
              # maxVal = max(matrix2Plot, na.rm=TRUE)
              # # cat(aCDPnormType, ": ", minVal, " --> ", maxVal, "\n", sep="" )
              # 
              # if(aCDPnormType == "scaled") {
              #   cols2Use = heatCols2
              #   # minVal=-4; maxVal=4
              #   breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
              #   # heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
              # } else {
              #   cols2Use = heatColsV
              #   # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
              #   # minVal = -3; maxVal = -1.5 
              #   breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
              #   # heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
              # }
              # 
              # # if ( ! file.exists( heatOutdir ) ) {
              # #   if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
              # # }
              # 
              # # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
              # PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayMap.", 
              #                                            nameTextb, ".refValt.pdf", sep=""))
              # pdf(PDFfilename, height=6, width=12)
              # 
              # # heatmap.2(t(matrix2Plot), col=heatCols2,
              # #           Rowv=TRUE, Colv=FALSE,
              # #           dendrogram="none",
              # #           scale="none", key=TRUE, density.info="none",
              # #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
              # #           RowSideColors=kClustColbar,
              # #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
              # 
              # # https://www.biostars.org/p/18211/
              # # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
              # heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
              #              Rowv=FALSE, Colv=FALSE,
              #              dendrogram="none",
              #              scale="none", key=TRUE, density.info="none",
              #              trace="none", cexCol=0.95, labCol=binDistDesc, labRow=cellLabels, margins=c(5,0),
              #              RowSideColors=rColBars, RowSideColorsSize=4,
              #              symm=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
              #              symkey=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
              #              symbreaks=ifelse(aCDPnormType=="scaled", TRUE, FALSE),
              #              keysize=0.2, 
              #              KeyValueName=titleText,
              #              #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
              #              # key.par=list(mar=c(3,0,3,0)),
              #              # lmat -- added 2 lattice sections (5 and 6) for padding
              #              lmat=rbind(c(0, 0, 5, 4), c(0,1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
              # 
              # dev.off()
              # # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6)
              
            } # alleliclib
            
            if (is.null(matrix2Plot.combined)) next
            
            matrix2Plot = matrix2Plot.combined
            print(dim(matrix2Plot))
            
            # TPcolbar = matrixColTPcols[[paste(aLib, "ref", sep=".")]][covOrd2Use.combined]
            TPcolbar = matrixColTPcols[[paste(aLib, "ref", sep=".")]]
            # covcolbar = covCols[[lib]][covOrd[[aLib]]][timepointThenCovOrd[[aLib]]]
            # rColBars = t(cbind(TPcolbar, covcolbar))
            # rownames(rColBars)=c("sample", "coverage")
            rColBars = t(cbind(TPcolbar))
            rownames(rColBars)=c("sample")
            
            alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
            rownames(alleleBars)=c("allele")
            colnames(alleleBars)=rownames(matrix2Plot)
            
            titleTextb = paste(titleText,  " bothAlleles  (n = ", dim(matrix2Plot)[2], ")", sep="")
            nameTextb = paste(nameText, "bothAlleles", sep=".")
            
            minVal = min(matrix2Plot, na.rm=TRUE)
            maxVal = max(matrix2Plot, na.rm=TRUE)
            # cat(aCDPnormType, ": ", minVal, " --> ", maxVal, "\n", sep="" )
            
            if(aCDPnormType == "scaled") {
              cols2Use = heatCols2
              # minVal=-4; maxVal=4
              breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
              # heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
            } else {
              cols2Use = heatColsV
              # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
              # minVal = -3; maxVal = -1.5 
              breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
              # heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
            }
            
            # # 20210317
            # # Tried to make height proportional to number of rows.
            # # But doesn't really work b/c it also scales the color scale bar.
            # map_height = 4 + (4 * (dim(matrix2Plot)[2] / 1000))
            map_height = 6
            
            # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
            PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayMap.", 
                                                       nameTextb, ".combined.pdf", sep=""))
            
            pdf(PDFfilename, height=map_height, width=12)
            
            # heatmap.2(t(matrix2Plot), col=heatCols2,
            #           Rowv=TRUE, Colv=FALSE,
            #           dendrogram="none",
            #           scale="none", key=TRUE, density.info="none",
            #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
            #           RowSideColors=kClustColbar,
            #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
            
            # https://www.biostars.org/p/18211/
            # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
            heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                         Rowv=FALSE, Colv=FALSE,
                         dendrogram="none",
                         scale="none", key=TRUE, density.info="none",
                         trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), labRow=cellLabels, margins=c(5,0),
                         RowSideColors=rColBars, RowSideColorsSize=4,
                         ColSideColors=t(alleleBars), ColSideColorsSize=4,
                         symm=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                         symkey=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                         symbreaks=ifelse(aCDPnormType=="scaled", TRUE, FALSE),
                         keysize=0.2, 
                         KeyValueName=titleTextb,
                         #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                         # key.par=list(mar=c(3,0,3,0)),
                         # lmat -- added 2 lattice sections (5 and 6) for padding
                         # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                         lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(1, 0.2, 5), lwid=c(0.2, 0.2, 5, 1))
            dev.off()
            # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
            
          } # plotRefAlt
        } # plotXCI
      } # whichChr
      
      # } # aCDPnormType
      
    } # aLib
    
    # } # aCompType
    
  } # aBins2Use
} # aThresh










##################################################
# 20210317
# *** ALL CELL STATES IN ONE PLOT ***
#
# Decay MAPS for active vs inactive cells using differentiated threshold
# Using LMD to indicate allelic direction of XCI
#
# BUT NOT BINARIZED VERSIONS.
# AND ONLY FOR BOTH ALLELES ALONGSIDE ONE ANOTHER.
# AND ONLY FOR LMD (aka LMD) XCI status.
# AND ONLY FOR SCALED CDPs.

# *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
aChr = "chrX"

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", aChr, sep="")), sep=""))
# if ( ! file.exists( subsupermainOutputDir ) ) {
#   if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
# }

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
# if ( ! file.exists( mainOutputDir ) ) {
#   if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
# }

decayCurveMaindir = file.path(mainOutputDir, "decayMaps")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    # for( aCompType in c( "LMD" )) {
    aCompType = "LMD"
    
    for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
      # aLib = "TsixStop_NPC"
      
      # # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
      # for( aCDPnormType in c( "scaled" ) ) {  
      aCDPnormType = "scaled" 
      
      decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele.allXCIstates.V2", sep="."), aCompType)
      if ( ! file.exists( decayOutdir ) ) {
        if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
      }
      
      print(aLib)
      whatLib = sub("_$", "", 
                    ifelse(!grepl("_", aLib), aLib, 
                           sub("[^_]+$", "", 
                               sub("_rep.$", "", aLib) )))
      tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
      
      # maxDiff = 1
      # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
      #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
      #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
      # }
      
      tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
      anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
      if ( aCompType == "LMD" ) {
        XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
      } else if ( grepl("dist", aCompType) ) {
        XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
      } else {
        XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
      }
      # head(anAllelicCDPcomp)
      
      LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
      refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
      altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
      
      # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
      # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
      # print(c(tenFPRcellNum, tenFPRrate))
      
      # *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
      for (whichChr in c("chrX", "chr1")) {
        # whichChr = "chrX"
        
        matrix2Plot.all = c()
        TPcolbar = c()
        
        for (plotXCI in c(FALSE, TRUE)) {
          
          for  (plotRefAlt in c(TRUE, FALSE)) {
            # whichChr="chrX"; plotXCI=TRUE; plotRefAlt=TRUE;
            
            if (!plotXCI && plotRefAlt==FALSE) next;
            
            matrix2Plot.combined = c()
            
            for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
              # alleliclib = "TsixStop_NPC.ref"
              print(alleliclib)
              
              if(plotXCI) {
                if(plotRefAlt) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][refXCIboolIdx][ timepointThenCovOrd[[alleliclib]][refXCIboolIdx] ]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][altXCIboolIdx][ timepointThenCovOrd[[alleliclib]][altXCIboolIdx] ]
                }
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                # covOrd2Use = covOrd[[alleliclib]][!XCIboolIdx][ timepointThenCovOrd[[alleliclib]][!XCIboolIdx] ]
              }
              
              if ( is.null(matrix2Use) || length(matrix2Use) == 0 || length(matrix2Use) == 14 ) next
              # length(matrix2Use) == 14 ==> matrix has one row
              # ncol(matrix2Use) < 25  heatmaps look bad
              
              # Combine CDPs from both alleles
              matrix2Plot.combined = rbind( matrix2Plot.combined, matrix2Use )
              # if ( grepl( ".ref", alleliclib)) {
              #   covOrd2Use.combined = covOrd2Use
              # }
              
            } # alleliclib
            
            if (is.null(matrix2Plot.combined)) next
            
            matrix2Plot.all = cbind(matrix2Plot.all, matrix2Plot.combined)
            
            # titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
            # nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
            
            TPcol = "lightgreen" # active
            if (plotXCI && plotRefAlt) TPcol = "darkblue"
            if (plotXCI && !plotRefAlt) TPcol = "red4"
            
            TPcolbar = c( TPcolbar, rep( TPcol, dim(matrix2Plot.combined)[2] ))
            # covcolbar = covCols[[lib]][covOrd[[aLib]]][timepointThenCovOrd[[aLib]]]
            # rColBars = t(cbind(TPcolbar, covcolbar))
            # rownames(rColBars)=c("sample", "coverage")
            
          } # plotRefAlt
        } # plotXCI
        
        rColBars = t(cbind(TPcolbar))
        rownames(rColBars)=c("XCI state")
        
        matrix2Plot = matrix2Plot.all
        
        alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
        rownames(alleleBars)=c("allele")
        colnames(alleleBars)=rownames(matrix2Plot)
        
        titleText = paste(whichChr, aLib, aCDPnormType, sep=" ")
        nameText = paste(whichChr, aLib, aCDPnormType, sep=".")
        
        titleTextb = paste(titleText,  " bothAlleles  (n = ", dim(matrix2Plot)[2], ")", sep="")
        nameTextb = paste(nameText, "bothAlleles", sep=".")
        
        minVal = min(matrix2Plot, na.rm=TRUE)
        maxVal = max(matrix2Plot, na.rm=TRUE)
        # cat(aCDPnormType, ": ", minVal, " --> ", maxVal, "\n", sep="" )
        
        if(aCDPnormType == "scaled") {
          cols2Use = heatCols2
          # minVal=-4; maxVal=4
          breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
          # heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
        } else {
          cols2Use = heatColsV
          # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
          # minVal = -3; maxVal = -1.5 
          breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
          # heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
        }
        
        # # 20210317
        # # Tried to make height proportional to number of rows.
        # # But doesn't really work b/c it also scales the color scale bar.
        # map_height = 4 + (4 * (dim(matrix2Plot)[2] / 1000))
        map_height = 8
        
        # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
        PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayMap.", 
                                                   nameTextb, ".combined.pdf", sep=""))
        
        pdf(PDFfilename, height=map_height, width=12)
        
        # heatmap.2(t(matrix2Plot), col=heatCols2,
        #           Rowv=TRUE, Colv=FALSE,
        #           dendrogram="none",
        #           scale="none", key=TRUE, density.info="none",
        #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
        #           RowSideColors=kClustColbar,
        #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
        
        # https://www.biostars.org/p/18211/
        # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
        heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                     Rowv=FALSE, Colv=FALSE,
                     dendrogram="none",
                     scale="none", key=TRUE, density.info="none",
                     trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), margins=c(5,0), labRow="", # labRow=cellLabels, 
                     RowSideColors=rColBars, RowSideColorsSize=4,
                     ColSideColors=t(alleleBars), ColSideColorsSize=4,
                     symm=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                     symkey=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                     symbreaks=ifelse(aCDPnormType=="scaled", TRUE, FALSE),
                     keysize=0.2, 
                     KeyValueName=titleTextb,
                     #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                     # key.par=list(mar=c(3,0,3,0)),
                     # lmat -- added 2 lattice sections (5 and 6) for padding
                     # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                     lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(1, 0.2, 5), lwid=c(0.2, 0.2, 5, 1))
        dev.off()
        # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
        
      } # whichChr
      
      # } # aCDPnormType
      
    } # aLib
    
    # } # aCompType
    
  } # aBins2Use
} # aThresh










##################################################
# 20210414
# *** ALL TIME POINTS AND ALL CELL STATES IN ONE PLOT ***
#
# Decay MAPS for active vs inactive cells using differentiated threshold
# Using LMD to indicate allelic direction of XCI
#
# BUT NOT BINARIZED VERSIONS.
# AND ONLY FOR BOTH ALLELES ALONGSIDE ONE ANOTHER.
# AND ONLY FOR LMD (aka LMD) XCI status.
# AND ONLY FOR SCALED CDPs.

# *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
aChr = "chrX"

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", aChr, sep="")), sep=""))
# if ( ! file.exists( subsupermainOutputDir ) ) {
#   if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
# }

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
# if ( ! file.exists( mainOutputDir ) ) {
#   if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
# }

decayCurveMaindir = file.path(mainOutputDir, "decayMaps")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    # for( aCompType in c( "LMD")) {
    aCompType = "LMD"
    
    libFlag = TRUE
    for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
      # aLib = "TsixStop_NPC"
      
      # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
      # for( aCDPnormType in c( "scaled" ) ) {  
      aCDPnormType ="scaled" 
      
      decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele.allXCIstates.V3", sep="."), aCompType)
      if ( ! file.exists( decayOutdir ) ) {
        if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
      }
      
      print(aLib)
      whatLib = sub("_$", "", 
                    ifelse(!grepl("_", aLib), aLib, 
                           sub("[^_]+$", "", 
                               sub("_rep.$", "", aLib) )))
      tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
      
      # 20210414 all libs
      if (libFlag) {
        matrix2Plot.all = list()
        XCIcolbar = list()
        TPcolbar = list()
        for (whichChr in c("chrX", "chr1")) {
          # whichChr = "chrX"
          matrix2Plot.all[[whichChr]] = c()
          XCIcolbar[[whichChr]] = c()
          TPcolbar[[whichChr]] = c()
        }
        libFlag = FALSE
      }
      
      # maxDiff = 1
      # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
      #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
      #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
      # }
      
      tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
      anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
      if ( aCompType == "LMD" ) {
        XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
      } else if ( grepl("dist", aCompType) ) {
        XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
      } else {
        XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
      }
      # head(anAllelicCDPcomp)
      
      LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
      refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
      altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
      
      # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
      # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
      # print(c(tenFPRcellNum, tenFPRrate))
      
      # *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
      for (whichChr in c("chrX", "chr1")) {
        # whichChr = "chrX"

        for (plotXCI in c(FALSE, TRUE)) {
          
          for  (plotRefAlt in c(TRUE, FALSE)) {
            # whichChr="chrX"; plotXCI=TRUE; plotRefAlt=TRUE;
            
            if (!plotXCI && plotRefAlt==FALSE) next;

            matrix2Plot.combined = c()
            
            for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
              # alleliclib = "TsixStop_NPC.ref"
              print(alleliclib)
              
              if(plotXCI) {
                if(plotRefAlt) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][refXCIboolIdx][ timepointThenCovOrd[[alleliclib]][refXCIboolIdx] ]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][altXCIboolIdx][ timepointThenCovOrd[[alleliclib]][altXCIboolIdx] ]
                }
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                # covOrd2Use = covOrd[[alleliclib]][!XCIboolIdx][ timepointThenCovOrd[[alleliclib]][!XCIboolIdx] ]
              }
              
              if ( is.null(matrix2Use) || length(matrix2Use) == 0 || length(matrix2Use) == 14 ) next
              # length(matrix2Use) == 14 ==> matrix has one row
              # ncol(matrix2Use) < 25  heatmaps look bad
              
              # Reorder cells randomnly
              set.seed(123)
              matrix2Use = matrix2Use[, sample(ncol(matrix2Use))]
              
              # Combine CDPs from both alleles
              matrix2Plot.combined = rbind( matrix2Plot.combined, matrix2Use )
              # if ( grepl( ".ref", alleliclib)) {
              #   covOrd2Use.combined = covOrd2Use
              # }
              
            } # alleliclib
            
            if (is.null(matrix2Plot.combined)) next
            
            matrix2Plot.all[[whichChr]]  = cbind(matrix2Plot.all[[whichChr]] , matrix2Plot.combined)
            
            # titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
            # nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
            
            XCIcol = "lightgreen" # active
            if (plotXCI && plotRefAlt) XCIcol = "darkblue"
            if (plotXCI && !plotRefAlt) XCIcol = "red4"
            XCIcolbar[[whichChr]] = c( XCIcolbar[[whichChr]] , rep( XCIcol, dim(matrix2Plot.combined)[2] ))
            TPcol = tpCols[[aLib]]
            TPcolbar[[whichChr]] = c( TPcolbar[[whichChr]], rep( TPcol, dim(matrix2Plot.combined)[2] ))
            # covcolbar = covCols[[lib]][covOrd[[aLib]]][timepointThenCovOrd[[aLib]]]
            # rColBars = t(cbind(TPcolbar, covcolbar))
            # rownames(rColBars)=c("sample", "coverage")
            
          } # plotRefAlt
        } # plotXCI
        
      } # whichChr
      
      # } # aCDPnormType
    
      if (aLib %in% mostDifferentiatedSample) {
        
        for (whichChr in c("chrX", "chr1")) {
          
          rColBars = t(cbind(TPcolbar[[whichChr]], XCIcolbar[[whichChr]]))
          rownames(rColBars)=c("time point", "XCI state")
          
          matrix2Plot = matrix2Plot.all[[whichChr]] 
          
          alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
          rownames(alleleBars)=c("allele")
          colnames(alleleBars)=rownames(matrix2Plot)
          
          titleText = paste(whichChr, aLib, aCDPnormType, sep=" ")
          nameText = paste(whichChr, aLib, aCDPnormType, sep=".")
          
          titleTextb = paste(titleText,  " bothAlleles  (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameTextb = paste(nameText, "bothAlleles", sep=".")
          
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          # cat(aCDPnormType, ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(aCDPnormType == "scaled") {
            cols2Use = heatCols2
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            # heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = heatColsV
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            # heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          # # 20210317
          # # Tried to make height proportional to number of rows.
          # # But doesn't really work b/c it also scales the color scale bar.
          # map_height = 4 + (4 * (dim(matrix2Plot)[2] / 1000))
          map_height = 12
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayMap.", 
                                                     nameTextb, ".combined.pdf", sep=""))
          
          pdf(PDFfilename, height=map_height, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), margins=c(5,0), labRow="", # labRow=cellLabels, 
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       ColSideColors=t(alleleBars), ColSideColorsSize=4,
                       symm=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                       symkey=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                       symbreaks=ifelse(aCDPnormType=="scaled", TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleTextb,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                       lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(0.5, 0.2, 5), lwid=c(0.2, 0.4, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
        } # whichChr
        
        libFlag = TRUE 
      } # if (aLib %in% mostDifferentiatedSample)
      
    } # aLib
    
    # } # aCompType
    
  } # aBins2Use
} # aThresh










##################################################
# 20210414
# *** ALL TIME POINTS AND ALL CELL STATES IN ONE PLOT ***
#
##################################################
# 20210501
# Included cell cycle and coverage
#
# Decay MAPS for active vs inactive cells using differentiated threshold
# Using LMD to indicate allelic direction of XCI
#
# BUT NOT BINARIZED VERSIONS.
# AND ONLY FOR BOTH ALLELES ALONGSIDE ONE ANOTHER.
# AND ONLY FOR LMD (aka LMD) XCI status.
# AND ONLY FOR SCALED CDPs.
#
##################################################

# *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
aChr = "chrX"

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", aChr, sep="")), sep=""))
# if ( ! file.exists( subsupermainOutputDir ) ) {
#   if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
# }

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
# if ( ! file.exists( mainOutputDir ) ) {
#   if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
# }

decayCurveMaindir = file.path(mainOutputDir, "decayMaps")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    # for( aCompType in c( "LMD")) {
    aCompType = "LMD"
    
    libFlag = TRUE
    for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
      # aLib = "TsixStop_NPC"
      
      # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
      # for( aCDPnormType in c( "scaled" ) ) {  
      aCDPnormType ="scaled" 
      
      decayOutdir = file.path(decayCurveMaindir, paste( aCDPnormType, "XCIvsNonXCI.differentiatedThreshold.byAllele.allXCIstates.V4", sep="."), aCompType)
      if ( ! file.exists( decayOutdir ) ) {
        if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
      }
      
      print(aLib)
      whatLib = sub("_$", "", 
                    ifelse(!grepl("_", aLib), aLib, 
                           sub("[^_]+$", "", 
                               sub("_rep.$", "", aLib) )))
      tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
      
      # 20210414 all libs
      if (libFlag) {
        matrix2Plot.all = list()
        XCIcolbar = list()
        TPcolbar = list()
        # 20210501
        covColbar = list()
        cellCycleColbar = list()
        for (whichChr in c("chrX", "chr1")) {
          # whichChr = "chrX"
          matrix2Plot.all[[whichChr]] = c()
          XCIcolbar[[whichChr]] = c()
          TPcolbar[[whichChr]] = c()
          covColbar[[whichChr]] = c()
          cellCycleColbar[[whichChr]] = c()
        }
        libFlag = FALSE
      }
      
      # maxDiff = 1
      # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
      #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
      #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
      # }
      
      tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
      anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
      if ( aCompType == "LMD" ) {
        XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
      } else if ( grepl("dist", aCompType) ) {
        XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
      } else {
        XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
      }
      # head(anAllelicCDPcomp)
      
      LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
      refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
      altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
      
      # tenFPRrate = round(sum(XCIboolIdx)/length(XCIboolIdx) *100, 2)
      # tenFPRcellNum = paste0( "(", sum(XCIboolIdx), " / ", length(XCIboolIdx), ")" )
      # print(c(tenFPRcellNum, tenFPRrate))
      
      # *** NOTE: CDP maps for both chr1 and chrX plotted in same chrX directory ***
      for (whichChr in c("chrX", "chr1")) {
        # whichChr = "chrX"
        
        for (plotXCI in c(FALSE, TRUE)) {
          
          for  (plotRefAlt in c(TRUE, FALSE)) {
            # whichChr="chrX"; plotXCI=TRUE; plotRefAlt=TRUE;
            
            if (!plotXCI && plotRefAlt==FALSE) next;
            
            matrix2Plot.combined = c()
            
            for (alleliclib in paste(aLib, c("ref", "alt"), sep=".")) {
              # alleliclib = "TsixStop_NPC.ref"
              print(alleliclib)
              
              if(plotXCI) {
                if(plotRefAlt) {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,refXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][refXCIboolIdx][ timepointThenCovOrd[[alleliclib]][refXCIboolIdx] ]
                } else {
                  matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,altXCIboolIdx]
                  # covOrd2Use = covOrd[[alleliclib]][altXCIboolIdx][ timepointThenCovOrd[[alleliclib]][altXCIboolIdx] ]
                }
              } else {
                matrix2Use = dataMatrices[[aThresh]][[whichChr]][[bins2Use]][[aCDPnormType]][[alleliclib]][,!XCIboolIdx]
                # covOrd2Use = covOrd[[alleliclib]][!XCIboolIdx][ timepointThenCovOrd[[alleliclib]][!XCIboolIdx] ]
              }
              
              if ( is.null(matrix2Use) || length(matrix2Use) == 0 || length(matrix2Use) == 14 ) next
              # length(matrix2Use) == 14 ==> matrix has one row
              # ncol(matrix2Use) < 25  heatmaps look bad
              
              # Reorder cells randomnly
              set.seed(123)
              matrix2Use = matrix2Use[, sample(ncol(matrix2Use))]
              
              # Combine CDPs from both alleles
              matrix2Plot.combined = rbind( matrix2Plot.combined, matrix2Use )
              # if ( grepl( ".ref", alleliclib)) {
              #   covOrd2Use.combined = covOrd2Use
              # }
              
            } # alleliclib
            
            if (is.null(matrix2Plot.combined)) next
            
            matrix2Plot.all[[whichChr]]  = cbind(matrix2Plot.all[[whichChr]] , matrix2Plot.combined)
            
            # titleText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=" ")
            # nameText = paste(whichChr, aLib, aCDPnormType, ifelse(plotXCI, "XCI", "notXCI"), sep=".")
            
            XCIcol = "lightgreen" # active
            if (plotXCI && plotRefAlt) XCIcol = "darkblue"
            if (plotXCI && !plotRefAlt) XCIcol = "red4"
            XCIcolbar[[whichChr]] = c( XCIcolbar[[whichChr]] , rep( XCIcol, dim(matrix2Plot.combined)[2] ))
            
            TPcol = tpCols[[aLib]]
            TPcolbar[[whichChr]] = c( TPcolbar[[whichChr]], rep( TPcol, dim(matrix2Plot.combined)[2] ))
            
            # 20210501
            # coverage colors
            # covColIdx = match( colnames(matrix2Plot.combined), 
            #                    names(covCols.biallelic[[paste0(aLib, ".ref")]]) )
            # covColbar[[whichChr]] = c( covColbar[[whichChr]],  
            #                            covCols.biallelic[[paste0(aLib, ".ref")]][covColIdx] )
            bialllelic_total_counts = colSums(matrixCsums.biallelic[[aLib]], na.rm=TRUE)
            covColIdx = match( colnames(matrix2Plot.combined), names(bialllelic_total_counts) )
            covColbar[[whichChr]] = c( covColbar[[whichChr]],
                                       bialllelic_total_counts[covColIdx] )
            
            # 20210501
            # autosomal CDP cell cycle group
            CDPsInClustersIdx = match( sub( paste0(aLib, ".ref."), "", colnames(matrix2Plot.combined) ), 
                                       kClusters[[aLib]]$cellID )
            cellCycleColbar[[whichChr]] = c( cellCycleColbar[[whichChr]],  
                                             kClustColors[ kClusters[[aLib]]$clusterID[CDPsInClustersIdx] ]  )
            
            # rColBars = t(cbind(TPcolbar, covcolbar))
            # rownames(rColBars)=c("sample", "coverage")
            
          } # plotRefAlt
        } # plotXCI
        
      } # whichChr
      
      # } # aCDPnormType
      
      if (aLib %in% mostDifferentiatedSample) {
        
        for (whichChr in c("chrX", "chr1")) {
          
          covColOrd = order(covColbar[[whichChr]])
          covColbarCols = viridis(length(covColbar[[whichChr]]), alpha=0.5)[covColOrd]

          rColBars = t(cbind(TPcolbar[[whichChr]], 
                             XCIcolbar[[whichChr]], 
                             cellCycleColbar[[whichChr]],
                             covColbarCols) )
          rownames(rColBars) = c("time point", "XCI state", "cell cycle group", "total contacts")
          
          matrix2Plot = matrix2Plot.all[[whichChr]] 
          
          alleleBars = t( as.matrix( c( rep("blue", dim(matrix2Plot)[1]/2), rep("red", dim(matrix2Plot)[1]/2) ) ) )
          rownames(alleleBars)=c("allele")
          colnames(alleleBars)=rownames(matrix2Plot)
          
          titleText = paste(whichChr, aLib, aCDPnormType, sep=" ")
          nameText = paste(whichChr, aLib, aCDPnormType, sep=".")
          
          titleTextb = paste(titleText,  " bothAlleles  (n = ", dim(matrix2Plot)[2], ")", sep="")
          nameTextb = paste(nameText, "bothAlleles", sep=".")
          
          minVal = min(matrix2Plot, na.rm=TRUE)
          maxVal = max(matrix2Plot, na.rm=TRUE)
          # cat(aCDPnormType, ": ", minVal, " --> ", maxVal, "\n", sep="" )
          
          if(aCDPnormType == "scaled") {
            cols2Use = heatCols2
            # minVal=-4; maxVal=4
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            # heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
          } else {
            cols2Use = heatColsV
            # minVal=log10(1/dim(matrix2Plot)[1]) # 10^-2 = proportion of 1/100
            # minVal = -3; maxVal = -1.5 
            breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
            # heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
          }
          
          # # 20210317
          # # Tried to make height proportional to number of rows.
          # # But doesn't really work b/c it also scales the color scale bar.
          # map_height = 4 + (4 * (dim(matrix2Plot)[2] / 1000))
          map_height = 12
          
          # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
          PDFfilename = file.path(decayOutdir, paste(projectName, ".", aCompType, ".decayMap.", 
                                                     nameTextb, ".combined.pdf", sep=""))
          
          pdf(PDFfilename, height=map_height, width=12)
          
          # heatmap.2(t(matrix2Plot), col=heatCols2,
          #           Rowv=TRUE, Colv=FALSE,
          #           dendrogram="none",
          #           scale="none", key=TRUE, density.info="none",
          #           trace="none",cexCol=1.2, labRow=NA, labCol=NA,
          #           RowSideColors=kClustColbar,
          #           symm=TRUE,symkey=TRUE,symbreaks=TRUE)
          
          # https://www.biostars.org/p/18211/
          # source("~/proj/2019_sciHiC_analysis/results/gbonora/20180819_sciHiC_FUCCI/analyzeContactDecayRates/heatmap.3.GB.R")
          heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                       Rowv=FALSE, Colv=FALSE,
                       dendrogram="none",
                       scale="none", key=TRUE, density.info="none",
                       trace="none", cexCol=0.95, labCol=c(binDistDesc, binDistDesc), margins=c(5,0), labRow="", # labRow=cellLabels, 
                       RowSideColors=rColBars, RowSideColorsSize=4,
                       ColSideColors=t(alleleBars), ColSideColorsSize=4,
                       symm=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                       symkey=ifelse(aCDPnormType=="scaled",TRUE, FALSE),
                       symbreaks=ifelse(aCDPnormType=="scaled", TRUE, FALSE),
                       keysize=0.2, 
                       KeyValueName=titleTextb,
                       #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
                       # key.par=list(mar=c(3,0,3,0)),
                       # lmat -- added 2 lattice sections (5 and 6) for padding
                       # lmat=rbind(c(0, 0, 5, 4), c(0, 1, 2, 3)), lhei=c(1, 5), lwid=c(0.5, 0.5, 5, 1))
                       lmat=rbind(c(0, 0, 6, 5), c(0, 0, 2, 0), c(0, 1, 3, 4)), lhei=c(0.5, 0.2, 5), lwid=c(0.2, 0.4, 5, 1))
          dev.off()
          # layout(mat=rbind(c(0, 6, 5, 4), c(0, 1, 2, 3)), heights=c(1, 5), widths=c(0.2, 1, 5, 0.2)); layout.show(6) 
          
        } # whichChr
        
        libFlag = TRUE 
      } # if (aLib %in% mostDifferentiatedSample)
      
    } # aLib
    
    # } # aCompType
    
  } # aBins2Use
} # aThresh










##################################################
# 20200812
# Cell IDs for XCI and non-XCI cells
# Using LMD to indicate allelic direction of XCI
#
# 20201118
# Calculate overlap of LMD and Spearman XCI Patski cells

XCItableMaindir = file.path("../data/", "XCItables")
if ( ! file.exists( XCItableMaindir ) ) {
  if ( ! dir.create( XCItableMaindir, showWarnings=FALSE, recursive=TRUE )  ) XCItableMaindir = mainOutputDir
}

mostDifferentiatedSample = c("EBdiff_F121_d20",
                             "EBdiff_F123_d20",
                             "EBdiff_TsixStop_d20",
                             "TsixStop_NPC",
                             "F121_NPC",
                             "Patski",
                             "combined")
names(mostDifferentiatedSample) = c("EBdiff_F121",
                                    "EBdiff_F123",
                                    "EBdiff_TsixStop",
                                    "TsixStop",
                                    "F121",
                                    "Patski",
                                    "combined")

for ( aThresh in names(tenFPRrateThresholds) ) {
  
  for( aBins2Use in names(tenFPRrateThresholds[[aThresh]]) ) {
    
    # 20201118
    XCIcells = list()
    XCIcells.bool = list()
    
    # for( aCompType in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]]) ) {
    for( aCompType in c( "Spearman", "LMD")) {
      
      # # 2020425
      # propXCItable = list()
      # propXCItable[["scaled"]] = c()
      # propXCItable[["prop"]] = c()
      
      # 20201118
      XCIcells[[aCompType]] = list()
      XCIcells.bool[[aCompType]] = list()
      
      for( aLib in names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]]) ) {
        
        # for( aCDPnormType in rev(names(tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[aLib]])) ) {
        for( aCDPnormType in c( "scaled" ) ) {  
          
          # aLib = "F121_NPC"
          print(aLib)
          whatLib = sub("_$", "", 
                        ifelse(!grepl("_", aLib), aLib, 
                               sub("[^_]+$", "", 
                                   sub("_rep.$", "", aLib) )))
          tenFPRrateThreshLib = mostDifferentiatedSample[whatLib]
          
          # maxDiff = 1
          # for( bLib in grep(whatLib, names(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]]), value=TRUE) ) {
          #   thisMax = min(allelicCDPcompList[[aThresh]][["chrX"]][[aBins2Use]][[aCompType]][[bLib]][[aCDPnormType]], na.rm=TRUE)
          #   maxDiff = ifelse(thisMax < maxDiff, thisMax, maxDiff)
          # }
          
          tenFPRrateThresh = tenFPRrateThresholds[[aThresh]][[aBins2Use]][[aCompType]][[tenFPRrateThreshLib]][[aCDPnormType]]
          aChr = "chrX"
          anAllelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][[aCompType]][[aLib]][[aCDPnormType]]
          if ( aCompType == "LMD" ) {
            XCIboolIdx = abs(anAllelicCDPcomp) >= tenFPRrateThresh
          } else if ( grepl("dist", aCompType) ) {
            XCIboolIdx = anAllelicCDPcomp >= tenFPRrateThresh
          } else {
            XCIboolIdx = anAllelicCDPcomp <= tenFPRrateThresh
          }
          # head(anAllelicCDPcomp)
          
          LMDallelicCDPcomp = allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]
          refXCIboolIdx = (LMDallelicCDPcomp >= 0) & XCIboolIdx
          altXCIboolIdx = (LMDallelicCDPcomp < 0) & XCIboolIdx
          
          # 20201118
          XCIcells[[aCompType]][[aLib]] = names(XCIboolIdx[XCIboolIdx])
          XCIcells.bool[[aCompType]][[aLib]] = XCIboolIdx

          #--------------------------------------------------------------#
          # 20200812
          # Cell IDs for XCI and non-XCI cells
          
          XCItableOutdir = file.path(XCItableMaindir, paste( aCompType, ".", aCDPnormType, ".XCIvsNonXCI.differentiatedThreshold.byAllele", sep=""))
          if ( ! file.exists( XCItableOutdir ) ) {
            if ( ! dir.create( XCItableOutdir, showWarnings=FALSE, recursive=TRUE )  ) XCItableOutdir = XCItableMaindir
          }
          
          tableFileName = file.path(XCItableOutdir, paste(aLib, ".nonXCIcells.tsv", sep=""))
          write.table(x=names(XCIboolIdx[XCIboolIdx==FALSE]),
                      file=tableFileName, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
          
          tableFileName = file.path(XCItableOutdir, paste(aLib, ".refXCIcells.tsv", sep=""))
          write.table(x=names(XCIboolIdx[refXCIboolIdx==TRUE]),
                      file=tableFileName, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
          
          tableFileName = file.path(XCItableOutdir, paste(aLib, ".altXCIcells.tsv", sep=""))
          write.table(x=names(altXCIboolIdx[altXCIboolIdx==TRUE]),
                      file=tableFileName, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        } # aCDPnormType
      } # aLib

    } # aCompType
    
    # 20201118
    # Check cell IDs are in the same order.
    all(names(allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["LMD"]][[aLib]][[aCDPnormType]]) == 
          names(allelicCDPcompList[[aThresh]][[aChr]][[aBins2Use]][["Spearman"]][[aLib]][[aCDPnormType]]))
    # [1] TRUE
    
    # Calculate overlap of LMD and Spearman XCI Patski cells
    # ***** Note: *****
    # Discovered an issue in the calculuaton of 'tenFPRrateThresh':
    # Was previoulsy using the AUC table rownames which were rounded thresholds.
    # For consitent results (XCI calls) do not round.
     for ( aLib in names(XCIcells[[aCompType]]) ) {
        overlap = sum( !is.na( match(XCIcells[["Spearman"]][[aLib]], XCIcells[["LMD"]][[aLib]]) ) ) 
        cat(aLib, ":\n",
            "XCI cells by Spearman: ", length(XCIcells[["Spearman"]][[aLib]]), "\n",
            "XCI cells by LMD: ", length(XCIcells[["LMD"]][[aLib]]), "\n",
            "overlap: ", overlap, "\n",
            "% XCI cells by Spearman: ", round(overlap / length(XCIcells[["Spearman"]][[aLib]]) * 100, 2), "\n",  
            "% XCI cells by LMD: ", round(overlap / length(XCIcells[["LMD"]][[aLib]]) * 100, 2), "\n",
            sep=""
            )
     }
    # Patski:
    # XCI cells by Spearman: 917
    # XCI cells by LMD: 943
    # overlap: 815
    # % XCI cells by Spearman: 88.88
    # % XCI cells by LMD: 86.43
    # 
    # combined:
    # XCI cells by Spearman: 917
    # XCI cells by LMD: 943
    # overlap: 815
    # % XCI cells by Spearman: 88.88
    # % XCI cells by LMD: 86.43
    
    # Calculate overlap of LMD and Spearman XCI Patski cells using Boolean. Should be same as above and is.
    for ( aLib in names(XCIcells[[aCompType]]) ) {
      overlap = sum(XCIcells.bool[["Spearman"]][[aLib]] & XCIcells.bool[["LMD"]][[aLib]])
      cat(aLib, ":\n",
          "XCI cells by Spearman: ", sum(XCIcells.bool[["Spearman"]][[aLib]]), "\n",
          "XCI cells by LMD: ", sum(XCIcells.bool[["LMD"]][[aLib]]), "\n",
          "overlap: ", overlap, "\n",
          "% XCI cells by Spearman: ", round(overlap / sum(XCIcells.bool[["Spearman"]][[aLib]]) * 100, 2), "\n",  
          "% XCI cells by LMD: ", round(overlap / sum(XCIcells.bool[["LMD"]][[aLib]]) * 100, 2), "\n",
          sep=""
      )
    }
    # Patski:
    # XCI cells by Spearman: 917
    # XCI cells by LMD: 943
    # overlap: 815
    # % XCI cells by Spearman: 88.88
    # % XCI cells by LMD: 86.43
    # 
    # combined:
    # XCI cells by Spearman: 917
    # XCI cells by LMD: 943
    # overlap: 815
    # % XCI cells by Spearman: 88.88
    # % XCI cells by LMD: 86.43
    
    # So 89% of Patski cells called inactive using the Spearman correlation were called inactive using the LMD method.
    
  } # aBins2Use
} # aThresh










##################################################
# 20201002 
# Decay curves for distal F121 cells along pseudotime branches

whichTrajectory="genome"; DEqval="DEqval1e-11"
# whichTrajectory="autosomes"; DEqval="DEqval1e-11"
# whichTrajectory="chrX"; DEqval="DEqval1e-05"

# Load table of distal cell IDs
tablefilename = paste0( "~/proj/2019_sciHiC_analysis/results/gbonora/20191105_sciHiC_allSamples_Nmasked/analyzeAllelicSciHiCdata.monocle2.catenatedAllelicVectors/",
                        "allelicAnalysisOutput.monocle2.catenatedAllelicVectors.F121only.", whichTrajectory, 
                        "/trajectoryAnalysis.semisupervised/",
                        "20191105_allelic_sciHiC_allSamples_Nmasked.monocle2results.trajectoryAnalysis.minExprThreshold0.1.numCellsExpressed10.downSampleProp1.", 
                        DEqval, ".pseudotime_per_branch.distal_cells.tsv" )
                        
ptData2 = read.table(tablefilename, sep="\t", header = TRUE)
head(ptData2)

branch_labels = sort( unique ( ptData2$branch_labels ))


for ( nThresh in names(allelicCDPcompList) ) {
  
  for ( whichChr in names(allelicCDPcompList[[nThresh]]) ) {
    
    subsupermainOutputDir = file.path(supermainOutputDir, 
                                      paste( "countThresh", nThresh, ifelse(whichChr=="all", "", paste(".", whichChr, sep="")), sep=""))
    if ( ! file.exists( subsupermainOutputDir ) ) {
      if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
    }
    
    for( bins2Use in names(allelicCDPcompList[[nThresh]][[whichChr]]) ) {
      
      mainOutputDir = file.path(subsupermainOutputDir, 
                                bins2Use)
      if ( ! file.exists( mainOutputDir ) ) {
        if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
      }  
      
      decayCurveMaindir = file.path(mainOutputDir, 
                                    paste0( "decayCurves.pseudotime_distal.", whichTrajectory, "_trajectory" ) )
      if ( ! file.exists( decayCurveMaindir ) ) {
        if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
      }
      
      lwd2Use = 0.5
      alpha2Use=0.01
      
      for (useScaled in c(FALSE, TRUE)) {
        
        decayOutdir = file.path(decayCurveMaindir, ifelse(useScaled, "scaled", "prop"))
        if ( ! file.exists( decayOutdir ) ) {
          if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
        }
        
        if (!useScaled) {
          ylims.prop = range(unlist(dataMatrix.prop), na.rm=TRUE)
        }
        
        for (lib in c( whichLibs, "F121") ) {
          # lib = "F121_d3"
          print(lib)
          
          titleText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), sep=" ")
          nameText = paste(whichChr, lib, ifelse(useScaled, "scaled", "prop"), "pseudotime_distal_vs_nondistal", 
                           paste0( whichTrajectory, "_trajectory" ),
                           sep=".")
          
          PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.overlaid.v2.", nameText, ".refValt.pdf", sep=""))
          pdf(PDFfilename, height=length(branch_labels) * 3, width=9)
          
          par(mfrow=c(length(branch_labels),2))
          
          for (aBranch in branch_labels) {
            
            for (isDistal in c(TRUE, FALSE) ) {
              
              meanCDP = list()
              nuPlotflag=TRUE  
              
              for (allele in c("ref", "alt")) {
                alleliclib = paste(lib, allele, sep=".")
                print(paste(aBranch, isDistal, alleliclib, nuPlotflag))
                
                if (lib == "F121") {
                  matrix2Use = c()
                  for (alib in whichLibs) {
                    if (useScaled) {
                      matrix2Use = cbind(matrix2Use, 
                                         dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[ paste(alib, allele, sep=".")]])
                    } else {
                      matrix2Use = cbind(matrix2Use, 
                                         dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[ paste(alib, allele, sep=".")]])
                    }
                  }
                } else { 
                  # if (!useScaled) {
                  #   matrix2Use = dataMatrix.prop[[alleliclib]]
                  # } else {
                  #   matrix2Use = dataMatrix.scaled[[alleliclib]]
                  # }
                  if (useScaled) {
                    matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
                  } else {
                    matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
                  }
                }
                
                distal_cell_ids = ptData2$cell[ ptData2$branch_labels == aBranch & ptData2$distal==isDistal ] 
                isDistal_boolidx = !is.na( match( sub( ".alt", "", sub( ".ref", "", colnames(matrix2Use) ) ), 
                                                  distal_cell_ids ) )
                matrix2Use = matrix2Use[,isDistal_boolidx]
                
                nCells = ncol(matrix2Use)
                print(paste( sum(isDistal_boolidx), nCells ))
                
                if ( is.null(nCells) || nCells==0 ) {
                  if ( nuPlotflag ) {
                    nCells = 0
                    if (useScaled) {
                      ylabText = "Scaled Cis counts"
                      ylims = c(-3, +3)
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = yat
                    } else {
                      ylabText = "Proportion of Cis contact (log10)"
                      # ylims = c(-3, -1.5)
                      ylims = ylims.prop
                      yat = seq(ylims[1], ylims[2], 0.5)
                      ylabels = signif(10^yat, 2)
                    }
                    plot( 0, type="l", col=col2Use, lwd=lwd2Use,
                          main=paste(titleText, " - ", aBranch, " branch\n",
                                     ifelse(isDistal, "pseudotime distal", "pseudotime non-distal"), 
                                     " (n=", nCells, ")", sep=""),
                          # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                          xlab="", ylab="",
                          ylim=ylims, xaxt="n", yaxt="n")
                    axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                    axis(2, at=yat, labels=ylabels)
                    nuPlotflag = FALSE
                  }
                  
                } else {
                  
                  # alaph2Use =  1/nCells
                  col2Use = setAlpha(alleleCols[sub("^.*\\.", "", alleliclib)], alpha2Use)
                  
                  for (iCell in 1:nCells) {
                    y = matrix2Use[,iCell]
                    if (iCell==1 && nuPlotflag){
                      if (useScaled) {
                        ylabText = "Scaled Cis counts"
                        ylims = c(-3, +3)
                        yat = seq(ylims[1], ylims[2], 0.5)
                        ylabels = yat
                      } else {
                        ylabText = "Proportion of Cis contact (log10)"
                        # ylims = c(-3, -1.5)
                        ylims = ylims.prop
                        yat = seq(ylims[1], ylims[2], 0.5)
                        ylabels = signif(10^yat, 2)
                      }
                      plot( y, type="l", col=col2Use, lwd=lwd2Use,
                            main=paste(titleText, " - ", aBranch, " branch\n",
                                       ifelse(isDistal, "pseudotime distal", "pseudotime non-distal"), 
                                       " (n=", nCells, ")", sep=""),
                            # xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                            xlab="", ylab="",
                            ylim=ylims, xaxt="n", yaxt="n")
                      axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
                      axis(2, at=yat, labels=ylabels)
                      nuPlotflag = FALSE
                    } else {
                      lines(y, col=col2Use, lwd=lwd2Use)
                    }
                  }
                  y = rowMeans(matrix2Use, na.rm = TRUE)
                  lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                  meanCDP[[alleliclib]] = y
                } # alleliclib
                
                # alleliclib = sub("\\.alt$", ".ref", alleliclib)
                # if (!useScaled) {
                #   matrix2Use = dataMatrix.prop[[alleliclib]]
                # } else {
                #   matrix2Use = dataMatrix.scaled[[alleliclib]]
                # }
                # if (useScaled) {
                #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[alleliclib]]
                # } else {
                #   matrix2Use = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[alleliclib]]
                # }
                # distal_cell_ids = ptData2$cell[ ptData2$branch_labels == aBranch & ptData2$distal==isDistal ] 
                # isDistal_boolidx = !is.na( match( sub( ".alt", "", sub( ".ref", "", colnames(matrix2Use) ) ), 
                #                                   distal_cell_ids ) )
                # matrix2Use = matrix2Use[,isDistal_boolidx]
                # y = rowMeans(matrix2Use, na.rm = TRUE)
                # lines(y, col=alleleCols[sub("^.*\\.", "", alleliclib)], lwd=2.5)
                # legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
                # interAllelicCDPcorr = cor( meanCDP[[1]], meanCDP[[2]], method='spearman', use='pairwise.complete.obs' )
                
                if (length(meanCDP) == 2) {
                  
                  if (lib == "F121") {
                    matrix2Use.ref = c()
                    matrix2Use.alt = c()
                    for (alib in whichLibs) {
                      if (useScaled) {
                        matrix2Use.ref = cbind(matrix2Use.ref, 
                                               dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(alib, "ref", sep=".")]])
                        aalleliclib = paste(lib, "alt", sep=".")
                        matrix2Use.alt = cbind(matrix2Use.ref, 
                                               dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(alib, "alt", sep=".")]])
                      } else {
                        matrix2Use.ref = cbind(matrix2Use.ref, 
                                               dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(alib, "ref", sep=".")]])
                        aalleliclib = paste(lib, "alt", sep=".")
                        matrix2Use.alt = cbind(matrix2Use.ref, 
                                               dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(alib, "alt", sep=".")]])
                      }
                    }
                  } else { 
                    # if (!useScaled) {
                    #   matrix2Use.ref = dataMatrix.prop[[paste(lib, "ref", sep=".")]]
                    #   matrix2Use.alt = dataMatrix.prop[[paste(lib, "alt", sep=".")]]
                    # } else {
                    #   matrix2Use.ref = dataMatrix.scaled[[paste(lib, "ref", sep=".")]]
                    #   matrix2Use.alt = dataMatrix.scaled[[paste(lib, "alt", sep=".")]]
                    # }
                    if (useScaled) {
                      matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "ref", sep=".")]]
                      matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["scaled"]][[paste(lib, "alt", sep=".")]]
                    } else {
                      matrix2Use.ref = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "ref", sep=".")]]
                      matrix2Use.alt = dataMatrices[[nThresh]][[whichChr]][[bins2Use]][["prop"]][[paste(lib, "alt", sep=".")]]
                    }
                  }
                  
                  isDistal_boolidx = !is.na( match( sub( ".alt", "", sub( ".ref", "", colnames(matrix2Use.ref) ) ), 
                                                    distal_cell_ids ) )
                  matrix2Use.ref = matrix2Use.ref[,isDistal_boolidx]
                  matrix2Use.alt = matrix2Use.alt[,isDistal_boolidx]
                  
                  CDPspearmans = c()
                  nCells = ncol(matrix2Use.ref)
                  for(iCell in 1:nCells) {
                    CDPspearman = cor(matrix2Use.ref[,iCell], matrix2Use.alt[,iCell],
                                      method='spearman', use='pairwise.complete.obs')
                    CDPspearmans = c(CDPspearmans, CDPspearman)
                  }
                  
                  meanInterAllelicCDPcorr = mean(CDPspearmans, na.rm=TRUE)
                  medianInterAllelicCDPcorr = median(CDPspearmans, na.rm=TRUE)
                  
                  legend("bottom", legend=c(
                    paste("spearman between mean CDPs = ", round(interAllelicCDPcorr, 2), sep=""),
                    paste("median spearman (individual CDPs) = ", round(medianInterAllelicCDPcorr, 2), sep="") ), 
                    bty='n')
                }
                
              } # nCell>0
              
            } # isDistal
            
          } # aBranch
          
          legend("topleft", legend=paste(lib, c("ref", "alt"), sep="."), fill=alleleCols)
          
          dev.off() 
          
        } # libs
        
      } # useScaled
      
    } # bins2Use
    
  } # whichChr
  
} # nThresh

# cd ~/proj/2019_sciHiC_analysis/results/gbonora/20191105_sciHiC_allSamples_Nmasked/analyzeAllelicContactDecayRates/analysisOutput.loadAllChromData 
# mkdir decayCurves.scaled.pseudotime_distal_vs_nondistal
# cd decayCurves.scaled.pseudotime_distal_vs_nondistal/
# find .. -type f -name "*.scaled.pseudotime_distal_vs_nondistal.genome_trajecotry.refValt.pdf" | while read FILE; do cp -pr ${FILE} .; done










##################################################
# 20200506
# Save concatenated allelic broad-bined CDPs for F121s

alleleCDPmatrixDir = "../data/allelic_catenated_contact_decay_profiles/"
if ( ! file.exists( alleleCDPmatrixDir ) ) {
  if ( ! dir.create( alleleCDPmatrixDir, showWarnings=FALSE, recursive=TRUE )  ) alleleCDPmatrixDir = supermainOutputDir
}


theCellType="F121"
theLibs=c(
  "F121_d0",
  "F121_d3",
  "F121_d7",
  "F121_d11",
  "F121_NPC"
)

# 20200921
theCellType="TsixStop"
theLibs=c(
  "TsixStop_d0",
  "TsixStop_d3",
  "TsixStop_d7",
  "TsixStop_d11",
  "TsixStop_NPC"
)

# 20200921
theCellType="F123"
theLibs=c(
  "F123_d0",
  "F123_NPC"
)

# 20200922
theCellType="F121andF123"
theLibs=c(
  "F121_d0",
  "F121_d3",
  "F121_d7",
  "F121_d11",
  "F121_NPC",
  "F123_d0",
  "F123_NPC"
)

# 20210323
theCellType="oldF121andTsixStop"
theLibs=c(
  "EBdiff_F121_d0",   
  "EBdiff_F121_d7",     
  "EBdiff_F121_d20",     
  "EBdiff_TsixStop_d0", 
  "EBdiff_TsixStop_d7", 
  "EBdiff_TsixStop_d20"
)

for ( theChrom in names(dataMatrix.PF) ) {
  
  theAlelleicCDPmatrix = c()
  
  for( theLib in theLibs ) {
    
    # theChrom="chrX"; theLib="F121_d0"
    allelicCDPvectors = c()
    for ( theAllele in c("ref", "alt") ) {
      allelicCDPvectors = rbind( allelicCDPvectors,
                                 dataMatrix.PF[[theChrom]][[paste(theLib, theAllele, sep=".")]] )
    }
    colnames(allelicCDPvectors) = sub( ".ref", "", colnames(allelicCDPvectors) )
    theAlelleicCDPmatrix = cbind(theAlelleicCDPmatrix, allelicCDPvectors)
  }
  rownames(theAlelleicCDPmatrix) = paste( c( rep( "ref", nrow(allelicCDPvectors)/2 ),
                                             rep( "alt", nrow(allelicCDPvectors)/2 ) ),
                                          rownames(allelicCDPvectors), 
                                          sep="_" )
  
  print( table( sub("\\..*$", "", colnames(theAlelleicCDPmatrix)) ) )
  
  # sciHi-C - filtered by chr1 & chrX allelic threshold
  # F121_d0 F121_d11  F121_d3  F121_d7 F121_NPC
  # 986      164      435      243      325
  #
  # 2020921 - filtered by chr1 & chrX allelic threshold
  # TsixStop_d0 TsixStop_d11  TsixStop_d3  TsixStop_d7 TsixStop_NPC 
  # 89           83           86          159          975 
  #
  # 2020921 - filtered by chr1 allelic threshold only
  # # F123_d0 F123_NPC 
  # # 326      216 
  # After resolving some filtering bug:
  # F123_d0 F123_NPC 
  # 1494      876 
  #
  # 2020922 - Both F121 and F123 filtered by chr1 allelic threshold only
  # F121_d0 F121_d11  F121_d3  F121_d7 F121_NPC  F123_d0 F123_NPC 
  # 1932      496      735      594      841     1494      876 
  #
  # 20210323
  # EBdiff_F121_d0     EBdiff_F121_d20      EBdiff_F121_d7 
  # 393                  89                  47 
  # EBdiff_TsixStop_d0 EBdiff_TsixStop_d20  EBdiff_TsixStop_d7 
  # 291                  80                  51 
  
  nameText = paste0("catenatedAllelicCDPs.", theCellType, ".", theChrom)
  tableFileName = file.path(alleleCDPmatrixDir, paste0(projectName, ".", nameText, ".tsv"))
  write.table(x=theAlelleicCDPmatrix,
              file=tableFileName, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  # system(paste0("gzip -f ", tableFileName))
  try(system(paste0("gzip -f \"", sub("~/", "$HOME/", tableFileName), "\"")))

}
