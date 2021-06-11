##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200106 non-allelic sciHic contact decay analysis
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



##################################################
# Functions

# smoothVector
smoothVector <- function(x, bins=3 ){
  return( shift(smth(x, "sma", bins), places=floor(bins/2), dir="left") )
}

parseIBC <- function(x){ return(strsplit(x=as.character(x), split="-")[[1]][1]) }
convertIbc2tp <- function(x){ return(ibc2tp[ibc2tp$ibc == as.character(x), 2]) }


# https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# https://www.biostars.org/p/18211/
# Load latest version of heatmap.3 function
# require(devtools)
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
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










##################################################
##################################################
##################################################
# Inputs and options

projectName = "non-allelic_sci-Hi-C"

inputDir = "../data/non-allelic_contact_decay_profiles/"

supermainOutputDir = "../results/non-allelic_sci-Hi-C/"
if ( ! file.exists( supermainOutputDir ) ) {
  if ( ! dir.create( supermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) supermainOutputDir = getwd()
}
setwd(supermainOutputDir)
getwd()

# # 20210310
# combined_label="combined"

# ---------------------------------------------------------------------------
# 20210226
# In the interest of time, just the differentiated cells for cell cycle analysis
libs=c(
"Patski",
"EBdiff_TsixStop_d20",
# "mES_EBdiff_D20",
"EBdiff_F121_d20",
"F121_NPC",
"TsixStop_NPC"
)

# timepointIDfile = "~/refdata/sciHiC_barcodeData/inner_barcode.timepoint.txt"
# ibc2tp = read.table(timepointIDfile)
# names(ibc2tp) = c("ibc", "tp")

# # Use this to color Patski libraries instead.
# tpCols=gray((1:4)/6)
# names(tpCols) = libs
tpCols = c(
    # Patski
    "gray",
    
    # oldTsix,
    "darkblue",
    "darkblue",
    
    # Tsix
    "black",
    
    # F121
    "black"
  )
names(tpCols) = libs
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 20210306
# Run for ESCs too... won't see the Xi but will see the Xa stucture re-emerge.
# And will see mitotic bands for group 1 cells.
libs=c(
  "EBdiff_TsixStop_d0",
  "EBdiff_F121_d0",
  "F121_d0",
  "TsixStop_d0"
)

# timepointIDfile = "~/refdata/sciHiC_barcodeData/inner_barcode.timepoint.txt"
# ibc2tp = read.table(timepointIDfile)
# names(ibc2tp) = c("ibc", "tp")

# # Use this to color Patski libraries instead.
# tpCols=gray((1:4)/6)
# names(tpCols) = libs
tpCols = c(
  # oldTsix,
  "red",
  "red",
  
  # Tsix
  "red",
  
  # F121
  "red"
)
names(tpCols) = libs
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 20210308
# Run for all F121 time points -- to be able to exclude mitotic cells from traj analysis.
libs=c(
  "F121_d0",
  "F121_d3",
  "F121_d7",
  "F121_d11",
  "F121_NPC"
)

# timepointIDfile = "~/refdata/sciHiC_barcodeData/inner_barcode.timepoint.txt"
# ibc2tp = read.table(timepointIDfile)
# names(ibc2tp) = c("ibc", "tp")

# # Use this to color Patski libraries instead.
# tpCols=gray((1:4)/6)
# names(tpCols) = libs
tpCols = c(
  "red",
  "orange",
  "green",
  "blue",
  "black"
)
names(tpCols) = libs

# 20210310
combined_label="F121_combined"

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 20210501
# Run for all TsixStop time points -- to be able to exclude mitotic cells from traj analysis.
libs=c(
  "TsixStop_d0",
  "TsixStop_d3",
  "TsixStop_d7",
  "TsixStop_d11",
  "TsixStop_NPC"
)

# timepointIDfile = "~/refdata/sciHiC_barcodeData/inner_barcode.timepoint.txt"
# ibc2tp = read.table(timepointIDfile)
# names(ibc2tp) = c("ibc", "tp")

# # Use this to color Patski libraries instead.
# tpCols=gray((1:4)/6)
# names(tpCols) = libs
tpCols = c(
  "red",
  "orange",
  "green",
  "blue",
  "black"
)
names(tpCols) = libs

# 20210501
combined_label="TsixStop_combined"

# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 20210502
# Patski alone.
# Note: For Patski cell cycle plots colored by XCI status 
#       use XCI analysis based on 100 contacts per chr1 and chrX allele rather than 50.

libs=c(
  "Patski"
)

# timepointIDfile = "~/refdata/sciHiC_barcodeData/inner_barcode.timepoint.txt"
# ibc2tp = read.table(timepointIDfile)
# names(ibc2tp) = c("ibc", "tp")

# # Use this to color Patski libraries instead.
# tpCols=gray((1:4)/6)
# names(tpCols) = libs
tpCols = c(
  # Patski
  "gray"
)
names(tpCols) = libs

# ---------------------------------------------------------------------------


# In-silico cell phasing over the cell cycle
# We counted the number of cis (intra-chromosomal) contacts per cell, binning contacts by distance into logarithmic bins 
# (143 bins, first one for contacts distanced < 1Kb, then each bin covers an exponent step of 0.125, using base 2). 
# log2(1000) ~ 9.96 ~10 
# length(2^seq(10, 27.625, 0.125)/(10^6)) + 1 (trans) = 143
# 2^27.625 = 206,992,033
# 2^9.875 = 939.0121
#
# Actually saved data to 2^27.750 = 225,726,413 Mb, or 143 bins + 1 trans = 144
binDist = c( 2^seq(10, 27.750, 0.125), Inf)
length(binDist)
# 144

plot3Dscatters=FALSE
# plot3Dscatters=TRUE

heatCols1 = heat.colors(256)
heatCols2 <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                    "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                    "#4393C3", "#2166AC", "#053061"))(256))
heatColsM = viridis(100, option="magma")
heatColsV = viridis(100, option="D")

# 20200108
# autosomalCDPs=FALSE
autosomalCDPs=TRUE










##################################################
##################################################
##################################################
# Load decay data

matrixFileName = list()
matrixData = list()

for(lib in libs) {
  # lib = "TsixStop_d3"
  
  # Counts file
  matrixFileName[[lib]] <- paste(inputDir,
                                 lib, 
                                 ifelse(autosomalCDPs, ".autosomes", ""), 
                                 ".logbinContactDecay", 
                                 ".mm10.txt.gz", sep="")
  matrixData[[lib]] <- read.table(gzfile(matrixFileName[[lib]]), header=TRUE)
  rownames(matrixData[[lib]]) = matrixData[[lib]][,1]
  matrixData[[lib]] = matrixData[[lib]][,-1]
  colnames(matrixData[[lib]]) = sub("\\.", "-", colnames(matrixData[[lib]]))
  print(c( lib, dim(matrixData[[lib]])))
}

# [1] "Patski" "144"    "4493"  
# [1] "EBdiff_TsixStop_d20" "144"                 "404"                
# [1] "EBdiff_F121_d20" "144"             "301"            
# [1] "F121_NPC" "144"      "1054"    
# [1] "TsixStop_NPC" "144"          "2535"   

# [1] "EBdiff_TsixStop_d0" "144"                "795"               
# [1] "EBdiff_F121_d0" "144"            "1479"          
# [1] "F121_d0" "144"     "2203"   
# [1] "TsixStop_d0" "144"         "342"  

# [1] "F121_d0" "144"     "2203"   
# [1] "F121_d3" "144"     "894"    
# [1] "F121_d7" "144"     "732"    
# [1] "F121_d11" "144"      "608"     
# [1] "F121_NPC" "144"      "1054"  

# [1] "TsixStop_d0" "144"         "342"        
# [1] "TsixStop_d3" "144"         "839"        
# [1] "TsixStop_d7" "144"         "822"        
# [1] "TsixStop_d11" "144"          "471"         
# [1] "TsixStop_NPC" "144"          "2535"  

# matrixData[[lib]][1:5,1:5]










##################################################
##################################################
##################################################
# Analyze

# Minimum contact counts per cell
nThresh=1000 

subsupermainOutputDir = file.path(supermainOutputDir, 
                                  paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""))
if ( ! file.exists( subsupermainOutputDir ) ) {
  if ( ! dir.create( subsupermainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) subsupermainOutputDir = supermainOutputDir
}
setwd(subsupermainOutputDir)
getwd()










##################################################
##################################################
##################################################
# Select cells w/ sufficient coverage
# Filter for cells in both libraries that pass filter and are common to both alleles.

# Coverage
matrixCsums = list()
dataMatrix.OG = list()

for(lib in libs) {
  dataMatrix.OG[[lib]] = data.matrix(matrixData[[lib]])
  matrixCsums[[lib]] = colSums(dataMatrix.OG[[lib]], na.rm=TRUE)
  names(matrixCsums[[lib]]) = colnames(dataMatrix.OG[[lib]])
  
  passFilterIdx = matrixCsums[[lib]] >= nThresh
  dataMatrix.OG[[lib]] = dataMatrix.OG[[lib]][,passFilterIdx]
  matrixCsums[[lib]] = matrixCsums[[lib]][passFilterIdx]
}


# 20210310 - to jointly analyze the non-allelic autosomal CDPS of F211 cells
# # 20210226 Just for coverage plots, I think, so skip
# Combine cells from all libraries
matrixCsums[[combined_label]] = c()
dataMatrix.OG[[combined_label]] = c()
cellIDs.combo = c()
cell.libOfOrigin = c()
for(lib in libs) {
  matrixCsums[[combined_label]] = c( matrixCsums[[combined_label]], matrixCsums[[lib]])
  dataMatrix.OG[[combined_label]] = cbind(dataMatrix.OG[[combined_label]], dataMatrix.OG[[lib]])
  # all(colnames(dataMatrix.OG[[lib]]) == names(matrixCsums[[lib]]))
  cellIDs.combo = c(cellIDs.combo, paste(lib, ".", names(matrixCsums[[lib]]), sep=""))
  cell.libOfOrigin = c(cell.libOfOrigin, rep(lib, length(matrixCsums[[lib]])))
}

names(matrixCsums[[combined_label]]) = cellIDs.combo
colnames(dataMatrix.OG[[combined_label]]) = cellIDs.combo

# Don't rename here, b/c complicates date assignment
# names(matrixCsums[[combined_label]]) = cellIDs.combo
# colnames(dataMatrix.OG[[combined_label]]) = cellIDs.combo

whichLibs = c(libs)
# whichLibs = c(libs, combined_label)
# whichLibs = c(combined_label)
# whichLibs = c("EBdiff_F121_d20")










##################################################
##################################################
##################################################
# Coverage level ordering and coloring
# Colors ordered by cell coverage

covCols = list()
covOrd = list()

for(lib in whichLibs) {
  covOrd[[lib]] = order(matrixCsums[[lib]])
  names(covOrd[[lib]]) = names(matrixCsums[[lib]])[ covOrd[[lib]] ]
  covCols[[lib]] = viridis(length(matrixCsums[[lib]]), alpha=0.5)
  covCols[[lib]][covOrd[[lib]]] = viridis(length(matrixCsums[[lib]]), alpha=0.5) 
}

# # check
# head(names(covOrd[[lib]]))
# head(matrixCsums[[lib]][ covOrd[[lib]] ])
# 
# barplot(rep(1,length( covCols[[lib]] )), col=covCols[[lib]] )
# barplot(rep(1,length( covCols[[lib]] )), col=covCols[[lib]][covOrd[[lib]]] )
# 
# barplot(matrixCsums[[lib]], col= covCols[[lib]], las=2, cex.names=0.5 )
# barplot(matrixCsums[[lib]][ covOrd[[lib]] ], col= covCols[[lib]][covOrd[[lib]]], las=2, cex.names=0.5 )










##################################################
##################################################
##################################################
# Preportion of cis reads in contact range

propDecayRate50kbTo8Mb = list()
for(lib in whichLibs) {
  propDecayRate50kbTo8Mb[[lib]] = (log10(dataMatrix.OG[[lib]][which.min(abs(binDist-50000)),]) 
                                   - log10(dataMatrix.OG[[lib]][which.min(abs(binDist-8000000)),])) /
    (log10(8000000)-log10(500000))  
}










##################################################
##################################################
##################################################
# Library labeling

matrixColTPs = list()
matrixColTPcols = list()
timepointThenCovOrd = list()

for (lib in whichLibs) {
  if (lib == combined_label ) {
    matrixColTPs[[lib]] = cell.libOfOrigin
    matrixColTPcols[[lib]] = rep("black", length(matrixColTPs[[lib]]))
    for ( a_sample in unique(cell.libOfOrigin) ) {
      matrixColTPcols[[lib]][matrixColTPs[[lib]]== a_sample] = tpCols[a_sample]
    }
    names(matrixColTPcols[[lib]]) = matrixColTPs[[lib]]
    
    # Order by library then by coverage
    # timepointThenCovOrd[[lib]] = order(cell.libOfOrigin)
    # timepointThenCovOrd[[lib]] = seq(1, length(matrixCsums[[lib]]) )
    timepointThenCovOrd[[lib]] = order(as.numeric(as.character(sub("NPC", "100",
                                                                   sub("d", "", 
                                                                       sub("TsixStop_", "", 
                                                                           sub("F121_", "", 
                                                                               matrixColTPs[[lib]][ covOrd[[lib]] ])))) )))
    names(timepointThenCovOrd[[lib]]) = matrixColTPs[[lib]][ covOrd[[lib]] ][timepointThenCovOrd[[lib]]]
  } else {
    matrixColTPs[[lib]] = rep( lib, length(matrixCsums[[lib]]) )
    matrixColTPcols[[lib]] = rep(tpCols[[lib]], length(matrixColTPs[[lib]]))
    names(matrixColTPcols[[lib]]) = matrixColTPs[[lib]]
    
    # Order by library then by coverage
    # timepointThenCovOrd[[lib]] = order(cell.libOfOrigin)
    timepointThenCovOrd[[lib]] = seq(1, length(matrixCsums[[lib]]) )
    names(timepointThenCovOrd[[lib]]) = matrixColTPs[[lib]][ covOrd[[lib]] ][timepointThenCovOrd[[lib]]]
  }
}










##################################################
##################################################
##################################################
# Coverage plots

covplotsMaindir = file.path(subsupermainOutputDir, "coveragePlots")
if ( ! file.exists( covplotsMaindir ) ) {
  if ( ! dir.create( covplotsMaindir, showWarnings=FALSE, recursive=TRUE )  ) covplotsMaindir = subsupermainOutputDir
}

for(lib in whichLibs) {

  titleText = paste(lib, sep=" ")
  nameText = paste(lib, sep=".")

  # Ranked by timepoint, then coverage, colored by timepoint
  PDFfilename = file.path(covplotsMaindir, paste(projectName, ".coverage.barplots.rankAndColByTimepoint.", nameText, ".pdf", sep=""))
  pdf(PDFfilename, height=6, width=12)
  # # unsorted
  # barplot(matrixCsums[[lib]],
  #         col=matrixColTPcols[[lib]], border=NA, names.arg=NA,
  # barplot(matrixCsums[[lib]][timepointThenCovOrd[[lib]]],
  #         col=matrixColTPcols[[lib]][timepointThenCovOrd[[lib]]], border=NA, names.arg=NA,
  barplot(matrixCsums[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]],
           col=matrixColTPcols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]], border=NA, names.arg=NA,
          xlab=paste(ncol(dataMatrix.OG[[lib]])),
          ylab="total cis coverage",
          main=paste("ranked coverage: ", titleText, sep=""))
          # sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "))
  dev.off()

}










##################################################
##################################################
##################################################
# Preprocess data

# Contact length range to use.
bins2Use="50kb-8Mb" # Nagano et al.

# Also work:
# bins2Use="1kb-80Mb"
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

mainOutputDir = file.path(subsupermainOutputDir, 
                          bins2Use)
if ( ! file.exists( mainOutputDir ) ) {
  if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = subsupermainOutputDir
}
setwd(mainOutputDir)
getwd()

matrixRvar = list()
dataMatrix = list()
dataMatrix.scaled = list()
dataMatrix.prop = list()

for(lib in whichLibs) {
  
  if (bins2Use=="50kb-8Mb") {
    # Only bins from 50kb to 8Mb
    binBoolIdx = binDist >= 50000 & binDist <= 8000000
    sum(binBoolIdx)
    dataMatrix[[lib]] = dataMatrix.OG[[lib]][binBoolIdx, ]
  } else if (bins2Use=="50kb-8Mb+Trans") {
    # Only bins from 50kb to 8Mb
    binBoolIdx = binDist >= 50000 & binDist <= 8000000
    sum(binBoolIdx)
    dataMatrix[[lib]] = dataMatrix.OG[[lib]][binBoolIdx | c(rep(FALSE, length(binBoolIdx)-1), TRUE), ]
  } else if (bins2Use=="1kb-80Mb") {
    binBoolIdx = binDist > 1024 & binDist <= 80000000
    sum(binBoolIdx)
    dataMatrix[[lib]] = dataMatrix.OG[[lib]][binBoolIdx, ]
  } else if (bins2Use=="20kb-Cis") {
    # Only bins from 20kb to 20Mb
    binBoolIdx = binDist >= 20000
    sum(binBoolIdx)
    dataMatrix[[lib]] = dataMatrix.OG[[lib]][binBoolIdx & c(rep(TRUE, length(binBoolIdx)-1), FALSE), ]
  } else if (bins2Use=="allCis") {
    binBoolIdx = rep(TRUE, nrow(dataMatrix.OG[[lib]])-1 )
    dataMatrix[[lib]] = dataMatrix.OG[[lib]][-nrow(dataMatrix.OG[[lib]]), ]
  }
  
  # PCA doesn't work on zero variance rows
  # Unlikely to have zero variance columns
  matrixRvar[[lib]] = apply(dataMatrix[[lib]], 1, var, na.rm=TRUE)
  zeroVarRowIdx = matrixRvar[[lib]] == 0 | is.na(matrixRvar[[lib]])
  dataMatrix[[lib]] = dataMatrix[[lib]][!zeroVarRowIdx,]
  
  # columnwise scaling
  dataMatrix.scaled[[lib]] = scale(dataMatrix[[lib]]) 
  
  # columnwise proportions
  dataMatrix.prop[[lib]] = apply(dataMatrix[[lib]], MARGIN=2, FUN=getProportions)
  # to avoid log10(0)
  dataMatrix.prop[[lib]] = apply(dataMatrix.prop[[lib]], MARGIN=2, FUN=addNZval)
  dataMatrix.prop[[lib]] = log10(dataMatrix.prop[[lib]])
  
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

binDist2Use = binDist[binBoolIdx]










##################################################
##################################################
##################################################
# Heatmaps of contact decay

heatMaindir = file.path(mainOutputDir, "decayMaps")
if ( ! file.exists( heatMaindir ) ) {
  if ( ! dir.create( heatMaindir, showWarnings=FALSE, recursive=TRUE )  ) heatMaindir = mainOutputDir
}
for (lib in whichLibs) {
  print(lib)
  
  for (useScaled in c(FALSE, TRUE)) {
    
    if (!useScaled) {
      matrix2Use = dataMatrix.prop[[lib]]
    } else {
      matrix2Use = dataMatrix.scaled[[lib]]
    }
    
    
    ##################################################
    # Sort cells by timepoint and coverage
    
    # head(matrixCsums[[lib]][ covOrd[[lib]] ])
    # tail(matrixCsums[[lib]][ covOrd[[lib]] ])
    
    matrix2Plot = matrix2Use[ ,covOrd[[lib]][timepointThenCovOrd[[lib]]] ]
    cellLabels = rep("", ncol(matrix2Plot))
    cellLabels[1] = colnames(matrix2Plot)[1]
    cellLabels[length(cellLabels)] = colnames(matrix2Plot)[length(cellLabels)]
    
    TPcolbar = matrixColTPcols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
    covcolbar = covCols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
    rColBars = t(cbind(TPcolbar, covcolbar))
    rownames(rColBars)=c("sample", "coverage" )
    
    titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), sep=".")
    
    minVal = min(matrix2Plot, na.rm=TRUE)
    maxVal = max(matrix2Plot, na.rm=TRUE)

    if(useScaled) {
      cols2Use = heatCols2
      breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
    } else {
      cols2Use = heatColsV
      breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
    }
    if ( ! file.exists( heatOutdir ) ) {
      if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
    }
    
    PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=6, width=12)
    
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

  }
}










##################################################
##################################################
##################################################
# Sorted heatmaps of contact decay

heatMaindir = file.path(mainOutputDir, "decayMaps")
if ( ! file.exists( heatMaindir ) ) {
  if ( ! dir.create( heatMaindir, showWarnings=FALSE, recursive=TRUE )  ) heatMaindir = mainOutputDir
}
for (lib in whichLibs) {
  print(lib)
  
  for (useScaled in c(FALSE, TRUE)) {
    
    if (!useScaled) {
      matrix2Use = dataMatrix.prop[[lib]]
    } else {
      matrix2Use = dataMatrix.scaled[[lib]]
    }
    
    
    ##################################################
    # Sort cells by timepoint and coverage
    
    # head(matrixCsums[[lib]][ covOrd[[lib]] ])
    # tail(matrixCsums[[lib]][ covOrd[[lib]] ])
    
    matrix2Plot = matrix2Use[ ,covOrd[[lib]][timepointThenCovOrd[[lib]]] ]
    cellLabels = rep("", ncol(matrix2Plot))
    cellLabels[1] = colnames(matrix2Plot)[1]
    cellLabels[length(cellLabels)] = colnames(matrix2Plot)[length(cellLabels)]
    
    TPcolbar = matrixColTPcols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
    covcolbar = covCols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]
    rColBars = t(cbind(TPcolbar, covcolbar))
    rownames(rColBars)=c("sample", "coverage" )
    
    titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), " (n = ", dim(matrix2Plot)[2], ")", sep="")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), sep=".")
    
    minVal = min(matrix2Plot, na.rm=TRUE)
    maxVal = max(matrix2Plot, na.rm=TRUE)
    
    if(useScaled) {
      cols2Use = heatCols2
      breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      heatOutdir = file.path(heatMaindir, paste("scaled", sep=""))
    } else {
      cols2Use = heatColsV
      breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      heatOutdir = file.path(heatMaindir, paste("prop", sep=""))
    }
    if ( ! file.exists( heatOutdir ) ) {
      if ( ! dir.create( heatOutdir, showWarnings=FALSE, recursive=TRUE )  ) heatOutdir = heatMaindir
    }
    
    
    ##################################################
    # Sort by PC1 score
    
    PC1scores <- prcomp(t(matrix2Use), scale = TRUE)$x[, 1]
  
    matrix2Plot = matrix2Use[,order(PC1scores, decreasing=TRUE)]
    
    TPcolbar = matrixColTPcols[[lib]][order(PC1scores, decreasing=TRUE)]
    covcolbar = covCols[[lib]][order(PC1scores, decreasing=TRUE)]
    rColBars = t(cbind(TPcolbar, covcolbar))
    rownames(rColBars)=c("sample", "coverage")
    
    titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), " sortedByPC1", " (n = ", dim(matrix2Plot)[2], ")", sep="")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), "sortedByPC1", sep=".")
    
    minVal = min(matrix2Plot, na.rm=TRUE)
    maxVal = max(matrix2Plot, na.rm=TRUE)
    if(useScaled) {
      cols2Use = heatCols2
      breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
    } else {
      cols2Use = heatColsV
      breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
    }
    
    PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=6, width=12)
    
    heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                 Rowv=FALSE, Colv=FALSE,
                 dendrogram="none",
                 scale="none", key=TRUE, density.info="none",
                 trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
    
    
    
    ##################################################
    # Sort by 50kb-80Mb contact decay rate

    matrix2Plot = matrix2Use[,order(propDecayRate50kbTo8Mb[[lib]])]
    
    TPcolbar = matrixColTPcols[[lib]][order(propDecayRate50kbTo8Mb[[lib]])]
    covcolbar = covCols[[lib]][order(propDecayRate50kbTo8Mb[[lib]])]
    rColBars = t(cbind(TPcolbar, covcolbar))
    rownames(rColBars)=c("sample", "coverage")
    
    titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), " sortedByPropDecayRate50kbTo8Mb", " (n = ", dim(matrix2Plot)[2], ")", sep="")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), "sortedByPropDecayRate50kbTo8Mb", sep=".")
    
    minVal = min(matrix2Plot, na.rm=TRUE)
    maxVal = max(matrix2Plot, na.rm=TRUE)
    if(useScaled) {
      cols2Use = heatCols2
      breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
    } else {
      cols2Use = heatColsV
      breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
    }
    
    PDFfilename = file.path(heatOutdir, paste(projectName, ".decaymap.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=6, width=12)
    
    heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                 Rowv=FALSE, Colv=FALSE,
                 dendrogram="none",
                 scale="none", key=TRUE, density.info="none",
                 trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
    
  }
}










##################################################
##################################################
##################################################
# Save matrices of CDPs

matrixMaindir = file.path(mainOutputDir, "decayMatrices")
if ( ! file.exists( matrixMaindir ) ) {
  if ( ! dir.create( matrixMaindir, showWarnings=FALSE, recursive=TRUE )  ) matrixMaindir = mainOutputDir
}
for (lib in whichLibs) {
  print(lib)
  
  for (matrixType in c("counts", "scaled", "prop")) {
    
    if (matrixType == "prop") {
      matrix2Use = dataMatrix.prop[[lib]]
    } else if (matrixType == "scaled") {
      matrix2Use = dataMatrix.scaled[[lib]]
    } else if (matrixType == "counts") {
      matrix2Use = dataMatrix[[lib]]
    }
    
    ##################################################
    # Sort cells by timepoint and coverage
    
    # head(matrixCsums[[lib]][ covOrd[[lib]] ])
    # tail(matrixCsums[[lib]][ covOrd[[lib]] ])
    
    matrix2Save = matrix2Use[ ,covOrd[[lib]][timepointThenCovOrd[[lib]]] ]
    # colnames(matrix2Save) = paste( colnames(matrix2Save),
    #                                names(matrixColTPcols[[lib]][covOrd[[lib]]][timepointThenCovOrd[[lib]]]),
    #                                sep=".")
    rownames(matrix2Save) = binDist2Use

    matrixOutdir = file.path(matrixMaindir, paste(matrixType, sep=""))
    if ( ! file.exists( matrixOutdir ) ) {
      if ( ! dir.create( matrixOutdir, showWarnings=FALSE, recursive=TRUE )  ) matrixOutdir = heatMaindir
    }

    nameText = paste(lib, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), bins2Use, matrixType, sep=".")
    matrixFileName = file.path(matrixOutdir, paste(projectName, ".", nameText, ".tsv", sep=""))
    write.table(x=matrix2Save, file=matrixFileName, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
    try(system(paste("gzip -f ", matrixFileName, sep="")))
    
    
    ##################################################
    # Sort by PC1 score
    # 20200109
    
    PC1scores <- prcomp(t(matrix2Use), scale = TRUE)$x[, 1]
    matrix2Save = matrix2Use[,order(PC1scores, decreasing=TRUE)]

    rownames(matrix2Save) = binDist2Use
    
    matrixOutdir = file.path(matrixMaindir, paste(matrixType, sep=""))
    if ( ! file.exists( matrixOutdir ) ) {
      if ( ! dir.create( matrixOutdir, showWarnings=FALSE, recursive=TRUE )  ) matrixOutdir = heatMaindir
    }
    
    nameText = paste(lib, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), bins2Use, matrixType, "sortedByPC1", sep=".")
    matrixFileName = file.path(matrixOutdir, paste(projectName, ".", nameText, ".tsv", sep=""))
    write.table(x=matrix2Save, file=matrixFileName, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
    try(system(paste("gzip -f ", matrixFileName, sep="")))
    
    
    ##################################################
    # Sort by 50kb-80Mb contact decay rate
    
    matrix2Save = matrix2Use[,order(propDecayRate50kbTo8Mb[[lib]])]

    rownames(matrix2Save) = binDist2Use
    
    matrixOutdir = file.path(matrixMaindir, paste(matrixType, sep=""))
    if ( ! file.exists( matrixOutdir ) ) {
      if ( ! dir.create( matrixOutdir, showWarnings=FALSE, recursive=TRUE )  ) matrixOutdir = heatMaindir
    }
    
    nameText = paste(lib, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), bins2Use, matrixType, "sortedByPropDecayRate50kbTo8Mb", sep=".")
    matrixFileName = file.path(matrixOutdir, paste(projectName, ".", nameText, ".tsv", sep=""))
    write.table(x=matrix2Save, file=matrixFileName, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
    try(system(paste("gzip -f ", matrixFileName, sep="")))
    
  }
}










##################################################
##################################################
##################################################
# Kmeans

kCenters = 4
# kCenters = 8

kclustMaindir = file.path(mainOutputDir, paste( "kclust.k=", kCenters, sep="" ) )
if ( ! file.exists( kclustMaindir ) ) {
  if ( ! dir.create( kclustMaindir, showWarnings=FALSE, recursive=TRUE )  ) kclustMaindir = mainOutputDir
}

# set.seed(1234); kClustColors =  sample(color, kCenters)
# kClustColors = c("#9ac15e",
#                   "#8d50a9",
#                   "#b35948",
#                   "#758a8a")
# kClustColors =  rainbow(kCenters)
# kClustColors = c("#006666",
#                  "#009999",
#                  "#00CCCC",
#                  "#00FFFF")
kClustColors = colorRampPalette(c("#00ffff", "#001e1e"))(kCenters)
names(kClustColors) = paste("group", 1:kCenters)

kClusters.OG = list()
kClusters = list()

for (lib in whichLibs) {
  print(lib)
  
  kClusters.OG[[lib]] = list()
  kClusters[[lib]] = list()
  
  for (useScaled in c(FALSE, TRUE)) {
    
    kClusters.OG[[lib]][[ifelse(useScaled, "scaled", "prop")]] = list()
    kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]] = list()
    
    for (useMethod in c("pearson", "spearman", "euclidean")) {
      
      if (!useScaled) {
        matrix2Use = dataMatrix.prop[[lib]]
        kclustOutdir = file.path(kclustMaindir, paste("prop.", useMethod, sep=""))
      } else {
        matrix2Use = dataMatrix.scaled[[lib]]
        kclustOutdir = file.path(kclustMaindir, paste("scaled.", useMethod, sep=""))
      }
      if ( ! file.exists( kclustOutdir ) ) {
        if ( ! dir.create( kclustOutdir, showWarnings=FALSE, recursive=TRUE )  ) kclustOutdir = kclustMaindir
      }
      
      if (useMethod == "euclidean") {
        theDistMatrix <- dist(t(matrix2Use))
      } else {
        theDistMatrix <- as.dist(1-cor(matrix2Use, method=useMethod))
      }
      
      set.seed(1234)
      kClusters.OG[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]] <- kmeans(theDistMatrix, centers=kCenters, nstart=1)$cluster
      
      
      ##################################################
      # Re-sort Kids by mean PC1 score
      # Doesn't work b/c PC score sign is arbitrary
      
      PC1scores <- prcomp(t(matrix2Use), scale = TRUE)$x[, 1]      

      
      ##################################################
      # Re-sort Kids by mean scaled score across half the region

      halfbins = ceiling(nrow(matrix2Use)/2)
      halfScores <- colSums(matrix2Use[1:halfbins,], na.rm=TRUE)

      meanmeanHalfScores = c()
      for (theKid in 1:kCenters ) {
        # theKid = 1
        meanmeanHalfScores[theKid] = mean(halfScores[kClusters.OG[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]] == theKid], trim=0.025, na.rm=TRUE)
      }
      meanHalfScoreOrder = order(meanmeanHalfScores)

      # Re-order
      kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]] = kClusters.OG[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      for (theKid in 1:kCenters ) {
        # cat( meanHalfScoreOrder[theKid], " --> ", theKid, "\n", sep="")
        kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]][kClusters.OG[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]] == meanHalfScoreOrder[theKid]] = theKid
      }

      # Check
      newmeanmeanHalfScores= c()
      for (theKid in 1:kCenters ) {
        # theKid = 1
        newmeanmeanHalfScores[theKid] = mean(halfScores[kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]] == theKid], trim=0.025, na.rm=TRUE)
      }
      meanmeanHalfScores[meanHalfScoreOrder]
      newmeanmeanHalfScores
      
      
      ##################################################
      # Sort cells by hierarchical clustering
      
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      
      kClustColbar = kClustColbar[order(kIds)]
      # names(kClustColbar) == colnames(matrix2Use)  
      
      # matrix2Plot = log10(dataMatrix[[lib]] + 1)
      # matrix2Plot = matrix2Plot[,order(kIds) ]
      matrix2Plot = matrix2Use[,order(kIds) ]
      # names(kClustColbar) == colnames(matrix2Plot)
      
      TPcolbar = matrixColTPcols[[lib]][ order(kIds) ]
      covcolbar = covCols[[lib]][ order(kIds) ]
      rColBars = t(cbind(TPcolbar, covcolbar, kClustColbar))
      rownames(rColBars)=c("sample", "coverage", paste(kCenters, "clusters") )
      
      titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), useMethod, " hclustSorted", " (n = ", dim(matrix2Plot)[2], ")", sep="")
      nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, "hclustSorted", sep=".")
      
      minVal = min(matrix2Plot, na.rm=TRUE)
      maxVal = max(matrix2Plot, na.rm=TRUE)
      if(useScaled) {
        cols2Use = heatCols2
        breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      } else {
        cols2Use = heatColsV
        breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      }
      
      PDFfilename = file.path(kclustOutdir, paste(projectName, ".kclust.", nameText, ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)

      heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                   Rowv=TRUE, Colv=FALSE,
                   dendrogram="none",
                   scale="none", key=TRUE, density.info="none",
                   trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
      
      
      
      ##################################################
      # Sort cells by 50kb-80Mb contact decay rate
      
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      
      kClustColbar = kClustColbar[order(propDecayRate50kbTo8Mb[[lib]])]
      # names(kClustColbar) == colnames(matrix2Use)  
      
      # matrix2Plot = log10(dataMatrix[[lib]] + 1)
      # matrix2Plot = matrix2Plot[,order(kIds) ]
      matrix2Plot = matrix2Use[,order(propDecayRate50kbTo8Mb[[lib]]) ]
      # names(kClustColbar) == colnames(matrix2Plot)
      
      TPcolbar = matrixColTPcols[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ]
      covcolbar = covCols[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ]
      rColBars = t(cbind(TPcolbar, covcolbar, kClustColbar))
      rownames(rColBars)=c("sample", "coverage", paste(kCenters, "clusters") )
      
      titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), useMethod, " sortedByDecayRate", " (n = ", dim(matrix2Plot)[2], ")", sep="")
      nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, "sortedByDecayRate", sep=".")
      
      minVal = min(matrix2Plot, na.rm=TRUE)
      maxVal = max(matrix2Plot, na.rm=TRUE)
      if(useScaled) {
        cols2Use = heatCols2
        breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      } else {
        cols2Use = heatColsV
        breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      }
      
      PDFfilename = file.path(kclustOutdir, paste(projectName, ".kclust.", nameText, ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)
      
      heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                   Rowv=FALSE, Colv=FALSE,
                   dendrogram="none",
                   scale="none", key=TRUE, density.info="none",
                   trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
      
      
      
      ##################################################
      # Sort cells by cluster
      
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      
      # Order by kIdx then by decay
      # timepointThenDecayOrd[[lib]] = order(as.numeric(as.character(sub("d", "", matrixColTPs[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ]))))
      # names(timepointThenDecayOrd[[lib]]) = matrixColTPs[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ][timepointThenDecayOrd[[lib]]]
      kClustOrder= order(kIds)
      
      kClustColbar = kClustColbar[kClustOrder]
      # names(kClustColbar) == colnames(matrix2Use)  
      
      # matrix2Plot = log10(dataMatrix[[lib]] + 1)
      # matrix2Plot = matrix2Plot[,order(kIds) ]
      matrix2Plot = matrix2Use[,kClustOrder]
      # names(kClustColbar) == colnames(matrix2Plot)
      
      TPcolbar = matrixColTPcols[[lib]][kClustOrder]
      covcolbar = covCols[[lib]][kClustOrder]
      rColBars = t(cbind(TPcolbar, covcolbar, kClustColbar))
      rownames(rColBars)=c("sample", "coverage", paste(kCenters, "clusters") )
      
      titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), useMethod, " sortedByCluster", " (n = ", dim(matrix2Plot)[2], ")", sep="")
      nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, "sortedByCluster", sep=".")
      
      minVal = min(matrix2Plot, na.rm=TRUE)
      maxVal = max(matrix2Plot, na.rm=TRUE)
      if(useScaled) {
        cols2Use = heatCols2
        breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      } else {
        cols2Use = heatColsV
        breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      }
      
      PDFfilename = file.path(kclustOutdir, paste(projectName, ".kclust.", nameText, ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)
      
      heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                   Rowv=FALSE, Colv=FALSE,
                   dendrogram="none",
                   scale="none", key=TRUE, density.info="none",
                   trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
      
      
      
      ##################################################
      # Sort cells by cluster and then PC1scores
      
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      
      # Order by kIdx then by PC1
      kClustOrderThenPC1ord = order(kIds[order(PC1scores)])
      
      kClustColbar = kClustColbar[order(PC1scores)][kClustOrderThenPC1ord]
      # names(kClustColbar) == colnames(matrix2Use)  
      
      # matrix2Plot = log10(dataMatrix[[lib]] + 1)
      # matrix2Plot = matrix2Plot[,order(kIds) ]
      matrix2Plot = matrix2Use[,order(PC1scores)][,kClustOrderThenPC1ord]
      # names(kClustColbar) == colnames(matrix2Plot)
      
      TPcolbar = matrixColTPcols[[lib]][order(PC1scores)][kClustOrderThenPC1ord]
      covcolbar = covCols[[lib]][order(PC1scores)][kClustOrderThenPC1ord]
      rColBars = t(cbind(TPcolbar, covcolbar, kClustColbar))
      rownames(rColBars)=c("sample", "coverage", paste(kCenters, "clusters") )
      
      titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), useMethod, " sortedByClusterThenPC1", " (n = ", dim(matrix2Plot)[2], ")", sep="")
      nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, "sortedByClusterThenPC1", sep=".")
      
      minVal = min(matrix2Plot, na.rm=TRUE)
      maxVal = max(matrix2Plot, na.rm=TRUE)
      if(useScaled) {
        cols2Use = heatCols2
        breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      } else {
        cols2Use = heatColsV
        breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      }
      
      PDFfilename = file.path(kclustOutdir, paste(projectName, ".kclust.", nameText, ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)
      
      heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                   Rowv=FALSE, Colv=FALSE,
                   dendrogram="none",
                   scale="none", key=TRUE, density.info="none",
                   trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
      
      
      
      ##################################################
      # Sort cells by cluster and then 50kb-80Mb contact decay rate 
      
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      
      # Order by kIdx then by decay
      # timepointThenDecayOrd[[lib]] = order(as.numeric(as.character(sub("d", "", matrixColTPs[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ]))))
      # names(timepointThenDecayOrd[[lib]]) = matrixColTPs[[lib]][ order(propDecayRate50kbTo8Mb[[lib]]) ][timepointThenDecayOrd[[lib]]]
      kClustOrderThenDecayOrd = order(kIds[order(propDecayRate50kbTo8Mb[[lib]])])
      
      # # barplot(propDecayRate50kbTo8Mb[[lib]][order(propDecayRate50kbTo8Mb[[lib]])][timepointThenDecayOrd[[lib]]],
      # #         col=matrixColTPcols[[lib]][order(propDecayRate50kbTo8Mb[[lib]])][timepointThenDecayOrd[[lib]]] )
      # barplot(propDecayRate50kbTo8Mb[[lib]][order(propDecayRate50kbTo8Mb[[lib]])][kClustOrderThenDecayOrd],
      #         col=kClustColbar[order(propDecayRate50kbTo8Mb[[lib]])][kClustOrderThenDecayOrd] )
      
      kClustColbar = kClustColbar[order(propDecayRate50kbTo8Mb[[lib]])][kClustOrderThenDecayOrd]
      # names(kClustColbar) == colnames(matrix2Use)  
      
      # matrix2Plot = log10(dataMatrix[[lib]] + 1)
      # matrix2Plot = matrix2Plot[,order(kIds) ]
      matrix2Plot = matrix2Use[,order(propDecayRate50kbTo8Mb[[lib]])][,kClustOrderThenDecayOrd]
      # names(kClustColbar) == colnames(matrix2Plot)
      
      TPcolbar = matrixColTPcols[[lib]][order(propDecayRate50kbTo8Mb[[lib]])][kClustOrderThenDecayOrd]
      covcolbar = covCols[[lib]][order(propDecayRate50kbTo8Mb[[lib]])][kClustOrderThenDecayOrd]
      rColBars = t(cbind(TPcolbar, covcolbar, kClustColbar))
      rownames(rColBars)=c("sample", "coverage", paste(kCenters, "clusters") )
      
      titleText = paste(lib, " ", ifelse(useScaled, "scaled", "prop"), useMethod, " sortedByClusterThenDecayRate", " (n = ", dim(matrix2Plot)[2], ")", sep="")
      nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, "sortedByClusterThenDecayRate", sep=".")
      
      minVal = min(matrix2Plot, na.rm=TRUE)
      maxVal = max(matrix2Plot, na.rm=TRUE)
      if(useScaled) {
        cols2Use = heatCols2
        breaks2Use = seq(minVal, maxVal, (maxVal-minVal)/length(cols2Use))
      } else {
        cols2Use = heatColsV
        breaks2Use = seq(-3, maxVal, (maxVal-(-3))/length(cols2Use))
      }
      
      PDFfilename = file.path(kclustOutdir, paste(projectName, ".kclust.", nameText, ".pdf", sep=""))
      pdf(PDFfilename, height=6, width=12)
      
      heatmap.3.GB(t(matrix2Plot), col=cols2Use, breaks=breaks2Use,
                   Rowv=FALSE, Colv=FALSE,
                   dendrogram="none",
                   scale="none", key=TRUE, density.info="none",
                   trace="none", cexCol=0.95, labCol=binDistDesc, labRow=NA, margins=c(5,0),
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
      
    }
  }
}










##################################################
##################################################
##################################################
# Decay curves

decayCurveMaindir = file.path(mainOutputDir, "decayCurves")
if ( ! file.exists( decayCurveMaindir ) ) {
  if ( ! dir.create( decayCurveMaindir, showWarnings=FALSE, recursive=TRUE )  ) decayCurveMaindir = mainOutputDir
}

lwd2Use = 0.5
alpha2Use=1

for (lib in whichLibs) {
  print(lib)
  
  for (useScaled in c(FALSE, TRUE)) {
    
    if (!useScaled) {
      matrix2Use = dataMatrix.prop[[lib]]
    } else {
      matrix2Use = dataMatrix.scaled[[lib]]
    }
    
    ######################
    # Color by Sample Type
    
    decayOutdir = file.path(decayCurveMaindir, "BySampleType", ifelse(useScaled, "scaled", "prop"))
    if ( ! file.exists( decayOutdir ) ) {
      if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
    }
    
    titleText = paste(lib, ifelse(useScaled, "scaled", "prop"), "colBySampleType", sep=" ")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), "colBySampleType", sep=".")
    
    PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=12, width=24)
    
    par(mfrow=c(2,2))
    
    for (theTP in names(tpCols)) {
      nCells = sum(matrixColTPs[[lib]]==theTP)
      if (nCells==0) next 
      # alaph2Use =  1/nCells
      col2Use = setAlpha(tpCols[theTP], alpha2Use)
      for (iCell in 1:nCells) {
        if (useScaled) {
          y = matrix2Use[, matrixColTPs[[lib]]==theTP][,iCell]
          ylabText = "Scaled Cis counts"
          ylims = c(-3, +3)
          yat = seq(-3, +3, 0.5)
          ylabels = yat
        } else {
          y = matrix2Use[, matrixColTPs[[lib]]==theTP][,iCell]
          ylabText = "Proportion of Cis contact (log10)"
          ylims = c(-3, -1.5)
          yat = seq(-3, -1.5 , 0.5)
          ylabels = signif(10^yat, 2)
        }
        if (iCell==1){
          plot( y, type="l", col=col2Use, lwd=lwd2Use,
                main=theTP, xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                ylim=ylims, xaxt="n", yaxt="n")
          axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
          axis(2, at=yat, labels=ylabels)
          if (theTP == names(tpCols)[1]) {
            legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
          }
        } else {
          lines(y, col=col2Use, lwd=lwd2Use)
        }
      }
    }
    title(main=titleText, outer=TRUE, cex.main=1.5, line=-2)
    
    dev.off()
    
    
    ######################
    # Color by Klust group
    
    kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
    # kIds = kClusters[[lib]][["scaled"]][["euclidean"]]
    kClustColbar =  rep("black", length(kIds))
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    
    kIdColTable = cbind(names(kIds), kIds, kClustColbar)
    colnames(kIdColTable) = c("cellID", "group", "color")
    
    decayOutdir = file.path(decayCurveMaindir, paste( "ByClusterGroup", ".k=", kCenters, sep="" ), ifelse(useScaled, "scaled", "prop"))
    if ( ! file.exists( decayOutdir ) ) {
      if ( ! dir.create( decayOutdir, showWarnings=FALSE, recursive=TRUE )  ) decayOutdir = decayCurveMaindir
    }
    
    titleText = paste(lib, ifelse(useScaled, "scaled", "prop"), "colByClusterGroup", sep=" ")
    nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), "colByClusterGroup", sep=".")
    
    PDFfilename = file.path(decayOutdir, paste(projectName, ".decayCurve.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=12, width=24)
    
    par(mfrow=c(2,2))
    
    for (theTP in names(tpCols)) {
      nCells = sum(matrixColTPs[[lib]]==theTP)
      if (nCells==0) next 
      # alaph2Use =  1/nCells
      for (iCell in 1:nCells) {
        col2Use = setAlpha(kIdColTable[ kIdColTable[,"cellID"] == colnames(matrix2Use[, matrixColTPs[[lib]]==theTP])[iCell], 
                                        "color" ], alpha2Use)
        if (useScaled) {
          y = matrix2Use[, matrixColTPs[[lib]]==theTP][,iCell]
          ylabText = "Scaled Cis counts"
          ylims = c(-3, +3)
          yat = seq(-3, +3, 0.5)
          ylabels = yat
        } else {
          y = matrix2Use[, matrixColTPs[[lib]]==theTP][,iCell]
          ylabText = "Proportion of Cis contacts (log10)"
          ylims = c(-3, -1.5)
          yat = seq(-3, -1.5 , 0.5)
          ylabels = signif(10^yat, 2)
        }
        if (iCell==1){
          plot( y, type="l", col=col2Use, lwd=lwd2Use,
                main=theTP, xlab="contact distance\n(bins of 2^x in increments of 0.125)", ylab=ylabText,
                ylim=ylims, xaxt="n", yaxt="n")
          axis(1, at=which(!is.na(binDistDesc)), labels=binDistDesc)
          axis(2, at=yat, labels=ylabels)
          if (theTP == names(tpCols)[1]) {
            legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
          }
        } else {
          lines(y, col=col2Use, lwd=lwd2Use)
        }
      }
    }
    title(main=titleText, outer=TRUE, cex.main=1.5, line=-2)
    
    dev.off()
  }
  
}










##################################################
##################################################
##################################################
# 20210310
# samples per cluster for combined analysis 

for (lib in whichLibs) {
  print(lib)
  if( lib == combined_label) {
    for (useScaled in c(TRUE)) {
      for (useMethod in c("spearman")) {
      
        cluster_df = as.data.frame( cbind(  
          sub("\\..*", "", names(kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]])),
          kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]]
        ) )
        
        colnames(cluster_df) = c("sample", "cluster_group")
        # Set factor order:
        cluster_df$sample <- factor(cluster_df$sample , levels = libs)
        cluster_df$cluster_group <- factor(cluster_df$cluster_group , levels = 1:kCenters)
        
        print(table(cluster_df))
      }
    }
  } 
}

# cluster_group
# sample        1    2    3    4
# F121_d0    89   32  114 1968
# F121_d3    82  123  291  398
# F121_d7    56  113  243  320
# F121_d11   26   65  190  327
# F121_NPC   28   24  222  780











##################################################
##################################################
##################################################
# Save clusters

clusteringOutputDir = file.path("../data/autosomal_non-allelic_cell_cycle_cluster_ids")

nThresh=1000 # min total contact counts per cell

# bins2Use="50kb-8Mb"
# # bins2Use="1kb-80Mb"
# 
# useScaled=TRUE
# # useMethod="euclidean"
# useMethod="spearman"
# 
# kCenters = 4
# # kCenters = 8
# 
# autosomalCDPs=TRUE

kclustFiledir = file.path(clusteringOutputDir, 
	paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), 
	bins2Use, paste( "kclust", ".k=", kCenters, sep="" ), "clusterIDs")
if ( ! file.exists( kclustFiledir ) ) {
  if ( ! dir.create( kclustFiledir, showWarnings=FALSE, recursive=TRUE )  ) kclustFiledir = kclustMaindir
}


for (lib in whichLibs) {
  print(lib)
  for (useScaled in c(FALSE, TRUE)) {
    for (useMethod in c("pearson", "spearman", "euclidean")) {

      nameText = paste(lib, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""),sep=""), bins2Use, 
                       ifelse(useScaled, "scaled", "prop"), 
                       paste("kclust.dist", useMethod, sep=""),
                       paste("kCenters", kCenters, sep=""),
                       sep=".")
      
      kclustFileName = file.path(kclustFiledir, paste(projectName, ".", nameText, ".tsv", sep=""))
      write.table(x=kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[useMethod]], 
                  file=kclustFileName, sep="\t",
                  row.names=TRUE, col.names=FALSE)
    }
  }
}










# ##################################################
# ##################################################
# ##################################################
# # hclust analysis
# 
# hclustMaindir = file.path(mainOutputDir, "hclust")
# if ( ! file.exists( hclustMaindir ) ) {
#   if ( ! dir.create( hclustMaindir, showWarnings=FALSE, recursive=TRUE )  ) hclustMaindir = mainOutputDir
# }
# 
# for (lib in whichLibs) {
#   print(lib)
#   
#   for (useScaled in c(FALSE, TRUE)) {
#     
#     for (useMethod in c("pearson", "spearman")) {
#       
#       if (!useScaled) {
#         prcorrRes <- cor(dataMatrix.prop[[lib]], method=useMethod)
#         hclustOutdir = file.path(hclustMaindir, paste("prop.", useMethod, sep=""))
#       } else {
#         prcorrRes <- cor(dataMatrix.scaled[[lib]], method=useMethod)
#         hclustOutdir = file.path(hclustMaindir, paste("scaled.", useMethod, sep=""))
#       }  
#       if ( ! file.exists( hclustOutdir ) ) {
#         if ( ! dir.create( hclustOutdir, showWarnings=FALSE, recursive=TRUE )  ) hclustOutdir = hclustMaindir
#       }
#       
#       theData = prcorrRes
#       
#       # heatCols1 = heat.colors(256)
#       heatCols2 <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
#                                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
#                                           "#4393C3", "#2166AC", "#053061"))(256))
#       
#       covcolbar = covCols[[lib]]
#       TPcolbar = matrixColTPcols[[lib]]
#       
#       rColBars = rbind(TPcolbar)
#       
#       # Kmeans clusters
#       for (distMethod in c("pearson", "spearman", "euclidean")) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][[distMethod]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         kClustColbar = kClustColbar
#         #  names(kClustColbar) == colnames(theData)
#         rColBars = rbind( rColBars, kClustColbar)
#       }
#       rownames(rColBars)=c("sample", 
#                            paste( kCenters, " ", c("pearson", "spearman", "euclidean"), " clusters", sep="") )
#       # rColBars = rColBars[,ncol(rColBars):1]
#       
#       titleText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, sep=" ")
#       nameText = paste(lib, ifelse(useScaled, "scaled", "prop"), useMethod, sep=".")
#       
#       PDFfilename = file.path(hclustOutdir, paste(projectName, ".hclust.", nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=12, width=12)
#       
#       heatmap.3.GB(theData, col=heatCols2, 
#                    distfun = function(x) as.dist(1-x), dendrogram="row",
#                    scale="none", key=TRUE, keysize=0.5, density.info="none",
#                    trace="none",cexCol=1.2, labRow=NA, labCol=NA, margins=c(0,0),
#                    RowSideColors=covcolbar, # RowSideColorsSize=10,
#                    ColSideColors=t(rColBars), ColSideColorsSize=10, # If used as RowSideColors then flipped
#                    symm=TRUE,symkey=TRUE,symbreaks=TRUE, 
#                    KeyValueName=titleText,
#                    lmat=rbind(c(0,0,6,0),c(0,0,2,0),c(4,1,3,0),c(0,0,5,0)), 
#                    lhei=c(0.5, 0.8, 5, 0.2),
#                    lwid=c(1.5, 0.2, 5, 0.2))
#       
#       dev.off()
#       
#       if (FALSE) {
#         ##################################################
#         # Testing:
#         
#         heatmap(theData, Rowv=NA, Colv=NA, col=heatCols1, scale="none") 
#         heatmap(theData, Rowv=NA, Colv=NA, col=heatCols2, scale="none") 
#         
#         corrplot(theData, col=heatCols2,  order="hclust",
#                  method="square", type = "upper", tl.col = "black", tl.srt = 45, tl.cex=0.5 )
#         
#         # require(gplots)
#         heatmap.2(theData, col=heatCols2, 
#                   distfun = function(x) as.dist(1-x), dendrogram="both",
#                   scale="none", key=TRUE, keysize=0.5, density.info="none",
#                   trace="none", cexCol=1.2, labRow=NA, labCol=NA, 
#                   RowSideColors=TPcolbar,
#                   ColSideColors=covcolbar,
#                   symm=TRUE,symkey=TRUE,symbreaks=TRUE,  
#                   lmat=rbind(c(6,0,5,0),c(0,0,2,0),c(4,1,3,0),c(0,0,0,0)), 
#                   lhei=c(1.25,0.2,3.5,0),
#                   lwid=c(1.25,0.2,3.5,0))
#         
#         hclust(d=as.dist(1-theData), method="average")
#       } # FALSE
#       
#     }
#   } 
# }










##################################################
##################################################
##################################################
# PCA analysis
#
# install.packages("factoextra")
# require(factoextra)

PCAmaindir = file.path(mainOutputDir, "PCA")
if ( ! file.exists( PCAmaindir ) ) {
  if ( ! dir.create( PCAmaindir, showWarnings=FALSE, recursive=TRUE )  ) PCAmaindir = mainOutputDir
}

PCAdims = 3

pchVal=20
pchcexVal=1.25

for (lib in whichLibs) {
  print(lib)
  
  for (useScaled in c(FALSE, TRUE)) {
    
    if (!useScaled) {
      PCAoutdir = file.path(PCAmaindir, "prop")
      titleText = paste(lib, "prop", sep=" ")
      nameText = paste(lib, "prop", sep=".")
      prcompRes <- prcomp(t(dataMatrix.prop[[lib]]), scale = TRUE)
    } else {
      PCAoutdir = file.path(PCAmaindir, "scaled")
      titleText = paste(lib, "scaled", sep=" ")
      nameText = paste(lib, "scaled", sep=".")
      prcompRes <- prcomp(t(dataMatrix.scaled[[lib]]), scale = TRUE)
    }  
    
    if ( ! file.exists( PCAoutdir ) ) {
      if ( ! dir.create( PCAoutdir, showWarnings=FALSE, recursive=TRUE )  ) PCAoutdir = PCAmaindir
    }
    
    length(prcompRes$x[,1])
    
    # https://stackoverflow.com/questions/29053223/plots-not-working-in-for-loop?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PCA.skree.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=6, width=12)
    fviz_plot = fviz_eig(prcompRes, main=lib)
    print(fviz_plot)
    dev.off()
    
    # fviz_pca_var(prcompRes,
    #              col.var = "contrib", # Color by contributions to the PC
    #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    #              repel = TRUE     # Avoid text overlapping
    # )
    
    
    # Plot indexed dims: color and order by coverage
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PCA.DimsVsIndices.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=4, width=4*PCAdims)
    par(mfrow=c(1,PCAdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDim in seq(1,PCAdims)) {
      PCAaLabel = paste("PC", theDim, sep="")
      plot(prcompRes$x[, PCAaLabel][ covOrd[[lib]] ], 
           col=covCols[[lib]][ covOrd[[lib]] ], 
           pch=pchVal, cex=pchcexVal, 
           xlab="index", ylab=paste("dim", theDim), main=NA)
      # text(rtsne_out$Y, labels=rownames(mydata))
    }
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE)
    title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
    dev.off()
    
    
    # Plot indexed dims: ranked by timepoint then coverage; colored by timepoint
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PCA.DimsVsIndices.rankAndColByTimepoint.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=4, width=4*PCAdims) 
    par(mfrow=c(1,PCAdims),
        mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDim in seq(1,PCAdims)) { 
      PCAaLabel = paste("PC", theDim, sep="") 
      plot(prcompRes$x[,  PCAaLabel][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
           col=matrixColTPcols[[lib]][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
           pch=pchVal, cex=pchcexVal, 
           xlab="index", ylab=paste("dim", theDim), main=NA)
      #  text(rtsne_out$Y, labels=rownames(mydata)) 
      if(theDim==1) {
        legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
      }
    } 
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE) 
    title(sub=paste("Ranked by timepoint then coverage; colored by timepoint"), line=1, outer=TRUE) 
    dev.off()
    
    
    # Plot indexed dims: Ranked by timepoint then coverage; colored by euclidean
    # Kmeans clusters
    # for (distMethod in c("pearson", "spearman",)) {
    kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
    kClustColbar =  rep("black", length(kIds))
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    # names(kClustColbar) == colnames(theData)
    # names(kClustColbar) == names(prcompRes$x[timepointThenCovOrd[[lib]], PCAaLabel])
    
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PCA.DimsVsIndices.rankAndColByEuclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=4, width=4*PCAdims) 
    par(mfrow=c(1,PCAdims),
        mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDim in seq(1,PCAdims)) { 
      PCAaLabel = paste("PC", theDim, sep="") 
      plot(prcompRes$x[, PCAaLabel][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
           col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
           pch=pchVal, cex=pchcexVal, 
           xlab="index", ylab=paste("dim", theDim), main=NA)
      #  text(rtsne_out$Y, labels=rownames(mydata)) 
      if(theDim==1) {
        legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
      }
    } 
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE) 
    title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Euclidean distance", paste( " k=", kCenters, sep="" )), line=1, outer=TRUE) 
    dev.off()

    
    # Plot indexed dims: Ranked by timepoint then coverage; colored by spearman
    # Kmeans clusters
    # for (distMethod in c("pearson", "spearman",)) {
    kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
    kClustColbar =  rep("black", length(kIds))
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    # names(kClustColbar) == colnames(theData)
    # names(kClustColbar) == names(prcompRes$x[timepointThenCovOrd[[lib]], PCAaLabel])
    
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PCA.DimsVsIndices.rankAndColBySpearclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=4, width=4*PCAdims) 
    par(mfrow=c(1,PCAdims),
        mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDim in seq(1,PCAdims)) { 
      PCAaLabel = paste("PC", theDim, sep="") 
      plot(prcompRes$x[, PCAaLabel][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
           col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
           pch=pchVal, cex=pchcexVal, 
           xlab="index", ylab=paste("dim", theDim), main=NA)
      #  text(rtsne_out$Y, labels=rownames(mydata)) 
      if(theDim==1) {
        legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
      }
    } 
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE) 
    title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Spearman distance", paste( "k=", kCenters, sep="" )), line=1, outer=TRUE) 
    dev.off()
    
    
    # Plot dim vs dim: color by coverage       
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PC.DimVsDim.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=3*(PCAdims-1), width=3*(PCAdims-1))
    par(mfrow=c(PCAdims-1, PCAdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDimA in seq(1,PCAdims-1)) {
      for (theDimB in seq(1,PCAdims)) {
        if (theDimB == theDimA) next
        if (theDimB > theDimA) {
          PCAaLabel = paste("PC", theDimA, sep="")
          PCAbLabel = paste("PC", theDimB, sep="")
          plot(prcompRes$x[, PCAaLabel], 
               prcompRes$x[, PCAbLabel], 
               col=covCols[[lib]], 
               pch=pchVal, cex=pchcexVal, 
               xlab=PCAaLabel,
               ylab=PCAbLabel,
               main=NA)
        } else {
          plot.new()
        }
      }
    }
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE)
    title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
    dev.off()
    
    
    # Plot dim vs dim: color by timepoint
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PC.DimVsDim.colByTimepoint.", nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=3*(PCAdims-1), width=3*(PCAdims-1))
    par(mfrow=c(PCAdims-1, PCAdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDimA in seq(1,PCAdims-1)) {
      for (theDimB in seq(1,PCAdims)) {
        if (theDimB == theDimA) next
        if (theDimB > theDimA) {
          PCAaLabel = paste("PC", theDimA, sep="")
          PCAbLabel = paste("PC", theDimB, sep="")
          plot(prcompRes$x[, PCAaLabel],
               prcompRes$x[, PCAbLabel],
               col=matrixColTPcols[[lib]], 
               pch=pchVal, cex=pchcexVal,
               xlab=PCAaLabel,
               ylab=PCAbLabel,
               main=NA)
        } else {
          plot.new()
          legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
        }
      }
    }
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE)
    # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
    dev.off()
    
    
    # Plot dim vs dim: color by euclidean clusters
    # Kmeans clusters
    # for (distMethod in c("pearson", "spearman",)) {
    kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
    kClustColbar =  rep("black", length(kIds))
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    # names(kClustColbar) == colnames(theData)
    # names(kClustColbar) == names(prcompRes$x[timepointThenCovOrd[[lib]], PCAaLabel])
    
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PC.DimVsDim.colByEuclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=3*(PCAdims-1), width=3*(PCAdims-1))
    par(mfrow=c(PCAdims-1, PCAdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDimA in seq(1,PCAdims-1)) {
      for (theDimB in seq(1,PCAdims)) {
        if (theDimB == theDimA) next
        if (theDimB > theDimA) {
          PCAaLabel = paste("PC", theDimA, sep="")
          PCAbLabel = paste("PC", theDimB, sep="")
          plot(prcompRes$x[, PCAaLabel],
               prcompRes$x[, PCAbLabel],
               col=kClustColbar, 
               pch=pchVal, cex=pchcexVal,
               xlab=PCAaLabel,
               ylab=PCAbLabel,
               main=NA)
        } else {
          plot.new()
          legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
          # #  require(scatterplot3d)
          # scatterplot3d(x=prcompRes$x[timepointThenCovOrd[[lib]], "PC1"], 
          #               y=prcompRes$x[timepointThenCovOrd[[lib]], "PC2"], 
          #               z=prcompRes$x[timepointThenCovOrd[[lib]], "PC3"], 
          #               xlab="PC1", ylab="PC2", zlab="PC3",
          #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
          #               axis=TRUE, grid=TRUE, box=TRUE,
          #               main=NA, type="p", angle=45) 
        }
      }
    }
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE)
    # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
    dev.off()
    
    
    # Plot dim vs dim: color by Spearman clusters
    # Kmeans clusters
    # for (distMethod in c("pearson", "spearman",)) {
    kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
    kClustColbar =  rep("black", length(kIds))
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    # names(kClustColbar) == colnames(theData)
    # names(kClustColbar) == names(prcompRes$x[timepointThenCovOrd[[lib]], PCAaLabel])
    
    PDFfilename = file.path(PCAoutdir, paste(projectName, ".PC.DimVsDim.colBySpearclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=3*(PCAdims-1), width=3*(PCAdims-1))
    par(mfrow=c(PCAdims-1, PCAdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
    for (theDimA in seq(1,PCAdims-1)) {
      for (theDimB in seq(1,PCAdims)) {
        if (theDimB == theDimA) next
        if (theDimB > theDimA) {
          PCAaLabel = paste("PC", theDimA, sep="")
          PCAbLabel = paste("PC", theDimB, sep="")
          plot(prcompRes$x[, PCAaLabel],
               prcompRes$x[, PCAbLabel],
               col=kClustColbar, 
               pch=pchVal, cex=pchcexVal,
               xlab=PCAaLabel,
               ylab=PCAbLabel,
               main=NA)
        } else {
          plot.new()
          legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
          # #  require(scatterplot3d)
          # scatterplot3d(x=prcompRes$x[timepointThenCovOrd[[lib]], "PC1"], 
          #               y=prcompRes$x[timepointThenCovOrd[[lib]], "PC2"], 
          #               z=prcompRes$x[timepointThenCovOrd[[lib]], "PC3"], 
          #               xlab="PC1", ylab="PC2", zlab="PC3",
          #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
          #               axis=TRUE, grid=TRUE, box=TRUE,
          #               main=NA, type="p", angle=45) 
        }
      }
    }
    title(main=paste("PCA ", titleText, sep=""), line=1, outer=TRUE)
    # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
    dev.off()
    

    if (plot3Dscatters) {
      # Plot 3D scattrplot: color by euclidean clusters
      PCA3dir = file.path(PCAoutdir, "3Dscatters", lib, "pngs")
      if ( ! file.exists( PCA3dir ) ) {
        if ( ! dir.create( PCA3dir, showWarnings=FALSE, recursive=TRUE )  ) PCA3dir = PCAoutdir
      }
      
      # Kmeans clusters
      # for (distMethod in c("pearson", "spearman",)) {
      kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
      kClustColbar =  rep("black", length(kIds))
      kClustColbar = kClustColors[kIds]
      names(kClustColbar) = names(kIds)
      # names(kClustColbar) == colnames(theData)
      # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
      
      someAngles = seq(5, 360, 5)
      for (i in 1:length(someAngles)) {
        PNGfilename = file.path(PCA3dir, paste(sprintf("%03d", i), ".png", sep=""))
        png(PNGfilename, height=8, width=8, units="in", res=150)
        
        # require(scatterplot3d)
        scatterplot3d(x=prcompRes$x[, "PC1"], 
                      y=prcompRes$x[, "PC2"], 
                      z=prcompRes$x[, "PC3"], 
                      xlab="PC1", ylab="PC2", zlab="PC3", type="p", 
                      color=kClustColbar, 
                      pch=pchVal, cex.symbols=pchcexVal*1.25,
                      axis=TRUE, grid=TRUE, box=FALSE,
                      main=someAngles[i], angle=someAngles[i]) 
        
        dev.off()
      }
    }
    
  }
}










# ##################################################
# ##################################################
# ##################################################
# # MDS
# 
# MDSdims = 3
# 
# MDSmaindir = file.path(mainOutputDir, "MDS" )
# if ( ! file.exists( MDSmaindir ) ) {
#   if ( ! dir.create( MDSmaindir, showWarnings=FALSE, recursive=TRUE )  ) MDSmaindir = mainOutputDir
# }
# 
# pchVal=20
# pchcexVal=1.25
# 
# for (lib in whichLibs) {
#   print(lib)
#   
#   for (useScaled in c(FALSE, TRUE)) {
#     
#     for (useEuclidean in c(FALSE, TRUE)) {
#       
#       if (!useScaled) {
#         MDSoutdir = file.path(MDSmaindir, "prop")
#         titleText = paste(lib, "prop", sep=" ")
#         nameText = paste(lib, "prop", sep=".")
#         if (useEuclidean) {
#           MDSoutdir = file.path(MDSoutdir, "euclidean")
#           titleText = paste(titleText, "euclidean", sep=" ")
#           nameText = paste(nameText, "euclidean", sep=".")
#           dmDist <- dist(t(dataMatrix.prop[[lib]])) # euclidean distances between the rows = scaler products of differences
#         } else {
#           MDSoutdir = file.path(MDSoutdir, "correlation")
#           titleText = paste(titleText, "correlation", sep=" ")
#           nameText = paste(nameText, "correlation", sep=".")
#           dmDist <- as.dist(1-cor(dataMatrix.prop[[lib]]))
#         }
#       } else {
#         MDSoutdir = file.path(MDSmaindir, "scaled")
#         titleText = paste(lib, "scaled", sep=" ")
#         nameText = paste(lib, "scaled", sep=".")
#         if (useEuclidean) {
#           MDSoutdir = file.path(MDSoutdir, "euclidean")
#           titleText = paste(titleText, "euclidean", sep=" ")
#           nameText = paste(nameText, "euclidean", sep=".")
#           dmDist <- dist(t(dataMatrix.scaled[[lib]])) # euclidean distances between the rows = scaler products of differences
#         } else {
#           MDSoutdir = file.path(MDSoutdir, "correlation")
#           titleText = paste(titleText, "correlation", sep=" ")
#           nameText = paste(nameText, "correlation", sep=".")
#           dmDist <- as.dist(1-cor(dataMatrix.scaled[[lib]]))
#         }
#       }
#       
#       if ( ! file.exists( MDSoutdir ) ) {
#         if ( ! dir.create( MDSoutdir, showWarnings=FALSE, recursive=TRUE )  ) MDSoutdir = MDSmaindir
#       }
#       
#       ##################################################
#       # MDS nD
#       set.seed(1234)
#       MDSfit <- cmdscale(dmDist, eig=TRUE, k=MDSdims) # k is the number of dim
#       
#       
#       # Plot indexed dims - ranked and color by coverage
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimsVsIndices.", nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=4, width=4*MDSdims)
#       par(mfrow=c(1,MDSdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDim in seq(1,MDSdims)) {
#         plot(MDSfit$points[,theDim][ covOrd[[lib]] ],
#              col=covCols[[lib]][ covOrd[[lib]] ],
#              pch=pchVal, cex=pchcexVal,
#              xlab="index", ylab=paste("dim", theDim), main=NA)
#         # text(rtsne_out$Y, labels=rownames(mydata))
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot indexed dims - ranked and color by timepoint
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimsVsIndices.rankAndColByTimepoint.", nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=4, width=4*MDSdims)
#       par(mfrow=c(1,MDSdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDim in seq(1,MDSdims)) {
#         plot(MDSfit$points[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              col=matrixColTPcols[[lib]][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              pch=pchVal, cex=pchcexVal,
#              xlab="index", ylab=paste("dim", theDim), main=NA)
#         # text(rtsne_out$Y, labels=rownames(mydata))
#         if(theDim==1){
#           legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       title(sub=paste("Ranked by timepoint then coverage; colored by timepoint"), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot indexed dims - Ranked by timepoint then coverage; colored by k-means cluster by Euclidean distance
#       # Kmeans clusters
#       # for (distMethod in c("pearson", "spearman",)) {
#       kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#       kClustColbar =  rep("black", length(kIds))
#       kClustColbar = kClustColors[kIds]
#       names(kClustColbar) = names(kIds)
#       # names(kClustColbar) == colnames(theData)
#       # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#       
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimsVsIndices.rankAndColByEuclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=4, width=4*MDSdims)
#       par(mfrow=c(1,MDSdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDim in seq(1,MDSdims)) {
#         plot(MDSfit$points[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              pch=pchVal, cex=pchcexVal,
#              xlab="index", ylab=paste("dim", theDim), main=NA)
#         # text(rtsne_out$Y, labels=rownames(mydata))
#         if(theDim==1){
#           legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Euclidean distance", paste( "k=", kCenters, sep="" )), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot indexed dims - Ranked by timepoint then coverage; colored by k-means cluster by Spearman distance
#       # Kmeans clusters
#       # for (distMethod in c("pearson", "spearman",)) {
#       kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
#       kClustColbar =  rep("black", length(kIds))
#       kClustColbar = kClustColors[kIds]
#       names(kClustColbar) = names(kIds)
#       # names(kClustColbar) == colnames(theData)
#       # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#       
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimsVsIndices.rankAndColBySpearclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=4, width=4*MDSdims)
#       par(mfrow=c(1,MDSdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDim in seq(1,MDSdims)) {
#         plot(MDSfit$points[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ],
#              pch=pchVal, cex=pchcexVal,
#              xlab="index", ylab=paste("dim", theDim), main=NA)
#         # text(rtsne_out$Y, labels=rownames(mydata))
#         if(theDim==1){
#           legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Spearman distance", paste( "k=", kCenters, sep="" )), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot dim vs dim - color by coverage
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimVsDim.", nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=3*(MDSdims-1), width=3*(MDSdims-1))
#       par(mfrow=c(MDSdims-1, MDSdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDimA in seq(1,MDSdims-1)) {
#         for (theDimB in seq(1,MDSdims)) {
#           if (theDimB == theDimA) next
#           if (theDimB > theDimA) {
#             plot(MDSfit$points[,theDimA],
#                  MDSfit$points[,theDimB],
#                  col=covCols[[lib]],
#                  pch=pchVal, cex=pchcexVal,
#                  xlab=paste("dim", theDimA),
#                  ylab=paste("dim", theDimB),
#                  main=NA)
#           } else {
#             plot.new()
#           }
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot dim vs dim - color by timepoint
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimVsDim.colByTimepoint.", nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=3*(MDSdims-1), width=3*(MDSdims-1))
#       par(mfrow=c(MDSdims-1, MDSdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDimA in seq(1,MDSdims-1)) {
#         for (theDimB in seq(1,MDSdims)) {
#           if (theDimB == theDimA) next
#           if (theDimB > theDimA) {
#             plot(MDSfit$points[,theDimA],
#                  MDSfit$points[,theDimB],
#                  col=matrixColTPcols[[lib]],
#                  pch=pchVal, cex=pchcexVal,
#                  xlab=paste("dim", theDimA),
#                  ylab=paste("dim", theDimB),
#                  main=NA)
#           } else {
#             plot.new()
#             legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
#           }
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot dim vs dim - color by Euclusters
#       # Kmeans clusters
#       # for (distMethod in c("pearson", "spearman",)) {
#       kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#       kClustColbar =  rep("black", length(kIds))
#       kClustColbar = kClustColors[kIds]
#       names(kClustColbar) = names(kIds)
#       # names(kClustColbar) == colnames(theData)
#       # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#       
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimVsDim.colByEuclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=3*(MDSdims-1), width=3*(MDSdims-1))
#       par(mfrow=c(MDSdims-1, MDSdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDimA in seq(1,MDSdims-1)) {
#         for (theDimB in seq(1,MDSdims)) {
#           if (theDimB == theDimA) next
#           if (theDimB > theDimA) {
#             plot(MDSfit$points[,theDimA],
#                  MDSfit$points[,theDimB],
#                  col=kClustColbar,
#                  pch=pchVal, cex=pchcexVal,
#                  xlab=paste("dim", theDimA),
#                  ylab=paste("dim", theDimB),
#                  main=NA)
#           } else {
#             plot.new()
#             legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#             # # require(scatterplot3d)
#             # scatterplot3d(x=MDSfit$points[timepointThenCovOrd[[lib]], 1],
#             #               y=MDSfit$points[timepointThenCovOrd[[lib]], 2],
#             #               z=MDSfit$points[timepointThenCovOrd[[lib]], 3],
#             #               xlab="dim1", ylab="dim2", zlab="dim3", type="p",
#             #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
#             #               axis=TRUE, grid=TRUE, box=FALSE,
#             #               main=NA, angle=45)
#           }
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       # Plot dim vs dim - color by Spearclusters
#       # Kmeans clusters
#       # for (distMethod in c("pearson", "spearman",)) {
#       kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
#       kClustColbar =  rep("black", length(kIds))
#       kClustColbar = kClustColors[kIds]
#       names(kClustColbar) = names(kIds)
#       # names(kClustColbar) == colnames(theData)
#       # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#       
#       PDFfilename = file.path(MDSoutdir, paste(projectName, ".MDS.DimVsDim.colBySpearmanclusters.", paste( "k=", kCenters, ".", sep="" ), nameText, ".pdf", sep=""))
#       pdf(PDFfilename, height=3*(MDSdims-1), width=3*(MDSdims-1))
#       par(mfrow=c(MDSdims-1, MDSdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#       for (theDimA in seq(1,MDSdims-1)) {
#         for (theDimB in seq(1,MDSdims)) {
#           if (theDimB == theDimA) next
#           if (theDimB > theDimA) {
#             plot(MDSfit$points[,theDimA],
#                  MDSfit$points[,theDimB],
#                  col=kClustColbar,
#                  pch=pchVal, cex=pchcexVal,
#                  xlab=paste("dim", theDimA),
#                  ylab=paste("dim", theDimB),
#                  main=NA)
#           } else {
#             plot.new()
#             legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#             # # require(scatterplot3d)
#             # scatterplot3d(x=MDSfit$points[timepointThenCovOrd[[lib]], 1],
#             #               y=MDSfit$points[timepointThenCovOrd[[lib]], 2],
#             #               z=MDSfit$points[timepointThenCovOrd[[lib]], 3],
#             #               xlab="dim1", ylab="dim2", zlab="dim3", type="p",
#             #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
#             #               axis=TRUE, grid=TRUE, box=FALSE,
#             #               main=NA, angle=45)
#           }
#         }
#       }
#       title(main=paste("MDS: ", titleText, sep=""), line=1, outer=TRUE)
#       # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#       dev.off()
#       
#       
#       if (plot3Dscatters) {
#         # Plot 3D scattrplot: color by euclidean clusters
#         MDS3dir = file.path(MDSoutdir, "3Dscatters", lib, "pngs")
#         if ( ! file.exists( MDS3dir ) ) {
#           if ( ! dir.create( MDS3dir, showWarnings=FALSE, recursive=TRUE )  ) MDS3dir = MDSoutdir
#         }
#         
#         # Kmeans clusters
#         # for (distMethod in c("pearson", "spearman",)) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         # names(kClustColbar) == colnames(theData)
#         # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#         
#         someAngles = seq(5, 360, 5)
#         for (i in 1:length(someAngles)) {
#           PNGfilename = file.path(MDS3dir, paste(sprintf("%03d", i), ".png", sep=""))
#           png(PNGfilename, height=8, width=8, units="in", res=150)
#           
#           # require(scatterplot3d)
#           scatterplot3d(x=MDSfit$points[, 1],
#                         y=MDSfit$points[, 2],
#                         z=MDSfit$points[, 3],
#                         xlab="dim1", ylab="dim2", zlab="dim3", type="p",
#                         color=kClustColbar,
#                         pch=pchVal, cex.symbols=pchcexVal*1.25,
#                         axis=TRUE, grid=TRUE, box=FALSE,
#                         main=someAngles[i], angle=someAngles[i])
#           
#           dev.off()
#         }
#       }
#       
#     }
#   }
# }










# ##################################################
# ##################################################
# ##################################################
# # TSNE
# 
# TSNEmaindir = file.path(mainOutputDir, "TSNE")
# if ( ! file.exists( TSNEmaindir ) ) {
#   if ( ! dir.create( TSNEmaindir, showWarnings=FALSE, recursive=TRUE )  ) TSNEmaindir = mainOutputDir
# }
# 
# TSNEdims = 3
# 
# pchVal=20
# pchcexVal=1.25
# 
# 
# for ( perplexVal in c(5, 25, 50) ) {
#   for (lib in whichLibs) {
#     print(lib)
#     
#     for (useScaled in c(FALSE, TRUE)) {
#       
#       for (useEuclidean in c(FALSE, TRUE)) {
#         
#         if (!useScaled) {
#           TSNEoutdir = file.path(TSNEmaindir, "prop")
#           titleText = paste(lib, "prop", sep=" ")
#           nameText = paste(lib, "prop", sep=".")
#           if (useEuclidean) {
#             TSNEoutdir = file.path(TSNEoutdir, "euclidean")
#             titleText = paste(titleText, "euclidean", sep=" ")
#             nameText = paste(nameText, "euclidean", sep=".")
#             dmDist <- dist(t(dataMatrix.prop[[lib]])) # euclidean distances between the rows = scaler products of differences
#           } else {
#             TSNEoutdir = file.path(TSNEoutdir, "correlation")
#             titleText = paste(titleText, "correlation", sep=" ")
#             nameText = paste(nameText, "correlation", sep=".")
#             dmDist <- as.dist(1-cor(dataMatrix.prop[[lib]]))
#           }
#         } else {
#           TSNEoutdir = file.path(TSNEmaindir, "scaled")
#           titleText = paste(lib, "scaled", sep=" ")
#           nameText = paste(lib, "scaled", sep=".")
#           if (useEuclidean) {
#             TSNEoutdir = file.path(TSNEoutdir, "euclidean")
#             titleText = paste(titleText, "euclidean", sep=" ")
#             nameText = paste(nameText, "euclidean", sep=".")
#             dmDist <- dist(t(dataMatrix.scaled[[lib]])) # euclidean distances between the rows = scaler products of differences
#           } else {
#             TSNEoutdir = file.path(TSNEoutdir, "correlation")
#             titleText = paste(titleText, "correlation", sep=" ")
#             nameText = paste(nameText, "correlation", sep=".")
#             dmDist <- as.dist(1-cor(dataMatrix.scaled[[lib]]))
#           }
#         }  
#         
#         
#         ##################################################
#         # tsne
#         # run Rtsne with default parameters
#         set.seed(1234)
#         tryRet = try(rtsne_out  <- Rtsne(dmDist, is_distance=TRUE, dims=TSNEdims, check_duplicates=FALSE, perplexity=perplexVal))
#         if (class(tryRet) == "try-error" ) next
#         
#         
#         TSNEoutdir = file.path(TSNEoutdir, paste("perplexity", perplexVal, sep=""))
#         if ( ! file.exists( TSNEoutdir ) ) {
#           if ( ! dir.create( TSNEoutdir, showWarnings=FALSE, recursive=TRUE )  ) TSNEoutdir = TSNEmaindir
#         }
#         
#         
#         # Plot indexed dims - color by coverage 
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimsVsIndices.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=4, width=4*TSNEdims)
#         par(mfrow=c(1,TSNEdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDim in seq(1,TSNEdims)) {
#           plot(rtsne_out$Y[,theDim], 
#                col=covCols[[lib]], pch=pchVal, cex=pchcexVal, 
#                xlab="index", ylab=paste("dim", theDim), main=NA)
#           # text(rtsne_out$Y, labels=rownames(mydata))
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
#         dev.off()
#         
#         
#         # Plot indexed dims - ranked by timepoint then coverage; colored by timepoint
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimsVsIndices.rankAndColByTimepoint.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=4, width=4*TSNEdims)
#         par(mfrow=c(1,TSNEdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDim in seq(1,TSNEdims)) {
#           plot(rtsne_out$Y[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                col=matrixColTPcols[[lib]][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                pch=pchVal, cex=pchcexVal, 
#                xlab="index", ylab=paste("dim", theDim), main=NA)
#           # text(rtsne_out$Y, labels=rownames(mydata))
#           if(theDim==1) {
#             legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub=paste("Ranked by timepoint then coverage; colored by timepoint"), line=1, outer=TRUE) 
#         dev.off()       
#         
#         
#         # Plot indexed dims - ranked by timepoint then coverage; colored by k-means cluster by Euclidean distance
#         # Kmeans clusters
#         # for (distMethod in c("pearson", "spearman",)) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         # names(kClustColbar) == colnames(theData)
#         # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#         
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimsVsIndices.rankAndColByEuclusters.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=4, width=4*TSNEdims)
#         par(mfrow=c(1,TSNEdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDim in seq(1,TSNEdims)) {
#           plot(rtsne_out$Y[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                pch=pchVal, cex=pchcexVal, 
#                xlab="index", ylab=paste("dim", theDim), main=NA)
#           # text(rtsne_out$Y, labels=rownames(mydata))
#           if(theDim==1) {
#             legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Euclidean distance"), line=1, outer=TRUE) 
#         dev.off()       
#         
#         
#         # Plot indexed dims - ranked by timepoint then coverage; colored by k-means cluster by Spearman distance
#         # Kmeans clusters
#         # for (distMethod in c("pearson", "spearman",)) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         # names(kClustColbar) == colnames(theData)
#         # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#         
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimsVsIndices.rankAndColBySpearclusters.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=4, width=4*TSNEdims)
#         par(mfrow=c(1,TSNEdims), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDim in seq(1,TSNEdims)) {
#           plot(rtsne_out$Y[,theDim][ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                col=kClustColbar[ covOrd[[lib]][ timepointThenCovOrd[[lib]] ] ], 
#                pch=pchVal, cex=pchcexVal, 
#                xlab="index", ylab=paste("dim", theDim), main=NA)
#           # text(rtsne_out$Y, labels=rownames(mydata))
#           if(theDim==1) {
#             legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub=paste("Ranked by timepoint then coverage; colored by k-means cluster by Spearman distance"), line=1, outer=TRUE) 
#         dev.off()           
#         
#         
#         # Plot dim vs dim - color by coverage 
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimVsDim.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=3*(TSNEdims-1), width=3*(TSNEdims-1))
#         par(mfrow=c(TSNEdims-1, TSNEdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDimA in seq(1,TSNEdims-1)) {
#           for (theDimB in seq(1,TSNEdims)) {
#             if (theDimB == theDimA) next
#             if (theDimB > theDimA) {
#               plot(rtsne_out$Y[,theDimA],
#                    rtsne_out$Y[,theDimB],
#                    col=covCols[[lib]], 
#                    pch=pchVal, cex=pchcexVal, 
#                    xlab=paste("dim", theDimA), 
#                    ylab=paste("dim", theDimB),
#                    main=NA)
#             } else {
#               plot.new()
#             }
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub="dark to light = lower to higher coverage", line=1, outer=TRUE)
#         dev.off()
#         
#         
#         # Plot dim vs dim - color by timepoint 
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimVsDim.colByTimepoint.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=3*(TSNEdims-1), width=3*(TSNEdims-1))
#         par(mfrow=c(TSNEdims-1, TSNEdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDimA in seq(1,TSNEdims-1)) {
#           for (theDimB in seq(1,TSNEdims)) {
#             if (theDimB == theDimA) next
#             if (theDimB > theDimA) {
#               plot(rtsne_out$Y[,theDimA],
#                    rtsne_out$Y[,theDimB],
#                    col=matrixColTPcols[[lib]], 
#                    pch=pchVal, cex=pchcexVal, 
#                    xlab=paste("dim", theDimA), 
#                    ylab=paste("dim", theDimB),
#                    main=NA)
#             } else {
#               plot.new()
#               legend("topleft", bty="n", legend=names(tpCols), fill=tpCols)
#             }
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         # title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#         dev.off()
#         
#         
#         # Plot dim vs dim - color by Euclusters
#         # Kmeans clusters
#         # for (distMethod in c("pearson", "spearman",)) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         # names(kClustColbar) == colnames(theData)
#         # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#         
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimVsDim.colByEuclusters.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=3*(TSNEdims-1), width=3*(TSNEdims-1))
#         par(mfrow=c(TSNEdims-1, TSNEdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDimA in seq(1,TSNEdims-1)) {
#           for (theDimB in seq(1,TSNEdims)) {
#             if (theDimB == theDimA) next
#             if (theDimB > theDimA) {
#               plot(rtsne_out$Y[,theDimA],
#                    rtsne_out$Y[,theDimB],
#                    col=kClustColbar, 
#                    pch=pchVal, cex=pchcexVal, 
#                    xlab=paste("dim", theDimA), 
#                    ylab=paste("dim", theDimB),
#                    main=NA)
#             } else {
#               plot.new()
#               legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#               # require(scatterplot3d)
#               # scatterplot3d(x=rtsne_out$Y[timepointThenCovOrd[[lib]], 1], 
#               #               y=rtsne_out$Y[timepointThenCovOrd[[lib]], 2], 
#               #               z=rtsne_out$Y[timepointThenCovOrd[[lib]], 3], 
#               #               xlab="dim1", ylab="dim2", zlab="dim3", type="p", 
#               #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
#               #               axis=TRUE, grid=TRUE, box=FALSE,
#               #               main=NA, angle=45) 
#             }
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#         dev.off()
#         
#         
#         # Plot dim vs dim - color by Spearclusters
#         # Kmeans clusters
#         # for (distMethod in c("pearson", "spearman",)) {
#         kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["spearman"]]
#         kClustColbar =  rep("black", length(kIds))
#         kClustColbar = kClustColors[kIds]
#         names(kClustColbar) = names(kIds)
#         # names(kClustColbar) == colnames(theData)
#         # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#         
#         PDFfilename = file.path(TSNEoutdir, paste(projectName, ".tsne.DimVsDim.colBySpearclusters.", nameText, ".pdf", sep=""))
#         pdf(PDFfilename, height=3*(TSNEdims-1), width=3*(TSNEdims-1))
#         par(mfrow=c(TSNEdims-1, TSNEdims-1), mar=c(4,4,2,2), oma=c(2,2,4,2))
#         for (theDimA in seq(1,TSNEdims-1)) {
#           for (theDimB in seq(1,TSNEdims)) {
#             if (theDimB == theDimA) next
#             if (theDimB > theDimA) {
#               plot(rtsne_out$Y[,theDimA],
#                    rtsne_out$Y[,theDimB],
#                    col=kClustColbar, 
#                    pch=pchVal, cex=pchcexVal, 
#                    xlab=paste("dim", theDimA), 
#                    ylab=paste("dim", theDimB),
#                    main=NA)
#             } else {
#               plot.new()
#               legend("topleft", bty="n", legend=names(kClustColors), fill=kClustColors)
#               # require(scatterplot3d)
#               # scatterplot3d(x=rtsne_out$Y[timepointThenCovOrd[[lib]], 1], 
#               #               y=rtsne_out$Y[timepointThenCovOrd[[lib]], 2], 
#               #               z=rtsne_out$Y[timepointThenCovOrd[[lib]], 3], 
#               #               xlab="dim1", ylab="dim2", zlab="dim3", type="p", 
#               #               color=kClustColbar, pch=pchVal, cex.symbols=pchcexVal*0.85,
#               #               axis=TRUE, grid=TRUE, box=FALSE,
#               #               main=NA, angle=45) 
#             }
#           }
#         }
#         title(main=paste("tSNE: ", titleText, sep=""), line=1, outer=TRUE)
#         title(sub=paste(names(tpCols), ":", tpCols,  sep="", collapse=", "), line=1, outer=TRUE)
#         dev.off()        
#         
#         if (plot3Dscatters) {
#           
#           # Plot 3D scattrplot: color by euclidean clusters
#           TNSE3dir = file.path(TSNEoutdir, "3Dscatters", lib, "pngs")
#           if ( ! file.exists( TNSE3dir ) ) {
#             if ( ! dir.create( TNSE3dir, showWarnings=FALSE, recursive=TRUE )  ) TNSE3dir = TSNEoutdir
#           }
#           
#           # Kmeans clusters
#           # for (distMethod in c("pearson", "spearman",)) {
#           kIds = kClusters[[lib]][[ifelse(useScaled, "scaled", "prop")]][["euclidean"]]
#           kClustColbar =  rep("black", length(kIds))
#           kClustColbar = kClustColors[kIds]
#           names(kClustColbar) = names(kIds)
#           # names(kClustColbar) == colnames(theData)
#           # names(kClustColbar) == names(prcompRes$x[ timepointThenCovOrd[[lib]], PCAaLabel])
#           
#           someAngles = seq(5, 360, 5)
#           for (i in 1:length(someAngles)) {
#             PNGfilename = file.path(TNSE3dir, paste(sprintf("%03d", i), ".png", sep=""))
#             png(PNGfilename, height=8, width=8, units="in", res=150)
#             
#             # require(scatterplot3d)
#             scatterplot3d(x=rtsne_out$Y[, 1], 
#                           y=rtsne_out$Y[, 2], 
#                           z=rtsne_out$Y[, 3], 
#                           xlab="dim1", ylab="dim2", zlab="dim3", type="p", 
#                           color=kClustColbar, 
#                           pch=pchVal, cex.symbols=pchcexVal*1.25,
#                           axis=TRUE, grid=TRUE, box=FALSE,
#                           main=someAngles[i], angle=someAngles[i]) 
#             
#             dev.off()
#           }
#         } 
#         
#       }
#     }
#   }
# }
# 
# # }










##################################################
##################################################
##################################################
# Cell cycle analysis



##################################################
# Load kmeans group IDs

otherProjectName="non-allelic_sci-Hi-C"

inputDir = "../data/autosomal_non-allelic_cell_cycle_cluster_ids/"

nThresh=1000

bins2Use="50kb-8Mb"
# bins2Use="1kb-80Mb"

useScaled=TRUE

# useMethod="euclidean"
useMethod="spearman"

kCenters = 4
# kCenters = 8

# autosomalCDPs=FALSE
autosomalCDPs=TRUE

kclustFiledir = file.path(inputDir, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), bins2Use, paste( "kclust", ".k=", kCenters, sep="" ), "clusterIDs")
if ( ! file.exists( kclustFiledir ) ) {
                          if ( ! dir.create( kclustFiledir, showWarnings=FALSE, recursive=TRUE )  ) kclustFiledir = kclustMaindir
                          }
                          
kClustersFromFile = list()

for (lib in whichLibs) {
  print(lib)
  nameText = paste(lib, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), 
                   bins2Use,
                   ifelse(useScaled, "scaled", "unscaled"),
                   paste("kclust.dist", useMethod, sep=""),
                   paste("kCenters", kCenters, sep=""),
                   sep=".")
  kclustFileName = file.path(kclustFiledir, paste(otherProjectName, ".", nameText, ".tsv", sep=""))
  kClustersFromFile[[lib]] =  read.table(file=kclustFileName, sep="\t", header=FALSE)
}


# # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# # kClustColors =  rainbow(kCenters)
# set.seed(1234)
# kClustColors =  sample(color, kCenters)
kClustColors = colorRampPalette(c("#00ffff", "#001e1e"))(kCenters)
names(kClustColors) = paste("group", 1:kCenters)

kClustDesc = paste(bins2Use,
                   ".", ifelse(useScaled, "scaled", "unscaled"), 
                   ".", paste("kclust.dist", useMethod, sep=""),
                   ".", paste("kCenters", kCenters, sep=""), 
                   sep="")



##################################################
# 20210501
# Load XCI state

# 20210502
# Note: For Patski cell cycle plots colored by XCI status 
#       use XCI analysis based on 100 contacts per chr1 and chrX allele rather than 50.
patski_along_100 = FALSE
# patski_along_100 = TRUE

if (patski_along_100) {
  XCIstatusFileDir = paste0("..data/LMD_XCI_cells_countthresh100")
} else {
  XCIstatusFileDir = paste0("..data/LMD_XCI_cells_countthresh50")
}

XCIstatusFromFile = list()

for (lib in whichLibs) {
  
  XCIstatusFromFile[[lib]] = list()
                    
  for (XCIstatus in c("altXCIcells", "refXCIcells", "nonXCIcells")) {

    cat(lib, "\t", XCIstatus, "\n", sep="")
    
    nameText = paste(lib, XCIstatus, sep=".")
    XCIstatusFileName = file.path(XCIstatusFileDir, paste(nameText, ".tsv", sep=""))
    if ( ! file.exists( XCIstatusFileName ) ) {
      XCIstatusFromFile[[lib]][[XCIstatus]] = NULL
    } else {
      XCIstatusFromFile[[lib]][[XCIstatus]] = try( read.table(file=XCIstatusFileName, sep="\t", header=FALSE) )
     if ( class(XCIstatusFromFile[[lib]][[XCIstatus]]) == "try-error" ) {
       XCIstatusFromFile[[lib]][[XCIstatus]] = NULL 
     }
    }
    
  }
  
}



##################################################
# 20181022
# Short range and mitotic contact proportions

propContactsShort = list()
propContactsMitotic = list()
propContactsLong = list()
propNonMitotic = list()

mitoticRanges = c("2to12Mb", "1to20Mb")

for (mitoticRange in mitoticRanges){
  propContactsShort[[mitoticRange]] = list()
  propContactsMitotic[[mitoticRange]] = list()
  propContactsLong[[mitoticRange]] = list()
  propNonMitotic[[mitoticRange]] = list()
  
  for(lib in whichLibs) {
    
    if (mitoticRange == "2to12Mb") {
      propContactsShort[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist<=2000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propContactsMitotic[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist>2000000 & binDist<=12000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propContactsLong[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist>12000000 & binDist<=80000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propNonMitotic[[mitoticRange]][[lib]] = propContactsShort[[mitoticRange]][[lib]] + propContactsLong[[mitoticRange]][[lib]]
      summary(propContactsMitotic[[mitoticRange]][[lib]])
      summary(propNonMitotic[[mitoticRange]][[lib]])
    } else if (mitoticRange == "1to20Mb") {
      propContactsShort[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist<=1000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propContactsMitotic[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist>1000000 & binDist<=20000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propContactsLong[[mitoticRange]][[lib]] = colSums(dataMatrix.OG[[lib]][binDist>20000000 & binDist<=80000000,]) / colSums(dataMatrix.OG[[lib]][binDist<Inf,]) 
      propNonMitotic[[mitoticRange]][[lib]] = propContactsShort[[mitoticRange]][[lib]] + propContactsLong[[mitoticRange]][[lib]]
      summary(propContactsMitotic[[mitoticRange]][[lib]])
      summary(propNonMitotic[[mitoticRange]][[lib]])
    }
  }
}



##################################################
# 20181022
# Cell cycle plot based on short range and mitotic contact proportions
# These plots will be the same regarless of the distance range used for the contact decay analysis
# except for the coloring by cluster group derived from the contact decay analysis. 

ccplotsMaindir = file.path(supermainOutputDir, paste( "countThresh", nThresh, ifelse(autosomalCDPs, ".autosomalCDPs", ""), sep=""), "cellCyclePlots")
if ( ! file.exists( ccplotsMaindir ) ) {
  if ( ! dir.create( ccplotsMaindir, showWarnings=FALSE, recursive=TRUE )  ) ccplotsMaindir = subsupermainOutputDir
}

axisLims=list()
axisLims[["short"]] = c(25, 95)
axisLims[["long"]] = c(0, 35)
if (patski_along_100) {
  axisLims[["mitotic"]] = c(0, 60)
} else {
  axisLims[["mitotic"]] = c(0, 70)
}

axisLabels=list()

for (mitoticRange in mitoticRanges){
  
  axisLabels[[mitoticRange]] = list()
  axisLabels[[mitoticRange]][["short"]] = paste("% short-range contacts (0-", ifelse(mitoticRange=="2to12Mb", "2", "1"),"Mb)", sep="")
  axisLabels[[mitoticRange]][["long"]] = paste("% long-range contacts (", ifelse(mitoticRange=="2to12Mb", "12", "20"), "-80Mb)", sep="")
  axisLabels[[mitoticRange]][["mitotic"]] = paste("% mitotic contacts (", sub("to", "-", mitoticRange), ")", sep="")
  
  for(lib in whichLibs) {
    
    titleText = paste(lib, sep=" ")
    nameText = paste(lib, sep=".")
    
    #################
    # Cluster ID colors
    
    ccplotsSubdir = file.path(ccplotsMaindir, 
                              paste("mitoticRange", mitoticRange, sep=""),
                              paste("ColorByClusterIDs", ".", kClustDesc, sep=""))
    if ( ! file.exists( ccplotsSubdir ) ) {
      if ( ! dir.create( ccplotsSubdir, showWarnings=FALSE, recursive=TRUE )  ) ccplotsSubdir = ccplotsMaindir
    }
    
    PDFfilename = file.path(ccplotsSubdir, paste(projectName, ".cellCyclePlots",
                                                 ".mitoticRange", mitoticRange, 
                                                 ".clusterIDs",  ".", 
                                                 kClustDesc, ".", 
                                                 nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=8, width=8)
    par(mfrow=c(2,2), mar=par()$mar-c(0,0,2,2), oma=par()$oma+c(0,0,2,0))
    
    kIds = kClustersFromFile[[lib]][,2]
    kClustColbar = kClustColors[kIds]
    names(kClustColbar) = names(kIds)
    
    plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
         propContactsShort[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["short"]],
         xlab=axisLabels[[mitoticRange]][["mitotic"]], 
         ylab=axisLabels[[mitoticRange]][["short"]], 
         main="")
    points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
           propContactsShort[[mitoticRange]][[lib]]*100, 
           pch=20, col=kClustColbar)
    
    y = propContactsShort[[mitoticRange]][[lib]]*100
    x = propContactsMitotic[[mitoticRange]][[lib]]*100
    y = y[x<30]; x = x[x<30]
    # lnmod <- lm(y~x)
    # # abline(lnmod, lty="solid", xlim=c(0,30)) # xlim doesn't work for abline
    # xs <- range(0, 30)
    # ys <- predict(lnmod, newdata=data.frame(x=xs))
    # lines(xs, ys, lty="dashed", lwd=2) 
    
    clusterCounts = table(kIds)
    legendText = c(paste("all", " (n = ", length(matrixColTPs[[lib]]), ")", sep=""))
    for (i in 1:length(clusterCounts)) {
      legendText = c(legendText, paste( "group ", i, " (n = ", clusterCounts[i], ")", sep=""))
      
    }
    legend("topright", legend=legendText, cex=0.7,
           # paste("slope (x<30%): ", round(lnmod$coefficients[2], 2), sep="")),
           fill=c("black", kClustColors, NA), border=NA, bty='n')
    
    plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
         propContactsLong[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["long"]],
         xlab=axisLabels[[mitoticRange]][["mitotic"]], 
         ylab=axisLabels[[mitoticRange]][["long"]], 
         main="")
    points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
           propContactsLong[[mitoticRange]][[lib]]*100, 
           pch=20, col=kClustColbar)
    
    plot(propContactsShort[[mitoticRange]][[lib]]*100, 
         propContactsLong[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["short"]] , ylim=axisLims[["long"]],
         xlab=axisLabels[[mitoticRange]][["short"]], 
         ylab=axisLabels[[mitoticRange]][["long"]], 
         main="")
    points(propContactsShort[[mitoticRange]][[lib]]*100, 
           propContactsLong[[mitoticRange]][[lib]]*100, 
           pch=20, col=kClustColbar)
    
    sp3d <- scatterplot3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
                          propContactsShort[[mitoticRange]][[lib]]*100, 
                          propContactsLong[[mitoticRange]][[lib]]*100, 
                          xlim=axisLims[["mitotic"]],
                          ylim=axisLims[["short"]],
                          zlim=axisLims[["long"]],
                          xlab=axisLabels[[mitoticRange]][["mitotic"]], 
                          ylab=axisLabels[[mitoticRange]][["short"]], 
                          zlab=axisLabels[[mitoticRange]][["long"]],
                          type="p", pch=21,
                          axis=TRUE, grid=TRUE, box=TRUE, 
                          angle=150, scale.y=1.25)
    # addgrids3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #            propContactsShort[[mitoticRange]][[lib]]*100, 
    #            propContactsLong[[mitoticRange]][[lib]]*100, 
    #            grid = c("xy", "xz", "yz"), 
    #            angle=150, scale.y=1.25)
    sp3d$points3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
                  propContactsShort[[mitoticRange]][[lib]]*100, 
                  propContactsLong[[mitoticRange]][[lib]]*100, 
                  type="p", pch=21, bg=kClustColbar)
    
    title( main=paste("Cell cycle inferred from short range and mitotic contact proportions\n", 
                      titleText, "; mitoticRange: ", mitoticRange, sep=""), outer=TRUE)
    
    dev.off()
    
    
    
    #################
    # XCI status colors
    
    library(plotrix)
    
    ccplotsSubdir = file.path(ccplotsMaindir, 
                              paste("mitoticRange", mitoticRange, sep=""),
                              paste("ColorByXCIstatus", 
                                    ifelse(allelic_CDP_analysis_thresh!=50, 
                                           paste0( ".allelic_count_thresh", allelic_CDP_analysis_thresh ), "" ),
                                    sep=""))
    if ( ! file.exists( ccplotsSubdir ) ) {
      if ( ! dir.create( ccplotsSubdir, showWarnings=FALSE, recursive=TRUE )  ) ccplotsSubdir = ccplotsMaindir
    }
    
    XCIstatus_col = c( rgb(150, 239, 147, maxColorValue = 255), # not XCI
                       rgb(17, 24, 147, maxColorValue = 255), # 129 XCI
                       rgb(170, 68, 68, maxColorValue = 255) # cast XCI
                       )
    names(XCIstatus_col) = c( "nonXCIcells", "refXCIcells", "altXCIcells" )
    
    XCIstatus_per_cell = rep( "NA", length( propContactsShort[[mitoticRange]][[lib]] ) )
    XCIstatusCol_per_cell = rep( "white", length( propContactsShort[[mitoticRange]][[lib]] ) )
    names(XCIstatus_per_cell) = names(propContactsShort[[mitoticRange]][[lib]])
    names(XCIstatus_per_cell)
    for ( XCIstatus in names(XCIstatusFromFile[[lib]]) ) {
      XCIstatus.cellIDs = XCIstatusFromFile[[lib]][[XCIstatus]]
      currStatus_idx = match( XCIstatusFromFile[[lib]][[XCIstatus]][,1], paste( lib, names(XCIstatus_per_cell), sep="." ) )
      XCIstatus_per_cell[ currStatus_idx ] = XCIstatus
      XCIstatusCol_per_cell[ currStatus_idx ] = XCIstatus_col[XCIstatus]
    }
    XCIstatus_per_cell = factor(XCIstatus_per_cell, 
                             levels=c("NA", "nonXCIcells", "refXCIcells", "altXCIcells"))
    
    XCI_status_counts = table( XCIstatus_per_cell )
    legendText = c(paste("all", " (n = ", length(XCIstatus_per_cell), ")", sep=""))
    for (i in 1:length(XCI_status_counts)) {
      legendText = c(legendText, paste( names(XCI_status_counts)[i], "(n = ", XCI_status_counts[i], ")", sep=""))
      
    }
    
    XCI_status_group_counts = table( XCIstatus_per_cell, kIds ) 
    XCI_status_group_proportions = XCI_status_group_counts / rowSums( XCI_status_group_counts )
    XCI_status_group_percent = round(100 * XCI_status_group_proportions)

    cluster_counts = table(kIds)
    cluster_counts_proportions = cluster_counts / sum( cluster_counts )
    cluster_counts_percent = round(100 * cluster_counts_proportions)
    
    XCI_status_group_enrichment = c()
    for (i in 1:nrow(XCI_status_group_proportions) ) {
      XCI_status_group_enrichment = rbind(XCI_status_group_enrichment,
                                          XCI_status_group_proportions[i, ] / cluster_counts_proportions )
    }
    rownames(XCI_status_group_enrichment) =  rownames(XCI_status_group_proportions)
    colnames(XCI_status_group_enrichment) =  colnames(XCI_status_group_proportions)
    XCI_status_group_enrichment = round(XCI_status_group_enrichment, 1)
    
    XCI_status_group_counts
    # Patski
    # kIds
    # XCIstatus_per_cell    1    2    3    4
    # NA            43  504 1263 1522
    # nonXCIcells   16    3   13  186
    # refXCIcells    5    5   26  907
    # altXCIcells    0    0    0    0
  
    rowSums(XCI_status_group_counts)
    # NA nonXCIcells refXCIcells altXCIcells 
    # 3332         218         943           0 
    # 218 + 943 = 1161
    
    PDFfilename = file.path(ccplotsSubdir, paste(projectName, ".cellCyclePlots",
                                                 ".mitoticRange", mitoticRange, 
                                                 ".XCIstatus", 
                                                 nameText, ".pdf", sep=""))
    pdf(PDFfilename, height=10, width=8)
    par(mfrow=c(3,2), mar=par()$mar-c(0,0,2,2), oma=par()$oma+c(0,0,2,0))

    plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
         propContactsShort[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["short"]],
         xlab=axisLabels[[mitoticRange]][["mitotic"]], 
         ylab=axisLabels[[mitoticRange]][["short"]], 
         main="")
    points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
           propContactsShort[[mitoticRange]][[lib]]*100, 
           pch=20, col=XCIstatusCol_per_cell)
    
    y = propContactsShort[[mitoticRange]][[lib]]*100
    x = propContactsMitotic[[mitoticRange]][[lib]]*100
    y = y[x<30]; x = x[x<30]
    # lnmod <- lm(y~x)
    # # abline(lnmod, lty="solid", xlim=c(0,30)) # xlim doesn't work for abline
    # xs <- range(0, 30)
    # ys <- predict(lnmod, newdata=data.frame(x=xs))
    # lines(xs, ys, lty="dashed", lwd=2) 
    
    legend("topright", legend=legendText, cex=0.7, 
           # paste("slope (x<30%): ", round(lnmod$coefficients[2], 2), sep="")),
           fill=c("black", "white", XCIstatus_col, NA), border="black", bty='n')
   
    plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
         propContactsLong[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["long"]],
         xlab=axisLabels[[mitoticRange]][["mitotic"]], 
         ylab=axisLabels[[mitoticRange]][["long"]], 
         main="")
    points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
           propContactsLong[[mitoticRange]][[lib]]*100, 
           pch=20, col=XCIstatusCol_per_cell)
    
    plot(propContactsShort[[mitoticRange]][[lib]]*100, 
         propContactsLong[[mitoticRange]][[lib]]*100, 
         pch=21, xlim=axisLims[["short"]] , ylim=axisLims[["long"]],
         xlab=axisLabels[[mitoticRange]][["short"]], 
         ylab=axisLabels[[mitoticRange]][["long"]], 
         main="")
    points(propContactsShort[[mitoticRange]][[lib]]*100, 
           propContactsLong[[mitoticRange]][[lib]]*100, 
           pch=20, col=XCIstatusCol_per_cell)
    
    sp3d <- scatterplot3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
                          propContactsShort[[mitoticRange]][[lib]]*100, 
                          propContactsLong[[mitoticRange]][[lib]]*100, 
                          xlim=axisLims[["mitotic"]],
                          ylim=axisLims[["short"]],
                          zlim=axisLims[["long"]],
                          xlab=axisLabels[[mitoticRange]][["mitotic"]], 
                          ylab=axisLabels[[mitoticRange]][["short"]], 
                          zlab=axisLabels[[mitoticRange]][["long"]],
                          type="p", pch=21,
                          axis=TRUE, grid=TRUE, box=TRUE, 
                          angle=150, scale.y=1.25)
    # addgrids3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #            propContactsShort[[mitoticRange]][[lib]]*100, 
    #            propContactsLong[[mitoticRange]][[lib]]*100, 
    #            grid = c("xy", "xz", "yz"), 
    #            angle=150, scale.y=1.25)
    sp3d$points3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
                  propContactsShort[[mitoticRange]][[lib]]*100, 
                  propContactsLong[[mitoticRange]][[lib]]*100, 
                  type="p", pch=21, bg=XCIstatusCol_per_cell)
    
    plot(c(0,100), c(0,100), type="n", axes=FALSE, xlab="", ylab="")
    addtable2plot(50, 0, round( XCI_status_group_proportions, 2), bty="o",
                  xjust=0.5, yjust=1, 
                  display.rownames=TRUE,
                  hlines=TRUE, vlines=TRUE, 
                  title=paste0( lib, " proportions per group by XCI status"), cex=1.25)
    addtable2plot(50, 70, round( t(cluster_counts_proportions), 2), bty="o",
                  xjust=0.5, yjust=1, 
                  display.rownames=FALSE,
                  hlines=TRUE, vlines=TRUE, 
                  title=paste0( lib, " overall proportions per group"), cex=1.25)
    
    plot(c(0,100), c(0,100), type="n", axes=FALSE, xlab="", ylab="")
    addtable2plot(50, 0, XCI_status_group_enrichment, bty="o",
                  xjust=0.5, yjust=1, 
                  display.rownames=TRUE,
                  hlines=TRUE, vlines=TRUE, 
                  title=paste0( lib, " XCI status enrichment per group"), cex=1.25)
    
    title( main=paste("Cell cycle inferred from short range and mitotic contact proportions\n", 
                      titleText, "; mitoticRange: ", mitoticRange, sep=""), outer=TRUE)
    
    dev.off()
    
    
    
    # #################
    # # Timepoint colors
    # 
    # ccplotsSubdir = file.path(ccplotsMaindir, 
    #                           paste("mitoticRange", mitoticRange, sep=""),
    #                           paste("ColorByTimepoints", sep=""))
    # if ( ! file.exists( ccplotsSubdir ) ) {
    #   if ( ! dir.create( ccplotsSubdir, showWarnings=FALSE, recursive=TRUE )  ) ccplotsSubdir = ccplotsMaindir
    # }
    # 
    # PDFfilename = file.path(ccplotsSubdir, paste(projectName, ".cellCyclePlots.mitoticRange", mitoticRange, ".", nameText, ".pdf", sep=""))
    # pdf(PDFfilename, height=8, width=8)
    # par(mfrow=c(2,2), mar=par()$mar-c(0,0,2,2), oma=par()$oma+c(0,0,2,0))
    # 
    # plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #      propContactsShort[[mitoticRange]][[lib]]*100, 
    #      pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["short"]],
    #      xlab=axisLabels[[mitoticRange]][["mitotic"]], 
    #      ylab=axisLabels[[mitoticRange]][["short"]], 
    #      main="")
    # points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #        propContactsShort[[mitoticRange]][[lib]]*100, 
    #        pch=20, col=matrixColTPcols[[lib]])
    # 
    # y = propContactsShort[[mitoticRange]][[lib]]*100
    # x = propContactsMitotic[[mitoticRange]][[lib]]*100
    # y = y[x<30]; x = x[x<30]
    # # lnmod <- lm(y~x)
    # # # abline(lnmod, lty="solid", xlim=c(0,30)) # xlim doesn't work for abline
    # # xs <- range(0, 30)
    # # ys <- predict(lnmod, newdata=data.frame(x=xs))
    # # lines(xs, ys, lty="dashed", lwd=2) 
    # 
    # nML1 = sum(matrixColTPs[[lib]]=="ML1")
    # nML2= sum(matrixColTPs[[lib]]=="ML2")
    # nPL1 = sum(matrixColTPs[[lib]]=="PL1")
    # nPL2 = sum(matrixColTPs[[lib]]=="PL2")
    # legend("topright", legend=c(paste("all", " (n = ", length(matrixColTPs[[lib]]), ")", sep=""),
    #                             paste( "ML1 (n = ", nML1, ")", sep=""),
    #                             paste( "ML2 (n = ", nML2, ")", sep=""),
    #                             paste( "PL1 (n = ", nPL1, ")", sep=""),
    #                             paste( "PL2 (n = ", nPL2, ")", sep="")),
    #                             # paste("slope (x<30%): ", round(lnmod$coefficients[2], 2), sep="")),
    #        fill=c("black", tpCols, NA), border=NA, bty='n')
    # 
    # plot(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #      propContactsLong[[mitoticRange]][[lib]]*100, 
    #      pch=21, xlim=axisLims[["mitotic"]], ylim=axisLims[["long"]],
    #      xlab=axisLabels[[mitoticRange]][["mitotic"]], 
    #      ylab=axisLabels[[mitoticRange]][["long"]], 
    #      main="")
    # points(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #        propContactsLong[[mitoticRange]][[lib]]*100, 
    #        pch=20, col=matrixColTPcols[[lib]])
    # 
    # plot(propContactsShort[[mitoticRange]][[lib]]*100, 
    #      propContactsLong[[mitoticRange]][[lib]]*100, 
    #      pch=21, xlim=axisLims[["short"]], ylim=axisLims[["long"]],
    #      xlab=axisLabels[[mitoticRange]][["short"]],
    #      ylab=axisLabels[[mitoticRange]][["long"]], 
    #      main="")
    # points(propContactsShort[[mitoticRange]][[lib]]*100, 
    #        propContactsLong[[mitoticRange]][[lib]]*100, 
    #        pch=20, col=matrixColTPcols[[lib]])
    # 
    # sp3d <- scatterplot3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #                       propContactsShort[[mitoticRange]][[lib]]*100, 
    #                       propContactsLong[[mitoticRange]][[lib]]*100, 
    #                       xlim=axisLims[["mitotic"]],
    #                       ylim=axisLims[["short"]],
    #                       zlim=axisLims[["long"]],
    #                       xlab=axisLabels[[mitoticRange]][["mitotic"]], 
    #                       ylab=axisLabels[[mitoticRange]][["short"]], 
    #                       zlab=axisLabels[[mitoticRange]][["long"]], 
    #                       type="p", pch=21,
    #                       axis=TRUE, grid=TRUE, box=TRUE, 
    #                       angle=150, scale.y=1.25)
    # # addgrids3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    # #            propContactsShort[[mitoticRange]][[lib]]*100, 
    # #            propContactsLong[[mitoticRange]][[lib]]*100, 
    # #            grid = c("xy", "xz", "yz"), 
    # #            angle=150, scale.y=1.25)
    # sp3d$points3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #               propContactsShort[[mitoticRange]][[lib]]*100, 
    #               propContactsLong[[mitoticRange]][[lib]]*100, 
    #               type="p", pch=21, bg=matrixColTPcols[[lib]])
    # 
    # title( main=paste("Cell cycle inferred from short range and mitotic contact proportions\n", 
    #                   titleText, "; mitoticRange: ", mitoticRange, sep=""), outer=TRUE)
    # 
    # dev.off()
    
    
    
    # # Separate timepoints
    # if (lib==combined_label) {
    #   
    #   #################
    #   # Different cols
    #   for (plotTP in names(tpCols)) {
    #     colsToUse=tpCols[plotTP]
    #     
    #     PDFfilename = file.path(ccplotsSubdir, paste(projectName, ".cellCyclePlots.mitoticRange", mitoticRange, ".", nameText, ".", plotTP, ".pdf", sep=""))
    #     pdf(PDFfilename, height=8, width=8)
    #     par(mfrow=c(2,2), mar=par()$mar-c(0,0,2,2), oma=par()$oma+c(0,0,2,0))
    #     
    #     plot(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["short"]],
    #          xlab=axisLabels[[mitoticRange]][["mitotic"]],
    #          ylab=axisLabels[[mitoticRange]][["short"]],
    #          main="")
    #     points(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            pch=20, col=colsToUse)
    #     
    #     y = propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100
    #     x = propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100
    #     y = y[x<30]; x = x[x<30]
    #     lnmod <- lm(y~x)
    #     # abline(lnmod, lty="solid", xlim=c(0,30)) # xlim doesn't work for abline
    #     xs <- range(0, 30)
    #     ys <- predict(lnmod, newdata=data.frame(x=xs))
    #     lines(xs, ys, lty="dashed", lwd=2) 
    #     
    #     nd = sum(matrixColTPs[[lib]]==plotTP)
    #     legend("topright", legend=c(paste("all", " (n = ", length(matrixColTPs[[lib]]), ")", sep=""),
    #                                 paste( plotTP, " (n = ", nd, ")", sep=""),
    #                                 paste("slope (x<30%): ", round(lnmod$coefficients[2], 2), sep="")),
    #            fill=c("black", colsToUse, NA), border=NA, bty='n')
    #     
    #     plot(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          pch=21, xlim=axisLims[["mitotic"]] , ylim=axisLims[["long"]],
    #          xlab=axisLabels[[mitoticRange]][["mitotic"]],
    #          ylab=axisLabels[[mitoticRange]][["long"]],
    #          main="")
    #     points(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            pch=20, col=colsToUse)
    #     
    #     plot(propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #          pch=21, xlim=axisLims[["short"]] , ylim=axisLims[["long"]],
    #          xlab=axisLabels[[mitoticRange]][["short"]], 
    #          ylab=axisLabels[[mitoticRange]][["long"]],
    #          main="")
    #     points(propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #            pch=20, col=colsToUse)
    #     
    #     sp3d <- scatterplot3d(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                           propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                           propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                           xlim=axisLims[["mitotic"]],
    #                           ylim=axisLims[["short"]],
    #                           zlim=axisLims[["long"]],
    #                           xlab=axisLabels[[mitoticRange]][["mitotic"]],
    #                           ylab=axisLabels[[mitoticRange]][["short"]],
    #                           zlab=axisLabels[[mitoticRange]][["long"]],
    #                           type="p", pch=21,
    #                           axis=TRUE, grid=TRUE, box=TRUE, 
    #                           angle=150, scale.y=1.25)
    #     # addgrids3d(propContactsMitotic[[mitoticRange]][[lib]]*100, 
    #     #            propContactsShort[[mitoticRange]][[lib]]*100, 
    #     #            propContactsLong[[mitoticRange]][[lib]]*100, 
    #     #            grid = c("xy", "xz", "yz"), 
    #     #            angle=150, scale.y=1.25)
    #     sp3d$points3d(propContactsMitotic[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                   propContactsShort[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                   propContactsLong[[mitoticRange]][[lib]][matrixColTPs[[lib]] == plotTP]*100, 
    #                   type="p", pch=21, bg=colsToUse)
    #     
    #     title( main=paste(plotTP, " cell cycle inferred from short range and mitotic contact proportions\n", 
    #                       titleText, "; mitoticRange: ", mitoticRange, sep=""), outer=TRUE)
    #     
    #     dev.off()
    #     
    #   }      
    
  }
}
