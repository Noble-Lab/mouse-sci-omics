##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200720 allelic sci-RNA-seq analysis using monocle 2 
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

require(monocle)
require(data.table)



##################################################
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}










##################################################
##################################################
##################################################
# Input

projectName = "allelic_sci-ATAC-seq"

inputDir = "../data/allelic_sci-ATAC-seq_count_matrices/"

superMainOutputDir = paste("../results/allelic_sci-ATAC-seq_monocle2/")
if ( ! file.exists( superMainOutputDir ) ) {
    if ( ! dir.create( superMainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) superMainOutputDir = getwd()
}
setwd(superMainOutputDir)
getwd()

tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")

sexCols=c("darkred", "midnightblue")
names(sexCols) = c("female", "male")

alleleCols=c("pink", "lightblue")
names(alleleCols) = c("129", "Cast")

pchVal=20
pchcexVal=1.25

theChroms = paste( "chr", c( c(1:19), "X"), sep="" )










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
barcodeInfo$gender = "female" 
barcodeInfo$gender[grepl("M", barcodeInfo$cellType)] = "male"

# Fix time points
barcodeInfo$timepoint[barcodeInfo$timepoint=='d4'] = 'd3' 
barcodeInfo$timepoint[barcodeInfo$timepoint=='d10'] = 'd11' 

head(barcodeInfo)
# barcode wellID cellType timepoint gender
# 1 CTGAAGCTCGCTTATCGTGGCGGAGGTCCAGGACGT     7G     Fd10       d11 female
# 2 ATTCAGAACATTATTATTTTAATTAGCCAGGCGAAG     5E    F_NPC       NPC female
# 3 CTGAAGCTAGCCGTAGTTTACGATCATCTAATCTTA     7F      Fd4        d3 female
# 4 TCTCGCGCGATCTTCGCATAGTTCTAGATAATCTTA    11F     Md10       d11   male
# 5 GAGATTCCTATTGCTGGATAGTTCTAGAGTACTGAC     4H      Fd7        d7 female
# 6 CGCTCATTAAGTAGCTCAGTAATGATCGTAATCTTA     3F      Fd4        d3 female










##################################################
##################################################
##################################################
# Load counts data and transform into CDS

countsFile = file.path(inputDir, "sciATAC_mouseDiff.allelic.counts.biallelicPeaks.txt")

count.table.OG = as.data.frame(fread(countsFile, header=FALSE, sep="\t"))
dim(count.table.OG)
# 8583610       3

colnames(count.table.OG) = c("locus", "barcode", "counts")
head(count.table.OG)










##################################################
##################################################
##################################################
# Transform table into count matrix - unfiltered

uniqLoci.OG = unique(count.table.OG$locus)
length(uniqLoci.OG)
# 156,929

count.matrix.OG = data.table::dcast(count.table.OG, locus ~ barcode, value.var="counts", 
                                    drop=FALSE, fill=0)
rownames(count.matrix.OG) = count.matrix.OG[,1]
count.matrix.OG = count.matrix.OG[,-1]

dim(count.matrix.OG)
# [1] 156,929   6914

# *********************************************************
# Later after filtering, removing singlet haplocells and 
# some removing allelically skewed regions (see below):
dim(count.matrix.OG)
# [1] 156,895   3902
# *********************************************************

# Rename loci
locIDs.OG <- rownames(count.matrix.OG) 
locIDs.OG = sub("_", "-", sub("_", ":", locIDs.OG))
rownames(count.matrix.OG) = locIDs.OG

# cell IDs
cell_list.OG <- colnames(count.matrix.OG)

# 20190313
# Remove singlets
i = 1
singletList = c()
while( i < length(cell_list.OG) ) {
    if(sub("_...$", "", cell_list.OG[i]) != sub("_...$", "", cell_list.OG[i+1])) {
        singletList = c(singletList, i)
        i = i + 1
    } else {
        i = i + 2
    }
}
length(singletList)
#  146

cell_list.OG = cell_list.OG[-singletList]
count.matrix.OG <- count.matrix.OG[,-singletList]
dim(count.matrix.OG)
# 156,929   6768










##################################################
##################################################
##################################################
# Counts tables on unfiltered data
# V2 only

bcIdx = match( sub("_.*$", "", colnames(count.matrix.OG)), barcodeInfo$barcode )
all( sub("_.*$", "", colnames(count.matrix.OG)) == barcodeInfo$barcode[bcIdx] )
# TRUE

# time-points
timepoints = as.matrix(barcodeInfo$timepoint[bcIdx])
rownames(timepoints) <- colnames(count.matrix.OG)
colnames(timepoints) <- "timepoints"
table(timepoints)
# timepoints
# d0  d11   d3   d7  NPC 
# 860 1286 2626 1474  522 

# gender
genders = as.matrix(barcodeInfo$gender[bcIdx])
rownames(genders) <- colnames(count.matrix.OG)
colnames(genders) <- "genders"
table(genders)
# genders
# female   male 
# 4340   2428 

table(timepoints, genders)
# genders
# timepoints female male
# d0     546  314
# d11    760  526
# d3    2204  422
# d7     308 1166
# NPC    522    0

# alleles
alleles = as.matrix(sub("^.*_", "", colnames(count.matrix.OG)))
alleles[alleles=="alt"] = "Cast"
alleles[alleles=="ref"] = "129"
rownames(alleles) <- colnames(count.matrix.OG)
colnames(alleles) <- "alleles"

table(alleles)
# alleles
# 129 Cast 
# 3384 3384 

table(timepoints, genders, alleles)
# , , alleles = 129
# genders
# timepoints female male
# d0     273  157
# d11    380  263
# d3    1102  211
# d7     154  583
# NPC    261    0
# 
# , , alleles = Cast
# genders
# timepoints female male
# d0     273  157
# d11    380  263
# d3    1102  211
# d7     154  583
# NPC    261    0










##################################################
##################################################
##################################################
# 20200516 
# Filter for low counts

totalCoverage = colSums(count.matrix.OG, na.rm=TRUE)

length(totalCoverage)
# 6768

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.0   179.8   793.0  1691.1  2332.2 27954.0 

# Find biased Genes adjusted for bias
totCovDensity <- density(log10(totalCoverage))
maxTCdensity = max(totCovDensity$y)

totCovLHDensity <- density(log10(totalCoverage[totalCoverage < 500]))
lhMaxTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(500) & totCovLHDensity$x>0][which.max(totCovLHDensity$y[totCovLHDensity$x<log10(500) & totCovLHDensity$x>0])]
lhMinTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(500) & totCovLHDensity$x>lhMaxTotCovDensity][which.min(totCovLHDensity$y[totCovLHDensity$x<log10(500) & totCovLHDensity$x>lhMaxTotCovDensity])]
print(c(lhMinTotCovDensity, 10^lhMinTotCovDensity))

PDFfilename = file.path(superMainOutputDir, paste(projectName, ".monocle2results.coverageDistribution.unfiltered.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
plot(totCovDensity, xlab="log10(total UMIs)")
abline(v=lhMinTotCovDensity, lty="dashed", col="red")
text(lhMinTotCovDensity, 0, round(10^lhMinTotCovDensity, 2), pos=4)
dev.off()

# Cells to exclude due to extreme allelic bias
# altLowCellsBoolIdx = totalCoverage > lhMinTotCovDensity
altLowCellsBoolIdx = totalCoverage > 500

count.matrix.OG = count.matrix.OG[,altLowCellsBoolIdx]

dim(count.matrix.OG)
# 156929   3961

totalCoverage = colSums(as.matrix(count.matrix.OG))

length(totalCoverage)
# 3961

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501    1080    2006    2765    3487   27954 

# cell IDs
cell_list.OG <- colnames(count.matrix.OG)

# 20190313
# Remove singlets again
i = 1
singletList = c()
while( i < length(cell_list.OG) ) {
    if(sub("_...$", "", cell_list.OG[i]) != sub("_...$", "", cell_list.OG[i+1])) {
        singletList = c(singletList, i)
        i = i + 1
    } else {
        i = i + 2
    }
}
length(singletList)
#  146

cell_list.OG = cell_list.OG[-singletList]
count.matrix.OG <- count.matrix.OG[,-singletList]
dim(count.matrix.OG)
# 156,929   3902










##################################################
##################################################
##################################################
# Find monoallelic peaks across all/most cells
# *** AND FILTER THEM OUT ***

refBoolIdx = grepl("_ref$", cell_list.OG)
altBoolIdx = grepl("_alt$", cell_list.OG)

if( sum(altBoolIdx) != sum(refBoolIdx) ) {
    cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
    next
}

# Binarize accessibility
count.matrix.OG.B = count.matrix.OG
count.matrix.OG.B[count.matrix.OG.B>0] = 1 

# Peaks coverage
biallelicPeakCoverage = rowSums(count.matrix.OG.B)

zeroCovBoolIdx = biallelicPeakCoverage == 0 
which(zeroCovBoolIdx)
# chr2:31865907-31866195 
# 85217 

# *********************** #
# *** FILTER THEM OUT *** #
count.matrix.OG = count.matrix.OG[!zeroCovBoolIdx,]
# *** FILTER THEM OUT *** #
# *********************** #

# *********************************************************
# After filtering singlet haplocells and some allelically skewed regions:
dim(count.matrix.OG)
# 156895   3902
# *********************************************************


##################################################
# Plot D-score distribution

count.matrix.OG.B = count.matrix.OG.B[!zeroCovBoolIdx,]

# # Peaks covered in only a small proportion of cells
# lowCovBoolIdx = biallelicPeakCoverage / ncol(count.matrix.OG.B) < 0.001 # 0.1 % 
# sum(lowCovBoolIdx);  round(sum(lowCovBoolIdx) / length(lowCovBoolIdx) * 100, 2)
# # [1] 4142
# # [1] 2.65 %

# Allelic totals per peak
refPeakTotals = rowSums( count.matrix.OG.B[, refBoolIdx] )
altPeakTotals = rowSums( count.matrix.OG.B[, altBoolIdx] )

# d-scores
dScores = (refPeakTotals / (refPeakTotals+altPeakTotals)) - 0.5

round(summary(dScores), 3)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.500  -0.105  -0.020  -0.028   0.054   0.500 
# Note slight bias towards Cast

# Find biased peaks adjusted for bias
dScoreDensity <- density(dScores)
maxDSdensity = max(dScoreDensity$y)

dScoreLHDensity <- density(dScores[dScores < 0])
dScoreRHDensity <- density(dScores[dScores > 0])
lhMinDscoreDensity = dScoreLHDensity$x[dScoreLHDensity$x<0 & dScoreLHDensity$x>-0.5][which.min(dScoreLHDensity$y[dScoreLHDensity$x<0 & dScoreLHDensity$x>-0.5])]
rhMinDscoreDensity = dScoreRHDensity$x[dScoreRHDensity$x>0 & dScoreRHDensity$x<0.5][which.min(dScoreRHDensity$y[dScoreRHDensity$x>0 & dScoreRHDensity$x<0.5])]
print(c(lhMinDscoreDensity, rhMinDscoreDensity))
# [1] -0.4757382  0.4739299

# Peaks to exclude due to extreme allelic bias
altBiasedPeakBoolIdx = dScores < lhMinDscoreDensity
refBiasedPeakBoolIdx = dScores > rhMinDscoreDensity

sum(altBiasedPeakBoolIdx); round(sum(altBiasedPeakBoolIdx) / length(altBiasedPeakBoolIdx) * 100, 2)
# [1] 775
# [1] 0.5

sum(refBiasedPeakBoolIdx); round(sum(refBiasedPeakBoolIdx) / length(refBiasedPeakBoolIdx) * 100, 2)
# [1] 387
# [1] 0.25

# plot dScore distribution 
thisOutputDir = file.path(superMainOutputDir, "allelicBiasAnalysis")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = superMainOutputDir
}

PDFfilename = file.path(thisOutputDir, paste(projectName, ".TotalAccessibilityDscoreDistribution.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
ggplot() + aes(dScores) + geom_density(alpha = 0.4, size=0.75) + 
    ggtitle( paste("D-score distribution based on total accessibility per peak across haplocells ( 129/(129+Cast) - 0.5; ", length(dScores), " genes)", sep="") ) +
    # geom_vline(xintercept=mean(dScores), linetype = "dashed", colour = "blue") + 
    # geom_text(aes(x=mean(dScores), label=round(mean(dScores), 2), y=0.25*maxDSdensity), colour="blue", angle=90, vjust = +1) +
    geom_vline(xintercept=lhMinDscoreDensity, linetype = "dotted", colour = "red") + 
    geom_text(aes(x=lhMinDscoreDensity, 
                  label=paste("d-score = ", round(lhMinDscoreDensity, 2),
                              "; # Cast-biased peaks = ", sum(altBiasedPeakBoolIdx), 
                              " (", round(sum(altBiasedPeakBoolIdx) / length(altBiasedPeakBoolIdx) * 100, 2), "%)", 
                              sep=""), y=0.25*maxDSdensity, hjust="left"), 
              colour="red", angle=90, vjust = -1) +
    geom_vline(xintercept=rhMinDscoreDensity, linetype = "dotted", colour = "red") + 
    geom_text(aes(x=rhMinDscoreDensity, 
                  label=paste("d-score = ", round(rhMinDscoreDensity, 2),  
                              "; # 129-biased peaks = ", sum(refBiasedPeakBoolIdx), 
                              " (", round(sum(refBiasedPeakBoolIdx) / length(refBiasedPeakBoolIdx) * 100, 2), "%)", 
                              sep=""), y=0.25*maxDSdensity, hjust="left"), 
              colour="red", angle=90, vjust = -1)
dev.off()










##################################################
##################################################
##################################################
# Generate CDS data

sapply(count.matrix.OG, class)
# as.matrix(count.matrix.OG[1:3,1:3])

# require(Matrix)
count.matrix.OG.sparse = Matrix::Matrix(as.matrix(count.matrix.OG), sparse=TRUE)

locus.annotations = as.data.frame( cbind( rownames(count.matrix.OG), rownames(count.matrix.OG)), stringsAsFactors = FALSE )
colnames(locus.annotations) = c("id", "gene_short_name")
rownames(locus.annotations) = locus.annotations$id

cell.annotations = as.data.frame( cbind( colnames(count.matrix.OG), "sciATAC"), stringsAsFactors = FALSE )
colnames(cell.annotations) = c("cell", "sample")
rownames(cell.annotations) = cell.annotations$cell

pd = new("AnnotatedDataFrame", data = cell.annotations)
fd = new("AnnotatedDataFrame", data = locus.annotations)

# pData(pd)
# pData(fd)

cds.OG = newCellDataSet(count.matrix.OG.sparse, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
# cds.OG = new("CellDataSet", exprs = count.matrix.OG.sparse, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
pData(cds.OG)$n.umi = apply(Biobase::exprs(cds.OG), 2, sum)










##################################################
##################################################
##################################################
# Continue exploring data

# Total expression per barcode
pData(cds.OG)$totalExpr <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!
pData(cds.OG)[1:5,]

# Total expression distributions 
summary(pData(cds.OG)$totalExpr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501    1120    2051    2798    3513   27954 

summary(log10(pData(cds.OG)$totalExpr))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.700   3.049   3.312   3.312   3.546   4.446 










##################################################
##################################################
##################################################
# Add metadata

bcIdx = match( sub("_.*$", "", colnames(count.matrix.OG)), barcodeInfo$barcode )
all( sub("_.*$", "", colnames(count.matrix.OG)) == barcodeInfo$barcode[bcIdx] )
# TRUE

# time-points
timepoints = as.matrix(barcodeInfo$timepoint[bcIdx])
rownames(timepoints) <- colnames(count.matrix.OG)
colnames(timepoints) <- "timepoints"
table(timepoints)
# timepoints
# d0  d11   d3   d7  NPC 
# 560  498 1724  756  364 

# gender
genders = as.matrix(barcodeInfo$gender[bcIdx])
rownames(genders) <- colnames(count.matrix.OG)
colnames(genders) <- "genders"
table(genders)
# female   male 
# 2692   1210 

table(timepoints, genders)
# timepoints female male
# d0     330  230
# d3    1480  244
# d7     160  596
# d11    358  140
# NPC    364    0

all( colnames(cds.OG) == rownames(timepoints) )
# TRUE
all( colnames(cds.OG) == rownames(genders) )
# TRUE

pData(cds.OG)$gender <- factor( genders, levels=c("female", "male") )
pData(cds.OG)$timepoint <- factor( timepoints, levels=c("d0", "d3", "d7", "d11", "NPC") )

# 20200709
# To more easily filter Xist regions
fData(cds.OG)$ensChr <- sub(":.*$", "", fData(cds.OG)$id )
fData(cds.OG)$ensStart <- as.numeric( as.character( sub("-.*$", "", sub("^.*:", "", fData(cds.OG)$id ) ) ) )
fData(cds.OG)$ensEnd <- as.numeric( as.character( sub("^.*-", "", sub("^.*:", "", fData(cds.OG)$id ) ) ) )

# 20200709
# And this too while I'm at it.
pData(cds.OG)$sex <- pData(cds.OG)$gender

# 20200729
# For allelic distrubution plots
cellAlleles = sub(".*_", "", rownames(pData(cds.OG)))
cellAlleles[cellAlleles=="alt"]="Cast"
cellAlleles[cellAlleles=="ref"]="129"
pData(cds.OG)$allele <- factor(cellAlleles, levels=c("129", "Cast"))

dim(cds.OG)
# Features  Samples 
# 156895     3902 

# rm(count.matrix.OG)

# save.image(file=file.path(superMainOutputDir,".RData")) # To .RData
# load(file=file.path(superMainOutputDir,".RData"))











##################################################
##################################################
##################################################
# Normalize and pre-process the data

cds.OG <- estimateSizeFactors(cds.OG)
pData(cds.OG)[1:5,]










##################################################
##################################################
##################################################
# Filter out extreme cells

# Total expression per barcode
pData(cds.OG)$totalCoverage <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!
pData(cds.OG)[1:5,]

# rm(cds)
upper_bound <- 10^(mean(log10(pData(cds.OG)$totalCoverage)) +
                       2*sd(log10(pData(cds.OG)$totalCoverage)))
lower_bound <- 10^(mean(log10(pData(cds.OG)$totalCoverage)) -
                       2*sd(log10(pData(cds.OG)$totalCoverage)))

# Filter
cds <- cds.OG[,pData(cds.OG)$totalCoverage > lower_bound &
                  pData(cds.OG)$totalCoverage < upper_bound]

dim(cds.OG)
# Features  Samples 
# 156895     3902 

# 20190313
# Remove singlets
cell_list = pData(cds)$cell
i = 1
singletList = c()
while( i < length(cell_list) ) {
    if(sub("_...$", "", cell_list[i]) != sub("_...$", "", cell_list[i+1])) {
        singletList = c(singletList, i)
        i = i + 1
    } else {
        i = i + 2
    }
}
length(singletList)
# 11

# cell_list.OG[singletList]
# cell_list.OG[singletList-1]
# cell_list.OG[singletList-2]
# cell_list.OG[singletList+1]
# cell_list.OG[singletList+2]

cds <- cds[,-singletList]

# 20200709
# To more easily filter Xist regions
fData(cds)$ensChr <- sub(":.*$", "", fData(cds)$id )
fData(cds)$ensStart <- as.numeric( as.character( sub("-.*$", "", sub("^.*:", "", fData(cds)$id ) ) ) )
fData(cds)$ensEnd <- as.numeric( as.character( sub("^.*-", "", sub("^.*:", "", fData(cds)$id ) ) ) )
#
# And this too while I'm at it.
pData(cds)$sex <- pData(cds)$gender

dim(cds)
# Features  Samples 
# 156895     3800 










##################################################
##################################################
##################################################
# Once you've excluded cells that do not pass your quality control filters, 
# you should verify that the expression values stored in your CellDataSet 
# follow a distribution that is roughly lognormal: 
p1 <- qplot(totalCoverage, data=pData(cds.OG), geom="density") +
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound)
p2 <- qplot(log10(totalCoverage), data=pData(cds.OG), geom="density") +
    geom_vline(xintercept = log10(lower_bound)) +
    geom_vline(xintercept = log10(upper_bound))
# p3 <- qplot(log10(totalCoverage), data = pData(cds), geom = "density") +
#     geom_vline(xintercept = log10(lower_bound)) +
#     geom_vline(xintercept = log10(upper_bound))
p3 <- ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(fill=gender), alpha = 0.4) 
p4 <- ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(color=timepoint), size=0.75) + scale_color_manual(values = tpCols)

PDFfilename = file.path(superMainOutputDir, paste(projectName, ".monocle2results.coverageDistributions.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
multiplot(p1, p2, p3, p4, layout=matrix(c(1,1,2,2,3,3,4,4), nrow=2, byrow=TRUE))
dev.off()

dim(cds)
# Features  Samples 
# 156895     3800 

# Note: many more d0 cells
table(cds$timepoint, cds$gender) /2 
# female male
# d0     165  114
# d3     705  119
# d7      77  293
# d11    177   70
# NPC    180    0

# 20200914 summary stats
# @ 500 UMIs
cbind( tapply(pData(cds)$totalExpr[pData(cds)$sex == "female" & pData(cds)$allele == "129"], 
              pData(cds)$timepoint[pData(cds)$sex == "female" & pData(cds)$allele == "129"], mean),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "female" & pData(cds)$allele == "Cast"], 
              pData(cds)$timepoint[pData(cds)$sex == "female" & pData(cds)$allele == "Cast"], mean),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "male" & pData(cds)$allele == "129"], 
              pData(cds)$timepoint[pData(cds)$sex == "male" & pData(cds)$allele == "129"], mean),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "male" & pData(cds)$allele == "Cast"], 
              pData(cds)$timepoint[pData(cds)$sex == "male" & pData(cds)$allele == "Cast"], mean),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "female" & pData(cds)$allele == "129"], 
              pData(cds)$timepoint[pData(cds)$sex == "female" & pData(cds)$allele == "129"], median),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "female" & pData(cds)$allele == "Cast"], 
              pData(cds)$timepoint[pData(cds)$sex == "female" & pData(cds)$allele == "Cast"], median),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "male" & pData(cds)$allele == "129"], 
              pData(cds)$timepoint[pData(cds)$sex == "male" & pData(cds)$allele == "129"], median),
       tapply(pData(cds)$totalExpr[pData(cds)$sex == "male" & pData(cds)$allele == "Cast"], 
              pData(cds)$timepoint[pData(cds)$sex == "male" & pData(cds)$allele == "Cast"], median)
)
# [,1]     [,2]     [,3]     [,4] [,5] [,6]   [,7]   [,8]
# d0  2277.012 2309.430 2579.132 2492.044 1734 1706 2286.5 2212.5
# d3  2923.847 2916.567 3415.143 3352.765 2359 2392 2913.0 2891.0
# d7  2580.649 2634.662 2006.287 1960.464 1958 1917 1485.0 1488.0
# d11 2096.079 2103.186 1874.014 1776.257 1794 1699 1179.0 1183.0
# NPC 2186.200 2190.467       NA       NA 1779 1758     NA     NA

# # 20200729
# # For allelic distrubution plots
# cellAlleles.cds = sub(".*_", "", rownames(pData(cds)))
# cellAlleles.cds[cellAlleles.cds=="alt"]="Cast"
# cellAlleles.cds[cellAlleles.cds=="ref"]="129"
# pData(cds)$allele <- factor(cellAlleles.cds, levels=c("129", "Cast"))
#
# Already applied to cds.OG so not necessary here.










##################################################
##################################################
##################################################
# Distribution of cell counts with expression

min_expr = 0.1  # == 0 i.e. >= 1 # see below
# min_expr = 1.1 # i.e. >= 2
# min_expr = 2.1 # i.e. >= 3

cds <- detectGenes(cds, min_expr = min_expr)



##################################################
# Distribution of cell counts region is accessible in:
# # Distribution of cell counts gene is expressed in:
summary(fData(cds)$num_cells_expressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   13.00   22.00   46.27   45.00 1819.00 


# min_expr = 0 == 0.1
summary(fData(cds)$num_cells_expressed) 
summary( Matrix::rowSums(Biobase::exprs(cds) > 0) )
summary( Matrix::rowSums(Biobase::exprs(cds) > 0.1) )
summary( Matrix::rowSums(Biobase::exprs(cds) >= 1) )

summary( Matrix::rowSums(Biobase::exprs(cds) > 1) )
summary( Matrix::rowSums(Biobase::exprs(cds) >=2) )

quantile(p=seq(0.75, 0.95, 0.05), fData(cds)$num_cells_expressed)
# 75% 80% 85% 90% 95% 
# 45  55  73 107 190

cellCounts=c()
for (i in c(2, 5, 10, 20, 25, 50, 100)) {
    cellCounts[as.character(i)] = sprintf( "%s%% of regions accessible in %s+ (%s%%) haplocells",
                                           round( sum(fData(cds)$num_cells_expressed > i) / length(fData(cds)$num_cells_expressed) * 100, 2 ), 
                                           i, round(i / nrow(pData(cds)) * 100, 2) )
}
cellCounts
# 2 
# "99.15% of regions accessible in 2+ (0.05%) haplocells" 
# 5 
# "95.95% of regions accessible in 5+ (0.13%) haplocells" 
# 10 
# "83.14% of regions accessible in 10+ (0.26%) haplocells" 
# 20 
# "53.9% of regions accessible in 20+ (0.53%) haplocells" 
# 25 
# "44.31% of regions accessible in 25+ (0.66%) haplocells" 
# 50 
# "22.09% of regions accessible in 50+ (1.32%) haplocells" 
# 100 
# "10.65% of regions accessible in 100+ (2.63%) haplocells" 



##################################################
# Distribution of proportions of cells gene is expressed in:
summary(fData(cds)$num_cells_expressed) / nrow(pData(cds))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.003421 0.005790 0.012176 0.011842 0.478684 

meanPropotion = summary(fData(cds)$num_cells_expressed)["Mean"]  / nrow(pData(cds))

# Expressed genes in top half of cells:
expressedGenes <-  row.names(subset(fData(cds), 
                                    num_cells_expressed >= ( meanPropotion * nrow(pData(cds)) ) ) )
length(expressedGenes)
# 37835

# top 5% most expressed genes in terms of cells:
expressedGenes <-  row.names(subset(fData(cds), 
                                    num_cells_expressed >= ( 0.05 * nrow(pData(cds)) ) ) )
length(expressedGenes)
# 7905

# 10 cells
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))



##################################################
# Distributions of num_cells_expressed

PDFfilename = file.path(superMainOutputDir, paste(projectName, ".monocle2results.numCellsAccessibleDistributions.pdf", sep=""))
pdf(PDFfilename, height=4*2, width=4*4)
par(mfrow=c(2, 4))

plot( density(log10(fData(cds)$num_cells_expressed + 1)), 
      main="log10 num_cells_accessible + 1", xlab="num cells", ylab="Region count density" )
legend("topright", legend=cellCounts, cex=1.0, bty='n')
legend("right", legend=paste0("No cells: ", sum(fData(cds)$num_cells_expressed<=0), " / ",
                              length(fData(cds)$num_cells_expressed), 
                              " (", round( sum(fData(cds)$num_cells_expressed<=0) / 
                                               length(fData(cds)$num_cells_expressed) * 100, 2), 
                              "% of genes)"), bty='n')
plot( density(log10(fData(cds)$num_cells_expressed)), 
      main="log10 num_cells_expressed", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds)$num_cells_expressed[fData(cds)$num_cells_expressed>0])), 
      main="log10 num_cells_expressed > 0", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds)$num_cells_expressed[fData(cds)$num_cells_expressed>0] + 1)), 
      main="log10 num_cells_expressed > 0 + 1", xlab="num cells", ylab="gene count density" )

plot( density(log10(fData(cds[expressed_genes,])$num_cells_expressed + 1)), 
      main="log10 num_cells_expressed + 1 for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds[expressed_genes,])$num_cells_expressed)), 
      main="log10 num_cells_expressed for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds[expressed_genes,])$num_cells_expressed[fData(cds[expressed_genes,])$num_cells_expressed>0])), 
      main="log10 num_cells_expressed > 0 for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds[expressed_genes,])$num_cells_expressed[fData(cds[expressed_genes,])$num_cells_expressed>0] + 1)), 
      main="log10 num_cells_expressed > 0 + 1 for expressed genes", xlab="num cells", ylab="gene count density" )

dev.off()










# ##################################################
# ##################################################
# ##################################################
# # Filter for expressed genes
# 
# min_expr = 0.1  # == 0 i.e. >= 1 # see below
# # min_expr = 1.1 # i.e. >= 2
# # min_expr = 2.1 # i.e. >= 3
# 
# min_cells = 5 # See analysis above
# # min_cells = 10 # See analysis above
# 
# mainOutputDir = file.path(superMainOutputDir, paste0("minExpr", min_expr, "_minCells", min_cells) )
# if ( ! file.exists( mainOutputDir ) ) {
#     if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
# }
# setwd(mainOutputDir)
# getwd()
# 
# save.image(file=file.path(mainOutputDir,".RData")) # To .RData
# # load(file=file.path(mainOutputDir,".RData"))










##################################################
##################################################
##################################################
# 20200415
# Plot and analyze distributions of absolute differential allelic total expression: *** per chromosome ***
#
# And facet by timepoint and sex
# Filter for expressed genes by timepoint and sex
# Use same ylims
#
# 20200417
# Save XCI status of all cells per time point
#
# 20200709
# Apply thresholds for minimum counts per chromosome 
# And exclude Xist

# Applied above:
# ##################################################
# # Code to make expression code work in this context
# 
# fData(cds.OG)$ensChr <- sub(":.*$", "", fData(cds.OG)$id)
# fData(cds)$ensStart <- as.numeric( as.character( sub("-.*$", "", sub("^.*:", "", fData(cds)$id) ) ) )
# fData(cds)$ensEnd <- as.numeric( as.character(sub("^.*-", "", fData(cds)$id) ) ) 
# 
# pData(cds.OG)$sex <- pData(cds.OG)$gender
# 
# # cds = cds.OG
# 
# ##################################################

sub_min_expr = 0.1; sub_min_cells = 5;
# sub_min_expr = 0.1; sub_min_cells = 10;
# sub_min_expr = 1.1; sub_min_cells = 10; 
# sub_min_expr = 1.1; sub_min_cells = 5; 

whichChroms2analze = list()
whichChroms2analze[["chr1andchrX"]] = c(1,20)
whichChroms2analze[["genome"]] = 1:20

whichChroms = "chr1andchrX"
# whichChroms = "genome" # 20200803 Not run yet

FPRthresh = 10

useMeanThreshold=FALSE
# useMeanThreshold=TRUE

# --------------------------------------------------
# 20200709
# Apply thresholds for minimum counts per chromosome 
min_counts_per_chromosome = 10
# min_counts_per_chromosome = 0

# Exclude Xist
rmXist = TRUE
# rmXist = FALSE
# --------------------------------------------------

library("ggplotify")
library(gridExtra)

modes <- function(d){
    i <- which(diff(sign(diff(d$y))) < 0) + 1
    data.frame(x = d$x[i], y = d$y[i])
}

thisOutputDir = file.path(superMainOutputDir, "differential_allelic_accessibility_distributions",  
                          paste0( whichChroms, 
                                  ifelse(useMeanThreshold, "", ".chromosomewideFPR"), 
                                  ifelse(min_counts_per_chromosome==0, "", 
                                         paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                  ifelse(rmXist, ".rmXist", "") 
                                  ) )
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = currOutputDir
}

alldaLogPlots.hist.scaled.na = list()
alldaLogPlots.hist.scaled.na.new = list()

autosomal.differential.FPRthreshold.vector = c()
autosomal.log2Ratio.FPRthreshold.vector = c()

XCIresultsTable.diff = c()
XCIresultsTable.l2r = c()

XCIcellsPerTimepoint = list()
XCIcellsPerTimepoint[[whichChroms]] = list()

autosomal.differential.FPRthreshold.perChrom = list()
autosomal.log2Ratio.FPRthreshold.perChrom = list()


############
# Genome-wide allelic counts
# For limits

refBoolIdx = grepl("_ref$", rownames(pData(cds)))
altBoolIdx = grepl("_alt$", rownames(pData(cds)))

if( sum(altBoolIdx) != sum(refBoolIdx) ) {
    cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
    next
}

all(pData(cds)$n.umi == pData(cds)$totalExpr)
# [1] TRUE

refExpr = pData(cds)$totalExpr
altExpr = pData(cds)$totalExpr

maxGenomeCount = max(c(refExpr, altExpr), na.rm=TRUE)

YlimsPerChrom.log2 = list()
YlimsPerChrom.log2[["chr1"]] = c(0, 0.25)
YlimsPerChrom.log2[["chrX"]] = c(0, 0.175)
if (whichChroms == "genome" ) {
    for( aChrom in paste0( "chr", 2:19 ) )
        YlimsPerChrom.log2[[aChrom]] = YlimsPerChrom.log2[["chr1"]]
}

XlimsPerChrom.log2 = list()
XlimsPerChrom.log2[["chr1"]] = c(-10, 10)
XlimsPerChrom.log2[["chrX"]] = c(-10, 10)
if (whichChroms == "genome" ) {
    for( aChrom in paste0( "chr", 2:19 ) )
        XlimsPerChrom.log2[[aChrom]] = XlimsPerChrom.log2[["chr1"]]
}

# --------------------------------------------------
# 20200709
# Apply thresholds for minimum counts per chromosome 
# For each chromosome 
#
# And exclude Xist

if (min_counts_per_chromosome > 0) {
    
    ###
    # Get cells showing minimum total chromosomal expression per chromosome
    
    min_exprC_boolidx_per_chrom = list()
    
    for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
        # for (theChrom in theChroms) {
        
        if (rmXist) {
            
            # Exclude Xist
            # From sciRNA-seq CDS:
            # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
            # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
            # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
            
            Xist_window = 5000
            Xist_region_start = 103460375 - Xist_window
            Xist_region_end = 103483217 + Xist_window
            # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
            # [1] 8
            
            chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
        } else {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
        }
        cdsC.OG <- cds[chromgenesNames,]
        print(dim(cdsC.OG))
        
        ############
        # chromosome-wide allelic counts
        
        pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
        
        refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
        altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
        
        if( sum(altBoolIdx) != sum(refBoolIdx) ) {
            cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
            next
        }
        
        refExprC = pData(cdsC.OG)$Cexpr[refBoolIdx]
        altExprC = pData(cdsC.OG)$Cexpr[altBoolIdx]
        
        totalExprC = refExprC + altExprC
        
        min_exprC_boolidx_per_chrom[[ theChrom ]] = min_exprC_boolidx = totalExprC >= min_counts_per_chromosome
        
        cat( theChrom, ": ", sum(min_exprC_boolidx_per_chrom[[ theChrom ]]), " / ", 
             length(min_exprC_boolidx_per_chrom[[ theChrom ]]),
             " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
             sep="" )
    }
    # w/o Xist:
    # Features  Samples 
    # 4301     3800 
    # chrX: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)

    # 20200803   
    # Only use chr1 and chrX counts as threshold
    if ( whichChroms == "genome" ){
        min_exprC_boolidx = min_exprC_boolidx_per_chrom[[ "chr1" ]] # intialize
        for (theChrom in "chrX") {
            # for (theChrom in theChroms) {
            print(theChrom)
            
            min_exprC_boolidx = min_exprC_boolidx & min_exprC_boolidx_per_chrom[[ theChrom ]]
        }
    } else {
        
        ###
        # Now consider cells with sufficient chromosomal counts for all chromosomes under consideration
        min_exprC_boolidx = min_exprC_boolidx_per_chrom[[ 1 ]] # intialize
        for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
            # for (theChrom in theChroms) {
            print(theChrom)
            
            min_exprC_boolidx = min_exprC_boolidx & min_exprC_boolidx_per_chrom[[ theChrom ]] 
        }
    }
    
    cat( "All chromosomes: ", sum(min_exprC_boolidx), " / ", 
         length(min_exprC_boolidx), 
         " cells pass filter", 
         paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
         sep="" )
    # w/o Xist:
    # All chromosomes: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)
    
    ###
    # Check min_exprC_boolidx how many female cells per time point?
    
    # all( grepl("_ref$", rownames(pData(cdsC.OG))) == grepl("_ref$", rownames(pData(cds))) )
    # TRUE
    refBoolIdx = grepl("_ref$", rownames(pData(cds)))
    altBoolIdx = grepl("_alt$", rownames(pData(cds)))
    
    for (aSex in c("female", "male")) {
        # aSex = "female"
        
        subsexBoolidx = pData(cds)$sex == aSex
        # all( pData(cdsC.OG)$sex == pData(cds)$sex )
        # [1] TRUE
        
        #for (aTP in c("all", names(tpCols))) {
        for (aTP in c(names(tpCols))) {
            # aTP = "NPC"
            
            print(c(aSex, aTP))
            
            if (aTP == "all") {
                tpSubsetBoolidx = rep(TRUE, length(pData(cds)$timepoint))
            } else {
                tpSubsetBoolidx = pData(cds)$timepoint == aTP
            }
            
            subsetBoolidx = tpSubsetBoolidx & subsexBoolidx
            
            # cds.TP.Sex <- cds[, subsetBoolidx]
            # print(dim(cdsC.TP.Sex))
            
            # For cells not haplocells
            # all( subsetBoolidx[ refBoolIdx ] == subsetBoolidx[ altBoolIdx ] )
            # TRUE
            
            cat( aTP, ": ", sum( min_exprC_boolidx & subsetBoolidx[ refBoolIdx ] ), " / ", 
                 sum( subsetBoolidx[ refBoolIdx ] ), 
                 " ", aTP, " ", aSex, " cells pass filter", 
                 paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
                 sep="" )
        }
        
        cat( sum( min_exprC_boolidx & subsexBoolidx[ refBoolIdx ] ), " / ", 
             sum( subsexBoolidx[ refBoolIdx ] ), 
             " ", aSex, " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
             sep="" )
    }
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 164 / 165 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 702 / 705 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 77 / 77 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 174 / 177 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 178 / 180 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1295 / 1304 female cells pass filter (min_counts_per_chromosome >= 10)
    # 
    # [1] "male" "d0"  
    # d0: 113 / 114 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 116 / 119 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 287 / 293 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 67 / 70 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 0 / 0 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 583 / 596 male cells pass filter (min_counts_per_chromosome >= 10)
}

# --------------------------------------------------

############
# For each chromosome 
for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
    # for (theChrom in theChroms) {
    print(theChrom)
    
    daLogPlots.hist.scaled.na = list()
    daLogPlots.hist.scaled.na.new = list()
    
    plotally = 0
    
    # --------------------------------------------------
    # 20200709
    # Exclude Xist

    if (rmXist) {
        
        # Exclude Xist
        # From sciRNA-seq CDS:
        # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
        # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
        # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
        
        Xist_window = 5000
        Xist_region_start = 103460375 - Xist_window
        Xist_region_end = 103483217 + Xist_window
        # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
        # [1] 8
        
        chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
    } else {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
    }
    cdsC.OG <- cds[chromgenesNames,]
    print(dim(cdsC.OG))
    
    # --------------------------------------------------
    
    ############
    # chromosome-wide allelic counts
    # Not currently used.
    #
    # 20200519 Now chromosome-wide FPR threshold is an option.
    
    pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
    
    refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
    altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
    
    if( sum(altBoolIdx) != sum(refBoolIdx) ) {
        cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
        next
    }
    
    refExprC = pData(cdsC.OG)$Cexpr[refBoolIdx]
    altExprC = pData(cdsC.OG)$Cexpr[altBoolIdx]
    
    maxChromCount = max(c(refExprC, altExprC), na.rm=TRUE)
    
    diffAllelicExp = refExprC - altExprC
    absDiffAllelicExp = abs(diffAllelicExp)
    diffDensity <- density(absDiffAllelicExp)
    maxX.diff= max(diffDensity$x, na.rm=TRUE)
    maxY.diff.density = max(diffDensity$y, na.rm=TRUE) * 1.25
    diffHist <- hist(absDiffAllelicExp, plot=FALSE)
    maxY.diff.hist = max(diffHist$counts, na.rm=TRUE) / 5 
    
    log2AllelicExp = log2(refExprC + 0.1) - log2(altExprC + 0.1)
    absLog2AllelicExp = abs(log2AllelicExp)
    l2rDensity <- density(absLog2AllelicExp)
    maxX.l2r= max(l2rDensity$x, na.rm=TRUE)
    maxY.l2r.density = max(l2rDensity$y, na.rm=TRUE) * 1.25
    l2rHist <- hist(absLog2AllelicExp, plot=FALSE)
    maxY.l2r.hist = max(l2rHist$counts, na.rm=TRUE) / 5 
    
    # 20200519 Now chromosome-wide FPR threshold is an option.
    autosomal.differential.FPRthreshold.perChrom[[theChrom]] =  quantile(absDiffAllelicExp, (100-FPRthresh)/100, na.rm=TRUE)
    autosomal.log2Ratio.FPRthreshold.perChrom[[theChrom]] =  quantile(absLog2AllelicExp, (100-FPRthresh)/100, na.rm=TRUE)
    
    for (aSex in c("female", "male")) {
        # aSex = "female"
        
        if(theChrom =="chrX" && aSex == "male") {
            next
        }
        
        subsexBoolidx = pData(cdsC.OG)$sex == aSex
        
        #for (aTP in c("all", names(tpCols))) {
        for (aTP in c(names(tpCols))) {
            # aTP = "NPC"
            
            print(c(theChrom, aSex, aTP))
            
            plotally = plotally + 1
            
            if (aTP == "all") {
                tpSubsetBoolidx = rep(TRUE, length(pData(cdsC.OG)$timepoint))
            } else {
                tpSubsetBoolidx = pData(cdsC.OG)$timepoint == aTP
            }
            
            subsetBoolidx = tpSubsetBoolidx & subsexBoolidx
            
            cdsC.TP.Sex <- cdsC.OG[, subsetBoolidx]
            print(dim(cdsC.TP.Sex))
            
            cdsC.TP.Sex <- detectGenes(cdsC.TP.Sex, min_expr = sub_min_expr)
            
            # The vector expressed_genes now holds the identifiers for genes expressed in at least n cells of the data set.
            expressed_genes <- row.names(subset(fData(cdsC.TP.Sex),
                                                num_cells_expressed >= sub_min_cells))
            print(length(expressed_genes))
            
            cdsC.TP.Sex <- cdsC.TP.Sex[expressed_genes,]
            print(dim(cdsC.TP.Sex))
            
            pData(cdsC.TP.Sex)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.TP.Sex))  
            
            refBoolIdx = grepl("_ref$", rownames(pData(cdsC.TP.Sex)))
            altBoolIdx = grepl("_alt$", rownames(pData(cdsC.TP.Sex)))
            
            if( sum(altBoolIdx) != sum(refBoolIdx) ) {
                cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
                next
            }

            # --------------------------------------------------
            # 20200709
            # Apply thresholds for minimum counts per chromosome 
            # For each chromosome 
            
            if ( min_counts_per_chromosome > 0 && sum(refBoolIdx) > 0 ) {
                
                all.refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
                all.altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
                
                if( sum(all.altBoolIdx) != sum(all.refBoolIdx) ) {
                    cat( sum(all.altBoolIdx), " != ", sum(all.refBoolIdx), "\n", sep="")
                    next
                }
                
                cat( aTP, ": ", sum( min_exprC_boolidx & subsetBoolidx[ all.refBoolIdx ] ), " / ", 
                     sum( subsetBoolidx[ all.refBoolIdx ] ), 
                     " ", aTP, " ", aSex, " cells pass filter", 
                     paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
                     sep="" )
                
                # length(all.refBoolIdx)
                # [1] 5392
                #
                # length(subsetBoolidx)
                # [1] 5392
                # 
                # all( subsetBoolidx[ all.refBoolIdx ] == subsetBoolidx[ all.altBoolIdx ] )
                # TRUE
                #
                # length(refBoolIdx); length(grepl("_ref$", rownames(pData(cdsC.TP.Sex))))
                # [1] 1080
                #
                # length(refBoolIdx[ seq(2, length(refBoolIdx), 2) ])
                # [1] 540
                #
                # length(min_exprC_boolidx[ subsetBoolidx[all.refBoolIdx] ])
                # [1] 540
                #
                # all(min_exprC_boolidx[ subsetBoolidx[all.refBoolIdx] ]) == all( min_exprC_boolidx[ subsetBoolidx[all.altBoolIdx] ])
                # [1] TRUE
                # 
                # all( (min_exprC_boolidx & subsetBoolidx[ all.refBoolIdx ])[ subsetBoolidx[ all.refBoolIdx ] ]
                # == min_exprC_boolidx[ subsetBoolidx[ all.refBoolIdx ] ] )
                # [1] TRUE
                #
                # min_exprC_boolidx is based on all cells regardless of sex and time point,
                # but applies to both ref and alt haplocells,
                # so must select min_exprC_boolidx entries that are sex and time point specific 
                # and apply to every alternate refBoolIdx / altBoolIdx entry.
                refBoolIdx[ seq(2, length(refBoolIdx), 2) ][ ! min_exprC_boolidx[ subsetBoolidx[ all.refBoolIdx ] ] ] = FALSE
                altBoolIdx[ seq(1, length(altBoolIdx), 2) ][ ! min_exprC_boolidx[ subsetBoolidx[ all.altBoolIdx ] ] ] = FALSE
                if( sum(altBoolIdx) != sum(refBoolIdx) ) {
                    cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
                    next
                }
                
            }
            
            # --------------------------------------------------

            refExprC = pData(cdsC.TP.Sex)$Cexpr[refBoolIdx]
            altExprC = pData(cdsC.TP.Sex)$Cexpr[altBoolIdx]
            
            pData(cdsC.TP.Sex)$totalExpr[refBoolIdx] = refExprC
            pData(cdsC.TP.Sex)$totalExpr[altBoolIdx] = altExprC
            
            
            #####################################
            # Absolute differential allelic accessibility
            
            diffAllelicExp = refExprC - altExprC
            absDiffAllelicExp = abs(diffAllelicExp)
            
            # Save absolute differential scores
            # When chrX: obtain mean value and select XCI cells
            if (theChrom != "chrX") {
                autosomal.differential.FPRthreshold.vector = c( autosomal.differential.FPRthreshold.vector,
                                                                quantile(absDiffAllelicExp, 0.9, na.rm=TRUE) )
            } else {
                # hist(autosomal.differential.FPRthreshold.vector)
                # summary(autosomal.differential.FPRthreshold.vector)
                autosomal.differential.FPRthreshold.mean = mean(autosomal.differential.FPRthreshold.vector, na.rm=TRUE)
                
                # 20200519 Now chromosome-wide FPR threshold is an option.
                if (useMeanThreshold) {
                    XCIboolIdx = absDiffAllelicExp > autosomal.differential.FPRthreshold.mean
                } else {
                    XCIboolIdx = absDiffAllelicExp > autosomal.differential.FPRthreshold.perChrom[["chr1"]]
                }
                
                theLabel = paste(FPRthresh, "% autosomal FPR = ", round(autosomal.differential.FPRthreshold.mean, 2),
                                 "\n# inactive cells =  ", sum(XCIboolIdx),
                                 " (", round(sum(XCIboolIdx) / length(XCIboolIdx) * 100, 2), "%)",
                                 sep="")
                print(c("delta: ", theLabel))
                
                XCIresultsRow.diff = c( autosomal.differential.FPRthreshold.mean, 
                                        sum(XCIboolIdx),
                                        length(XCIboolIdx),
                                        sum(XCIboolIdx) / length(XCIboolIdx))
                XCIresultsTable.diff = cbind(XCIresultsTable.diff,
                                             XCIresultsRow.diff)
            }
            
            pData(cdsC.TP.Sex)$DAexpr = rep(0, length(refBoolIdx))
            pData(cdsC.TP.Sex)$DAexpr[refBoolIdx] = diffAllelicExp
            pData(cdsC.TP.Sex)$DAexpr[altBoolIdx] = diffAllelicExp
            
            pData(cdsC.TP.Sex)$ADAexpr = rep(0, length(refBoolIdx))
            pData(cdsC.TP.Sex)$ADAexpr[refBoolIdx] = absDiffAllelicExp
            pData(cdsC.TP.Sex)$ADAexpr[altBoolIdx] = absDiffAllelicExp
            
            theXlims = c(0, maxX.diff)
            # theXlims = NULL
            
            # theYlims = c(0, maxY.diff.density)
            theYlims = NULL
            maxY.diff.density = 0
            maxY.diff.hist = 0
            
            
            ####################################
            ## Absolute log2 ratio of allelic accessibility
            
            log2AllelicExp = log2(refExprC + 0.1) - log2(altExprC + 0.1)
            absLog2AllelicExp = abs(log2AllelicExp)
            
            # Save absolute differential scores
            # When chrX: obtain mean value and select XCI cells
            if (theChrom != "chrX") {
                autosomal.log2Ratio.FPRthreshold.vector = c( autosomal.log2Ratio.FPRthreshold.vector,
                                                             quantile(absLog2AllelicExp, (100-FPRthresh)/100, na.rm=TRUE) )
            } else {
                # hist(autosomal.log2Ratio.FPRthreshold.vector)
                # summary(autosomal.log2Ratio.FPRthreshold.vector)
                autosomal.log2Ratio.FPRthreshold.mean = mean(autosomal.log2Ratio.FPRthreshold.vector, na.rm=TRUE)
                
                # 20200519 Now chromosome-wide FPR threshold is an option.
                if (useMeanThreshold) {
                    XCIboolIdx = absLog2AllelicExp > autosomal.log2Ratio.FPRthreshold.mean
                } else {
                    XCIboolIdx = absLog2AllelicExp > autosomal.log2Ratio.FPRthreshold.perChrom[["chr1"]]
                }
                
                theLabel = paste(FPRthresh, "% autosomal FPR = ", round(autosomal.log2Ratio.FPRthreshold.mean, 2),
                                 "\n# inactive cells =  ", sum(XCIboolIdx),
                                 " (", round(sum(XCIboolIdx) / length(XCIboolIdx) * 100, 2), "%)",
                                 sep="")
                print(c("l2r: ", theLabel))
                
                XCIresultsRow.l2r = c( autosomal.log2Ratio.FPRthreshold.mean, 
                                       sum(XCIboolIdx),
                                       length(XCIboolIdx),
                                       sum(XCIboolIdx) / length(XCIboolIdx))
                XCIresultsTable.l2r = cbind(XCIresultsTable.l2r,
                                            XCIresultsRow.l2r)
                
                # > head(pData(cdsC.TP.Sex)$Cexpr[refBoolIdx][XCIboolIdx])
                # [1] 2 1 1 1 0 0
                # > head(refExprC[XCIboolIdx])
                # [1] 2 1 1 1 0 0
                
                XCIcellsPerTimepoint[[whichChroms]][[aTP]] = cbind( refExprC,
                                                                    altExprC,
                                                                    log2AllelicExp,
                                                                    autosomal.log2Ratio.FPRthreshold.mean,
                                                                    0)
                rownames(XCIcellsPerTimepoint[[whichChroms]][[aTP]]) = sub("_ref", "", rownames(pData(cdsC.TP.Sex))[refBoolIdx])
                colnames(XCIcellsPerTimepoint[[whichChroms]][[aTP]]) = c("refExpr", "altExpr", "log2Expr", "auto10FPRthreshold", "XCI")
                XCIcellsPerTimepoint[[whichChroms]][[aTP]][,"XCI"][XCIboolIdx] = 1
            }
            
            pData(cdsC.TP.Sex)$LAexpr = rep(0, length(refBoolIdx))
            pData(cdsC.TP.Sex)$LAexpr[refBoolIdx] = log2AllelicExp
            pData(cdsC.TP.Sex)$LAexpr[altBoolIdx] = log2AllelicExp
            
            pData(cdsC.TP.Sex)$ALAexpr = rep(0, length(refBoolIdx))
            pData(cdsC.TP.Sex)$ALAexpr[refBoolIdx] = absLog2AllelicExp
            pData(cdsC.TP.Sex)$ALAexpr[altBoolIdx] = absLog2AllelicExp
            
            # Same for all plots
            # theXlims = c(0, 5)
            # theYlims = c(0, maxY.l2r.density)
            
            # Different for all plots
            # theXlims = NULL
            theXlims = c(0, maxX.l2r)
            
            theYlims = NULL
            maxY.l2r.density = 0
            maxY.l2r.hist = 0
            
            theXlims = c(-maxX.l2r, maxX.l2r)
            
            ##########            
            # Scaled non-absolute histograms
            # https://stackoverflow.com/questions/11766856/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion
            # https://stackoverflow.com/questions/16339204/normalizing-faceted-histograms-separately-in-ggplot2
            # https://stackoverflow.com/questions/4725339/percentage-on-y-lab-in-a-faceted-ggplot-barchart
            
            # Tried to get max values for plot text:
            # # https://stackoverflow.com/questions/58785930/r-find-maximum-of-density-plot
            # d <- hist(pData(cdsC.TP.Sex)[refBoolIdx,"ALAexpr"], plot=FALSE)
            # maxY.l2r.hist = max(d$density, na.rm=TRUE)
            # theYlims = c(0, maxY.l2r.hist)
            
            daLogPlot.hist.scaled.na <- ggplot(pData(cdsC.TP.Sex)[refBoolIdx,], aes(x=LAexpr)) +
                # geom_histogram(aes(y=..count../sum(..count..)), binwidth = 0.25, center=0, alpha = 0.4, size=0.75) +
                geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth = 0.25, center=0, alpha = 0.4, size=0.75) +
                ggtitle( paste(theChrom, ": ", aSex, " filtered cells (n = ", dim(cdsC.TP.Sex[,refBoolIdx])[2], "; ",
                               " genes = ", dim(cdsC.TP.Sex[,refBoolIdx])[1], ")", sep="") ) +
                # ylab("Proportion of cells") +
                # xlab("log2 ratio of allelic accessibility") +
                theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.y=element_blank()) +
                coord_cartesian(xlim = theXlims, ylim = theYlims)
            
            if (theChrom == "chrX") {
                daLogPlot.hist.scaled.na <- daLogPlot.hist.scaled.na +
                    geom_vline(xintercept=-autosomal.log2Ratio.FPRthreshold.mean, linetype = "dotted", colour = "red") +
                    geom_vline(xintercept=autosomal.log2Ratio.FPRthreshold.mean, linetype = "dotted", colour = "red") +
                    geom_text(aes(x=autosomal.log2Ratio.FPRthreshold.mean, y=0.90*maxY.l2r.hist, hjust="left",
                                  label=theLabel),
                              # colour="red", angle=90, vjust = -1)
                              colour="red", angle=0, vjust = -1)
            }
            
            daLogPlots.hist.scaled.na[[paste0("p", as.character(plotally))]] = as.grob(daLogPlot.hist.scaled.na)
            
            
            ##########            
            # Scaled non-absolute histograms - new fancy version
            # https://stackoverflow.com/questions/11766856/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion
            # https://stackoverflow.com/questions/16339204/normalizing-faceted-histograms-separately-in-ggplot2
            # https://stackoverflow.com/questions/4725339/percentage-on-y-lab-in-a-faceted-ggplot-barchart
            
            # Tried to get max values for plot text:
            # # https://stackoverflow.com/questions/58785930/r-find-maximum-of-density-plot
            # d <- hist(pData(cdsC.TP.Sex)[refBoolIdx,"ALAexpr"], plot=FALSE)
            # maxY.l2r.hist = max(d$density, na.rm=TRUE)
            # theYlims = c(0, maxY.l2r.hist)
            theYlims = YlimsPerChrom.log2[[theChrom]]
            
            # theXlims = c(-maxX.l2r, maxX.l2r)
            theXlims = XlimsPerChrom.log2[[theChrom]]
            
            daLogPlot.hist.scaled.na.new <- ggplot(pData(cdsC.TP.Sex)[refBoolIdx,], aes(x=LAexpr)) +
                # geom_histogram(aes(y=..count../sum(..count..)), binwidth = 0.25, center=0, alpha = 0.4, size=0.75) +
                geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth = 0.25, 
                               center=0, alpha = 0.8, size=0.75, fill = tpCols[aTP]) +
                theme_minimal() +
                ggtitle( paste(theChrom, ": ", aSex, " filtered cells (n = ", dim(cdsC.TP.Sex[,refBoolIdx])[2], "; ",
                               " genes = ", dim(cdsC.TP.Sex[,refBoolIdx])[1], ")", sep="") ) +
                # ylab("Proportion of cells") +
                # xlab("log2 ratio of allelic accessibility") +
                theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.y=element_blank()) +
                coord_cartesian(xlim = theXlims, ylim = theYlims)
            
            if (theChrom == "chrX") {
                daLogPlot.hist.scaled.na.new <- daLogPlot.hist.scaled.na.new +
                    geom_vline(xintercept=-autosomal.log2Ratio.FPRthreshold.mean, linetype = "dashed", colour = "gray", lwd=1.5) +
                    geom_vline(xintercept=autosomal.log2Ratio.FPRthreshold.mean, linetype = "dashed", colour = "gray", lwd=1.5)
                # geom_text(aes(x=autosomal.log2Ratio.FPRthreshold.mean, y=0.90*maxY.l2r.hist, hjust="left",
                #               label=theLabel),
                #           # colour="red", angle=90, vjust = -1)
                #           colour="red", angle=0, vjust = -1)
            }
            
            daLogPlots.hist.scaled.na.new[[paste0("p", as.character(plotally))]] = as.grob(daLogPlot.hist.scaled.na.new)
            
            # rm(cdsC.TP.Sex)
            
        } # timepoint
        
    } # sex
    
    alldaLogPlots.hist.scaled.na[[theChrom]] <- marrangeGrob(grobs = daLogPlots.hist.scaled.na, layout_matrix=matrix(1:10, nrow=2, byrow=TRUE),
                                                             top = paste(theChrom, " scaled histograms of log2 differential accessibility based on total accessibility per cell across ", sum(subsetBoolidx), " cells for expressed genes" ) )
    
    alldaLogPlots.hist.scaled.na.new[[theChrom]] <- marrangeGrob(grobs = daLogPlots.hist.scaled.na.new, layout_matrix=matrix(1:10, nrow=2, byrow=TRUE),
                                                                 top = paste(theChrom, " scaled histograms of log2 differential expression based on total expression per cell across ", sum(subsetBoolidx), " cells for expressed genes" ) )
    
} # chrom


#########
# Tabulate results

colnames(XCIresultsTable.diff) = names(tpCols)
colnames(XCIresultsTable.l2r) = names(tpCols)
rownames(XCIresultsTable.diff) = c("threshold", "XCI_counts", "total_count", "proportion")
rownames(XCIresultsTable.l2r) = c("threshold", "XCI_counts", "total_count", "proportion")

tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCItally.byDifferentialAllelicTotalAccessibility.", FPRthresh, ".byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                               ".", whichChroms, ".tsv", sep=""))
write.table(x = XCIresultsTable.diff, file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCItally.byLogRatioAllelicTotalAccessibility.", FPRthresh, ".byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                               ".", whichChroms, ".tsv", sep=""))
write.table(x = XCIresultsTable.l2r, file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


#########
# XCI cell table

for (aTP in names(XCIcellsPerTimepoint[[whichChroms]])) {
    tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCIcells.", aTP, ".byLog2AllelicAccessibility.", FPRthresh, ".byTimePoint", 
                                                   ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                   ".", whichChroms, ".tsv", sep=""))
    write.table(x = XCIcellsPerTimepoint[[whichChroms]][[aTP]], file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
}


#########
# XCI cell table - all cells

XCIcellsAll = list()
XCIcellsAll[[whichChroms]] = c()
for (aTP in names(XCIcellsPerTimepoint[[whichChroms]])) {
    XCIcellsAll[[whichChroms]] = rbind( XCIcellsAll[[whichChroms]],
                                        cbind( XCIcellsPerTimepoint[[whichChroms]][[aTP]], aTP ) )
}
colnames(XCIcellsAll[[whichChroms]]) = c("refExpr", "altExpr", "log2Expr", "auto10FPRthreshold", "XCI", "timepoint" )  
dim(XCIcellsAll[[whichChroms]])
#  1304    6

# save
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCIcells.all.byLog2AllelicAccessibilityn.", FPRthresh, "pctFPR.byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                               ".", whichChroms, ".tsv", sep=""))
write.table(x = XCIcellsAll[[whichChroms]], file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


#########
# plot

require(gridExtra)
require(grid)


# Scaled histograms of log2 ratio of  allelic accessibility per chromosome
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.V2.scaled.histogramsOfLog2AllelicAccessibilityPerChromosome.faceted",
                                             ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                             ".", whichChroms, ".pdf", sep=""))
pdf(PDFfilename, height=2*5, width=5*5)
for (aGrob in names(alldaLogPlots.hist.scaled.na)) {
    print(aGrob)
    print(alldaLogPlots.hist.scaled.na[[aGrob]])
}
dev.off()

# Fancy scaled histograms of log2 ratio of  allelic expression per chromosome
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.V2.scaled.fancyHistogramsOfLog2AllelicAccessibilityPerChromosome.faceted",
                                             ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                             ".", whichChroms, ".pdf", sep=""))
pdf(PDFfilename, height=2*5, width=5*5)
for (aGrob in names(alldaLogPlots.hist.scaled.na.new)) {
    print(aGrob)
    print(alldaLogPlots.hist.scaled.na.new[[aGrob]])
}
dev.off()










##################################################
##################################################
##################################################
# 20200506
# Save concatenated allelic UMIs for F121s
#
# Filter for expressed genes across timepoints and sexes
#
# 20200709
# Apply thresholds for minimum counts per chromosome 
# And allow for exclusion of Xist... but actually keep it in this context.

alleleCDPmatrixDir = file.path("../data", "allelic_sci-ATAC-seq_catenated_count_matrices")
if ( ! file.exists( alleleCDPmatrixDir ) ) {
    if ( ! dir.create( alleleCDPmatrixDir, showWarnings=FALSE, recursive=TRUE )  ) alleleCDPmatrixDir = superMainOutputDir
}

# F121libs=c(
#     "F121_d0",
#     "F121_d3",
#     "F121_d7",
#     "F121_d11",
#     "F121_NPC"
# )

# Applied above:
# ##################################################
# # Code to make expression code work in this context
# 
# fData(cds.OG)$ensChr <- sub(":.*$", "", fData(cds.OG)$id)
# fData(cds)$ensStart <- as.numeric( as.character( sub("-.*$", "", sub("^.*:", "", fData(cds)$id) ) ) )
# fData(cds)$ensEnd <- as.numeric( as.character(sub("^.*-", "", fData(cds)$id) ) ) 
# 
# pData(cds.OG)$sex <- pData(cds.OG)$gender
# 
# # cds = cds.OG
# 
# ##################################################

sub_min_expr = 0.1; sub_min_cells = 5;
# sub_min_expr = 0.1; sub_min_cells = 10;
# sub_min_expr = 1.1; sub_min_cells = 10; 
# sub_min_expr = 1.1; sub_min_cells = 5; 

whichChroms2analze = list()
whichChroms2analze[["chr1andchrX"]] = c(1,20)
whichChroms2analze[["genome"]] = 1:20

whichChroms = "chr1andchrX"
# whichChroms = "genome" 

# --------------------------------------------------
# 20200709
# Apply thresholds for minimum counts per chromosome 
# For each chromosome 
#
# And exclude Xist

min_counts_per_chromosome = 10
# min_counts_per_chromosome = 0

# Exclude Xist
# rmXist = TRUE
rmXist = FALSE

if (min_counts_per_chromosome > 0) {
    
    ###
    # Get cells showing minimum total chromosomal expression per chromosome
    
    min_exprC_boolidx_per_chrom = list()
    
    for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
        # for (theChrom in theChroms) {
        
        if (rmXist) {
            
            # Exclude Xist
            # From sciRNA-seq CDS:
            # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
            # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
            # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
            
            Xist_window = 5000
            Xist_region_start = 103460375 - Xist_window
            Xist_region_end = 103483217 + Xist_window
            # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
            # [1] 8
            
            chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
        } else {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
        }
        cdsC.OG <- cds[chromgenesNames,]
        print(dim(cdsC.OG))
        
        ############
        # chromosome-wide allelic counts
        
        pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
        
        refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
        altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
        
        if( sum(altBoolIdx) != sum(refBoolIdx) ) {
            cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
            next
        }
        
        refExprC = pData(cdsC.OG)$Cexpr[refBoolIdx]
        altExprC = pData(cdsC.OG)$Cexpr[altBoolIdx]
        
        totalExprC = refExprC + altExprC
        
        min_exprC_boolidx_per_chrom[[ theChrom ]] = min_exprC_boolidx = totalExprC >= min_counts_per_chromosome
        
        cat( theChrom, ": ", sum(min_exprC_boolidx_per_chrom[[ theChrom ]]), " / ", 
             length(min_exprC_boolidx_per_chrom[[ theChrom ]]),
             " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
             sep="" )
    }
    # w/o Xist:
    # Features  Samples 
    # 4301     3800 
    # chrX: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)
    
    # 20200803   
    # Only use chr1 and chrX counts as threshold
    if ( whichChroms == "genome" ){
        min_exprC_boolidx = min_exprC_boolidx_per_chrom[[ "chr1" ]] # intialize
        for (theChrom in "chrX") {
            # for (theChrom in theChroms) {
            print(theChrom)
            
            min_exprC_boolidx = min_exprC_boolidx & min_exprC_boolidx_per_chrom[[ theChrom ]]
        }
    } else {
        
        ###
        # Now consider cells with sufficient chromosomal counts for all chromosomes under consideration
        min_exprC_boolidx = min_exprC_boolidx_per_chrom[[ 1 ]] # intialize
        for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
            # for (theChrom in theChroms) {
            print(theChrom)
            
            min_exprC_boolidx = min_exprC_boolidx & min_exprC_boolidx_per_chrom[[ theChrom ]] 
        }
    }
    
    cat( "All chromosomes: ", sum(min_exprC_boolidx), " / ", 
         length(min_exprC_boolidx), 
         " cells pass filter", 
         paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
         sep="" )
    # w/o Xist:
    # All chromosomes: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)
    
    ###
    # Check min_exprC_boolidx how many female cells per time point?
    
    # all( grepl("_ref$", rownames(pData(cdsC.OG))) == grepl("_ref$", rownames(pData(cds))) )
    # TRUE
    refBoolIdx = grepl("_ref$", rownames(pData(cds)))
    altBoolIdx = grepl("_alt$", rownames(pData(cds)))
    
    for (aSex in c("female", "male")) {
        # aSex = "female"
        
        subsexBoolidx = pData(cds)$sex == aSex
        # all( pData(cdsC.OG)$sex == pData(cds)$sex )
        # [1] TRUE
        
        #for (aTP in c("all", names(tpCols))) {
        for (aTP in c(names(tpCols))) {
            # aTP = "NPC"
            
            print(c(aSex, aTP))
            
            if (aTP == "all") {
                tpSubsetBoolidx = rep(TRUE, length(pData(cds)$timepoint))
            } else {
                tpSubsetBoolidx = pData(cds)$timepoint == aTP
            }
            
            subsetBoolidx = tpSubsetBoolidx & subsexBoolidx
            
            # cds.TP.Sex <- cds[, subsetBoolidx]
            # print(dim(cdsC.TP.Sex))
            
            # For cells not haplocells
            # all( subsetBoolidx[ refBoolIdx ] == subsetBoolidx[ altBoolIdx ] )
            # TRUE
            
            cat( aTP, ": ", sum( min_exprC_boolidx & subsetBoolidx[ refBoolIdx ] ), " / ", 
                 sum( subsetBoolidx[ refBoolIdx ] ), 
                 " ", aTP, " ", aSex, " cells pass filter", 
                 paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
                 sep="" )
        }
        
        cat( sum( min_exprC_boolidx & subsexBoolidx[ refBoolIdx ] ), " / ", 
             sum( subsexBoolidx[ refBoolIdx ] ), 
             " ", aSex, " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
             sep="" )
    }
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 164 / 165 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 702 / 705 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 77 / 77 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 174 / 177 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 178 / 180 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1295 / 1304 female cells pass filter (min_counts_per_chromosome >= 10)
    # 
    # [1] "male" "d0"  
    # d0: 113 / 114 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 116 / 119 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 287 / 293 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 67 / 70 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 0 / 0 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 583 / 596 male cells pass filter (min_counts_per_chromosome >= 10)
}

# --------------------------------------------------

############
# For each chromosome 
for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
    # for (theChrom in theChroms) {
    # theChrom =  "chrX"
    print(theChrom)
    
    # --------------------------------------------------
    # 20200709
    # Allow for exclusion of Xist - but don't exclude it.
    
    if (rmXist) {
        
        # Exclude Xist
        # From sciRNA-seq CDS:
        # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
        # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
        # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
        
        Xist_window = 5000
        Xist_region_start = 103460375 - Xist_window
        Xist_region_end = 103483217 + Xist_window
        # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
        # [1] 8
        
        chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
    } else {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
    }
    cdsC.OG <- cds[chromgenesNames,]
    print(dim(cdsC.OG))
    
    # --------------------------------------------------

    # F121 cells only
    for (aSex in c("female")) {
        # aSex = "female"
        
        if(theChrom =="chrX" && aSex == "male") {
            next
        }
        
        subsexBoolidx = pData(cdsC.OG)$sex == aSex
        
        # filter across time points so that vectors are the same length
        for (aTP in c("all")) { 
            # for (aTP in c("all", names(tpCols))) {
            # for (aTP in c(names(tpCols))) {
            # aTP = "all"
            
            print(c(theChrom, aSex, aTP))
            
            if (aTP == "all") {
                tpSubsetBoolidx = rep(TRUE, length(pData(cdsC.OG)$timepoint))
            } else {
                tpSubsetBoolidx = pData(cdsC.OG)$timepoint == aTP
            }
            
            subsetBoolidx = tpSubsetBoolidx & subsexBoolidx
            
            cdsC.TP.Sex <- cdsC.OG[, subsetBoolidx]
            print(dim(cdsC.TP.Sex))
            
            cdsC.TP.Sex <- detectGenes(cdsC.TP.Sex, min_expr = sub_min_expr)
            
            # The vector expressed_genes now holds the identifiers for genes expressed in at least n cells of the data set.
            expressed_genes <- row.names(subset(fData(cdsC.TP.Sex),
                                                num_cells_expressed >= sub_min_cells))
            print(length(expressed_genes))
            
            cdsC.TP.Sex <- cdsC.TP.Sex[expressed_genes,]
            print(dim(cdsC.TP.Sex))
            
            refBoolIdx = grepl("_ref$", rownames(pData(cdsC.TP.Sex)))
            altBoolIdx = grepl("_alt$", rownames(pData(cdsC.TP.Sex)))
            
            if( sum(altBoolIdx) != sum(refBoolIdx) ) {
                cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
                next
            }
            
            # --------------------------------------------------
            # 20200709
            # Apply thresholds for minimum counts per chromosome 
            # For each chromosome 
            
            if ( min_counts_per_chromosome > 0 && sum(refBoolIdx) > 0 ) {
                
                all.refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
                all.altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
                
                if( sum(all.altBoolIdx) != sum(all.refBoolIdx) ) {
                    cat( sum(all.altBoolIdx), " != ", sum(all.refBoolIdx), "\n", sep="")
                    next
                }
                
                cat( aTP, ": ", sum( min_exprC_boolidx & subsetBoolidx[ all.refBoolIdx ] ), " / ", 
                     sum( subsetBoolidx[ all.refBoolIdx ] ), 
                     " ", aTP, " ", aSex, " cells pass filter", 
                     paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
                     sep="" )
                
                # min_exprC_boolidx is based on all cells regardless of sex and time point,
                # but applies to both ref and alt haplocells,
                # so must select min_exprC_boolidx entries that are sex and time point specific 
                # and apply to every alternate refBoolIdx / altBoolIdx entry.
                refBoolIdx[ seq(2, length(refBoolIdx), 2) ][ ! min_exprC_boolidx[ subsetBoolidx[ all.refBoolIdx ] ] ] = FALSE
                altBoolIdx[ seq(1, length(altBoolIdx), 2) ][ ! min_exprC_boolidx[ subsetBoolidx[ all.altBoolIdx ] ] ] = FALSE
                if( sum(altBoolIdx) != sum(refBoolIdx) ) {
                    cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
                    next
                }
                
            }
            
            # --------------------------------------------------
            
            ##################################################
            # Save filtered CDS to file
            
            CDSmatrix = as.matrix(Biobase::exprs(cdsC.TP.Sex))
            
            dim(CDSmatrix)
            # 3801 4340
            
            CDSmatrix[1:2, 1:2]
            # d0_A01_B02_RT_16_alt d0_A01_B02_RT_16_ref
            # ENSMUSG00000000037.16                    0                    0
            # ENSMUSG00000000134.17                    0                    0
            
            CDSmatrix[, refBoolIdx][1:2, 1:2]
            # ENSMUSG00000000037.16                    0                    0
            # ENSMUSG00000000134.17                    0                    0
            CDSmatrix[, altBoolIdx][1:2, 1:2]
            # ENSMUSG00000000037.16                    0                    0
            # ENSMUSG00000000134.17                    0                    0
            
            # # # pData(cds)$sex # factor( sampleSex, levels=c("female", "male") )
            # # # pData(cds)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
            # # 
            # # length(pData(cds)$sex) # factor( sampleSex, levels=c("female", "male") )
            # # 
            # # length(pData(cds)$timepoint) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
            # 
            # colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cdsC.TP.Sex)$timepoint, ".", pData(cdsC.TP.Sex)$sex, sep="")
            
            # Concatenate UMIs for two alleles
            CDSmatrix.catenated = rbind(CDSmatrix[, refBoolIdx], CDSmatrix[, altBoolIdx])
            rownames(CDSmatrix.catenated) = paste( c( rep( "ref", nrow(CDSmatrix) ),
                                                      rep( "alt", nrow(CDSmatrix) ) ),
                                                   rownames(CDSmatrix.catenated), 
                                                   sep="_" )
            colnames(CDSmatrix.catenated) = sub("_ref$", "", colnames(CDSmatrix.catenated))
            
            dim(CDSmatrix.catenated)
            # [1] 7248 1304
            
            length(cdsC.TP.Sex$timepoint[refBoolIdx])
            # [1] 1304
            
            colnames(CDSmatrix.catenated) =  paste0("F121_", cdsC.TP.Sex$timepoint[refBoolIdx], ".", colnames(CDSmatrix.catenated))
            
            matrixfilename = file.path(alleleCDPmatrixDir, paste0(projectName, ".UMImatrix.filtered", 
                                                                  ".minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                                                  ".F121.", theChrom, 
                                                                  ifelse(min_counts_per_chromosome==0, "", 
                                                                         paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                                                  ifelse(rmXist, ".rmXist", ""),
                                                                  ".tsv"))
            write.table(CDSmatrix.catenated, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
            try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
            
            print( table(cdsC.TP.Sex$timepoint[refBoolIdx]) )
            # sciATAC-seq
            # # # d0   d3   d7  d11  NPC 
            # # # 273 1102  154  380  261 
            # # d0  d3  d7 d11 NPC 
            # # 165 705  77 177 180 
            # min_counts_per_chromosome > = 10 (chr1 & chrX only)
            # d0  d3  d7 d11 NPC 
            # 164 702  77 174 178 
            
        } # timepoint
        
    } # sex
    
} # chrom










##################################################
##################################################
##################################################
# 20200615
# Heatmaps and clustering of allelic accessibility
#
# 20200629
# Only XCI cells
#
# 20200709
# Apply thresholds for minimum counts per chromosome 
# And potentially exclude Xist

# install.packages("tidyHeatmap")
require("ComplexHeatmap")

# 20200709 *** CHANGE THIS ***
# # XCIthresh = 1.2239 # 10% FPR thresholds based on chr1 distributions of allelic ratio of total chromosomal expression
XCIthresh = 0.98080 # 10% FPR thresholds based on chr1 distributions of allelic ratio of total chromosomal expression - after filterint by min_counts_per_chromosome >= 10

# Applied above:
# ##################################################
# # Code to make expression code work in this context
# 
# fData(cds.OG)$ensChr <- sub(":.*$", "", fData(cds.OG)$id)
# fData(cds)$ensStart <- as.numeric( as.character( sub("-.*$", "", sub("^.*:", "", fData(cds)$id) ) ) )
# fData(cds)$ensEnd <- as.numeric( as.character(sub("^.*-", "", fData(cds)$id) ) ) 
# 
# pData(cds.OG)$sex <- pData(cds.OG)$gender
# 
# # cds = cds.OG
# 
# ##################################################

sub_min_expr = 0.1; sub_min_cells = 5

whichChroms2analze = list()
whichChroms2analze[["chr1andchrX"]] = c(1,20)
# whichChroms2analze[["all"]] = 1:20
whichChroms2analze[["chrX"]] = 20
whichChroms2analze[["chr1"]] = 1

whichChroms = "chr1andchrX"
# whichChroms = "all"
# whichChroms = "chrX"
# whichChroms = "chr1"

plotAll=FALSE
# plotAll=TRUE

# --------------------------------------------------
# 20200709
# Apply thresholds for minimum counts per chromosome 
# For each chromosome 
#
# And exclude Xist

min_counts_per_chromosome = 10
# min_counts_per_chromosome = 0

# Exclude Xist
# rmXist = TRUE
rmXist = FALSE

if (min_counts_per_chromosome > 0) {
    
    ###
    # Get cells showing minimum total chromosomal expression per chromosome
    
    min_exprC_boolidx_per_chrom = list()
    
    for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
        # for (theChrom in theChroms) {
        
        if (rmXist) {
            
            # Exclude Xist
            # From sciRNA-seq CDS:
            # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
            # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
            # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
            
            Xist_window = 5000
            Xist_region_start = 103460375 - Xist_window
            Xist_region_end = 103483217 + Xist_window
            # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
            # [1] 8
            
            chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
        } else {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
        }
        cdsC.OG <- cds[chromgenesNames,]
        print(dim(cdsC.OG))
        
        ############
        # chromosome-wide allelic counts
        
        pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
        
        refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
        altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
        
        if( sum(altBoolIdx) != sum(refBoolIdx) ) {
            cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
            next
        }
        
        refExprC = pData(cdsC.OG)$Cexpr[refBoolIdx]
        altExprC = pData(cdsC.OG)$Cexpr[altBoolIdx]
        
        totalExprC = refExprC + altExprC
        
        min_exprC_boolidx_per_chrom[[ theChrom ]] = min_exprC_boolidx = totalExprC >= min_counts_per_chromosome
        
        cat( theChrom, ": ", sum(min_exprC_boolidx_per_chrom[[ theChrom ]]), " / ", 
             length(min_exprC_boolidx_per_chrom[[ theChrom ]]),
             " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
             sep="" )
    }
    # w/o Xist:
    # Features  Samples 
    # 4301     3800 
    # chrX: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)
    
    ###
    # Now consider cells with sufficient chromosomal counts for all chromosomes under consideration
    min_exprC_boolidx = min_exprC_boolidx_per_chrom[[ 1 ]] # intialize
    for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
        # for (theChrom in theChroms) {
        print(theChrom)
        
        min_exprC_boolidx = min_exprC_boolidx & min_exprC_boolidx_per_chrom[[ theChrom ]] 
    }
    
    cat( "All chromosomes: ", sum(min_exprC_boolidx), " / ", 
         length(min_exprC_boolidx), 
         " cells pass filter", 
         paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome), ")\n",
         sep="" )
    # w/o Xist:
    # All chromosomes: 1878 / 1900 cells pass filter (min_counts_per_chromosome >= 10)
    
    ###
    # Check min_exprC_boolidx how many female cells per time point?
    
    # all( grepl("_ref$", rownames(pData(cdsC.OG))) == grepl("_ref$", rownames(pData(cds))) )
    # TRUE
    refBoolIdx = grepl("_ref$", rownames(pData(cds)))
    altBoolIdx = grepl("_alt$", rownames(pData(cds)))
    
    for (aSex in c("female", "male")) {
        # aSex = "female"
        
        subsexBoolidx = pData(cds)$sex == aSex
        # all( pData(cdsC.OG)$sex == pData(cds)$sex )
        # [1] TRUE
        
        #for (aTP in c("all", names(tpCols))) {
        for (aTP in c(names(tpCols))) {
            # aTP = "NPC"
            
            print(c(aSex, aTP))
            
            if (aTP == "all") {
                tpSubsetBoolidx = rep(TRUE, length(pData(cds)$timepoint))
            } else {
                tpSubsetBoolidx = pData(cds)$timepoint == aTP
            }
            
            subsetBoolidx = tpSubsetBoolidx & subsexBoolidx
            
            # cds.TP.Sex <- cds[, subsetBoolidx]
            # print(dim(cdsC.TP.Sex))
            
            # For cells not haplocells
            # all( subsetBoolidx[ refBoolIdx ] == subsetBoolidx[ altBoolIdx ] )
            # TRUE
            
            cat( aTP, ": ", sum( min_exprC_boolidx & subsetBoolidx[ refBoolIdx ] ), " / ", 
                 sum( subsetBoolidx[ refBoolIdx ] ), 
                 " ", aTP, " ", aSex, " cells pass filter", 
                 paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
                 sep="" )
        }
        
        cat( sum( min_exprC_boolidx & subsexBoolidx[ refBoolIdx ] ), " / ", 
             sum( subsexBoolidx[ refBoolIdx ] ), 
             " ", aSex, " cells pass filter", 
             paste0(" (min_counts_per_chromosome >= ", min_counts_per_chromosome ), ")\n",
             sep="" )
    }
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 164 / 165 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 702 / 705 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 77 / 77 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 174 / 177 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 178 / 180 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1295 / 1304 female cells pass filter (min_counts_per_chromosome >= 10)
    # 
    # [1] "male" "d0"  
    # d0: 113 / 114 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 116 / 119 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 287 / 293 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 67 / 70 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 0 / 0 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 583 / 596 male cells pass filter (min_counts_per_chromosome >= 10)
}

# --------------------------------------------------

for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
    # theChrom = "chrX"
    # theChrom = "chr1"
    print(theChrom)
    
    thisOutputDir = file.path( superMainOutputDir, 
                               paste0( "heatmaps",
                                       ifelse(min_counts_per_chromosome==0, "", 
                                              paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                       ifelse(rmXist, ".rmXist", "") ),
                               theChrom )
    if ( ! file.exists( thisOutputDir ) ) {
        if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = superMainOutputDir
    }
    
    chromgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == theChrom) ]
    
    # --------------------------------------------------
    # 20200709
    # Allow for exclusion of Xist - but don't exclude it.
    
    if (rmXist) {
        
        # Exclude Xist
        # From sciRNA-seq CDS:
        # fData(cds)[ fData(cds)$gene_short_name == "Xist" , ]
        # gene_short_name                   id ensChr  ensStart    ensEnd ensStrand
        # # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3   chrX 103460375 103483217         -
        
        Xist_window = 5000
        Xist_region_start = 103460375 - Xist_window
        Xist_region_end = 103483217 + Xist_window
        # sum(fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end)
        # [1] 8
        
        chromgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr == theChrom & !( fData(cds)$ensChr == theChrom & fData(cds)$ensStart >= Xist_region_start & fData(cds)$ensEnd <= Xist_region_end ) ]
    } else {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
    }
    cdsC.OG <- cds[chromgenesNames,]
    print(dim(cdsC.OG))
    
    # --------------------------------------------------
    
    # --------------------------------------------------
    # 20200709
    # Apply thresholds for minimum counts per chromosome 
    # For each chromosome 
    
    if ( min_counts_per_chromosome > 0 ) {
        cdsC.OG <- cdsC.OG[,rep(min_exprC_boolidx, each=2)]
        print(dim(cdsC.OG))
    }
    
    # --------------------------------------------------
    
    # 20200617 Not sure the is necessary.
    # The vector expressed_genes now holds the identifiers for genes expressed in at least n cells of the data set.
    cdsC.OG <- detectGenes(cdsC.OG, min_expr = sub_min_expr)
    expressed_genes <- row.names(subset(fData(cdsC.OG),
                                        num_cells_expressed >= sub_min_cells))
    print(length(expressed_genes))
    
    cdsC.OG <- cdsC.OG[expressed_genes,]
    print( "After filtering for expressed genes: " )
    print( dim(cdsC.OG) ) 
    
    # Gene coverage
    CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
    biallelicGeneCoverage = rowSums(CDSmatrix)
    
    zeroCovBoolIdx = biallelicGeneCoverage == 0 
    print( c( sum(zeroCovBoolIdx),
              round(sum(zeroCovBoolIdx) / length(zeroCovBoolIdx) * 100, 2) ) )
    
    cdsC.OG = cdsC.OG[!zeroCovBoolIdx,]
    print(dim(cdsC.OG))
    
    # Allelic totals per Gene
    refBoolIdx = grepl("_ref$", rownames(pData(cdsC.OG)))
    altBoolIdx = grepl("_alt$", rownames(pData(cdsC.OG)))
    if( sum(altBoolIdx) != sum(refBoolIdx) ) {
        cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
        next
    }
    
    # sex
    femaleBoolIdx = pData(cdsC.OG)$sex == "female"
    maleBoolIdx = pData(cdsC.OG)$sex == "male"
    
    # rownames(pData(cdsC.OG))[1:3]
    # # [1] "d0_A01_B02_RT_17_alt"
    # # [2] "d0_A01_B02_RT_17_ref"
    # # [3] "d0_A01_B02_RT_28_alt"
    # 
    # rownames(pData(cdsC.OG))[refBoolIdx][1:3]
    # # [1] "d0_A01_B02_RT_17_ref" "d0_A01_B02_RT_28_ref" "d0_A01_B03_RT_27_ref"
    
    
    ################
    # Get cell XCI status
    
    # dim(XCIcellsAll[[whichChroms]])[1]
    # 1344 cells 
    
    # XCIcellsAll[[whichChroms]][1:3,]
    # # refExpr altExpr   log2Expr auto10FPRthreshold XCI
    # # d0_A01_B02_RT_17       5       9 -0.8353693           1.223988   0
    # # d0_A01_B03_RT_6       14       6  1.2088140           1.223988   1
    # # d0_A02_B02_RT_17       0       1 -3.4594316           1.223988   1
    
    theseCellIDs = sub("_ref", "", rownames(pData(cdsC.OG))[refBoolIdx])
    length(theseCellIDs)
    # [1] 1878
    
    cellsWithXCIstatus = match( theseCellIDs, rownames(XCIcellsAll[[whichChroms]]) )
    # all( rownames(theseCellIDs)[!is.na(cellsWithXCIstatus)] == 
    #          rownames(XCIcellsAll[[whichChroms]])[na.omit(cellsWithXCIstatus)] )
    
    theseCellIDsXCIstatus = rep(NaN, length(theseCellIDs))
    theseCellIDsXCIstatus[ !is.na(cellsWithXCIstatus) ] = XCIcellsAll[[whichChroms]][,"log2Expr"][ na.omit(cellsWithXCIstatus) ]
    pData(cdsC.OG)$XCIstatus[refBoolIdx] = as.numeric( as.character( theseCellIDsXCIstatus ))
    pData(cdsC.OG)$XCIstatus[altBoolIdx] = as.numeric( as.character( theseCellIDsXCIstatus ))
    
    # # Sort cells by XCI status for cells the have it
    # cellXCIstatusOrder = order( as.numeric( as.character( pData(cdsC.OG)$XCIstatus[refBoolIdx] ) ) )
    
    XCIboolIdx = abs(pData(cdsC.OG)$XCIstatus) > XCIthresh
    
    # Sort genes by location
    genePosOrder = order( as.numeric( as.character( fData(cdsC.OG)$ensStart ) ) )
    cdsC.OG = cdsC.OG[ genePosOrder, ] 
    
    # color scales
    logUMIcols = circlize::colorRamp2(breaks = range(log10(pData(cdsC.OG)$n.umi[refBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx] + 0.1), na.rm=TRUE),
                                      colors = c("white", "green4"),
                                      transparency = .2)
    
    XCIcolors = circlize::colorRamp2(breaks = c( min(pData(cdsC.OG)$XCIstatus[refBoolIdx], na.rm=TRUE), 0 , max(pData(cdsC.OG)$XCIstatus[refBoolIdx], na.rm=TRUE) ),
                                     colors = c("blue", "white", "red"),
                                     transparency = .2)
    
    
    
    ##########
    # Allelic accessibility
    
    CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
    # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
    rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
    
    CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
    # # 20200617 Use mad and IQR instead of var but always 0
    # CDS.var = apply(CDSmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
    # CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
    topQuantThresh = 0.75
    topQuantValue = quantile(CDS.var, p=topQuantThresh, na.rm=TRUE)
    topQuantGeneBoolIdx = CDS.var >= topQuantValue
    CDSmatrix.topQuantVar = CDSmatrix[ topQuantGeneBoolIdx,  ]
    # dim(CDSmatrix.topQuantVar)
    
    # Print top most variable genes 
    topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
    toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
    toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.topQuantVar)[topvarGenesBoolIdx]
    rownames(CDSmatrix.topQuantVar) = toprownames
    
    # Stack allelic genes one on top of the other
    CDSmatrix.diploid.topQuantVar = rbind(CDSmatrix.topQuantVar[,refBoolIdx],
                                          CDSmatrix.topQuantVar[,altBoolIdx])
    # dim(CDSmatrix.diploid.topQuantVar)
    
    # sort cells by sex then time point
    sexThenTpOrder = order( pData(cdsC.OG)$sex[refBoolIdx], pData(cdsC.OG)$timepoint[refBoolIdx] )
    # pData(cdsC.OG)[refBoolIdx,c("sex", "timepoint")][sexThenTpOrder,]
    CDSmatrix.diploid.topQuantVar = CDSmatrix.diploid.topQuantVar[,sexThenTpOrder]
    # dim(CDSmatrix.diploid.topQuantVar)
    
    col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx][sexThenTpOrder], 
                               sex = pData(cdsC.OG)$sex[refBoolIdx][sexThenTpOrder],
                               log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx] + 0.1)[sexThenTpOrder],
                               XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx][sexThenTpOrder], 
                               col = list( sex = sexCols, timepoint = tpCols,
                                           log10totalUMIs = logUMIcols,
                                           XCIstatus = XCIcolors ))
    
    # Can't annotate rows in this version 
    # row_ha = HeatmapAnnotation(allele = c(rep("129", sum(refBoolIdx)),
    #                                       rep("Cast", sum(altBoolIdx)) ),
    #                            col = list( allele = alleleCols))
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicAccessibility",
                                                 ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                 ".top", round((1-topQuantThresh)*100), "pctVar",
                                                 ".", theChrom, ".pdf", sep=""))
    pdf(PDFfilename, width=10, height=10)
    
    Heatmap(log10(CDSmatrix.diploid.topQuantVar+0.1), 
            cluster_rows = FALSE, 
            cluster_columns = FALSE, 
            top_annotation = col_ha, 
            show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
            show_column_names=FALSE,
            split = c( rep("129", dim(CDSmatrix.topQuantVar+0.1)[1] ),
                       rep("Cast", dim(CDSmatrix.topQuantVar+0.1)[1]) ),
            column_title = paste0( "Cells (n=", dim(CDSmatrix.topQuantVar+0.1)[2] / 2, ")" ), 
            row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " genes (n=", dim(CDSmatrix.topQuantVar+0.1)[1], " per allele)" ),
            name = "log10 UMIs", 
            use_raster = TRUE, raster_device = "png")
    
    dev.off()
    
    
    
    if(plotAll) {
        
        ##########
        # Allelic gene accessibility - female cells only 
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        CDSmatrix = CDSmatrix[,femaleBoolIdx]
        
        CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        # CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.topQuantVar = CDSmatrix[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.topQuantVar)
        
        # Print top most variable genes 
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.topQuantVar) = toprownames
        
        # Stack allelic genes one on top of the other
        CDSmatrix.diploid.topQuantVar = rbind(CDSmatrix.topQuantVar[ ,refBoolIdx[femaleBoolIdx] ],
                                              CDSmatrix.topQuantVar[ ,altBoolIdx[femaleBoolIdx] ])
        # dim(CDSmatrix.diploid.topQuantVar)
        
        # sort cells by time point
        tpOrder = order( pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx] )
        # pData(cdsC.OG)[refBoolIdx & femaleBoolIdx, c("timepoint")][tpOrder]
        CDSmatrix.diploid.topQuantVar = CDSmatrix.diploid.topQuantVar[,tpOrder]
        # dim(CDSmatrix.diploid.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][tpOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[tpOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][tpOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        # Can't annotate rows in this version 
        # row_ha = HeatmapAnnotation(allele = c(rep("129", sum(refBoolIdx)),
        #                                       rep("Cast", sum(altBoolIdx)) ),
        #                            col = list( allele = alleleCols))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicAccessibility",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, ".femaleCells.pdf", sep=""))
        pdf(PDFfilename, width=10, height=10)
        
        Heatmap(log10(CDSmatrix.diploid.topQuantVar+0.1), 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                split = c( rep("129", dim(CDSmatrix.topQuantVar+0.1)[1] ),
                           rep("Cast", dim(CDSmatrix.topQuantVar+0.1)[1]) ),
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.topQuantVar+0.1)[2] / 2, ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " genes (n=", dim(CDSmatrix.topQuantVar+0.1)[1], " per allele)" ),
                name = "log10 UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # Allelic gene accessibility - female cells only 
        # with XCI status > 10% FPR threshold base on chr1 XCI distribution
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        CDSmatrix = CDSmatrix[,femaleBoolIdx]
        # dim(CDSmatrix)
        
        CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        # CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.topQuantVar = CDSmatrix[ topQuantGeneBoolIdx,  ]
        dim(CDSmatrix.topQuantVar)
        
        # Print top most variable genes 
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.topQuantVar) = toprownames
        
        # Stack allelic genes one on top of the other
        CDSmatrix.diploid.topQuantVar = rbind(CDSmatrix.topQuantVar[ ,refBoolIdx[XCIboolIdx & femaleBoolIdx] ],
                                              CDSmatrix.topQuantVar[ ,altBoolIdx[XCIboolIdx & femaleBoolIdx] ])
        # dim(CDSmatrix.diploid.topQuantVar)
        
        # sort cells by time point
        tpOrder = order( pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx & XCIboolIdx] )
        # pData(cdsC.OG)[refBoolIdx & femaleBoolIdx, c("timepoint")][tpOrder]
        CDSmatrix.diploid.topQuantVar = CDSmatrix.diploid.topQuantVar[,tpOrder]
        # dim(CDSmatrix.diploid.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx & XCIboolIdx][tpOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx & XCIboolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[tpOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx][tpOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        # Can't annotate rows in this version 
        # row_ha = HeatmapAnnotation(allele = c(rep("129", sum(refBoolIdx)),
        #                                       rep("Cast", sum(altBoolIdx)) ),
        #                            col = list( allele = alleleCols))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicAccessibility",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, ".femaleXCIcells.pdf", sep=""))
        pdf(PDFfilename, width=10, height=10)
        
        Heatmap(log10(CDSmatrix.diploid.topQuantVar+0.1), 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                split = c( rep("129", dim(CDSmatrix.topQuantVar+0.1)[1] ),
                           rep("Cast", dim(CDSmatrix.topQuantVar+0.1)[1]) ),
                column_title = paste0( "Female XCI cells (n=", dim(CDSmatrix.topQuantVar+0.1)[2] / 2, ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " genes (n=", dim(CDSmatrix.topQuantVar+0.1)[1], " per allele)" ),
                name = "log10 UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[,refBoolIdx] + 0.1) /
                                          (CDSmatrix[,altBoolIdx] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by sex then time point
        sexThenTpOrder = order( pData(cdsC.OG)$sex[refBoolIdx], pData(cdsC.OG)$timepoint[refBoolIdx] )
        # pData(cdsC.OG)[refBoolIdx,c("sex", "timepoint")][sexThenTpOrder,]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,sexThenTpOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx][sexThenTpOrder], 
                                   sex = pData(cdsC.OG)$sex[refBoolIdx][sexThenTpOrder],
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx] + 0.1)[sexThenTpOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx][sexThenTpOrder], 
                                   col = list( sex = sexCols, timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, ".pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by time point
        tpOrder = order( pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx] )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,tpOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][tpOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[tpOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][tpOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleCells.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off() 
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by ABSOLUTE XCI log2 ratio
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx& femaleBoolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( abs( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx] ), decreasing=FALSE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleCells.absXCIorder.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleCells.XCIorder.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio
        # with XCI status > 10% FPR threshold base on chr1 XCI distribution
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx & XCIboolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx & XCIboolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleXCIcells.XCIorder.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female XCI cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio - cluster genes
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleCells.XCIorder.clusteredGenes.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = TRUE, 
                clustering_distance_rows = "spearman",
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio - cluster genes
        # with XCI status > 10% FPR threshold base on chr1 XCI distribution
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        topQuantThresh = 0.75
        topQuantValue = quantile(CDS.var, p=topQuantThresh)
        topQuantGeneBoolIdx = CDS.var >= topQuantValue
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio[ topQuantGeneBoolIdx,  ]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        # Print top most variable genes    
        topvarGenesBoolIdx = order( CDS.var[ topQuantGeneBoolIdx ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(topQuantGeneBoolIdx) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.topQuantVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.topQuantVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.topQuantVar = CDSmatrix.allelicRatio.topQuantVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.topQuantVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx & XCIboolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".top", round((1-topQuantThresh)*100), "pctVar",
                                                     ".", theChrom, "femaleXCIcells.XCIorder.clusteredGenes.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.topQuantVar, 
                cluster_rows = TRUE, 
                clustering_distance_rows = "spearman",
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female XCI cells (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[2], ")" ), 
                row_title = paste0( "top ", round((1-topQuantThresh)*100), "% most variable ", theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.topQuantVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio - cluster genes - all genes
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        
        # Filter zero variance genes for clustering
        # which(CDS.var == 0) 
        # # Rhox5    Tsx Tceal8 
        # # 69    223    270 
        CDSmatrix.allelicRatio.allVar = CDSmatrix.allelicRatio[ CDS.var != 0, ]
        # dim(CDSmatrix.allelicRatio.allVar)
        
        # Print top genes    
        topvarGenesBoolIdx = order( CDS.var[ CDS.var != 0 ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(CDS.var != 0) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.allVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.allVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.allVar = CDSmatrix.allelicRatio.allVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.allVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".allGenes",
                                                     ".", theChrom, "femaleCells.XCIorder.clusteredGenes.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.allVar, 
                cluster_rows = TRUE, 
                clustering_distance_rows = "spearman",
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female cells (n=", dim(CDSmatrix.allelicRatio.allVar+0.1)[2], ")" ), 
                row_title = paste0( theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.allVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
        
        
        ##########
        # log2 ratio allelic gene accessibility - female cells - ranked by XCI log2 ratio - cluster genes - all genes
        # with XCI status > 10% FPR threshold base on chr1 XCI distribution
        
        CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
        # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
        rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
        
        # log2 allelic ratio for all genes
        CDSmatrix.allelicRatio= log2( (CDSmatrix[ ,refBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) /
                                          (CDSmatrix[ ,altBoolIdx & femaleBoolIdx & XCIboolIdx ] + 0.1) )
        # dim(CDSmatrix.allelicRatio)
        
        CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, var, na.rm=TRUE)
        # # 20200617 Use mad and IQR instead of var but always 0
        # CDS.var = apply(CDSmatrix.allelicRatio, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
        
        # Filter zero variance genes for clustering
        # which(CDS.var == 0) 
        # # Rhox5    Tsx Tceal8 
        # # 69    223    270 
        CDSmatrix.allelicRatio.allVar = CDSmatrix.allelicRatio[ CDS.var != 0, ]
        # dim(CDSmatrix.allelicRatio.allVar)
        
        # Print top genes    
        topvarGenesBoolIdx = order( CDS.var[ CDS.var != 0 ], decreasing=TRUE ) <= 50
        toprownames = rep("", sum(CDS.var != 0) ) 
        toprownames[topvarGenesBoolIdx] = rownames(CDSmatrix.allelicRatio.allVar)[topvarGenesBoolIdx]
        rownames(CDSmatrix.allelicRatio.allVar) = toprownames
        
        # sort cells by XCI status
        xciOrder = order( pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx], decreasing=TRUE )
        # pData(cdsC.OG)[refBoolIdx,c(timepoint")][tpOrder]
        CDSmatrix.allelicRatio.allVar = CDSmatrix.allelicRatio.allVar[,xciOrder]
        # dim(CDSmatrix.allelicRatio.allVar)
        
        col_ha = HeatmapAnnotation(timepoint = pData(cdsC.OG)$timepoint[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx & XCIboolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx & XCIboolIdx] + 0.1)[xciOrder],
                                   XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx][xciOrder], 
                                   col = list( timepoint = tpCols,
                                               log10totalUMIs = logUMIcols,
                                               XCIstatus = XCIcolors ))
        
        PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.log2allelicRatio",
                                                     ".minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                     ".allGenes",
                                                     ".", theChrom, "femaleXCIcells.XCIorder.clusteredGenes.pdf", sep=""))
        pdf(PDFfilename, width=10, height=6)
        
        Heatmap(CDSmatrix.allelicRatio.allVar, 
                cluster_rows = TRUE, 
                clustering_distance_rows = "spearman",
                cluster_columns = FALSE, 
                top_annotation = col_ha, 
                show_row_names=TRUE, row_names_gp = gpar(fontsize = 5),
                show_column_names=FALSE,
                column_title = paste0( "Female XCI cells (n=", dim(CDSmatrix.allelicRatio.allVar+0.1)[2], ")" ), 
                row_title = paste0( theChrom, " log2 ratio genes (n=", dim(CDSmatrix.allelicRatio.allVar+0.1)[1], ")" ),
                name = "log2 ratio allelic UMIs", 
                use_raster = TRUE, raster_device = "png")
        
        dev.off()
        
    } # plotAll
    
} # chrom










##################################################
##################################################
##################################################
# Partion genes into chromosome-specific groups

GgenesNames <- row.names(fData(cds))[ ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
cdsG.OG = cds[GgenesNames,]
dim(cdsG.OG)

# 20191112 Exclude Y genes from autosomes
AgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr != "chrX" & fData(cds)$ensChr != "chrY" & ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
cdsA.OG = cds[AgenesNames,]
dim(cdsA.OG)

XgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ]
cdsX.OG <- cds[XgenesNames,]
dim(cdsX.OG)

# 20200421 
XgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ]
femaleCells <- row.names(pData(cds))[ which(pData(cds)$sex == "female") ]
cdsX.female.OG <- cds[XgenesNames, femaleCells]
dim(cdsX.female.OG)

# 20200421 
GgenesNames <- row.names(fData(cds))[ ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
femaleCells <- row.names(pData(cds))[ which(pData(cds)$sex == "female") ]
cdsG.female.OG <- cds[GgenesNames, femaleCells]
dim(cdsG.female.OG)

# 20200421 
GgenesNames <- row.names(fData(cds))[ ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
maleCells <- row.names(pData(cds))[ which(pData(cds)$sex == "male") ]
cdsG.male.OG <- cds[GgenesNames, maleCells]
dim(cdsG.male.OG)

# 20200421 
AgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr != "chrX" & fData(cds)$ensChr != "chrY" & ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
femaleCells <- row.names(pData(cds))[ which(pData(cds)$sex == "female") ]
cdsA.female.OG <- cds[AgenesNames, femaleCells]
dim(cdsA.female.OG)

# 20200421 
AgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr != "chrX" & fData(cds)$ensChr != "chrY" & ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
maleCells <- row.names(pData(cds))[ which(pData(cds)$sex == "male") ]
cdsA.male.OG <- cds[AgenesNames, maleCells]
dim(cdsA.male.OG)


#pData(cds)$totalExprV2 <- Matrix::colSums(Biobase::exprs(cds))  
pData(cdsG.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsG.OG))  
pData(cdsA.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsA.OG))  
pData(cdsX.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsX.OG))  
pData(cdsX.female.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsX.female.OG))  
pData(cdsG.female.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsG.female.OG))  
pData(cdsG.male.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsG.male.OG))  
pData(cdsA.female.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsA.female.OG))  
pData(cdsA.male.OG)$Texpr <- Matrix::colSums(Biobase::exprs(cdsA.male.OG))  



##################################################
# List of CDSs per region

CDSs = list()
CDSs[["all"]] <- cds
CDSs[["genome"]] <- cdsG.OG
CDSs[["autosomes"]] <- cdsA.OG
CDSs[["chrX"]] <- cdsX.OG
CDSs[["chrX.female"]] <- cdsX.female.OG
CDSs[["genome.female"]] <- cdsG.female.OG
CDSs[["genome.male"]] <- cdsG.male.OG
CDSs[["autosomes.female"]] <- cdsA.female.OG
CDSs[["autosomes.male"]] <- cdsA.male.OG

names(CDSs)

rm(cdsG.OG, cdsA.OG, cdsX.OG, cdsX.female.OG, cdsG.female.OG,  cdsG.male.OG,  cdsA.female.OG,  cdsA.male.OG, 
   # cdsG, cdsA, cdsX, cdsX.female,
   GgenesNames, AgenesNames, XgenesNames)
gc()



##################################################
# Remove singlets again b/c of filtering.
# i.e. cells with reads from only one parental allele
#
# What about for chrX? 
# No, data not filtered for chrX b/c haplocells could be Xa alleles.

for (whichChr in names(CDSs) ){
    # whichChr="autosomes"
    print(whichChr)
    
    if (whichChr %in% c("chrX") ) {
        next 
    }
    
    # Remove singlets
    i = 1
    singletList = c()
    while( i < nrow(pData(CDSs[[whichChr]])) ) {
        if(sub("_...$", "", pData(CDSs[[whichChr]])$cell[i]) != sub("_...$", "", pData(CDSs[[whichChr]])$cell[i+1])) {
            singletList = c(singletList, i)
            i = i + 1
        } else {
            i = i + 2
        }
    }
    print(length(singletList))
    pData(CDSs[[whichChr]])$cell[singletList]
    pData(CDSs[[whichChr]])$cell[singletList-1]
    pData(CDSs[[whichChr]])$cell[singletList-2]
    pData(CDSs[[whichChr]])$cell[singletList+1]
    pData(CDSs[[whichChr]])$cell[singletList+2]
    
    if (length(singletList) > 0 ){
        CDSs[[whichChr]] <- CDSs[[whichChr]][,-singletList]
    }
}  



##################################################
# 2019723
# Save images for later in current folder

gc()
sort( sapply(mget(ls()),object.size) )
format(object.size(x=lapply(ls(), get)), unit='auto')

save.image(file=file.path(superMainOutputDir,".RData")) # To .RData
# # load(file=file.path(superMainOutputDir,".RData"))
