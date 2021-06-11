##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20190505 non-allelic sci-ATAC-seq analysis using monocle 2 
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

projectName = "non-allelic_sci-ATAC-seq""

inputDir = "../data/non-allelic_sci-ATAC-seq_count_matrices/"

mainOutputDir = "../results/non-allelic_sci-RNA-seq_monocle2/"
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










##################################################
##################################################
##################################################
# Save unfiltered data w/ meta data
# 20200516 Used to be saved in ~cisTopic.R script

dim(count.matrix.OG)
# 157060   5586

# Rename cells
cell_list.OG <- colnames(count.matrix.OG)

# Setup cell metadata
bcIdx = match( colnames(count.matrix.OG), barcodeInfo$barcode )
all( colnames(count.matrix.OG) == barcodeInfo$barcode[bcIdx] )
# TRUE

# time-points
timepoints.OG = as.matrix(barcodeInfo$timepoint[bcIdx])
rownames(timepoints.OG) <- colnames(count.matrix.OG)
colnames(timepoints.OG) <- "timepoints"
table(timepoints.OG)
# timepoints.OG
# d0  d11   d3   d7  NPC 
# 632 1154 2096 1297  407 

# sex
sex.OG = as.matrix(barcodeInfo$sex[bcIdx])
rownames(sex.OG) <- colnames(count.matrix.OG)
colnames(sex.OG) <- "sex"
table(sex.OG)
# female   male 
# 3613   1973 

table(timepoints.OG, sex.OG)
# timepoints.OG female male
# d0     382  250
# d3    1673  423
# d7     441  856
# d11    710  444
# NPC    407    0

# cell_list.OG = sub("\\..*", "", cell_list.OG)
cell_list.OG.annotated = paste(cell_list.OG, timepoints.OG, sex.OG, sep=".")
names(cell_list.OG.annotated) = cell_list.OG

count.matrix.OG.annotated = count.matrix.OG
colnames(count.matrix.OG.annotated) = cell_list.OG.annotated

# 20200516 Used to be saved in ~cisTopic.R script
matrixfilename = file.path(mainOutputDir, paste(projectName, ".countsMatrix.unfiltered.tsv", sep=""))
write.table(count.matrix.OG.annotated, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))










##################################################
##################################################
##################################################
# 20200727
# Check total umis per timepoint and sex
# For comparison against 

for (sex in c("female", "male")) {
    for (timepoint in names(tpCols)) {
        subsetBoolidx = sex.OG == sex & timepoints.OG == timepoint 
        cat(sex, "\t", timepoint, "\t", 
            format( trim=TRUE, scientific=FALSE, big.mark=",",
                    sum( count.matrix.OG.annotated[, subsetBoolidx] ) ), 
            "\n", sep="")
    }
}
cat( "total\t", 
     format( trim=TRUE, scientific=FALSE, big.mark=",",
             sum( count.matrix.OG.annotated ) ), 
     "\n", sep="")

# female  d0      1,695,880
# female  d3      11,335,101
# female  d7      1,101,693
# female  d11     1,947,718
# female  NPC     1,964,760
# male    d0      1,325,327
# male    d3      2,015,546
# male    d7      3,119,345
# male    d11     719,196

# total 25,224,566










##################################################
##################################################
##################################################
# 20200516 
# Filter for low counts

totalCoverage = colSums(count.matrix.OG, na.rm=TRUE)

length(totalCoverage)
# 5586

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5      69     457    4516    5296  124568 


# ***************************************************
# 20200516 
#
# Not sure why this there are cells with counts < 100.
# These should have been filtered out during processing by
# /Users/ajs/Documents/Gian's Stuff/NobleLab/proj/2019_sciATAC_analysis/results/gbonora/ATAC-Pipeline-Dev/mergeall.GB.py
# And the outputs suggests that they were:
# 
# $ more processReads/readProcessing/sciATAC_mouseDiff/qc_info/QC_stats.txt
#
# Mergeall QC info for ['/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDiff_Nmasked/readPro
#                       cessing/sciATAC_mouseDiff.split.q10.sorted.bam']Mergeall QC info for ['/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data
#                                                                                             /data_20191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff.split.q10.sort.bam']
# Total reads after deduplication: 93253924
# 
# Total reads after removing blacklist: 0 /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouseDif
# f_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.clean.bed
# 
# Total reads after removing bad barcodes: Mergeall QC info for ['/net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_2
#                                                                0191105_sciATAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff.split.q10.sort.bam']
# Total reads after removing blacklist: 91582452 /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_m
# ouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.clean.bed
# 
# Total reads after removing bad barcodes: 89851658 /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATA
# C_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.clean.bed
# 
# Read per cell cutoff: 100
# 
# Total reads after removing low read cells: 89496592 /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciA
# TAC_mouseDiff_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.for_macs.bed
# 
# Total cells with more than 100 reads: 5586 /net/noble/vol8/gbonora/sciOmics/2019_sciATAC_analysis/data/data_20191105_sciATAC_mouse
# Diff_Nmasked/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.high_read_cells.txt
#
# ***************************************************
#
# And it seems this is due to the input file produced by mergeall.GB.py, not my processing thereof in R:
# $ cat ~/proj/2019_sciATAC_analysis/results/gbonora/20191105_sciATAC_mouseDiff_Nmasked/processReads/readProcessing/sciATAC_mouseDiff/sciATAC_mouseDiff.counts.txt | \
# awk '{cells[$2]+=$3} END{for(cellID in cells) {print cells[cellID]}}' | sort -nr | tail -n 10
# 13
# 13
# 12
# 12
# 12
# 10
# 8
# 7
# 5
#
# ***************************************************


# Find biased Genes adjusted for bias
totCovDensity <- density(log10(totalCoverage))
maxTCdensity = max(totCovDensity$y)

totCovLHDensity <- density(log10(totalCoverage[totalCoverage < 500]))
lhMaxTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(500) & totCovLHDensity$x>0][which.max(totCovLHDensity$y[totCovLHDensity$x<log10(500) & totCovLHDensity$x>0])]
lhMinTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(500) & totCovLHDensity$x>lhMaxTotCovDensity][which.min(totCovLHDensity$y[totCovLHDensity$x<log10(500) & totCovLHDensity$x>lhMaxTotCovDensity])]
print(c(lhMinTotCovDensity, 10^lhMinTotCovDensity))

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.coverageDistribution.unfiltered.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
plot(totCovDensity, xlab="log10(total UMIs)")
abline(v=lhMinTotCovDensity, lty="dashed", col="red")
text(lhMinTotCovDensity, 0, round(10^lhMinTotCovDensity, 2), pos=4)
dev.off()

# altLowCellsBoolIdx = totalCoverage > lhMinTotCovDensity
altLowCellsBoolIdx = totalCoverage > 500

count.matrix.OG = count.matrix.OG[,altLowCellsBoolIdx]

dim(count.matrix.OG)
# # 7021 3107 : chr12 - unfiltered
# # 156,369 3109 
# 
# 157,060   5,586
#
# 157060   2750

totalCoverage = colSums(as.matrix(count.matrix.OG, na.rm=TRUE))

length(totalCoverage)
# 2750

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 501    1851    5508    9054   12272  124568 


# cell IDs
cell_list.OG <- colnames(count.matrix.OG)










##################################################
##################################################
##################################################
# Save filtered data w/ meta data

# Setup cell metadata
bcIdx = match( colnames(count.matrix.OG), barcodeInfo$barcode )
all( colnames(count.matrix.OG) == barcodeInfo$barcode[bcIdx] )
# TRUE

# time-points
timepoints.OG = as.matrix(barcodeInfo$timepoint[bcIdx])
rownames(timepoints.OG) <- colnames(count.matrix.OG)
colnames(timepoints.OG) <- "timepoints"
table(timepoints.OG)
# d0  d11   d3   d7  NPC 
# 367  465 1116  568  234 

# sex
sex.OG = as.matrix(barcodeInfo$sex[bcIdx])
rownames(sex.OG) <- colnames(count.matrix.OG)
colnames(sex.OG) <- "sex"
table(sex.OG)
# female   male 
# 1828    922 

table(timepoints.OG, sex.OG)
# timepoints.OG female male
# d0     227  140
# d3     961  155
# d7     108  460
# d11    298  167
# NPC    234    0

# cell_list.OG = sub("\\..*", "", cell_list.OG)
cell_list.OG.annotated = paste(cell_list.OG, timepoints.OG, sex.OG, sep=".")
names(cell_list.OG.annotated) = cell_list.OG

count.matrix.OG.annotated = count.matrix.OG
colnames(count.matrix.OG.annotated) = cell_list.OG.annotated

matrixfilename = file.path("../data", paste(projectName, ".countsMatrix.filtered.500UMIs.tsv", sep=""))
write.table(count.matrix.OG.annotated, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))










##################################################
##################################################
##################################################
# Generate CDS data

sapply(count.matrix.OG, class)
# as.matrix(count.matrix.OG[1:3,1:3])

# require(Matrix)
count.matrix.OG.sparse = Matrix::Matrix(as.matrix(count.matrix.OG), sparse=TRUE)

# count.matrix.OG = Matrix::Matrix(count.matrix.OG.sparse, sparse=FALSE)
# dim(count.matrix.OG)
# # [1] 157060   2836

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

all( pData(cds.OG)$n.umi == pData(cds.OG)$totalCoverage )
# TRUE










##################################################
##################################################
##################################################
# Explore data

dim(pData(cds.OG))
# [1] 2750    4

dim(fData(cds.OG))
# [1] 157060      2

dim(Biobase::exprs(cds.OG))
# 1] 157060   2750

# as.matrix(Biobase::exprs(cds.OG))










##################################################
##################################################
##################################################
# Add metadata

bcIdx = match( colnames(count.matrix.OG), barcodeInfo$barcode )
all( colnames(count.matrix.OG) == barcodeInfo$barcode[bcIdx] )
# TRUE

# time-points
timepoints = as.matrix(barcodeInfo$timepoint[bcIdx])
rownames(timepoints) <- colnames(count.matrix.OG)
colnames(timepoints) <- "timepoints"
table(timepoints)
# d0  d11   d3   d7  NPC 
# 367  465 1116  568  234 

# sex
sex = as.matrix(barcodeInfo$sex[bcIdx])
rownames(sex) <- colnames(count.matrix.OG)
colnames(sex) <- "sex"
table(sex)
# female   male 
# 1828    922 

table(timepoints, sex)
# timepoints.OG female male
# d0     227  140
# d3     961  155
# d7     108  460
# d11    298  167
# NPC    234    0

all( colnames(cds.OG) == rownames(timepoints) )
# TRUE
all( colnames(cds.OG) == rownames(sex) )
# TRUE

pData(cds.OG)$sex <- factor( sex, levels=c("female", "male") )
# pData(cds.OG)$sex <- factor( pData(cds.OG)$sex , levels=c("male", "female", ) )
pData(cds.OG)$timepoint <- factor( timepoints, levels=c("d0", "d3", "d7", "d11", "NPC") )

# save.image() # To .RData
# load(file = ".RData")










##################################################
##################################################
##################################################
# Normalize and pre-process the data

cds.OG <- estimateSizeFactors(cds.OG)
pData(cds.OG)[1:5,]
# cell
# AGCGATAGAACGAATTCGGGCGGAGGTCGTACTGAC AGCGATAGAACGAATTCGGGCGGAGGTCGTACTGAC
# AGCGATAGAACGAATTCGGGCGGAGGTCTATAGCCT AGCGATAGAACGAATTCGGGCGGAGGTCTATAGCCT
# AGCGATAGAACGAATTCGGGTCTCGCCGCAGGACGT AGCGATAGAACGAATTCGGGTCTCGCCGCAGGACGT
# AGCGATAGAACGAATTCGGGTCTCGCCGGGCTCTGA AGCGATAGAACGAATTCGGGTCTCGCCGGGCTCTGA
# AGCGATAGAACGAATTCGGTAATGATCGCCTATCCT AGCGATAGAACGAATTCGGTAATGATCGCCTATCCT
# sample Size_Factor n.umi sex timepoint
# AGCGATAGAACGAATTCGGGCGGAGGTCGTACTGAC sciATAC   2.4012234  8522 female       NPC
# AGCGATAGAACGAATTCGGGCGGAGGTCTATAGCCT sciATAC   3.1935539 11334 female       NPC
# AGCGATAGAACGAATTCGGGTCTCGCCGCAGGACGT sciATAC   1.9351798  6868 female       NPC
# AGCGATAGAACGAATTCGGGTCTCGCCGGGCTCTGA sciATAC   0.6689163  2374 female       NPC
# AGCGATAGAACGAATTCGGTAATGATCGCCTATCCT sciATAC   3.9960280 14182 female       NPC
# num_genes_expressed totalCoverage
# AGCGATAGAACGAATTCGGGCGGAGGTCGTACTGAC                4184          8522
# AGCGATAGAACGAATTCGGGCGGAGGTCTATAGCCT                5833         11334
# AGCGATAGAACGAATTCGGGTCTCGCCGCAGGACGT                3510          6868
# AGCGATAGAACGAATTCGGGTCTCGCCGGGCTCTGA                1269          2374
# AGCGATAGAACGAATTCGGTAATGATCGCCTATCCT                6789         14182



# how many express a particular gene, or how many genes are expressed by a given cell?
cds.OG <- detectGenes(cds.OG, min_expr = 0.1)
# function (cds, min_expr = NULL) 
# {
#     if (is.null(min_expr)) {
#         min_expr <- cds@lowerDetectionLimit
#     }
#     fData(cds)$num_cells_expressed <- Matrix::rowSums(exprs(cds) > 
#                                                           min_expr)
#     pData(cds)$num_genes_expressed <- Matrix::colSums(exprs(cds) > 
#                                                           min_expr)
#     cds
# }



print(head(fData(cds.OG)))
# id           gene_short_name
# chr10:100002601-100002817 chr10:100002601-100002817 chr10:100002601-100002817
# chr10:100004264-100004846 chr10:100004264-100004846 chr10:100004264-100004846
# chr10:100010337-100010585 chr10:100010337-100010585 chr10:100010337-100010585
# chr10:100015382-100016722 chr10:100015382-100016722 chr10:100015382-100016722
# chr10:100022027-100022672 chr10:100022027-100022672 chr10:100022027-100022672
# chr10:100031148-100031491 chr10:100031148-100031491 chr10:100031148-100031491
# num_cells_expressed
# chr10:100002601-100002817                  16
# chr10:100004264-100004846                  43
# chr10:100010337-100010585                  24
# chr10:100015382-100016722                 358
# chr10:100022027-100022672                  73
# chr10:100031148-100031491                  25

# The vector expressed_genes now holds the identifiers for genes expressed in at least 10 cells of the data set.
accessible_loci <- row.names(subset(fData(cds.OG),
                                    num_cells_expressed >= 10))
length(accessible_loci)
# 157052

print(head(pData(cds.OG)))










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
# 157060     2750 

dim(cds)
# Features  Samples 
# 157060     2727 










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
p3 <- ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(fill=sex), alpha = 0.4) 
p4 <- ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(color=timepoint), size=0.75) + scale_color_manual(values = tpCols)

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.coverageDistributions.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
multiplot(p1, p2, p3, p4, layout=matrix(c(1,1,2,2,3,3,4,4), nrow=2, byrow=TRUE))
dev.off()










##################################################
##################################################
##################################################
# Save filtered CDS to file

CDSmatrix = as.matrix(Biobase::exprs(cds))

dim(CDSmatrix)
# # [1] 52636  7769
# 157060   2727

CDSmatrix[1:2, 1:2]
# RowBcol2_A01_A01_RT_0 RowBcol2_A01_A01_RT_12
# ENSMUSG00000000001.4                      0                      0
# ENSMUSG00000000003.15                     0                      0

pData(cds)$sex # factor( sampleSex, levels=c("female", "male") )
pData(cds)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )

length(pData(cds)$sex) # factor( sampleSex, levels=c("female", "male") )
length(pData(cds)$timepoint) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
# 2727

colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cds)$timepoint, ".", pData(cds)$sex, sep="")

matrixfilename = file.path("../data", paste(projectName, ".countsMatrix.filtered.500UMIs.final.tsv", sep=""))
write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))

table(pData(cds)$timepoint, pData(cds)$sex)
# 500 UMIs
# female male
# d0     227  140
# d3     945  153
# d7     107  457
# d11    297  167
# NPC    234    0

summary(pData(cds)$num_genes_expressed)
# 500 UMIs
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 268     982    2788    4198    6014   22243 

# 20200914 summary stats
# @ 500 UMIs
tapply(pData(cds)$n.umi[pData(cds)$sex == "female"], pData(cds)$timepoint[pData(cds)$sex == "female"], summary)
tapply(pData(cds)$n.umi[pData(cds)$sex == "male"], pData(cds)$timepoint[pData(cds)$sex == "male"], summary)

cbind( tapply(pData(cds)$n.umi[pData(cds)$sex == "female"], pData(cds)$timepoint[pData(cds)$sex == "female"], mean),
       tapply(pData(cds)$n.umi[pData(cds)$sex == "male"], pData(cds)$timepoint[pData(cds)$sex == "male"], mean),
       tapply(pData(cds)$n.umi[pData(cds)$sex == "female"], pData(cds)$timepoint[pData(cds)$sex == "female"], median),
       tapply(pData(cds)$n.umi[pData(cds)$sex == "male"], pData(cds)$timepoint[pData(cds)$sex == "male"], median)
)

# [,1]      [,2] [,3] [,4]
# d0   7382.872  9382.986 4811 7091
# d3  10672.233 12122.549 7494 9756
# d7   9304.252  6320.425 6644 3947
# d11  6214.441  4054.174 3425 1777
# NPC  8326.620        NA 6143   NA










##################################################
##################################################
##################################################
# Clean up

sort( sapply(ls(),function(x){object.size(get(x))}) ) 

print(object.size(x=lapply(ls(), get)), units="Mb")

rm( count.matrix.OG.B, theSamplePlots, altAccess, refAccess, hitPeaksBoolIdx, count.table.OG, CDSmatrix, locus.annotations, cell.annotations )

gc()

print(object.size(x=lapply(ls(), get)), units="Mb")

# save.image() # To .RData
# load(file = ".RData")










##################################################
##################################################
##################################################
# Partion genes into chromosome-specific groups

dim(cds)
# Features  Samples 
# 157060     2727 

fData(cds)$ensChr = sub(":.*$", "", row.names(fData(cds)))
loci = sub("^.*:", "", row.names(fData(cds)))
fData(cds)$ensStart = sub("-.*$", "", loci)
fData(cds)$ensEnd = sub("^.*-", "", loci)

AgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr != "chrX" & fData(cds)$ensChr != "chrY" & ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("MT", fData(cds)$ensChr) ]
cdsA.OG = cds[AgenesNames,]
dim(cdsA.OG)
# Features  Samples 
# 152660     2727

XgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ]
cdsX.OG <- cds[XgenesNames,]
dim(cdsX.OG)
# Features  Samples
#   4364     2727 

#pData(cds)$totalCoverageV2 <- Matrix::colSums(Biobase::exprs(cds))
pData(cdsA.OG)$Aexpr <- Matrix::colSums(Biobase::exprs(cdsA.OG))
pData(cdsX.OG)$Xexpr <- Matrix::colSums(Biobase::exprs(cdsX.OG))

# Filter cdsA
upper_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) +
                       2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) -
                       2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
cdsA <- cdsA.OG[,pData(cdsA.OG)$Aexpr > lower_bound &
                    pData(cdsA.OG)$Aexpr < upper_bound]
dim(cdsA.OG)
# Features  Samples
# 152660     2727 
dim(cdsA)
# Features  Samples
# 152660     2721

# Filter cdsX
upper_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) +
                       2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) -
                       2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
cdsX <- cdsX.OG[, which(pData(cdsX.OG)$Xexpr > lower_bound &
                            pData(cdsX.OG)$Xexpr < upper_bound)]
dim(cdsX.OG)
# Features  Samples
# 4364     2727
dim(cdsX)
# Features  Samples
# 4364     2642 

table(pData(cds)$timepoint, pData(cds)$sex)
# female male
# d0     227  140
# d3     945  153
# d7     107  457
# d11    297  167
# NPC    234    0

table(pData(cdsA)$timepoint, pData(cdsA)$sex)
# female male
# d0     227  140
# d3     939  153
# d7     107  457
# d11    297  167
# NPC    234    0

table(pData(cdsX)$timepoint, pData(cdsX)$sex)
# female male
# d0     224  138
# d3     928  149
# d7     102  438
# d11    284  152
# NPC    227    0



# 20200916
cds.noNPCs = cds[,pData(cds)$timepoint != "NPC"]
dim(cds.noNPCs)
# Features  Samples 
# 157060     2493



# 20201010
cdsA.noNPCs = cdsA[,pData(cdsA)$timepoint != "NPC"]
dim(cdsA.noNPCs)
# Features  Samples 
# 152696     2487 



# 20201010
cdsX.noNPCs = cdsX[,pData(cdsX)$timepoint != "NPC"]
dim(cdsX.noNPCs)
# Features  Samples 
# 4364     2415 










##################################################
##################################################
##################################################
# Accessibility distributions faceted by genome-wide, autosome and chrX

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.accessibilityDistributions.faceted.pdf", sep=""))
pdf(PDFfilename, height=4, width=12)
print( ggplot(pData(cds.OG), aes(x=log10(totalCoverage))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("All features (", dim(cds.OG)[1], ") and unfiltered cells (n = ", dim(cds.OG)[2], ")", sep="") ) )
print( ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("All features (", dim(cds)[1], ") and filtered cells (n = ", dim(cds)[2], ")", sep="") ) )
print( ggplot(pData(cdsA.OG), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("Autosomal features (", dim(cdsA.OG)[1], ") and unfiltered cells (n = ", dim(cdsA.OG)[2], ")", sep="") ) )
print( ggplot(pData(cdsX.OG), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("chrX features (", dim(cdsX.OG)[1], ") and unfiltered cells (n = ", dim(cdsX.OG)[2], ")", sep="") ) )
print( ggplot(pData(cdsA), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("Autosomal features (", dim(cdsA)[1], ") and filtered cells (n = ", dim(cdsA)[2], ")", sep="") ) )
print( ggplot(pData(cdsX), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, nrow = 1) + ggtitle( paste("chrX features (", dim(cdsX)[1], ") and filtered cells (n = ", dim(cdsX)[2], ")", sep="") ) )
dev.off()










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- background
#
# Monocle introduced the strategy of using RNA-Seq for single cell trajectory analysis. 
# Rather than purifying cells into discrete states experimentally, 
# Monocle uses an algorithm to learn the sequence of gene expression changes each cell 
# must go through as part of a dynamic biological process. 
#
# Once it has learned the overall "trajectory" of gene expression changes, 
# Monocle can place each cell at its proper position in the trajectory. 
#
# You can then use Monocle's differential analysis toolkit to find genes 
# regulated over the course of the trajectory, as described in the section 
# Finding Genes that Change as a Function of Pseudotime. 
#
# If there are multiple outcome for the process, 
# Monocle will reconstruct a "branched" trajectory. 
# These branches correspond to cellular "decisions", 
# and Monocle provides powerful tools for identifying the genes affected by them 
# and involved in making them. You can see how to analyze branches in the section 
# Analyzing Branches in Single-Cell Trajectories. 
#
# Monocle relies on a machine learning technique called reversed graph embedding 
# to construct single-cell trajectories. You can read more about the theoretical 
# foundations of Monocle's approach in the section Theory Behind Monocle, 
# or consult the references shown at the end of the vignette. 
#
# The ordering workflow has three main steps, each of which involve a significant machine learning task.
# 
# Step 1: choosing genes that define progress
# Inferring a single-cell trajectory is a machine learning problem. 
# The first step is to select the genes Monocle will use as input for its machine learning approach. 
# This is called feature selection, and it has a major impact in the shape of the trajectory. 
#
# In single-cell RNA-Seq, genes expressed at low levels are often very noisy, 
# but some contain important information regarding the state of the cell. 
# Monocle orders cells by examining the pattern of expression of these genes across the cell population. 
# Monocle looks for genes that vary in "interesting" (i.e. not just noisy) ways, 
# and uses these to structure the data.
# Monocle provides you with a variety of tools to select genes that will yield a robust, 
# accurate, and biologically meaningful trajectory. 
# You can use these tools to either perform a completely "unsupervised" analysis, 
# in which Monocle has no forehand knowledge of which gene you consider important. 
# Alternatively, you can make use of expert knowledge in the form of genes that are already known to define 
# biolgical progress to shape Monocle's trajectory. 
# We consider this mode "semi-supervised", because Monocle will augment the markers you provide with other, 
# related genes.
#
# Step 2: reducing the dimensionality of the data
# Once we have selected the genes we will use to order the cells, 
# Monocle applies a dimensionality reduction to the data. 
# Monocle uses a recently developed algoithm called
# Reversed Graph Embedding to reduce the data's dimensionality.
#
# Step 3: ordering the cells in pseudotime
# With the expression data projected into a lower dimensional space, 
# Monocle is ready to learn the trajectory that 
# describes how cells transition from one state into another. 
# Monocle assumes that the trajectory has a tree structure, 
# with one end of it the "root", and the others the "leaves". 
# Monocle's job is to fit the best tree it can to the data. 
# This task is called manifold learning.
# A cell at the beginning of the biological process starts at the root and progresses
# along the trunk until it reaches the first branch, if there is one. 
# That cell must then choose a path, and moves further and further along the tree until it reaches a leaf. 
# A cell's pseudotime value is the distance it would have to travel to get back to the root. 










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- Semi-supervised mode

thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis",
                          ".semisupervised") )
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress
# through myogenesis. We ultimately want a set of genes that increase (or decrease)
# in expression as a function of progress through the process we're studying.
#
# Ideally, we'd like to use as little prior knowledge of the biology of the system
# under study as possible. We'd like to discover the important ordering genes from the data,
# rather than relying on literature and textbooks, because that might introduce bias in the ordering.
# We'll start with one of the simpler ways to do this, but we generally
# recommend a somewhat more sophisticated approach called "dpFeature".
#
# One effective way to isolate a set of ordering genes is to simply compare the cells collected
# at the beginning of the process to those at the end and find the differentially expressed genes,
# as described above.
# The command below will find all genes that are differentially expressed in response to the switch
# from growth medium to differentiation medium:

minExprThreshold=0.1
# minExprThreshold=1

numCellsExpressed=10
# numCellsExpressed=5

# 20200721
# Slect most variable loci
# topVarProp=0
topVarProp=0.75

cds <- detectGenes(cds, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cds), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)
# 20200924 - genome subsampled
# 101321 @ minExprThreshold=0.1 and num_cells_expressed >= 10

# 20200518
set.seed(123)
coveredLoci = sample( coveredLoci, size=floor(dsProp * length(coveredLoci)) )
length(coveredLoci)

# 20200721
# Get top most variable loci. Code from 'analyzeAllelicSciATAC.monocle2.R'
# counts matrix
CDSmatrix = as.matrix(Biobase::exprs(cds[coveredLoci]))
# all(rownames(pData(cds[coveredLoci])) == colnames(CDSmatrix))
rownames(CDSmatrix) = fData(cds[coveredLoci])$gene_short_name
CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
# # 20200617 Use mad and IQR instead of var but always 0
# CDS.var = apply(CDSmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
# CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
topQuantValue = quantile(CDS.var, p=topVarProp)
topQuantGeneBoolIdx = CDS.var >= topQuantValue
# CDSmatrix.topQuantVar = CDSmatrix[ topQuantGeneBoolIdx,  ]
sum(topQuantGeneBoolIdx)
# 39255 # top 25%
#
# 20200924 - genome subsampled
# 25331 @ minExprThreshold=0.1 and num_cells_expressed >= 10

# diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci <- differentialGeneTest(cds[coveredLoci[topQuantGeneBoolIdx],],
                                                     fullModelFormulaStr = '~timepoint',
                                                     cores = 4)

# 20200721
# Different DEqvals
# Increase stringency for DE loci
DEqval=10^-21 # ***

ordering_timepoint_DE_loci<- row.names(subset(timepoint_DE_loci, qval < DEqval))
length(ordering_timepoint_DE_loci)
# 2996 @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21 ***

# 20210331
# Save DE bins
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".DEbins.txt", sep=""))
write.table(sort(ordering_timepoint_DE_loci), 
            tablefilename, row.names = FALSE, col.names = FALSE, quote=FALSE)

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cds <- setOrderingFilter(cds, ordering_timepoint_DE_loci)

# rm(timepoint_DE_loci); gc()
print(object.size(x=lapply(ls(), get)), units="Mb")



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree', verbose=TRUE)



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cds <- orderCells(cds)



####################
####################
# "State" is just Monocle's term for the segment of the tree.
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
cds <- orderCells(cds, root_state = GM_state(cds))



####################
####################
# Tabulate results
table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]
# 2996  @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21 ***
#     d0   d3   d7  d11  NPC
# 1  367 1092  165   78    0
# 3    0    6  388  315    0
# 2    0    0   11   71  234



####################
####################
# Table of cells per branch
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.tsv", sep=""))
table_res <- rbind(
    table(pData(cds[, pData(cds)$sex == "female"])$State, pData(cds[, pData(cds)$sex == "female"])$timepoint),
    table(pData(cds[, pData(cds)$sex == "male"])$State, pData(cds[, pData(cds)$sex == "male"])$timepoint)
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.tsv", sep=""))
table_res <- rbind(
    t( t( table(pData(cds[, pData(cds)$sex == "female"])$State, pData(cds[, pData(cds)$sex == "female"])$timepoint) ) / 
        colSums( table(pData(cds[, pData(cds)$sex == "female"])$State, pData(cds[, pData(cds)$sex == "female"])$timepoint) ) ),
    t( t( table(pData(cds[, pData(cds)$sex == "male"])$State, pData(cds[, pData(cds)$sex == "male"])$timepoint) ) / 
        colSums( table(pData(cds[, pData(cds)$sex == "male"])$State, pData(cds[, pData(cds)$sex == "male"])$timepoint) ) )
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportion log 2 ratio
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.log2ratio.tsv", sep=""))
table_res <- 
    log2( ( t( t( table(pData(cds[, pData(cds)$sex == "female"])$State, pData(cds[, pData(cds)$sex == "female"])$timepoint) ) / 
                 colSums( table(pData(cds[, pData(cds)$sex == "female"])$State, pData(cds[, pData(cds)$sex == "female"])$timepoint) ) ) + 0.01 )  / 
            ( t( t( table(pData(cds[, pData(cds)$sex == "male"])$State, pData(cds[, pData(cds)$sex == "male"])$timepoint) ) / 
                     colSums( table(pData(cds[, pData(cds)$sex == "male"])$State, pData(cds[, pData(cds)$sex == "male"])$timepoint) ) ) + 0.01 ) )
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cds, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cds, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cds, color_by = "sex") + geom_point(size=2, aes(color = sex)) )
# print( plot_cell_trajectory(cds, color_by = "sample") + geom_point(size=2, aes(color = sample)) )
# print( plot_cell_trajectory(cds, color_by = "batchNum") + geom_point(size=2, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space - faceted.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# 20201005
# CDFs of pseudotime density

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cds)[, c("cell", "timepoint", "Pseudotime", "State") ]

ptData$relativePseudotime = NaN
for ( aState in ptData$State ) {
    ptData$branch_relative_Pseudotime[ ptData$State == aState ] =  zero2one( ptData$Pseudotime[ ptData$State == aState ] )
    if( min( ptData$Pseudotime[ ptData$State == aState ] ) == min(  ptData$Pseudotime ) ) {
        ptData$branch_relative_Pseudotime[ ptData$State == aState ] = 1-  ptData$branch_relative_Pseudotime[ ptData$State == aState ]
    }
}

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".pseudotime_CDFs_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(ptData, aes(branch_relative_Pseudotime, colour = State)) + 
    stat_ecdf(geom = "step") + 
    labs(title = "sci-RNA pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") 

print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()



tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".pseudotime_per_branch.tsv", sep=""))
write.table(ptData, 
            tablefilename,sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



####################
####################
# Save pData/pseudotime

head(pData(cds))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.",
                                               paste0(".minExprThreshold", minExprThreshold), 
                                               paste0(".numCellsExpressed", numCellsExpressed),
                                               ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".tsv", sep=""))
write.table(pData(cds), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# Compare to sciRNA pData/pseudotime

sciRNApDatafilename="~/proj/2019_sciRNA_analysis/results/gbonora/20190118_sciRNA_d3-d11andNPCs/analyzeSciRNAdata/analysisOutput.monocle2/trajectoryAnalysis.semisupervised/20190118_sciRNA_d3-d11andNPCs.pData.pseudotime.tsv"
sciRNApData = read.table(sciRNApDatafilename, sep="\t", header=TRUE)

head(sciRNApData)



####################
####################
# Save environment

# downsample=FALSE; thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis", ifelse(downsample, ".downsampled", ""), ".semisupervised") )
# minExprThreshold=0.1; numCellsExpressed=10; dsProp=NA; topVarProp=0.75; DEqval=0.01; DEqval=10^-18; DEqval=10^-21 
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genome.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed),
                                                 ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                                 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                                 ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ) 
) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )










##################################################
##################################################
##################################################
# 20200916 Genome-wide trajectory analysis *** EXCLUDING NPCs ***
# Constructing Single Cell Trajectories -- Semi-supervised mode

thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.noNPCs.semisupervised")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress
# through myogenesis. We ultimately want a set of genes that increase (or decrease)
# in expression as a function of progress through the process we're studying.
#
# Ideally, we'd like to use as little prior knowledge of the biology of the system
# under study as possible. We'd like to discover the important ordering genes from the data,
# rather than relying on literature and textbooks, because that might introduce bias in the ordering.
# We'll start with one of the simpler ways to do this, but we generally
# recommend a somewhat more sophisticated approach called "dpFeature".
#
# One effective way to isolate a set of ordering genes is to simply compare the cells collected
# at the beginning of the process to those at the end and find the differentially expressed genes,
# as described above.
# The command below will find all genes that are differentially expressed in response to the switch
# from growth medium to differentiation medium:

minExprThreshold=0.1
# minExprThreshold=1

numCellsExpressed=10
# numCellsExpressed=5

# 20200721
# Pre-select most variable loci
# topVarProp=0
topVarProp=0.75

cds.noNPCs <- detectGenes(cds.noNPCs, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cds.noNPCs), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)
# 153489

# 20200518
set.seed(123)
coveredLoci = sample( coveredLoci, size=floor(dsProp * length(coveredLoci)) )
length(coveredLoci)

# 20200721
# Get top most variable loci. Code from 'analyzeAllelicSciATAC.monocle2.R'
# counts matrix
cds.noNPCsmatrix = as.matrix(Biobase::exprs(cds.noNPCs[coveredLoci]))
# all(rownames(pData(cds.noNPCs[coveredLoci])) == colnames(cds.noNPCsmatrix))
rownames(cds.noNPCsmatrix) = fData(cds.noNPCs[coveredLoci])$gene_short_name
cds.noNPCs.var = apply(log10(cds.noNPCsmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
# # 20200617 Use mad and IQR instead of var but always 0
# cds.noNPCs.var = apply(cds.noNPCsmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
# cds.noNPCs.var = apply(log10(cds.noNPCsmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
topQuantValue = quantile(cds.noNPCs.var, p=topVarProp)
topQuantGeneBoolIdx = cds.noNPCs.var >= topQuantValue
# cds.noNPCsmatrix.topQuantVar = cds.noNPCsmatrix[ topQuantGeneBoolIdx,  ]
sum(topQuantGeneBoolIdx)
# 38373 # top 25%

# diff_test_res <- differentialGeneTest(cds.noNPCs[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci <- differentialGeneTest(cds.noNPCs[coveredLoci[topQuantGeneBoolIdx],],
                                          fullModelFormulaStr = '~timepoint',
                                          cores = 4)

# 20200721
# Different DEqvals
# Increase stringency for DE loci
DEqval=10^-21

ordering_timepoint_DE_loci<- row.names(subset(timepoint_DE_loci, qval < DEqval))
length(ordering_timepoint_DE_loci)
# 1248 @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21

# 20210331
# Save DE bins
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".DEbins.txt", sep=""))
write.table(sort(ordering_timepoint_DE_loci), 
            tablefilename, row.names = FALSE, col.names = FALSE, quote=FALSE)

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cds.noNPCs <- setOrderingFilter(cds.noNPCs, ordering_timepoint_DE_loci)

# rm(timepoint_DE_loci); gc()
print(object.size(x=lapply(ls(), get)), units="Mb")



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cds.noNPCs <- reduceDimension(cds.noNPCs, max_components = 2, method = 'DDRTree', verbose=TRUE)



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cds.noNPCs <- orderCells(cds.noNPCs)



####################
####################
# "State" is just Monocle's term for the segment of the tree.
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
cds.noNPCs <- orderCells(cds.noNPCs, root_state = GM_state(cds.noNPCs))



####################
####################
# Tabulate results

table_res <- table(pData(cds.noNPCs)$State, pData(cds.noNPCs)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]

# d0   d3   d7  d11  NPC
# 367 1098  564  464    0

table(pData(cds.noNPCs)$State, pData(cds.noNPCs)$timepoint, pData(cds.noNPCs)$sex)
# , ,  = female
#    d0  d3  d7 d11 NPC
# 1 227 945 107 297   0
# , ,  = male
#    d0  d3  d7 d11 NPC
# 1 140 153 457 167   0

head(pData(cds.noNPCs))



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cds.noNPCs, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cds.noNPCs, color_by = "State") + geom_point(size=2, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds.noNPCs, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cds.noNPCs, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cds.noNPCs, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cds.noNPCs, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cds.noNPCs, color_by = "sex") + geom_point(size=2, aes(color = sex)) )
# print( plot_cell_trajectory(cds.noNPCs, color_by = "sample") + geom_point(size=2, aes(color = sample)) )
# print( plot_cell_trajectory(cds.noNPCs, color_by = "batchNum") + geom_point(size=2, aes(color = batchNum)) )
dev.off()



####################
####################
# Save pData/pseudotime

head(pData(cds.noNPCs))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.",
                                               paste0(".minExprThreshold", minExprThreshold), 
                                               paste0(".numCellsExpressed", numCellsExpressed),
                                               ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".tsv", sep=""))
write.table(pData(cds.noNPCs), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# Compare to sciRNA pData/pseudotime

sciRNApDatafilename="~/proj/2019_sciRNA_analysis/results/gbonora/20190118_sciRNA_d3-d11andNPCs/analyzeSciRNAdata/analysisOutput.monocle2/trajectoryAnalysis.semisupervised/20190118_sciRNA_d3-d11andNPCs.pData.pseudotime.tsv"
sciRNApData = read.table(sciRNApDatafilename, sep="\t", header=TRUE)

head(sciRNApData)



####################
####################
# 20200928
# CDFs of pseudotime density

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cds.noNPCs)[, c("cell", "timepoint", "Pseudotime", "State", "sex") ]

ptData$branch_relative_Pseudotime = NaN
for ( aState in ptData$State ) {
    ptData$branch_relative_Pseudotime[ ptData$State == aState ] =  zero2one( ptData$Pseudotime[ ptData$State == aState ] )
    if( min( ptData$Pseudotime[ ptData$State == aState ] ) == min(  ptData$Pseudotime ) ) {
        ptData$branch_relative_Pseudotime[ ptData$State == aState ] = 1-  ptData$branch_relative_Pseudotime[ ptData$State == aState ]
    }
}

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".pseudotime_CDFs_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(ptData, aes(branch_relative_Pseudotime, colour = State)) + 
    stat_ecdf(geom = "step") + 
    labs(title = "sci-ATAC pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") 

ggplot(ptData, aes(branch_relative_Pseudotime, colour = sex)) + 
    stat_ecdf(geom = "step") + 
    labs(title = "sci-ATAC pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "sex\n") 

print( plot_cell_trajectory(cds.noNPCs, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()



tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ".pseudotime_per_branch.tsv", sep=""))
write.table(ptData, 
            tablefilename,sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



####################
####################
# Save environment

# minExprThreshold=0.1; numCellsExpressed=10; dsProp=NA; topVarProp=0.75; DEqval=10^-21
# thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.noNPCs.semisupervised")
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genome.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed),
                                                 ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                                 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                                 ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ) 
) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories using Autosomal genes -- Semi-supervised mode

thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.autosomes.semisupervised")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress
# through myogenesis. We ultimately want a set of genes that increase (or decrease)
# in expression as a function of progress through the process we're studying.
#
# Ideally, we'd like to use as little prior knowledge of the biology of the system
# under study as possible. We'd like to discover the important ordering genes from the data,
# rather than relying on literature and textbooks, because that might introduce bias in the ordering.
# We'll start with one of the simpler ways to do this, but we generally
# recommend a somewhat more sophisticated approach called "dpFeature".
#
# One effective way to isolate a set of ordering genes is to simply compare the cells collected
# at the beginning of the process to those at the end and find the differentially expressed genes,
# as described above.
# The command below will find all genes that are differentially expressed in response to the switch
# from growth medium to differentiation medium:

minExprThreshold=0.1
# minExprThreshold=1
numCellsExpressed=10

# 20200721
# Rather than downsampling, pre-select most variable loci
# topVarProp=0
topVarProp=0.75

cdsA <- detectGenes(cdsA, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cdsA), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)

# 20200518
set.seed(123)
coveredLoci = sample( coveredLoci, size=floor(dsProp * length(coveredLoci)) )
length(coveredLoci)

# 20200721
# Get top most variable loci. Code from 'analyzeAllelicSciATAC.monocle2.R'
# counts matrix
CDSmatrix = as.matrix(Biobase::exprs(cdsA[coveredLoci]))
# all(rownames(pData(cds[coveredLoci])) == colnames(CDSmatrix))
rownames(CDSmatrix) = fData(cds[coveredLoci])$gene_short_name
CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, var, na.rm=TRUE)
# # 20200617 Use mad and IQR instead of var but always 0
# CDS.var = apply(CDSmatrix, MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
# CDS.var = apply(log10(CDSmatrix + 0.1), MARGIN = 1, mad, constant = 1, high = TRUE, na.rm=TRUE)
topQuantValue = quantile(CDS.var, p=topVarProp)
topQuantGeneBoolIdx = CDS.var >= topQuantValue
# CDSmatrix.topQuantVar = CDSmatrix[ topQuantGeneBoolIdx,  ]
sum(topQuantGeneBoolIdx)
# 38159 # top 25%

# diff_test_res <- differentialGeneTest(cdsA[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci_A <- differentialGeneTest(cdsA[coveredLoci[topQuantGeneBoolIdx],],
                                            fullModelFormulaStr = '~timepoint',
                                            cores = 4)

# 20200721
# Different DEqvals
# Increase stringency for DE loci
DEqval=10^-21 ***

ordering_timepoint_DE_loci_A<- row.names(subset(timepoint_DE_loci_A, qval < DEqval))
length(ordering_timepoint_DE_loci_A)
# 2919 @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21

# 20210331
# Save DE bins
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".DEbins.txt", sep=""))
write.table(sort(ordering_timepoint_DE_loci_A), 
            tablefilename, row.names = FALSE, col.names = FALSE, quote=FALSE)

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsA <- setOrderingFilter(cdsA, ordering_timepoint_DE_loci_A)

# rm(timepoint_DE_loci_A); gc()
print(object.size(x=lapply(ls(), get)), units="Mb")



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsA <- reduceDimension(cdsA, max_components = 2, method = 'DDRTree', verbose=TRUE)



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cdsA <- orderCells(cdsA)



####################
####################
# "State" is just Monocle's term for the segment of the tree.
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}

cdsA <- orderCells(cdsA, root_state = GM_state(cdsA))



####################
####################
# Tabulate results

table_res <- table(pData(cdsA)$State, pData(cdsA)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]
# # 2919 @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21
#     d0   d3   d7  d11  NPC
# 1  364 1051   73   32    0
# 3    3   41  417  319    0
# 2    0    0   74  113  234


table(pData(cdsA)$State, pData(cdsA)$timepoint, pData(cdsA)$sex)
# 2919 @ minExprThreshold=0.1 and num_cells_expressed >= 10 and topVarProp=25% and DEqval=10^-21
# , ,  = female
# d0  d3  d7 d11 NPC
# 1 226 914  39  31   0
# 2   0   0   2  53 234
# 3   1  25  66 213   0
# 
# , ,  = male
# d0  d3  d7 d11 NPC
# 1 138 137  34   1   0
# 2   0   0  72  60   0
# 3   2  16 351 106   0

head(pData(cdsA))



####################
####################
# Table of cells per branch
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cdsA)$State, pData(cdsA)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.tsv", sep=""))
table_res <- rbind(
    table(pData(cdsA[, pData(cdsA)$sex == "female"])$State, pData(cdsA[, pData(cdsA)$sex == "female"])$timepoint),
    table(pData(cdsA[, pData(cdsA)$sex == "male"])$State, pData(cdsA[, pData(cdsA)$sex == "male"])$timepoint)
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.tsv", sep=""))
table_res <- rbind(
    t( t( table(pData(cdsA[, pData(cdsA)$sex == "female"])$State, pData(cdsA[, pData(cdsA)$sex == "female"])$timepoint) ) / 
           colSums( table(pData(cdsA[, pData(cdsA)$sex == "female"])$State, pData(cdsA[, pData(cdsA)$sex == "female"])$timepoint) ) ),
    t( t( table(pData(cdsA[, pData(cdsA)$sex == "male"])$State, pData(cdsA[, pData(cdsA)$sex == "male"])$timepoint) ) / 
           colSums( table(pData(cdsA[, pData(cdsA)$sex == "male"])$State, pData(cdsA[, pData(cdsA)$sex == "male"])$timepoint) ) )
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportion log 2 ratio
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.log2ratio.tsv", sep=""))
table_res <- 
    log2( ( t( t( table(pData(cdsA[, pData(cdsA)$sex == "female"])$State, pData(cdsA[, pData(cdsA)$sex == "female"])$timepoint) ) / 
                   colSums( table(pData(cdsA[, pData(cdsA)$sex == "female"])$State, pData(cdsA[, pData(cdsA)$sex == "female"])$timepoint) ) ) + 0.01 )  / 
              ( t( t( table(pData(cdsA[, pData(cdsA)$sex == "male"])$State, pData(cdsA[, pData(cdsA)$sex == "male"])$timepoint) ) / 
                       colSums( table(pData(cdsA[, pData(cdsA)$sex == "male"])$State, pData(cdsA[, pData(cdsA)$sex == "male"])$timepoint) ) ) + 0.01 ) )
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cdsA, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=2, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "sex") + geom_point(size=2, aes(color = sex)) )
# print( plot_cell_trajectory(cdsA, color_by = "sample") + geom_point(size=2, aes(color = sample)) )
# print( plot_cell_trajectory(cdsA, color_by = "batchNum") + geom_point(size=2, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                             ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                             ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# Save pData/pseudotime

head(pData(cdsA))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.autosomes.",
                                               paste0(".minExprThreshold", minExprThreshold), 
                                               paste0(".numCellsExpressed", numCellsExpressed),
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".tsv", sep=""))
write.table(pData(cdsA), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# Save environment

# thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.autosomes.semisupervised")
# minExprThreshold=0.1; numCellsExpressed=10; dsProp=NA; topVarProp=0.75; DEqval=0.01; DEqval=10^-21 
RdataFileName = file.path(thisOutputDir, paste0( ".RData.autosomes.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed),
                                                 ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                                 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                                 ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" )
                                                 ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories using X-linked genes -- Semi-supervised mode

thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.chrX.semisupervised")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress
# through myogenesis. We ultimately want a set of genes that increase (or decrease)
# in expression as a function of progress through the process we're studying.
#
# Ideally, we'd like to use as little prior knowledge of the biology of the system
# under study as possible. We'd like to discover the important ordering genes from the data,
# rather than relying on literature and textbooks, because that might introduce bias in the ordering.
# We'll start with one of the simpler ways to do this, but we generally
# recommend a somewhat more sophisticated approach called "dpFeature".
#
# One effective way to isolate a set of ordering genes is to simply compare the cells collected
# at the beginning of the process to those at the end and find the differentially expressed genes,
# as described above.
# The command below will find all genes that are differentially expressed in response to the switch
# from growth medium to differentiation medium:

dim(cdsX)
# Features  Samples 
# 4364     2642 

cdsX <- detectGenes(cdsX, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cdsX), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)


# diff_test_res <- differentialGeneTest(cdsX[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci_X <- differentialGeneTest(cdsX[coveredLoci,],
                                          fullModelFormulaStr = '~timepoint',
                                          cores = 4)

ordering_timepoint_DE_loci_X <- row.names(subset(timepoint_DE_loci_X, qval < 0.01))

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsX <- setOrderingFilter(cdsX, ordering_timepoint_DE_loci_X)

rm(timepoint_DE_loci_X); gc()
print(object.size(x=lapply(ls(), get)), units="Mb")



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsX <- reduceDimension(cdsX, max_components = 2, method = 'DDRTree', verbose=TRUE)



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cdsX <- orderCells(cdsX)



####################
####################
# "State" is just Monocle's term for the segment of the tree.
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}

cdsX <- orderCells(cdsX, root_state = GM_state(cdsX))



####################
####################
# Tabulate results

table_res <- table(pData(cdsX)$State, pData(cdsX)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]

# 4359 @ minExprThreshold=0.1 and numCellsExpressed >= 10
#    d0  d3  d7 d11 NPC
# 1 300 911 162  71   0
# 2  59 155 290 267  10
# 3   3  11  88  98 217


table(pData(cdsX)$State, pData(cdsX)$timepoint, pData(cdsX)$sex)

# 4359 @ minExprThreshold=0.1 and numCellsExpressed >= 10
# female
#    d0  d3  d7 d11 NPC
# 1 194 799  53  66   0
# 2  29 119  43 169  10
# 3   1  10   6  49 217
# 
# male
#    d0  d3  d7 d11 NPC
# 1 106 112 109   5   0
# 2  30  36 247  98   0
# 3   2   1  82  49   0

head(pData(cdsX))



####################
####################
# Table of cells per branch
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                              ".byTimePoint.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cdsX)$State, pData(cdsX)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.tsv", sep=""))
table_res <- rbind(
    table(pData(cdsX[, pData(cdsX)$sex == "female"])$State, pData(cdsX[, pData(cdsX)$sex == "female"])$timepoint),
    table(pData(cdsX[, pData(cdsX)$sex == "male"])$State, pData(cdsX[, pData(cdsX)$sex == "male"])$timepoint)
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.tsv", sep=""))
table_res <- rbind(
    t( t( table(pData(cdsX[, pData(cdsX)$sex == "female"])$State, pData(cdsX[, pData(cdsX)$sex == "female"])$timepoint) ) / 
           colSums( table(pData(cdsX[, pData(cdsX)$sex == "female"])$State, pData(cdsX[, pData(cdsX)$sex == "female"])$timepoint) ) ),
    t( t( table(pData(cdsX[, pData(cdsX)$sex == "male"])$State, pData(cdsX[, pData(cdsX)$sex == "male"])$timepoint) ) / 
           colSums( table(pData(cdsX[, pData(cdsX)$sex == "male"])$State, pData(cdsX[, pData(cdsX)$sex == "male"])$timepoint) ) )
)
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportion log 2 ratio
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ifelse( !is.na(dsProp), paste0(".downSampleProp", dsProp), ""), 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "" ),
                                               ifelse(topVarProp!=1, paste0(".topVarProp", round(topVarProp*100)), "" ),
                                               ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.log2ratio.tsv", sep=""))
table_res <- 
    log2( ( t( t( table(pData(cdsX[, pData(cdsX)$sex == "female"])$State, pData(cdsX[, pData(cdsX)$sex == "female"])$timepoint) ) / 
                   colSums( table(pData(cdsX[, pData(cdsX)$sex == "female"])$State, pData(cdsX[, pData(cdsX)$sex == "female"])$timepoint) ) ) + 0.01 )  / 
              ( t( t( table(pData(cdsX[, pData(cdsX)$sex == "male"])$State, pData(cdsX[, pData(cdsX)$sex == "male"])$timepoint) ) / 
                       colSums( table(pData(cdsX[, pData(cdsX)$sex == "male"])$State, pData(cdsX[, pData(cdsX)$sex == "male"])$timepoint) ) ) + 0.01 ) )
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.chrX", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cdsX, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=2, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "sex") + geom_point(size=2, aes(color = sex)) )
# print( plot_cell_trajectory(cdsX, color_by = "sample") + geom_point(size=2, aes(color = sample)) )
# print( plot_cell_trajectory(cdsX, color_by = "batchNum") + geom_point(size=2, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.chrX", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# Save environment

# minExprThreshold=0.1; numCellsExpressed=10; # No downsampling
RdataFileName = file.path(thisOutputDir, paste0( ".RData.chrX.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed)  ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )

