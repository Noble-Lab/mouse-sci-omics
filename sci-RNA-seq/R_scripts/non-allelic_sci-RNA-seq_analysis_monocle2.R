##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200325 non-allelic sci-RNA-seq analysis using monocle 2 
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

projectName = "non-allelic_sci-RNA-seq""

inputDir = "../data/non-allelic_sci-RNA-seq_count_matrices/"

mainOutputDir = "../results/non-allelic_sci-RNA-seq_monocle2/"
if ( ! file.exists( mainOutputDir ) ) {
    if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()
    
tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")

pchVal=20
pchcexVal=1.25

libs=c(
	"d0",
	"d3-11andNPCs"
)

CDSpaths = list()
CDSpaths[["d0"]] = file.path(inputDir, "d0")
CDSpaths[["d3-11andNPCs"]] = file.path(inputDir, "d3-11andNPCs")    










##################################################
##################################################
##################################################
# Load counts data and transform into CDS

cat("\n\n# ***********\nLoad and merge CDSs\n", sep="")

# Adapted from Jonathan Packer
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
	df = read.table(
		mat.path,
		col.names = c("gene.idx", "cell.idx", "count"),
		colClasses = c("integer", "integer", "integer"))
	
	gene.annotations = read.table(
		gene.annotation.path,
		col.names = c("id", "gene_short_name"),
		colClasses = c("character", "character"))
	
	cell.annotations = read.table(
		cell.annotation.path,
		col.names = c("cell", "sample"),
		colClasses = c("character", "factor"))
	
	rownames(gene.annotations) = gene.annotations$id
	rownames(cell.annotations) = cell.annotations$cell
	
	# add a dummy cell to ensure that all genes are included in the matrix
	# even if a gene isn't expressed in any cell
	df = rbind(df, data.frame(
		gene.idx = c(1, nrow(gene.annotations)),
		cell.idx = rep(nrow(cell.annotations)+1, 2),
		count = c(1, 1)))
	
	mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
	mat = mat[, 1:(ncol(mat)-1)]
	
	rownames(mat) = gene.annotations$id
	colnames(mat) = cell.annotations$cell
	
	pd = new("AnnotatedDataFrame", data = cell.annotations)
	fd = new("AnnotatedDataFrame", data = gene.annotations)
	
	cds = newCellDataSet(mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
	
	return(cds)
}

CDSperLib = list()

for (theLib in names(CDSpaths) ) {
	print(theLib)
	CDSperLib[[theLib]] = load.cds(mat.path=file.path(CDSpaths[[theLib]], "UMI.count.matrix"),
								   gene.annotation.path=file.path(CDSpaths[[theLib]], "gene.annotations"),
								   cell.annotation.path=file.path(CDSpaths[[theLib]], "cell.annotations"))
	# print(head(CDSperLib[[theLib]]))
	# print(tail(CDSperLib[[theLib]]))
	print(dim(CDSperLib[[theLib]]))
}
# [1] "d0"
# Features  Samples 
# 52636     3939 
# [1] "d3-11andNPCs"
# Features  Samples 
# 52636     6679 

for (theLib in names(CDSpaths) ) {
	print(theLib)
	totalExpr <- Matrix::colSums(exprs(CDSperLib[[theLib]]))  
	print(sum(totalExpr))
	print(summary(totalExpr))
	cat("Cells > 500: ", sum(totalExpr>500), "\n", sep="")
}
# [1] "d0"
# [1] 7654857
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     666    1379    1943    2578   23448 
# Cells > 750: 3,174
# 
# [1] "d3-11andNPCs"
# [1] 10715422
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     445    1012    1604    2050   39144 
# Cells > 750: 4,811
# 
# 4,811 + 3,174 = 7,985 ~ 7.999
        
        
        
        
        
        
        
        
        
        
##################################################
##################################################
##################################################
# Use Seurat to merge CDSs
# https://davetang.org/muse/2018/01/24/merging-two-10x-single-cell-datasets/
#
# Other approaches did not work:
#
# cds = Biobase::combine(as(cdsA,"ExpressionSet"), as(cdsB,"ExpressionSet"))
# # library(sos);findFn("combine")
#
# require("a4Base")
# Elist <- list(cdsA,cdsB)
# do.call(BiocGenerics::combine, Elist)

require(Seurat)
find.package("Seurat")
# [1] "/net/gs/vol3/software/modules-sw/R/3.5.0/Linux/RHEL6/x86_64/lib64/R/library/Seurat"

# Merge first two libraries
theLibA = names(CDSperLib)[1]
theLibB = names(CDSperLib)[2]
cdsA.Seurat <-  exportCDS(CDSperLib[[theLibA]], export_to="Seurat", export_all = TRUE)
cdsB.Seurat <-  exportCDS(CDSperLib[[theLibB]], export_to="Seurat", export_all = TRUE)

cds.Merged.Seurat <- MergeSeurat(object1 = cdsA.Seurat,
								 object2 = cdsB.Seurat,
								 add.cell.id1 = theLibA,
								 add.cell.id2 = theLibB)

# Merge other libraries
for (theLib in names(CDSperLib)[-c(1:2)] ) {
	# for (theLib in names(CDSperLib)[3] ) { # For testing
	print(theLib)
	cds.Seurat = exportCDS(CDSperLib[[theLib]], export_to="Seurat", export_all = TRUE)
	cds.Merged.Seurat.Tmp <- MergeSeurat(object1 = cds.Merged.Seurat,
										 object2 = cds.Seurat,
										 add.cell.id2 = theLib)
	# slotNames(cds.Merged.Seurat.Tmp)
	print( unique(sub("_.*", "", cds.Merged.Seurat.Tmp@cell.names)) )
	cds.Merged.Seurat = cds.Merged.Seurat.Tmp
}

unique(sub("_.*", "", cds.Merged.Seurat@cell.names))
# [1] "RowBcol2" "RowBcol3" "RowDcol5" "RowEcol6" "RowFcol7"

# Convert back to Monocle CDS
cds.OG <- importCDS(cds.Merged.Seurat, import_all=TRUE)

rm("cdsA.Seurat", "cdsB.Seurat", "cds.Seurat", "cds.Merged.Seurat.Tmp", "cds.Merged.Seurat")
gc()
        
        
        
##################################################
# Explore data

dim(pData(cds.OG))
# 10618     3
dim(fData(cds.OG))
# [1] 52636  1
dim(experimentData(cds.OG))
# NULL
dim(Biobase::exprs(cds.OG))
# 52636 10618










##################################################
##################################################
##################################################
# Fix fData ids and short names

fData(cds.OG)$id <- fData(CDSperLib[[theLibA]])$id
fData(cds.OG)$gene_short_name <- fData(CDSperLib[[theLibA]])$gene_short_name

head(fData(cds.OG))
# gene_short_name                    id
# ENSMUSG00000000001.4            Gnai3  ENSMUSG00000000001.4
# ENSMUSG00000000003.15            Pbsn ENSMUSG00000000003.15
# ENSMUSG00000000028.14           Cdc45 ENSMUSG00000000028.14
# ENSMUSG00000000031.15             H19 ENSMUSG00000000031.15
# ENSMUSG00000000037.16           Scml2 ENSMUSG00000000037.16
# ENSMUSG00000000049.11            Apoh ENSMUSG00000000049.11

dim(fData(cds.OG))
# [1] 52636     2










##################################################
##################################################
##################################################
# Add metadata

getBatches <- function(plate2IDs) {
	batchNum = rep(2, length(plate2IDs))
	batchNum[plate2IDs == "RowBcol2"] = 1
	batchNum[plate2IDs == "RowBcol3"] = 1
	return(batchNum)
}

getRTplateNumbers <- function (cellIDs, plate2IDs) {
	RTids = as.numeric(sub(".*_", "", cellIDs))
	# RTids[plate2IDs == "RowBcol2"] = 1
	# RTids[plate2IDs == "RowBcol3"] = 2
	RTids[RTids>=0 & RTids<=23] = 1
	RTids[RTids>=24 & RTids<=47] = 2
	RTids[RTids>=384 & RTids<=479] = 3
	RTids[RTids>=480 & RTids<=575] = 4
	RTids[RTids>=576 & RTids<=671] = 5
	RTids[RTids>=672 & RTids<=767] = 6
	return(RTids)
}

getRTwellNumbers <- function (cellIDs, plate2IDs) {
	RTids = as.numeric(sub(".*_", "", cellIDs))
	# RTids[plate2IDs == "RowBcol2"] = NaN
	# RTids[plate2IDs == "RowBcol3"] = NaN
	return( RTids %% 12 )
}

getSampleSex <- function (cellIDs, plate2IDs) {
	RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
	sampleSex = rep("NA", length(RTplateNumbers))
	sampleSex[RTplateNumbers==1]="female"
	sampleSex[RTplateNumbers==2]="male"
	sampleSex[RTplateNumbers>2]="female"
	sampleSex[RTplateNumbers>=5]="male"
	return(sampleSex)
}

getSampleTimePoints <- function (cellIDs, plate2IDs) {
	RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
	RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
	sampleTimePoints = rep("NA", length(RTplateNumbers))
	sampleTimePoints[RTplateNumbers==1]="d0"
	sampleTimePoints[RTplateNumbers==2]="d0"
	sampleTimePoints[RTplateNumbers==3 & RTwellNumbers<=5]="d3"
	sampleTimePoints[RTplateNumbers==3 & RTwellNumbers>5]="d7"
	sampleTimePoints[RTplateNumbers==4 & RTwellNumbers<=5]="d11"
	sampleTimePoints[RTplateNumbers==4 & RTwellNumbers>5]="NPC"
	sampleTimePoints[RTplateNumbers==5 & RTwellNumbers<=5]="d3"
	sampleTimePoints[RTplateNumbers==5 & RTwellNumbers>5]="d7"
	sampleTimePoints[RTplateNumbers==6 & RTwellNumbers<=5]="d11"
	sampleTimePoints[RTplateNumbers==6 & RTwellNumbers>5]="NPC"
	return(sampleTimePoints)
}

head(pData(cds.OG))

plate2IDs = pData(cds.OG)$sample
unique(plate2IDs)
# [1] "d0"           "d3-11andNPCs"

table(plate2IDs)
# d0 d3-11andNPCs 
# 3939         6679 

cellIDs = pData(cds.OG)$cell
cellCount = length(cellIDs)
# 10618

table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))
# d0 d3-11andNPCs 
# 3939         6679 

plate2cellIDs = rownames(pData(cds.OG))
plate2RTids = sub("_.*_", "_", plate2cellIDs)
unique(plate2RTids)

RTids = as.numeric(as.character(sub(".*RT_", "", cellIDs)))
sort(unique(RTids))
range(unique(RTids))
# 0 781

summary(RTids)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    30.0   468.0   375.7   639.0   781.0 

summary(RTids[RTids<384])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   10.00   23.00   22.94   35.00   64.00 

(781-384+1)/96
# [1] 4

RTplateBatches = getBatches(plate2IDs)
table(RTplateBatches)

RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
table(RTplateNumbers)

RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
table(RTwellNumbers)

sampleSex = getSampleSex(cellIDs, plate2IDs)
table(sampleSex)
# female   male 
# 4909   5709 

sampleTimePoints = getSampleTimePoints(cellIDs, plate2IDs)
table(sampleTimePoints)
# d0  d11   d3   d7   NA  NPC 
# 3939 1721 1779 1837   13 1329 

cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="d3",][1:3,]
# cellIDs          plate2IDs RTplateNumbers
# [1,] "A01_A01_RT_399" "He_5"    "3"
# [2,] "A01_A01_RT_458" "He_5"    "3"
# [3,] "A01_A01_RT_461" "He_5"    "3"

table(sampleTimePoints, sampleSex)
# sampleTimePoints female male
# d0    2047 1892
# d3     991  788
# d7     595 1242
# d11    676 1045
# NA       0   13
# NPC    600  729

cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="NA",]
# cellIDs         plate2IDs RTplateNumbers
# [1,] "A08_A08_RT_64" "He_6"    "64"
# [2,] "A10_A10_RT_64" "He_6"    "64"

########################################
# Remove those troublesome He_6 64 cells
########################################
cds.OG <- cds.OG[,sampleTimePoints!="NA"]

# Rerun stats
plate2IDs = pData(cds.OG)$sample
table(plate2IDs)
# d0 d3-11andNPCs 
# 3939         6666 

cellIDs = pData(cds.OG)$cell
cellCount = length(cellIDs)
# 10605

table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))
# d0 d3-11andNPCs 
# 3939         6666 

plate2cellIDs = rownames(pData(cds.OG))
plate2RTids = sub("_.*_", "_", plate2cellIDs)
unique(plate2RTids)

RTids = as.numeric(sub(".*_", "", cellIDs))
unique(RTids)
range(unique(RTids))
# [1]   0 767

RTplateBatches = getBatches(plate2IDs)
table(RTplateBatches)
# 2 
# 10605 

RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
table(RTplateNumbers)
# 1    2    3    4    5    6 
# 2047 1892 1586 1276 2030 1774 

RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
table(RTwellNumbers)
# 0    1    2    3    4    5    6    7    8    9   10   11 
# 742  961  707  971 1056 1087  883  905  834 1001  714  744 

sampleSex = getSampleSex(cellIDs, plate2IDs)
table(sampleSex)
# female   male 
# 4909   5696 

sampleTimePoints = getSampleTimePoints(cellIDs, plate2IDs)
table(sampleTimePoints)
# d0  d11   d3   d7  NPC 
# 3939 1721 1779 1837 1329  

table(sampleTimePoints, sampleSex)
# d0    2047 1892
# d3     991  788
# d7     595 1242
# d11    676 1045
# NPC    600  729

pData(cds.OG)$batchNum <- factor( RTplateBatches, levels=c("1", "2") )
pData(cds.OG)$sex <- factor( sampleSex, levels=c("female", "male") )
pData(cds.OG)$timepoint <- factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )










##################################################
##################################################
##################################################
# Add gencode chromosomal annotations

cat("\n\n# ***********\nGencode chromosomal annotations\n", sep="")

# 20200131
# Use annotations used in preprocessing pipeline
geneAnnotationsFile="../refdata/gene-models/mouse/gencode.vM15.genes.bed"
geneAnnotations = read.table(geneAnnotationsFile)
names(geneAnnotations) = c("chr", "start", "end", "id", "dummy", "strand")

geneIDs = geneAnnotations$id
length(geneIDs)
# 52636

cdsGeneIDs = fData(cds.OG)$id 
length(cdsGeneIDs)
# 52636

rangeGenesInCDS = match(cdsGeneIDs, geneIDs)
sum(!is.na(rangeGenesInCDS)) / length(rangeGenesInCDS) * 100
# 100

fData(cds.OG)$ensChr <- geneAnnotations$chr[rangeGenesInCDS]
fData(cds.OG)$ensStart <- geneAnnotations$start[rangeGenesInCDS]
fData(cds.OG)$ensEnd <- geneAnnotations$end[rangeGenesInCDS]
fData(cds.OG)$ensStrand <- geneAnnotations$strand[rangeGenesInCDS]










##################################################
##################################################
##################################################
# Explore data

# Total expression per barcode
pData(cds.OG)$totalExpr <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!
pData(cds.OG)[1:5,]

# Total expression distributions
summary(pData(cds.OG)$totalExpr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     506    1147    1732    2248   39144 

summary(log10(pData(cds.OG)$totalExpr))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.000   2.704   3.060   3.011   3.352   4.593 

# Total expression per secondary plate
for (plate2ID in unique(plate2IDs) ) {
	cat(plate2ID, "\t", sum(pData(cds.OG)$totalExpr[plate2IDs == plate2ID], na.rm=TRUE), "\n", sep="")
}
# d0	7654857
# d3-11andNPCs	10710784

# Total expression per RT (primaty) plate
for (plate1ID in unique(RTplateNumbers) ) {
	cat(plate1ID, "\t", sum(pData(cds.OG)$totalExpr[RTplateNumbers == plate1ID], na.rm=TRUE), "\n", sep="")
}
# 1	3852609
# 2	3802248
# 3	2064587
# 4	2928101
# 5	2154817
# 6	3563279

for (timePoint in unique(pData(cds.OG)$timepoint) ) {
	cat(timePoint, "\t", sum(pData(cds.OG)$totalExpr[ pData(cds.OG)$timepoint == timePoint ], na.rm=TRUE), "\n", sep="")
}
# d0	7654857
# d3	1845685
# d7	2373719
# d11	2778291
# NPC	3713089










##################################################
##################################################
##################################################
# Filter for female (F121) or male (F123)

cellTypeFilter="allSamples"
# cellTypeFilter="F121"
# cellTypeFilter="F123"

mainOutputDir =file.path(mainOutputDir, cellTypeFilter)
if ( ! file.exists( mainOutputDir ) ) {
    if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

if (cellTypeFilter == "allSamples") actualCellTypeFilter=NA
if (cellTypeFilter == "F121") actualCellTypeFilter="female"
if (cellTypeFilter == "F123") actualCellTypeFilter="male"

if(!is.na(actualCellTypeFilter)) {
    cat("\n\n# ***********\nFilter by sex", "\n", sep="")
    cds.OG = cds.OG[, pData(cds.OG)$sex==actualCellTypeFilter]
}

cellTypeFilter; dim(cds.OG)
# Features  Samples 
# 52636    10605 

table(    pData(cds.OG)$timepoint)
# d0   d3   d7  d11  NPC 
# 3939 1779 1837 1721 1329 










##################################################
##################################################
##################################################
# 20200516 
# Filter for low counts

totalCoverage = pData(cds.OG)$totalExpr 

length(totalCoverage)
# 10605

# Find biased Genes adjusted for bias
totCovDensity <- density(log10(totalCoverage))
maxTCdensity = max(totCovDensity$y)

# Allelic:
# maxTotCovDensity = totCovDensity$x[totCovDensity$y == maxTCdensity]
# print(c(maxTotCovDensity, 10^maxTotCovDensity))
# maxTotCovDensity = 160 # *** NB ***

totCovLHDensity <- density(log10(totalCoverage[totalCoverage < 200]))
lhMaxTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(200) & totCovLHDensity$x>0][which.max(totCovLHDensity$y[totCovLHDensity$x<log10(200) & totCovLHDensity$x>0])]
lhMinTotCovDensity = totCovLHDensity$x[totCovLHDensity$x<log10(200) & totCovLHDensity$x>lhMaxTotCovDensity][which.min(totCovLHDensity$y[totCovLHDensity$x<log10(200) & totCovLHDensity$x>lhMaxTotCovDensity])]
print(c(lhMinTotCovDensity, 10^lhMinTotCovDensity))

lhMaxTotCovDensity = 10^lhMaxTotCovDensity
lhMinTotCovDensity = 10^lhMinTotCovDensity

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.coverageDistribution.unfiltered.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
plot(totCovDensity, xlab="log10(total UMIs)")
abline(v=log10(lhMinTotCovDensity), lty="dashed", col="red")
# abline(v=log10(lhMaxTotCovDensity), lty="dashed", col="red")
text(log10(lhMinTotCovDensity), 0, lhMinTotCovDensity, pos=4)
dev.off()

# Cells to exclude due to low coverage
maxTotCovDensity = 200 # *** NB ***
altLowCellsBoolIdx = totalCoverage > maxTotCovDensity
sum(altLowCellsBoolIdx)

cds.OG = cds.OG[,altLowCellsBoolIdx]

dim(cds.OG)
# Features  Samples 
# 52636     9505 

table(pData(cds.OG)$timepoint)
# d0   d3   d7  d11  NPC 
# 3492 1585 1678 1565 1185 

totalCoverage = pData(cds.OG)$totalExpr 

length(totalCoverage)
# 9505

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201     676    1317    1916    2435   39144 










##################################################
##################################################
##################################################
# Normalize data

cat("\n\n# ***********\nNormalize data", "\n", sep="")

cds.OG <- estimateSizeFactors(cds.OG)
pData(cds.OG)[1:5,]
# Cell   sample Size_Factor n.umi sex
# RowBcol2_A01_A01_RT_0   A01_A01_RT_0 RowBcol2   0.7530199  1272 female
# RowBcol2_A01_A01_RT_12 A01_A01_RT_12 RowBcol2   0.3759179   635 female
# RowBcol2_A01_A01_RT_16 A01_A01_RT_16 RowBcol2   0.9418669  1591 female
# RowBcol2_A01_A01_RT_17 A01_A01_RT_17 RowBcol2   1.5883273  2683 female
# RowBcol2_A01_A01_RT_24 A01_A01_RT_24 RowBcol2   1.0608582  1792 female
# timepoint totalExpr num_genes_expressed
# RowBcol2_A01_A01_RT_0         d3      1272                 952
# RowBcol2_A01_A01_RT_12        d3       635                 528
# RowBcol2_A01_A01_RT_16        d3      1591                1154
# RowBcol2_A01_A01_RT_17        d3      2683                1720
# RowBcol2_A01_A01_RT_24        d3      1792                1235

sum(exprs(cds.OG) < 1 & exprs(cds.OG) > 0)
# [1] 0

# how many express a particular gene, or how many genes are expressed by a given cell?
cds.OG <- detectGenes(cds.OG, min_expr = 0.1) # 0.1 means exprs(cds.OG) > 0
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
# id gene_short_name num_cells_expressed
# ENSMUSG00000000001.4   ENSMUSG00000000001.4           Gnai3                 188
# ENSMUSG00000000003.15 ENSMUSG00000000003.15            Pbsn                   1
# ENSMUSG00000000028.14 ENSMUSG00000000028.14           Cdc45                 143
# ENSMUSG00000000031.15 ENSMUSG00000000031.15             H19                   6
# ENSMUSG00000000037.16 ENSMUSG00000000037.16           Scml2                 142
# ENSMUSG00000000049.11 ENSMUSG00000000049.11            Apoh                   7

# The vector expressed_genes now holds the identifiers for genes expressed
# in at least 10 cells of the data set.
expressed_genes <- row.names(subset(fData(cds.OG),
									num_cells_expressed >= 10))
length(expressed_genes)
# 19392

print(head(pData(cds.OG)))








##################################################
##################################################
##################################################
# Filter out extreme cells: cds.OG --> cds

cat("\n\n# ***********\nFilter out extreme cells: cds.OG --> cds", "\n", sep="")

# rm(cds)
upper_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) +
					   2*sd(log10(pData(cds.OG)$totalExpr)))
lower_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) -
					   2*sd(log10(pData(cds.OG)$totalExpr)))
cds <- cds.OG[,pData(cds.OG)$totalExpr > lower_bound &
				  pData(cds.OG)$totalExpr < upper_bound]

dim(cds.OG)
# Stringent:
# Features  Samples 
# 52636     9505 

dim(cds)
# Stringent:
# Features  Samples 
# 52636     9206 

head(pData(cds))

# NO! Use global size factors.
# # Probably better to run this again:
# cds <- estimateSizeFactors(cds)

pData(cds)$totalExprV2 <- Matrix::colSums(Biobase::exprs(cds))  










##################################################
##################################################
##################################################
# Once you've excluded cells that do not pass your quality control filters,
# you should verify that the expression values stored in your CellDataSet
# follow a distribution that is roughly lognormal:

p1 <- qplot(totalExpr, data=pData(cds.OG), geom="density") +
	geom_vline(xintercept = lower_bound) +
	geom_vline(xintercept = upper_bound)
p2 <- qplot(log10(totalExpr), data=pData(cds.OG), geom="density") +
	geom_vline(xintercept = log10(lower_bound)) +
	geom_vline(xintercept = log10(upper_bound))
# p3 <- qplot(log10(totalExpr), data = pData(cds), geom = "density") +
#     geom_vline(xintercept = log10(lower_bound)) +
#     geom_vline(xintercept = log10(upper_bound))
p3 <- ggplot(pData(cds), aes(x=log10(totalExpr))) + 
	geom_density(aes(fill=sex), alpha = 0.4)
p4 <- ggplot(pData(cds), aes(x=log10(totalExpr))) + 
	geom_density(aes(color=timepoint), size=0.75) + 
	scale_color_manual(values = tpCols)
# Same as above two plots:
# p5 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExprV2))) + 
#     geom_density(aes(fill=sex), alpha = 0.4) 
# p6 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExprV2))) + 
#     geom_density(aes(color=timepoint), size=0.75) + 
#     scale_color_manual(values = tpCols)

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.expressionDistributions.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
multiplot(p1, p2, p3, p4, layout=matrix(c(1,1,2,2,3,3,4,4), nrow=2, byrow=TRUE))
dev.off()

rm(p1, p2, p3, p4)










##################################################
##################################################
##################################################
# Save filtered CDS to file

if (cellTypeFilter == "allSamples") {
	CDSmatrix = as.matrix(Biobase::exprs(cds))
	
	dim(CDSmatrix)
	# 52636 10185
	
	CDSmatrix[1:2, 1:2]
	# RowBcol2_A01_A01_RT_0 RowBcol2_A01_A01_RT_12
	# ENSMUSG00000000001.4                      0                      0
	# ENSMUSG00000000003.15                     0                      0
	
	pData(cds)$sex # factor( sampleSex, levels=c("female", "male") )
	pData(cds)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	
	length(pData(cds)$sex) # factor( sampleSex, levels=c("female", "male") )
	length(pData(cds)$timepoint) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	# 9206 - stringent
	
	colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cds)$timepoint, ".", pData(cds)$sex, sep="")
	rownames(CDSmatrix) = paste(fData(cds)$ensChr, ":", fData(cds)$ensStart, "-", fData(cds)$ensEnd, "|", fData(cds)$ensStrand, "|", fData(cds)$id, "|", fData(cds)$gene_short_name, sep="")
	
	matrixfilename = file.path("../data", paste(projectName, ".UMImatrix.filtered.200UMIs.tsv", sep=""))
	# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	fwrite(as.data.frame(CDSmatrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
	
	table(pData(cds)$timepoint, pData(cds)$sex)
	# @ 200UMIs Nmasked
	# female male
	# d0    1714 1670
	# d3     881  659
	# d7     520 1126
	# d11    586  942
	# NPC    514  594
	
	summary(pData(cds)$num_genes_expressed)
	# @ 100UMIs Nmasked
	# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	# 33     377     742     910    1273    3976 
	#
	# @ 200UMIs Nmasked
	#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
	# 85.0   463.0   811.5   969.6  1323.0  3717.0 
}










##################################################
##################################################
##################################################
# Partition genes into chromosome-specific groups

cat("\n\n# ***********\nPartition cells: genome, autosomes, and chrX", "\n", sep="")

GgenesNames <- row.names(fData(cds))[ ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
cdsG.OG = cds[GgenesNames,]
dim(cdsG.OG)
# @ 200 UMIs
# Features  Samples 
# 52513     9206 

# 20191112 Exclude Y genes from autosomes
AgenesNames <- row.names(fData(cds))[ fData(cds)$ensChr != "chrX" & fData(cds)$ensChr != "chrY" & ! is.na(fData(cds)$ensChr) & ! grepl("\\.1", fData(cds)$ensChr) & ! grepl("chrM", fData(cds)$ensChr) ]
cdsA.OG = cds[AgenesNames,]
dim(cdsA.OG)
# @ 200 UMIs
# Features  Samples 
# 48329     9206 

XgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ]
cdsX.OG <- cds[XgenesNames,]
dim(cdsX.OG)
# @ 200 UMIs
# Features  Samples 
# 2614     9206 

pData(cdsG.OG)$Gexpr <- Matrix::colSums(Biobase::exprs(cdsG.OG))  
pData(cdsA.OG)$Aexpr <- Matrix::colSums(Biobase::exprs(cdsA.OG))  
pData(cdsX.OG)$Xexpr <- Matrix::colSums(Biobase::exprs(cdsX.OG))  










##################################################
##################################################
##################################################
# Filter cdsG
upper_bound <- 10^(mean(log10(pData(cdsG.OG)$Gexpr)) +
					   2*sd(log10(pData(cdsG.OG)$Gexpr)))
lower_bound <- 10^(mean(log10(pData(cdsG.OG)$Gexpr)) -
					   2*sd(log10(pData(cdsG.OG)$Gexpr)))
cdsG <- cdsG.OG[,pData(cdsG.OG)$Gexpr > lower_bound &
					pData(cdsG.OG)$Gexpr < upper_bound]

dim(cdsG)
# Features  Samples 
# 52513     9825 
#
# @ 200 UMIs
# Features  Samples 
# 52513     9021 

summary(pData(cds)$num_genes_expressed)

table(pData(cdsG)$timepoint, pData(cdsG)$sex)
# @ 200 UMIs
# female male
# d0    1682 1654
# d3     868  646
# d7     501 1102
# d11    574  925
# NPC    497  572


# 20200914 summary stats
# @ 200 UMIs
tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "female"], pData(cdsG)$timepoint[pData(cdsG)$sex == "female"], summary)
tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "male"], pData(cdsG)$timepoint[pData(cdsG)$sex == "male"], summary)

cbind( tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "female"], pData(cdsG)$timepoint[pData(cdsG)$sex == "female"], mean),
	   tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "male"], pData(cdsG)$timepoint[pData(cdsG)$sex == "male"], mean),
	   tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "female"], pData(cdsG)$timepoint[pData(cdsG)$sex == "female"], median),
	   tapply(pData(cdsG)$Gexpr[pData(cdsG)$sex == "male"], pData(cdsG)$timepoint[pData(cdsG)$sex == "male"], median)
)
# [,1]      [,2]   [,3]   [,4]
# d0  1941.279 1923.2551 1561.5 1537.5
# d3  1179.256  902.1935  909.5  625.5
# d7  1324.072 1210.4673  947.0  915.0
# d11 1793.383 1414.5297 1421.0 1028.0
# NPC 2508.909 2536.1888 2158.0 2253.5


# Filter cdsA
upper_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) +
					   2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) -
					   2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
cdsA <- cdsA.OG[,pData(cdsA.OG)$Aexpr > lower_bound &
					pData(cdsA.OG)$Aexpr < upper_bound]

dim(cdsA)
# Features  Samples 
# 48329     9824 
#
# @ 200 UMIs
# Features  Samples 
# 48329     9012 

table(pData(cdsA)$timepoint, pData(cdsA)$sex)
# female male
# d0    1844 1753
# d3     935  741
# d7     556 1190
# d11    628  994
# NPC    541  642
#
# @ 200 UMIs
# female male
# d0    1686 1652
# d3     862  646
# d7     500 1102
# d11    573  924
# NPC    496  571

# Filter cdsX
upper_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) +
					   2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) -
					   2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
cdsX <- cdsX.OG[, which(pData(cdsX.OG)$Xexpr > lower_bound &
							pData(cdsX.OG)$Xexpr < upper_bound)]

dim(cdsX)
# Features  Samples 
# 2614    9818 
#
# @ 200 UMIs
# Features  Samples 
# 2614     8867 

table(pData(cdsX)$timepoint, pData(cdsX)$sex)
# female male
# d0    1848 1772
# d3     944  723
# d7     547 1184
# d11    629  982
# NPC    551  638
#
# @ 200 UMIs
# female male
# d0    1621 1650
# d3     863  607
# d7     500 1082
# d11    570  900
# NPC    501  573










##################################################
##################################################
##################################################
# 20210520
# Save filtered CDS to file

if (cellTypeFilter == "allSamples") {
	CDSmatrix = as.matrix(Biobase::exprs(cdsG))
	
	dim(CDSmatrix)
	# 52513  9021
	
	CDSmatrix[1:2, 1:2]
	# d0_A01_B02_RT_0 d0_A01_B02_RT_12
	# ENSMUSG00000000001.4                0                0
	# ENSMUSG00000000003.15               0                0
	
	pData(cdsG)$sex # factor( sampleSex, levels=c("female", "male") )
	pData(cdsG)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	
	length(pData(cdsG)$sex) # factor( sampleSex, levels=c("female", "male") )
	length(pData(cdsG)$timepoint) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	# # 9206 - stringent @ 200UMIs
	# 9021 - stringer @ 200UMIs and after log extreme removal
	
	colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cdsG)$timepoint, ".", pData(cdsG)$sex, sep="")
	rownames(CDSmatrix) = paste(fData(cdsG)$ensChr, ":", fData(cdsG)$ensStart, "-", fData(cdsG)$ensEnd, "|", fData(cdsG)$ensStrand, "|", fData(cdsG)$id, "|", fData(cdsG)$gene_short_name, sep="")
	
	matrixfilename = file.path("../data", paste(projectName, ".UMImatrix.filtered.200UMIs.final.tsv", sep=""))
	# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	fwrite(as.data.frame(CDSmatrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
	
}   




# 20200916
cdsG.noNPCs = cdsG[,pData(cdsG)$timepoint != "NPC"]
dim(cdsG.noNPCs)
# Features  Samples 
# 52513     7952 



# 20201010
cdsA.noNPCs = cdsA[,pData(cdsA)$timepoint != "NPC"]
dim(cdsA.noNPCs)
# Features  Samples 
# 48329     7945 



# 20201010
cdsX.noNPCs = cdsX[,pData(cdsX)$timepoint != "NPC"]
dim(cdsX.noNPCs)
# Features  Samples 
#     2614     7793



# save.image() # To .RData
# load(file = ".RData")










##################################################
##################################################
##################################################
# Expression distributions faceted by genome-wide, autosome and chrX

# ncol5
# mainOutputDir = paste("~/proj/2019_sciRNA_analysis/results/gbonora/20190118_sciRNA_d3-d11andNPCs/analyzeSciRNAdata/analysisOutput.monocle2")
PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle2results.expressionDistributions.faceted.ncol5.pdf", sep=""))
pdf(PDFfilename, height=3, width=10)
ggplot(pData(cds.OG), aes(x=log10(totalExpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("All features (", dim(cds.OG)[1], ") and unfiltered cells (n = ", dim(cds.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cds), aes(x=log10(totalExpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("All features (", dim(cds)[1], ") and filtered cells (n = ", dim(cds)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
# ggplot(pData(cds), aes(x=log10(totalExprV2))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint) + ggtitle("All features again")
ggplot(pData(cdsG.OG), aes(x=log10(Gexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode genome-wide genes (", dim(cdsG.OG)[1], ") and unfiltered cells (n = ", dim(cdsG.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cdsA.OG), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode autosomal genes (", dim(cdsA.OG)[1], ") and unfiltered cells (n = ", dim(cdsA.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cdsX.OG), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode chrX genes (", dim(cdsX.OG)[1], ") and unfiltered cells (n = ", dim(cdsX.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cdsG), aes(x=log10(Gexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode genome-wide genes (", dim(cdsG)[1], ") and filtered cells (n = ", dim(cdsG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cdsA), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode autosomal genes (", dim(cdsA)[1], ") and filtered cells (n = ", dim(cdsA)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
ggplot(pData(cdsX), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode chrX genes (", dim(cdsX)[1], ") and filtered cells (n = ", dim(cdsX)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25))
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
# Constructing Single Cell Trajectories -- Semi-supervised mode for Gencode genome-wide genes only -- filtered

cat("\n\n# ***********\nSemi-supervised clustering - genome-wide", "\n", sep="")
thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis.genomewide", 
												 ".gencode.filtered.semisupervised" ) )
if ( ! file.exists( thisOutputDir ) ) {
	if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

minExprThreshold=0.1 # 18953 genes; 9825 cells; # 0.1 means exprs(cds.OG) > 0

# Stringent 200 UMIs
# 14,037 / 18,744 genes; 9021 cells; minExprThreshold=0.1



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress.

dim(cdsG)
# Stringent 200 UMIs
# Features  Samples
# 52513     9021
   
# What limit?
cdsG <- detectGenes(cdsG, min_expr = minExprThreshold)
sum( fData(cdsG)$num_cells_expressed >= 10 )
sum( pData(cdsG)$num_genes_expressed >= 10 )

expressedGenes <-  row.names(subset(fData(cdsG), num_cells_expressed >= 10))
length(expressedGenes)
# 16589 - genome - downsampled

# diff_test_res <- differentialGeneTest(cdsG[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DEG_genes <- differentialGeneTest(cdsG[expressedGenes,],
											fullModelFormulaStr = '~timepoint',
											cores = 1)
ordering_timepoint_DEG_genes<- row.names(subset(timepoint_DEG_genes, qval < 0.01))
length(ordering_timepoint_DEG_genes)
# Stringent 200 UMIs
# 14,037 / 18,744; minExprThreshold=0.1
#
# No NPCs:
# 12,200; minExprThreshold=0.1
#
# Downsampled:
# 11,922; minExprThreshold=0.1


# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsG <- setOrderingFilter(cdsG, ordering_timepoint_DEG_genes)



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsG <- reduceDimension(cdsG, max_components = 2, method = 'DDRTree')



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cdsG <- orderCells(cdsG)



####################
####################
# "State" is just Monocle's term for the segment of the tree.
GM_state <- function(cdsG){
	if (length(unique(pData(cdsG)$State)) > 1){
		T0_counts <- table(pData(cdsG)$State, pData(cdsG)$timepoint)[,"d0"]
		return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
	} else {
		return (1)
	}
}
cdsG <- orderCells(cdsG, root_state = GM_state(cdsG))



####################
####################
# Tabulate results

table_res <- table(pData(cdsG)$State, pData(cdsG)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]
# Stringent 200 UMIs ***
# minExprThreshold=0.1
# d0   d3   d7  d11  NPC
# 1 3328  781   48   35    0
# 2    0  656 1516 1346    0
# 3    8   77   39  118 1069


table(pData(cdsG)$State, pData(cdsG)$timepoint, pData(cdsG)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# , ,  = female
# d0   d3   d7  d11  NPC
# 1 1678  645   36   35    0
# 2    0  213  456  459    0
# 3    4   10    9   80  497
# 
# , ,  = male
# d0   d3   d7  d11  NPC
# 1 1650  136   12    0    0
# 2    0  443 1060  887    0
# 3    4   67   30   38  572

table(pData(cdsG)$timepoint, pData(cdsG)$sex)
# minExprThreshold=0.1
# female male
# d0    1844 1752
# d3     937  739
# d7     556 1190
# d11    628  995
# NPC    542  642



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".byTimePoint.cellCountsPerBranch.femaleVsMale.tsv", sep=""))
table_res <- rbind(
	table(pData(cdsG[, pData(cdsG)$sex == "female"])$State, pData(cdsG[, pData(cdsG)$sex == "female"])$timepoint),
	table(pData(cdsG[, pData(cdsG)$sex == "male"])$State, pData(cdsG[, pData(cdsG)$sex == "male"])$timepoint)
)
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.tsv", sep=""))
table_res <- rbind(
	t( t( table(pData(cdsG[, pData(cdsG)$sex == "female"])$State, pData(cdsG[, pData(cdsG)$sex == "female"])$timepoint) ) / 
		   colSums( table(pData(cdsG[, pData(cdsG)$sex == "female"])$State, pData(cdsG[, pData(cdsG)$sex == "female"])$timepoint) ) ),
	t( t( table(pData(cdsG[, pData(cdsG)$sex == "male"])$State, pData(cdsG[, pData(cdsG)$sex == "male"])$timepoint) ) / 
		   colSums( table(pData(cdsG[, pData(cdsG)$sex == "male"])$State, pData(cdsG[, pData(cdsG)$sex == "male"])$timepoint) ) )
)
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Table of cells per branch - male vs female proportion log 2 ratio
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".byTimePoint.cellCountsPerBranch.femaleVsMale.proportions.log2ratio.tsv", sep=""))
table_res <- 
	log2( ( t( t( table(pData(cdsG[, pData(cdsG)$sex == "female"])$State, pData(cdsG[, pData(cdsG)$sex == "female"])$timepoint) ) / 
				   colSums( table(pData(cdsG[, pData(cdsG)$sex == "female"])$State, pData(cdsG[, pData(cdsG)$sex == "female"])$timepoint) ) ) + 0.01 )  / 
			  ( t( t( table(pData(cdsG[, pData(cdsG)$sex == "male"])$State, pData(cdsG[, pData(cdsG)$sex == "male"])$timepoint) ) / 
					   colSums( table(pData(cdsG[, pData(cdsG)$sex == "male"])$State, pData(cdsG[, pData(cdsG)$sex == "male"])$timepoint) ) ) + 0.01 ) )
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
print( plot_cell_trajectory(cdsG, color_by = "State") + geom_point(size=1, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsG, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsG, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsG, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "sex") + geom_point(size=1, aes(color = sex)) )
print( plot_cell_trajectory(cdsG, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsG, color_by = "sample") + geom_point(size=1, aes(color = sample)) )
print( plot_cell_trajectory(cdsG, color_by = "batchNum") + geom_point(size=1, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsG, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsG, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# 20201005
# CDFs of pseudotime density

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cdsG)[, c("cell", "timepoint", "Pseudotime", "State") ]

ptData$relativePseudotime = NaN
for ( aState in ptData$State ) {
	ptData$branch_relative_Pseudotime[ ptData$State == aState ] =  zero2one( ptData$Pseudotime[ ptData$State == aState ] )
	if( min( ptData$Pseudotime[ ptData$State == aState ] ) == min(  ptData$Pseudotime ) ) {
		ptData$branch_relative_Pseudotime[ ptData$State == aState ] = 1-  ptData$branch_relative_Pseudotime[ ptData$State == aState ]
	}
}


PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".pseudotime_CDFs_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(ptData, aes(branch_relative_Pseudotime, colour = State)) + 
	stat_ecdf(geom = "step") + 
	labs(title = "sci-Hi-C CDP pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") 

print( plot_cell_trajectory(cdsG, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()

	
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
												".pseudotime_per_branch.tsv", sep=""))
write.table(ptData, 
			tablefilename,sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



####################
####################
# And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, 
# using your biological knowledge of the system. 
# For example, in this experiment, a highly proliferative population of progenitor cells are generating two types of post-mitotic cells. 
# So the root should have cells that express high levels of proliferation markers.
# We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:

markerGenes = c("Pou5f1", "Sox2", "Nanog", "Klf4",
				# "Sox1",  # Sox1 NPC marker!
				"Otx2", # Transiently expressed during EB differentiation
				"Sema5a", # Upreguglated during EB differentiation
				"Snhg18", # SNHG18 was demonstrated to regulate the progression of epithelial-mesenchymal transition and cytoskeleton remodeling, thereby affecting cell motility.
				"Eomes", # a NPC intermediate marker, expressed highly in the intermediate progenitor stage of the developing neuron.
				"Pax6", # responsible for instructing cells in embryonic development to differentiate into cells of the brain.
				"Map2", # NPC markers
				"Xist") # XCI regulation
				
markerGenesNames <- row.names(subset(fData(cdsA),  gene_short_name %in% markerGenes))
 
if (length(markerGenesNames) > 1) {
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByState.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	print( plot_genes_jitter(cdsG[markerGenesNames,],
							 grouping = "State",
							 min_expr = 0.1) )
	dev.off()
	
	# To confirm that the ordering is correct we can select a couple of markers of myogenic progress. 
	# Plotting these genes demonstrates that ordering looks good:
	expressedGenes <-  row.names(subset(fData(cdsG), num_cells_expressed >= 10))
	cdsG.expressed <- cdsG[expressedGenes,]
	markerGenesNames.expressed <- row.names(subset(fData(cdsG.expressed),
												   gene_short_name %in% markerGenes))
	cdsG.expressed.subset <- cdsG.expressed[markerGenesNames.expressed,]
	fData(cdsG.expressed.subset)
	# gene_short_name                   id num_cells_expressed
	# ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3                 726
	# use_for_ordering
	# ENSMUSG00000086503.3             TRUE
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	print( plot_genes_in_pseudotime(cdsG.expressed.subset, color_by = "State") + geom_point(size=0.25, aes(color = State)) )
	print( plot_genes_in_pseudotime(cdsG.expressed.subset, color_by = "timepoint") + geom_point(size=0.25, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
	print( plot_genes_in_pseudotime(cdsG.expressed.subset, color_by = "sex") + geom_point(size=0.25, aes(color = sex)) )
	print( plot_genes_in_pseudotime(cdsG.expressed.subset, color_by = "sample") + geom_point(size=0.25, aes(color = sample)) )
	dev.off()
	
}



####################
####################
# Expression kinetics plots for a smaller subset of Genes
markerGeneSubset = c("Klf4",
					 "Otx2", # Transiently expressed during EB differentiation
					 "Sema5a", # Upreguglated during EB differentiation
					 "Map2",
					 "Xist") # XCI regulation)  # NPC markers
markerGeneSubsetNames <- row.names(subset(fData(cdsG),  gene_short_name %in% markerGeneSubset))
genesInExperimentIdx = match( markerGeneSubset, fData(cdsG)$gene_short_name )
markerGeneSubsetNames = rownames(fData(cdsG))[na.omit(genesInExperimentIdx)]

markerGeneSubset = as.data.frame(markerGeneSubset)
colnames(markerGeneSubset) = "geneID"
rownames(markerGeneSubset) = markerGeneSubsetNames
markerGeneSubset = factor(markerGeneSubset$geneID, levels = c("Klf4",
															  "Otx2", # Transiently expressed during EB differentiation
															  "Sema5a", # Upreguglated during EB differentiation
															  "Map2",
															  "Xist") ) # XCI regulation)  # NPC markers

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.timepoints.geneSubset.pdf", sep=""))
pdf(PDFfilename, height=10, width=5)
print( plot_genes_in_pseudotime(cdsG[markerGeneSubsetNames,], color_by = "timepoint", panel_order=markerGeneSubset) + 
		   geom_point(size=0.25, aes(color = timepoint)) + 
		   scale_color_manual(values = tpCols) )
dev.off()

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.sex.geneSubset.pdf", sep=""))
pdf(PDFfilename, height=10, width=5)
print( plot_genes_in_pseudotime(cdsG[markerGeneSubsetNames,], color_by = "sex", panel_order=markerGeneSubset) + 
		   geom_point(size=0.25, aes(color = sex)) )
dev.off()



####################
####################
# Clustering Genes by Pseudotemporal Expression Pattern
# A common question that arises when studying time-series gene expression studies is:
# "which genes follow similar kinetic trends"?
# Monocle can help you answer this question by grouping genes that have similar trends,
# so you can analyze these groups to see what they have in common.
# Monocle provides a convenient way to visualize all pseudotime-dependent genes.
# The function plot_pseudotime_heatmap takes a CellDataSet object
# (usually containing a only subset of significant genes)
# and generates smooth expression curves much like plot_genes_in_pseudotime.
# Then, it clusters these genes and plots them using the pheatmap package.
# This allows you to visualize modules of genes that co-vary across pseudotime.

if (length(markerGenesNames) > 1) {
	
	# Takes an age if cores used b/c can't allocate memory!
	pseudoTime_diff_test_res <- differentialGeneTest(cdsG[markerGenesNames,],
													 # pseudoTime_diff_test_res <- differentialGeneTest(cdsG.expressed.subset,
													 fullModelFormulaStr = "~ sm.ns(Pseudotime)",
													 cores = 1)
	pseudoTime_sigGenesNames <- row.names(subset(pseudoTime_diff_test_res, qval < 0.1))
	
	if (length(pseudoTime_sigGenesNames) > 1) {
		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsG[pseudoTime_sigGenesNames,],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsG[markerGenesNames,],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()

		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.femaleCells.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsG[pseudoTime_sigGenesNames, pData(cdsG)$sex=="female"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsG[markerGenesNames, pData(cdsG)$sex=="female"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()

		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.femaleCells.noNPCs.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsG[pseudoTime_sigGenesNames, pData(cdsG)$sex=="female" & pData(cdsG)$timepoint != "NPC"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsG[markerGenesNames, pData(cdsG)$sex=="female" & pData(cdsG)$timepoint != "NPC"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()
		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.maleCells.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsG[pseudoTime_sigGenesNames, pData(cdsG)$sex=="male"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsG[markerGenesNames, pData(cdsG)$sex=="male"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()
		
		
	}
	
} 



####################
####################
# Save trajectory information

head(pData(cdsG))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.genome.semisupervised.tsv", sep=""))
write.table(pData(cdsG), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# save environment

# thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.genomewide.gencode.filtered.semisupervised"); minExprThreshold=0.1
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genomewide.gencode.filtered.semisupervised", 
												 ifelse(minExprThreshold!=0.1, paste0(".minExprThreshold", minExprThreshold), "" ) ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )




















##################################################
##################################################
##################################################
# 20200916 Genome-wide trajectory analysis *** EXCLUDING NPCs ***
# Constructing Single Cell Trajectories -- Semi-supervised mode for Gencode genome-wide genes only -- filtered

cat("\n\n# ***********\nSemi-supervised clustering - genome-wide - excluding NPCs", "\n", sep="")

thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.noNPCs.genomewide.gencode.filtered.semisupervised")
if ( ! file.exists( thisOutputDir ) ) {
	if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}


minExprThreshold=0.1 # 18953 genes; 9825 cells; # 0.1 means exprs(cds.OG) > 0



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress.

dim(cdsG.noNPCs)
# Features  Samples 
# 52513     7952 

# What limit?
cdsG.noNPCs <- detectGenes(cdsG.noNPCs, min_expr = minExprThreshold)
sum( fData(cdsG.noNPCs)$num_cells_expressed >= 10 )
sum( pData(cdsG.noNPCs)$num_genes_expressed >= 10 )
# 19376 7952

expressedGenes <-  row.names(subset(fData(cdsG.noNPCs), num_cells_expressed >= 10))
length(expressedGenes)

# diff_test_res <- differentialGeneTest(cdsG.noNPCs[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DEG_genes <- differentialGeneTest(cdsG.noNPCs[expressedGenes,],
											fullModelFormulaStr = '~timepoint',
											cores = 1)
ordering_timepoint_DEG_genes<- row.names(subset(timepoint_DEG_genes, qval < 0.01))
length(ordering_timepoint_DEG_genes)
# Stringent 200 UMIs
# 14,037 / 18,744; minExprThreshold=0.1
# 12,200 / 19,376; minExprThreshold=0.1 - No NPCs

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsG.noNPCs <- setOrderingFilter(cdsG.noNPCs, ordering_timepoint_DEG_genes)



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsG.noNPCs <- reduceDimension(cdsG.noNPCs, max_components = 2, method = 'DDRTree')



####################
####################
# Trajectory step 3: order cells along the trajectory
# Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
cdsG.noNPCs <- orderCells(cdsG.noNPCs)



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
cdsG.noNPCs <- orderCells(cdsG.noNPCs, root_state = GM_state(cdsG.noNPCs))



####################
####################
# Tabulate results

table_res <- table(pData(cdsG.noNPCs)$State, pData(cdsG.noNPCs)$timepoint)
table_res[order(table_res[,1], table_res[,2], decreasing=TRUE), ]

table(pData(cdsG.noNPCs)$State, pData(cdsG.noNPCs)$timepoint)
# Stringent 200 UMIs
# minExprThreshold=0.1
# No NPCs
# d0   d3   d7  d11  NPC
# 3336 1514 1603 1499    0

table(pData(cdsG.noNPCs)$State, pData(cdsG.noNPCs)$timepoint, pData(cdsG.noNPCs)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# No NPCs
# , ,  = female
# d0   d3   d7  d11  NPC
# 1 1682  868  501  574    0
# , ,  = male
# 1 1654  646 1102  925    0

table(pData(cdsG.noNPCs)$timepoint, pData(cdsG.noNPCs)$sex)
# minExprThreshold=0.1
# No NPCs
# female male
# d0    1682 1654
# d3     868  646
# d7     501 1102
# d11    574  925
# NPC      0    0



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.genomewide.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "State") + geom_point(size=1, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "sex") + geom_point(size=1, aes(color = sex)) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "sample") + geom_point(size=1, aes(color = sample)) )
print( plot_cell_trajectory(cdsG.noNPCs, color_by = "batchNum") + geom_point(size=1, aes(color = batchNum)) )
dev.off()



####################
####################
# Save trajectory information

head(pData(cdsG.noNPCs))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.genome.noNPCs.semisupervised.tsv", sep=""))
write.table(pData(cdsG.noNPCs), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# 20200928
# CDFs of pseudotime density

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cdsG.noNPCs)[, c("cell", "timepoint", "Pseudotime", "State", "sex") ]

ptData$branch_relative_Pseudotime = NaN
for ( aState in ptData$State ) {
	ptData$branch_relative_Pseudotime[ ptData$State == aState ] =  zero2one( ptData$Pseudotime[ ptData$State == aState ] )
	if( min( ptData$Pseudotime[ ptData$State == aState ] ) == min(  ptData$Pseudotime ) ) {
		ptData$branch_relative_Pseudotime[ ptData$State == aState ] = 1 -  ptData$branch_relative_Pseudotime[ ptData$State == aState ]
	}
}

if (sexCorrected) {
	ptData$branch_relative_Pseudotime = 1 -  ptData$branch_relative_Pseudotime
}

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".pseudotime_CDFs_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(ptData, aes(branch_relative_Pseudotime, colour = State)) + 
	stat_ecdf(geom = "step") + 
	labs(title = "sci-RNA pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") 

ggplot(ptData, aes(branch_relative_Pseudotime, colour = sex)) + 
	stat_ecdf(geom = "step") + 
	labs(title = "sci-RNA pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "sex\n") 

print( plot_cell_trajectory(cdsG.noNPCs, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()



tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
												".pseudotime_per_branch.tsv", sep=""))
write.table(ptData, 
			tablefilename, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



####################
####################
# save environment

# thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.noNPCs.genomewide.gencode.filtered.semisupervised"); minExprThreshold=0.1
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genomewide.gencode.noNPCs.filtered.semisupervised", 
												 ifelse(minExprThreshold!=0.1, paste0(".minExprThreshold", minExprThreshold), "" ) ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )    










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- Semi-supervised mode for Gencode autosomal genes only -- filtered

cat("\n\n# ***********\nSemi-supervised clustering - autosomes", "\n", sep="")

thisOutputDir = file.path(mainOutputDir, "trajectoryAnalysis.autosomes.gencode.filtered.semisupervised")
if ( ! file.exists( thisOutputDir ) ) {
	if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

minExprThreshold=0.1 # 18321 genes; 9824 cells; # 0.1 means exprs(cds.OG) > 0



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress.

dim(cdsA)
# Stringent 200 UMIs
# Features  Samples 
# 48329     9012 

# What limit?
cdsA <- detectGenes(cdsA, min_expr = minExprThreshold)
sum( fData(cdsA)$num_cells_expressed >= 10 )
sum( pData(cdsA)$num_genes_expressed >= 10 )

expressedGenes <-  row.names(subset(fData(cdsA), num_cells_expressed >= 10))
length(expressedGenes)

# diff_test_res <- differentialGeneTest(cdsA[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DEG_genes <- differentialGeneTest(cdsA[expressedGenes,],
											fullModelFormulaStr = '~timepoint',
											cores = 1)
ordering_timepoint_DEG_genes<- row.names(subset(timepoint_DEG_genes, qval < 0.01))
length(ordering_timepoint_DEG_genes)
# Stringent 200 UMIs
# 13,532 / 18,115 genes; minExprThreshold=0.1

# 20210331
# Save DE bins
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".DEbins.txt", sep=""))

write.table(sort(ordering_timepoint_DEG_genes), 
			tablefilename, row.names = FALSE, col.names = FALSE, quote=FALSE)

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsA <- setOrderingFilter(cdsA, ordering_timepoint_DEG_genes)



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsA <- reduceDimension(cdsA, max_components = 2, method = 'DDRTree')



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

table(pData(cdsA)$State, pData(cdsA)$timepoint)
# Stringent 200 UMIs
# minExprThreshold=0.1
# d0   d3   d7  d11  NPC
# 2 3328  830   55   36    0
# 3    0  531 1485 1308    0
# 1   10  147   62  153 1067

table(pData(cdsA)$State, pData(cdsA)$timepoint, pData(cdsA)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# female
# d0   d3   d7  d11  NPC
# 2 1681  664   40   36    0
# 3    0  175  440  434    0
# 1    5   23   20  103  496
# 
# male
# d0   d3   d7  d11  NPC
# 2 1647  166   15    0    0
# 3    0  356 1045  874    0
# 1    5  124   42   50  571

table(pData(cdsA)$timepoint, pData(cdsA)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# female male
# d0    1686 1652
# d3     862  646
# d7     500 1102
# d11    573  924
# NPC    496  571



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
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
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=1, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsA, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "sex") + geom_point(size=1, aes(color = sex)) )
print( plot_cell_trajectory(cdsA, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsA, color_by = "sample") + geom_point(size=1, aes(color = sample)) )
print( plot_cell_trajectory(cdsA, color_by = "batchNum") + geom_point(size=1, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsA, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsA, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, 
# using your biological knowledge of the system. 
# For example, in this experiment, a highly proliferative population of progenitor cells are generating two types of post-mitotic cells. 
# So the root should have cells that express high levels of proliferation markers.
# We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:

# markerGenes = c("Pou5f1", "Sox2", "Nanog", "Klf4",
#                 "T", "Eomes", "Mixl1", "Evx1",      # T NPC intermediate marker!
#                 "Cap1", "Mpped1", "Sox1", "Smtn",   # Sox1 NPC marker!
#                 "Otx2", "Sema5a", "Snhg18", "Xist",
#                 "Plzf", "Tuj1", "Pax6", "Map2" )  # NPC
markerGenes = c( "Nanog", "Klf4", "Sox2", "Pou5f1",
				# "Sox1",  # Sox1 NPC marker!
				"Otx2", # Transiently expressed during EB differentiation
				"Snhg18", # SNHG18 was demonstrated to regulate the progression of epithelial-mesenchymal transition and cytoskeleton remodeling, thereby affecting cell motility.
				"Sema5a", # Upreguglated during EB differentiation
				"Eomes", # a NPC intermediate marker, expressed highly in the intermediate progenitor stage of the developing neuron.
				"Pax6", # responsible for instructing cells in embryonic development to differentiate into cells of the brain.
				"Map2" )  # NPC markers
# "Xist", # XCI regulation

# markerGenesNames <- row.names(subset(fData(cdsA),  gene_short_name %in% markerGenes))
genesInExperimentIdx = match( markerGenes, fData(cdsA)$gene_short_name )

markerGenesNames = rownames(fData(cdsA))[na.omit(genesInExperimentIdx)]
length(markerGenesNames)

markerGeneAnnotation  <- as.data.frame( c( rep("pluripotency factors", 4), 
										c( rep("EB diff markers", 3) ),
										c( rep("NPC diff markers", 3) )
) )
colnames(markerGeneAnnotation) = c("groups")
rownames(markerGeneAnnotation) = markerGenesNames
markerGeneAnnotation$groups = factor(markerGeneAnnotation$groups, levels = c("pluripotency factors", "EB diff markers", "NPC diff markers"))

if (length(markerGenesNames) > 1) {
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByState.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	print( plot_genes_jitter(cdsA[markerGenesNames,],
							 grouping = "State",
							 min_expr = 0.1) )
	dev.off()
	
	# To confirm that the ordering is correct we can select a couple of markers of myogenic progress. 
	# Plotting these genes demonstrates that ordering looks good:
	expressedGenes <-  row.names(subset(fData(cdsA), num_cells_expressed >= 10))
	cdsA.expressed <- cdsA[expressedGenes,]
	markerGenesNames.expressed <- row.names(subset(fData(cdsA.expressed),
												   gene_short_name %in% markerGenes))
	cdsA.expressed.subset <- cdsA.expressed[markerGenesNames.expressed,]
	fData(cdsA.expressed.subset)
	# gene_short_name                   id num_cells_expressed
	# ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3                 726
	# use_for_ordering
	# ENSMUSG00000086503.3             TRUE
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	print( plot_genes_in_pseudotime(cdsA.expressed.subset, color_by = "State") + geom_point(size=0.25, aes(color = State)) )
	print( plot_genes_in_pseudotime(cdsA.expressed.subset, color_by = "timepoint") + geom_point(size=0.25, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
	print( plot_genes_in_pseudotime(cdsA.expressed.subset, color_by = "sex") + geom_point(size=0.25, aes(color = sex)) )
	print( plot_genes_in_pseudotime(cdsA.expressed.subset, color_by = "sample") + geom_point(size=0.25, aes(color = sample)) )
	dev.off()
	
}



####################
####################
# Expression kinetics plots for a smaller subset of Genes
markerGeneSubset = c("Klf4",
					 "Sox2",
					 "Otx2", # Transiently expressed during EB differentiation
					 "Sema5a", # Upreguglated during EB differentiation
					 "Map2" )  # NPC markers
markerGeneSubsetNames <- row.names(subset(fData(cdsA),  gene_short_name %in% markerGeneSubset))
genesInExperimentIdx = match( markerGeneSubset, fData(cdsA)$gene_short_name )
markerGeneSubsetNames = rownames(fData(cdsA))[na.omit(genesInExperimentIdx)]

markerGeneSubset = as.data.frame(markerGeneSubset)
colnames(markerGeneSubset) = "geneID"
rownames(markerGeneSubset) = markerGeneSubsetNames
markerGeneSubset = factor(markerGeneSubset$geneID, levels = c("Klf4",
															  "Sox2",
															  "Otx2", # Transiently expressed during EB differentiation
															  "Sema5a", # Upreguglated during EB differentiation
															  "Map2",
															  "Xist") ) # XCI regulation)  # NPC markers

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.timepoints.geneSubset.pdf", sep=""))
pdf(PDFfilename, height=10, width=5)
print( plot_genes_in_pseudotime(cdsA[markerGeneSubsetNames,], color_by = "timepoint", panel_order=markerGeneSubset) + 
		   geom_point(size=0.25, aes(color = timepoint)) + 
		   scale_color_manual(values = tpCols) )
dev.off()



####################
####################
# Clustering Genes by Pseudotemporal Expression Pattern
# A common question that arises when studying time-series gene expression studies is:
# "which genes follow similar kinetic trends"?
# Monocle can help you answer this question by grouping genes that have similar trends,
# so you can analyze these groups to see what they have in common.
# Monocle provides a convenient way to visualize all pseudotime-dependent genes.
# The function plot_pseudotime_heatmap takes a CellDataSet object
# (usually containing a only subset of significant genes)
# and generates smooth expression curves much like plot_genes_in_pseudotime.
# Then, it clusters these genes and plots them using the pheatmap package.
# This allows you to visualize modules of genes that co-vary across pseudotime.

source(file.path(mainOutputDir, "../..", "plot_pseudotime_heatmap.GB.R" ))

if (length(markerGenesNames) > 1) {
	
	# Takes an age if cores used b/c can't allocate memory!
	pseudoTime_diff_test_res <- differentialGeneTest(cdsA[markerGenesNames,],
													 # pseudoTime_diff_test_res <- differentialGeneTest(cdsA.expressed.subset,
													 fullModelFormulaStr = "~ sm.ns(Pseudotime)",
													 cores = 1)
	pseudoTime_sigGenesNames <- row.names(subset(pseudoTime_diff_test_res, qval < 0.1))
	
	if (length(pseudoTime_sigGenesNames) > 1) {
		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsA[pseudoTime_sigGenesNames,],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsA[markerGenesNames,],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()

		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.maleCells.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsA[pseudoTime_sigGenesNames, pData(cdsA)$sex=="male"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsA[markerGenesNames, pData(cdsA)$sex=="male"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap.GB(cdsA[markerGenesNames, pData(cdsA)$sex=="male"],
								cluster_rows = FALSE,
								add_annotation_row = markerGeneAnnotation,
								cores = 1,
								show_rownames = TRUE)
		dev.off()

		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.femaleCells.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsA[pseudoTime_sigGenesNames, pData(cdsA)$sex=="female"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsA[markerGenesNames, pData(cdsA)$sex=="female"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap.GB(cdsA[markerGenesNames, pData(cdsA)$sex=="female"],
								cluster_rows = FALSE,
								add_annotation_row = markerGeneAnnotation,
								cores = 1,
								show_rownames = TRUE)
		dev.off()
		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.autosomes.gencode.minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.femaleCells.noNPCs.pdf", sep=""))
		pdf(PDFfilename, height=10, width=10)
		plot_pseudotime_heatmap(cdsA[pseudoTime_sigGenesNames, pData(cdsA)$sex=="female" & pData(cdsA)$timepoint != "NPC"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		grid::grid.newpage()
		plot_pseudotime_heatmap(cdsA[markerGenesNames, pData(cdsA)$sex=="female" & pData(cdsA)$timepoint != "NPC"],
								hclust_method = "average",
								num_clusters = 4,
								cores = 1,
								show_rownames = TRUE)
		dev.off()

	}
	
}



####################
####################
# Save trajectory information

head(pData(cdsA))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.autosomes.semisupervised.tsv", sep=""))
write.table(pData(cdsA), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# Save environment

# minExprThreshold=0.1
RdataFileName = file.path(thisOutputDir, paste0( ".RData.autosomes.gencode.filtered.semisupervised", ifelse(minExprThreshold!=0.1, paste0(".minExprThreshold", minExprThreshold), "" ) ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- Semi-supervised mode for Gencode X chromosome genes only -- filtered

cat("\n\n# ***********\nSemi-supervised clustering - X-linked", "\n", sep="")

thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis.chrX.gencode.filtered.semisupervised" ) ) 
if ( ! file.exists( thisOutputDir ) ) {
	if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

minExprThreshold=0.1 # 624 genes; 8276 cells; # 0.1 means exprs(cds.OG) > 0



####################
####################
# Trajectory step 1: choose genes that define a cell's progress
# First, we must decide which genes we will use to define a cell's progress.

dim(cdsX)
# Features  Samples 
# 2614     9818 

# What limit?
cdsX <- detectGenes(cdsX, min_expr = minExprThreshold)
sum( fData(cdsX)$num_cells_expressed >= 10 )
sum( pData(cdsX)$num_genes_expressed >= 10 )

expressedGenes <-  row.names(subset(fData(cdsX), num_cells_expressed >= 10))
length(expressedGenes)

# diff_test_res <- differentialGeneTest(cdsX[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DEG_genes <- differentialGeneTest(cdsX[expressedGenes,],
											fullModelFormulaStr = '~timepoint',
											cores = 1)

ordering_timepoint_DEG_genes<- row.names(subset(timepoint_DEG_genes, qval < 0.01))
length(ordering_timepoint_DEG_genes)
# # 263 / 275 genes
# 404 / 624 genes
# 361 / 423 genes

# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cdsX <- setOrderingFilter(cdsX, ordering_timepoint_DEG_genes)

# # 20190507
# # Error in kmeans(t(Z), K, centers = centers) : 
# #     initial centers are not distinct
# # https://github.com/cole-trapnell-lab/monocle-release/issues/26



####################
####################
# Trajectory step 2: reduce data dimensionality
# Next, we will reduce the space down to one with two dimensions,
# which we will be able to easily visualize and interpret while Monocle is ordering the cells.
cdsX <- reduceDimension(cdsX, max_components = 2, method = 'DDRTree')



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

table(pData(cdsX)$State, pData(cdsX)$timepoint)
# Stringent 200 UMIs
# minExprThreshold=0.1
# d0   d3   d7  d11  NPC
# 1 3131  763   84   60    1
# 2   99  612 1474 1364   64
# 3   41   95   24   46 1009

table(pData(cdsX)$State, pData(cdsX)$timepoint, pData(cdsX)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# , ,  = female
# d0   d3   d7  d11  NPC
# 1 1588  692   44   50    0
# 2   19  128  442  479   23
# 3   14   43   14   41  478
# 
# , ,  = male
# d0   d3   d7  d11  NPC
# 1 1543   71   40   10    1
# 2   80  484 1032  885   41
# 3   27   52   10    5  531

table(pData(cdsX)$timepoint, pData(cdsX)$sex)
# Stringent 200 UMIs
# minExprThreshold=0.1
# female male
# d0    1621 1650
# d3     863  607
# d7     500 1082
# d11    570  900
# NPC    501  573



####################
####################
# Table of cells per branch - male vs female
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
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
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.chrX.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=1, aes(color = State)) )
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cdsX, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "sex") + geom_point(size=1, aes(color = sex)) )
print( plot_cell_trajectory(cdsX, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~sex, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~State, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~timepoint, nrow = 2) )
print( plot_cell_trajectory(cdsX, color_by = "sample") + geom_point(size=1, aes(color = sample)) )
print( plot_cell_trajectory(cdsX, color_by = "batchNum") + geom_point(size=1, aes(color = batchNum)) )
dev.off()



####################
####################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, 
											 ".monocle2results.trajectoryAnalysis.chrX.gencode.minExprThreshold", minExprThreshold, 
											 ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cdsX, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cdsX, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()



####################
####################
# Save trajectory information

head(pData(cdsX))

pDatafilename = file.path(thisOutputDir, paste(projectName, ".pData.pseudotime.chrX.semisupervised.tsv", sep=""))
write.table(pData(cdsX), pDatafilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
# try(system(paste("gzip -f \"", sub("~/", "$HOME/", pDatafilename), "\"", sep = "")))



####################
####################
# Save environment

# minExprThreshold=0.1 # 624 genes; 8276 cells; # 0.1 means exprs(cds.OG) > 0
RdataFileName = file.path(thisOutputDir, paste0( ".RData.chrX.gencode.filtered.semisupervised", ifelse(minExprThreshold!=0.1, paste0(".minExprThreshold", minExprThreshold), "" ) ) )
# save.image(file =  RdataFileName )
# load(file = RdataFileName )
