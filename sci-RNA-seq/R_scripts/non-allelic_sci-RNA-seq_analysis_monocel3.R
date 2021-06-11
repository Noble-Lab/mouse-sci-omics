##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200323 non-allelic sci-RNA-seq analysis using monocle 3 
#
# See https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/
#
##################################################
##################################################
##################################################



##################################################
# Set up
# 
# R/3.6.1
#
##################################################



rm(list=ls())
gc()
sessionInfo()
    









##################################################
##################################################
##################################################
# Packages and functions

require(monocle3)
require(ggplot2)
require(data.table)
require(Matrix)
# https://blog.revolutionanalytics.com/2014/07/magrittr-simplifying-r-code-with-pipes.html
require(magrittr)
require(dplyr)


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

mainOutputDir = "../results/non-allelic_sci-RNA-seq_monocle3/"
if ( ! file.exists( mainOutputDir ) ) {
    if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

# cellTypeFilter="allSamples"
# cellTypeFilter="F121"
# cellTypeFilter="F123"
mainOutputDir =file.path(mainOutputDir, cellTypeFilter)
if ( ! file.exists( mainOutputDir ) ) {
    if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()










##################################################
##################################################
##################################################
# Steps:
# 1. Load CDSs & merge.
# 2. Meta data & and cell type assertion.
# 3. Gene annotations.
# 4. Cell filtering: very high and very low.
# 5. Gene partitioning: all, Ensemble genome-wide, autosomes, and chrX and assign to CDS list.
# 6. Distribution plots.
# 7. Save environment.
    
cat("\n\n# ##########\n# ##########\nPROCESSING STARTS HERE\n# ##########\n", sep="")

tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")

tpJitterCols=c("black", "black", "black", "black", grey(0.8))
names(tpJitterCols) = c("d0", "d3", "d7", "d11", "NPC")

sexCols=c("red", "blue")
names(sexCols) = c("female", "male")

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
	
	# pd = new("AnnotatedDataFrame", data = cell.annotations)
	# fd = new("AnnotatedDataFrame", data = gene.annotations)
	
	# cds = newCellDataSet(mat, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
	cds = new_cell_data_set(expression_data = mat, cell_metadata = cell.annotations, gene_metadata = gene.annotations)
	
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
# Cells > 500: 3,174
# 
# [1] "d3-11andNPCs"
# [1] 10715422
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     445    1012    1604    2050   39144 
# Cells > 500: 4,811
# 
# 4,811 + 3,174 = 7,985 ~ 7.999










##################################################
##################################################
##################################################
# Merge CDSs

cds.OG <- combine_cds(CDSperLib, keep_all_genes=TRUE, cell_names_unique = FALSE)
# print(head(cds.OG))
# print(tail(cds.OG))
print(dim(cds.OG))

totalExpr <- Matrix::colSums(exprs(cds.OG))  
print(sum(totalExpr))
print(summary(totalExpr))

# [1] 52636 10618
# [1] 18370279
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100.0   504.2  1144.0  1730.1  2247.0 39144.0 

# pData(cds.OG)
# fData(cds.OG)
# dim(exprs(cds.OG))










##################################################
##################################################
##################################################
# Add metadata

cat("\n\n# ***********\nMeta data CDSs\n", sep="")

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
# unique(plate2IDs)
table(plate2IDs)

cellIDs = pData(cds.OG)$cell
cellCount = length(cellIDs)
cellCount

# table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))

plate2cellIDs = rownames(pData(cds.OG))
plate2RTids = pData(cds.OG)$sample
# unique(plate2RTids)

RTids = as.numeric(sub(".*_", "", cellIDs))
# unique(RTids)
range(unique(RTids))
summary(RTids)
summary(RTids[RTids<384])

(767-384+1)/96

RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
# table(RTplateNumbers)

RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
# table(RTwellNumbers)

sampleSex = getSampleSex(cellIDs, plate2IDs)
# table(sampleSex)

sampleTimePoints = getSampleTimePoints(cellIDs, plate2IDs)
# table(sampleTimePoints)

# cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="d3",][1:3,]

table(sampleTimePoints, sampleSex)
# sampleTimePoints female male
# d0    2047 1892
# d11    676 1045
# d3     991  788
# d7     595 1242
# NA       0   13
# NPC    600  729

# cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="NA",]

########################################
# Remove those troublesome He_6 64 cells
########################################
cds.OG <- cds.OG[,sampleTimePoints!="NA"]

# Rerun stats
plate2IDs = pData(cds.OG)$sample
# table(plate2IDs)

cellIDs = pData(cds.OG)$cell
cellCount = length(cellIDs)
cellCount

table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))

plate2cellIDs = rownames(pData(cds.OG))
plate2RTids = pData(cds.OG)$sample
# unique(plate2RTids)

RTids = as.numeric(sub(".*_", "", cellIDs))
# unique(RTids)
range(unique(RTids))
summary(RTids)
summary(RTids[RTids<384])

(767-384+1)/96

RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
table(RTplateNumbers)

RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
# table(RTwellNumbers)

sampleSex = getSampleSex(cellIDs, plate2IDs)
# table(sampleSex)

sampleTimePoints = getSampleTimePoints(cellIDs, plate2IDs)
# table(sampleTimePoints)

# cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="d3",][1:3,]

table(sampleTimePoints, sampleSex)
# sampleTimePoints female male
# d0    2047 1892
# d11    676 1045
# d3     991  788
# d7     595 1242
# NPC    600  729

pData(cds.OG)$sex <- factor( sampleSex, levels=c("female", "male") )
pData(cds.OG)$timepoint <- factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )

batchNums = rep(2, nrow(pData(cds.OG)))
batchNums[pData(cds.OG)$timepoint == "d0"] = 1
batchNums[pData(cds.OG)$timepoint %in% 
			  c("d3", "d7", "d11")] = 2
batchNums[pData(cds.OG)$timepoint == "NPC"] = 3

pData(cds.OG)$batch <- factor( batchNums, levels=c(1 , 2, 3) )










##################################################
##################################################
##################################################
# Add gencode chromosomal annotations

cat("\n\n# ***********\nGencode chromosomal annotations\n", sep="")

# 20200131
# Use annotations used in Packer pipeline
# geneAnnotationsFile="~/refdata/jspacker/gene-models/human-and-mouse/gencode.v27.vM15.genes.bed"
geneAnnotationsFile="~/refdata/jspacker/gene-models/mouse/gencode.vM15.genes.bed"
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
# 20210419 
# Filter for low counts

totalCoverage = pData(cds.OG)$totalExpr 

# Cells to exclude due to low coverage
maxTotCovDensity = 200 # *** NB ***
altLowCellsBoolIdx = totalCoverage > maxTotCovDensity
sum(altLowCellsBoolIdx)
# 9505

cds.OG = cds.OG[,altLowCellsBoolIdx]











##################################################
##################################################
#################################################
# Filter for female (F121) or male (F123)

cat("\n\n# ***********\nFilter by sex or not", "\n", sep="")

# cellTypeFilter="allSamples"
# cellTypeFilter="F121"
# cellTypeFilter="F123"

if (cellTypeFilter == "allSamples") actualCellTypeFilter=NA
if (cellTypeFilter == "F121") actualCellTypeFilter="female"
if (cellTypeFilter == "F123") actualCellTypeFilter="male"

if(!is.na(actualCellTypeFilter)) {
	cds.OG = cds.OG[, pData(cds.OG)$sex==actualCellTypeFilter]
}

cellTypeFilter; dim(cds.OG)
# "allSamples"
# [1] 52636  9505










##################################################
##################################################
##################################################
# Normalize the data

cat("\n\n# ***********\nNormalize the data\n", sep="")

cds.OG <- estimate_size_factors(cds.OG)

# *** No longer runs ***
# Takes an age
# cds.OG <- estimateDispersions(cds.OG)

# how many express a particular gene, or how many genes are expressed by a given cell?
cds.OG <- detect_genes(cds.OG, min_expr = 0.1)

print(head(fData(cds.OG)))
# id gene_short_name num_cells_expressed
# ENSMUSG00000000001.4   ENSMUSG00000000001.4           Gnai3                 188
# ENSMUSG00000000003.15 ENSMUSG00000000003.15            Pbsn                   1
# ENSMUSG00000000028.14 ENSMUSG00000000028.14           Cdc45                 143
# ENSMUSG00000000031.15 ENSMUSG00000000031.15             H19                   6
# ENSMUSG00000000037.16 ENSMUSG00000000037.16           Scml2                 142
# ENSMUSG00000000049.11 ENSMUSG00000000049.11            Apoh                   7

# The vector expressed_genes now holds the identifiers for genes expressed in at least 10 cells of the data set.
expressed_genes <- row.names(subset(fData(cds.OG),
									num_cells_expressed >= 10))
length(expressed_genes)
# # 19,392
# 19331

# print(head(pData(cds.OG)))










##################################################
##################################################
##################################################
# Filter out extreme cells

cat("\n\n# ***********\nFilter cells\n", sep="")

pData(cds.OG)$totalExpr <- Matrix::colSums(exprs(cds.OG))  

sum(pData(cds.OG)$totalExpr)
# 18,365,641
# UMIs for cells >= 500 UMIs

summary(pData(cds.OG)$totalExpr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100     506    1147    1732    2248   39144 

upper_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) +
					   2*sd(log10(pData(cds.OG)$totalExpr)))
lower_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) -
					   2*sd(log10(pData(cds.OG)$totalExpr)))
cds <- cds.OG[,pData(cds.OG)$totalExpr > lower_bound &
				  pData(cds.OG)$totalExpr < upper_bound]

dim(cds.OG)
# # Features  Samples 
# 52636  9505

dim(cds)
# # Features  Samples 
# 52636  9206

# NO! Use global size factors.
# # Probably better to run this again:
# cds <- estimate_size_factors(cds)

# 20191007: This should probabbly be rerun after filtering?
# how many express a particular gene, or how many genes are expressed by a given cell?
cds <- detect_genes(cds, min_expr = 0.1)

# pData(cds)$totalExprV2 <- Matrix::colSums(exprs(cds))  










##################################################
##################################################
##################################################
# Once you've excluded cells that do not pass your quality control filters, 
# you should verify that the expression values stored in your CellDataSet 
# follow a distribution that is roughly lognormal: 

p1 <- qplot(totalExpr, data=as.data.frame(pData(cds.OG)), geom="density") +
	geom_vline(xintercept = lower_bound) +
	geom_vline(xintercept = upper_bound)
p2 <- qplot(log10(totalExpr), data=as.data.frame(pData(cds.OG)), geom="density") +
	geom_vline(xintercept = log10(lower_bound)) +
	geom_vline(xintercept = log10(upper_bound))
p3 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExpr))) + 
	geom_density(aes(fill=sex), alpha = 0.4) 
p4 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExpr))) + 
	geom_density(aes(color=timepoint), size=0.75) + 
	scale_color_manual(values = tpCols)
# Same as above two plots:
# p5 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExprV2))) + 
#     geom_density(aes(fill=sex), alpha = 0.4) 
# p6 <- ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExprV2))) + 
#     geom_density(aes(color=timepoint), size=0.75) + 
#     scale_color_manual(values = tpCols)

PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle3results.expressionDistributions.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
multiplot(p1, p2, p3, p4, layout=matrix(c(1,1,2,2,3,3,4,4), nrow=2, byrow=TRUE))
dev.off()

rm(p1, p2, p3, p4)

# stop()









##################################################
##################################################
##################################################
# Partition from genome, autosomes and chrX using Ensemble annotation

cat("\n\n# ***********\nPartition cells: genome, autosomes, and chrX", "\n", sep="")

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

pData(cdsG.OG)$Gexpr <- Matrix::colSums(exprs(cdsG.OG))  
pData(cdsA.OG)$Aexpr <- Matrix::colSums(exprs(cdsA.OG))  
pData(cdsX.OG)$Xexpr <- Matrix::colSums(exprs(cdsX.OG))  

# Filter cdsG
upper_bound <- 10^(mean(log10(pData(cdsG.OG)$Gexpr)) +
					   2*sd(log10(pData(cdsG.OG)$Gexpr)))
lower_bound <- 10^(mean(log10(pData(cdsG.OG)$Gexpr)) -
					   2*sd(log10(pData(cdsG.OG)$Gexpr)))
cdsG <- cdsG.OG[,pData(cdsG.OG)$Gexpr > lower_bound &
					pData(cdsG.OG)$Gexpr < upper_bound]

dim(cdsG.OG)
dim(cdsG)
table(pData(cdsG)$timepoint, pData(cdsG)$sex)
# 52513  9021

# Filter cdsA
upper_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) +
					   2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])) -
					   2*sd(log10(pData(cdsA.OG)$Aexpr[ which(pData(cdsA.OG)$Aexpr > 0) ])))
cdsA <- cdsA.OG[,pData(cdsA.OG)$Aexpr > lower_bound &
					pData(cdsA.OG)$Aexpr < upper_bound]

dim(cdsA.OG)
dim(cdsA)
table(pData(cdsA)$timepoint, pData(cdsA)$sex)

# Filter cdsX
upper_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) +
					   2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
lower_bound <- 10^(mean(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])) -
					   2*sd(log10(pData(cdsX.OG)$Xexpr[ which(pData(cdsX.OG)$Xexpr > 0) ])))
cdsX <- cdsX.OG[, which(pData(cdsX.OG)$Xexpr > lower_bound &
							pData(cdsX.OG)$Xexpr < upper_bound)]

dim(cdsX.OG)
dim(cdsX)
table(pData(cdsX)$timepoint, pData(cdsX)$sex)

# how many express a particular gene, or how many genes are expressed by a given cell?
cdsG <- detect_genes(cdsG, min_expr = 0.1)
cdsA <- detect_genes(cdsA, min_expr = 0.1)
cdsX <- detect_genes(cdsX, min_expr = 0.1)










##################################################
##################################################
##################################################
# Expression distributions faceted by genome-wide, autosome and chrX

# ncol5
# mainOutputDir = paste("~/proj/2019_sciRNA_analysis/results/gbonora/20191105_sciRNA_mouseDiff_Nmasked/analyzeSciRNAdata/analysisOutput.monocle2")
PDFfilename = file.path(mainOutputDir, paste(projectName, ".monocle3results.expressionDistributions.faceted.ncol5.pdf", sep=""))
pdf(PDFfilename, height=3, width=10)
print( ggplot(as.data.frame(pData(cds.OG)), aes(x=log10(totalExpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("All features (", dim(cds.OG)[1], ") and unfiltered cells (n = ", dim(cds.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 2)) )
print( ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("All features (", dim(cds)[1], ") and filtered cells (n = ", dim(cds)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 2)) )
# print( ggplot(as.data.frame(pData(cds)), aes(x=log10(totalExprV2))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint) + ggtitle("All features again") )
print( ggplot(as.data.frame(pData(cdsG.OG)), aes(x=log10(Gexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode genome-wide genes (", dim(cdsG.OG)[1], ") and unfiltered cells (n = ", dim(cdsG.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
print( ggplot(as.data.frame(pData(cdsA.OG)), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode autosomal genes (", dim(cdsA.OG)[1], ") and unfiltered cells (n = ", dim(cdsA.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
print( ggplot(as.data.frame(pData(cdsX.OG)), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode chrX genes (", dim(cdsX.OG)[1], ") and unfiltered cells (n = ", dim(cdsX.OG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
print( ggplot(as.data.frame(pData(cdsG)), aes(x=log10(Gexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode genome-wide genes (", dim(cdsG)[1], ") and filtered cells (n = ", dim(cdsG)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
print( ggplot(as.data.frame(pData(cdsA)), aes(x=log10(Aexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode autosomal genes (", dim(cdsA)[1], ") and filtered cells (n = ", dim(cdsA)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
print( ggplot(as.data.frame(pData(cdsX)), aes(x=log10(Xexpr))) + geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + facet_wrap( ~ timepoint, ncol = 5) + ggtitle( paste("Gencode chrX genes (", dim(cdsX)[1], ") and filtered cells (n = ", dim(cdsX)[2], ")", sep="") ) + coord_cartesian(ylim=c(0, 1.25)) )
dev.off()










##################################################
##################################################
##################################################
# For all gene sets, analyze cells:

cat("\n\n# ***********\nFor all gene sets, analyze cells\n", sep="")

theCDSs = list()
theCDSs[["allGenes.genomewide"]] = cds
theCDSs[["gencodeGenes.genomewide"]] = cdsG
theCDSs[["gencodeGenes.autosomal"]] = cdsA
theCDSs[["gencodeGenes.chrX"]] = cdsX

rm( cdsG, cdsA, cdsX, cdsG.OG, cdsA.OG, cdsX.OG ); gc()










##################################################
##################################################
##################################################
# Save this 
save.image() # To .RData










##################################################
##################################################
##################################################
# Load saved environment

cat("\n\n# ***********\nLoad saved environment\n", sep="")

# Cell type subfolder
mainOutputDir =file.path(mainOutputDir, cellTypeFilter)
if ( ! file.exists( mainOutputDir ) ) {
	if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

load(file = ".RData")

print( names(theCDSs) )
cat( "\n" )

theCDSdescS = names(theCDSs)










##################################################
##################################################
##################################################
# Actual analysis:
# 1. Preprocessing (& batch correction) based on DE genes
# 2. UMAP dimensionality reduction only  based on DE genes
# 3. DE analysis across time course time points
# 4. Marker genes
# 5. Expression tables

    
cat("\n\n# ##########\n# ##########\nACTUAL ANALYSIS\n# ##########\n", sep="")

for (theCDSdesc in theCDSdescS) {
	# theCDSdesc = theCDSdescS 
	# theCDSdesc = "gencodeGenes.genomewide"
	
	cat("\n\n# ***********\n# ***********\nCDS: ", theCDSdesc, "\n# ***********\n", sep="")
	
	thisOutputDir = file.path( mainOutputDir, paste( "analysis.", theCDSdesc, ".filtered", sep="" ) )
	if ( ! file.exists( thisOutputDir ) ) {
		if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
	}
	# load( file =  file.path( thisOutputDir,  paste( ".RData.", theCDSdesc, ".postDE", sep="" ) ) )
	
	# if (exists("cds")) rm(cds)
	# if (exists("cds.DE")) rm(cds.DE)
	# gc()
	cds = theCDSs[[theCDSdesc]]
	print( dim(cds) )
	# 52513  9021
	
	
	
	##################################################
	# Preprocess
	
	cat("\n\n# ***********\nPreprocess cells\n", sep="")
	
	# Preprocessing
	cds <- preprocess_cds(cds, num_dim = nDims)
	# cds <- preprocess_cds(cds, num_dim = 50, method="LSI")
	
	
	
	##################################################
	# Remove batch effects

	# Try to remove batch effects
	if (is.na(batchCorrect)  || batchCorrect) {
		cat("\n\n# ***********\nBatch correct\n", sep="")
	
		if (is.na(batchCorrect)) {
			cds = align_cds(cds, num_dim = nDims, alignment_group = "sex", verbose=TRUE)
		} else {
			# Reset gene list
			cds = align_cds(cds, num_dim = nDims, alignment_group = "batch", verbose=TRUE)
		}
	
		# if you don't run preprocess first
		# # Aligning cells from different batches using Batchelor. 
		# # Please remember to cite:
		# #     Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). 'Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.' Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091
		# #
		# # Error in array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x),  : 
		# #                                                                      'data' must be of a vector type, was 'NULL'
		# traceback()
		# # 5: array(x, c(length(x), 1L), if (!is.null(names(x))) list(names(x), 
		# #                                                            NULL) else NULL)
		# # 4: as.matrix.default(preproc_res)
		# # 3: as.matrix(preproc_res)
		# # 2: batchelor::fastMNN(as.matrix(preproc_res), batch = colData(cds)[, 
		# #                                                                    alignment_group], k = alignment_k, cos.norm = FALSE, pc.input = TRUE)
		# # 1: align_cds(cds, num_dim = 50, alignment_group = "batch")
		
		
		
		###########################
		# Save filtered CDS to file

		if (cellTypeFilter == "allSamples") {
			CDSmatrix = as.matrix(exprs(cds))
		
			dim(CDSmatrix)
			# # 52636 10185
			# 52513  9021
		
			CDSmatrix[1:2, 1:2]
			# A01_B02_RT_0_d0 A01_B02_RT_12_d0
			# ENSMUSG00000000001.4                0                0
			# ENSMUSG00000000003.15               0                0
		
			pData(cds)$sex # factor( sampleSex, levels=c("female", "male") )
			pData(cds)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
		
			length(pData(cds)$sex) # factor( sampleSex, levels=c("female", "male") )
			length(pData(cds)$timepoint) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
			# 9021 
		
			colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cds)$timepoint, ".", pData(cds)$sex, sep="")
			rownames(CDSmatrix) = paste(fData(cds)$ensChr, ":", fData(cds)$ensStart, "-", fData(cds)$ensEnd, "|", fData(cds)$ensStrand, "|", fData(cds)$id, "|", fData(cds)$gene_short_name, sep="")
		
			matrixfilename = file.path(mainOutputDir, paste(projectName, ".UMImatrix.filtered.corrected.tsv", sep=""))
			# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
			fwrite(as.data.frame(CDSmatrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
			try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
		
			table(pData(cds)$timepoint, pData(cds)$sex)

			# # @ 200UMIs Nmasked Monocle2
			# # female male
			# # d0    1714 1670
			# # d3     881  659
			# # d7     520 1126
			# # d11    586  942
			# # NPC    514  594
			# 
			# # @ 200UMIs Nmasked Monocle3 - filtered
			# female male
			# d0    1682 1654
			# d3     868  646
			# d7     501 1102
			# d11    574  925
			# NPC    497  572
		
			summary(pData(cds)$num_genes_expressed)
			# # @ 200UMIs Nmasked Monocle2
			# #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
			# # 85.0   463.0   811.5   969.6  1323.0  3717.0
			# 
			# # @ 200UMIs Nmasked Monocle3 - batch corrected
			# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
			# 128.0   468.0   810.0   958.2  1309.0  3327.0 
		}
	
	}
	
	
	
	#################################
	# Run UMAP to get low dimension representation and visualize
	
	cat("\n\n# ***********\nUMAP and tSNE\n", sep="")
	
	set.seed(123)
	cds <- reduce_dimension(cds)
	# cds <- reduce_dimension(cds, reduction_method = "tSNE")
	# cds <- reduce_dimension(cds, reduction_method = "LSI")
	
	
	
	##################################################
	# Determine X:A ratios
	
	if(grepl("genomewide", theCDSdesc)) {
		
		theAutosomes = paste0("chr", as.character(1:19))
		AgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr %in% theAutosomes) ]
		XgenesNames <- row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ]
		cat(length(AgenesNames), " autosomal & ", length(XgenesNames), " X-linked genes\n", sep="")
		
		AmeanExp = Matrix::colMeans(exprs(cds[AgenesNames, ]))  
		XmeanExp = Matrix::colMeans(exprs(cds[XgenesNames, ]))  
		
		XAratio = XmeanExp/AmeanExp
		summary(XAratio)
		
		pData(cds)$XAratio.fullrange <- XAratio
		
		# truncate:
		XAratio.trunc = XAratio
		XAratio.trunc[XAratio.trunc<0.5] = 0.5
		# XAratio.trunc[XAratio.trunc>2] = 2
		XAratio.trunc[XAratio.trunc>1.75] = 1.75
		
		pData(cds)$XAratio <- XAratio.trunc
		
		
		###############
		# X:A ratio using expressed genes only
		
		# This should have already been done, but redo in case:
		cds <- detect_genes(cds, min_expr = 0.1)
		# genes expressed in 10% of cells:
		expressed_genes_10pct <- row.names(subset(fData(cds), num_cells_expressed >= 0.1 * nrow(pData(cds))))
		
		AgenesNames <- intersect( expressed_genes_10pct, row.names(fData(cds))[ which(fData(cds)$ensChr %in% theAutosomes) ])
		XgenesNames <- intersect( expressed_genes_10pct, row.names(fData(cds))[ which(fData(cds)$ensChr == "chrX") ] )
		cat(length(AgenesNames), " autosomal & ", length(XgenesNames), " X-linked genes\n", sep="")
		
		AmeanExp = Matrix::colMeans(exprs(cds[AgenesNames, ]))  
		XmeanExp = Matrix::colMeans(exprs(cds[XgenesNames, ]))  
		
		XAratio = XmeanExp/AmeanExp
		summary(XAratio)
		
		pData(cds)$XAratio.10pctGenes.fullrange <- XAratio
		
		# truncate:
		XAratio.trunc = XAratio
		# XAratio.trunc[XAratio.trunc<0.5] = 0.5
		XAratio.trunc[XAratio.trunc>2] = 2
		
		pData(cds)$XAratio.10pctGenes <- XAratio.trunc
		
	}
	
	
	
	##################################################
	# Plot results
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle3results.dimensionalityReduction.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	
	print( plot_pc_variance_explained(cds) )
	
	print( plot_cells(cds, color_cells_by = 'timepoint', reduction_method = "PCA", 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
			   geom_point(size=0.5, aes(color = timepoint)) + 
			   scale_color_manual(values = tpCols) )
	# print( plot_cells(cds, color_cells_by = 'timepoint', reduction_method = "tSNE", 
	#                   label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
	#            geom_point(size=0.5, aes(color = timepoint)) + 
	#            scale_color_manual(values = tpCols) )
	print( plot_cells(cds, color_cells_by = 'timepoint', 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
			   geom_point(size=0.5, aes(color = timepoint)) + 
			   scale_color_manual(values = tpCols) )
	# print( plot_cells(cds, color_cells_by = 'timepoint', reduction_method = "LSI", 
	#     label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
	#     geom_point(size=0.5, aes(color = timepoint)) + 
	#     scale_color_manual(values = tpCols) )
	
	if (cellTypeFilter == "allSamples") {
		print( plot_cells(cds, color_cells_by = 'sex', reduction_method = "PCA", 
						  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				   geom_point(size=0.5, aes(color = sex)) + 
				   scale_color_manual(values = sexCols) )
		# print( plot_cells(cds, color_cells_by = 'sex', reduction_method = "tSNE", 
		#                   label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
		#            geom_point(size=0.5, aes(color = sex)) + 
		#            scale_color_manual(values = sexCols) )
		print( plot_cells(cds, color_cells_by = 'sex',
						  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				   geom_point(size=0.5, aes(color = sex)) + 
				   scale_color_manual(values = sexCols) )
	}
	
	dev.off()
	
	
	
	if(grepl("genomewide", theCDSdesc)) {
		
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle3results.dimensionalityReduction.XAratios.pdf", sep=""))
		pdf(PDFfilename, height=10, width=20)
		
		print(plot_cells(cds, color_cells_by = 'XAratio', reduction_method = "PCA", 
						 label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				  geom_point(size=0.5, aes(color = XAratio)) +
				  facet_wrap(~sex, ncol=2) + ggtitle("XAratio"))
		# print(plot_cells(cds, color_cells_by = 'XAratio', reduction_method = "tSNE", 
		#                  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
		#           geom_point(size=0.5, aes(color = XAratio)) +
		#           facet_wrap(~sex, ncol=2) + ggtitle("XAratio"))
		print(plot_cells(cds, color_cells_by = 'XAratio', 
						 label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				  geom_point(size=0.5, aes(color = XAratio)) +
				  facet_wrap(~sex, ncol=2) + ggtitle("XAratio"))  
		
		print(plot_cells(cds, color_cells_by = 'XAratio.10pctGenes', reduction_method = "PCA", 
						 label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				  geom_point(size=0.5, aes(color = XAratio.10pctGenes)) +
				  facet_wrap(~sex, ncol=2) + ggtitle("XAratio (genes expressed in 10% of cells)"))
		# print(plot_cells(cds, color_cells_by = 'XAratio.10pctGenes', reduction_method = "tSNE", 
		#                  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
		#           geom_point(size=0.5, aes(color = XAratio.10pctGenes)) +
		#           facet_wrap(~sex, ncol=2) + ggtitle("XAratio (genes expressed in 10% of cells)"))
		print(plot_cells(cds, color_cells_by = 'XAratio.10pctGenes', 
						 label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				  geom_point(size=0.5, aes(color = XAratio.10pctGenes)) +
				  facet_wrap(~sex, ncol=2) + ggtitle("XAratio (genes expressed in 10% of cells)"))  
		
		dev.off()
	}

	
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle3results.dimensionalityReduction.geneExpression.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	
	if (theCDSdesc == "gencodeGenes.genomewide") {
		print( plot_cells(cds, 1, 2, color = "timepoint", genes = c("Klf4", "Sox2", "Xist", "Otx2", "Sema5a", "Map2"),
						  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) )
	} 
	
	if (theCDSdesc == "gencodeGenes.autosomal") {
		print( plot_cells(cds, 1, 2, color = "timepoint", genes = c("Klf4", "Sox2", "Otx2", "Sema5a", "Map2"),
						  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) )
	}
	
	if (theCDSdesc != "gencodeGenes.autosomal" && cellTypeFilter == "allSamples") {
		print( plot_cells(cds, 1, 2, color = "sex", genes = "Xist", 
						  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
				   facet_wrap(~sex, ncol=2) + ggtitle("Xist") )
	}
	
	dev.off()
	
	
	
	##################################################
	# Cluster cells and visualize the clusters and partitions
	# Note: Partitions are basically super-clusters
	
	cat("\n\n# ***********\nCluster\n", sep="")
	
	cds <- cluster_cells(cds)
	# cds <- cluster_cells(cds, reduction_method = "tSNE")
	# cds <- cluster_cells(cds, reduction_method = "LSI")
	
	# Error in dyn.load(file, DLLpath = DLLpath, ...) :
	#     unable to load shared object '/net/noble/vol2/home/gbonora/R/x86_64-pc-linux-gnu-library/3.6/leidenbase/libs/leidenbase.so':
	#     libiconv.so.2: cannot open shared object file: No such file or directory
	
	head(partitions(cds, reduction_method = "UMAP")) # partitions accessor to give list of partitions
	head(clusters(cds, reduction_method = "UMAP"))
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle3results.clusteringUsingUMAP.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	
	print( plot_cells(cds, color_cells_by = 'timepoint', 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) + 
			   geom_point(size=0.5, aes(color = timepoint)) + 
			   scale_color_manual(values = tpCols) )
	print( plot_cells(cds, color_cells_by="partition", group_cells_by="partition", 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) )
	print( plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster", 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=FALSE) )
	
	dev.off()
	
	
	
	##################################################
	# Learn the trajectory graph
	# Next, we will fit a principal graph within each parition using the learn_graph() function:
	
	cat("\n\n# ***********\nTrajectory analysis\n", sep="")

	# Learn A disjoint graph in each partition
	cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
	
	PDFfilename = file.path(thisOutputDir, 
							paste(projectName, 
								  ".monocle3results.trajectory.byPartition.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)
	
	print( plot_cells(cds, color_cells_by = 'timepoint', 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=TRUE) + 
			   geom_point(size=0.5, aes(color = timepoint)) + 
			   scale_color_manual(values = tpCols) )
	
	dev.off()

	# When use_partition = FALSE, a single graph is learned across all partitions.
	cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
	
	PDFfilename = file.path(thisOutputDir, 
							paste(projectName, 
								  ".monocle3results.trajectory.pdf", sep=""))
	pdf(PDFfilename, height=10, width=10)

	print( plot_cells(cds, color_cells_by = 'timepoint', 
					  label_cell_groups=FALSE, label_groups_by_cluster=FALSE, show_trajectory_graph=TRUE) + 
			   geom_point(size=0.5, aes(color = timepoint)) + 
			   scale_color_manual(values = tpCols) )
	
	dev.off()
	
	
	
	##################################################
	# Differential analysis on time-course only
	
	cat("\n\n# ***********\nDE analysis on time points only\n", sep="")
	
	tpCols.TP = tpCols[grep("d", names(tpCols))]
	cds.TP = cds[,pData(cds)$timepoint %in% names(tpCols.TP)]
	gene_fits.TP <- fit_models(cds.TP, model_formula_str = "~timepoint", cores = 1)
	
	# gene_fits.TP
	# gene_fits.TP$model[1]
	gene_fits.TP$model_summary[1]
	# $ENSMUSG00000000001.4
	# Generalized Linear Model of class 'speedglm':
	#     
	#     Call:  speedglm::speedglm(formula = model_formula, family = stats::quasipoisson(),      model = FALSE, y = FALSE, acc = 0.001) 
	# 
	# Coefficients:
	#     ------------------------------------------------------------------ 
	#     Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)    -3.111     0.1470 -21.168 5.32e-97 ***
	#     timepointd3     1.034     0.1951   5.300 1.18e-07 ***
	#     timepointd7     1.306     0.1834   7.120 1.16e-12 ***
	#     timepointd11    1.220     0.1890   6.458 1.12e-10 ***
	#     
	#     ------------------------------------------------------------------- 
	#     Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
	# 
	# ---
	#     null df: 8634; null deviance: 5683.07;
	# residuals df: 8631; residuals deviance: 5438.35;
	# # obs.: 8635; # non-zero weighted obs.: 8635;
	# AIC: NA; log Likelihood: NA;
	# RSS: 29827.7; dispersion: 3.455879; iterations: 6;
	# rank: 4; max tolerance: 2.3e-10; convergence: TRUE.
	
	# evaluate_fits(gene_fits.TP)
	
	fit_coefs.TP <- coefficient_table(gene_fits.TP)
	head(fit_coefs.TP)
	# # A tibble: 6 x 18
	# id    gene_short_name ensChr ensStart ensEnd ensStrand num_cells_expre… model
	# <chr> <chr>           <chr>     <int>  <int> <chr>                <int> <nam>
	#     1 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 2 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 3 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 4 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 5 ENSM… Pbsn            X        7.78e7 7.79e7 -                        5 <spe…
	# 6 ENSM… Pbsn            X        7.78e7 7.79e7 -                        5 <spe…
	# # … with 10 more variables: model_summary <named list>, status <chr>,
	# #   term <chr>, estimate <dbl>, std_err <dbl>, test_val <dbl>, p_value <dbl>,
	# #   normalized_effect <dbl>, model_component <chr>, q_value <dbl>
	
	tp_terms.TP <- fit_coefs.TP %>% filter(grepl("timepoint.*", term))
	head(tp_terms.TP)
	# # A tibble: 6 x 18
	# id    gene_short_name ensChr ensStart ensEnd ensStrand num_cells_expre… model
	# <chr> <chr>           <chr>     <int>  <int> <chr>                <int> <lis>
	#     1 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 2 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 3 ENSM… Gnai3           3        1.08e8 1.08e8 -                      921 <spe…
	# 4 ENSM… Pbsn            X        7.78e7 7.79e7 -                        5 <spe…
	# 5 ENSM… Pbsn            X        7.78e7 7.79e7 -                        5 <spe…
	# 6 ENSM… Pbsn            X        7.78e7 7.79e7 -                        5 <spe…
	# # … with 10 more variables: model_summary <list>, status <chr>, term <chr>,
	# #   estimate <dbl>, std_err <dbl>, test_val <dbl>, p_value <dbl>,
	# #   normalized_effect <dbl>, model_component <chr>, q_value <dbl>
	
	# range(tp_terms$q_value)
	# range(tp_terms.TP$q_value)
	# range(tp_terms.TP$p_value)
	
	# *** q_value < 0.01 ***
	tp_terms.sig.TP = tp_terms.TP %>% filter (q_value < 0.01) %>%
		dplyr::select(gene_short_name, term, q_value, estimate)
	tp_terms.sig.TP = tp_terms.sig.TP[order(tp_terms.sig.TP$q_value),]
	dim(tp_terms.sig.TP)
	# 5073    4
	
	genes.sig.TP = unique(tp_terms.sig.TP$gene_short_name)
	head(genes.sig.TP)
	# [1] "Cdk8"    "Camk1d"  "Airn"    "Gm42418" "Col4a1"  "Cmss1" 
	
	# Get min p-value per gene id
	cds_subset.TP <- cds.TP[rowData(cds.TP)$gene_short_name %in% genes.sig.TP,]
	# head(cds_subset.TP)
	# class: cell_data_set 
	# dim: 6 8635 
	# metadata(2): cds_version citations
	# assays(1): counts
	# rownames(6): ENSMUSG00000000001.4 ENSMUSG00000000031.15 ...
	# ENSMUSG00000000149.10 ENSMUSG00000000184.12
	# rowData names(7): id gene_short_name ... ensStrand num_cells_expressed
	# colnames(8635): A01_B02_RT_0_d0 A01_B02_RT_12_d0 ...
	# H12_F07_RT_735_d3-11andNPCs H12_F07_RT_759_d3-11andNPCs
	# colData names(10): cell sample ... totalExprV2 Gexpr
	# reducedDimNames(3): PCA UMAP tSNE
	# spikeNames(0):
	dim(cds_subset.TP)
	# 2656 8635
	
	# https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2#11889798
	# Let's create an additional data frame to hold the text annotations:
	geneQ.TP = c()
	for(theGene in genes.sig.TP) {
		minQ = min(tp_terms.sig.TP$q_value[tp_terms.sig.TP$gene_short_name == theGene])
		geneQ.TP = rbind(geneQ.TP, c(theGene, minQ))
	}
	geneQ.TP = as.data.frame(geneQ.TP)
	colnames(geneQ.TP) = c("gene_short_name", "q_value")
	geneQ.TP$q_value = paste("p-val=", signif(as.numeric(as.character(geneQ.TP$q_value)), 3), sep="")
	dim(geneQ.TP)
	# 2655    2
	
	# Plot DE gene distributions
	PDFfilename = file.path(thisOutputDir, 
							paste(projectName, 
								  ".monocle3results.DEgenes.TP.qval0.01.justViolins.pdf", sep=""))
	pdf(PDFfilename, height=10, width=15)
	
	i = 1
	while (i <= min(nrow(geneQ.TP), 500)) {
		# while (i <=  20) {
		toi = ifelse((i+9)> nrow(geneQ.TP),  nrow(geneQ.TP), i+9 )
		
		genes.sig.toPlot = genes.sig.TP[i:toi]
		cds_subset.toPlot <- cds_subset.TP[na.omit(match(genes.sig.toPlot, rowData(cds_subset.TP)$gene_short_name)), ]
		geneQ.TP.toPlot <- geneQ.TP[na.omit(match(genes.sig.toPlot, geneQ.TP$gene_short_name)), ]
		# rowData(cds_subset.toPlot); geneQ.TP.toPlot
		fData(cds_subset.toPlot)$q_value = geneQ.TP.toPlot$q_value
		# Order facets by q_value
		# https://stackoverflow.com/questions/49899417/how-to-change-the-order-of-ggplot2-facets
		# No need for a seperate geom_text call and concern yourself with font size, style and text position either
		q_valueLabels = paste(fData(cds_subset.toPlot)$gene_short_name, "\n", 
							  fData(cds_subset.toPlot)$q_value, "\n",
							  "n=", fData(cds_subset.toPlot)$num_cells_expressed, sep="")
		fData(cds_subset.toPlot)$gene_short_name = factor(q_valueLabels,
														  levels=q_valueLabels,
														  ordered=TRUE)
		fData(cds_subset.toPlot)$q_value = geneQ.TP.toPlot$q_value
		
		cat(i, "-", toi, ": ", paste(dim(cds_subset.toPlot), collapse="-"), "\n", sep="")
		# Error: cds_subset has more than 100 genes - pass only the subset of the CDS to be plotted.
		print( plot_genes_violin(cds_subset.toPlot, group_cells_by="timepoint", ncol=5, nrow=2) + 
				   geom_violin(scale = "width", aes(fill = timepoint)) + 
				   scale_fill_manual(values = tpCols) +
				   #http://www.sthda.com/english/forum/topic-26+ggplot2-changing-geom_jitter-in-specific-colors.php
				   # geom_jitter(height = 0, width = 0.1, size=0.05, aes(color = timepoint), show.legend=FALSE) +
				   # scale_color_manual(values = tpJitterCols) +
				   # # For some reason using a separate df cause issue w/ multiple over printing
				   # # geom_text( data = geneQ.TP.toPlot, 
				   # geom_text(
				   #     mapping = aes(x = Inf, y = Inf, label = q_value),
				   #     hjust = 0.5, vjust = 1, 
				   #     size = 1*ggplot2:::.pt, family=theme_get()$text$family
				   #     # https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
				   #     # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
			   #     # https://stackoverflow.com/questions/34610165/what-is-the-default-font-for-ggplot2
			   # ) +
			   # guides(colour = guide_legend(override.aes = list(size=2))) +
			   theme(axis.text.x=element_text(angle=45, hjust=1)) + 
				   ggtitle(paste("Top ", i, " - ", (i+9), " DE genes", sep="")) )
		i = i + 11
	}
	# rm("genes.sig.toPlot")
	
	dev.off()
	
	# Save 'DE' table
	geneQcount.TP = c()
	for(theGene in genes.sig.TP) {
		minQ = min(tp_terms.sig.TP$q_value[tp_terms.sig.TP$gene_short_name == theGene])
		cellCount = fData(cds_subset.TP[rowData(cds_subset.TP)$gene_short_name == theGene, ])$num_cells_expressed
		geneQcount.TP = rbind(geneQcount.TP, c(theGene, minQ, cellCount))
	}
	geneQcount.TP = as.data.frame(geneQcount.TP)
	colnames(geneQcount.TP) = c("gene_short_name", "q_value", "num_cells_expressed")
	dim(geneQcount.TP)
	# 2655    3
	
	outputFile =  file.path(thisOutputDir, 
							paste(projectName, 
								  ".monocle3results.DEgenes.TP.pval0.01.tsv", sep=""))
	write.table(geneQcount.TP, file=outputFile, sep="\t", row.names=FALSE, col.names=TRUE)
	system(paste('gzip -f ', outputFile, '\n', sep=""), intern=TRUE)
	
	
	
	##################################################
	# Tabulate expression summaries
	# Copied and adapted from 'analyzeSciRNA.monocle2.R'
	
	cat("\n\n# ***********\nTabulate expression summaries\n", sep="")
	
	# 20200120
	# https://rdrr.io/github/cole-trapnell-lab/monocle3/man/exprs-cell_data_set-method.html
	# CDSmatrix = as.matrix( Biobase::exprs(cds) )
	CDSmatrix = as.matrix( exprs(cds) )
	
	print( dim(CDSmatrix) )
	# # [1] 52636  7769
	# 52073  9817
	
	pData(cds)$sex # factor( sampleSex, levels=c("female", "male") )
	pData(cds)$timepoint #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	
	print( length(pData(cds)$sex) ) # factor( sampleSex, levels=c("female", "male") )
	# 9817
	
	print( length(pData(cds)$timepoint) ) #factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
	# 9817
	
	colnames(CDSmatrix) = paste(colnames(CDSmatrix), ".", pData(cds)$timepoint, ".", pData(cds)$sex, sep="")
	rownames(CDSmatrix) = paste("chr", fData(cds)$ensChr, ":", fData(cds)$ensStart, "-", fData(cds)$ensEnd, "|", fData(cds)$ensStrand, "|", fData(cds)$id, "|", fData(cds)$gene_short_name, sep="")
	
	meanMatrix = Matrix::rowMeans(CDSmatrix, na.rm=TRUE)   
	matrixfilename = file.path(thisOutputDir, paste(projectName, ".UMIsPerGene.mean.", theCDSdescS, ".tsv", sep=""))
	# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	fwrite(as.data.frame(meanMatrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
	
	medianMatrix = apply(X=CDSmatrix, MARGIN=1, FUN=median, na.rm=TRUE)   
	matrixfilename = file.path(thisOutputDir, paste(projectName, ".UMIsPerGene.median.", theCDSdescS, ".tsv", sep=""))
	# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	fwrite(as.data.frame(medianMatrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
	
	Q75Matrix = apply(X=CDSmatrix, MARGIN=1, FUN=quantile, probs=0.75, na.rm=TRUE)   
	matrixfilename = file.path(thisOutputDir, paste(projectName, ".UMIsPerGene.Q75.", theCDSdescS, ".tsv", sep=""))
	# write.table(CDSmatrix, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	fwrite(as.data.frame(Q75Matrix), matrixfilename, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
	try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))
	
	# table(pData(cds)$timepoint, pData(cds)$sex)
	
	
	
	##################################################
	# Find marker genes expressed by each timepoint (not by each cluster)
	# i.e. Differential analysis
	
	cat("\n\n# ***********\nFind marker genes\n", sep="")
	
	# The data frame marker_test_res contains a number of metrics for how specifically expressed each gene is in each group (time point or partition). 
	# We could group the cells according to cluster, partition, or any categorical variable in colData(cds). 
	marker_test_res = top_markers(cds, group_cells_by="timepoint", reference_cells=NULL, cores=1)
	head(marker_test_res)
	
	# colnames(marker_test_res)
	# cell_group fraction_expressing             gene_id     gene_short_name 
	# marker_score marker_test_p_value marker_test_q_value     mean_expression pseudo_R2         specificity 
	
	# You can rank the table according to one or more of the specificity metrics 
	# and take the top gene for each cluster. 
	# For example, pseudo_R2 is one such measure. We can rank markers according to pseudo_R2 like this:
	
	n_genes = 5
	
	# top_specific_markers = marker_test_res %>%
	#     filter(fraction_expressing >= 0.10) %>%
	#     group_by(cell_group) %>%
	#     top_n(n_genes, pseudo_R2) %>%
	# arrange(cell_group)
	
	top_specific_markers = marker_test_res %>%
		filter(fraction_expressing >= 0.10) %>%
		group_by(cell_group) %>%
		top_n(n_genes, specificity) %>%
		arrange(specificity) %>%
		arrange(cell_group)
	
	dim(top_specific_markers)
	top_specific_markers[1:(n_genes*6), c("gene_short_name", "cell_group","specificity")]
	
	# # Now, we can plot the expression and fraction of cells
	# # that express each marker in each group with the plot_genes_by_group function:
	# PDFfilename = file.path(thisOutputDir,
	#                         paste(projectName, ".monocle3results.top",
	#                               n_genes, "markerGenes.allTimepoints.pdf", sep=""))
	# pdf(PDFfilename, height=15, width=10)
	# 
	# # top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
	# top_specific_marker_ids = unique(top_specific_markers %>%
	#                                      pull(gene_short_name))
	# plot_genes_by_group(cds,
	#                     markers=top_specific_marker_ids,
	#                     group_cells_by="timepoint",
	#                     ordering_type="cluster_row_col",# ordering_type="maximal_on_diag", # "none",
	#                     max.size=3) +
	#     ggtitle(paste("Top ", n_genes, " marker genes", sep=""))
	# dev.off()
	
	for (theTP in names(tpCols)) {
		PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle3results.top", 
													 n_genes, "markerGenes.", theTP, ".pdf", sep=""))
		pdf(PDFfilename, height=5, width=10)
		
		top_specific_marker_ids = unique(top_specific_markers %>% 
											 filter(cell_group == theTP) %>% 
											 pull(gene_short_name))
		
		print( plot_genes_by_group(cds,
								   markers=top_specific_marker_ids,
								   group_cells_by="timepoint",
								   ordering_type="none",
								   max.size=5) + 
				   ggtitle(paste("Top ", n_genes, " ", theTP, " marker genes", sep="")) )
		
		# # Color clusters by marker genes
		# print( plot_cells(cds,
		#                   genes=ciliated_genes,
		#                   label_groups_by_cluster=FALSE,
		#                   label_leaves=FALSE,
		#                   label_branch_points=FALSE,
		#                   label_cell_groups=FALSE,
		#                   show_trajectory_graph=FALSE) +
		#            geom_point(size=0.5))
		
		dev.off()
		
	}  
	
	
	
	##################################################
	# Save this work
	
	save.image( file =  file.path( thisOutputDir, paste( ".RData.", theCDSdesc, ".postDE", sep="" ) ) )
	# load( file =  file.path( thisOutputDir,  paste( ".RData.", theCDSdesc, ".postDE", sep="" ) ) )
