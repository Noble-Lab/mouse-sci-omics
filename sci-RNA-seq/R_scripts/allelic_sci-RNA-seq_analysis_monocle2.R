##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200322 allelic sci-RNA-seq analysis using monocle 2 
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
##################################################
##################################################
# Input

projectName = "allelic_sci-RNA-seq"

inputDir = "../data/allelic_sci-RNA-seq_count_matrices/"

superMainOutputDir = paste("../results/allelic_sci-RNA-seq_monocle2/")
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

libs=c(
    "d0",
    "d3-11andNPCs"
)

CDSpaths = list()
CDSpaths[["d0"]] = file.path(inputDir, "d0", "final-output")
CDSpaths[["d3-11andNPCs"]] = file.path(inputDir, "d3-11andNPCs", "final-output")

theChroms = paste( "chr", c( c(1:19), "X"), sep="" )










##################################################
##################################################
##################################################
# Generate CDS data

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
	pData(cds)$n.umi = apply(Biobase::exprs(cds), 2, sum)
	
	return(cds)
}


CDSperLib = list()

for (theLib in names(CDSpaths) ) {
	print(theLib)
	CDSperLib[[theLib]] = load.cds(mat.path=file.path(CDSpaths[[theLib]], "UMI.count.matrix"), 
								   gene.annotation.path=file.path(CDSpaths[[theLib]], "gene.annotations"), 
								   cell.annotation.path=file.path(CDSpaths[[theLib]], "cell.annotations"))
}

head(pData(CDSperLib[[theLib]]))
# # cell   sample Size_Factor n.umi
# # A01_A01_RT_1_alt   A01_A01_RT_1_alt RowBcol3          NA   130
# # A01_A01_RT_1_ref   A01_A01_RT_1_ref RowBcol3          NA   148
# # A01_A01_RT_27_alt A01_A01_RT_27_alt RowBcol3          NA   528
# # A01_A01_RT_27_ref A01_A01_RT_27_ref RowBcol3          NA   605
# # A01_A01_RT_28_alt A01_A01_RT_28_alt RowBcol3          NA   107
# # A01_A01_RT_28_ref A01_A01_RT_28_ref RowBcol3          NA   161
# 
# cell       sample Size_Factor n.umi
# A01_D05_RT_458_alt A01_D05_RT_458_alt d3-11andNPCs          NA   106
# A01_D05_RT_458_ref A01_D05_RT_458_ref d3-11andNPCs          NA   111
# A01_D05_RT_484_alt A01_D05_RT_484_alt d3-11andNPCs          NA   107
# A01_D05_RT_484_ref A01_D05_RT_484_ref d3-11andNPCs          NA   104
# A01_D05_RT_500_alt A01_D05_RT_500_alt d3-11andNPCs          NA   111
# A01_D05_RT_500_ref A01_D05_RT_500_ref d3-11andNPCs          NA   124

totHaploCells = 0
for (theLib in names(CDSpaths) ) {
	cat( theLib, "\t", nrow(pData(CDSperLib[[theLib]])), "\n", sep="" )
	totalExpr <- Matrix::colSums(exprs(CDSperLib[[theLib]]))  
	print(sum(totalExpr))
	print(summary(totalExpr))
	totHaploCells = totHaploCells + nrow(pData(CDSperLib[[theLib]]))
}

# After 20200323 pipeline 'split -a 3' fix:
# 
# d0	3,564
# [1] 866383
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100.0   143.0   196.0   243.1   286.0  1667.0 
# 
# d3-11andNPCs	5,526
# [1] 1540970
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100.0   147.0   207.0   278.9   328.0  4086.0 

totHaploCells
# 3,564 + 5,526 = 9,090


# ##################################################
# # Load CDS data
# 
# CDSperLib = list()
# 
# for (theLib in names(CDSpaths) ) {
#     print(theLib)
#     CDSperLib[[theLib]] = readRDS(file=CDSpaths[[theLib]])
# }





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

dim(pData(cds.OG))/c(2,1)
# 4545    4

dim(fData(cds.OG))
# [1] 52636  1


dim(Biobase::exprs(cds.OG))/c(1,2)
# [1] 52636  4545

	








##################################################
##################################################
##################################################
# Fix fData ids and short names

# head(fData(cds.OG))
# # gene_short_name
# # ENSMUSG00000000001.4   ENSMUSG00000000001.4
# # ENSMUSG00000000003.15 ENSMUSG00000000003.15
# # ENSMUSG00000000028.14 ENSMUSG00000000028.14
# # ENSMUSG00000000031.15 ENSMUSG00000000031.15
# # ENSMUSG00000000037.16 ENSMUSG00000000037.16
# # ENSMUSG00000000049.11 ENSMUSG00000000049.11
# 
# head(fData(CDSperLib[[theLibA]]))
# # id gene_short_name
# # ENSMUSG00000000001.4   ENSMUSG00000000001.4           Gnai3
# # ENSMUSG00000000003.15 ENSMUSG00000000003.15            Pbsn
# # ENSMUSG00000000028.14 ENSMUSG00000000028.14           Cdc45
# # ENSMUSG00000000031.15 ENSMUSG00000000031.15             H19
# # ENSMUSG00000000037.16 ENSMUSG00000000037.16           Scml2
# # ENSMUSG00000000049.11 ENSMUSG00000000049.11            Apoh

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

getRTplateNumbers <- function (cellIDs, plate2IDs) {
	RTids = as.numeric(sub(".*_", "",  sub("_...$", "", cellIDs)))
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
	RTids = as.numeric(sub(".*_", "",  sub("_...$", "", cellIDs)))
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
# # [1] "RowBcol2" "RowBcol3" "He_5"     "He_6"     "He_7" 
# 
# "d0"           "d3-11andNPCs"

table(plate2IDs)/2
# d0 d3-11andNPCs 
# 1782         2763 

cellIDs = pData(cds.OG)$cell
cellCount = length(cellIDs)/2
# 4545

table(cellIDs, plate2IDs)
colSums(table(cellIDs, plate2IDs))/2
# d0 d3-11andNPCs 
# 1782         2763 

cellAlleles = sub(".*_", "", rownames(pData(cds.OG)))

plate2cellIDs = rownames(pData(cds.OG))
plate2RTids = sub("_.*_", "_", sub("_...$", "", plate2cellIDs))
unique(plate2RTids)

RTids = as.numeric(sub(".*_", "", plate2RTids))
unique(RTids)
range(unique(RTids))
# [1]   0 767

summary(RTids)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0    28.0   476.0   367.4   642.0   767.0 

summary(RTids[RTids<384])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00   10.00   23.00   22.85   34.00   47.00 

(767-384+1)/96
# [1] 4

# cellIDs[ pData(cds.OG)$sample =="He_5" ]

RTplateNumbers = getRTplateNumbers(cellIDs, plate2IDs)
table(RTplateNumbers)/2
# 1   2   3   4   5   6 
# 905 877 506 797 548 912 

RTwellNumbers = getRTwellNumbers(cellIDs, plate2IDs)
table(RTwellNumbers)/2
# 0   1   2   3   4   5   6   7   8   9  10  11 
# 280 352 262 402 471 464 444 448 380 438 288 316 

sampleSex = getSampleSex(cellIDs, plate2IDs)
table(sampleSex) /2
# female   male 
# 2208   2337 

sampleTimePoints = getSampleTimePoints(cellIDs, plate2IDs)
table(sampleTimePoints) /2
# d0  d11   d3   d7  NPC 
# 1782  722  469  585  987 

cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="d0",][1:4,]
# cellIDs             plate2IDs RTplateNumbers
# [1,] "A01_B02_RT_16_alt" "d0"      "1"           
# [2,] "A01_B02_RT_16_ref" "d0"      "1"           
# [3,] "A01_B02_RT_17_alt" "d0"      "1"           
# [4,] "A01_B02_RT_17_ref" "d0"      "1"  

table(sampleTimePoints, sampleSex, cellAlleles) 
# , , cellAlleles = alt
# 
# sampleSex
# sampleTimePoints female male
# d0     905  877
# d11    333  389
# d3     320  149
# d7     186  399
# NPC    464  523
# 
# , , cellAlleles = ref
# 
# sampleSex
# sampleTimePoints female male
# d0     905  877
# d11    333  389
# d3     320  149
# d7     186  399
# NPC    464  523

cbind(cellIDs, plate2IDs, RTplateNumbers)[sampleTimePoints=="NA",]
# NADA!

# pData(cds.OG)$batchNum <- factor( RTplateBatches, levels=c("1", "2") )
pData(cds.OG)$sex <- factor( sampleSex, levels=c("female", "male") )
pData(cds.OG)$timepoint <- factor( sampleTimePoints, levels=c("d0", "d3", "d7", "d11", "NPC") )
cellAlleles = sub(".*_", "", rownames(pData(cds.OG)))
cellAlleles[cellAlleles=="alt"]="Cast"
cellAlleles[cellAlleles=="ref"]="129"
pData(cds.OG)$allele <- factor(cellAlleles, levels=c("129", "Cast"))

# cell   sample Size_Factor n.umi
# RowBcol2_A01_A01_RT_0_alt   A01_A01_RT_0_alt RowBcol2          NA   116
# RowBcol2_A01_A01_RT_0_ref   A01_A01_RT_0_ref RowBcol2          NA   108
# RowBcol2_A01_A01_RT_16_alt A01_A01_RT_16_alt RowBcol2          NA   159
# RowBcol2_A01_A01_RT_16_ref A01_A01_RT_16_ref RowBcol2          NA   146
# RowBcol2_A01_A01_RT_17_alt A01_A01_RT_17_alt RowBcol2          NA   246
# RowBcol2_A01_A01_RT_17_ref A01_A01_RT_17_ref RowBcol2          NA   264
# batchNum sex timepoint allele
# RowBcol2_A01_A01_RT_0_alt         1 female        d0   Cast
# RowBcol2_A01_A01_RT_0_ref         1 female        d0    129
# RowBcol2_A01_A01_RT_16_alt        1 female        d0   Cast
# RowBcol2_A01_A01_RT_16_ref        1 female        d0    129
# RowBcol2_A01_A01_RT_17_alt        1 female        d0   Cast
# RowBcol2_A01_A01_RT_17_ref        1 female        d0    129










##################################################
##################################################
##################################################
# Add gencode chromosomal annotations

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
# Continue exploring data

# Total expression per barcode
pData(cds.OG)$totalExpr <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!
pData(cds.OG)[1:5,]

# Total expression distributions 
summary(pData(cds.OG)$totalExpr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 100.0   145.0   202.0   264.8   311.0  4086.0 

summary(log10(pData(cds.OG)$totalExpr))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.000   2.161   2.305   2.348   2.493   3.611 

# Total expression per secondary plate
for (plate2ID in unique(plate2IDs) ) {
	cat(plate2ID, "\t", sum(pData(cds.OG)$totalExpr[plate2IDs == plate2ID], na.rm=TRUE), "\n", sep="") 
}
# d0	866383
# d3-11andNPCs	1540970

# Total expression per RT (primaty) plate
for (plate1ID in unique(RTplateNumbers) ) {
	cat(plate1ID, "\t", sum(pData(cds.OG)$totalExpr[RTplateNumbers == plate1ID], na.rm=TRUE), "\n", sep="") 
}
# 1	433470
# 2	432913
# 3	220301
# 4	503767
# 5	218710
# 6	598192

for (timePoint in unique(pData(cds.OG)$timepoint) ) {
	cat(timePoint, "\t", sum(pData(cds.OG)$totalExpr[ pData(cds.OG)$timepoint == timePoint ], na.rm=TRUE), "\n", sep="") 
}
# d0	866383
# d3	196032
# d7	242979
# d11	350967
# NPC	750992










##################################################
##################################################
##################################################
# 20200516 
# Filter for low counts

totalCoverage = pData(cds.OG)$totalExpr 

length(totalCoverage)
# 9090

# Find biased Genes adjusted for bias
totCovDensity <- density(log10(totalCoverage))
maxTCdensity = max(totCovDensity$y)

maxTotCovDensity = totCovDensity$x[totCovDensity$y == maxTCdensity]
print(c(maxTotCovDensity, 10^maxTotCovDensity))

maxTotCovDensity = 160 # *** NB Based on distribution minima **** 

PDFfilename = file.path(superMainOutputDir, paste(projectName, ".monocle2results.coverageDistribution.unfiltered.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
plot(totCovDensity, xlab="log10(total UMIs)")
abline(v=log10(maxTotCovDensity), lty="dashed", col="red")
text(log10(maxTotCovDensity), 0, maxTotCovDensity, pos=4)
dev.off()

# Cells to exclude due to extreme allelic bias
altLowCellsBoolIdx = totalCoverage > maxTotCovDensity
sum(altLowCellsBoolIdx)

cds.OG = cds.OG[,altLowCellsBoolIdx]

dim(cds.OG)
# Features  Samples 
# 52636     6122 

totalCoverage = pData(cds.OG)$totalExpr 

length(totalCoverage)
# 6122

summary(totalCoverage)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 161.0   201.0   263.0   330.1   380.0  4086.0 

# Remove singlets
# i.e. cells with reads from only one parental allele
i = 1
singletList = c()
while( i < nrow(pData(cds.OG)) ) {
	if(sub("_...$", "", pData(cds.OG)$cell[i]) != sub("_...$", "", pData(cds.OG)$cell[i+1])) {
		singletList = c(singletList, i)
		i = i + 1
	} else {
		i = i + 2
	}
}
length(singletList)
# 438

pData(cds.OG)$cell[singletList]

cds.OG <- cds.OG[,-singletList]

# Check alt == ref
refBoolIdx = grepl("_ref$", rownames(pData(cds.OG)))
altBoolIdx = grepl("_alt$", rownames(pData(cds.OG)))
if( sum(altBoolIdx) != sum(refBoolIdx) ) {
	cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
}

dim(cds.OG)
# Features  Samples 
# 52636     5684 
	









##################################################
##################################################
##################################################
# Normalize and pre-process the data

cds.OG <- estimateSizeFactors(cds.OG)
pData(cds.OG)[1:5,]

# # 20190723
# 20200323
# Unfiltered by expressed genes
save.image(file=file.path(superMainOutputDir,".RData")) # To .RData
# load(file=file.path(superMainOutputDir,".RData"))

dim(cds.OG)
# Features  Samples 
# 52636     5684 

# Note: many more d0 cells
table(cds.OG$timepoint, cds.OG$sex)
# female male
# d0    1098 1064
# d3     290  148
# d7     216  362
# d11    436  420
# NPC    760  890



##################################################
# Distribution of cell counts with expression

min_expr = 0.1  # == 0 i.e. >= 1 # see below
# min_expr = 1.1 # i.e. >= 2
# min_expr = 2.1 # i.e. >= 3

cds.OG <- detectGenes(cds.OG, min_expr = min_expr)



##################################################
# Distribution of cell counts gene is expressed in:
summary(fData(cds.OG)$num_cells_expressed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0    34.7     5.0  6781.0 
#
# i.e. Only 25% of genes are expressed in 5 or more cells.

# min_expr = 0 == 0.1
summary(fData(cds.OG)$num_cells_expressed) 
summary( Matrix::rowSums(Biobase::exprs(cds.OG) > 0) )
summary( Matrix::rowSums(Biobase::exprs(cds.OG) > 0.1) )
summary( Matrix::rowSums(Biobase::exprs(cds.OG) >= 1) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     0.0    34.7     5.0  6781.0 

summary( Matrix::rowSums(Biobase::exprs(cds.OG) > 1) )
summary( Matrix::rowSums(Biobase::exprs(cds.OG) >=2) )
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.000    0.000    6.865    0.000 5138.000 

quantile(p=seq(0.75, 0.95, 0.05), fData(cds.OG)$num_cells_expressed)
# 75% 80% 85% 90% 95% 
# 5  13  31  72 190 

cellCounts=c()
for (i in c(2, 5, 10, 20, 25, 50, 100)) {
	cellCounts[as.character(i)] = sprintf( "%s%% of genes expressed in %s+ (%s%%) haplocells",
										   round( sum(fData(cds.OG)$num_cells_expressed > i) / length(fData(cds.OG)$num_cells_expressed) * 100, 2 ), 
										   i, round(i / nrow(pData(cds.OG)) * 100, 2) )
}
cellCounts
# 29.25% of genes expressed in 2+ (0.02%) haplocells
# 24.85% of genes expressed in 5+ (0.06%) haplocells
# 21.22% of genes expressed in 10+ (0.11%) haplocells
# 17.4% of genes expressed in 20+ (0.22%) haplocells
# 16.14% of genes expressed in 25+ (0.28%) haplocells
# 12.03% of genes expressed in 50+ (0.55%) haplocells
# 8.18% of genes expressed in 100+ (1.1%) haplocells



##################################################
# Distribution of proportions of cells gene is expressed in:
summary(fData(cds.OG)$num_cells_expressed) / nrow(pData(cds.OG))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000000 0.004004 0.000577 0.782483 

meanPropotion = summary(fData(cds.OG)$num_cells_expressed)["Mean"]  / nrow(pData(cds.OG))

# Expressed genes in top half of cells:
expressedGenes <-  row.names(subset(fData(cds.OG), 
									num_cells_expressed >= ( meanPropotion * nrow(pData(cds.OG)) ) ) )
length(expressedGenes)
# 7509

# top 5% most expressed genes in terms of cells:
expressedGenes <-  row.names(subset(fData(cds.OG), 
									num_cells_expressed >= ( 0.05 * nrow(pData(cds.OG)) ) ) )
length(expressedGenes)
# 19

# 10 cells
expressed_genes <- row.names(subset(fData(cds.OG),
									num_cells_expressed >= 10))



##################################################
# Distributions of num_cells_expressed

PDFfilename = file.path(superMainOutputDir, paste(projectName, ".monocle2results.numCellsExpressedDistributions.pdf", sep=""))
pdf(PDFfilename, height=4*2, width=4*4)
par(mfrow=c(2, 4))

plot( density(log10(fData(cds.OG)$num_cells_expressed + 1)), 
	  main="log10 num_cells_expressed + 1", xlab="num cells", ylab="gene count density" )
legend("topright", legend=cellCounts, cex=1.0, bty='n')
legend("right", legend=paste0("No cells: ", sum(fData(cds.OG)$num_cells_expressed<=0), " / ",
							  length(fData(cds.OG)$num_cells_expressed), 
							  " (", round( sum(fData(cds.OG)$num_cells_expressed<=0) / 
											   length(fData(cds.OG)$num_cells_expressed) * 100, 2), 
							  "% of genes)"), bty='n')
plot( density(log10(fData(cds.OG)$num_cells_expressed)), 
	  main="log10 num_cells_expressed", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds.OG)$num_cells_expressed[fData(cds.OG)$num_cells_expressed>0])), 
	  main="log10 num_cells_expressed > 0", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds.OG)$num_cells_expressed[fData(cds.OG)$num_cells_expressed>0] + 1)), 
	  main="log10 num_cells_expressed > 0 + 1", xlab="num cells", ylab="gene count density" )

plot( density(log10(fData(cds.OG[expressed_genes,])$num_cells_expressed + 1)), 
	  main="log10 num_cells_expressed + 1 for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds.OG[expressed_genes,])$num_cells_expressed)), 
	  main="log10 num_cells_expressed for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds.OG[expressed_genes,])$num_cells_expressed[fData(cds.OG[expressed_genes,])$num_cells_expressed>0])), 
	  main="log10 num_cells_expressed > 0 for expressed genes", xlab="num cells", ylab="gene count density" )
plot( density(log10(fData(cds.OG[expressed_genes,])$num_cells_expressed[fData(cds.OG[expressed_genes,])$num_cells_expressed>0] + 1)), 
	  main="log10 num_cells_expressed > 0 + 1 for expressed genes", xlab="num cells", ylab="gene count density" )

dev.off()


    







##################################################
##################################################
##################################################
# Filter for expressed genes

min_expr = 0.1  # == 0 i.e. >= 1 # see below
# min_expr = 1.1 # i.e. >= 2
# min_expr = 2.1 # i.e. >= 3

min_cells = 5
# min_cells = 10

mainOutputDir = file.path(superMainOutputDir, paste0("minExpr", min_expr, "_minCells", min_cells) )
if ( ! file.exists( mainOutputDir ) ) {
    if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

# save.image(file=file.path(mainOutputDir,".RData")) # To .RData
# load(file=file.path(mainOutputDir,".RData"))



##################################################
# Perform filtering for expressed genes

# how many express a particular gene, or how many genes are expressed by a given cell?
cds.OG <- detectGenes(cds.OG, min_expr = min_expr)
# function (cds, min_expr = NULL) 
# {
#     if (is.null(min_expr)) {
#         min_expr <- cds@lowerDetectionLimit
#     }
#     fData(cds)$num_cells_expressed <- Matrix::rowSums(Biobase::exprs(cds) > 
#                                                           min_expr)
#     pData(cds)$num_genes_expressed <- Matrix::colSums(Biobase::exprs(cds) > 
#                                                           min_expr)
#     cds
# }
print(head(fData(cds.OG)))
    
# The vector expressed_genes now holds the identifiers for genes expressed in at least n cells of the data set.
expressed_genes <- row.names(subset(fData(cds.OG),
									num_cells_expressed >= min_cells))
length(expressed_genes)
# 11,461 # 0.1 == 0
# 4,707 # 1

print(head(pData(cds.OG)))

cds.OG <- cds.OG[expressed_genes,]
print(dim(cds.OG))

# 0.1 in >= 5 cells
# Features  Samples 
# 12969     5684 

    
    







##################################################
##################################################
##################################################
# Filter out extreme cells

# rm(cds)
upper_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) +
					   2*sd(log10(pData(cds.OG)$totalExpr)))
lower_bound <- 10^(mean(log10(pData(cds.OG)$totalExpr)) -
					   2*sd(log10(pData(cds.OG)$totalExpr)))
cds <- cds.OG[,pData(cds.OG)$totalExpr > lower_bound &
				  pData(cds.OG)$totalExpr < upper_bound]

# Remove singlets
# i.e. cells with reads from only one parental allele
i = 1
singletList = c()
while( i < nrow(pData(cds)) ) {
	if(sub("_...$", "", pData(cds)$cell[i]) != sub("_...$", "", pData(cds)$cell[i+1])) {
		singletList = c(singletList, i)
		i = i + 1
	} else {
		i = i + 2
	}
}
length(singletList)
# 36

pData(cds)$cell[singletList]

cds <- cds[,-singletList]

# Check alt == ref
refBoolIdx = grepl("_ref$", rownames(pData(cds)))
altBoolIdx = grepl("_alt$", rownames(pData(cds)))
if( sum(altBoolIdx) != sum(refBoolIdx) ) {
	cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
}

dim(cds.OG)
# Features  Samples 
# 12969     5684 


dim(cds)
# Features  Samples 
# 12969     5392 

head(pData(cds))

# pData(cds)$sex <- factor( pData(cds)$sex, levels=c("female", "male") )
# pData(cds)$timepoint <- factor( pData(cds)$timepoint, levels=c("d0", "d3", "d7", "d11", "NPC") )
# pData(cds)$batchNum <- factor( pData(cds)$batchNum, levels=c("1", "2") )

# rm(CDSperLib)

fData(cds)[ grepl("Xist", fData(cds)[,"gene_short_name"]), ]
# gene_short_name                   id ensChr  ensStart ensEnd ensStrand num_cells_expressed
# ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3      X 103460366  103483254         -                 424










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

sub_min_expr = 0.1; sub_min_cells = 5;
# sub_min_expr = 0.1; sub_min_cells = 10;
# sub_min_expr = 1.1; sub_min_cells = 10; 
# sub_min_expr = 1.1; sub_min_cells = 5; 

whichChroms2analze = list()
whichChroms2analze[["chr1andchrX"]] = c(1,20)
whichChroms2analze[["genome"]] = 1:20

whichChroms = "chr1andchrX"
# whichChroms = "genome"

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

thisOutputDir = file.path( mainOutputDir, "differential_allelic_expression_distributions", 
                           paste0( whichChroms, 
                                   ifelse(useMeanThreshold, "", ".chromosomewideFPR"), 
                                   ifelse(min_counts_per_chromosome==0, "", 
                                          paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                   ifelse(rmXist, ".rmXist", "")
                                   ) )
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
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
YlimsPerChrom.log2[["chr1"]] = c(0, 0.20)
YlimsPerChrom.log2[["chrX"]] = c(0, 0.15)
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
        
        # Exclude Xist
        if (rmXist) {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
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
    # Features  Samples 
    # 925     5392 
    # chr1: 2696 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/ Xist:
    # Features  Samples 
    # 390     5392 
    # chrX: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # Features  Samples 
    # 389     5392 
    # chrX: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)

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
    # w/ Xist:
    # All chromosomes: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # All chromosomes: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    
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
    # w/ Xist:
    # [1] "female" "d0"    
    # d0: 489 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 78 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 145 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 209 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1048 / 1344 female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 486 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 76 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 143 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 187 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1019 / 1344female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # [1] "male" "d0"  
    # d0: 295 / 503 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 57 / 73 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 125 / 181 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 152 / 206 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 225 / 389 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 854 / 1352 male cells pass filter (min_counts_per_chromosome >= 10)
    
}

# --------------------------------------------------
    
############
# For each chromosome 
for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
    # for (theChrom in theChroms) {
    print(theChrom)
    
    if (plotOldPlots) {
        
        daScatterPlots = list()
        
        daPlots = list()
        daPlots.scaled = list()
        daPlots.hist = list()
        daPlots.hist.scaled = list()
        
        daLogPlots = list()
        daLogPlots.scaled = list()
        daLogPlots.hist = list()
        daLogPlots.hist.scaled = list()
        
    }
    
    daLogPlots.hist.scaled.na = list()
    daLogPlots.hist.scaled.na.new = list()
    
    plotally = 0
    
    # --------------------------------------------------
    # 20200709
    # Exclude Xist
    
    if (rmXist) {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
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
    autosomal.log2Ratio.FPRthreshold.perChrom[[theChrom]] =  quantile(absLog2AllelicExp, (100-FPRthresh)/100, na.rm=TRUE)

    for (aSex in c("female", "male")) {
        # aSex = "female"
        
        if(theChrom =="chrX" && aSex == "male") {
            next
        }
        
        subsexBoolidx = pData(cdsC.OG)$sex == aSex
        
        #for (aTP in c("all", names(tpCols))) {
        for (aTP in c(names(tpCols))) {
            # aTP = "d0"
            
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

            # length(refBoolIdx)
            # [1] 688
            
            # --------------------------------------------------
            # 20200709
            # Apply thresholds for minimum counts per chromosome 
            # For each chromosome 
            
            if ( min_counts_per_chromosome > 0) {
            
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


                 
            ##########            
            # Scatter plots
            
            df <- as.data.frame(
                cbind( log10(pData(cdsC.TP.Sex)$totalExpr[refBoolIdx] + 0.1),
                       log10(pData(cdsC.TP.Sex)$totalExpr[altBoolIdx] + 0.1) ) )
            colnames(df) = c("allele_129", "allele_Cast")
            
            daScatterPlot <- ggplot(df, aes(x=allele_129, y=allele_Cast)) +
                geom_point(size=1) +
                ggtitle( paste(theChrom, ": ", aSex, " filtered cells (n = ", dim(cdsC.TP.Sex[,refBoolIdx])[2], "; ",
                               " genes = ", dim(cdsC.TP.Sex[,refBoolIdx])[1], ")", sep="") ) +
                ylab("129 UMIs (log10)") +
                xlab("Cast UMIs (log10)") +
                coord_cartesian(xlim = c(0, log10(maxChromCount+0.1)), ylim = c(0, log10(maxChromCount+0.1)))  +
                geom_abline(intercept = 0, slope = 1, colour = 'red', linetype="dashed") 
            
            daScatterPlots[[paste0("p", as.character(plotally))]] = as.grob(daScatterPlot)
            
            
            
            #####################################
            # Absolute differential allelic expression
            
            diffAllelicExp = refExprC - altExprC
            absDiffAllelicExp = abs(diffAllelicExp)
            
            # Save absolute differential scores
            # When chrX: obtain mean value and select XCI cells
            if (theChrom != "chrX") {
                autosomal.differential.FPRthreshold.vector = c( autosomal.differential.FPRthreshold.vector,
                                                                quantile(absDiffAllelicExp, (100-FPRthresh)/100, na.rm=TRUE) )
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
            ## Absolute log2 ratio of allelic expression
            
            log2AllelicExp = log2(refExprC + 0.1) - log2(altExprC + 0.1)
            absLog2AllelicExp = abs(log2AllelicExp)

            # Save absolute differential scores
            # When chrX: obtain mean value and select XCI cells
            if (theChrom != "chrX") {
                autosomal.log2Ratio.FPRthreshold.vector = c( autosomal.log2Ratio.FPRthreshold.vector,
                                                         quantile(absLog2AllelicExp, 0.9, na.rm=TRUE) )
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
            
            theXlims = c(-maxX.l2r, maxX.l2r)
            
            daLogPlot.hist.scaled.na <- ggplot(pData(cdsC.TP.Sex)[refBoolIdx,], aes(x=LAexpr)) +
                # geom_histogram(aes(y=..count../sum(..count..)), binwidth = 0.25, center=0, alpha = 0.4, size=0.75) +
                geom_histogram(aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]), binwidth = 0.25, center=0, alpha = 0.4, size=0.75) +
                ggtitle( paste(theChrom, ": ", aSex, " filtered cells (n = ", dim(cdsC.TP.Sex[,refBoolIdx])[2], "; ",
                               " genes = ", dim(cdsC.TP.Sex[,refBoolIdx])[1], ")", sep="") ) +
                ylab("Proportion of cells") +
                xlab("log2 ratio of allelic expression") +
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
                # xlab("log2 ratio of allelic expression") +
                theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                      axis.title.y=element_blank(), axis.text.y=element_blank()) +
                coord_cartesian(xlim = theXlims, ylim = theYlims)
            
            if (theChrom == "chrX") {
                daLogPlot.hist.scaled.na.new <- daLogPlot.hist.scaled.na.new +
                    geom_vline(xintercept=-autosomal.log2Ratio.FPRthreshold.mean, linetype = "dotted", colour = "gray", lwd=1.5) +
                    geom_vline(xintercept=autosomal.log2Ratio.FPRthreshold.mean, linetype = "dotted", colour = "gray", lwd=1.5)
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
                                                             top = paste(theChrom, " scaled histograms of log2 differential expression based on total expression per cell across ", sum(subsetBoolidx), " cells for expressed genes" ) )
    
    alldaLogPlots.hist.scaled.na.new[[theChrom]] <- marrangeGrob(grobs = daLogPlots.hist.scaled.na.new, layout_matrix=matrix(1:10, nrow=2, byrow=TRUE),
                                                                 top = paste(theChrom, " scaled histograms of log2 differential expression based on total expression per cell across ", sum(subsetBoolidx), " cells for expressed genes" ) )
 
} # chrom


#########
# Tabulate results

colnames(XCIresultsTable.diff) = names(tpCols)
colnames(XCIresultsTable.l2r) = names(tpCols)
rownames(XCIresultsTable.diff) = c("threshold", "XCI_counts", "total_count", "proportion")
rownames(XCIresultsTable.l2r) = c("threshold", "XCI_counts", "total_count", "proportion")

tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCItally.byDifferentialAllelicTotalExpression.", FPRthresh, "pctFPR.byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                               ".", whichChroms, 
                                               ifelse(stringentFilter, ".stringentlyFiltered", ""), ".tsv", sep=""))
write.table(x = XCIresultsTable.diff, file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCItally.byLogRatioAllelicTotalExpression.", FPRthresh, "pctFPR.byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                               ".", whichChroms, 
                                               ifelse(stringentFilter, ".stringentlyFiltered", ""), ".tsv", sep=""))
write.table(x = XCIresultsTable.l2r, file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


#########
# XCI cell table

for (aTP in names(XCIcellsPerTimepoint[[whichChroms]])) {
    tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCIcells.", aTP, ".byLog2AllelicExpression.", FPRthresh, "pctFPR.byTimePoint", 
                                                   ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                                   ".", whichChroms, 
                                                   ifelse(stringentFilter, ".stringentlyFiltered", ""),  ".tsv", sep=""))
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
#  1344    6

# save
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.XCIcells.all.byLog2AllelicExpression.", FPRthresh, "pctFPR.byTimePoint", 
                                               ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells,
                                               ".", whichChroms, 
                                               ifelse(stringentFilter, ".stringentlyFiltered", ""),  ".tsv", sep=""))
write.table(x = XCIcellsAll[[whichChroms]], file = tablefilename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)


#########
# plot

require(gridExtra)
require(grid)

# Scaled histograms of log2 ratio of  allelic expression per chromosome
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.V2.scaled.histogramsOfLog2AllelicExpressionPerChromosome.faceted",
                                             ".numCellsExpressedPerSubset.minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                             ".", whichChroms, ".pdf", sep=""))
pdf(PDFfilename, height=2*5, width=5*5)
for (aGrob in names(alldaLogPlots.hist.scaled.na)) {
    print(aGrob)
    print(alldaLogPlots.hist.scaled.na[[aGrob]])
}
dev.off()

# Fancy scaled histograms of log2 ratio of  allelic expression per chromosome
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.V2.scaled.fancyHistogramsOfLog2AllelicExpressionPerChromosome.faceted",
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

alleleCDPmatrixDir = file.path("../data", "allelic_sci-RNA-seq_catenated_count_matrices")
if ( ! file.exists( alleleCDPmatrixDir ) ) {
    if ( ! dir.create( alleleCDPmatrixDir, showWarnings=FALSE, recursive=TRUE )  ) alleleCDPmatrixDir = mainOutputDir
}

# F121libs=c(
#     "F121_d0",
#     "F121_d3",
#     "F121_d7",
#     "F121_d11",
#     "F121_NPC"
# )

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
# And potentially exclude Xist

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
        
        # Exclude Xist
        if (rmXist) {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
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
    # Features  Samples 
    # 925     5392 
    # chr1: 2696 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/ Xist:
    # Features  Samples 
    # 390     5392 
    # chrX: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # Features  Samples 
    # 389     5392 
    # chrX: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    
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
    # w/ Xist:
    # All chromosomes: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # All chromosomes: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    
    
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
    # w/ Xist:
    # [1] "female" "d0"    
    # d0: 489 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 78 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 145 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 209 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1048 / 1344 female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 486 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 76 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 143 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 187 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1019 / 1344female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # [1] "male" "d0"  
    # d0: 295 / 503 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 57 / 73 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 125 / 181 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 152 / 206 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 225 / 389 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 854 / 1352 male cells pass filter (min_counts_per_chromosome >= 10)

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
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
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
            
            if ( min_counts_per_chromosome > 0) {
                
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
            
            ##################################################
            # Save filtered CDS to file

            CDSmatrix = as.matrix(Biobase::exprs(cdsC.TP.Sex))
            # 
            # dim(CDSmatrix)
            # # 349 4238
            # 
            # CDSmatrix[1:2, 1:2]
            # # d0_A01_B02_RT_16_alt d0_A01_B02_RT_16_ref
            # # ENSMUSG00000000037.16                    0                    0
            # # ENSMUSG00000000134.17                    0                    0
            # 
            # CDSmatrix[, refBoolIdx][1:2, 1:2]
            # # ENSMUSG00000000037.16                    0                    0
            # # ENSMUSG00000000134.17                    0                    0
            # CDSmatrix[, altBoolIdx][1:2, 1:2]
            # # ENSMUSG00000000037.16                    0                    0
            # # ENSMUSG00000000134.17                    0                    0
 
            # Concatenate UMIs for two alleles
            CDSmatrix.catenated = rbind(CDSmatrix[, refBoolIdx], CDSmatrix[, altBoolIdx])
            rownames(CDSmatrix.catenated) = paste( c( rep( "ref", nrow(CDSmatrix) ),
                                                      rep( "alt", nrow(CDSmatrix) ) ),
                                                   rownames(CDSmatrix.catenated), 
                                                   rep(paste0( fData(cdsC.TP.Sex)$ensChr, ":",
                                                           fData(cdsC.TP.Sex)$ensStart, "-",
                                                           fData(cdsC.TP.Sex)$ensEnd ), 
                                                       # fData(cdsC.TP.Sex)$ensStrand 
                                                       2 ),
                                                   sep="_" )

            dim(CDSmatrix.catenated)
            # [1] 698 2119
            
            length(cdsC.TP.Sex$timepoint[refBoolIdx])
            # [1] 2119
            
            colnames(CDSmatrix.catenated) =  sub("d3-11andNPCs_", "",
                                                 sub( "_ref", "", 
                                                      paste0( "F121_", cdsC.TP.Sex$timepoint[refBoolIdx], ".",
                                                            colnames(CDSmatrix.catenated) ) ) )
            
            
            matrixfilename = file.path(alleleCDPmatrixDir, paste0(projectName, ".UMImatrix.filtered", 
                                                                  ".minExpr", sub_min_expr, ".minCells", sub_min_cells, 
                                                                  ".F121.", theChrom, 
                                                                  ifelse(stringentFilter, ".stringentlyFiltered", ""), 
                                                                  ifelse(min_counts_per_chromosome==0, "", 
                                                                         paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                                                  ifelse(rmXist, ".rmXist", ""),
                                                                  ".tsv"))
            write.table(CDSmatrix.catenated, matrixfilename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
            try(system(paste("gzip -f \"", sub("~/", "$HOME/", matrixfilename), "\"", sep = "")))

            print( table(cdsC.TP.Sex$timepoint[refBoolIdx]) )

            # sciRNA-seq
            # # d0  d3  d7 d11 NPC 
            # # 540 142 104 214 344
            # min_counts_per_chromosome > = 10 (chr1 & chrX only)
            # d0  d3  d7 d11 NPC 
            # 489 127  78 145 209 
             
        } # timepoint
        
    } # sex
    
} # chrom










##################################################
##################################################
##################################################
# 20200615
# Heatmaps and clustering of allelic expression
#
# 20200629
# Only XCI cells
#
# 20200709
# Apply thresholds for minimum counts per chromosome 
# And potentially exclude Xist

# install.packages("tidyHeatmap")
require("ComplexHeatmap")

# 20200709
# # XCIthresh = 1.2239 # 10% FPR thresholds based on chr1 distributions of allelic ratio of total chromosomal expression
XCIthresh = 1.2050998 # 10% FPR thresholds based on chr1 distributions of allelic ratio of total chromosomal expression - after filterint by min_counts_per_chromosome >= 10

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

# --------------------------------------------------
# 20200709
# Apply thresholds for minimum counts per chromosome 
# For each chromosome 
# 
# And potentially exclude Xist

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
        
        # Allow for exclusion of Xist - but don't exclude it.
        if (rmXist) {
            chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
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
    # Features  Samples 
    # 925     5392 
    # chr1: 2696 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/ Xist:
    # Features  Samples 
    # 390     5392 
    # chrX: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # Features  Samples 
    # 389     5392 
    # chrX: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    
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
    # w/ Xist:
    # All chromosomes: 1902 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:
    # All chromosomes: 1870 / 2696 cells pass filter (min_counts_per_chromosome >= 10)
    
    
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
    # w/ Xist:
    # [1] "female" "d0"    
    # d0: 489 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 78 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 145 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 209 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1048 / 1344 female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # w/o Xist:    
    # [1] "female" "d0"    
    # d0: 486 / 540 d0 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d3"    
    # d3: 127 / 142 d3 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d7"    
    # d7: 76 / 104 d7 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "d11"   
    # d11: 143 / 214 d11 female cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "female" "NPC"   
    # NPC: 187 / 344 NPC female cells pass filter (min_counts_per_chromosome >= 10)
    # 1019 / 1344 female cells pass filter (min_counts_per_chromosome >= 10)
    #
    # [1] "male" "d0"  
    # d0: 295 / 503 d0 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d3"  
    # d3: 57 / 73 d3 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d7"  
    # d7: 125 / 181 d7 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "d11" 
    # d11: 152 / 206 d11 male cells pass filter (min_counts_per_chromosome >= 10)
    # [1] "male" "NPC" 
    # NPC: 225 / 389 NPC male cells pass filter (min_counts_per_chromosome >= 10)
    # 854 / 1352 male cells pass filter (min_counts_per_chromosome >= 10)
    
}

# --------------------------------------------------

for (theChrom in theChroms[ whichChroms2analze[[whichChroms]] ]) {
    # theChrom = "chrX"
    # theChrom = "chr1"
    print(theChrom)

    thisOutputDir = file.path( mainOutputDir, 
                               paste0( "heatmaps",
                                       ifelse(min_counts_per_chromosome==0, "", 
                                              paste0(".min_counts_per_chromosome", min_counts_per_chromosome)),
                                       ifelse(rmXist, ".rmXist", "") ),
                               theChrom )
    if ( ! file.exists( thisOutputDir ) ) {
        if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
    }
    
    # --------------------------------------------------
    # 20200709
    # Allow for exclusion of Xist - but don't.
    
    if (rmXist) {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom & fData(cds)$gene_short_name != "Xist" ) ]
    } else {
        chromgenesNames <- row.names(fData(cds))[ which( fData(cds)$ensChr == theChrom ) ]
    }
    cdsC.OG <- cds[chromgenesNames,]
    print(dim(cdsC.OG))
    # Features  Samples 
    # 390     5392 
    
    # --------------------------------------------------

    # --------------------------------------------------
    # 20200709
    # Apply thresholds for minimum counts per chromosome 
    # For each chromosome 
    
    if ( min_counts_per_chromosome > 0) {
        cdsC.OG <- cdsC.OG[,rep(min_exprC_boolidx, each=2)]
        print(dim(cdsC.OG))
        # Features  Samples 
        # 390     3804 
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
    # # 1344 cells 
    # 1019
    
    # XCIcellsAll[[whichChroms]][1:3,]
    # # refExpr altExpr   log2Expr auto10FPRthreshold XCI
    # # d0_A01_B02_RT_17       5       9 -0.8353693           1.223988   0
    # # d0_A01_B03_RT_6       14       6  1.2088140           1.223988   1
    # # d0_A02_B02_RT_17       0       1 -3.4594316           1.223988   1

    theseCellIDs = sub("_ref", "", rownames(pData(cdsC.OG))[refBoolIdx])
    length(theseCellIDs)
    # # [1] 2696
    # [1] 1902
    
    cellsWithXCIstatus = match( theseCellIDs, rownames(XCIcellsAll[[whichChroms]]) )
    # all( rownames(theseCellIDs)[!is.na(cellsWithXCIstatus)] == rownames(XCIcellsAll[[whichChroms]])[na.omit(cellsWithXCIstatus)] )
    
    # any(is.na(XCIcellsAll[[whichChroms]][,"log2Expr"]))
    # [1] FALSE
    #    
    # length(XCIcellsAll[[whichChroms]][,"log2Expr"])
    # [1] 1019
    #
    # !is.na(cellsWithXCIstatus)
    # [1] 1019
    # 
    # length(na.omit(cellsWithXCIstatus))
    # [1] 1019

    theseCellIDsXCIstatus = rep(NaN, length(theseCellIDs))
    theseCellIDsXCIstatus[ !is.na(cellsWithXCIstatus) ] = XCIcellsAll[[whichChroms]][,"log2Expr"][ na.omit(cellsWithXCIstatus) ]
    pData(cdsC.OG)$XCIstatus[refBoolIdx] = as.numeric( as.character( theseCellIDsXCIstatus ))
    pData(cdsC.OG)$XCIstatus[altBoolIdx] = as.numeric( as.character( theseCellIDsXCIstatus ))
    
    # # Sort cells by XCI status for cells the have it
    # cellXCIstatusOrder = order( as.numeric( as.character( pData(cdsC.OG)$XCIstatus[refBoolIdx] ) ) )
    
    sum( is.na(pData(cdsC.OG)$XCIstatus))
    # 1766
    sum( is.na(pData(cdsC.OG)$XCIstatus) & femaleBoolIdx )
    # 58
    
    XCIboolIdx = abs(pData(cdsC.OG)$XCIstatus) > XCIthresh
    XCIboolIdx[ is.na(XCIboolIdx) ] = FALSE
    # sum(XCIboolIdx)/2
    # [1] 511

    # Sort genes by location
    genePosOrder = order( as.numeric( as.character( fData(cdsC.OG)$ensStart ) ) )
    cdsC.OG = cdsC.OG[ genePosOrder, ] 
    dim(cdsC.OG)
    
    # color scales
    logUMIcols = circlize::colorRamp2(breaks = range(log10(pData(cdsC.OG)$n.umi[refBoolIdx] + pData(cdsC.OG)$n.umi[altBoolIdx] + 0.1), na.rm=TRUE),
                                      colors = c("white", "green4"),
                                      transparency = .2)
    
    XCIcolors = circlize::colorRamp2(breaks = c( min(pData(cdsC.OG)$XCIstatus[refBoolIdx], na.rm=TRUE), 0 , max(pData(cdsC.OG)$XCIstatus[refBoolIdx], na.rm=TRUE) ),
                                     colors = c("blue", "white", "red"),
                                     transparency = .2)
    
    
    
    ##########
    # Allelic gene expression
    
    # counts matrix
    CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
    # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
    rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
    
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

    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicExpression",
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
    
    
    
    ##########
    # Allelic gene expression - female cells only 
    
    # counts matrix
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
    dim(CDSmatrix.topQuantVar)
    
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
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicExpression",
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
    # Allelic gene expression - female cells only 
    # with XCI status > 10% FPR threshold base on chr1 XCI distribution
    
    # counts matrix
    CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
    # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
    rownames(CDSmatrix) = fData(cdsC.OG)$gene_short_name
    
    CDSmatrix = CDSmatrix[ ,XCIboolIdx & femaleBoolIdx ]
    dim(CDSmatrix)
    
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
                               log10totalUMIs = log10(pData(cdsC.OG)$n.umi[refBoolIdx & femaleBoolIdx & XCIboolIdx] 
                                                      + pData(cdsC.OG)$n.umi[altBoolIdx & femaleBoolIdx & XCIboolIdx] + 0.1)[tpOrder],
                               XCIstatus = pData(cdsC.OG)$XCIstatus[refBoolIdx & femaleBoolIdx & XCIboolIdx][tpOrder], 
                               col = list( timepoint = tpCols,
                                           log10totalUMIs = logUMIcols,
                                           XCIstatus = XCIcolors ))
    
    # Can't annotate rows in this version 
    # row_ha = HeatmapAnnotation(allele = c(rep("129", sum(refBoolIdx)),
    #                                       rep("Cast", sum(altBoolIdx)) ),
    #                            col = list( allele = alleleCols))
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.heatmap.allelicExpression",
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
    # log2 ratio allelic gene expression
    
    # counts matrix
    CDSmatrix = as.matrix(Biobase::exprs(cdsC.OG))
    # all(rownames(pData(cdsC.OG)) == colnames(CDSmatrix))
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
    # log2 ratio allelic gene expression - female cells
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleCells.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by ABSOLUTE XCI log2 ratio
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleCells.absXCIorder.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleCells.XCIorder.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio
    # with XCI status > 10% FPR threshold base on chr1 XCI distribution

    # counts matrix
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
                                                 ".", theChrom, ".femaleXCIcells.XCIorder.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio - cluster genes
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleCells.XCIorder.clusteredGenes.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio - cluster genes
    # with XCI status > 10% FPR threshold base on chr1 XCI distribution
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleXCIcells.XCIorder.clusteredGenes.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio - cluster genes - all genes
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleCells.XCIorder.clusteredGenes.pdf", sep=""))
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
    # log2 ratio allelic gene expression - female cells - ranked by XCI log2 ratio - cluster genes - all genes
    # with XCI status > 10% FPR threshold base on chr1 XCI distribution
    
    # counts matrix
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
                                                 ".", theChrom, ".femaleXCIcells.XCIorder.clusteredGenes.pdf", sep=""))
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
    
    
} # chrom










#################################################
##################################################
##################################################
# Partition genes into chromosome-specific groups

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

save.image(file=file.path(mainOutputDir,".RData")) # To .RData
# load(file=file.path(mainOutputDir,".RData"))










##################################################
##################################################
##################################################
# 2020415
# Log10 + 0.1 expression distributions faceted by genome-wide, autosome and chrX
# NOTE: THIS  NO LONGER ONLY PLOTS EXPRESSED GENES (b/c log10(Texpr) drops 0 cells)
#       THIS MIGHT AFFECT PLOTS OF INDIVIUAL CHROMS (e.g. chrX)

thisOutputDir = file.path( mainOutputDir, "total_expression_distributions")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.expressionDistributions.log+0.1.faceted.pdf", sep=""))
pdf(PDFfilename, height=3, width=10)

# print( ggplot(pData(cds.OG), aes(x=log10(totalExpr))) + 
#            geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + 
#            facet_wrap( ~ timepoint, ncol = 5) + 
#            ggtitle( paste("All features (", dim(cds.OG)[1], ") and unfiltered cells (n = ", dim(cds.OG)[2], ")", sep="") ) )
# print( ggplot(pData(CDSs[["all"]]), aes(x=log10(totalExpr))) + 
#            geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + 
#            facet_wrap( ~ timepoint, ncol = 5) + 
#            ggtitle( paste("All features (", dim(CDSs[["all"]])[1], ") and filtered cells (n = ", dim(CDSs[["all"]])[2], ")", sep="") ) )

for ( whichChr in c("genome", "autosomes", "chrX" )) {
    print( ggplot(pData(CDSs[[whichChr]]), aes(x=log10(Texpr + 0.1))) + 
               geom_density(aes(fill=sex), alpha = 0.4, size=0.75) + 
               facet_wrap( ~ timepoint, ncol = 5) + 
               ggtitle( paste("Gencode ", whichChr, " genes (", dim(CDSs[[whichChr]])[1], ") and filtered cells (n = ", dim(CDSs[[whichChr]])[2], ")", sep="") ) )
}

dev.off()

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.expressionDistributions.log+0.1.faceted.allelic.pdf", sep=""))
pdf(PDFfilename, height=6, width=10)

for ( whichChr in c("genome", "autosomes", "chrX" )) {
    print( ggplot(pData(CDSs[[whichChr]]), aes(x=log10(Texpr + 0.1))) + 
               geom_density(aes(fill=allele), alpha = 0.4, size=0.75) + 
               coord_cartesian(ylim = c(0, 3)) +
               facet_wrap( ~ sex + timepoint, ncol = 5) + 
               ggtitle( paste("Gencode ", whichChr, " genes (", dim(CDSs[[whichChr]])[1], ") and filtered cells (n = ", dim(CDSs[[whichChr]])[2], ")", sep="") ) )
}

dev.off()










##################################################
##################################################
##################################################
# Plot allelic unlogged expression distributions of total expression: *** per chromosome ***

thisOutputDir = file.path( mainOutputDir, "total_expression_distributions")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

daPlots = list()
for (theChrom in theChroms) {
    # print(theChrom)
    
    chromgenesNames <- row.names(fData(CDSs[["all"]]))[ which(fData(CDSs[["all"]])$ensChr == theChrom) ]
    cdsC.OG <- CDSs[["all"]][chromgenesNames,]
    dim(cdsC.OG)
    
    pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
    
    daPlots[[theChrom]] <- ggplot(pData(cdsC.OG), aes(x=Cexpr)) + 
        geom_density(aes(fill=allele), alpha = 0.4, size=0.75) + 
        facet_wrap( ~ sex + timepoint, ncol = 5) + 
        ggtitle( paste(theChrom, " genes (", dim(cdsC.OG)[1], ") and filtered cells (n = ", dim(cdsC.OG)[2], ")", sep="") ) 
    
    # Scaled: https://stats.stackexchange.com/questions/172311/r-geom-density-values-in-y-axis
    daPlots[[paste(theChrom, ".scaled", sep="")]] <- ggplot(pData(cdsC.OG), aes(x=Cexpr)) + 
        geom_density(aes(y=..scaled.., fill=allele), alpha = 0.4, size=0.75) + 
        facet_wrap( ~ sex + timepoint, ncol = 5) + 
        ggtitle( paste(theChrom, " genes (", dim(cdsC.OG)[1], ") and filtered cells (n = ", dim(cdsC.OG)[2], ")", sep="") ) 
    
    rm(cdsC.OG)
    
}


require(gridExtra)
require(grid)

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.expressionDistributionsPerChromosome.faceted.pdf", sep=""))
pdf(PDFfilename, height=6*20, width=10*2)
grid.arrange(grobs=daPlots, nrow=20, ncol=2,
             main=paste("Expression distributions per chromosome", sep=""))
dev.off()










##################################################
##################################################
##################################################
# 2020415
# Plot allelic log10 + 0.1  expression distributions of total expression: *** per chromosome ***

thisOutputDir = file.path( mainOutputDir, "total_expression_distributions")
if ( ! file.exists( thisOutputDir ) ) {
    if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
}

daPlots = list()
for (theChrom in theChroms) {
    # print(theChrom)
    
    chromgenesNames <- row.names(fData(CDSs[["all"]]))[ which(fData(CDSs[["all"]])$ensChr == theChrom) ]
    cdsC.OG <- CDSs[["all"]][chromgenesNames,]
    dim(cdsC.OG)
    
    pData(cdsC.OG)$Cexpr <- Matrix::colSums(Biobase::exprs(cdsC.OG))  
    
    daPlots[[theChrom]] <- ggplot(pData(cdsC.OG), aes(x=log10(Cexpr + 0.1))) + 
        geom_density(aes(fill=allele), alpha = 0.4, size=0.75) + 
        facet_wrap( ~ sex + timepoint, ncol = 5) + 
        ggtitle( paste(theChrom, " genes (", dim(cdsC.OG)[1], ") and filtered cells (n = ", dim(cdsC.OG)[2], ")", sep="") ) 
    
    # Scaled: https://stats.stackexchange.com/questions/172311/r-geom-density-values-in-y-axis
    daPlots[[paste(theChrom, ".scaled", sep="")]] <- ggplot(pData(cdsC.OG), aes(x=log10(Cexpr + 0.1))) + 
        geom_density(aes(y=..scaled.., fill=allele), alpha = 0.4, size=0.75) + 
        facet_wrap( ~ sex + timepoint, ncol = 5) + 
        ggtitle( paste(theChrom, " genes (", dim(cdsC.OG)[1], ") and filtered cells (n = ", dim(cdsC.OG)[2], ")", sep="") ) 
    
    rm(cdsC.OG)
    
}


require(gridExtra)
require(grid)

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.expressionDistributionsPerChromosome.log10+0.1.faceted.pdf", sep=""))
pdf(PDFfilename, height=6*20, width=10*2)
grid.arrange(grobs=daPlots, nrow=20, ncol=2,
             main=paste("Expression distributions per chromosome", sep=""))
dev.off()










##################################################
##################################################
##################################################
# Allelic stuff
# Values used in cluster plots


##################################################
# Overall allelic ratios & pairings

for (whichChr in names(CDSs) ){
    # whichChr="autosomes"
    print(whichChr)
    
    altBoolIdx = grepl("_alt$", rownames(pData(CDSs[[whichChr]])))
    refBoolIdx = grepl("_ref$", rownames(pData(CDSs[[whichChr]])))
    
    if( sum(altBoolIdx) != sum(refBoolIdx) ) {
        cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
        next
    }
    
    pData(CDSs[[whichChr]])$cellIndex[altBoolIdx] = seq(1, sum(altBoolIdx), 1)
    pData(CDSs[[whichChr]])$cellIndex[refBoolIdx] = seq(1, sum(refBoolIdx), 1)
    
    overallAltAllelicRatios = pData(CDSs[[whichChr]])$totalExpr[altBoolIdx] / (pData(CDSs[[whichChr]])$totalExpr[altBoolIdx] + pData(CDSs[[whichChr]])$totalExpr[refBoolIdx])
    overallRefAllelicRatios = pData(CDSs[[whichChr]])$totalExpr[refBoolIdx] / (pData(CDSs[[whichChr]])$totalExpr[altBoolIdx] + pData(CDSs[[whichChr]])$totalExpr[refBoolIdx])
    
    pData(CDSs[[whichChr]])$allelicRatio[altBoolIdx] <- overallAltAllelicRatios
    pData(CDSs[[whichChr]])$allelicRatio[refBoolIdx] <- overallRefAllelicRatios
    
    print(summary(pData(CDSs[[whichChr]])$allelicRatio[altBoolIdx]))
    print(summary(pData(CDSs[[whichChr]])$allelicRatio[refBoolIdx]))
    
    rm( altBoolIdx,
        refBoolIdx,
        overallAltAllelicRatios,
        overallRefAllelicRatios
    )
    
    gc()
    
}









    
##################################################
##################################################
##################################################
# Proportion of skewed/monoallelic genes 

# dim( Biobase::exprs(CDSs[[whichChr]]) )
# colnames( Biobase::exprs(CDSs[[whichChr]]) )

for (whichChr in names(CDSs) ){
    # whichChr="autosomes"
    print(whichChr)
    
    refBoolIdx = grepl("_ref$", rownames(pData(CDSs[[whichChr]])))
    altBoolIdx = grepl("_alt$", rownames(pData(CDSs[[whichChr]])))
    
    if( sum(altBoolIdx) != sum(refBoolIdx) ) {
        cat( sum(altBoolIdx), " != ", sum(refBoolIdx), "\n", sep="")
        next
    }
    
    refExpr = as.matrix(Biobase::exprs(CDSs[[whichChr]])[,refBoolIdx])
    altExpr = as.matrix(Biobase::exprs(CDSs[[whichChr]])[,altBoolIdx])
    
    # biallelic expression 
    refExprBoolIdx = refExpr > 0 
    altExprBoolIdx = altExpr > 0 
    biallelicExpr = refExprBoolIdx & altExprBoolIdx
    expressedGeneBoolIdx = refExprBoolIdx | altExprBoolIdx
    propRefExpr = colSums(refExprBoolIdx, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    propAltExpr = colSums(altExprBoolIdx, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    propBiallelic = colSums(biallelicExpr, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    propMonoallelic = 1 - propBiallelic
    
    pData(CDSs[[whichChr]])$propMonoAllelicGenes[altBoolIdx] <- propMonoallelic
    pData(CDSs[[whichChr]])$propMonoAllelicGenes[refBoolIdx] <- propMonoallelic
    pData(CDSs[[whichChr]])$propRefOnlyExprGenes[altBoolIdx] <- propRefExpr
    pData(CDSs[[whichChr]])$propRefOnlyExprGenes[refBoolIdx] <- propRefExpr
    pData(CDSs[[whichChr]])$propAltOnlyExprGenes[altBoolIdx] <- propAltExpr
    pData(CDSs[[whichChr]])$propAltOnlyExprGenes[refBoolIdx] <- propAltExpr
    
    # skewed expression
    dScore = (altExpr / (altExpr + refExpr)) - 0.5 # Range -0.5 -> 0.5 ~ ref skewed -> alt skewed
    refSkewedBoolIdx = dScore <= - 0.35
    altSkewedBoolIdx = dScore >= 0.35
    skewedGenesBoolIdx = refSkewedBoolIdx | altSkewedBoolIdx
    propRefSkewed = colSums(refSkewedBoolIdx, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    propAltSkewed = colSums(altSkewedBoolIdx, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    propSkewed = colSums(skewedGenesBoolIdx, na.rm=TRUE) / colSums(expressedGeneBoolIdx, na.rm=TRUE)
    
    pData(CDSs[[whichChr]])$propSkewedGenes[altBoolIdx] <- propSkewed
    pData(CDSs[[whichChr]])$propSkewedGenes[refBoolIdx] <- propSkewed
    pData(CDSs[[whichChr]])$propRefSkewedGenes[altBoolIdx] <- propRefSkewed
    pData(CDSs[[whichChr]])$propRefSkewedGenes[refBoolIdx] <- propRefSkewed
    pData(CDSs[[whichChr]])$propAltSkewedGenes[altBoolIdx] <- propAltSkewed
    pData(CDSs[[whichChr]])$propAltSkewedGenes[refBoolIdx] <- propAltSkewed
    
    # Percentage of expressed genes
    pData(CDSs[[whichChr]])$percentageExpressedGenes[altBoolIdx] <- colSums(expressedGeneBoolIdx) / nrow(expressedGeneBoolIdx) * 100
    pData(CDSs[[whichChr]])$percentageExpressedGenes[refBoolIdx] <- colSums(expressedGeneBoolIdx) / nrow(expressedGeneBoolIdx) * 100
    
    print(summary(pData(CDSs[[whichChr]])$percentageExpressedGenes[altBoolIdx] ))
    print(summary(pData(CDSs[[whichChr]])$propSkewedGenes[altBoolIdx]))
    print(summary(pData(CDSs[[whichChr]])$propRefSkewedGenes[altBoolIdx]))
    print(summary(pData(CDSs[[whichChr]])$propAltSkewedGenes[altBoolIdx]))
    
    rm(    refBoolIdx,
           altBoolIdx,
           refExpr,
           altExpr,
           refExprBoolIdx,
           altExprBoolIdx,
           biallelicExpr,
           expressedGeneBoolIdx,
           propRefExpr,
           propAltExpr,
           propBiallelic,
           propMonoallelic,
           dScore,
           refSkewedBoolIdx,
           altSkewedBoolIdx,
           skewedGenesBoolIdx,
           propRefSkewed,
           propAltSkewed,
           propSkewed)
    gc()
    
}










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
# Allelic CDSs

allelicCDSs = list()
for (whichChr in names(CDSs) ){
    
    if (whichChr %in% c("all") ) {
        next 
    }
    
    alleleBoolIdx = list()
    alleleBoolIdx[["ref"]] = grepl("_ref$", rownames(pData(CDSs[[whichChr]])))
    alleleBoolIdx[["alt"]] = grepl("_alt$", rownames(pData(CDSs[[whichChr]])))
    
    for(theAllele in names(alleleBoolIdx)) {
        # theAllele="ref"
        
        allelicCDSs[[paste(whichChr, ".", theAllele, sep="")]] = CDSs[[whichChr]][,alleleBoolIdx[[theAllele]] ]
    }
}










##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- Semi-supervised mode -- allelic
# Use the known cell state!

for ( whichChr in names(allelicCDSs) ){
    print(whichChr)
    
    thisOutputDir = file.path(mainOutputDir, paste("trajectoryAnalysis.semisupervised.", whichChr, sep="") )
    if ( ! file.exists( thisOutputDir ) ) {
        if ( ! dir.create( thisOutputDir, showWarnings=FALSE, recursive=TRUE )  ) thisOutputDir = mainOutputDir
    }
    
    # save.image(file =  file.path(thisOutputDir, paste( ".RData.semisupervised.", whichChr, sep="") )
    # load(file =  file.path(thisOutputDir, paste( ".RData.semisupervised.", whichChr, sep="") )
    
    ####################
    ####################
    # Trajectory step 1: choose genes that define a cell's progress
    # First, we must decide which genes we will use to define a cell's progress.
    
    # What limit?
    minExprThreshold=2 # 7079 genes; 7475 cells
    allelicCDSs[[whichChr]] <- detectGenes(allelicCDSs[[whichChr]], min_expr = min_expr)
    sum( fData(allelicCDSs[[whichChr]])$num_cells_expressed >= 10)
    sum( pData(allelicCDSs[[whichChr]])$num_genes_expressed >= 10)
    
    expressedGenes <-  row.names(subset(fData(allelicCDSs[[whichChr]]), num_cells_expressed >= min_cells))
    length(expressedGenes)
    
    # diff_test_res <- differentialGeneTest(allelicCDSs[[whichChr]][expressed_genes,], fullModelFormulaStr = "~Media")
    timepoint_DEG_genes <- differentialGeneTest(allelicCDSs[[whichChr]][expressedGenes,],
                                                fullModelFormulaStr = '~timepoint',
                                                cores = 1)
    ordering_timepoint_DEG_genes<- row.names(subset(timepoint_DEG_genes, qval < 0.01))
    length(ordering_timepoint_DEG_genes)
    
    # Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
    # because the next several functions will depend on them.
    allelicCDSs[[whichChr]] <- setOrderingFilter(allelicCDSs[[whichChr]], ordering_timepoint_DEG_genes)
    
    
    ####################
    ####################
    # Trajectory step 2: reduce data dimensionality
    # Next, we will reduce the space down to one with two dimensions,
    # which we will be able to easily visualize and interpret while Monocle is ordering the cells.
    allelicCDSs[[whichChr]] <- reduceDimension(allelicCDSs[[whichChr]], max_components = 2, method = 'DDRTree')
    
    
    ####################
    ####################
    # Trajectory step 3: order cells along the trajectory
    # Now that the space is reduced, it's time to order the cells using the orderCells function as shown below.
    allelicCDSs[[whichChr]] <- orderCells(allelicCDSs[[whichChr]])
    
    
    ##################################################
    # "State" is just Monocle's term for the segment of the tree.
    GM_state <- function(cds){
        if (length(unique(pData(cds)$State)) > 1){
            T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
            return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
        } else {
            return (1)
        }
    }
    allelicCDSs[[whichChr]] <- orderCells(allelicCDSs[[whichChr]], root_state = GM_state(allelicCDSs[[whichChr]]))
    
    
    ##################################################
    # Tabulate results
    
    tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".tabulatedResults.byTimePoint.tsv", sep=""))
    write.table( x = table(pData(allelicCDSs[[whichChr]])$State, pData(allelicCDSs[[whichChr]])$timepoint),
                 file = tablefilename, sep="\t" )
    
    tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".tabulatedResults.byTimePointAndGender.tsv", sep=""))
    write.table( x = table(pData(allelicCDSs[[whichChr]])$State, pData(allelicCDSs[[whichChr]])$timepoint, pData(allelicCDSs[[whichChr]])$sex),
                 file = tablefilename, sep="\t" )
    
    
    
    ##################################################
    # Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".byTimePoint.plotTrajectories.pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    # The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
    # near one of the tips of the tree, while the others are distributed amongst the two "branches". 
    # Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
    # so we often have to call orderCells again using the root_state argument to specify the beginning.
    # First, we plot the trajectory, this time coloring the cells by "State":
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "State") + geom_point(size=1, aes(color = State)) )
    # If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
    # Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~State, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~sex, nrow = 1) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "State") + geom_point(size=1, aes(color = State)) + facet_wrap(~timepoint, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~State, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "Pseudotime") + geom_point(size=1, aes(color = Pseudotime)) + facet_wrap(~sex, nrow = 1) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "timepoint") + geom_point(size=1, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "sex") + geom_point(size=1, aes(color = sex)) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~sex, nrow = 1) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~State, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "sex") + geom_point(size=1, aes(color = sex)) + facet_wrap(~timepoint, nrow = 2) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "sample") + geom_point(size=1, aes(color = sample)) )
    print( plot_cell_trajectory(allelicCDSs[[whichChr]], color_by = "batchNum") + geom_point(size=1, aes(color = batchNum)) )
    dev.off()
    
    
    
    
    
    ##################################################
    # And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, 
    # using your biological knowledge of the system. 
    # For example, in this experiment, a highly proliferative population of progenitor cells are generating two types of post-mitotic cells. 
    # So the root should have cells that express high levels of proliferation markers.
    # We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:
    
    markerGenes = c("Pou5f1", "Sox2", "Nanog", "Klf4",
                    "Eomes", # T NPC intermediate marker!
                    "Sox1",  # Sox1 NPC marker!
                    "Otx2", "Sema5a", "Snhg18", "Xist",
                    "Pax6", "Map2" )  # NPC
    markerGenesNames <- row.names(subset(fData(allelicCDSs[[whichChr]]),  gene_short_name %in% markerGenes))
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByState.pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    print( plot_genes_jitter(allelicCDSs[[whichChr]][markerGenesNames,],
                             grouping = "State",
                             min_expr = min_expr) )
    dev.off()
    
    # To confirm that the ordering is correct we can select a couple of markers of myogenic progress. 
    # Plotting these genes demonstrates that ordering looks good:
    expressedGenes <-  row.names(subset(fData(allelicCDSs[[whichChr]]), num_cells_expressed >= 10))
    allelicCDS.expressed <- allelicCDSs[[whichChr]][expressedGenes,]
    markerGenesNames.expressed <- row.names(subset(fData(allelicCDS.expressed),
                                                   gene_short_name %in% markerGenes))
    allelicCDS.expressed.subset <- allelicCDS.expressed[markerGenesNames.expressed,]
    fData(allelicCDS.expressed.subset )
    # gene_short_name                   id num_cells_expressed
    # ENSMUSG00000086503.3            Xist ENSMUSG00000086503.3                 726
    # use_for_ordering
    # ENSMUSG00000086503.3             TRUE
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".byTimePoint.expressionPatternsByPseudotime.pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    print( plot_genes_in_pseudotime(allelicCDS.expressed.subset, color_by = "State") + geom_point(size=0.25, aes(color = State)) )
    print( plot_genes_in_pseudotime(allelicCDS.expressed.subset, color_by = "timepoint") + geom_point(size=0.25, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
    print( plot_genes_in_pseudotime(allelicCDS.expressed.subset, color_by = "sex") + geom_point(size=0.25, aes(color = sex)) )
    print( plot_genes_in_pseudotime(allelicCDS.expressed.subset, color_by = "sample") + geom_point(size=0.25, aes(color = sample)) )
    dev.off()
    
    
    
    
    
    
    
    
    
    ##################################################
    ##################################################
    ##################################################
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
    
    # Takes an age if cores used b/c can't allocate memory!
    pseudoTime_diff_test_res <- differentialGeneTest(allelicCDSs[[whichChr]][markerGenesNames,],
                                                     # pseudoTime_diff_test_res <- differentialGeneTest(allelicCDS.expressed.subset,
                                                     fullModelFormulaStr = "~ sm.ns(Pseudotime)",
                                                     cores = 1)
    pseudoTime_sigGenesNames <- row.names(subset(pseudoTime_diff_test_res, qval < 0.1))
    
    PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis.", whichChr, ".minExprThreshold", minExprThreshold, ".byTimePoint.pseudotimeHeatmap.coVaryingGeneClusters.pdf", sep=""))
    pdf(PDFfilename, height=10, width=10)
    print( plot_pseudotime_heatmap(allelicCDSs[[whichChr]][pseudoTime_sigGenesNames,],
                                   num_clusters = 5,
                                   cores = 1,
                                   show_rownames = TRUE) )
    dev.off()
    
} # for ( hichChr in names(allelicCDSs) ){
