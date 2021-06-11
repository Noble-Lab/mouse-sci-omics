##################################################
##################################################
##################################################
#
# 20200720 allelic sci-Hi-C analysis using Monocle2
#          using concatened chromosome-specific CDP vectors.
#          Initially for F121 cells only.
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

projectName = "allelic_sci-Hi-C"

inputDir = "../data/allelic_catenated_contact_decay_profiles/"

whichChrom="genome"
# whichChrom="chr1_and_chrX"
# whichChrom="chrX"
# whichChrom="chr1"
# whichChrom="autosomes"

# 20200921
cellTypesIncluded="F121only"

# 20200916
# cellTypesIncluded="F121only.noNPCs"

mainOutputDir = paste0("../results/allelic_sci-Hi-C_monocle2_", cellTypesIncluded, "_", whichChrom )
if ( ! file.exists( mainOutputDir ) ) {
  if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

pchVal=20
pchcexVal=1.25
altpchVal=17

tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")
tpShapes = rep( pchVal, length(tpCols) )
names(tpShapes) = names(tpCols)
 
alleleCols=c("darkred", "darkblue")
names(alleleCols) = c("ref", "alt")










##################################################
##################################################
##################################################
# Load counts data

cellTypeFileString = cellTypesIncluded
if ( grepl( "F121", cellTypesIncluded ) ) cellTypeFileString = "F121"

# 20200721
if (whichChrom == "genome") {

    count.table.OG = c()
    
    for (theChrom in c( paste0( "chr", 1:19 ), "chrX" ) ) { 
        
        count.table.perSample = c()
        
        # 20210320
        # cellTypesIncluded="F121_and_TsixStop"
        for (aCellTypeFileString in cellTypeFileString) {
            
            inputFileName = paste0("allelicContactDecayRates.catenatedAllelicCDPs.", aCellTypeFileString, ".",
                                   theChrom, ".tsv.gz" )
            countsFile = file.path(inputDir, inputFileName)
            
            count.table.OG.tmp = as.data.frame(fread(countsFile, header=TRUE, sep="\t") )
            
            rownames(count.table.OG.tmp) = count.table.OG.tmp[,1]
            count.table.OG.tmp = count.table.OG.tmp[,-1]
            
            rownames(count.table.OG.tmp)  = paste0( sub( "_" , paste0( "_", theChrom, ":"), rownames(count.table.OG.tmp) ),
                                                    "-", sub("^.*_", "", rownames(count.table.OG.tmp) ) )
            # print( dim(count.table.OG.tmp) ) 
            
            if (is.null(count.table.perSample)) {
                count.table.perSample = count.table.OG.tmp
            } else {
                count.table.perSample = cbind( count.table.perSample, 
                                               count.table.OG.tmp )
            }
                
        }
        
        count.table.OG = rbind( count.table.OG,
                                count.table.perSample )
    }
         
# 20200821
} else if (whichChrom == "autosomes") {
    
	count.table.OG = c()
	
	for (theChrom in c( paste0( "chr", 1:19 ) ) ) { 
		
		count.table.perSample = c()
		
		# 20210320
		# cellTypesIncluded="F121_and_TsixStop"
		for (aCellTypeFileString in cellTypeFileString) {
			
			inputFileName = paste0("allelicContactDecayRates.catenatedAllelicCDPs.", aCellTypeFileString, ".",
								   theChrom, ".tsv.gz" )
			countsFile = file.path(inputDir, inputFileName)
			
			count.table.OG.tmp = as.data.frame(fread(countsFile, header=TRUE, sep="\t") )
			
			rownames(count.table.OG.tmp) = count.table.OG.tmp[,1]
			count.table.OG.tmp = count.table.OG.tmp[,-1]
			
			rownames(count.table.OG.tmp)  = paste0( sub( "_" , paste0( "_", theChrom, ":"), rownames(count.table.OG.tmp) ),
													"-", sub("^.*_", "", rownames(count.table.OG.tmp) ) )
			# print( dim(count.table.OG.tmp) ) 
			
			if (is.null(count.table.perSample)) {
				count.table.perSample = count.table.OG.tmp
			} else {
				count.table.perSample = cbind( count.table.perSample, 
											   count.table.OG.tmp )
			}
			
		}
		
		count.table.OG = rbind( count.table.OG,
								count.table.perSample )
	}
        
} else if (whichChrom == "chr1_and_chrX") {

    count.table.OG = c()
    for (theChrom in c( "chr1", "chrX" ) ) { 
        inputFileName = paste0("allelicContactDecayRates.catenatedAllelicCDPs.", cellTypeFileString, ".",
                               theChrom, ".tsv.gz" )
        countsFile = file.path(inputDir, inputFileName)
        
        count.table.OG.tmp = as.data.frame(fread(countsFile, header=TRUE, sep="\t") )
        
        rownames(count.table.OG.tmp) = count.table.OG.tmp[,1]
        count.table.OG.tmp = count.table.OG.tmp[,-1]
        
        rownames(count.table.OG.tmp)  = paste0( sub( "_" , paste0( "_", theChrom, ":"), rownames(count.table.OG.tmp) ),
                                                "-", sub("^.*_", "", rownames(count.table.OG.tmp) ) )
        # print( dim(count.table.OG.tmp) ) 
        
        count.table.OG = rbind( count.table.OG,
                                count.table.OG.tmp )
    }
    
} else {
    
    inputFileName = paste0("allelicContactDecayRates.catenatedAllelicCDPs.", cellTypeFileString, ".",
                           whichChrom, ".tsv.gz" )
    countsFile = file.path(inputDir, inputFileName)
    
    count.table.OG = as.data.frame(fread(countsFile, header=TRUE, sep="\t" ))
    # count.table.OG = as.data.frame(read.table(gzfile(countsFile), header=TRUE, sep="\t" ))
    
    rownames(count.table.OG) = count.table.OG[,1]
    count.table.OG = count.table.OG[,-1]
    
    rownames(count.table.OG)  = paste0( sub( "_" , paste0( "_", whichChrom, ":"), rownames(count.table.OG) ),
                                        "-", sub("^.*_", "", rownames(count.table.OG) ) )
}

# count.table.OG[1:3, 1:3]

dim(count.table.OG)
# 286 2153 - chrX
# 286 2153 - chr1
#
# Note: 286/2 = 143 cis log bins

# 572 2153 - chr1_and_chrX
# 5720 2153 - genome
# 5434 2153 - autosomes










##################################################
##################################################
##################################################
# Generate CDS data

sapply(count.table.OG, class)
# as.matrix(count.matrix.OG[1:3,1:3])

# require(Matrix)
count.matrix.OG.sparse = Matrix::Matrix(as.matrix(count.table.OG), sparse=TRUE)

# count.table.OG = Matrix::Matrix(count.matrix.OG.sparse, sparse=FALSE)
# dim(count.table.OG)
# # [1] 157060   2836

locus.annotations = as.data.frame( cbind( rownames(count.table.OG), rownames(count.table.OG)), stringsAsFactors = FALSE )
colnames(locus.annotations) = c("id", "gene_short_name")
rownames(locus.annotations) = locus.annotations$id

cell.annotations = as.data.frame( cbind( colnames(count.table.OG), "sciHiC"), stringsAsFactors = FALSE )
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
# 2153    4

dim(fData(cds.OG))
# 572   2

dim(Biobase::exprs(cds.OG))
#  572 2153










##################################################
##################################################
##################################################
# Add metadata

# time-points
timepoints = as.matrix( sub( "^.*_", "", sub( "\\..*$", "", rownames(pData(cds.OG)) ) ) )
rownames(timepoints) <- rownames(pData(cds.OG))
colnames(timepoints) <- "timepoints"
table(timepoints)
# timepoints
# d0 d11  d3  d7 NPC 
# 986 164 435 243 325 

all( colnames(cds.OG) == rownames(timepoints) )
# TRUE

pData(cds.OG)$timepoint <- factor( timepoints, levels=names(tpCols) )










##################################################
##################################################
##################################################
# Normalize and pre-process the data

cds.OG <- estimateSizeFactors(cds.OG)
pData(cds.OG)[1:5,]

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

# The vector expressed_genes now holds the identifiers for genes expressed in at least 10 cells of the data set.
accessible_loci <- row.names(subset(fData(cds.OG),
                                    num_cells_expressed >= 10))
length(accessible_loci)
# 566

print(head(pData(cds.OG)))










##################################################
##################################################
##################################################
# Filter out extreme cells

# Total expression per barcode
pData(cds.OG)$totalCoverage <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!

summary(pData(cds.OG)$totalCoverage )

# chr1_and_chrX
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 268.0   460.0   641.0   946.4   973.0 33159.0 
#
# chrX
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 101.0   144.0   196.0   281.6   289.0 10054.0 
#
# chr1
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 159.0   314.0   442.0   664.8   680.0 23105.0
#
# genome
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2233    4006    5661    8386    8612  281463 
#
# autsomes
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2127    3863    5463    8105    8310  271409 


set_lower_bound = 0 
# chr1_and_chrX
# set_lower_bound = 460 # Q1
#
# chrX 
# set_lower_bound = 196 # Q2
#
# chr1
# set_lower_bound = 442 # Q2
#
# genome
# set_lower_bound = 5661 # Q2

# rm(cds)
upper_bound <- 10^(mean(log10(pData(cds.OG)$totalCoverage)) +
                       2*sd(log10(pData(cds.OG)$totalCoverage)))
lower_bound <- 10^(mean(log10(pData(cds.OG)$totalCoverage)) -
                       2*sd(log10(pData(cds.OG)$totalCoverage)))

if (set_lower_bound > 0) lower_bound = set_lower_bound

# Filter
cds <- cds.OG[,pData(cds.OG)$totalCoverage > lower_bound &
                  pData(cds.OG)$totalCoverage < upper_bound]

dim(cds.OG)

dim(cds)
# chr1_and_chrX
# Features  Samples 
# 572     2056 
# 
# set_lower_bound = 460
# Features  Samples 
# 572     1513 
#
#
#
# chrX
# Features  Samples 
# 286     2056 
# 
# set_lower_bound = 196
# Features  Samples 
# 286      974 
#
#
#
# chr1
# Features  Samples 
# 286     2055
# 
# set_lower_bound = 196
# Features  Samples 
# 286      974 
#
#
#
# genome
# Features  Samples 
# 5720     2056 
#
# set_lower_bound = 5661
# Features  Samples 
# 5720      978 

# autsomes
# Features  Samples 
# 5434     2056 


# 20200916 Genome-wide trajectory analysis *** EXCLUDING NPCs ***
if (cellTypesIncluded == "F121only.noNPCs") {
    cds = cds[,pData(cds)$timepoint != "NPC"]
    dim(cds)
    # Features  Samples 
    # 5720     1739 










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
p3 <- ggplot(pData(cds), aes(x=log10(totalCoverage))) + geom_density(aes(color=timepoint), size=0.75) + scale_color_manual(values = tpCols)

PDFfilename = file.path(mainOutputDir, 
                        paste(projectName, ".monocle2results.coverageDistributions",
                              ifelse(set_lower_bound >0, paste0( ".set_lower_bound", set_lower_bound ), ""),
                              ".pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
multiplot(p1, p2, p3, layout=matrix(c(1,1,2,2,3,3,3,3), nrow=2, byrow=TRUE))
dev.off()









##################################################
##################################################
##################################################
# save.image() # To .RData
# load(file = ".RData")










##################################################
##################################################
##################################################
# Run monocle2


##################################################
##################################################
##################################################
# Constructing Single Cell Trajectories -- Semi-supervised mode

# set_lower_bound = 0
thisOutputDir = file.path(mainOutputDir, 
                          paste0( "trajectoryAnalysis.semisupervised",
                                  ifelse(set_lower_bound >0, paste0( ".set_lower_bound", set_lower_bound ), "")
                          ) )
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

# 20200516
# minExprThreshold=0.1; numCellsExpressed=5
#
# 20200517: Tried this to reduce number of loci prevent process getting killed. 
#
# 20200518: Now downsampling.


minExprThreshold=0.1
# minExprThreshold=1
numCellsExpressed=10

# 20200518
# Try downsampling:
# dsProp=NA
dsProp=1
# dsProp=0.5



##################################################
# Run analysis 

cds <- detectGenes(cds, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cds), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)

# 565 - F121 chr1_chrX
# 281 - F121 chrX
# 5491 - F121 genome
# 5210 - F121 autosome
# 5481 - F121 genome no NPCs

set.seed(123)
coveredLoci = sample( coveredLoci, size=floor(dsProp * length(coveredLoci)) )
length(coveredLoci)

# diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci <- differentialGeneTest(cds[coveredLoci,],
                                          fullModelFormulaStr = '~timepoint',
                                          cores = 1)

DEqval=10^-11 # ***
# DEqval=10^-5 # chrX

ordering_timepoint_DE_loci<- row.names(subset(timepoint_DE_loci, qval < DEqval))
length(ordering_timepoint_DE_loci)

###############
# F121 -chrX
# 20210329 - 3 branches!
# 183 # 10^-4.5 ***

###############
# F121 - genome
# 2040 # 10^-11 ***

###############
# F121 - autosomes
# 1924 # 10^-11 ***

###############
# F121 genome - noNPCs
# 2199 # 10^-9 ***



##################################################
# 2021031
# Save DE bins
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ".DEbins.txt", sep=""))
write.table(sort(ordering_timepoint_DE_loci), 
            tablefilename, row.names = FALSE, col.names = FALSE, quote=FALSE)
#
##################################################


# Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object,
# because the next several functions will depend on them.
cds <- setOrderingFilter(cds, ordering_timepoint_DE_loci)

# rm(timepoint_DE_loci); gc()
# print(object.size(x=lapply(ls(), get)), units="Mb")



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

##################################################
# "State" is just Monocle's term for the segment of the tree.
if ( grepl( "^F121_and_TsixStop", cellTypesIncluded ) || grepl( "^newF121_", cellTypesIncluded ) ) {
    # 20210322
    GM_state <- function(cds){
        if (length(unique(pData(cds)$State)) > 1){
            T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[, "F121_d0"]
            return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
        } else {
            return (1)
        }
    }
} else {
    GM_state <- function(cds){
        if (length(unique(pData(cds)$State)) > 1){
            T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
            return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
        } else {
            return (1)
        }
    }
}

cds <- orderCells(cds, root_state = GM_state(cds))


##################################################
# Tabulate results

table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]

###############
# F121 chrX
#
######
# 20210329 - 3 branches!
#
# 183 # 10^-4.5 ***
# d0  d3  d7 d11 NPC
# 1 668  24  22   6   4
# 2  92 280  66  25   5
# 3 212 107 122 115 308
#
# 157 # 10^-6.5 **
# d0  d3  d7 d11 NPC
# 2 306   2   4   1   1
# 3 257 313  91  36   9
# 1 409  96 115 109 307

###############
# F121 genome
#
# DEqval=10^-11 ***
# d0  d3  d7 d11 NPC
# 2 781  10  18   7  37 - early == 1
# 1 108 358  94  30  37 - intermediate == 2
# 3  82  50  96 105 243 - late == 3

###############
# F121 autosomes
#
# 10^-11 ***
# d0  d3  d7 d11 NPC
# 1 758   4  10   5  35
# 3 147 375 120  36  72
# 2  66  39  78 101 210

###############
# F121 genome - noNPCs
#
# DEqval=10^-9 ***
# d0  d3  d7 d11 NPC
# 4 812  13  20   7   0
# 3  72   4  12  13   0
# 5  62 217  55  18   0
# 1  23 182 106 100   0
# 2   2   2  15   4   0



##################################################
# Table of cell counts per branch
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".byTimePoint.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



##################################################
# Table of cells per branch - proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ".byTimePoint.cellCountsPerBranch.proportions.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res <- table_res[order(table_res[,1], decreasing=TRUE), ]
table_res <- t( t( table_res ) / colSums(table_res) )
table_res
write.table(table_res, 
            tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



##################################################
#  20210320 Determine F121 batches / replicates 

unique(sub("-.*$", "", pData(cds)$cell))
# [1] "F121_d0.F121_PCR14_R11_pooled"  "F121_d0.F121_PCR14_R7_pooled"  
# [3] "F121_d0.R3_lane1"               "F121_d0.R3_lane2"              
# [5] "F121_d0.R4_lane1"               "F121_d0.R4_lane2"              
# [7] "F121_d3.F121_PCR14_R11_pooled"  "F121_d3.F121_PCR14_R7_pooled"  
# [9] "F121_d7.F121_PCR14_R11_pooled"  "F121_d7.F121_PCR14_R7_pooled"  
# [11] "F121_d11.F121_PCR14_R11_pooled" "F121_d11.F121_PCR14_R7_pooled" 
# [13] "F121_NPC.R3_lane1"              "F121_NPC.R3_lane2"             
# [15] "F121_NPC.R5_lane1"              "F121_NPC.R5_lane2"

batch_number = rep( "run2", length(pData(cds)$cell) )
batch_number[ grep( "F121_PCR14", pData(cds)$cell) ] = "run1"

pData(cds)$batch_number <- factor( batch_number, levels=c("run1", "run2") )
pData(cds)$timepoint_and_batch =  factor( paste( pData(cds)$timepoint, pData(cds)$batch_number, sep="_" ),
										  levels = unique(paste( pData(cds)$timepoint, pData(cds)$batch_number, sep="_" )) )

batch_cols = c("black", "gray")
names( batch_cols ) = c("run1", "run2")

unique( pData(cds)$timepoint_and_batch )
# "d0_run1"  "d0_run2"  "d3_run1"  "d7_run1"  "d11_run1" "NPC_run2"

timepoint_and_batch_cols = c( "red4", tpCols )
names( timepoint_and_batch_cols ) = unique( pData(cds)$timepoint_and_batch )

d0_batch_cols = timepoint_and_batch_cols
d0_batch_cols[ ! grepl( "d0", names(d0_batch_cols) ) ]  = rgb(t(col2rgb(c("white"))), alpha=0, maxColorValue = 255)

batch_shapes = tpShapes
names(batch_shapes) = names( batch_cols )

timepoint_and_batch_shapes = c(pchVal, tpShapes)
names(timepoint_and_batch_shapes) = names(timepoint_and_batch_cols)

d0_batch_shapes = c(pchVal, tpShapes)
names(d0_batch_shapes) = names(d0_batch_cols)



##################################################
#  20210320 
# Branch counts and proportons for F121 batch / replicate 

##################################################
# Table of cell counts per branch with replicates 
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".byTimePoint_and_batch.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$timepoint_and_batch)
table_res[order(table_res[,1], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)


##################################################
# Table of cells per branch - proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".byTimePoint_and_batch.cellCountsPerBranch.proportions.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$timepoint_and_batch)
table_res <- table_res[order(table_res[,1], decreasing=TRUE), ]
table_res <- t( t( table_res ) / colSums(table_res) )
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



##################################################
# 20210329 
# Branch counts and proportons for all time points

##################################################
# Table of cell counts per branch with replicates 
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".byAllTimePoints.cellCountsPerBranch.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$alltimepoints)
table_res[order(table_res[,1], decreasing=TRUE), ]
write.table(table_res[order(table_res[,1], decreasing=TRUE), ], 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)

   
##################################################
# Table of cells per branch - proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".byAllTimePoints.cellCountsPerBranch.proportions.tsv", sep=""))
table_res <- table(pData(cds)$State, pData(cds)$alltimepoints)
table_res <- table_res[order(table_res[,1], decreasing=TRUE), ]
table_res <- t( t( table_res ) / colSums(table_res) )
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



##################################################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
print( plot_ordering_genes(cds) )
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cds, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) )
print( plot_cell_trajectory(cds, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + 
           geom_point(size=2, aes(color = timepoint, shape = timepoint)) + 
           scale_color_manual(values = tpCols) +
           scale_shape_manual(values = tpShapes) )
print( plot_cell_trajectory(cds, color_by = "totalCoverage") + geom_point(size=2, aes(color = totalCoverage)))

if ( cellTypesIncluded == "F121only" || grepl( "F121_and_TsixStop", cellTypesIncluded ) || grepl( "^newF121_", cellTypesIncluded ) ) {
    print( plot_cell_trajectory(cds, color_by = "batch_number") + 
               geom_point(size=2, aes(color = batch_number, shape = batch_number)) + 
               scale_color_manual(values = batch_cols) +
               scale_shape_manual(values = batch_shapes) )
    print( plot_cell_trajectory(cds, color_by = "timepoint_and_batch") + 
               geom_point(size=2, aes(color = timepoint_and_batch, shape = timepoint_and_batch)) + 
               scale_color_manual(values = timepoint_and_batch_cols) +
               scale_shape_manual(values = timepoint_and_batch_shapes) )
    print( plot_cell_trajectory(cds, color_by = "timepoint_and_batch") + 
               geom_point(size=2, aes(color = timepoint_and_batch, shape = timepoint_and_batch)) + 
               scale_color_manual(values = d0_batch_cols)+
               scale_shape_manual(values = d0_batch_shapes) )
}

dev.off()


##################################################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space - faceted.

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds, color_by = "State") + 
           geom_point(size=2, aes(color = State)) + 
           facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + 
           geom_point(size=2, aes(color = timepoint, shape = timepoint)) + 
           scale_color_manual(values = tpCols) +
           scale_shape_manual(values = tpShapes) +
           facet_wrap(~State, nrow = 1) )
dev.off()


##################################################
# 20210329
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space - faceted.
# All time points

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".byAllTimePoints.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds, color_by = "alltimepoints") + 
		   geom_point(size=2, aes(color = alltimepoints)) + 
		   scale_color_manual(values = allTPcols) +
		   facet_wrap(~alltimepoints, nrow = 2) )
dev.off()



##################################################
# 20210324
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space - faceted by time point

timepoint_facet_rows = NULL
timepoint_facet_cols = 5

if ( nchar( paste0(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
		   ".numCellsExpressed", numCellsExpressed,
		   ".downSampleProp", dsProp, 
		   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
		   ".byTimePoint.plotTrajectories.faceted_by_timepoint.pdf" ) ) > 200 ){ 
	PDFfilename = file.path(thisOutputDir, paste0(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".FacetedByTimePoint.plotTrajectories.pdf"))
} else {
	PDFfilename = file.path(thisOutputDir, paste0(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".byTimePoint.plotTrajectories.faceted_by_timepoint.pdf"))
}
pdf(PDFfilename, height=10, width=10)
# If there are a ton of time points, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the cells are located:
print( plot_cell_trajectory(cds, color_by = "timepoint") + 
		   geom_point(size=0.1, aes(color = timepoint, shape = timepoint)) + 
		   scale_color_manual(values = tpCols) +
		   scale_shape_manual(values = tpShapes) +
		   facet_wrap(~timepoint, nrow = timepoint_facet_rows, ncol = timepoint_facet_cols) )
dev.off()



##################################################
# Pseudotime CDFs per branch.

##################################################
# 20210308
# Relabel branches here rather than below.
#
# And realized I could try to use the CDFs to automatically determined distal thresholds:
# By using ecdf inflection point for distal thresholds.

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cds)[, c("cell", "timepoint", "Pseudotime", "State", "totalCoverage") ]

ptData$branch_relative_Pseudotime = NaN
for ( aState in sort( unique( ptData$State ) ) ) {
	ptData$branch_relative_Pseudotime[ ptData$State == aState ] =  zero2one( ptData$Pseudotime[ ptData$State == aState ] )
	if( min( ptData$Pseudotime[ ptData$State == aState ] ) == min(  ptData$Pseudotime ) ) {
		ptData$branch_relative_Pseudotime[ ptData$State == aState ] = 1 -  ptData$branch_relative_Pseudotime[ ptData$State == aState ]
	}
}

# 20210308
# Relabel branches 
ptData2 = ptData
    
# Must manually annotate states/branches so the branch labels match across modalities
# 20210320
if (whichChrom == "genome" && DEqval==10^-11) {
	ptData2$branch_labels = "NA"
	ptData2$branch_labels[ ptData2$State == 2 ] = "early"
	ptData2$branch_labels[ ptData2$State == 1 ] = "intermediate"
	ptData2$branch_labels[ ptData2$State == 3 ] = "late"
	branch_distal_limits = c( 0.91, 0.83, 0.73 )
	names(branch_distal_limits) = c( "early", "intermediate", "late")
} else
	if (whichChrom == "autosomes") {
		ptData2$branch_labels = "NA"
		ptData2$branch_labels[ ptData2$State == 2 ] = "early"
		ptData2$branch_labels[ ptData2$State == 1 ] = "intermediate"
		ptData2$branch_labels[ ptData2$State == 3 ] = "late"
		branch_distal_limits = c( 0.91, 0.75, 0.7 )
		names(branch_distal_limits) = c( "early", "intermediate", "late")
	} else if (whichChrom == "chrX" &&  DEqval=10^-5) {
		ptData2$branch_labels = "NA"
		ptData2$branch_labels[ ptData2$State == 4 ] = "early"
		ptData2$branch_labels[ ptData2$State == 3 ] = "intermediate"
		ptData2$branch_labels[ ptData2$State == 1 ] = "late"
		ptData2$branch_labels[ ptData2$State == 5 ] = "mystery"
		branch_distal_limits = c( 0.8, 0.8, 0.8, 0.25 )
		names(branch_distal_limits) = c( "early", "intermediate", "late", "mystery")
	} else if (whichChrom == "chrX" && DEqval==10^-4.5) {
		ptData2$branch_labels = "NA"
		ptData2$branch_labels[ ptData2$State == 1 ] = "early"
		ptData2$branch_labels[ ptData2$State == 2 ] = "intermediate"
		ptData2$branch_labels[ ptData2$State == 3 ] = "late"
		branch_distal_limits = c( 0.84, 0.77, 0.85 )
		names(branch_distal_limits) = c( "early", "intermediate", "late")
	}
}
    
    
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".pseudotime_CDFs_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

branch_distal_limits_df = as.data.frame(branch_distal_limits)
branch_distal_limits_df$branch_label = rownames(as.data.frame(branch_distal_limits_df))

print( ggplot(ptData2, aes(branch_relative_Pseudotime, color = branch_labels)) + 
		   stat_ecdf(geom = "step") + 
		   labs(title = "sci-Hi-C pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") +
		   geom_vline(data = branch_distal_limits_df, aes(xintercept = branch_distal_limits, 
														  color = branch_label), linetype = "dashed")) 

print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()


tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".pseudotime_per_branch.tsv", sep=""))
write.table(ptData2, 
			tablefilename,sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



##################################################
# 20200928
# Get and plot cells at distal ends of branches.

ptData2$distal=FALSE
for ( aBranch in sort( unique( ptData2$branch_labels ) ) ) {
	branch_distal_cell_boolidx = ptData2$branch_relative_Pseudotime[ ptData2$branch_labels == aBranch ] >= branch_distal_limits[aBranch]
	ptData2$distal[ ptData2$branch_labels == aBranch ][ branch_distal_cell_boolidx ] = TRUE
}

cds.distal_cells = cds
pData(cds.distal_cells)$distal_cells = ptData2$distal

# Show distal cells along trajectory:
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".plotTrajectories.distal_cells.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
# print( plot_cell_trajectory(cds.distal_cells, color_by = "timepoint") + 
#            geom_point(size=2, aes(color = timepoint)) + 
#            scale_color_manual(values = tpCols) + 
#            facet_wrap(~distal_cells, nrow = 1) )
print( plot_cell_trajectory(cds.distal_cells, color_by = "distal_cells") + 
		   geom_point(size=2, aes(color = distal_cells)) + 
		   scale_color_manual(values = c("grey", "red")) )
dev.off()

# Table of cells distal classification
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".pseudotime_per_branch.distal_cells.tsv", sep=""))
write.table(ptData2, 
			tablefilename, sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)

# Table of cells per branch - proportions
tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											   ".minExprThreshold", minExprThreshold, 
											   ".numCellsExpressed", numCellsExpressed,
											   ".downSampleProp", dsProp, 
											   ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											   ".pseudotime_per_branch.distal_cells.proportions.tsv", sep=""))
table_res <- table( ptData2$distal, ptData2$branch_labels )
table_res <- table_res[order(table_res[,1], decreasing=TRUE), ]
table_res <- t( t( table_res ) / colSums(table_res) )
table_res <- rbind( table_res, table( ptData2$branch_labels ))
if(rownames(table_res)[1]) {
	rownames(table_res) = c("distal", "non-distal", "cell_count")
} else {
	rownames(table_res) = c("non-distal", "distal", "cell_count")
}
table_res
write.table(table_res, 
			tablefilename,sep="\t", row.names = TRUE, col.names = NA, quote=FALSE)



##################################################
# 20200928
# Coverage per time point

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".total_coverage_distributions_per_timepoint.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(pData(cds), aes(totalCoverage, colour = timepoint)) + scale_color_manual(values = tpCols) + 
	stat_ecdf(geom = "step") + 
	labs(title = "sci-Hi-C coverage\n", y = "Proportion", color = "time point\n") 

dev.off()



##################################################
# 20210308
# Coverage per branch

PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".total_coverage_distributions_per_branch.pdf", sep=""))
pdf(PDFfilename, height=5, width=6)

ggplot(ptData2, aes(totalCoverage, colour = branch_labels)) + 
	stat_ecdf(geom = "step") + 
	labs(title = "sci-Hi-C coverage\n", y = "Proportion", color = "branches\n") 

dev.off()
    
    
    
##################################################
# 20201007
# Which are the most differential bins @ the distal ends of branches.    
# T-test each bin between branches.

# exprs(cds)[1:5,1:2]
# 
# dim(exprs(cds))
# [1] 5434 2056
# 5434/143
# [1] 38

chroms_to_test = c("chr1", "chr2", "chr3", "chrX")

theBranchLabels = sort( unique( ptData2$branch_labels ) )

binDist = c( 2^seq(10, 27.750, 0.125))
length(binDist)
# 143
# NOTE: THESE HAVE NOT BEEN REBINNED BY 10 BINS

prettyBPs <- function(x) {
	cx = paste(as.character(round(x,1)), "bp", sep="")
	cx[x>=10^3] = paste(as.character(round(x[x>=10^3]/10^3,2)), "kb", sep="")
	cx[x>=10^6] = paste(as.character(round(x[x>=10^6]/10^6,1)), "Mb", sep="")
	return(cx)
}

binDistDesc = prettyBPs(binDist)
if (length(binDistDesc) > 100) {
	binDistDesc[seq(2,length(binDistDesc),2)] = ""
}
binDistDesc

ttest_result_list = list()

# 20210428
for (testDistal in c("distal", "proximal", "all")) {
	
	ttest_result_list[[testDistal]] = list()
	
	for ( atestchrom in chroms_to_test) {
		# atestchrom = "chr1"
		
		print( atestchrom )
		
		ttest_table = c()
		
		for ( aBranch in theBranchLabels ) {
			
			for ( bBranch in theBranchLabels ) {
				if ( bBranch <= aBranch ) next
				print( paste( aBranch, bBranch, sep=" vs " ) )
				
				if (testDistal == "distal") {
					aBranch_distal_cell_boolidx = ptData2$distal & ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = ptData2$distal & ptData2$branch_labels == bBranch 
				} else  if (testDistal == "proximal") {
					aBranch_distal_cell_boolidx = (! ptData2$distal) & ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = (! ptData2$distal) & ptData2$branch_labels == bBranch 
				} else {
					aBranch_distal_cell_boolidx = ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = ptData2$branch_labels == bBranch 
				}
				
				atestchrom_text = paste0( "_", atestchrom, ":")
				chrom_bin_boolIdx = grepl( atestchrom_text, rownames(exprs(cds)) )
				sum(chrom_bin_boolIdx)
				
				ttests_per_bin = c()
				
				for ( bin_idx in 1:sum(chrom_bin_boolIdx) ) {
					
					aBranch_bin_scores = as.matrix( exprs(cds)[ chrom_bin_boolIdx, aBranch_distal_cell_boolidx ] )[ bin_idx, ]
					bBranch_bin_scores = as.matrix( exprs(cds)[ chrom_bin_boolIdx, bBranch_distal_cell_boolidx ] )[ bin_idx, ]
					
					t_result = t.test(aBranch_bin_scores, bBranch_bin_scores)
					
					ttests_per_bin = c(ttests_per_bin, t_result$p.value)
				}
				
				# ttests_per_bin = as.matrix( p.adjust(ttests_per_bin, method="BH") )
				ttests_per_bin = as.matrix( ttests_per_bin )
				
				rownames(ttests_per_bin) = rownames(exprs(cds)[ chrom_bin_boolIdx, aBranch_distal_cell_boolidx ])
				colnames(ttests_per_bin) =  paste( aBranch, bBranch, sep="-vs-" )
				
				ttest_table = cbind(ttest_table, ttests_per_bin)
			}
		}
		
		ttest_result_list[[testDistal]][[atestchrom]] = ttest_table
	}
	
} #testDistal
	

##################
# Bar plot results

# 20210428
for (testDistal in c("distal", "proximal", "all")) {
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".differential_bins_between_branches", 
												 ifelse(testDistal!="distal", paste0(".", testDistal), ""), 
												 ".pdf", sep=""))
	pdf(PDFfilename, height=3*3, width=4*2)
	par(mfrow=c(3,2))
	for ( atestchrom in names(ttest_result_list[[testDistal]]) ) {
		for ( acomparison in colnames(ttest_result_list[[testDistal]][[ atestchrom ]]) ) {
			# acomparison = "early-vs-intermediate"
			
			for (allele in c("ref", "alt") ) {
				# allele = "ref"
				row_boolidx = grepl( allele, names( ttest_result_list[[testDistal]][[ atestchrom ]][, acomparison] ) )
				barplot( -log10( p.adjust( ttest_result_list[[testDistal]][[ atestchrom ]][row_boolidx, acomparison], method="BH" ) ), 
						 names.arg=binDistDesc, las=2, cex.names=0.5,
						 ylab="Adjuste p-value (-log10)", 
						 main=paste( atestchrom, acomparison, allele, sep=" " ) )
			}
		}
	}
	dev.off()

} # testDistal



##################################################
# 20201014
# Plot log2 fold change in average CDP scores between branches of distal cells.

chroms_to_test = c("chr1", "chr2", "chr3", "chrX")

theBranchLabels = sort( unique( ptData2$branch_labels ) )

binDist = c( 2^seq(10, 27.750, 0.125))
length(binDist)
# 143
# NOTE: THESE HAVE NOT BEEN REBINNED BY 10 BINS

# 20210320
if(dropfirstbins) {
	binDist = c( 2^seq(10 + (dropNbins * 0.125), 27.750, 0.125))
	length(binDist)
	# 140 - drop first 3 bins
}

prettyBPs <- function(x) {
	cx = paste(as.character(round(x,1)), "bp", sep="")
	cx[x>=10^3] = paste(as.character(round(x[x>=10^3]/10^3,2)), "kb", sep="")
	cx[x>=10^6] = paste(as.character(round(x[x>=10^6]/10^6,1)), "Mb", sep="")
	return(cx)
}

binDistDesc = prettyBPs(binDist)
if (length(binDistDesc) > 100) {
	binDistDesc[seq(2,length(binDistDesc),2)] = ""
}
binDistDesc

log2ratio_list = list()

# 20210428
for (testDistal in c("distal", "proximal", "all")) {
	
	log2ratio_list[[testDistal]] = list()
	
	for ( atestchrom in chroms_to_test) {
		# atestchrom = "chr1"
		
		print( atestchrom )
		
		log2ratio_table = c()
		
		for ( aBranch in theBranchLabels ) {
			
			for ( bBranch in theBranchLabels ) {
				if ( bBranch <= aBranch ) next
				print( paste( aBranch, bBranch, sep=" vs " ) )
				
				if (testDistal == "distal") {
					aBranch_distal_cell_boolidx = ptData2$distal & ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = ptData2$distal & ptData2$branch_labels == bBranch 
				} else  if (testDistal == "proximal") {
					aBranch_distal_cell_boolidx = (! ptData2$distal) & ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = (! ptData2$distal) & ptData2$branch_labels == bBranch 
				} else {
					aBranch_distal_cell_boolidx = ptData2$branch_labels == aBranch 
					bBranch_distal_cell_boolidx = ptData2$branch_labels == bBranch 
				}
				
				atestchrom_text = paste0( "_", atestchrom, ":")
				chrom_bin_boolIdx = grepl( atestchrom_text, rownames(exprs(cds)) )
				sum(chrom_bin_boolIdx)
				
				aBranch_bin_meanscores = apply(as.matrix( exprs(cds)[ chrom_bin_boolIdx, aBranch_distal_cell_boolidx ] ), 1, mean, na.rm=TRUE )
				bBranch_bin_meanscores = apply(as.matrix( exprs(cds)[ chrom_bin_boolIdx, bBranch_distal_cell_boolidx ] ), 1, mean, na.rm=TRUE )
				
				log2ratio_branch_bin_meanscores = as.matrix( log2( ( (aBranch_bin_meanscores + 1) / (bBranch_bin_meanscores + 1) ) ) )
				colnames(log2ratio_branch_bin_meanscores) =  paste( aBranch, bBranch, sep="-vs-" )
				
				log2ratio_table = cbind(log2ratio_table, log2ratio_branch_bin_meanscores)
				
			}
		}
		
		log2ratio_list[[testDistal]][[atestchrom]] = log2ratio_table
	}
	
} # testDistal



##################################################
# Plot results

ymax = lapply( lapply( log2ratio_list, abs ), max, na.rm=TRUE )
if ( cellTypesIncluded == "F121only" ) qval_threshold = 10^-6
if ( cellTypesIncluded == "F121only.interphase" ) qval_threshold = 10^-3
highlight_color = rgb(0.85, 0.85, 0.25, 1/4)

# 20210428
for (testDistal in c("distal", "proximal", "all")) {
	
	####
	# V1
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".log2ratio_bins_between_branches",
												 ifelse(testDistal!="distal", paste0(".", testDistal), ""), 
												 ".pdf", sep=""))
	pdf(PDFfilename, height=3*3, width=4*2)
	par(mfrow=c(3,2))
	for ( atestchrom in names(log2ratio_list[[testDistal]]) ) {
		for ( acomparison in colnames(log2ratio_list[[testDistal]][[atestchrom]]) ) {
			# acomparison = "early-vs-intermediate"
			
			for (allele in c("ref", "alt") ) {
				# allele = "ref"
				row_boolidx = grepl( allele, names( log2ratio_list[[testDistal]][[atestchrom]][, acomparison] ) )
				plot( log2ratio_list[[testDistal]][[atestchrom]][row_boolidx, acomparison], 
					  col = alleleCols[ allele ], lwd=2, 
					  type="l", ylim=c(-ymax[[ atestchrom ]], +ymax[[ atestchrom ]]), 
					  xlab="bin", ylab="log2 ratio", 
					  main=paste( atestchrom, acomparison, allele, "\nadjusted p-val <", qval_threshold, sep=" " ) )
				
				# Add shading for significant regions
				significant_difference_idx = which( p.adjust( ttest_result_list[[testDistal]][[ atestchrom ]][row_boolidx, acomparison], method="BH" ) < qval_threshold )
				for(i in significant_difference_idx ) {
					print( i )
					rect( i, -ymax[[ atestchrom ]], 
						  i+1, +ymax[[ atestchrom ]], 
						  col = highlight_color, border=NA )
					i = i+1
				}
				
				# Re-overplot lines
				lines( log2ratio_list[[testDistal]][[atestchrom]][row_boolidx, acomparison], 
					   col = alleleCols[ allele ], lwd=2 )
				
				abline(h=0, lty="dashed", col="black", lwd=1)
			}
		}
	}
	dev.off()
	
	####
	# V2
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".log2ratio_bins_between_branches",
												 ifelse(testDistal!="distal", paste0(".", testDistal), ""), 
												 ".V2.pdf", sep=""))
	pdf(PDFfilename, height=3*3, width=4*2)
	par(mfrow=c(3,2))
	for ( atestchrom in names(log2ratio_list[[testDistal]]) ) {
		for ( acomparison in colnames(log2ratio_list[[testDistal]][[atestchrom]]) ) {
			# acomparison = "early-vs-intermediate"
			
			for (allele in c("ref", "alt") ) {
				# allele = "ref"
				row_boolidx = grepl( allele, names( log2ratio_list[[testDistal]][[atestchrom]][, acomparison] ) )
				plot( -log2ratio_list[[testDistal]][[atestchrom]][row_boolidx, acomparison], 
					  col = alleleCols[ allele ], lwd=2, 
					  type="l", ylim=c(-ymax[[ atestchrom ]], +ymax[[ atestchrom ]]), 
					  xlab="bin", ylab="log2 ratio", 
					  main=paste( atestchrom, acomparison, allele, "\nadjusted p-val <", qval_threshold, sep=" " ) )
				
				# Add shading for significant regions
				significant_difference_idx = which( p.adjust( ttest_result_list[[testDistal]][[ atestchrom ]][row_boolidx, acomparison], method="BH" ) < qval_threshold )
				for(i in significant_difference_idx ) {
					print( i )
					rect( i, -ymax[[ atestchrom ]], 
						  i+1, +ymax[[ atestchrom ]], 
						  col = highlight_color, border=NA )
					i = i+1
				}
				
				# Re-overplot lines
				lines( -log2ratio_list[[testDistal]][[atestchrom]][row_boolidx, acomparison], 
					   col = alleleCols[ allele ], lwd=2 )
				
				abline(h=0, lty="dashed", col="black", lwd=1)
			}
		}
	}
	dev.off()
	
} # testDistal











##################################################
# 20210503
# Projections of cells V2
# For all cells and distal cells

dim(cds)
# Features  Samples 
# 5720     2056 

for ( whichCells in c( "all", "distal", "proximal" ) ) {
	
	if ( whichCells == "all" ) {
		theTpCols = tpCols[pData(cds)$timepoint]
	} else if ( whichCells == "distal" ) {
		theTpCols = tpCols[pData(cds[, pData(cds.distal_cells)$distal_cells ])$timepoint]
	} else if ( whichCells == "proximal" ) {
		theTpCols = tpCols[pData(cds[, ! pData(cds.distal_cells)$distal_cells ])$timepoint]
	}
	
	if ( whichCells == "all" ) {
		cds_matrix = t( as.matrix( exprs( cds ) ) )
	} else if ( whichCells == "distal" ) {
		cds_matrix = t( as.matrix( exprs( cds[, pData(cds.distal_cells)$distal_cells ] ) ) )
	} else if ( whichCells == "proximal" ) {
		cds_matrix = t( as.matrix( exprs( cds[, ! pData(cds.distal_cells)$distal_cells ] ) ) )
	}
	
	dim(cds_matrix)
	
	# remove zero variance rows/columns
	cds_matrix = cds_matrix[ apply(cds_matrix, 1, var) != 0, apply(cds_matrix, 2, var) != 0 ]
	
	if ( whichCells == "all" ) {
		cds_matrix.DE = t( as.matrix( exprs(cds[ordering_timepoint_DE_loci, ]) ) )
	} else if ( whichCells == "distal" ) {
		cds_matrix.DE = t( as.matrix( exprs(cds[ordering_timepoint_DE_loci, pData(cds.distal_cells)$distal_cells ]) ) )
	} else if ( whichCells == "proximal" ) {
		cds_matrix.DE = t( as.matrix( exprs(cds[ordering_timepoint_DE_loci, !pData(cds.distal_cells)$distal_cells ]) ) )
	}        
	dim(cds_matrix.DE)
	
	cor_cds = cor( t( cds_matrix ) )
	cor_cds.DE = cor( t( cds_matrix.DE ) )
	
	dim(cor_cds)
	dim(cor_cds.DE)
	# 1253 1253
	
	distal_cds_pca = list()
	distal_cds_pca[[ "PCA" ]] = prcomp(cds_matrix, scale=TRUE)
	distal_cds_pca[[ "PCA_DE" ]] = prcomp(cds_matrix.DE, scale=TRUE)
	distal_cds_pca[[ "PCA_CORR" ]] = prcomp(cor_cds, scale=TRUE)
	distal_cds_pca[[ "PCA_CORR_DE" ]] = prcomp(cor_cds.DE, scale=TRUE)
	
	require(umap)
	distal_cds_umap = list()
	distal_cds_umap[[ "UMAP" ]] = umap(cds_matrix, scale=TRUE)
	distal_cds_umap[[ "UMAP_DE" ]] = umap(cds_matrix.DE, scale=TRUE)
	distal_cds_umap[[ "UMAP_CORR" ]] = umap(cor_cds, scale=TRUE)
	distal_cds_umap[[ "UMAP_CORR_DE" ]] = umap(cor_cds.DE, scale=TRUE)
	
	# PCA centroids per cell type
	pca_centroids = list()
	for ( projection_type in names(distal_cds_pca) ){
		projection_data = distal_cds_pca[[ projection_type ]]
		# PCA centroids per cell type
		the_centroids = c()
		for ( cell_type in names(tpCols) ) {
			# cell_type = "F121_d0"
			cds_boolidx = grepl( cell_type, rownames(projection_data$x) )
			the_centroid = apply( projection_data$x[cds_boolidx,c(1,2)], 2, mean, na.rm=TRUE )
			print( paste0( cell_type, ": " ) )
			print( the_centroid )
			the_centroids = rbind(the_centroids, the_centroid)
		}
		rownames(the_centroids) = names(tpCols)
		pca_centroids[[ projection_type ]] = the_centroids
	}
	
	# UMAP centroids per cell type
	umap_centroids = list()
	for ( projection_type in names(distal_cds_umap) ){
		projection_data = distal_cds_umap[[ projection_type ]]
		# UMAP centroids per cell type
		the_centroids = c()
		for ( cell_type in names(tpCols) ) {
			# cell_type = "F121_d0"
			cds_boolidx = grepl( cell_type, rownames(projection_data$layout) )
			the_centroid = apply( projection_data$layout[cds_boolidx,c(1,2)], 2, mean, na.rm=TRUE )
			print( paste0( cell_type, ": " ) )
			print( the_centroid )
			the_centroids = rbind(the_centroids, the_centroid)
		}
		rownames(the_centroids) = names(tpCols)
		umap_centroids[[ projection_type ]] = the_centroids
	}
	
	
	
	###################
	# projection  plots
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".CDP_projections.", whichCells, ".pdf", sep=""))
	pdf(PDFfilename, height=8, width=15)
	par(mfrow=c(2,4), oma=c(0,0,2,0))
	
	for ( projection_type in names(distal_cds_pca) ) {
		plot(distal_cds_pca[[projection_type]]$x[,1], 
			 distal_cds_pca[[projection_type]]$x[,2], 
			 col=theTpCols, cex=0.75, pch=19,
			 main=paste0(projection_type) , xlab="PC1", ylab="PC2")
		points( pca_centroids[[ projection_type ]][, 1], pca_centroids[[ projection_type ]][, 2], col="yellow", cex=3, pch=19 )
		points( pca_centroids[[ projection_type ]][, 1], pca_centroids[[ projection_type ]][, 2], col=tpCols, cex=2, pch=19 )
	}
	
	for ( projection_type in names(distal_cds_umap) ) {
		plot(distal_cds_umap[[projection_type]]$layout[,1], 
			 distal_cds_umap[[projection_type]]$layout[,2], 
			 col=theTpCols, cex=0.75, pch=19,
			 main=paste0( projection_type) , xlab="UMAP1", ylab="UMAP2")
		points( umap_centroids[[ projection_type ]][, 1], umap_centroids[[ projection_type ]][, 2], col="yellow", cex=3, pch=19 )
		points( umap_centroids[[ projection_type ]][, 1], umap_centroids[[ projection_type ]][, 2], col=tpCols, cex=2, pch=19 )
	}
	legend( "topright", legend=names(tpCols), fill=tpCols)
	
	title(main=paste0("CDP projections - ", whichCells, " cells"), outer=TRUE)
	dev.off()
	
	
	###################
	# projection  centroid dendograms
	
	PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
												 ".minExprThreshold", minExprThreshold, 
												 ".numCellsExpressed", numCellsExpressed,
												 ".downSampleProp", dsProp, 
												 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
												 ".CDP_projection_centroid_dendrograms.", whichCells, ".pdf", sep=""))
	pdf(PDFfilename, height=8, width=15)
	par(mfrow=c(2,4), oma=c(0,0,2,0))
	
	for ( projection_type in names(pca_centroids) ) {
		plot( hclust( dist( pca_centroids[[ projection_type ]] )), main=projection_type )
	}
	
	for ( projection_type in names(umap_centroids) ) {
		plot( hclust( dist( umap_centroids[[ projection_type ]] )), main=projection_type )
	}
	
	title(main=paste0("CDP projections centroid dendrograms - ", whichCells, " cells"), outer=TRUE)
	dev.off()
	
} # "all", "distal", "proximal"



##################################################
# 20210503
# Projections of cells at distal ends of branches
    
dim(cds)
# Features  Samples 
# 5720     2056 

distal_matrix = t( as.matrix( exprs( cds[, pData(cds.distal_cells)$distal_cells ] ) ) )
dim(distal_matrix)
# [1] 5720 1253
sum(pData(cds.distal_cells)$distal_cells)
# 1253

# remove zero variance rows/columns
distal_matrix = distal_matrix[ apply(distal_matrix, 1, var) != 0, apply(distal_matrix, 2, var) != 0 ]

distal_matrix.DE = t( as.matrix( exprs(cds[ordering_timepoint_DE_loci, pData(cds.distal_cells)$distal_cells ]) ) )
dim(distal_matrix.DE)
#  2040 1253

cor_cds = cor( t( distal_matrix ) )
cor_cds.DE = cor( t( distal_matrix.DE ) )

dim(cor_cds)
dim(cor_cds.DE)
# 1253 1253

distal_cds_pca = list()
distal_cds_pca[[ "PCA" ]] = prcomp(distal_matrix, scale=TRUE)
distal_cds_pca[[ "PCA_DE" ]] = prcomp(distal_matrix.DE, scale=TRUE)
distal_cds_pca[[ "PCA_CORR" ]] = prcomp(cor_cds, scale=TRUE)
distal_cds_pca[[ "PCA_CORR_DE" ]] = prcomp(cor_cds.DE, scale=TRUE)

require(umap)
distal_cds_umap = list()
distal_cds_umap[[ "UMAP" ]] = umap(distal_matrix, scale=TRUE)
distal_cds_umap[[ "UMAP_DE" ]] = umap(distal_matrix.DE, scale=TRUE)
distal_cds_umap[[ "UMAP_CORR" ]] = umap(cor_cds, scale=TRUE)
distal_cds_umap[[ "UMAP_CORR_DE" ]] = umap(cor_cds.DE, scale=TRUE)

# PCA centroids per cell type
pca_centroids = list()
for ( projection_type in names(distal_cds_pca) ){
	projection_data = distal_cds_pca[[ projection_type ]]
	# PCA centroids per cell type
	the_centroids = c()
	for ( cell_type in names(tpCols) ) {
		# cell_type = "F121_d0"
		cds_boolidx = grepl( cell_type, rownames(projection_data$x) )
		the_centroid = apply( projection_data$x[cds_boolidx,c(1,2)], 2, mean, na.rm=TRUE )
		print( paste0( cell_type, ": " ) )
		print( the_centroid )
		the_centroids = rbind(the_centroids, the_centroid)
	}
	rownames(the_centroids) = names(tpCols)
	pca_centroids[[ projection_type ]] = the_centroids
}

# UMAP centroids per cell type
umap_centroids = list()
for ( projection_type in names(distal_cds_umap) ){
	projection_data = distal_cds_umap[[ projection_type ]]
	# PCA centroids per cell type
	the_centroids = c()
	for ( cell_type in names(tpCols) ) {
		# cell_type = "F121_d0"
		cds_boolidx = grepl( cell_type, rownames(projection_data$layout) )
		the_centroid = apply( projection_data$x[cds_boolidx,c(1,2)], 2, mean, na.rm=TRUE )
		print( paste0( cell_type, ": " ) )
		print( the_centroid )
		the_centroids = rbind(the_centroids, the_centroid)
	}
	rownames(the_centroids) = names(tpCols)
	umap_centroids[[ projection_type ]] = the_centroids
}

# plots
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
											 ".minExprThreshold", minExprThreshold, 
											 ".numCellsExpressed", numCellsExpressed,
											 ".downSampleProp", dsProp, 
											 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
											 ".distal_cell_CDP_projections.pdf", sep=""))
pdf(PDFfilename, height=8, width=15)

par(mfrow=c(2,4))

for ( projection_type in names(distal_cds_pca) ) {
	plot(distal_cds_pca[[projection_type]]$x[,1], 
		 distal_cds_pca[[projection_type]]$x[,2], 
		 col=tpCols[pData(cds[, pData(cds.distal_cells)$distal_cells ])$timepoint], cex=0.75, pch=19,
		 main=paste0(projection_type) , xlab="PC1", ylab="PC2")
	points( pca_centroids[[ projection_type ]][, 1], pca_centroids[[ projection_type ]][, 2], col="yellow", cex=3, pch=19 )
	points( pca_centroids[[ projection_type ]][, 1], pca_centroids[[ projection_type ]][, 2], col=tpCols, cex=2, pch=19 )
	legend( "bottomright", legend=names(tpCols), fill=tpCols)
}

for ( projection_type in names(distal_cds_umap) ) {
	plot(distal_cds_umap[[projection_type]]$layout[,1], 
		 distal_cds_umap[[projection_type]]$layout[,2], 
		 col=tpCols[pData(cds[, pData(cds.distal_cells)$distal_cells ])$timepoint], cex=0.75, pch=19,
		 main=paste0( projection_type) , xlab="UMAP1", ylab="UMAP2")
	points( umap_centroids[[ projection_type ]][, 1], umap_centroids[[ projection_type ]][, 2], col="yellow", cex=3, pch=19 )
	points( umap_centroids[[ projection_type ]][, 1], umap_centroids[[ projection_type ]][, 2], col=tpCols, cex=2, pch=19 )
	legend( "bottomright", legend=names(tpCols), fill=tpCols)
	
}
legend( "topright", legend=names(tpCols), fill=tpCols)

title(main="CDP projections -- distal cells", outer=TRUE)

dev.off()
    
    
    
##################################################
# Save object
# set_lower_bound = 0; thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis.semisupervised", ifelse(set_lower_bound >0, paste0( ".set_lower_bound", set_lower_bound ), "") ) )
# minExprThreshold=0.1; numCellsExpressed=10; dsProp=1; DEqval=0.01; # DEqval=10^-11; # DEqval=10^-5
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genome.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed),
                                                 ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                                 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "") 
) )

# save.image(file =  RdataFileName )
# load(file = RdataFileName )
