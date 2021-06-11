##################################################
##################################################
##################################################
#
# Giancarlo Bonora
#
# 20200817 allelic sci-RNA-seq analysis using Monocle2
#          using concatened chromosome-specific vectors
#          for F121 cells only
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

# 20200720
# min_counts_per_chromosome=0
min_counts_per_chromosome=10

# 20200803
whichChrom="genome"

# 20200916
excludeNPCs=FALSE
# excludeNPCs=TRUE

projectName = "allelic_sci-RNA-seq"

inputDir = paste("../data/allelic_sci-RNA-seq_catenated_count_matrices/")
    
mainOutputDir = paste("../results/allelic_sci-RNA-seq_monocle2_catenated_allelic_vectors/")
if (min_counts_per_chromosome>0) mainOutputDir = paste0( mainOutputDir, ".min_counts_per_chromosome", min_counts_per_chromosome )
if (whichChrom != "genome") mainOutputDir = paste0( mainOutputDir, ".", whichChrom)

# 20200916
if (excludeNPCs) mainOutputDir = paste0( mainOutputDir, ".noNPCs")

if ( ! file.exists( mainOutputDir ) ) {
  if ( ! dir.create( mainOutputDir, showWarnings=FALSE, recursive=TRUE )  ) mainOutputDir = getwd()
}
setwd(mainOutputDir)
getwd()

tpCols=c("red", "orange", "green", "blue", "black")
names(tpCols) = c("d0", "d3", "d7", "d11", "NPC")
 
pchVal=20
pchcexVal=1.25










##################################################
##################################################
##################################################
# Load counts data
# https://stackoverflow.com/questions/53967385/r-libs-user-ignored-by-r


if (whichChrom == "genome") {
    
    count.table.OG = c()
    for (theChrom in c( paste0( "chr", 1:19 ), "chrX" ) ) { 
        inputFileName = paste0("allelic_sci-RNA-seq.UMImatrix.filtered.minExpr0.1.minCells5.F121.",
                               theChrom, 
                               ".stringentlyFiltered",
                               ifelse(min_counts_per_chromosome>0, paste0(".min_counts_per_chromosome", min_counts_per_chromosome), ""), 
                               ".tsv.gz" )
        countsFile = file.path(inputDir, inputFileName)
        count.table.OG = rbind( count.table.OG,
                                as.data.frame(fread(countsFile, header=TRUE, sep="\t") ) )
    }
    
} 

rownames(count.table.OG) = count.table.OG[,1]
count.table.OG = count.table.OG[,-1]

dim(count.table.OG)
#  20966  1048 = genome  - .min_counts_per_chromosome10










##################################################
##################################################
##################################################
# No need to transform table; already a count matrix 

uniqBarcodes.OG =  unique(colnames(count.table.OG))
length(uniqBarcodes.OG)
# 1048

uniqLoci.OG = unique(rownames(count.table.OG))
length(uniqLoci.OG)
# 20966

# Rename loci
locIDs.OG <- rownames(count.table.OG) 

# cell IDs
cell_list.OG <- colnames(count.table.OG)










##################################################
##################################################
##################################################
# Generate CDS data

sapply(count.table.OG, class)
# as.matrix(count.matrix.OG[1:3,1:3])

# require(Matrix)
count.matrix.OG.sparse = Matrix::Matrix(as.matrix(count.table.OG), sparse=TRUE)

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
# 1048    4

dim(fData(cds.OG))
# 20966   2

dim(Biobase::exprs(cds.OG))
#   20966  1048



##################################################
# Add metadata

# time-points
timepoints = as.matrix( sub( "^.*_", "", sub( "\\..*$", "", rownames(pData(cds.OG)) ) ) )
rownames(timepoints) <- rownames(pData(cds.OG))
colnames(timepoints) <- "timepoints"
table(timepoints)
# timepoints
# d0 d11  d3  d7 NPC 
# 489 145 127  78 209 

all( colnames(cds.OG) == rownames(timepoints) )
# TRUE

pData(cds.OG)$timepoint <- factor( timepoints, levels=c("d0", "d3", "d7", "d11", "NPC") )










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
# 11,157

print(head(pData(cds.OG)))










##################################################
##################################################
##################################################
# Filter out extreme cells

# Total expression per barcode
pData(cds.OG)$totalCoverage <- Matrix::colSums(Biobase::exprs(cds.OG))       # !NB!

summary(pData(cds.OG)$totalCoverage )

# genome
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 305.0   441.0   579.5   649.7   800.8  1472.0 

set_lower_bound = 0 

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
# genome
# Features  Samples 
# 20966     1048 

dim(cds)
# genome
# Features  Samples 
# 20966     1022

# 20200916 Genome-wide trajectory analysis *** EXCLUDING NPCs ***
if (excludeNPCs) {
    cds = cds[,pData(cds)$timepoint != "NPC"]
    dim(cds)
    # Features  Samples 
    # 20966      824 
}
    










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
# Use the known cell state!

# set_lower_bound=0
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

minExprThreshold=0.1
numCellsExpressed=10

# 20200518
# Try downsampling:
# dsProp=NA
dsProp=1
# dsProp=0.5

cds <- detectGenes(cds, min_expr = minExprThreshold)
coveredLoci <-  row.names(subset(fData(cds), num_cells_expressed >= numCellsExpressed))
length(coveredLoci)
# 10868 - genome

dim(cds)
# Features  Samples 
#    20966     1022 

set.seed(123)
coveredLoci = sample( coveredLoci, size=floor(dsProp * length(coveredLoci)) )
length(coveredLoci)

# diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~Media")
timepoint_DE_loci <- differentialGeneTest(cds[coveredLoci,],
                                          fullModelFormulaStr = '~timepoint',
                                          cores = 4)

DEqval=0.01

ordering_timepoint_DE_loci<- row.names(subset(timepoint_DE_loci, qval < DEqval))
length(ordering_timepoint_DE_loci)

###############
# Genome
# 4230 # 0.01

###############
# Genome - No NPCs
# 2261 # 0.01



##################################################
# 20210331
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
GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$timepoint)[,"d0"]
        return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
cds <- orderCells(cds, root_state = GM_state(cds))



##################################################
# Tabulate results

table_res <- table(pData(cds)$State, pData(cds)$timepoint)
table_res[order(table_res[,1], decreasing=TRUE), ]

###############
# Genome
# DEqval=0.01
# d0  d3  d7 d11 NPC
# 1 477  82   1   3   0
# 3   3  43  76 131   0
# 2   1   1   0   6 198


###############
# Genome - No NPCs
# DEqval=0.01
# d0  d3  d7 d11 NPC
# 1 480  42   1   2   0
# 3   1  57   3   8   0
# 2   0  27  73 130   0



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
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".byTimePoint.plotTrajectories.pdf", sep=""))
pdf(PDFfilename, height=10, width=10)
# The trajectory has a tree-like structure. We can see that the cells collected at time zero are located 
# near one of the tips of the tree, while the others are distributed amongst the two "branches". 
# Monocle doesn't know a priori which of the trajectory of the tree to call the "beginning",
# so we often have to call orderCells again using the root_state argument to specify the beginning.
# First, we plot the trajectory, this time coloring the cells by "State":
# print( plot_cell_trajectory(cds, color_by = "Cluster") + geom_point(size=2, aes(color = Cluster)) )
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) )
print( plot_cell_trajectory(cds, color_by = "Pseudotime") + geom_point(size=2, aes(color = Pseudotime)) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) )
dev.off()



##################################################
# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space - faceted
PDFfilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                             ".minExprThreshold", minExprThreshold, 
                                             ".numCellsExpressed", numCellsExpressed,
                                             ".downSampleProp", dsProp, 
                                             ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                             ".byTimePoint.plotTrajectories.faceted.pdf", sep=""))
pdf(PDFfilename, height=10, width=20)
# If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. 
# Sometimes it can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located:
print( plot_cell_trajectory(cds, color_by = "State") + geom_point(size=2, aes(color = State)) + facet_wrap(~State, nrow = 1) )
print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )
dev.off()

# table_res[order(table_res[,1], decreasing=TRUE), ]



##################################################
# 20200928
# CDFs of pseudotime density

zero2one <-function(x) { return( (x-min(x))/(max(x)-min(x)) ) }

ptData = pData(cds)[, c("cell", "timepoint", "Pseudotime", "State") ]

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
    labs(title = "sci-RNA pseudotime CDFs\n", x = "renormalized pseudotime along branch", y = "Proportion", color = "Branch\n") 

print( plot_cell_trajectory(cds, color_by = "timepoint") + geom_point(size=2, aes(color = timepoint)) + scale_color_manual(values = tpCols) + facet_wrap(~State, nrow = 1) )

dev.off()



tablefilename = file.path(thisOutputDir, paste(projectName, ".monocle2results.trajectoryAnalysis", 
                                               ".minExprThreshold", minExprThreshold, 
                                               ".numCellsExpressed", numCellsExpressed,
                                               ".downSampleProp", dsProp, 
                                               ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), ""), 
                                               ".pseudotime_per_branch.tsv", sep=""))
write.table(ptData, 
            tablefilename,sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)



##################################################
# Save object
# set_lower_bound=0; thisOutputDir = file.path(mainOutputDir, paste0( "trajectoryAnalysis.semisupervised", ifelse(set_lower_bound >0, paste0( ".set_lower_bound", set_lower_bound ), "") ) )
# minExprThreshold=0.1; numCellsExpressed=10; dsProp=1; DEqval=0.01
RdataFileName = file.path(thisOutputDir, paste0( ".RData.genome.filtered.semisupervised", 
                                                 paste0(".minExprThreshold", minExprThreshold), 
                                                 paste0(".numCellsExpressed", numCellsExpressed),
                                                 ifelse(!is.na(dsProp), paste0(".dsProp", dsProp), "" ),
                                                 ifelse(DEqval!=0.01, paste0(".DEqval", DEqval), "") 
) )

# save.image(file =  RdataFileName )
# load(file = RdataFileName )
