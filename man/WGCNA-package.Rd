\name{WGCNA-package}
\alias{WGCNA-package}
\alias{WGCNA}
\docType{package}
\title{
Weighted Gene Co-Expression Network Analysis
}
\description{
Functions necessary to perform Weighted Correlation Network Analysis. WGCNA is also known as weighted
gene co-expression network analysis when dealing with gene expression data. Many functions of WGCNA can
also be used for general association networks specified by a symmetric adjacency matrix.
}
\details{
\tabular{ll}{
Package: \tab WGCNA\cr
Version: \tab 1.00\cr
Date: \tab 2011-02-05\cr
Depends: \tab R (>= 2.3.0), stats, impute, grDevices, dynamicTreeCut (>= 1.20), utils, flashClust, qvalue,
Hmisc, splines \cr
Suggests: AnnotationDbi, GO.db, org.*.eg.db, infotheo, entropy, minet
ZipData: \tab no\cr
License: \tab GPL (>= 2)\cr
URL: \tab http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/\cr
}

Index:
\preformatted{
GTOMdist                Generalized Topological Overlap Measure
TOMdist                 Topological overlap matrix dissimilarity
TOMplot                 Graphical representation of the Topological Overlap Matrix
TOMsimilarity           Topological overlap matrix similarity
TOMsimilarityFromExpr   Topological overlap matrix similarity
WGCNA-package           Weighted Gene Co-Expression Network Analysis
addErrorBars            Add error bars to a barplot.
addGrid                 Add grid lines to an existing plot.
addGuideLines           Add vertical "guide lines" to a dendrogram plot
addTraitToMEs           Add trait information to multi-set module
                        eigengene structure
adjacency               Calculate network adjacency
alignExpr               Align expression data with given vector
automaticNetworkScreening
                        One-step automatic network gene screening
automaticNetworkScreeningGS
                        One-step automatic network gene screening with
                        external gene significance
AUV2predicted           Prediction of weighted mutual information adjacency matrix by correlation
bicor                   Biweight Midcorrelation
bicorAndPvalue          Biweight Midcorrelation and the associated p-value
blockwiseConsensusModules
                        Find consensus modules across several datasets.
blockwiseModules        Automatic network construction and module
                        detection
checkAdjMat             Check adjacency matrix
checkSets               Check structure and retrieve sizes of a group
                        of datasets
clusterCoef             Clustering coefficient calculation
collapseRows            Collapse Rows in Numeric Matrix
collectGarbage          Iterative garbage collection
colQuantileC            Fast colunm-wise quantile of a matrix
conformityBasedNetworkConcepts
                        Calculation of conformity-based network concepts
consensusMEDissimilarity
                        Consensus dissimilarity of module eigengenes.
consensusOrderMEs       Put close eigenvectors next to each other in
                        several sets.
consensusProjectiveKMeans
                        Consensus projective K-means (pre-)clustering
                        of expression data
cor                     Faster calculation of Pearson correlations
corAndPvalue            Correlation and the associated p-value
cor1                    Faster calculation of column correlations of a matrix
corFast                 Faster calculation of Pearson correlations
corPredictionSuccess    ~~function to do ... ~~
corPvalueFisher         Fisher's asymptotic p-value for correlation
corPvalueStudent        Student asymptotic p-value for correlation
correlationPreservation
                        Preservation of eigengene correlations
cutreeStatic            Constant height tree cut
cutreeStaticColor       Constant height tree cut using color labels
displayColors           Show colors used to label modules
dynamicMergeCut         Threshold for module merging
exportNetworkToVisANT   Export network data in format readable by VisANT
exportNetworkToCytoscape
                        Export network data in format readable by Cytoscape
fixDataStructure        Put single-set data into a form useful for multiset calculations
fundamentalNetworkConcepts
			Calculation of fundamental network concepts
GOenrichmentAnalysis    Calculate enrichment p-values of clusters in GO terms
goodGenes               Filter genes with too many missing entries
goodGenesMS             Filter genes with too many missing entries across multiple sets
goodSamples             Filter samples with too many missing entries
goodSamplesGenes        Iterative filtering of samples and genes with
                        too many missing entries
goodSamplesGenesMS      Iterative filtering of samples and genes with too many missing entries across
                        multiple data sets
goodSamplesMS           Filter samples with too many missing entries across multiple data sets
greenBlackRed           Green-black-red color sequence
greenWhiteRed           Green-white-red color sequence
hubGeneSignificance     Hubgene significance
initProgInd             Inline display of progress
intramodularConnectivity
                        Calculation of intramodular connectivity
keepCommonProbes        Keep probes that are shared among given data sets
labeledBarplot          Barplot with text or color labels
labeledHeatmap          Produce a labeled heatmap plot
labelPoints		Attempt to intelligently label points in a scatterplot
labels2colors           Convert numerical labels to colors
matchLabels             Relabel modules to best approximate a reference labeling
mergeCloseModules       Merge close modules in gene expression data
metaZfunction           Meta-analysis Z statistic
moduleColor.getMEprefix
                        Get the prefix used to label module eigengenes
moduleEigengenes        Calculate module eigengenes
moduleNumber            Fixed-height cut of a dendrogram
modulePreservation	Calculation of module preservation statistics
multiSetMEs             Calculate module eigengenes
mutualInfoAdjacency     Calculate weighted adjacency matrices based on mutual information
nPresent                Number of present data entries
nearestNeighborConnectivity
                        Connectivity to a constant number of nearest neighbors
nearestNeighborConnectivityMS
                        Connectivity to a constant number of nearest
                        neighbors across multiple data sets
networkConcepts         Calculations of network concepts
networkScreening        Network screening 
networkScreeningGS      Network screening with external gene significance
normalizeLabels         Transform numerical labels into normal order
numbers2colors          Color representation for a numeric variable
orderMEs                Put close eigenvectors next to each other
overlapTable            Overlap counts and Fisher exact tests for two sets of module labels
pickHardThreshold       Analysis of scale free topology for
                        hard-thresholding.
pickSoftThreshold       Analysis of scale free topology for
                        soft-thresholding
plotClusterTreeSamples
                        Annotated clustering dendrogram of microarray samples
plotColorUnderTree      Plot color rows under a dendrogram
plotDendroAndColors     Dendrogram plot with color annotation of objects
plotEigengeneNetworks   Eigengene network plot
plotMEpairs             Pairwise scatterplots of eigengenes
plotModuleSignificance
                        Barplot of module significance
plotNetworkHeatmap      Network heatmap plot
pmean                   Parallel mean
pmedian                 Parallel median
pquantile               Parallel quantile
preservationNetworkConnectivity
                        Network preservation calculations
projectiveKMeans        Projective K-means (pre-)clustering of
                        expression data
propVarExplained        Proportion of variance explained by eigengenes
randIndex               Calculation of (adjusted) Rand index
rankPvalue              Estimate the p-value for ranking consistently high (or low) on multiple lists
recutBlockwiseTrees     Repeat blockwise module detection from
                        pre-calculated data
recutConsensusTrees     Repeat blockwise consensus module detection
                        from pre-calculated data
redWhiteGreen           Red-white-green color sequence
relativeCorPredictionSuccess
                        Compare prediction success
removeGreyME            Removes the grey eigengene from a given
                        collection of eigengenes.
scaleFreeFitIndex	Calculation of fitting statistics for evaluating scale free topology fit.
scaleFreePlot           Visual check of scale-free topology
setCorrelationPreservation
                        Summary correlation preservation measure
sigmoidAdjacencyFunction
                        Sigmoid-type adacency function
signedKME               Signed eigengene-based connectivity
signumAdjacencyFunction
                        Hard-thresholding adjacency function
simulateDatExpr         Simulation of expression data
simulateDatExpr5Modules
                        
simulateEigengeneNetwork
                        Simulate eigengene network from a causal model
simulateModule          Simulate a gene co-expression module
simulateMultiExpr       Simulate multi-set expression data
simulateSmallLayer      Simulate small modules
sizeGrWindow            Open a graphics window of given width and height
softConnectivity        Calculation of soft (weighted) connectevity
spaste                  Space-less paste
standardColors          Colors this library uses for labeling modules
standardScreeningBinaryTrait
                        Standard screening for a binary trait
standardScreeningCensoredTime
                        Standard screening with regard to a Censored Time Variable
stdErr                  Standard error
unsignedAdjacency       Calculation of unsigned adjacency
vectorTOM               Topological overlap for a subset of the whole set of genes
vectorizeMatrix         Turn a matrix into a vector of non-redundant components
verboseBarplot          Barplot with error bars, annotated by Kruskal-Wallis p-value
verboseBoxplot          Boxplot annotated by a Kruskal-Wallis p-value
verboseScatterplot      Scatterplot annotated by regression line and p-value
}
}
\author{
Peter Langfelder <Peter.Langfelder@gmail.com> and Steve Horvath <SHorvath@mednet.ucla.edu>, with
contributions by Jun Dong, Jeremy Miller, Lin Song, Andy Yip, and Bin Zhang

Maintainer: Peter Langfelder <Peter.Langfelder@gmail.com>
}
\references{
Peter Langfelder and Steve Horvath (2008) WGCNA: an R package for weighted correlation network analysis. BMC
Bioinformatics 2008, 9:559

Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24

Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression Network Analysis. PLoS Comput Biol
4(8): e1000117

Yip A, Horvath S (2007) Gene network interconnectedness and the generalized topological overlap measure.
BMC Bioinformatics 2007, 8:22

Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a hierarchical cluster tree: the Dynamic
Tree Cut library for R. Bioinformatics. November/btm563

Langfelder P, Horvath S (2007) Eigengene networks for studying the relationships between co-expression
modules. BMC Systems Biology 2007, 1:54
}
\keyword{ package }
