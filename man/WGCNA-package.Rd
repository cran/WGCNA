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
Version: \tab 1.51\cr
Date: \tab 2016-03-08\cr
Depends: \tab R (>= 3.0), dynamicTreeCut (>= 1.62), fastcluster, Hmisc \cr
Imports: \tab stats, impute, grDevices, utils, splines, reshape, foreach, doParallel, matrixStats (>= 0.8.1),
GO.db, AnnotationDbi \cr 
Suggests: \tab org.Hs.eg.db, org.Mm.eg.db, infotheo, entropy, minet, survival \cr 
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
accuracyMeasures        Accuracy measures for a 2x2 confusion matrix
addErrorBars            Add error bars to a barplot.
addGrid                 Add grid lines to an existing plot.
addGuideLines           Add vertical "guide lines" to a dendrogram plot
addTraitToMEs           Add trait information to multi-set module
                        eigengene structure
adjacency               Calculate network adjacency
adjacency.fromSimilarity
                        Calculate network adjacency from a similarity matrix
adjacency.polyReg       Adjacency based on polynomial regression
adjacency.splineReg     Adjacency based on natural cubic spline regression
alignExpr               Align expression data with given vector
automaticNetworkScreening
                        One-step automatic network gene screening
automaticNetworkScreeningGS
                        One-step automatic network gene screening with
                        external gene significance
AFcorMI                 Prediction of weighted mutual information adjacency matrix 
                        by correlation
bicor                   Biweight Midcorrelation
bicorAndPvalue          Biweight Midcorrelation and the associated p-value
blockwiseConsensusModules
                        Find consensus modules across several datasets.
blockwiseIndividualTOMs
                        Calculate individual topological overlaps across 
                        multi-set data
blockwiseModules        Automatic network construction and module
                        detection
BloodLists              (data) Gene sets for user enrichment analysis 
blueWhiteRed            Blue-white-red color sequence
BrainLists              (data) Gene sets for user enrichment analysis 
BrainRegionMarkers      (data) Gene Markers for Regions of the Human Brain
checkAdjMat             Check adjacency matrix
checkSets               Check structure and retrieve sizes of a group
                        of datasets
checkSimilarity         Check a similarity matrix
chooseOneHubInEachModule
                        Choose a hub gene in each module
chooseTopHubInEachModule
                        Choose the top hub gene in each module
clusterCoef             Clustering coefficient calculation
coClustering            Cluster preservation based on co-clustering
coClustering.permutationTest
                        Permutation test for co-clustering 
collapseRows            Collapse Rows in Numeric Matrix
collapseRowsUsingKME    Selects one representative row per group based on kM
collectGarbage          Iterative garbage collection
colQuantileC            Fast colunm-wise quantile of a matrix
conformityBasedNetworkConcepts
                        Calculation of conformity-based network concepts
conformityDecomposition
                        Conformity vector(s) and factorizability measure(s) 
                        of a network
consensusDissTOMandTree
                        Consensus TOM-based dissimilarity and clustering tree
consensusKME            Consensus eigengene-based connectivity
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
coxRegressionResiduals
                        Deviance- and martingale residuals from a Cox regression model
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
ImmunePathwayLists      (data) Immune Pathways with Corresponding Gene Markers
initProgInd             Inline display of progress
intramodularConnectivity
intramodularConnectivity.fromExpr
                        Calculation of intramodular connectivity
keepCommonProbes        Keep probes that are shared among given data sets
kMEcomparisonScatterplot
                        Scatterplots of eigengene-based connectivity
labeledBarplot          Barplot with text or color labels
labeledHeatmap          Produce a labeled heatmap plot
labelPoints		Attempt to intelligently label points in a scatterplot
labels2colors           Convert numerical labels to colors
lowerTri2matrix         Reconstruct a symmetric matrix from a distance 
                        (lower-triangular) representation
matchLabels             Relabel modules to best approximate a reference labeling
mergeCloseModules       Merge close modules in gene expression data
metaAnalysis            Meta-analysis significance statistics
metaZfunction           Meta-analysis Z statistic
moduleColor.getMEprefix
                        Get the prefix used to label module eigengenes
moduleEigengenes        Calculate module eigengenes
moduleMergeUsingKME     Merge modules and reassign genes using kME
moduleNumber            Fixed-height cut of a dendrogram
modulePreservation	Calculation of module preservation statistics
multiSetMEs             Calculate module eigengenes
multiData.eigengeneSignificance
                        Calculate eigengene significance for multiple data sets
mutualInfoAdjacency     Calculate weighted adjacency matrices based on mutual information
nPresent                Number of present data entries
nearestNeighborConnectivity
                        Connectivity to a constant number of nearest neighbors
nearestNeighborConnectivityMS
                        Connectivity to a constant number of nearest
                        neighbors across multiple data sets
nearestCentroidPredictor
                        Nearest centroid predictor for two-class problems
networkConcepts         Calculations of network concepts
networkScreening        Network screening 
networkScreeningGS      Network screening with external gene significance
normalizeLabels         Transform numerical labels into normal order
numbers2colors          Color representation for a numeric variable
orderBranchesUsingHubGenes
                        Optimize dendrogram using branch swaps and reflections
orderMEs                Put close eigenvectors next to each other
overlapTable            Overlap counts and Fisher exact tests for two sets of module labels
overlapTableUsingKME    Determines significant overlap between modules in two networks based on kME tables
pickHardThreshold       Analysis of scale free topology for hard-thresholding
pickHardThreshold.fromSimilarity
                        Analysis of scale free topology for hard-thresholding
pickSoftThreshold       Analysis of scale free topology for soft-thresholding
pickSoftThreshold.fromSimilarity
                        Analysis of scale free topology for soft-thresholding
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
populationMeansInAdmixture
                        Estimation of population-specific mean values in an admixed population
pquantile               Parallel quantile
preservationNetworkConnectivity
                        Network preservation calculations
projectiveKMeans        Projective K-means (pre-)clustering of
                        expression data
propVarExplained        Proportion of variance explained by eigengenes
proportionsInAdmixture  
                        Estimation of proportion of pure populations in an admixture
qvalue                  q-value calculation from package qvalue 
randIndex               Calculation of (adjusted) Rand index
randomGLMpredictor	Ensemble predictor based on bagging of generalized linear models 
rankPvalue              Estimate the p-value for ranking consistently high (or low) on multiple lists
recutBlockwiseTrees     Repeat blockwise module detection from
                        pre-calculated data
recutConsensusTrees     Repeat blockwise consensus module detection
                        from pre-calculated data
redWhiteGreen           Red-white-green color sequence
reflectTwoBranches      Reflect branches in a dendrogram
relativeCorPredictionSuccess
                        Compare prediction success
removeGreyME            Removes the grey eigengene from a given
                        collection of eigengenes.
removePrincipalComponents
                        Remove leading principal components from data
returnGeneSetsAsLists   Return pre-defined gene lists in several biomedical categories.
scaleFreeFitIndex	Calculation of fitting statistics for evaluating scale free topology fit.
scaleFreePlot           Visual check of scale-free topology
SCsLists                (data) Stem Cell-Related Genes with Corresponding Gene Markers
selectBranch            Find a branch in a dendrogram
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
softConnectivity.fromSimilarity
                        Calculation of soft (weighted) connectevity
spaste                  Space-less paste
standardColors          Colors this library uses for labeling modules
standardScreeningBinaryTrait
                        Standard screening for a binary trait
standardScreeningCensoredTime
                        Standard screening with regard to a Censored Time Variable
stdErr                  Standard error
stratifiedBarplot	Bar plots of data across two splitting parameters
swapTwoBranches         Swap branches in a dendrogram
TrueTrait               Estimate the true trait underlying a list of surrogate markers
transposeBigData        Block-by-block transpose of large matrices
unsignedAdjacency       Calculation of unsigned adjacency
userListEnrichment      Measure enrichment between inputted and user-defined lists
vectorTOM               Topological overlap for a subset of the whole set of genes
vectorizeMatrix         Turn a matrix into a vector of non-redundant components
verboseBarplot          Barplot with error bars, annotated by Kruskal-Wallis p-value
verboseBoxplot          Boxplot annotated by a Kruskal-Wallis p-value
verboseScatterplot      Scatterplot annotated by regression line and p-value
verboseIplot            Scatterplot annotated by regression line, p-value, and color for density
votingLinearPredictor   Voting linear predictor
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

Peter Langfelder, Steve Horvath (2012)
Fast R Functions for Robust Correlations and Hierarchical Clustering.
Journal of Statistical Software, 46(11), 1-17.
\url{http://www.jstatsoft.org/v46/i11/}


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
