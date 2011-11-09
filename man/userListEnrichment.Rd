\name{userListEnrichment}
\alias{userListEnrichment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Measure enrichment between inputted and user-defined lists
}
\description{
This function measures list enrichment between inputted lists of genes and files containing user-defined lists of genes.  Significant enrichment is measured using a hypergeometric test.  A pre-made collection of brain-related lists can also be loaded.  The function writes the significant enrichments to a file, but also returns all overlapping genes across all comparisons.
}
\usage{
userListEnrichment(geneR, labelR, fnIn=NULL, catNmIn = fnIn, nameOut = "enrichment.csv", 
                   useBrainLists = FALSE, useBloodAtlases=FALSE, omitCategories = "grey", 
                   outputCorrectedPvalues = TRUE, useStemCellLists = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneR}{
A vector of gene (or other) identifiers.  This vector should include ALL genes in your analysis (i.e., the genes correspoding to your labeled lists AND the remaining background reference genes).  
}
  \item{labelR}{
A vector of labels (for example, module assignments) corresponding to the geneR list.  NOTE: For all background reference genes that have no corresponding label, use the label "background" (or any label included in the omitCategories parameter).
}
  \item{fnIn}{
A vector of file names containing user-defined lists.  These files must be in one of three specific formats (see details section).  The default (NULL) may only be used if one of the "use_____" parameters is TRUE.
}
  \item{catNmIn}{
A vector of category names corresponding to each fnIn.  This name will be appended to each overlap corresponding to that filename.  The default sets the category names as the corresponding file names.
}
  \item{nameOut}{
Name of the file where the output enrichment information will be written.  (Note that this file includes only a subset of what is returned by the function.)
}
  \item{useBrainLists}{
If TRUE, a pre-made set of brain-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
  \item{useBloodAtlases}{
If TRUE, a pre-made set of blood-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
  \item{omitCategories}{
Any labelR entries corresponding to these categories will be ignored.  The default ("grey") will ignore unassigned genes in a standard WGCNA network.
}
  \item{outputCorrectedPvalues}{
If TRUE (default) only pvalues that are significant after correcting for multiple comparisons (using Bonferroni method) will be outputted to nameOut.  Otherwise the uncorrected p-values will be outputted to the file.  Note that both sets of p-values for all comparisons are reported in the returned "pValues" parameter.
}
  \item{useStemCellLists}{
If TRUE, a pre-made set of stem cell (SC)-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
}
\details{
User-inputted files for fnIn can be in one of three formats:

1) Text files (must end in ".txt") with one list per file, where the first line is the list descriptor and the remaining lines are gene names corresponding to that list, with one gene per line.  For example
Ribosome
RPS4
RPS8
...

2) Gene / category files (must be csv files), where the first line is the column headers corresponding to Genes and Lists, and the remaining lines correspond to the genes in each list, for any number of genes and lists.  For example:
Gene, Category
RPS4, Ribosome
RPS8, Ribosome
...
NDUF1, Mitohcondria
NDUF3, Mitochondria
...
MAPT, AlzheimersDisease
PSEN1, AlzheimersDisease
PSEN2, AlzheimersDisease
...

3) Module membership (kME) table in csv format.  Currently, the module assignment is the only thing that is used, so as long as the Gene column is 2nd and the Module column is 3rd, it doesn't matter what is in the other columns.  For example,
PSID, Gene, Module, <other columns>
<psid>, RPS4, blue, <other columns>
<psid>, NDUF1, red, <other columns>
<psid>, RPS8, blue, <other columns>
<psid>, NDUF3, red, <other columns>
<psid>, MAPT, green, <other columns>
...
}
\value{
  \item{pValues}{
A matrix showing the number of overlapping genes and both the uncorrected and Bonferroni corrected p-values for every pair of list overlaps tested.
}
  \item{ovGenes}{
A list of character vectors corresponding to the overlapping genes for every pair of list overlaps tested.  Specific overlaps can be found by typing <variableName>$ovGenes$'<labelR> -- <comparisonCategory>'.  See example below.
}
  \item{sigOverlaps}{
Identical information that is written to nameOut.  P-values (corrected or uncorrected, depending on outputCorrectedPvalues) corresponding to all significant enrichments.
}
}
\references{
The primary reference for this function is: Miller JA, Cai C, Langfelder P, Geschwind DH, Kurian SM, Salomon DR, Horvath S. (2011) Strategies for aggregating gene expression data: the collapseRows R function. BMC Bioinformatics 12:322.

If you have any suggestions for lists to add to this function, please e-mail Jeremy Miller at jeremyinla@gmail.com

-------------------------------------
References for the pre-defined brain lists (useBrainLists=TRUE, in alphabetical order by category descriptor) are as follows:

ABA ==> Cell type markers from: Lein ES, et al. (2007) Genome-wide atlas of gene expression in the adult mouse brain. Nature 445:168-176.

ADvsCT_inCA1  ==> Lists of genes found to be increasing or decreasing with Alzheimer's disease in 3 studies:
1. Blalock => Blalock E, Geddes J, Chen K, Porter N, Markesbery W, Landfield P (2004) Incipient Alzheimer's disease: microarray correlation analyses reveal major transcriptional and tumor suppressor responses. PNAS 101:2173-2178.
2. Colangelo => Colangelo V, Schurr J, Ball M, Pelaez R, Bazan N, Lukiw W (2002) Gene expression profiling of 12633 genes in Alzheimer hippocampal CA1: transcription and neurotrophic factor down-regulation and up-regulation of apoptotic and pro-inflammatory signaling. J Neurosci Res 70:462-473. 
3. Liang => Liang WS, et al (2008) Altered neuronal gene expression in brain regions differentially affected by Alzheimer's disease: a reference data set. Physiological genomics 33:240-56.

Blalock_AD ==> Modules from a network using the data from: Blalock E, Geddes J, Chen K, Porter N, Markesbery W, Landfield P (2004) Incipient Alzheimer's disease: microarray correlation analyses reveal major transcriptional and tumor suppressor responses. PNAS 101:2173-2178. 

CA1vsCA3 ==> Lists of genes enriched in CA1 and CA3 relative to other each and to other areas of the brain, from several studies:
1. Ginsberg => Ginsberg SD, Che S (2005) Expression profile analysis within the human hippocampus: comparison of CA1 and CA3 pyramidal neurons. J Comp Neurol 487:107-118.
2. Lein => Lein E, Zhao X, Gage F (2004) Defining a molecular atlas of the hippocampus using DNA microarrays and high-throughput in situ hybridization. J Neurosci 24:3879-3889.
3. Newrzella => Newrzella D, et al (2007) The functional genome of CA1 and CA3 neurons under native conditions and in response to ischemia. BMC Genomics 8:370.
4. Torres => Torres-Munoz JE, Van Waveren C, Keegan MG, Bookman RJ, Petito CK (2004) Gene expression profiles in microdissected neurons from human hippocampal subregions. Brain Res Mol Brain Res 127:105-114.
5. GorLorT => In either Ginsberg or Lein or Torres list.

Cahoy ==> Definite (10+ fold) and probable (1.5+ fold) enrichment from: Cahoy JD, et al. (2008) A transcriptome database for astrocytes, neurons, and oligodendrocytes: A new resource for understanding brain development and function. J Neurosci 28:264-278.

CTX ==> Modules from the CTX (cortex) network from: Oldham MC, et al. (2008) Functional organization of the transcriptome in human brain. Nat Neurosci 11:1271-1282.

DiseaseGenes ==> Probable (C or better rating as of 16 Mar 2011) and possible (all genes in database as of ~2008) genetics-based disease genes from: http://www.alzforum.org/

EarlyAD ==> Genes whose expression is related to cognitive markers of early Alzheimer's disease vs. non-demented controls with AD pathology, from: Parachikova, A., et al (2007) Inflammatory changes parallel the early stages of Alzheimer disease. Neurobiology of Aging 28:1821-1833.

HumanChimp ==> Modules showing region-specificity in both human and chimp from: Oldham MC, Horvath S, Geschwind DH (2006) Conservation and evolution of gene coexpression networks in human and chimpanzee brains. Proc Natl Acad Sci USA 103: 17973-17978.

HumanMeta ==> Modules from the human network from: Miller J, Horvath S, Geschwind D (2010) Divergence of human and mouse brain transcriptome highlights Alzheimer disease pathways. Proc Natl Acad Sci 107:12698-12703.

JAXdiseaseGene ==> Genes where mutations in mouse and/or human are known to cause any disease.  WARNING: this list represents an oversimplification of data!  This list was created from the Jackson Laboratory: Bult CJ, Eppig JT, Kadin JA, Richardson JE, Blake JA; Mouse Genome Database Group (2008) The Mouse Genome Database (MGD): Mouse biology and model systems. Nucleic Acids Res 36 (database issue):D724-D728.

Lu_Aging ==> Modules from a network using the data from: Lu T, Pan Y, Kao S-Y, Li C, Kohane I, Chan J, Yankner B (2004) Gene regulation and DNA damage in the ageing human brain. Nature 429:883-891.

MicroglialMarkers ==> Markers for microglia and macrophages from several studies: 
1. GSE772 => Gan L, et al. (2004) Identification of cathepsin B as a mediator of neuronal death induced by Abeta-activated microglial cells using a functional genomics approach. J Biol Chem 279:5565-5572.
2. GSE1910 => Albright AV, Gonzalez-Scarano F (2004) Microarray analysis of activated mixed glial (microglia) and monocyte-derived macrophage gene expression. J Neuroimmunol 157:27-38.
3. AitGhezala => Ait-Ghezala G, Mathura VS, Laporte V, Quadros A, Paris D, Patel N, et al. Genomic regulation after CD40 stimulation in microglia: relevance to Alzheimer's disease. Brain Res Mol Brain Res 2005;140(1-2):73-85.
4. 3treatments_Thomas => Thomas, DM, Francescutti-Verbeem, DM, Kuhn, DM (2006) Gene expression profile of activated microglia under conditions associated with dopamine neuronal damage. The FASEB Journal 20:515-517.

MitochondrialType ==> Mitochondrial genes from the somatic vs. synaptic fraction of mouse cells from: Winden KD, et al. (2009) The organization of the transcriptional network in specific neuronal classes. Mol Syst Biol 5:291.

MO ==> Markers for many different things provided to my by Mike Oldham.  These were originally from several sources:
1. 2+_26Mar08 => Genetics-based disease genes in two or more studies from http://www.alzforum.org/ (compiled by Mike Oldham).
2. Bachoo => Bachoo, R.M. et al. (2004) Molecular diversity of astrocytes with implications for neurological disorders. PNAS 101, 8384-8389.
3. Foster => Foster, LJ, de Hoog, CL, Zhang, Y, Zhang, Y, Xie, X, Mootha, VK, Mann, M. (2006) A Mammalian Organelle Map by Protein Correlation Profiling. Cell 125(1): 187-199.
4. Morciano => Morciano, M. et al. Immunoisolation of two synaptic vesicle pools from synaptosomes: a proteomics analysis. J. Neurochem. 95, 1732-1745 (2005).
5. Sugino => Sugino, K. et al. Molecular taxonomy of major neuronal classes in the adult mouse forebrain. Nat. Neurosci. 9, 99-107 (2006).

MouseMeta ==> Modules from the mouse network from: Miller J, Horvath S, Geschwind D (2010) Divergence of human and mouse brain transcriptome highlights Alzheimer disease pathways. Proc Natl Acad Sci 107:12698-12703.

Sugino/Winden ==> Conservative list of genes in modules from the network from: Winden K, Oldham M, Mirnics K, Ebert P, Swan C, Levitt P, Rubenstein J, Horvath S, Geschwind D (2009). The organization of the transcriptional network in specific neuronal classes. Molecular systems biology 5.
NOTE: Original data came from this neuronal-cell-type-selection experiment in mouse: Sugino K, Hempel C, Miller M, Hattox A, Shapiro P, Wu C, Huang J, Nelson S (2006). Molecular taxonomy of major neuronal classes in the adult mouse forebrain. Nat Neurosci 9:99-107

Voineagu ==> Several Autism-related gene categories from: Voineagu I, Wang X, Johnston P, Lowe JK, Tian Y, Horvath S, Mill J, Cantor RM, Blencowe BJ, Geschwind DH. (2011). Transcriptomic analysis of autistic brain reveals convergent molecular pathology. Nature 474(7351):380-4

-------------------------------------
References for the pre-defined blood atlases (useBloodAtlases=TRUE, in alphabetical order by category descriptor) are as follows:

Grigoryev ==> Grigoryev YA, Kurian SM, Avnur Z, Borie D, Deng J, et al. (2010).  Deconvoluting post-transplant immunity: cell subset-specific mapping reveals pathways for activation and expansion of memory T, monocytes and B cells.  PLoS One. 5(10):e13358.

HAEMatlas ==> Watkins NA, Gusnanto A, de Bono B, De S, Miranda-Saavedra D, et al. (2009). A HaemAtlas: characterizing gene expression in differentiated human blood cells. Blood. 113(19):e1-9.

IRIS ==> Abbas AB, Baldwin D, Ma Y, Ouyang W, Gurney A, et al. (2005). Immune response in silico (IRIS): immune-specific genes identified from a compendium of microarray expression data.  Genes Immun. 6(4):319-31.

Willinger ==> Willinger T, Freeman T, Hasegawa H, McMichael AJ, Callan MF (2005). Molecular signatures distinguish human central memory from effector memory CD8 T cell subsets. J Immunol. 175(9):5895-903. 

-------------------------------------
References for the pre-defined stem cell (SC) lists (useStemCellLists=TRUE, in alphabetical order by category descriptor) are as follows:

Cui ==> genes differentiating erythrocyte precursors (CD36+ cells) from multipotent human primary
hematopoietic stem cells/progenitor cells (CD133+ cells), from: Cui K, Zang C, Roh TY, Schones DE, Childs
RW, Peng W, Zhao K. (2009). Chromatin signatures in multipotent human hematopoietic stem cells indicate the
fate of bivalent genes during differentiation. Cell Stem Cell 4:80-93

Lee ==> gene lists related to Polycomb proteins in human embryonic SCs, from (a highly-cited paper!): Lee
TI, Jenner RG, Boyer LA, Guenther MG, Levine SS, Kumar RM, Chevalier B, Johnstone SE, Cole MF, Isono K, et
al. (2006) Control of developmental regulators by polycomb in human embryonic stem cells. Cell 125:301-313
}
\author{
Jeremy Miller
}
\examples{
# Example: first, read in some gene names and split them into categories
data(BrainLists);
listGenes = unique(as.character(BrainLists[,1]))
set.seed(100)
geneR = sort(sample(listGenes,10000))
categories = sort(rep(standardColors(10),1000))
categories[sample(1:10000,1000)] = "grey"
write(c("TESTLIST1",geneR[700:1500], sep="\n"),"TESTLIST1.txt")
write(c("TESTLIST2",geneR[1600:2400],sep="\n"),"TESTLIST2.txt")

# Now run the function!
testResults = userListEnrichment(geneR, labelR=categories, fnIn=c("TESTLIST1.txt","TESTLIST2.txt"), catNmIn=c("TEST1","TEST2"), 
nameOut = "testEnrichment.csv",useBrainLists=TRUE, omitCategories ="grey")

# To see a list of all significant enrichments, either open the file "testEnrichments.csv" in the current directory, or type:
testResults$sigOverlaps

# To see all of the overlapping genes between two categories (whether or not the p-value is significant), type restResults$ovGenes$'<labelR> -- <comparisonCategory>'.  For example:
testResults$ovGenes$"black -- TESTLIST1__TEST1"
testResults$ovGenes$"red -- salmon_M12_Ribosome__HumanMeta"

# More detailed overlap information is in the pValue output.  For example:
head(testResults$pValue)
}

\keyword{misc}
