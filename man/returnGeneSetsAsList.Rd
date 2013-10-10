\name{returnGeneSetsAsList}
\alias{returnGeneSetsAsList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return pre-defined gene lists in several biomedical categories.
}
\description{
This function returns gene sets for use with other R functions.  These gene sets can include inputted lists of genes and files containing user-defined lists of genes, as well as a pre-made collection of brain, blood, and other biological lists.  The function returns gene lists associated with each category for use with other enrichment strategies (i.e., GSVA).
}
\usage{
returnGeneSetsAsList(
   fnIn = NULL, catNmIn = fnIn, 
   useBrainLists = FALSE, useBloodAtlases = FALSE, 
   useStemCellLists = FALSE, useBrainRegionMarkers = FALSE, 
   useImmunePathwayLists = FALSE, geneSubset=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnIn}{
A vector of file names containing user-defined lists.  These files must be in one of three specific formats (see details section).  The default (NULL) may only be used if one of the "use_____" parameters is TRUE.
}
  \item{catNmIn}{
A vector of category names corresponding to each fnIn.  This name will be appended to each overlap corresponding to that filename.  The default sets the category names as the corresponding file names.
}
  \item{useBrainLists}{
If TRUE, a pre-made set of brain-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
  \item{useBloodAtlases}{
If TRUE, a pre-made set of blood-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
  \item{useStemCellLists}{
If TRUE, a pre-made set of stem cell (SC)-derived enrichment lists will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  See references section for related references.
}
  \item{useBrainRegionMarkers}{
If TRUE, a pre-made set of enrichment lists for human brain regions will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  These lists are derived from data from the Allen Human Brain Atlas (http://human.brain-map.org/).  See references section for more details.
}
  \item{useImmunePathwayLists}{
If TRUE, a pre-made set of enrichment lists for immune system pathways will be added to any user-defined lists for enrichment comparison.  The default is FALSE.  These lists are derived from the lab of Daniel R Saloman.  See references section for more details.
}
  \item{geneSubset}{
A vector of gene (or other) identifiers.  If entered, only genes in this list will be returned in the output, otherwise all genes in each category will be returned (default, geneSubset=NULL).  
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
  \item{geneSets}{
A list of categories in alphabetical order, where each compnent of the list is a character vector of all genes corresponding to the named category.  For example: geneSets = list(category1=c("gene1","gene2"),category2=c("gene3","gene4","gene5"))
}
}
\references{
Please see the help file for userListEnrichment in the WGCNA library for references for the pre-defined lists.
}
\author{
Jeremy Miller
}
\examples{
# Example: Return a list of genes for various immune pathways
geneSets   = returnGeneSetsAsList(useImmunePathwayLists=TRUE)
geneSets[7:8]
}

\keyword{misc}
