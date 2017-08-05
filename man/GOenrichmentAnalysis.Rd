\name{GOenrichmentAnalysis}
\alias{GOenrichmentAnalysis}
\title{ Calculation of GO enrichment (experimental)}
\description{
NOTE: GOenrichmentAnalysis is deprecated. Please use function
          enrichmentAnalysis from R package anRichment, available from
          https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/

WARNING: This function should be considered experimental. The arguments and resulting values (in
particular, the enrichment p-values) are not
yet finalized and may change in the future. The function should only be used to get a quick and rough
overview of GO enrichment in the modules in a data set; for a publication-quality analysis, please use an
established tool.

  Using Bioconductor's annotation packages, this function calculates enrichments and returns terms with
best enrichment values.
}
\usage{
GOenrichmentAnalysis(labels, 
                     entrezCodes, 
                     yeastORFs = NULL,
                     organism = "human", 
                     ontologies = c("BP", "CC", "MF"), 
                     evidence = "all",
                     includeOffspring = TRUE, 
                     backgroundType = "givenInGO",
                     removeDuplicates = TRUE,
                     leaveOutLabel = NULL,
                     nBestP = 10, pCut = NULL, 
                     nBiggest = 0, 
                     getTermDetails = TRUE,
                     verbose = 2, indent = 0)
}
\arguments{
  \item{labels}{ cluster (module, group) labels of genes to be analyzed. Either a single vector, or a
matrix. In the matrix case, each column will be analyzed separately; analyzing a collection of module
assignments in one function call will be faster than calling the function several tinmes. For each row,
the labels in all columns must correspond to the same gene specified in \code{entrezCodes}.  }

  \item{entrezCodes}{ Entrez (a.k.a. LocusLink) codes of the genes whose labels are given in
                      \code{labels}. A single vector; the i-th entry corresponds to row i of the matrix
\code{labels} (or to the i-the entry if \code{labels} is a vector). }

  \item{yeastORFs}{ if \code{organism=="yeast"} (below), this argument can be used to input yeast open
reading frame (ORF) identifiers instead of Entrez codes. Since the GO mappings for yeast are provided in
terms of ORF identifiers, this may lead to a more accurate GO enrichment analysis. If given, the argument
\code{entrezCodes} is ignored. }

  \item{organism}{ character string specifying the organism for which to perform the analysis. Recognized
values are (unique abbreviations of) \code{"human", "mouse", "rat", "malaria", "yeast", "fly", "bovine",
"worm", "canine", "zebrafish", "chicken"}. }

  \item{ontologies}{ vector of character strings specifying GO ontologies to be included in the analysis. 
Can be any subset of \code{"BP", "CC", "MF"}. The result will contain the terms with highest enrichment
in each specified category, plus a separate list of terms with best enrichment in all ontologies
combined. }

  \item{evidence}{ vector of character strings specifying admissible evidence for each gene in its
specific term, or "all" for all evidence codes. See Details or http://www.geneontology.org/GO.evidence.shtml for
available evidence codes and their meaning.}

  \item{includeOffspring}{ logical: should genes belonging to the offspring of each term be included in
the term? As a default, only genes belonging directly to each term are associated with the term. Note
that the calculation of enrichments with offspring included can be quite slow for large data sets.}

  \item{backgroundType}{specification of the background to use. Recognized values are (unique
abbreviations of) \code{"allGiven", "allInGO", "givenInGO"}, meaning that the functions will take all
genes given in \code{labels} as backround (\code{"allGiven"}), all genes present in any of the GO
categories (\code{"allInGO"}), or the intersection of given genes and genes present in GO
(\code{"givenInGO"}). The default is recommended for genome-wide enrichment studies. }

  \item{removeDuplicates}{logical: should duplicate entries in \code{entrezCodes} be removed? If
\code{TRUE}, only the first occurence of each unique Entrez code will be kept. The cluster labels
\code{labels} will be adjusted accordingly.}

  \item{leaveOutLabel}{optional specifications of module labels for which enrichment calculation is not
desired. Can be a single label or a vector of labels to be ignored. However, if in any of the sets no
labels are left to calculate enrichment of, the function will stop with an error.}

  \item{nBestP}{ specifies the number of terms with highest enrichment whose detailed information will be
returned. }

  \item{pCut}{ alternative specification of terms to be returned: all terms whose enrichment p-value is
more significant than \code{pCut} will be returned. If \code{pCut} is given, \code{nBestP} is ignored. }

  \item{nBiggest}{ in addition to returning terms with highest enrichment, terms that contain most of the
genes in each cluster can be returned by specifying the number of biggest terms per cluster to be
returned. This may be useful for development and testing purposes. }

  \item{getTermDetails}{ logical indicating whether detailed information on the most enriched terms should
be returned. }

  \item{verbose}{ integer specifying the verbosity of the function. Zero means silent, positive values
will cause the function to print progress reports.}

  \item{indent}{ integer specifying indentation of the diagnostic messages. Zero means no indentation,
each unit adds two spaces.}

}
\details{
  This function is basically a wrapper for the annotation packages available from Bioconductor. It
requires the packages GO.db, AnnotationDbi, and org.xx.eg.db, where xx is the code corresponding to the
organism that the user wishes to analyze (e.g., Hs for human Homo Sapiens, Mm for mouse Mus Musculus
etc). For each cluster specified in the input, the function calculates all enrichments in the specified
ontologies, and collects information about the terms with highest enrichment. The enrichment p-value is
calculated using Fisher exact test. As background we use all of the supplied genes that are present in at
least one term in GO (in any of the ontologies). 

For best results, the newest annotation libraries should be used. Because of the way
Bioconductor is set up, to get the newest annotation libraries you may have to use the current version of
R. 

According to http://www.geneontology.org/GO.evidence.shtml, the following codes are used by GO:
\preformatted{
  Experimental Evidence Codes
      EXP: Inferred from Experiment
      IDA: Inferred from Direct Assay
      IPI: Inferred from Physical Interaction
      IMP: Inferred from Mutant Phenotype
      IGI: Inferred from Genetic Interaction
      IEP: Inferred from Expression Pattern

  Computational Analysis Evidence Codes
      ISS: Inferred from Sequence or Structural Similarity
      ISO: Inferred from Sequence Orthology
      ISA: Inferred from Sequence Alignment
      ISM: Inferred from Sequence Model
      IGC: Inferred from Genomic Context
      IBA: Inferred from Biological aspect of Ancestor
      IBD: Inferred from Biological aspect of Descendant
      IKR: Inferred from Key Residues
      IRD: Inferred from Rapid Divergence
      RCA: inferred from Reviewed Computational Analysis

  Author Statement Evidence Codes
      TAS: Traceable Author Statement
      NAS: Non-traceable Author Statement

  Curator Statement Evidence Codes
      IC: Inferred by Curator
      ND: No biological Data available

  Automatically-assigned Evidence Codes
      IEA: Inferred from Electronic Annotation

  Obsolete Evidence Codes
      NR: Not Recorded 
}
}
\value{
  A list with the following components:
  \item{keptForAnalysis }{ logical vector with one entry per given gene. \code{TRUE} if the entry was
used for enrichment analysis. Depending on the setting of \code{removeDuplicates} above, only a single
entry per gene may be used. }

  \item{inGO }{ logical vector with one entry per given gene. \code{TRUE} if the gene belongs to any GO
term, \code{FALSE} otherwise. Also \code{FALSE} for genes not used for the analysis because of
duplication. }

If input \code{labels} contained only one vector of labels, the following components:

  \item{countsInTerms }{ a matrix whose rows correspond to given cluster, and whose columns correspond to
GO terms, contaning number of genes in the intersection of the corresponding module and GO term. Row and
column names are set appropriately.}

  \item{enrichmentP}{a matrix whose rows correspond to given cluster, and whose columns correspond to 
GO terms, contaning enrichment p-values of each term in each cluster. Row and
column names are set appropriately.}

  \item{bestPTerms}{a list of lists with each inner list corresponding to an ontology given in
\code{ontologies} in input, plus one component corresponding to all given ontologies combined. 
The name of each component is set appropriately. Each inner list contains two components: 
\code{enrichment} is
a data frame containing the highest enriched terms for each module; and \code{forModule} is a list of
lists with one inner list per module, appropriately named. Each inner list contains one component per term.
If input \code{getTermDeyails} is \code{TRUE}, 
this component is yet another list and contains components \code{termName} (term name), 
\code{enrichmentP} (enrichment P value), \code{termDefinition} (GO term definition), 
\code{termOntology} (GO term ontology), \code{geneCodes} (Entrez codes of module genes in this term),
\code{genePositions} (indices of the genes listed in \code{geneCodes} within the given \code{labels}). 
Thus, to obtain information on say the second term of the 5th module in ontology BP, 
one can look at the appropriate row of \code{bestPTerms$BP$enrichment}, or one can reference
\code{bestPTerms$BP$forModule[[5]][[2]]}. The author of the function apologizes for any confusion this
structure of the output may cause. }

  \item{biggestTerms}{a list of the same format as \code{bestPTerms}, containing information about the
terms with most genes in the module for each supplied ontology. }

If input \code{labels} contained more than one vector, instead of the above components the return value
contains a list named \code{setResults} that has one component per given set; each component is a list
containing the above components for the corresponding set.

}
\author{ Peter Langfelder }
\seealso{ Bioconductor's annotation packages such as GO.db and organism-specific annotation
packages such as org.Hs.eg.db. }

\keyword{ misc }% __ONLY ONE__ keyword per line
