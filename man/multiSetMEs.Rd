\name{multiSetMEs}
\alias{multiSetMEs}
\title{Calculate module eigengenes. }
\description{
  Calculates module eigengenes for several sets.
}
\usage{
multiSetMEs(exprData, 
            colors, 
            universalColors = NULL, 
            useSets = NULL, 
            useGenes = NULL,
            impute = TRUE, 
            nPC = 1, 
            align = "along average", 
            excludeGrey = FALSE,
            grey = if (is.null(universalColors)) {
                       if (is.numeric(colors)) 0 else "grey"
                   } else
                       if (is.numeric(universalColors)) 0 else "grey",
            subHubs = TRUE,
            trapErrors = FALSE, 
            returnValidOnly = trapErrors,
            softPower = 6,
            verbose = 1, indent = 0)
}
\arguments{
  \item{exprData}{Expression data in a multi-set format (see \code{\link{checkSets}}). A vector of lists,
with each list corresponding to one microarray dataset and expression data in the component \code{data},
that is \code{expr[[set]]$data[sample, probe]} is the expression of probe \code{probe} in sample
\code{sample} in dataset \code{set}. The number of samples can be different between the sets, but the
probes must be the same. }

  \item{colors}{A matrix of dimensions (number of probes, number of sets) giving the module assignment of
each gene in each set. The color "grey" is interpreted as unassigned.}

  \item{universalColors}{Alternative specification of module assignment. A single vector of length
(number of probes) giving the module assignment of each gene in all sets (that is the modules are common
to all sets). If given, takes precedence over \code{color}.}

  \item{useSets}{If calculations are requested in (a) selected set(s) only, the set(s) can be specified
here. Defaults to all sets.}

  \item{useGenes}{Can be used to restrict calculation to a subset of genes (the same subset in all
sets). If given, \code{validColors} in the returned list will only contain colors for the genes
specified in \code{useGenes}.}

  \item{impute}{Logical. If \code{TRUE}, expression data will be checked for the presence of \code{NA}
entries and if the latter are present, numerical data will be imputed, using function \code{impute.knn}
and probes from the same module as the missing datum. The function \code{impute.knn} uses a fixed random
seed giving repeatable results.} 

  \item{nPC}{Number of principal components to be calculated. If only eigengenes are needed, it is best
to set it to 1 (default). If variance explained is needed as well, use value \code{NULL}. This will cause
all principal components to be computed, which is slower.}

  \item{align}{Controls whether eigengenes, whose orientation is undetermined, should be aligned with
average expression (\code{align = "along average"}, the default) or left as they are (\code{align = ""}).
Any other value will trigger an error.}

  \item{excludeGrey}{Should the improper module consisting of 'grey' genes be excluded from the
eigengenes?}

  \item{grey}{Value of \code{colors} or \code{universalColors} (whichever applies)
designating the improper module. Note that if the appropriate colors argument is a
factor of numbers, the default value will be incorrect.}

  \item{subHubs}{Controls whether hub genes should be substituted for missing eigengenes. If
\code{TRUE}, each missing eigengene (i.e., eigengene whose calculation failed and the error
was trapped) will be replaced by a weighted average of the most connected hub genes in the
corresponding module. If this calculation fails, or if \code{subHubs==FALSE}, the value of
\code{trapErrors} will determine whether the offending module
will be removed or whether the function will issue an error and stop.}

  \item{trapErrors}{Controls handling of errors from that may arise when there are too many
\code{NA} entries in expression data. If \code{TRUE}, errors from calling these functions will be
trapped without abnormal exit.
If \code{FALSE}, errors will cause the function to stop. Note, however, that \code{subHubs} takes
precedence in the sense that if \code{subHubs==TRUE} and \code{trapErrors==FALSE}, an error will be
issued only if both the principal component and the hubgene calculations have failed. } 

  \item{returnValidOnly}{Boolean. Controls whether the returned data frames of module eigengenes 
contain columns
corresponding only to modules whose eigengenes or hub genes could be calculated correctly in every
set (\code{TRUE}), or whether
the data frame should have columns for each of the input color labels (\code{FALSE}).}

  \item{softPower}{The power used in soft-thresholding the adjacency matrix. Only used when the
hubgene approximation is necessary because the principal component calculation failed. It must be
non-negative. The default
value should only be changed if there is a clear indication that it leads to incorrect results.}

  \item{verbose}{Controls verbosity of printed progress messages. 0 means silent, up to (about) 5 the
verbosity gradually increases.}

  \item{indent}{A single non-negative integer controlling indentation of printed messages. 0 means no
indentation, each unit above that adds two spaces. }
}
\details{
This function calls \code{\link{moduleEigengenes}} for each set in \code{exprData}.

Module eigengene is defined as the first principal component of the expression matrix of the
corresponding module. The calculation may fail if the expression data has too many missing entries.
Handling of such errors is controlled by the arguments \code{subHubs} and
\code{trapErrors}.
If \code{subHubs==TRUE}, errors in principal component calculation will be trapped and a substitute
calculation of hubgenes will be attempted. If this fails as well, behaviour depends on
\code{trapErrors}: if \code{TRUE}, the offending
module will be ignored and the return value will allow the user to remove the module from further
analysis; if \code{FALSE}, the function will stop.
If \code{universalColors} is given, any offending
module will be removed from all sets (see \code{validMEs} in return value below). 

From the user's point of view, setting \code{trapErrors=FALSE} ensures that if the function returns
normally, there will be a valid eigengene (principal component or hubgene) for each of the input
colors. If the user sets \code{trapErrors=TRUE}, all calculational (but not input) errors will be
trapped, but the user should check the output (see below) to make sure all modules have a valid
returned eigengene.

While the principal component calculation can fail even on relatively sound data
(it does not take all that many "well-placed" \code{NA} to torpedo the
calculation),
it takes many more irregularities in the data for the hubgene calculation to
fail. In fact such a failure signals there likely is something seriously wrong with the data.

}
\value{
  A vector of lists similar in spirit to the input \code{exprData}. For each set there is a list with the
following components:
  \item{data}{Module eigengenes in a data frame, with each column corresponding to one eigengene.
The columns are named by the corresponding color with an \code{"ME"} prepended, e.g., \code{MEturquoise}
etc. Note that, when \code{trapErrors == TRUE} and \code{returnValidOnly==FALSE}, 
this data frame also contains entries corresponding to 
removed modules, if any. (\code{validMEs} below indicates which eigengenes are valid and \code{allOK}
whether all module eigengens were successfully calculated.) }
  \item{averageExpr}{If \code{align == "along average"}, a dataframe containing average normalized
expression in each module. The columns are named by the corresponding color with an \code{"AE"}
prepended, e.g., \code{AEturquoise} etc.}
  \item{varExplained}{A dataframe in which each column corresponds to a module, with the component
\code{varExplained[PC, module]} giving the variance of module \code{module} explained by the principal
component no. \code{PC}. This is only accurate if all principal components have been computed (input
\code{nPC = NULL}). At most 5 principal components are recorded in this dataframe.}
  \item{nPC}{A copy of the input \code{nPC}.}
  \item{validMEs}{A boolean vector. Each component (corresponding to the columns in \code{data})
is \code{TRUE} if the corresponding eigengene is valid, and \code{FALSE}
if it is invalid. Valid eigengenes include both principal components and their hubgene
approximations.
When \code{returnValidOnly==FALSE}, by definition all returned eigengenes are valid and the
entries of \code{validMEs} are all \code{TRUE}. }
  \item{validColors}{A copy of the input colors (\code{universalColors} if set, otherwise 
\code{colors[, set]}) with entries corresponding to invalid modules set to
\code{grey} if given, otherwise 0 if the appropriate input colors are numeric and "grey" otherwise.}
  \item{allOK}{Boolean flag signalling whether all eigengenes have been calculated correctly, either
as principal components or as the hubgene approximation. If \code{universalColors} is set, this flag
signals whether all eigengenes are valid in all sets.}
  \item{allPC}{Boolean flag signalling whether all returned eigengenes are principal components.
This flag (as well as the subsequent ones) is set independently for each set.}
  \item{isPC}{Boolean vector. Each component (corresponding to the columns in \code{eigengenes}) is
\code{TRUE} if the corresponding eigengene is the first principal component and \code{FALSE} if it
is the hubgene approximation or is invalid. }
  \item{isHub}{Boolean vector. Each component (corresponding to the columns in \code{eigengenes}) is
\code{TRUE} if the corresponding eigengene is the hubgene approximation and \code{FALSE} if it
is the first principal component or is invalid.}
  \item{validAEs}{Boolean vector. Each component (corresponding to the columns in \code{eigengenes}) is
\code{TRUE} if the corresponding module average expression is valid.}
  \item{allAEOK}{Boolean flag signalling whether all returned module average expressions contain
valid data. Note that \code{returnValidOnly==TRUE} does not imply \code{allAEOK==TRUE}: 
some invalid average expressions may be
returned if their corresponding eigengenes have been calculated correctly.}


}
\author{
  Peter Langfelder, \email{Peter.Langfelder@gmail.com}
}
\seealso{\code{\link{moduleEigengenes}}}
\keyword{misc}
