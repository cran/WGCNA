\name{moduleEigengenes}
\alias{moduleEigengenes}
\title{Calculate module eigengenes.}
\description{
  Calculates module eigengenes (1st principal component) of modules in a given single dataset.
}
\usage{
moduleEigengenes(expr, 
                 colors, 
                 impute = TRUE, 
                 nPC = 1, 
                 align = "along average", 
                 excludeGrey = FALSE, 
                 grey = if (is.numeric(colors)) 0 else "grey",
                 subHubs = TRUE,
                 trapErrors = FALSE, 
                 returnValidOnly = trapErrors, 
                 softPower = 6,
                 scale = TRUE,
                 verbose = 0, indent = 0)
}
\arguments{
  \item{expr}{Expression data for a single set in the form of a data frame where rows are samples and
columns are genes (probes).}

  \item{colors}{A vector of the same length as the number of probes in \code{expr}, giving module color
for all probes (genes). Color \code{"grey"} is reserved for unassigned genes. } 

  \item{impute}{If \code{TRUE}, expression data will be checked for the presence of \code{NA} entries and
if the latter are present, numerical data will be imputed, using function \code{impute.knn} and probes from
the same module as the missing datum. The function \code{impute.knn} uses a fixed random seed giving
repeatable results.}

  \item{nPC}{Number of principal components and variance explained entries to be calculated. Note
that only the first principal component is returned; the rest are used only for the calculation of
proportion of variance explained. The number of returned variance explained entries is
currently \code{min(nPC, 10)}. If given \code{nPC} is greater than 10, a warning is issued.}

  \item{align}{Controls whether eigengenes, whose orientation is undetermined, should be aligned with
average expression (\code{align = "along average"}, the default) or left as they are (\code{align = ""}).
Any other value will trigger an error.}

  \item{excludeGrey}{Should the improper module consisting of 'grey' genes be excluded from the
eigengenes?}

  \item{grey}{Value of \code{colors} designating the improper module. Note that if \code{colors} is a
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

  \item{returnValidOnly}{logical; controls whether the returned data frame of module eigengenes 
contains columns
corresponding only to modules whose eigengenes or hub genes could be calculated correctly 
(\code{TRUE}), or whether
the data frame should have columns for each of the input color labels (\code{FALSE}).}

  \item{softPower}{The power used in soft-thresholding the adjacency matrix. Only used when the
hubgene approximation is necessary because the principal component calculation failed. It must be
non-negative. The default
value should only be changed if there is a clear indication that it leads to incorrect results.}

  \item{scale}{logical; can be used to turn off scaling of the expression data before calculating the
singular value decomposition. The scaling should only be turned off if the data has been scaled
previously, in which case the function can run a bit faster. Note however that the function first
imputes, then scales the expression data in each module. If the expression contain missing data, scaling
outside of the function and letting the function impute missing data may lead to slightly different
results than if the data is scaled within the function.}

  \item{verbose}{Controls verbosity of printed progress messages. 0 means silent, up to (about) 5 the
verbosity gradually increases.}

  \item{indent}{A single non-negative integer controlling indentation of printed messages. 0 means no
indentation, each unit above that adds two spaces. }

}
\details{
  Module eigengene is defined as the first principal component of the expression matrix of the
corresponding module. The calculation may fail if the expression data has too many missing entries.
Handling of such errors is controlled by the arguments \code{subHubs} and 
\code{trapErrors}. 
If \code{subHubs==TRUE}, errors in principal component calculation will be trapped and a substitute
calculation of hubgenes will be attempted. If this fails as well, behaviour depends on
\code{trapErrors}: if \code{TRUE}, the offending 
module will be ignored and the return value will allow the user to remove the module from further
analysis; if \code{FALSE}, the function will stop. 

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
  A list with the following components:
  \item{eigengenes}{Module eigengenes in a dataframe, with each column corresponding to one eigengene.
The columns are named by the corresponding color with an \code{"ME"} prepended, e.g., \code{MEturquoise}
etc. If \code{returnValidOnly==FALSE}, module eigengenes whose calculation failed have 
all components set to \code{NA}.}
  \item{averageExpr}{If \code{align == "along average"}, a dataframe containing average normalized
expression in each module. The columns are named by the corresponding color with an \code{"AE"}
prepended, e.g., \code{AEturquoise} etc.}
  \item{varExplained}{A dataframe in which each column corresponds to a module, with the component
\code{varExplained[PC, module]} giving the variance of module \code{module} explained by the principal
component no. \code{PC}. The calculation is exact irrespective of the number of computed principal
components. At most 10 variance explained values are recorded in this dataframe.}
  \item{nPC}{A copy of the input \code{nPC}.} 
  \item{validMEs}{A boolean vector. Each component (corresponding to the columns in \code{data})
is \code{TRUE} if the corresponding eigengene is valid, and \code{FALSE}
if it is invalid. Valid eigengenes include both principal components and their hubgene
approximations.
When \code{returnValidOnly==FALSE}, by definition all returned eigengenes are valid and the
entries of \code{validMEs} are all \code{TRUE}. }
  \item{validColors}{A copy of the input colors with entries corresponding to invalid modules set to
\code{grey} if given, otherwise 0 if \code{colors} is numeric and "grey" otherwise.}
  \item{allOK}{Boolean flag signalling whether all eigengenes have been calculated correctly, either
as principal components or as the hubgene average approximation.}
  \item{allPC}{Boolean flag signalling whether all returned eigengenes are principal components.}
  \item{isPC}{Boolean vector. Each component (corresponding to the columns in \code{eigengenes}) is
\code{TRUE} if the corresponding eigengene is the first principal component and \code{FALSE} if it
is the hubgene approximation or is invalid.}
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
\references{ 
  Zhang, B. and Horvath, S. (2005), "A General Framework for Weighted Gene Co-Expression Network
Analysis", Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17}
\author{ 
Steve Horvath \email{SHorvath@mednet.ucla.edu}, Peter Langfelder
\email{Peter.Langfelder@gmail.com} }
\seealso{\code{\link{svd}}, \code{impute.knn}} 
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
