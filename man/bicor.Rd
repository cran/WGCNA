\name{bicor}
\alias{bicor}
\title{ Biweight Midcorrelation }
\description{
  Calculate biweight midcorrelation efficiently for matrices.
}
\usage{
bicor(x, y = NULL, 
      robustX = TRUE, robustY = TRUE, 
      use = "all.obs", 
      maxPOutliers = 1,
      quick = 0,
      pearsonFallback = "individual",
      cosine = FALSE, 
      cosineX = cosine,
      cosineY = cosine,
      nThreads = 0, 
      verbose = 0, indent = 0)
}
\arguments{

  \item{x}{ a vector or matrix-like numeric object }

  \item{y}{ a vector or matrix-like numeric object }

  \item{robustX}{ use robust calculation for \code{x}?}

  \item{robustY}{ use robust calculation for \code{y}?}

  \item{use}{ specifies handling of \code{NA}s. One of (unique abbreviations of) "all.obs",
"pairwise.complete.obs". }

  \item{maxPOutliers}{ specifies the maximum percentile of data that can be considered outliers on either
side of the median separately. For each side of the median, if
higher percentile than \code{maxPOutliers} is considered an outlier by the weight function based on
\code{9*mad(x)}, the width of the weight function is increased such that the percentile of outliers on
that side of the median equals \code{maxPOutliers}. Using \code{maxPOutliers=1} will effectively disable
all weight function broadening; using \code{maxPOutliers=0} will give results that are quite similar (but
not equal to) Pearson correlation. }  

  \item{quick}{ real number between 0 and 1 that controls the handling of missing data in the
calculation of correlations. See details. }

  \item{pearsonFallback}{Specifies whether the bicor calculation should revert to Pearson when median
absolute deviation (mad) is zero. Recongnized values are (abbreviations of) 
\code{"none", "individual", "all"}. If set to
\code{"none"}, zero mad will result in \code{NA} for the corresponding correlation. 
If set to \code{"individual"}, Pearson calculation will be used only for columns that have zero mad.
If set to \code{"all"}, the presence of a single zero mad will cause the whole variable to be treated in
Pearson correlation manner (as if the corresponding \code{robust} option was set to \code{FALSE}). }

  \item{cosine}{ logical: calculate cosine biweight midcorrelation? 
Cosine bicorrelation is similar to standard bicorrelation but the median subtraction is not performed. }

  \item{cosineX}{ logical: use the cosine calculation for \code{x}? This setting does not affect \code{y} 
and can be used to give a hybrid cosine-standard bicorrelation. }

  \item{cosineY}{ logical: use the cosine calculation for \code{y}? This setting does not affect \code{x}
and can be used to give a hybrid cosine-standard bicorrelation. }

  \item{nThreads}{ non-negative integer specifying the number of parallel threads to be used by certain
parts of correlation calculations. This option only has an effect on systems on which a POSIX thread
library is available (which currently includes Linux and Mac OSX, but excludes Windows).
If zero, the number of online processors will be used if it can be determined dynamically, otherwise
correlation calculations will use 2 threads. 
Note that this option does not affect what is usually the most expensive part of the
calculation, namely the matrix multiplication. The matrix multiplication is carried out by BLAS routines provided by R;
these can be sped up by installing a fast BLAS and making R use it.}


  \item{verbose}{ if non-zero, the underlying C function will print some diagnostics.}

  \item{indent}{ indentation for diagnostic messages. Zero means no indentation, each unit adds
two spaces. }

}
\details{
  This function implements biweight midcorrelation calculation (see references). If \code{y} is not
supplied, midcorrelation of columns of \code{x} will be calculated; otherwise, the midcorrelation between
columns of \code{x} and \code{y} will be calculated. Thus, \code{bicor(x)} is equivalent to
\code{bicor(x,x)} but is more efficient. 

The options \code{robustX}, \code{robustY} allow the user to revert the calculation to standard
correlation calculation. This is important, for example, if any of the variables is binary
(or, more generally, discrete) as in such cases the robust methods produce meaningless results. 
If both \code{robustX}, \code{robustY} are set to \code{FALSE}, the function calculates the
standard Pearson correlation (but is slower than the function \code{\link{cor}}).

 The argument \code{quick} specifies the precision of handling of missing data in the correlation
calculations. Value \code{quick = 0} will cause all
calculations to be executed accurately, which may be significantly slower than calculations without
missing data. Progressively higher values will speed up the
calculations but introduce progressively larger errors. Without missing data, all column meadians and
median absolute deviations (MADs) can be pre-calculated before the covariances are calculated. When
missing data are present, 
exact calculations require the column medians and MADs to be calculated for each covariance. The
approximate calculation uses the pre-calculated median and MAD and simply ignores missing data in the
covariance calculation. If the number of missing data is high, the pre-calculated medians and MADs may
be very different from the actual ones, thus potentially introducing large errors.
The \code{quick} value times the
number of rows specifies the maximum difference in the
number of missing entries for median and MAD calculations on the one hand and covariance on the other
hand that will be tolerated before a recalculation is triggered. The hope is that if only a few missing
data are treated approximately, the error introduced will be small but the potential speedup can be
significant.

The choice \code{"all"} for \code{pearsonFallback} is not fully implemented in the sense that there are
rare but possible cases in which the calculation is equivalent to \code{"individual"}. This may happen if
the \code{use} option is set to \code{"pairwise.complete.obs"} and
the missing data are arranged such that each individual mad is non-zero, but when two columns are analyzed
together, the missing data from both columns may make a mad zero. In such a case, the calculation is treated
as Pearson, but other columns will be treated as bicor. 

}
\value{
  A matrix of biweight midcorrelations. Dimnames on the result are set appropriately.
}
\references{ 

Peter Langfelder, Steve Horvath (2012)
Fast R Functions for Robust Correlations and Hierarchical Clustering.
Journal of Statistical Software, 46(11), 1-17.
\url{http://www.jstatsoft.org/v46/i11/}

"Dealing with Outliers in Bivariate Data: Robust Correlation", Rich Herrington, 
http://www.unt.edu/benchmarks/archives/2001/december01/rss.htm

"Introduction to Robust Estimation and Hypothesis Testing", Rand Wilcox, Academic Press, 1997.

"Data Analysis and Regression: A Second Course in Statistics", Mosteller and Tukey, Addison-Wesley,
 1977, pp. 203-209. 
}
\author{ Peter Langfelder}
\keyword{ robust }
