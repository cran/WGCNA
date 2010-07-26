\name{metaZfunction}
\alias{metaZfunction}
\title{
  Meta-analysis Z statistic
}
\description{
The function calculates a meta analysis Z statistic based on an input data frame of Z statistics.
}
\usage{
metaZfunction(datZ, columnweights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datZ}{
Matrix or data frame of Z statistics (assuming standard normal distribution under the null hypothesis). Rows
correspond to genes, columns to independent data sets.  }

  \item{columnweights}{
optional vector of non-negative numbers for weighing the columns of datZ.  }
}
\details{
For example, if datZ has 3 columns whose columns are labelled Z1,Z2,Z3  then ZMeta= (Z1+Z2+Z3)/sqrt(3).
Under the null hypothesis (where all Z statistics follow a standard normal distribution and the Z statistics
are independent), ZMeta also follows a standard normal distribution.  To calculate a 2 sided p-value, one an
use the following code
pvalue=2*pnorm(-abs(ZMeta) ) 
}
\value{
Vector of meta analysis Z statistic. Under the null hypothesis this should follow a standard normal
distribution. 
}
\author{
Steve Horvath
}
\keyword{misc}
