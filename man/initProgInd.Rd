\name{Inline display of progress}
\alias{initProgInd}
\alias{updateProgInd}
\title{ Inline display of progress }
\description{
  These functions provide an inline display of pregress. 
}
\usage{
initProgInd(leadStr = "..", trailStr = "", quiet = !interactive())
updateProgInd(newFrac, progInd, quiet = !interactive())
}
\arguments{
  \item{leadStr}{ character string that will be printed before the actual progress number. }
  \item{trailStr}{ character string that will be printed after the actual progress number. }
  \item{quiet}{ can be used to silence the indicator for non-interactive sessions whose output is
typically redirected to a file. }
  \item{newFrac}{ new fraction of progress to be displayed. }
  \item{progInd}{ an object of class \code{progressIndicator} that encodes previously printed message. }
}
\details{
  A progress indicator is a simple inline display of progress intended to satisfy impatient users during
lengthy operations. The function \code{initProgInd} initializes a progress indicator (at zero);
\code{updateProgInd} updates it to a specified fraction. 

Note that excessive use of \code{updateProgInd} may lead to a performance penalty (see examples). 
}
\value{
  Both functions return an object of class \code{progressIndicator} that holds information on the last
printed value and should be used for subsequent updates of the indicator. 
}
\author{ Peter Langfelder }
\examples{
max = 10;
prog = initProgInd("Counting: ", "done");
for (c in 1:max)
{
  Sys.sleep(0.10);
  prog = updateProgInd(c/max, prog);
}
printFlush("");

printFlush("Example 2:");
prog = initProgInd();
for (c in 1:max)
{
  Sys.sleep(0.10);
  prog = updateProgInd(c/max, prog);
}
printFlush("");

## Example of a significant slowdown:

## Without progress indicator:

system.time( {a = 0; for (i in 1:10000) a = a+i; } )

## With progress indicator, some 50 times slower:

system.time( 
  {
    prog = initProgInd("Counting: ", "done");
    a = 0; 
    for (i in 1:10000) 
    {
      a = a+i; 
      prog = updateProgInd(i/10000, prog);
    }
  }   
)
}

\keyword{misc}% __ONLY ONE__ keyword per line
