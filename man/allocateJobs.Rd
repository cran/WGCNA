\name{allocateJobs}
\alias{allocateJobs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Divide tasks among workers
}
\description{
This function calculates an even splitting of a given number of tasks among a given number of workers
(threads).
}
\usage{
allocateJobs(nTasks, nWorkers)
}
\arguments{
  \item{nTasks}{
number of tasks to be divided
}
  \item{nWorkers}{
number of workers
}
}
\details{
Tasks are labeled consecutively 1,2,..., \code{nTasks}. The tasks are split in contiguous blocks as evenly
as possible. 
}
\value{
A list with one component per worker giving the task indices to be worked on by each worker. If there are
more workers than tasks, the tasks for the extra workers are 0-length numeric vectors. 
}
\author{
Peter Langfelder
}
\examples{
allocateJobs(10, 3);
allocateJobs(2,4);
}
\keyword{misc}% __ONLY ONE__ keyword per line
