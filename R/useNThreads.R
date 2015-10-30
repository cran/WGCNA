# Function to control the number of threads to use in threaded calculations.

.useNThreads = function(nThreads = 0)
{
  if (nThreads==0)
  {
    nt.env = Sys.getenv(.threadAllowVar, unset = NA);
    if (is.na(nt.env)) return(1);
    if (nt.env=="") return(1);

    if (nt.env=="ALL_PROCESSORS") return (.nProcessorsOnline());

    nt = suppressWarnings(as.numeric(nt.env));

    if (!is.finite(nt)) return(2);

    return(nt);
  } else
    return (nThreads);
}

.nProcessorsOnline = function()
{
  n = detectCores();
  if (!is.numeric(n)) n = 2;
  if (!is.finite(n)) n = 2;
  if (n<1) n = 2;
  n;
}

allowWGCNAThreads = function(nThreads = NULL)
{
  # Stop any clusters that may be still running
  disableWGCNAThreads()
  # Enable WGCNA threads
  if (is.null(nThreads)) nThreads = .nProcessorsOnline();
  if (!is.numeric(nThreads) || nThreads < 2)
    stop("nThreads must be numeric and at least 2.");

  if (nThreads > .nProcessorsOnline())
    printFlush(paste("Warning in allowWGCNAThreads: Requested number of threads is higher than number\n",
                     "of available processors (or cores). Using too many threads may degrade code",
                     "performance. It is recommended that the number of threads is no more than number\n",
                     "of available processors.\n"))

  printFlush(paste("Allowing multi-threading with up to", nThreads, "threads."));

  pars = list(nThreads);
  names(pars) = .threadAllowVar;
  do.call(Sys.setenv, pars);
  invisible(nThreads);
}

disableWGCNAThreads = function()
{
  Sys.unsetenv(.threadAllowVar);
  pars = list(1)
  names(pars) = .threadAllowVar
  do.call(Sys.setenv, pars)
  if (exists(".revoDoParCluster", where = ".GlobalEnv")) 
  {
    stopCluster(get(".revoDoParCluster", pos = ".GlobalEnv"));
  }
  registerDoSEQ();
}

.checkAvailableMemory = function()
{
  size = 0;
  res = .C("checkAvailableMemoryForR", size = as.double(size), PACKAGE = "WGCNA")
  res$size;
}

# Function to calculate an appropriate blocksize

blockSize = function(matrixSize, rectangularBlocks = TRUE, maxMemoryAllocation = NULL, overheadFactor = 3)
{
  if (is.null(maxMemoryAllocation))
  {
    maxAlloc = .checkAvailableMemory();
  } else {
    maxAlloc = maxMemoryAllocation/8;
  }
  maxAlloc = maxAlloc/overheadFactor;

  if (rectangularBlocks)
  {
     blockSz = floor(maxAlloc/matrixSize);
  } else
     blockSz = floor(sqrt(maxAlloc));

  return( min (matrixSize, blockSz) )
}

#======================================================================================================
#
# enableWGCNAThreads
#
#======================================================================================================

enableWGCNAThreads = function(nThreads = NULL)
{
  nCores = detectCores();
  if (is.null(nThreads))
  {
    if (nCores < 4) nThreads = nCores else nThreads = nCores - 1;
  }
  if (!is.numeric(nThreads) || nThreads < 2)
        stop("nThreads must be numeric and at least 2.")
  if (nThreads > nCores)
     printFlush(paste("Warning in allowWGCNAThreads: Requested number of threads is higher than number\n",
          "of available processors (or cores). Using too many threads may degrade code",
          "performance. It is recommended that the number of threads is no more than number\n",
          "of available processors.\n"))
  printFlush(paste("Allowing parallel execution with up to", nThreads, "working processes."))
  pars = list(nThreads)
  names(pars) = .threadAllowVar
  do.call(Sys.setenv, pars)

  # Register a parallel backend for foreach
  registerDoParallel(nThreads);

  # Return the number of threads invisibly
  invisible(nThreads)
}

WGCNAnThreads = function()
{
  n = suppressWarnings(as.numeric(as.character(Sys.getenv(.threadAllowVar, unset = 1))));
  if (is.na(n)) n = 1;
  if (length(n)==0) n = 1;
  n;
}

#========================================================================================================
#
# allocateJobs
#
#========================================================================================================

# Facilitates multi-threading by producing an even allocation of jobs 
# Works even when number of jobs is less than number of threads in which case some components of the
# returned allocation will have length 0.

allocateJobs = function(nTasks, nWorkers)
{
  if (is.na(nWorkers))
  {
    warning("In function allocateJobs: 'nWorkers' is NA. Will use 1 worker.");
    nWorkers = 1;
  }
  n1 = floor(nTasks/nWorkers);
  n2 = nTasks - nWorkers*n1;
  allocation = list();
  start = 1;
  for (t in 1:nWorkers)
  {
    end = start + n1 - 1 + as.numeric(t<=n2);
    if (start > end)
    {
      allocation[[t]] = numeric(0);
    } else allocation[[t]] = c(start:end);
    start = end+1;
  }

  allocation;
}


  
