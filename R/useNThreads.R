# Function to control the number of threads to use in threaded calculations.

.threadAllowVar = "ALLOW_WGCNA_THREADS"

.useNThreads = function(nThreads = 0)
{
  if (nThreads==0)
  {
    nt.env = Sys.getenv(.threadAllowVar, unset = NA);
    if (is.na(nt.env)) return(1);
    if (nt.env=="") return(1);

    if (nt.env=="ALL_PROCESSORS") return (.nProcessorsOnline());

    nt = suppressWarnings(as.numeric(nt.env));

    if (is.na(nt)) return(1);

    return(nt);
  } else
    return (nThreads);
}

.nProcessorsOnline = function()
{
  n = 0;
  res = .C("nProcessorsForR", n = as.integer(n));
  res$n;
}

allowWGCNAThreads = function(nThreads = NULL)
{
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
}
