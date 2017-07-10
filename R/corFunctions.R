# slight re-definition of the bicor function

bicor = function(x, y = NULL, robustX = TRUE, robustY = TRUE, use = 'all.obs', maxPOutliers = 1, quick = 0,
                 pearsonFallback = "individual", 
                 cosine = FALSE,
                 cosineX = cosine, cosineY = cosine,
                 nThreads = 0, verbose = 0, indent = 0)
{
  Cerrors = c("Memory allocation error")
  nKnownErrors = length(Cerrors);
  na.method = pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method))
      stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
          "'all.obs', 'pairwise.complete.obs'"))
  if (na.method==1)
  {
    if (sum(is.na(x))> 0)
      stop("Missing values present in input variable 'x'. Consider using use = 'pairwise.complete.obs'.");
    if (!is.null(y))
    {
      if (sum(is.na(y)) > 0)
        stop("Missing values present in input variable 'y'. Consider using use = 'pairwise.complete.obs'.");
    }
  }

  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(na.method))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))

  if (quick < 0) stop("quick must be non-negative.");
  if (nThreads < 0) stop("nThreads must be non-negative.");
  if (is.null(nThreads) || (nThreads==0)) nThreads = .useNThreads();

  x = as.matrix(x);
  if (prod(dim(x))==0) stop("'x' has a zero dimension."); 
  storage.mode(x) = "double";
  nNA = 0L;
  err = 0L;
  warnX = 0L;
  warnY = 0L;
  quick = as.double(quick);
  maxPOutliers = as.double(maxPOutliers);
  fallback = as.integer(fallback);
  cosineX = as.integer(cosineX);
  robustX = as.integer(robustX);
  nThreads = as.integer(nThreads);
  verbose = as.integer(verbose); indent = as.integer(indent)
  if (is.null(y))
  {
    if (!robustX)
    {
      res = cor(x, use = use)
    } else {
      res = .Call("bicor1_call", x, 
               maxPOutliers, 
               quick, 
               fallback,
               cosineX, 
               nNA, err, warnX, 
               nThreads, verbose, indent,
               PACKAGE = "WGCNA");
    }
    if (!is.null(colnames(x))) dimnames(res) = list(colnames(x),  colnames(x));
    if (warnX > 0)
    {
      # For now have only one warning
      warning(paste("bicor: zero MAD in variable 'x'.", .zeroMADWarnings[fallback]));
    }
  } else {
    y = as.matrix(y);
    storage.mode(y) = "double";
    if (prod(dim(y))==0) stop("'y' has a zero dimension."); 
    if (nrow(x)!=nrow(y))
      stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
    cosineY = as.integer(cosineY);
    robustY = as.integer(robustY);
    res = .Call("bicor2_call", x, y,
             robustX, robustY,
             maxPOutliers, 
             quick, 
             fallback,
             cosineX,
             cosineY,
             nNA, err,
             warnX, warnY, 
             nThreads,
             verbose, indent,
             PACKAGE = "WGCNA");
    if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
        dimnames(res) = list(dimnames(x)[[2]], dimnames(y)[[2]]);
    if (warnX > 0)
      warning(paste("bicor: zero MAD in variable 'x'.", .zeroMADWarnings[fallback]));
    if (warnY > 0)
      warning(paste("bicor: zero MAD in variable 'y'.", .zeroMADWarnings[fallback]));
  }
  if (err > 0)
  {
    if (err > nKnownErrors)
    {
      stop(paste("An error occurred in compiled code. Error code is", err));
    } else {
      stop(paste(Cerrors[err], "occurred in compiled code. "));
    }
  }
  if (nNA > 0)
  {
    warning(paste("Missing values generated in calculation of bicor.",
                  "Likely cause: too many missing entries, zero median absolute deviation, or zero variance."));
  }
  res;
}

# Code to call my implementation of correlation
# For less than 100 correlations, use stats::cor since that is usually faster, particularly when no missing
# data are present, likely due to the complicated threading I do in the WGCNA correlations.  

cor = function(x, y = NULL, use = "all.obs", method = c("pearson", "kendall", "spearman"),
               quick = 0, 
               cosine = FALSE, 
               cosineX = cosine, cosineY = cosine,
               drop = FALSE,
               nThreads = 0, verbose = 0, indent = 0)
{
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
        "everything", "na.or.complete"), nomatch = 0)
    method <- match.arg(method)

    x = as.matrix(x);
    nx = ncol(x);
    if (!is.null(y)) 
    {
      y = as.matrix(y);
      ny = ncol(y);
    } else ny = nx;

    if ((method=="pearson") && ( (na.method==1) || (na.method==3) ))
    {
      Cerrors = c("Memory allocation error")
      nKnownErrors = length(Cerrors);
      na.method = pmatch(use, c("all.obs", "pairwise.complete.obs"))
      if (is.na(na.method))
          stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
              "'all.obs', 'pairwise.complete.obs'"))
      if (na.method==1)
      {
         if (sum(is.na(x))> 0)
           stop("Missing values present in input variable 'x'. Consider using use = 'pairwise.complete.obs'.");
         if (!is.null(y))
         {
           if (sum(is.na(y)) > 0)
             stop("Missing values present in input variable 'y'. Consider using use = 'pairwise.complete.obs'.");
         }
      }
      
      if (quick < 0) stop("quick must be non-negative.");
      if (nThreads < 0) stop("nThreads must be non-negative.");
      if (is.null(nThreads) || (nThreads==0)) nThreads = .useNThreads();
 
      if (prod(dim(x))==0) stop("'x' has a zero dimension."); 
      storage.mode(x)= "double";
      nNA = 0L
      err = 0L
      cosineX = as.integer(cosineX);
      nThreads = as.integer(nThreads);
      verbose = as.integer(verbose);
      indent = as.integer(indent);
      if (is.null(y))
      {
         res = .Call("cor1Fast_call", x,
                   quick, cosine,
                   nNA, err, nThreads,
                   verbose, indent, PACKAGE = "WGCNA");
         if (!is.null(dimnames(x)[[2]])) dimnames(res) = list(dimnames(x)[[2]],  dimnames(x)[[2]] );
      } else {
         y = as.matrix(y);
         storage.mode(y)= "double";
         cosineY = as.integer(cosineY);
         if (prod(dim(y))==0) stop("'y' has a zero dimension."); 
         if (nrow(x)!=nrow(y))
            stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
         res = .Call("corFast_call", x, y, 
                 quick, 
                 cosineX, 
                 cosineY,
                 nNA, err,
                 nThreads,
                 verbose, indent, 
                 PACKAGE = "WGCNA");
         if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
            dimnames(res) = list(dimnames(x)[[2]], dimnames(y)[[2]]);
      }
      if (err > 0)
      {
        if (err > nKnownErrors)
        {
          stop(paste("An error occurred in compiled code. Error code is", err));
        } else {
          stop(paste(Cerrors[err], "occurred in compiled code. "));
        }
      }
      if (nNA > 0)
      {
        warning(paste("Missing values generated in calculation of cor.",
                      "Likely cause: too many missing entries or zero variance."));
      }
      if (drop) res[, , drop = TRUE] else res;
    } else {
      stats::cor(x,y, use, method);
    }
}


# Wrappers for compatibility with older scripts

cor1 = function(x, use = "all.obs", verbose = 0, indent = 0) 
{
   cor(x, use = use, verbose = verbose, indent = indent);
}

corFast = function(x, y = NULL, use = "all.obs",
                quick = 0, nThreads = 0, verbose = 0, indent = 0)
{
  cor(x,y, use, method = "pearson", quick, nThreads, verbose, indent)
}


      
