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
  nNA = 0;
  err = 0;
  warnX = 0;
  warnY = 0;
  if (is.null(y))
  {
    if (!robustX)
    {
      bi = cor(x, use = use)
    } else {
      bi = matrix(0, ncol(x), ncol(x));
      res = .C("bicor1Fast", x = as.double(x), nrow = as.integer(nrow(x)), ncol = as.integer(ncol(x)),
               maxPOutliers = as.double(maxPOutliers), 
               quick = as.double(quick), 
               fallback = as.integer(fallback),
               cosine = as.integer(cosineX), 
               res = as.double(bi), nNA = as.integer(nNA),
               err = as.integer(err), 
               warn = as.integer(warnX), nThreads = as.integer(nThreads),
               verbose = as.integer(verbose), indent = as.integer(indent),
               NAOK = TRUE, PACKAGE = "WGCNA");
    }
    dim(res$res) = dim(bi);
    if (!is.null(dimnames(x)[[2]])) dimnames(res$res) = list(dimnames(x)[[2]],  dimnames(x)[[2]] );
    if (res$warn > 0)
    {
      # For now have only one warning
      warning(paste("bicor: zero MAD in variable 'x'.", .zeroMADWarnings[fallback]));
    }
  } else {
    y = as.matrix(y);
    if (prod(dim(y))==0) stop("'y' has a zero dimension."); 
    if (nrow(x)!=nrow(y))
      stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
    bi = matrix(0, ncol(x), ncol(y));
    res = .C("bicorFast", x = as.double(x), nrow = as.integer(nrow(x)), ncolx = as.integer(ncol(x)),
             y = as.double(y), ncoly = as.integer(ncol(y)),
             robustX = as.integer(robustX), robustY = as.integer(robustY),
             maxPOutliers = as.double(maxPOutliers), 
             quick = as.double(quick), 
             fallback = as.integer(fallback),
             cosineX = as.integer(cosineX),
             cosineY = as.integer(cosineY),
             res = as.double(bi), nNA = as.integer(nNA), err = as.integer(err),
             warnX = as.integer(warnX), 
             warnY = as.integer(warnY), 
             nThreads = as.integer(nThreads),
             verbose = as.integer(verbose), indent = as.integer(indent), NAOK = TRUE,
             PACKAGE = "WGCNA");
    dim(res$res) = dim(bi);
    if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
        dimnames(res$res) = list(dimnames(x)[[2]], dimnames(y)[[2]]);
    if (res$warnX > 0)
      warning(paste("bicor: zero MAD in variable 'x'.", .zeroMADWarnings[fallback]));
    if (res$warnY > 0)
      warning(paste("bicor: zero MAD in variable 'y'.", .zeroMADWarnings[fallback]));
  }
  if (res$err > 0)
  {
    if (err > nKnownErrors)
    {
      stop(paste("An error occurred in compiled code. Error code is", err));
    } else {
      stop(paste(Cerrors[err], "occurred in compiled code. "));
    }
  }
  if (res$nNA > 0)
  {
    warning(paste("Missing values generated in calculation of bicor.",
                  "Likely cause: too many missing entries, zero median absolute deviation, or zero variance."));
  }
  res$res;
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
      nNA = as.integer(0);
      err = as.integer(0);
      cosine = as.integer(cosine);
      nThreads = as.integer(nThreads);
      verbose = as.integer(verbose);
      indent = as.integer(indent);
      if (is.null(y))
      {
         res = .Call("cor1Fast_call", x,
                   quick, cosine,
                   nNA, err, nThreads,
                   verbose, indent, package = "WGCNA");
         if (!is.null(dimnames(x)[[2]])) dimnames(res) = list(dimnames(x)[[2]],  dimnames(x)[[2]] );
      } else {
         y = as.matrix(y);
         if (prod(dim(y))==0) stop("'y' has a zero dimension."); 
         if (nrow(x)!=nrow(y))
            stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
         bi = matrix(0, ncol(x), ncol(y));
         res = .C("corFast", x = as.double(x), nrow = as.integer(nrow(x)), ncolx = as.integer(ncol(x)),
                 y = as.double(y), ncoly = as.integer(ncol(y)),
                 quick = as.double(quick), 
                 cosineX = as.integer(cosineX), 
                 cosineY = as.integer(cosineY),
                 res = as.double(bi), nNA = as.integer(nNA), err = as.integer(err),
                 nThreads = as.integer(nThreads),
                 verbose = as.integer(verbose), indent = as.integer(indent), NAOK = TRUE,
                 PACKAGE = "WGCNA");
         res = res$res;
         dim(res) = dim(bi);
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


      
