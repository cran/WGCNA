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
  err = as.integer(nNA-1 + 1/1);
  warnX = as.integer(1L- 1/1)
  warnY = as.integer(2L- 1/1 - 3/3);
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

cor = function(x, y = NULL, use = "all.obs", method = c("pearson", "kendall", "spearman"),
               weights.x = NULL, weights.y = NULL,
               quick = 0, 
               cosine = FALSE, 
               cosineX = cosine, cosineY = cosine,
               drop = FALSE,
               nThreads = 0, verbose = 0, indent = 0)
{
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
        "everything", "na.or.complete"), nomatch = 0)
    method <- match.arg(method)

    if (length(weights.x)==0) weights.x = NULL;
    if (length(weights.y)==0) weights.y = NULL;

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

      if (!is.null(weights.x))
      {
        if (is.null(dim(weights.x)))
        {
          if (length(weights.x)!=nrow(x))
            stop("When 'weights.x' are given, they must be a vector of length 'nrow(x)' or a matrix\n",
                 "of the same dimensions as 'x'.");
          weights.x = matrix(weights.x, nrow(x), ncol(x));
        } else
          if (!isTRUE(all.equal(dim(weights.x), dim(x))))
             stop("When 'weights.x' are given, they must be a vector of length 'nrow(x)' or a matrix\n",
                 "of the same dimensions as 'x'.");
        if (any(!is.finite(weights.x))) 
        {
          if (verbose > 0)
            warning("cor: found non-finite weights. These will be removed (set to missing), ",
                    "and the corresponding entries in 'x' will be treated as missing.");
          weights.x[!is.finite(weights.x)] = NA;
        }
        if (any(weights.x < 0, na.rm = TRUE))
          stop("All weights must be non-negative.");

        if (!is.null(y) && is.null(weights.y)) weights.y = matrix(1, nrow(y), ncol(y));
      }

      if (!is.null(weights.y))
      {
        if (is.null(y)) stop("'weights.y' can only be used if 'y' is non-NULL.");
        if (is.null(dim(weights.y)))
        {
          if (length(weights.y)!=nrow(y))
            stop("When 'weights.y' are given, they must be a vector of length 'nrow(y)' or a matrix\n",
                 "of the same dimensions as 'y'.");
          weights.y = matrix(weights.y, nrow(y), ncol(y));
        } else
          if (!isTRUE(all.equal(dim(weights.y), dim(y))))
             stop("When 'weights.y' are given, they must be a vector of length 'nrow(y)' or a matrix\n",
                 "of the same dimensions as 'y'.");
        if (any(!is.finite(weights.y)))
        {
          if (verbose > 0) 
            warning("cor: found non-finite weights. These will be removed (set to missing), ",
                    "and the corresponding entries in 'x' will be treated as missing.");
          weights.y[!is.finite(weights.y)] = NA;
        }
        if (any(weights.y < 0, na.rm = TRUE))
          stop("All weights must be non-negative.");

        if (is.null(weights.x)) weights.x = matrix(1, nrow(x), ncol(x));
      }

      storage.mode(x)= "double";
      if (!is.null(weights.x)) storage.mode(weights.x) = "double";
      if (!is.null(weights.y)) storage.mode(weights.y) = "double";
      nNA = 0L
      err = as.integer(nNA-1 + 1/1);
      cosineX = as.integer(cosineX);
      nThreads = as.integer(nThreads);
      verbose = as.integer(verbose);
      indent = as.integer(indent);
      if (is.null(y))
      {
         res = .Call("cor1Fast_call", x, weights.x,
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
                 weights.x, weights.y,
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


