# Parallel version of quantile, mean and median

.dimensions = function(x)
{
   if (is.null(dim(x))) return(length(x))
   return(dim(x))
}

.shiftList = function(c0, lst)
{
  if (length(lst)==0) return(list(c0));
  ll = length(lst)
  out = list();
  out[[1]] = c0;
  for (i in 1:ll)
    out[[i+1]] = lst[[i]];
  out;
}

pquantile = function(prob, ...)
{
   pars = list(...)
   nPars = length(pars)
   dn = NULL
   for (p in 1:nPars)
   {
       if (mode(pars[[p]])!="numeric") 
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dim)))
             stop("Argument dimensions are not consistent.");
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."));
       if (is.null(dn)) dn = dimnames(pars[[p]]);
       pars[[p]] = as.numeric(pars[[p]]);
   }
   x = as.matrix(as.data.frame(pars))
   if (any(is.na(x))) 
      warning("The input contains missing data that will be removed.")
   q = apply(x, 1, quantile, prob = prob, na.rm = TRUE);
   rnq = rownames(q);
   if (length(dim) > 1) dim(q) = (if (length(prob)==1) dim else c(length(prob), dim));
   if (!is.null(dn)) dimnames(q) = (if (length(prob)==1) dn else .shiftList(rnq, dn));
   q
}

pmedian = function(...) { pquantile(0.5, ...)}

pmean = function(...)
{
   pars = list(...)
   nPars = length(pars)
   dn = NULL;
   for (p in 1:nPars)
   {
       if (mode(pars[[p]])!="numeric") 
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dim)))
             stop("Argument dimensions are not consistent.");
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."));
       if (is.null(dn)) dn = dimnames(pars[[p]]);
       pars[[p]] = as.numeric(pars[[p]]);
   }
   x = as.matrix(as.data.frame(pars))
   if (any(is.na(x))) 
      warning("The input contains missing data that will be removed.")
   q = rowMeans(x, na.rm = TRUE);
   if (length(dim) > 1) dim(q) = dim;
   if (!is.null(dn)) dimnames(q) = dn;
   q
}

