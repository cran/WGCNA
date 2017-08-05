# This function calls the C++ implementation of column quantile.

pmedian = function(...) { pquantile(prob = 0.5, ...)}

colQuantileC = function(data, p)
{  
  data = as.matrix(data)
  storage.mode(data) = "double";
  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!");
  p = as.double(as.character(p));
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!");
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."));

  .Call("quantileC_call", data, p, PACKAGE = "WGCNA");
}

rowQuantileC = function(data, p)
{  
  data = as.matrix(data)
  storage.mode(data) = "double";
  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!");
  ncol = ncol(data);
  nrow = nrow(data);
  quantiles = rep(0, nrow);

  p = as.double(as.character(p));
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!");
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."));

  .Call("rowQuantileC_call", data, p, PACKAGE = "WGCNA");
}

pquantile = function(prob, ...)
{
   pars = list(...)
   pquantile.fromList(pars, prob);
}


pquantile.fromList = function(dataList, prob)
{
   dn = .checkListDimConsistencyAndGetDimnames(dataList);
   if (length(prob) > 1) warning("pquantile2: only the first element of 'prob' will be used.");
   q = .Call("parallelQuantile", dataList, as.numeric(prob[1]));
   dimnames(q) = dn;
   q
}


pmean = function(..., weights = NULL)
{
  pmean.fromList(dataList = list(...), weights = weights)
}

pmean.fromList = function(dataList, weights = NULL)
{
   dn = .checkListDimConsistencyAndGetDimnames(dataList);
   if (is.null(weights)) weights = rep(1, length(dataList))
   q = .Call("parallelMean", dataList, as.numeric(weights));
   dimnames(q) = dn;
   q
}

#pmin.wgcna = function(...)
#{
#  pminWhich.fromList(dataList = list(...))$min
#}

pminWhich.fromList = function(dataList)
{
   dn = .checkListDimConsistencyAndGetDimnames(dataList);
   q = .Call("parallelMin", dataList);
   dimnames(q$min) = dimnames(q$which) = dn;
   q
}

minWhichMin = function(x, byRow = FALSE, dims = 1)
{
  d = dim(x);
  if (length(d) <= 2 && dims==1)
  {
    x = as.matrix(x);
    .Call("minWhich_call", x, as.integer(byRow), PACKAGE = "WGCNA")
  } else {
    if (dims < 1 || dims >= length(d)) stop("Invalid 'dims'. Must be between 1 and length(dim(x))-1.");
    d1 = d[1:dims];
    d2 = d[(dims+1):length(d)];
    dim(x) = c(prod(d1), prod(d2));
    out = .Call("minWhich_call", x, as.integer(byRow));
    if (byRow && length(d1) > 1)
    {
      dim(out$min) = d1;
      dim(out$which) = d1;
    } else if (!byRow && length (d2) > 1)
    {
      dim(out$min) = d2;
      dim(out$which) = d2;
    }
    out;
  }
}


