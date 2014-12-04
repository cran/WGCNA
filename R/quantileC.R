# This function calls the C++ implementation of column quantile.

colQuantileC = function(data, p)
{  
  data = as.matrix(data)

  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!");
  ncol = ncol(data);
  nrow = nrow(data);
  quantiles = rep(0, ncol);

  p = as.numeric(as.character(p));
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!");
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."));

  res = .C("quantileC", data = as.double(data), nrow = as.integer(nrow), ncol = as.integer(ncol), 
           p = as.double(p), quantiles = as.double(quantiles), NAOK = TRUE);

  res$quantiles;
}

rowQuantileC = function(data, p)
{  
  data = as.matrix(data)

  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!");
  ncol = ncol(data);
  nrow = nrow(data);
  quantiles = rep(0, nrow);

  p = as.numeric(as.character(p));
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!");
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."));

  res = .C("rowQuantileC", data = as.double(data), nrow = as.integer(nrow), ncol = as.integer(ncol), 
           p = as.double(p), quantiles = as.double(quantiles), NAOK = TRUE);

  res$quantiles;
}
