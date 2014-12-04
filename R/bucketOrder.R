.bucketOrder = function(data, min = min(data), max = max(data), nIntervals, exact = FALSE)
{
  
  if (any(!is.finite(data))) stop("'data' cannot contain missing or infinite values.");
  .Call("bucketOrder_R", as.numeric(data), as.double(min), as.double(max), as.integer(nIntervals),
        as.integer(exact));
  
}
