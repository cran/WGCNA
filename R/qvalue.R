# qvalue function by John D. Storey <email: jstorey@u.washington.edu>, modified by Peter Langfelder 
# for use in WGCNA


qvalue <- function(p, lambda=seq(0,0.90,0.05), pi0.method="smoother", fdr.level=NULL, robust=FALSE, 
  smooth.df = 3, smooth.log.pi0 = FALSE) {
#Input
#=============================================================================
#p: a vector of p-values (only necessary input)
#fdr.level: a level at which to control the FDR (optional)
#lambda: the value of the tuning parameter to estimate pi0 (optional)
#pi0.method: either "smoother" or "bootstrap"; the method for automatically
#           choosing tuning parameter in the estimation of pi0, the proportion
#           of true null hypotheses
#robust: an indicator of whether it is desired to make the estimate more robust
#        for small p-values and a direct finite sample estimate of pFDR (optional)
#gui: A flag to indicate to 'qvalue' that it should communicate with the gui.  ## change by Alan
#     Should not be specified on command line.
#smooth.df: degrees of freedom to use in smoother (optional)
#smooth.log.pi0: should smoothing be done on log scale? (optional)
#
#Output
#=============================================================================
#call: gives the function call
#pi0: an estimate of the proportion of null p-values
#qvalues: a vector of the estimated q-values (the main quantity of interest)
#pvalues: a vector of the original p-values
#significant: if fdr.level is specified, an indicator of whether the q-value
#    fell below fdr.level (taking all such q-values to be significant controls
#    FDR at level fdr.level)

#This is just some pre-processing

    if(min(p)<0 || max(p)>1) 
        stop("qvalue: p-values not in valid range.")
    if(length(lambda)>1 && length(lambda)<4) 
        stop("qvalue: If length of lambda greater than 1, you need at least 4 values.")
    if(length(lambda)>1 && (min(lambda) < 0 || max(lambda) >= 1)) 
        stop("qvalue: Lambda must be within [0, 1).")
    m <- length(p)
#These next few functions are the various ways to estimate pi0
    if(length(lambda)==1) {
        if(lambda<0 || lambda>=1) 
            stop("qvalue: Lambda must be within [0, 1).")

        pi0 <- mean(p >= lambda)/(1-lambda)
        pi0 <- min(pi0,1)
    }
    else {
        pi0 <- rep(0,length(lambda))
        for(i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
        }

        if(pi0.method=="smoother") {
            if(smooth.log.pi0)
              pi0 <- log(pi0)

            spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
            pi0 <- predict(spi0,x=max(lambda))$y

            if(smooth.log.pi0)
              pi0 <- exp(pi0)
            pi0 <- min(pi0,1)
        }
        else if(pi0.method=="bootstrap") {
            minpi0 <- min(pi0)
            mse <- rep(0,length(lambda))
            pi0.boot <- rep(0,length(lambda))
            for(i in 1:100) {
                p.boot <- sample(p,size=m,replace=TRUE)
                for(i in 1:length(lambda)) {
                    pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
                }
                mse <- mse + (pi0.boot-minpi0)^2
            }
            pi0 <- min(pi0[mse==min(mse)])
            pi0 <- min(pi0,1)
        }
        else {  ## change by Alan: check for valid choice of 'pi0.method' (only necessary on command line)
            stop("qvalue:: 'pi0.method' must be one of 'smoother' or 'bootstrap'.")
            return(0)
        }
    }
    if(pi0 <= 0) 
        stop("qvalue:: The estimated pi0 <= 0. Check that you have valid p-values or use another lambda method.")
    if(!is.null(fdr.level) && (fdr.level<=0 || fdr.level>1))  ## change by Alan:  check for valid fdr.level
        stop("qvalue:: 'fdr.level' must be within (0, 1].")
#The estimated q-values calculated here
    u <- order(p)

    # change by Alan
    # ranking function which returns number of observations less than or equal
    qvalue.rank <- function(x) {
      idx <- sort.list(x)

      fc <- factor(x)
      nl <- length(levels(fc))
      bin <- as.integer(fc)
      tbl <- tabulate(bin)
      cs <- cumsum(tbl)
 
      tbl <- rep(cs, tbl)
      tbl[idx] <- tbl

      return(tbl)
    }

    v <- qvalue.rank(p)
    
    qvalue <- pi0*m*p/v
    if(robust) {
        qvalue <- pi0*m*p/(v*(1-(1-p)^m))
    }
    qvalue[u[m]] <- min(qvalue[u[m]],1)
    for(i in (m-1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]],qvalue[u[i+1]],1)
    }
#The results are returned
    if(!is.null(fdr.level)) {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, fdr.level=fdr.level, ## change by Alan
          significant=(qvalue <= fdr.level), lambda=lambda)
    }
    else {
        retval <- list(call=match.call(), pi0=pi0, qvalues=qvalue, pvalues=p, lambda=lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
}
