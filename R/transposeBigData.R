transposeBigData = function (x, blocksize = 20000) 
{
    isdataframe = is.data.frame(x)
    ismatrix = is.matrix(x)
    if (!(isdataframe | ismatrix)) 
        stop("Input is neither a data frame nor a matrix")
    if (blocksize < 2) 
        stop("This blocksize makes no sense. It should be a positive integer>1.")
    nrow1 = nrow(x)
    ncol1 = ncol(x)
    xTranspose = matrix(NA, nrow = ncol1, ncol = nrow1)
    if (nrow1 <= ncol1) {
        no.blocks = as.integer(ncol1/blocksize)
        if (no.blocks >= 1) {
            for (i in 1:no.blocks) {
                blockIndex = (i - 1) * blocksize + 1:blocksize
                xTranspose[blockIndex, ] = t(x[, blockIndex])
            }
        }
        if (ncol1 - no.blocks * blocksize == 1) {
            xTranspose[ncol1, ] = t(x[, ncol1])
        }
        if (ncol1 - no.blocks * blocksize > 1) {
            finalindex = (no.blocks * blocksize + 1):ncol1
            xTranspose[finalindex, ] = t(x[, finalindex])
        }
    }
    if (nrow1 > ncol1) {
        no.blocks = as.integer(nrow1/blocksize)
        if (no.blocks >= 1) {
            for (i in 1:no.blocks) {
                blockIndex = (i - 1) * blocksize + 1:blocksize
                xTranspose[, blockIndex] = t(x[blockIndex, ])
            }
        }
        if (nrow1 - no.blocks * blocksize == 1) {
            xTranspose[, nrow1] = t(x[nrow1, ])
        }
        if (nrow1 - no.blocks * blocksize > 1) {
            finalindex = (no.blocks * blocksize + 1):nrow1
            xTranspose[, finalindex] = t(x[finalindex, ])
        }
    }
    if (isdataframe) {
        xTranspose = data.frame(xTranspose)
        dimnames(xTranspose)[[1]] = dimnames(x)[[2]]
        dimnames(xTranspose)[[2]] = dimnames(x)[[1]]
    }
    if (ismatrix) {
        dimnames(xTranspose)[[1]] = dimnames(x)[[2]]
        dimnames(xTranspose)[[2]] = dimnames(x)[[1]]
    }
    xTranspose
}
