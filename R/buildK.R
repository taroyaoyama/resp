
## build stiffness Matrix
buildK <- function(k) {
    nNodes <- length(k)
    Mat <- matrix(0, nNodes, nNodes)
    mat <- matrix(c(1, -1, -1, 1), 2, 2)
    Mat[1, 1] <- k[1]
    if (nNodes != 1) {
        for (i in 2:nNodes) {
            Mat[(i-1):i, (i-1):i] <- Mat[(i-1):i, (i-1):i] + k[i] * mat
        }
    }
    return (Mat)
}
