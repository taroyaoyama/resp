
#' Building stiffness Matrix from stiffnesses of elements
#' for MDOF response calculation
#'
#' @param ks vector of stiffnesses of elements
#'
#' @export

build_k <- function(k) {

    nNodes <- length(k)

    Mat <- matrix(0, nNodes, nNodes)
    submat <- matrix(c(1, -1, -1, 1), 2, 2)

    Mat[1, 1] <- k[1]

    if (nNodes != 1) {

        for (i in 2:nNodes) {

            Mat[(i-1):i, (i-1):i] <- Mat[(i-1):i, (i-1):i] + k[i] * submat
        }
    }

    return (Mat)
}
