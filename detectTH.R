#' moran.I
#'
#' calculate the Moran's I
#'
#' @param x the variable upon which the Moran's I is calculated
#' @param weight the weight matrix
#' @return a list that contains the Moran's I and its p-value
#' @noRd


moran.I <- function(x, weight, na.rm = FALSE){


    if (dim(weight)[1] != dim(weight)[2])
        stop("'weight' must be a square matrix")
    n <- length(x)
    if (dim(weight)[1] != n)
        stop("'weight' must have as many rows as observations in 'x'")
    ei <- -1/(n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        }
        else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(list(observed = NA, expected = ei, sd = NA,
                        p.value = NA))
        }
    }
    ROWSUM <- rowSums(weight)
    ROWSUM[ROWSUM == 0] <- 1
    weight <- weight/ROWSUM
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    obs <- (n/s) * (cv/v)


    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    s.sq <- s^2
    k <- (sum(y^4)/n)/(v/n)^2

    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) -
                     k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/
                    ((n - 1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))


    pv <- pnorm(obs, mean = ei, sd = sdi, lower.tail = FALSE)

    list(obs = obs, p.value = pv)
}




#' detect.pv
#'
#' Returns the potential partitioning variable(s)
#'
#' @param membership a length-n vector of membership labels
#' @param X.p A subset of covariates to be examined
#' @param size the maximum number of threshold variables allowed in the model; default size = the total number of covariates provided
#' @param No.return number of results to be returned; default is the total number of threshold variables explored
#' @return a list consisting of two items: the first item is the p values of Moran's I, ranked in an ascending order; the second item contains the corresponding covariates that yield those p values
#' @examples
#'
#' set.seed(121)
#'
#' n <- 50
#' p <- 3
#'
#' X <- matrix(rnorm(n * p), nrow = n)
#' Xj <- X[,1] # the threshold variable
#' beta1 <- rep(3,p)
#' beta2 <- rep(-3,p)
#'
#' index.g1 <- which(Xj <= 0)
#' index.g2 <- which(Xj > 0)
#'
#' y.g1 <- X[index.g1,] %*% beta1
#' y.g2 <- X[index.g2,] %*% beta2
#'
#' y <- rep(0,n)
#' y[index.g1] <- y.g1
#' y[index.g2] <- y.g2
#'
#' y <- y + rnorm(n = n, sd = 0.5)
#' m <- HP(X,y, method = "latent", max.no.cluster = 2)$membership
#' detect.pv(m, X.p = X)
#'
#'@export detect.pv


detect.pv <- function(membership, X.p, size = 1, No.return){

    index.list <- list()
    P.val <- c()
    k <- 1 # list object counter

    X.p <- as.matrix(X.p)

    subsize <- size

    while (subsize >= 1){

        index <- utils::combn(ncol(X.p), subsize)

        no.subsets <- ncol(index)

        for ( i in 1:no.subsets){

            x.p <- X.p[ ,index[,i] ]

            dist.xp <- dist(x.p)

            weight <- mst(dist.xp)

            res <- moran.I(membership, weight)

            P.val <- c(P.val, res$p.value)

            index.list[[k]]<- index[,i]

            k <- k + 1

        }
        subsize <- subsize - 1
    }

    index.sorted <- sort(P.val, index.return = T)$ix

    if (missing(No.return)){

        No.return <- k-1
    }

    index.max <- index.sorted[1:No.return]

    index.list.sorted <- index.list[index.max]
    P.val.sorted <- P.val[index.max]
    res.list <- list(index = index.list.sorted,
                     P.value = P.val.sorted)

    return(res.list)

}






