#' MLE-Based Test
#'
#' Method for testing difference of means in partially paired data using
#' maximum likelihood estimation under either the assumption of
#' heteroscedasticity or homoscedasticity as proposed by Lin and Stivers [1974]
#' or Ekbohm [1976], respectively.
#' 
#' @name mle.test
#'
#' @param x A numeric vector, or 2xn (or nx2) matrix. If \code{x} is a vector,
#' then data should be paired with \code{y} by index. If \code{x} is a matrix,
#' then pairing is done within column. The data must contain at least two
#' pairwise complete observations, and at least four pairwise incomplete
#' \code{(NA, !NA)} observations with at least two of those incomplete pairs
#' having NA in \code{x}.
#' @param y A numeric vector. Must contain at least two NAs paired with
#' non-missing values in \code{x}. For more information on pairwise
#' requirements, see description of \code{x}. Ignored if \code{x} is a matrix.
#' @param var.equal Logical indicating whether to use Lin and Stivers' test
#' under heteroscedasticity or Ekbohm's test under homoscedasticity.
#' @param alternative The alternative hypothesis, one of \code{"two.sided"},
#' \code{"greater"}, or \code{"less"}
#' 
#' @return
#' A list with class \code{"htest"} used for internal purposes to produce
#' printout.
#' 
#' @examples
#' # growth in pine trees between 3 and 5 years of age
#' x <- Loblolly[Loblolly$age == 3, "height"]
#' y <- Loblolly[Loblolly$age == 5, "height"]
#' # assign NA to first five obs of age 3, x, and last three obs of age 5, y
#' x[1:5] <- NA
#' y[11:14] <- NA
#' 
#' # we expect a difference but test to make sure it exists, under the
#' # assumption of heteroscedasticity
#' mle.test(x, y)
#' # and under the assumption of homoscedasticity
#' mle.test(x, y, var.equal=TRUE)
#' 
#' # we can also pass x and y as a matrix
#' X <- cbind(x, y)
#' mle.test(X)
#' 
#' @export
mle.test <- function (x, y=NULL, var.equal=FALSE,
                      alternative="two.sided") {
    # type checking
    e1 <- paste("Either x must be an n by 2 matrix, or x and y must both be",
                "length n vectors.")
    if (is.matrix(x)) {
        if (length(dim(x)) != 2 || !(2 %in% dim(x))) {
            stop(paste("x has bad dimensionality.", e1))
        }
        
        if (dim(x)[1] == 2) {
            X <- t(x)
        } else {
            X <- x
        }
    } else if (is.vector(x) && is.vector(y)) {
        if (length(x) != length(y)) {
            stop(paste("x and y have incompatible lengths.", e1))
        }
        X <- cbind(x, y)
    } else {
        stop(paste("x or y has incorrect type.", e1))
    }
    
    if  (!(alternative %in% c("two.sided", "less", "greater"))) {
        stop(paste("Incorrectly specified alternative, must be one of",
                   "two.sided, less, or greater."))
    }
    
    x.na <- is.na(X[,1])
    y.na <- is.na(X[,2])
    n.obs <- list(
        paired = sum(!(x.na | y.na)),
        unpaired.x = sum(!x.na & y.na),
        unpaired.y = sum(x.na & !y.na)
    )
    # observation size check
    for (n in n.obs) {
        if (n < 2) {
            stop(paste("For each observation type: paired, unpaired", 
                       "type 1, and unpaired type 2; there must be at least",
                       "2 observations."))
        }
    }
    
    if (any(is.na(x) & is.na(y))){
        warning("Warning: Pairwise NA values are dropped.")
    }
    
    out <- list()
    # tests
    x1 <- X[!(is.na(X[,1]) | is.na(X[,2])),]
    n <-  colSums(is.na(X)) # == c(n3, n2)
    n1 <- dim(x1)[1]
    r <- cor(x1[,1], x1[,2])
    
    d <- prod(n1+n)-prod(n)*r^2
    m <- c(mean(X[y.na, 1]), mean(X[x.na, 2]))
    m1 <- apply(x1, 2, mean)
    v1 <- apply(x1, 2, var)
    
    if (var.equal == TRUE) {
        fg <- n1 * (n1 + cor(x1)%*%n) / d
        v <- c(var(X[y.na, 1]), var(X[x.na, 2]))
        
        V1 <- sum(v1*(n1-1)) + (1+r^2)*v%*%(n-1)[2:1]
        V1 <- V1 / (2*(n1-1) + (1+r^2)*(sum(n)-2))
        V1 <- V1*(1-r) * (sum(n)*(1+r) + 2*n1) / d
        stat <- sum(c(1, -1)*cbind((1-fg), fg)*cbind(m, m1)) / sqrt(V1)
        
        out$method <- "Ekbohm's MLE-Based Test Under Homoscedasticity"
    } else {
        C <- cov(x1)
        # note that fg[2]=-g, with g as presented in Kuan and Huang [2013]
        fg <- n1 * c(1, -1) * (n1 + (C/v1)%*%n) / d
        V1 <- t(fg)%*%C%*%fg/n1 + (v1/n[2:1])%*%(c(1, -1)-fg)^2
        stat <- (m%*%(c(1, -1)-fg) + m1%*%fg) / sqrt(V1)
        
        out$method <- "Lin and Stivers' MLE-Based Test Under Heteroscedasticity"
    }
    
    if (alternative == "less") {
        p <- pt(stat, n1)
    } else if (alternative == "greater") {
        p <- pt(stat, n1, lower.tail=FALSE)
    } else {
        p <- 2*pt(-abs(stat), n1)
    }
    
    # create test object
    out$null.value <- 0
    names(out$null.value) <- "difference in means"
    out$alternative <- alternative
    out$estimate <- m
    names(out$estimate) <- c("mean x", "mean y")
    if (is.null(y)) {
        out$data.name <- deparse(substitute(x)) 
    } else {
        out$data.name <- paste(deparse(substitute(x)),
                               "and", deparse(substitute(y)))
    }
    out$statistic <- stat
    names(out$statistic) <- "Z"
    out$parameter <- n1
    names(out$parameter) <- c("df")
    out$p.value <- p
    
    class(out) <- "htest"
    return(out)
}