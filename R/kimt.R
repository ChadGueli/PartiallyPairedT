#' Kim et al.’s Modified t-Based Test
#'
#' Method for testing difference of means in partially paired data with a
#' t-based test as proposed by Kim et al. [2004].
#' 
#' @name kim.mod.t
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
#' # we expect a difference but test to make sure it exists
#' kim.mod.t(x, y)
#' 
#' # we can also pass x and y as a matrix
#' X <- cbind(x, y)
#' kim.mod.t(X)
#' 
#' @export
kim.mod.t <- function (x, y=NULL, alternative="two.sided") {
    # type checking
    e1 <- paste("Either x must be an n by 2 matrix, or x and y must both be",
                "length n vectors.")
    if (is.matrix(x)) {
        if (length(dim(x)) != 2 || !(2 %in% dim(x))) {
            stop(paste("Error: x has bad dimensionality.", e1))
        }
        
        if (dim(x)[1] == 2) {
            x_ <- x[1,]
            y_ <- x[2,]
        } else {
            x_ <- x[,1]
            y_ <- x[,2]
        }
    } else if (is.vector(x) && is.vector(y)) {
        if (length(x) != length(y)) {
            stop(paste("Error: x and y have incompatible lengths.", e1))
        }
        x_ <- x
        y_ <- y
    } else {
        stop(paste("Error: x or y has incorrect type.", e1))
    }
    
    if  (!(alternative %in% c("two.sided", "less", "greater"))) {
        stop(paste("Error: incorrectly specified alternative, must be one of",
                   "two.sided, less, or greater."))
    }
    
    x.na <- is.na(x_)
    y.na <- is.na(y_)
    obs <- list(
        d = x_[!(x.na | y.na)] - y_[!(x.na | y.na)],
        x = x_[!x.na & y.na],
        y = -y_[x.na & !y.na]
    )
    if (any(is.na(x_) & is.na(y_))){
        warning("Warning: Pairwise NA values are dropped.")
    }
    
    # observation size check
    for (vec in obs) {
        if (length(vec) < 2) {
            stop(paste("Error: For each observation type: paired, unpaired", 
                       "type 1, and unpaired type 2; there must be at least",
                       "2 observations."))
        }
    }
    
    # tests
    n_ <- c(length(obs$x), length(obs$y))
    n <- c(length(obs$d), rep(1/mean(1/n_), 2)) # n to be used in numerator
    n_ <- c(n[1], n[2]^2/n_) # n to be used in denominator
    
    means <- sapply(obs, mean) # note that y is negated
    vars <- sapply(obs, var)
    stat <- means%*%n / sqrt(vars%*%n_)
    
    if (alternative == "less") {
        p <- pnorm(stat)
    } else if (alternative == "greater") {
        p <- pnorm(stat, lower.tail=FALSE)
    } else {
        p <- 2*pnorm(-abs(stat))
    }
    
    # create test object
    out <- list()
    out$null.value <- 0
    names(out$null.value) <- "difference in means"
    out$alternative <- alternative
    out$method <- "Kim et al.'s Modified t-Based Test"
    out$estimate <- c(means[1],
                      means[2]+means[3])
    names(out$estimate) <- c("mean difference", "difference of means")
    if (is.null(y)) {
        out$data.name <- deparse(substitute(x)) 
    } else {
        out$data.name <- paste(deparse(substitute(x)),
                               "and", deparse(substitute(y)))
    }
    out$statistic <- stat
    names(out$statistic) <- "t-Based Statistic"
    out$p.value <- p
    
    class(out) <- "htest"
    return(out)
}