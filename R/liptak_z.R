
#' Liptak's Weighted Z
#'
#' Method for testing difference of means in partially paired data based on
#' p-value pooling meta-analysis as proposed by Liptak [1958], with p-values
#' obtained by t-test.
#' 
#' @name liptak.Z
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
#' @param null.mu Number indicating the true difference of means.
#' @param alternative The alternative hypothesis, one of \code{"two.sided"},
#' \code{"greater"}, or \code{"less"}
#' @param var.equal Logical indicating whether the variance is equal for
#'  testing difference of means. Default \code{NULL} will run a test for
#'  equality of variances.
#'  @param alpha Level at which to reject the assumption of equal variance.
#'  Ignored if \code{var.equal} is not \code{NULL}.
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
#' liptak.Z(x, y)
#' 
#' # we can also pass x and y as a matrix
#' X <- cbind(x, y)
#' liptak.Z(X)
#' 
#' @export
liptak.Z <- function (x, y=NULL, null.mu=0, alternative="two.sided",
                      var.equal=NULL, var.alpha=0.2) {
    # type checking
    e1 <- paste("Either x must be an n by 2 matrix, or x and y must both be",
                "length n vectors.")
    if (is.matrix(x)) {
        if (length(dim(x)) != 2 || !(2 %in% dim(x))) {
            stop(paste("x has bad dimensionality.", e1))
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
            stop(paste("x and y have incompatible lengths.", e1))
        }
        x_ <- x
        y_ <- y
    } else {
        stop(paste("x or y has incorrect type.", e1))
    }
    
    if  (!(alternative %in% c("two.sided", "less", "greater"))) {
        stop(paste("Incorrectly specified alternative, must be one of",
                   "two.sided, less, or greater."))
    }
    
    x.na <- is.na(x_)
    y.na <- is.na(y_)
    obs <- list(
        paired.x = x_[!(x.na | y.na)],
        paired.y = y_[!(x.na | y.na)],
        unpaired.x = x_[!x.na & y.na],
        unpaired.y = y_[x.na & !y.na]
    )
    if (any(is.na(x_) & is.na(y_))){
        warning("Warning: Pairwise NA values are dropped.")
    }
    
    # observation size check
    for (vec in obs) {
        if (length(vec) < 2) {
            stop(paste("For each observation type: paired, unpaired", 
                       "type 1, and unpaired type 2; there must be at least",
                       "2 observations."))
        }
    }
    
    # tests
    if (is.null(var.equal)){
        if (var.test(obs$unpaired.x, obs$unpaired.y)$p.value > var.alpha) {
            var.equal <- TRUE
        } else {
            var.equal <- FALSE
        }
    }
    
    if (alternative == "two.sided") {
        paired.t <- t.test(obs$paired.x, obs$paired.y, alternative="greater",
                           paired=TRUE)
        unpaired.t <- t.test(obs$unpaired.x, obs$unpaired.y,
                             alternative="greater",
                             var.equal=var.equal)
    } else {
        paired.t <- t.test(obs$paired.x, obs$paired.y, alternative=alternative,
                           paired=TRUE)
        unpaired.t <- t.test(obs$unpaired.x, obs$unpaired.y,
                             alternative=alternative,
                             var.equal=var.equal)
    }
    p <- c(paired.t$p.value, unpaired.t$p.value)
    
    w <- c(2*length(obs$paired.x),
           length(obs$unpaired.x)+length(obs$unpaired.y))
    w <- sqrt(w/sum(w))
    z <- qnorm(1-p)
    stat <- w%*%z
    p.c <- 1-pnorm(stat)
    if (alternative == "two.sided") {
        p.c <- if (p.c < 1/2) 2*p.c else 2*(1-p.c)
    }
    
    # create test object
    out <- list()
    out$null.value <- null.mu
    names(out$null.value) <- "difference in means"
    out$alternative <- alternative
    out$method <- paste("Liptak's Weighted Z-test with",
                        if (var.equal) "Equal" else "Enequal",
                        "Variance")
    out$estimate <- c(paired.t$estimate,
                      unpaired.t$estimate[1]-unpaired.t$estimate[2])
    names(out$estimate) <- c("mean difference", "difference of means")
    if (is.null(y)) {
        out$data.name <- deparse(substitute(x)) 
    } else {
        out$data.name <- paste(deparse(substitute(x)),
                               "and", deparse(substitute(y)))
    }
    out$statistic <- stat
    names(out$statistic) <- "weighted mean of Zs"
    out$parameter <- c(paired.t$parameter, unpaired.t$parameter)
    names(out$parameter) <- c("paired df", "unpaired df")
    out$p.value <- p.c
    
    class(out) <- "htest"
    return(out)
}