#' Looney and Jones’s Corrected Z-Test
#'
#' Method for testing difference of means in partially paired data using
#' correlation between paired observations to correct variance as proposed by
#' Looney and Jones [2003].
#' 
#' @name LJ.Zcorr
#'
#' @param x A numeric vector or 2xn matrix. If \code{x} is a vector, then
#' data should be paired with \code{y} by index. If \code{x} is a matrix,
#' then pairing is done within column.
#' @param y A numeric vector. Ignored if \code{x} is a matrix.
#' @param alternative The alternative hypothesis, one of \code{"two.sided"},
#' \code{"greater"}, or \code{"less"}
#'  
#' @return
#' A list with class \code{"htest"} used for internal purposes to produce
#' printout.
#' 
#' @examples
#' # growth in orange trees between 118 and 664 days of age
#' x <- Orange[Orange$age == 118, "circumference"]
#' y <- Orange[Orange$age == 664, "circumference"]
#' # assign NA to first tree at age 118, x, and last tree at age 664, y
#' x[1] <- NA
#' y[5] <- NA
#' 
#' # we expect a difference but test to make sure it exists
#' LJ.Zcorr(x, y)
#' 
#' # we can also pass x and y as a matrix
#' X <- cbind(x, y)
#' LJ.Zcorr(X)
#' 
#' # Finally, we don't require every group be represented. For example,
#' # if there are no unpaired x obs:
#' x <- Orange[Orange$age == 118, "circumference"]
#' # then the test is still valid:
#' LJ.Zcorr(x, y)
#' 
#' @export
LJ.Zcorr <- function (x, y=NULL, alternative="two.sided") {
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
    
    if (any(is.na(x_) & is.na(y_))){
        warning("Warning: Pairwise NA values are dropped.")
    }
    
    # tests
    n1 <- sum(!(is.na(x_) | is.na(y_)))
    n <- c(sum(!is.na(x_)), sum(!is.na(y_)))
    for (i in n) {
        if (i < 2) {
            stop(paste("There must be at least two of both type 1",
                        "and type 2 observation."))
        }
    }
    C <- cov(cbind(x_, y_), use="pairwise.complete.obs")
    # adjust C[1, 2] and C[2, 1] while eliminating NA
    C[2:3] <- ifelse(is.na(C[2:3]), 0, -n1*C[2:3]/n) 
    stat <- (mean(x_, na.rm=TRUE)-mean(y_, na.rm=TRUE))/sqrt(sum(C %*% (1/n)))
    
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
    out$method <- "Looney and Jones' Corrected Z"
    out$estimate <- c(mean(x_, na.rm=TRUE), mean(y_, na.rm=TRUE))
    names(out$estimate) <- c("mean x", "mean y")
    if (is.null(y)) {
        out$data.name <- deparse(substitute(x)) 
    } else {
        out$data.name <- paste(deparse(substitute(x)),
                               "and", deparse(substitute(y)))
    }
    out$statistic <- stat
    names(out$statistic) <- "Zcorr"
    out$p.value <- p
    
    class(out) <- "htest"
    return(out)
}