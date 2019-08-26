#' Utility functions
#' @name utils
NULL


#' @rdname utils
#' @description \code{rmse} computes the root mean squared error (RMSE) 
#' of two numeric vectors/matrices
#' @param x1 a numeric vector or matrix
#' @param x2 a numeric vector or matrix
#' @examples 
#' rmse(rnorm(10), rnorm(10))
#' @export
rmse <- function(x1, x2){
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  if(any(dim(x1) != dim(x2)))
    stop("inputs are of different dimensions")
  sqrt(colMeans((x1 - x2)^2))
}

#' @rdname utils
#' @description \code{freq} computes the frequency counts of 
#' a numeric or character vector
#' @param x a numeric or character vector
#' @param vals valid values, \code{NULL} to include all values
#' @param decimal round results to n-th decimal places
#' @return \code{freq} returns the frequency counts and percentages in a data.frame
#' @examples 
#' freq(round(runif(100, 1, 5)))
#' @export
freq <- function(x, vals=NULL, decimal=NULL){
  if(is.null(vals)) 
    vals <- sort(unique(x))
  rs <- data.frame(table(factor(x, levels=vals, labels=vals)), stringsAsFactors=FALSE)
  colnames(rs) <- c("value", "freq")
  
  rs$perc <- rs$freq / sum(rs$freq)
  rs$cum_freq <- cumsum(rs$freq)
  rs$cum_perc <- cumsum(rs$perc)
  
  if(!is.null(decimal)){
    rs$perc <- round(rs$perc, decimal)
    rs$cum_perc <- round(rs$cum_perc, decimal)
  }
  rs
}

#' @rdname utils
#' @description \code{cronbach_alpha} computes the Cronbach's alpha 
#' internal consistency reliability index
#' @param u oberved responses, 2d matrix
#' @examples
#' cronbach_alpha(model_3pl_gendata(1000, 20)$u)
#' @importFrom stats var
#' @export
cronbach_alpha <- function(u){
  k <- ncol(u)
  total_var <- var(rowSums(u, na.rm=TRUE))
  item_var <- sum(apply(u, 2, var, na.rm=TRUE))
  k / (k - 1) * (1 - item_var / total_var)
}

#' @rdname utils 
#' @description \code{spearman_brown} predicts the reliability when the 
#' current test is extended to n times longer
#' @param rho the reliability of the current test
#' @param n_len extend the test to n times longer
#' @examples
#' spearman_brown(.70, 2)
#' @export
spearman_brown <- function(rho, n_len){
  n_len * rho / (1 + (n_len - 1) * rho)
}

#' @rdname utils
#' @description \code{spearman_brown_reverse} computes how many times
#' the current test length needs to be extended in order to reach targeted
#' reliability
#' @param target_rho the targeted reliability
#' @examples 
#' spearman_brown_reverse(.70, .85)
#' @export
spearman_brown_reverse <- function(rho, target_rho){
  target_rho * (1 - rho) / rho / (1 - target_rho)
}

#' @rdname utils
#' @description \code{quadratic kappa} computes the quadratic weighted kappa
#' of two numeric vectors
#' @examples 
#' quadratic_kappa(round(runif(100, 1, 5)), round(runif(100, 1, 5)))
#' @importFrom reshape2 acast
#' @export 
quadratic_kappa <- function(x1, x2) { 
  if(length(x1) != length(x2))
    stop('inputs are of different size')
  n_data <- length(x1)
  vals <- sort(unique(c(x1, x2)))
  n_vals <- length(vals)
  
  w <- abs(outer(vals, vals, '-'))^2 / (n_vals - 1)^2
  obs <- expand.grid(x1=vals, x2=vals)
  obs$n <- 0
  for(i in 1:nrow(obs))
    obs$n[i] <- sum(x1==obs$x1[i] & x2==obs$x2[i], na.rm=TRUE)
  obs <- acast(obs, x1 ~ x2, value.var='n')
  obs <- obs / sum(obs, na.rm=TRUE)
  obs <- obs[match(rownames(obs), vals), match(colnames(obs), vals)]
  exp <- outer(rowSums(obs, na.rm=TRUE), colSums(obs, na.rm=TRUE), "*")
  1 - sum(obs*w, na.rm=TRUE) / sum(exp*w, na.rm=TRUE)
}



