#' Generalized Partial Credit Model
#' @description Common computations and operatoins for the GPCM
#' @name model_gpcm
NULL

#' @rdname model_gpcm
#' @param t ability parameters, 1d vector
#' @param a discrimination parameters, 1d vector
#' @param b item location parameters, 1d vector
#' @param d item category parameters, 2d vector
#' @param D the scaling constant, default=1.702
#' @param d0 insert an initial category value
#' @return \code{model_gpcm_prob} returns the resulting probabilities in a 3d array
#' @details 
#' Use \code{NA} to represent unused category.
#' @examples 
#' with(model_gpcm_gendata(10, 5, 3), model_gpcm_prob(t, a, b, d))
#' @export 
model_gpcm_prob <- function(t, a, b, d, D=1.702, d0=NULL){
  d <- as.matrix(d)
  if(!is.null(d0))
    d <- cbind(d0, d)
  n_p <- length(t)
  n_i <- length(b)
  n_c <- ncol(d)
  p <- model_gpcm_probC(t, a, b, d, D)
  aperm(array(unlist(p), dim=c(n_p, n_c, n_i)), c(1, 3, 2))
}

#' @rdname model_gpcm
#' @return \code{model_gpcm_info} returns the resulting information in a 3d array
#' @examples 
#' with(model_gpcm_gendata(10, 5, 3), model_gpcm_info(t, a, b, d))
#' @export 
model_gpcm_info <- function(t, a, b, d, D=1.702, d0=NULL){
  p <- model_gpcm_prob(t, a, b, d, D, d0)
  n_i <- dim(p)[2]
  n_c <- dim(p)[3]
  if(length(a) == 1)
    a <- rep(a, n_i)
  rs <- array(NA, dim=dim(p))
  for(j in 1:n_i) {
    p_j <- if(length(t) == 1) t(p[,j,]) else p[,j,]
    rs[,j,] <- (D * a[j])^2 * p_j * colSums(1:n_c * t(p_j) * outer(1:n_c, colSums(t(p_j) * 1:n_c, na.rm=T), '-'), na.rm=T)
  }
  rs
}

#' @rdname model_gpcm
#' @param u observed scores (starting from 0), 2d matrix
#' @param log TRUE to return log-likelihood
#' @return \code{model_gpcm_lh} returns the resulting likelihood in a matrix
#' @examples 
#' with(model_gpcm_gendata(10, 5, 3), model_gpcm_lh(u, t, a, b, d))
#' @export
model_gpcm_lh <- function(u, t, a, b, d, D=1.702, d0=NULL, log=FALSE){
  p <- model_gpcm_prob(t, a, b, d, D, d0)
  ix <- model_polytomous_3dindex(u)
  lh <- array(p[ix], dim=dim(u))
  if(log) lh <- log(lh)
  lh
}

#' @rdname model_gpcm
#' @param n_p the number of people to be generated
#' @param n_i the number of items to be generated
#' @param n_c the number of score categories
#' @param sort_d \code{TRUE} to sort d parameters for each item
#' @param t_dist parameters of the normal distribution used to generate t-parameters
#' @param a_dist parameters of the lognormal distribution parameters of a-parameters
#' @param b_dist parameters of the normal distribution used to generate b-parameters
#' @param d_dist parameters of the normal distribution used to generate d-parameters
#' @param t_bounds the bounds of the ability parameters
#' @param a_bounds the bounds of the discrimination parameters
#' @param b_bounds the bounds of the difficulty parameters
#' @param d_bounds the bounds of the category parameters
#' @param missing the proportion or number of missing responses
#' @param ... additional arguments
#' @return \code{model_gpcm_gendata} returns the generated response matrix and parameters
#' @examples
#' model_gpcm_gendata(10, 5, 3)
#' model_gpcm_gendata(10, 5, 3, missing=.1)
#' @importFrom stats rnorm rlnorm rmultinom
#' @export
model_gpcm_gendata <- function(n_p, n_i, n_c, t=NULL, a=NULL, b=NULL, d=NULL, D=1.702, sort_d=FALSE, 
                               t_dist=c(0, 1), a_dist=c(-.1, .2), b_dist=c(0, .8), d_dist=c(0, 1),
                               t_bounds=c(-3, 3), a_bounds=c(.01, 2.5), b_bounds=c(-3, 3), d_bounds=c(-3, 3),
                               missing=NULL, ...){
  if(is.null(t)){
    t <- rnorm(n_p, mean=t_dist[1], sd=t_dist[2])
    t[t < t_bounds[1]] <- t_bounds[1]
    t[t > t_bounds[2]] <- t_bounds[2]
  }
  if(is.null(a)){
    a <- rlnorm(n_i, meanlog=a_dist[1], sdlog=a_dist[2])
    a[a < a_bounds[1]] <- a_bounds[1]
    a[a > a_bounds[2]] <- a_bounds[2]
  }
  if(is.null(b)){
    b <- rnorm(n_i, mean=b_dist[1], sd=b_dist[2])
    b[b < b_bounds[1]] <- b_bounds[1]
    b[b > b_bounds[2]] <- b_bounds[2]
  }
  if(is.null(d)) {
    d <- matrix(rnorm(n_i * n_c, mean=d_dist[1], sd=d_dist[2]), nrow=n_i, ncol=n_c)
    d[d < d_bounds[1]] <- d_bounds[1]
    d[d > d_bounds[2]] <- d_bounds[2]
    d[,  1] <- 0
    d[, -1] <- d[, -1] - rowMeans(d[, -1, drop=FALSE])
    if(sort_d)
      d[, -1] <- t(apply(d[, -1, drop=FALSE], 1, sort))
  }
  
  if(length(t) == 1)
    t <- rep(t, n_p)
  if(length(a) == 1)
    a <- rep(a, n_i)
  if(length(t) != n_p)
    stop('wrong dimensions for t')
  if(length(a) != n_i)
    stop('wrong dimensions for a')
  if(length(b) != n_i)
    stop('wrong dimensions for b')
  if(any(dim(d) != c(n_i, n_c)))
    stop('wrong dimensions for d')
  
  p <- model_gpcm_prob(t, a, b, d, D, NULL)
  u <- apply(p, 1:2, function(x) which.max(rmultinom(1, 1, x)[,1]) - 1)
  #u <-apply(p, 2, function(x) rowSums(runif(n_p) >= t(apply(x, 1, cumsum))))
  if(!is.null(missing)){
    missing <- floor(ifelse(missing < 1, missing * n_p * n_i, missing))
    idx <- sample(length(u), missing)
    u[cbind(ceiling(idx/n_i), (idx-1)%%n_i+1)] <- NA
  }
  
  list(u=u, t=t, a=a, b=b, d=d)
}


#' @rdname model_gpcm
#' @param scale the scale, 't' for theta or 'b' for b-parameters
#' @param mean the mean of the new scale
#' @param sd the standard deviation of the new scale
#' @return \code{model_gpcm_rescale} returns t, a, b, d parameters on the new scale
#' @importFrom stats sd
#' @export
model_gpcm_rescale <- function(t, a, b, d, scale=c("t", "b"), mean=0, sd=1){
  scale <- switch(match.arg(scale), "t"=t, "b"=b)
  slope <- sd / sd(scale)
  intercept <- mean - slope * mean(scale)
  t <- slope * t + intercept
  b <- slope * b + intercept
  a <- a / slope
  d <- d * slope
  list(t=t, a=a, b=b, d=d)
}

#' @rdname model_gpcm
#' @param type the type of plot, prob for ICC and info for IIFC
#' @param total TRUE to sum values over items
#' @param item_level TRUE to add item level data
#' @param xaxis the values of x-axis
#' @return \code{model_gpcm_plot} returns a \code{ggplot} object
#' @examples
#' # Figure 1 in Muraki, 1992 (APM)
#' b <- matrix(c(-2,0,2,-.5,0,2,-.5,0,2), nrow=3, byrow=TRUE)
#' model_gpcm_plot(a=c(1,1,.7), b=rowMeans(b), d=rowMeans(b)-b, D=1.0, d0=0)
#' # Figure 2 in Muraki, 1992 (APM)
#' b <- matrix(c(.5,0,NA,0,0,0), nrow=2, byrow=TRUE)
#' model_gpcm_plot(a=.7, b=rowMeans(b, na.rm=TRUE), d=rowMeans(b, na.rm=TRUE)-b, D=1.0, d0=0)
#' # Figure 3 in Muraki, 1992 (APM)
#' b <- matrix(c(1.759,-1.643,3.970,-2.764), nrow=2, byrow=TRUE)
#' model_gpcm_plot(a=c(.778,.946), b=rowMeans(b), d=rowMeans(b)-b, D=1.0, d0=0)
#' # Figure 1 in Muraki, 1993 (APM)
#' b <- matrix(c(0,-2,4,0,-2,2,0,-2,0,0,-2,-2,0,-2,-4), nrow=5, byrow=TRUE)
#' model_gpcm_plot(a=1, b=rowMeans(b), d=rowMeans(b)-b, D=1.0)
#' # Figure 2 in Muraki, 1993 (APM)
#' b <- matrix(c(0,-2,4,0,-2,2,0,-2,0,0,-2,-2,0,-2,-4), nrow=5, byrow=TRUE)
#' model_gpcm_plot(a=1, b=rowMeans(b), d=rowMeans(b)-b, D=1.0, type='info', item_level=TRUE)
#' @import ggplot2
#' @importFrom stats aggregate
#' @export
model_gpcm_plot <- function(a, b, d, D=1.702, d0=NULL, type=c('prob', 'info'), item_level=FALSE, total=FALSE, xaxis=seq(-6, 6, .1)){
  rs <- switch(match.arg(type), "prob"=model_gpcm_prob, "info"=model_gpcm_info)(xaxis, a, b, d, D, d0)
  n_p <- dim(rs)[1]
  n_i <- dim(rs)[2]
  n_c <- dim(rs)[3]
  y <- NULL
  for(i in 1:n_i)
    y <- rbind(y, data.frame(theta=rep(xaxis, n_c), item=paste('Item', i), category=paste('Category', rep(1:n_c, each=n_p)), x=as.vector(rs[, i, ])))
  if(item_level)
    y <- rbind(y, cbind(aggregate(y$x, by=list(theta=y$theta, item=y$item), sum), category='Total'))
  if(total)
    y <- cbind(aggregate(y$x, by=list(theta=y$theta, category=y$category), sum), item='Total')
  
  y <- y[!is.na(y$x),]
  ggplot(y, aes_string(x="theta", y="x", color="category")) +
    geom_line() + facet_wrap(~item, scales='free') +
    xlab(expression(theta)) + ylab(type) +
    guides(color=FALSE) + theme_bw() + theme(legend.key=element_blank())
}

#' @rdname model_gpcm
#' @param verbose TRUE to print rough maximum likelihood values
#' @return \code{model_gpcm_plot_loglh} returns a \code{ggplot} object
#' @examples
#' with(model_gpcm_gendata(5, 50, 3), model_gpcm_plot_loglh(u, a, b, d))
#' @import ggplot2
#' @export
model_gpcm_plot_loglh <- function(u, a, b, d, D=1.702, d0=NULL, xaxis=seq(-6, 6, .1), verbose=FALSE){
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_t <- length(xaxis)
  rs <- array(NA, dim=c(n_p, n_t))
  for(i in 1:n_t)
    rs[, i] <- rowSums(model_gpcm_lh(u, rep(xaxis[i], n_p), a, b, d, D, d0, log=TRUE))
  if(verbose) 
    print(apply(rs, 1, function(x){xaxis[which.max(x)]}))
 
  rs <- data.frame(theta=rep(xaxis, each=n_p), people=rep(1:n_p, n_t), value=as.vector(rs))
  rs$people <- factor(rs$people)
  ggplot(rs, aes_string(x="theta", y="value", color="people")) +
    geom_line() + xlab(expression(theta)) + ylab("Log-likelihood") +
    guides(color=FALSE) + theme_bw()
}

