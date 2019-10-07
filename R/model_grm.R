#' Graded Response Model
#' @description Common computations and operations for the GRM
#' @name model_grm
NULL

#' @rdname model_grm
#' @param t ability parameters, 1d vector
#' @param a discrimination parameters, 1d vector
#' @param b item location parameters, 2d matrix
#' @param D the scaling constant, default=1.702
#' @param raw TRUE to return P*
#' @return \code{model_grm_prob} returns the resulting probabilities in a 3d array
#' @examples 
#' with(model_grm_gendata(10, 5, 3), model_grm_prob(t, a, b))
#' @export
model_grm_prob <- function(t, a, b, D=1.702, raw=FALSE){
  b <- as.matrix(b)
  n_p <- length(t)
  n_i <- nrow(b)
  n_c <- ncol(b)
  if(raw) {
    p <- model_grm_prob_rawC(t, a, b, D)
    p <- aperm(array(unlist(p), dim=c(n_p, n_c+2, n_i)), c(1, 3, 2))
  } else {
    p <- model_grm_probC(t, a, b, D)
    p <- aperm(array(unlist(p), dim=c(n_p, n_c+1, n_i)), c(1, 3, 2))
  }
  p
}

#' @rdname model_grm
#' @return \code{model_grm_info} returns the resulting information in a 3d array
#' @examples 
#' with(model_grm_gendata(10, 5, 3), model_grm_info(t, a, b))
#' @export 
model_grm_info <- function(t, a, b, D=1.702){
  p  <- p_ <- Rirt::model_grm_prob(t, a, b, D)
  p_[is.na(p_)] <- 0
  p_ <- aperm(apply(p_, c(1, 2), function(x) rev(cumsum(c(0, rev(x))))), c(2, 3, 1))
  n_c <- dim(p)[3]
  dv1_p_ <- aperm(p_ * (1 - p_), c(2, 3, 1)) * D * a
  dv2_p_ <- aperm((1 - 2 * p_) * p_ * (1 - p_), c(2, 3, 1)) * (D * a)^2
  dv1_p <- dv1_p_[,1:n_c,,drop=FALSE] - dv1_p_[,-1,,drop=FALSE]
  dv1_p <- aperm(dv1_p, c(3, 1, 2))
  dv2_p <- dv2_p_[,1:n_c,,drop=FALSE] - dv2_p_[,-1,,drop=FALSE]
  dv2_p <- aperm(dv2_p, c(3, 1, 2))
  1 / p * dv1_p^2 - dv2_p
}

#' @rdname model_grm
#' @param u observed scores (starting from 0), 2d matrix
#' @param log TRUE to return log-likelihood
#' @return \code{model_grm_lh} returns the resulting likelihood in a matrix
#' @examples 
#' with(model_grm_gendata(10, 5, 3), model_grm_lh(u, t, a, b))
#' @export
model_grm_lh <- function(u, t, a, b, D=1.702, log=FALSE){
  p <- model_grm_prob(t, a, b, D)
  ix <- model_polytomous_3dindex(u)
  lh <- array(p[ix], dim=dim(u))
  if(log) lh <- log(lh)
  lh
}

#' @rdname model_grm
#' @param n_p the number of people to be generated
#' @param n_i the number of items to be generated
#' @param n_c the number of score categories
#' @param t_dist parameters of the normal distribution used to generate t-parameters
#' @param a_dist parameters of the lognormal distribution used to generate a-parameters
#' @param b_dist parameters of the normal distribution used to generate b-parameters
#' @param t_bounds the bounds of the ability parameters
#' @param a_bounds the bounds of the discrimination parameters
#' @param b_bounds the bounds of the difficulty parameters
#' @param missing the proportion or number of missing responses
#' @return \code{model_grm_gendata} returns the generated response data and parameters in a list
#' @examples
#' model_grm_gendata(10, 5, 3)
#' model_grm_gendata(10, 5, 3, missing=.1)
#' @importFrom stats rnorm rlnorm rmultinom
#' @export
model_grm_gendata <- function(n_p, n_i, n_c, t=NULL, a=NULL, b=NULL, D=1.702, 
                              t_dist=c(0, 1), a_dist=c(-.1, .2), b_dist=c(0, .8), 
                              t_bounds=c(-3, 3), a_bounds=c(.01, 2.5), b_bounds=c(-3, 3),
                              missing=NULL){
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
  if(is.null(b)) {
    b <- matrix(rnorm(n_i * (n_c - 1), mean=b_dist[1], sd=b_dist[2]), nrow=n_i)
    b <- t(apply(b, 1, sort))
    b <- matrix(b, nrow=n_i, ncol=n_c-1)
  }
  
  if(length(t) == 1)
    t <- rep(t, n_p)
  if(length(a) == 1)
    a <- rep(a, n_i)
  if(length(t) != n_p)
    stop('wrong dimensions for t')
  if(length(a) != n_i)
    stop('wrong dimensions for a')
  if(any(dim(b) != c(n_i, n_c - 1)))
    stop('wrong dimensions for b')
  
  p <- model_grm_prob(t, a, b, D)
  u <- apply(p, 1:2, function(x) which.max(rmultinom(1, 1, x)[,1]) - 1)
  #u <- apply(p, 2, function(x) rowSums(runif(n_p) >= t(apply(x, 1, cumsum))))
  if(!is.null(missing)){
    missing <- floor(ifelse(missing < 1, missing * n_p * n_i, missing))
    idx <- sample(length(u), missing)
    u[cbind(ceiling(idx/n_i), (idx-1)%%n_i+1)] <- NA
  }
  
  list(u=u, t=t, a=a, b=b)
}

#' @rdname model_grm
#' @param scale the scale, 't' for theta or 'b' for b-parameters
#' @param mean the mean of the new scale
#' @param sd the standard deviation of the new scale
#' @return \code{model_grm_rescale} returns t, a, b parameters on the new scale
#' @importFrom stats sd
#' @export
model_grm_rescale <- function(t, a, b, scale=c("t", "b"), mean=0, sd=1){
  scale <- switch(match.arg(scale), "t"=t, "b"=b)
  slope <- sd / sd(scale)
  intercept <- mean - slope * mean(scale)
  t <- slope * t + intercept
  b <- slope * b + intercept
  a <- a / slope
  list(t=t, a=a, b=b)
}

#' @rdname model_grm
#' @param type the type of plot, prob for ICC and info for IIFC
#' @param total TRUE to sum values over items
#' @param item_level TRUE to combine categories
#' @param xaxis the values of x-axis
#' @return \code{model_grm_plot} returns a \code{ggplot} object
#' @examples
#' with(model_grm_gendata(10, 5, 3), model_grm_plot(a, b, type='prob'))
#' with(model_grm_gendata(10, 5, 3), model_grm_plot(a, b, type='info', item_level=TRUE))
#' @import ggplot2
#' @importFrom stats aggregate
#' @export
model_grm_plot <- function(a, b, D=1.702, type=c('prob', 'info'), item_level=FALSE, total=FALSE, xaxis=seq(-6, 6, .1), raw=FALSE){
  rs <- switch(match.arg(type), "prob"=model_grm_prob(xaxis, a, b, D, raw), "info"=model_grm_info(xaxis, a, b, D))
  n_p <- dim(rs)[1]
  n_i <- dim(rs)[2]
  n_c <- dim(rs)[3]
  y <- NULL
  for(i in 1:n_i)
    y <- rbind(y, data.frame(theta=rep(xaxis, n_c), item=paste('Item', i), category=paste('Category', rep(1:n_c, each=n_p)), x=as.vector(rs[,i,])))
  if(item_level) y <- rbind(y, cbind(aggregate(y$x, by=list(theta=y$theta, item=y$item), sum), category='Total'))
  if(total) y <- cbind(aggregate(y$x, by=list(theta=y$theta, category=y$category), sum), item='Total')
  
  y <- y[!is.na(y$x),]
  ggplot(y, aes_string(x="theta", y="x", color="category")) +
    geom_line() + facet_wrap(~item, scales='free') +
    xlab(expression(theta)) + ylab(type) +
    guides(color=FALSE) + theme_bw() + theme(legend.key=element_blank())
}

#' @rdname model_grm
#' @param verbose TRUE to print rough maximum likelihood values
#' @return \code{model_grm_plot_loglh} returns a \code{ggplot} object
#' @examples
#' with(model_grm_gendata(5, 50, 3), model_grm_plot_loglh(u, a, b))
#' @import ggplot2
#' @export
model_grm_plot_loglh <- function(u, a, b, D=1.702, xaxis=seq(-6, 6, .1), verbose=FALSE){
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_t <- length(xaxis)
  rs <- array(NA, dim=c(n_p, n_t))
  for(i in 1:n_t)
    rs[, i] <- rowSums(model_grm_lh(u, rep(xaxis[i], n_p), a, b, D, log=TRUE))
  if(verbose)
    print(apply(rs, 1, function(x){xaxis[which.max(x)]}))
  
  rs <- data.frame(theta=rep(xaxis, each=n_p), people=rep(1:n_p, n_t), value=as.vector(rs))
  rs$people <- factor(rs$people)
  ggplot(rs, aes_string(x="theta", y="value", color="people")) +
    geom_line() + xlab(expression(theta)) + ylab("Log-likelihood") +
    guides(color=FALSE) + theme_bw()
}
