#' Estimation of the 3PL model
#' @description Estimate the 3PL model using the joint or marginal 
#' maximum likelihood estimation methods
#' @name estimate_3pl 
NULL

#' @rdname estimate_3pl
#' @description \code{model_3pl_eap} scores response vectors using the EAP method
#' @return \code{model_3pl_eap} returns theta estimates and standard errors in a list
#' @examples 
#' with(model_3pl_gendata(10, 40), 
#'     cbind(true=t, est=model_3pl_eap(u, a, b, c)$t))
#' @importFrom stats dnorm 
#' @export
model_3pl_eap <- function(u, a, b, c, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4)){
  quad <- hermite_gauss('11')
  quad$w <- quad$w * exp(quad$t^2) * dnorm(quad$t, priors[1], priors[2])
  p <- model_3pl_prob(quad$t, a, b, c, D)
  lh <- exp(ifelse(is.na(u), 0, u) %*% t(log(p)) + ifelse(is.na(1 - u), 0, 1 - u) %*% t(log(1 - p)))
  t <- colSums(t(lh) * quad$w * quad$t) / colSums(t(lh) * quad$w)
  t[t < bounds_t[1]] <- bounds_t[1]
  t[t > bounds_t[2]] <- bounds_t[2]
  t_sd <- colSums(t(lh) * quad$w * outer(quad$t, t, '-')^2) / colSums(t(lh) * quad$w)
  t_sd <- sqrt(t_sd)
  list(t=t, sd=t_sd)
}

#' @rdname estimate_3pl
#' @description \code{model_3pl_map} scores response vectors using the MAP method
#' @return \code{model_3pl_map} returns theta estimates in a list
#' @examples 
#' with(model_3pl_gendata(10, 40), 
#'      cbind(true=t, est=model_3pl_map(u, a, b, c)$t))
#' @export
model_3pl_map <- function(u, a, b, c, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4), iter=30, conv=1e-3){
  t <- rnorm(dim(u)[1], 0, .01)
  t_free <- rep(T, length(t))
  for(m in 1:iter){
    dv_t <- model_3pl_dv_jmle(model_3pl_dv_Pt, u, t, a, b, c, D)
    dv_t$dv1 <- colSums(dv_t$dv1, na.rm=TRUE)
    dv_t$dv2 <- colSums(dv_t$dv2, na.rm=TRUE)
    if(!is.null(priors)){
      dv_t$dv1 <- dv_t$dv1 - (t - priors[1]) / priors[2]^2
      dv_t$dv2 <- dv_t$dv2 - 1 / priors[2]^2
    }
    nr_t <- nr_iteration(t, t_free, dv_t, 1.0, 1.0, bounds_t)
    t <- nr_t$param
    if(max(abs(nr_t$h[t_free])) < conv) break
  }
  list(t=t)
}

#' @rdname estimate_3pl
#' @keywords internal
model_3pl_dv_Pt <- function(t, a, b, c, D){
  p <- t(model_3pl_prob(t, a, b, c, D))
  dv1 <- D * a * (p - c) * (1 - p) / (1 - c)
  dv2 <- (D * a / (1 - c))^2 * (p - c) * (1 - p) * (1 + c - 2*p)
  list(dv1=dv1, dv2=dv2, p=p)
}

#' @rdname estimate_3pl
#' @keywords internal
model_3pl_dv_Pa <- function(t, a, b, c, D){
  p <- t(model_3pl_prob(t, a, b, c, D))
  dv1 <- D * t(outer(t, b, '-')) * (p - c) * (1 - p) / (1 - c)
  dv2 <- (D * t(outer(t, b, '-')) / (1 - c))^2 * (p - c) * (1 - p) * (1 + c - 2*p)
  list(dv1=dv1, dv2=dv2, p=p)
}

#' @rdname estimate_3pl
#' @keywords internal
model_3pl_dv_Pb <- function(t, a, b, c, D){
  p <- t(model_3pl_prob(t, a, b, c, D))
  dv1 <- - D * a * (p - c) * (1 - p) / (1 - c)
  dv2 <- (D * a / (1 - c))^2 * (p - c) * (1 - p) * (1 + c - 2*p)
  list(dv1=dv1, dv2=dv2, p=p)
}

#' @rdname estimate_3pl
#' @keywords internal
model_3pl_dv_Pc <- function(t, a, b, c, D){
  p <- t(model_3pl_prob(t, a, b, c, D))
  dv1 <- (1 - p) / (1 - c)
  dv2 <- array(0, dim=dim(dv1))
  list(dv1=dv1, dv2=dv2, p=p)
}

#' @rdname estimate_3pl
#' @keywords internal
model_3pl_dv_jmle <- function(pdv_fn, u, t, a, b, c, D){
  pdv <- pdv_fn(t, a, b, c, D)
  dv0 <- (t(u) - pdv$p) / pdv$p / (1 - pdv$p)
  dv1 <- dv0 * pdv$dv1
  dv2 <- dv0 * pdv$dv2 - dv1^2
  list(dv1=dv1, dv2=dv2)
}


#' @rdname estimate_3pl
#' @description \code{model_3pl_jmle} estimates the item and ability parameters
#' using the joint maximum likelihood estimation (JMLE) method
#' @param u observed response matrix, 2d matrix
#' @param t ability parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param a discrimination parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param b difficulty parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param c pseudo-guessing parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param D the scaling constant, 1.702 by default
#' @param iter the maximum iterations, default=100
#' @param conv the convergence criterion
#' @param nr_iter the maximum newton-raphson iterations, default=10
#' @param scale the mean and SD of the theta scale, default=\code{c(0,1)} in JMLE
#' @param bounds_t the bounds of ability parameters
#' @param bounds_a the bounds of discrimination parameters
#' @param bounds_b the bounds of difficulty parameters
#' @param bounds_c the bounds of guessing parameters
#' @param priors prior distributions, a list
#' @param decay decay rate, default=1
#' @param verbose TRUE to print details for debugging
#' @param true_params a list of true parameters for evaluating the parameter recovery
#' @return \code{model_3pl_jmle} returns estimated t, a, b, c parameters in a list
#' @examples
#' \donttest{
#' # generate data
#' x <- model_3pl_gendata(2000, 40)
#' # free estimation, 40 iterations
#' y <- model_3pl_jmle(x$u, true_params=x, iter=40, verbose=TRUE)
#' # fix c-parameters, 40 iterations
#' y <- model_3pl_jmle(x$u, c=0, true_params=x, iter=40)
#' }
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_3pl_jmle <- function(u, t=NA, a=NA, b=NA, c=NA, D=1.702, iter=100, conv=1e-3, nr_iter=10, scale=c(0, 1), bounds_t=c(-4, 4), bounds_a=c(.01, 2.5), bounds_b=c(-4, 4), bounds_c=c(0, .4), priors=list(t=c(0, 1)), decay=1, verbose=FALSE, true_params=NULL){
  # internal config
  h_max <- 1
  tracking <- list(fit=rep(NA, iter), t=rep(NA, iter), a=rep(NA, iter), b=rep(NA, iter), c=rep(NA, iter))
  
  # initial values
  n_p <- nrow(u)
  n_i <- ncol(u)
  if(length(t) == 1) t <- rep(t, n_p)
  t[t_free <- is.na(t)] <- rnorm(sum(is.na(t)), 0, .01)
  if(length(a) == 1) a <- rep(a, n_i)
  a[a_free <- is.na(a)] <- rlnorm(sum(is.na(a)), -.1, .01)
  if(length(b) == 1) b <- rep(b, n_i)
  b[b_free <- is.na(b)] <- rnorm(sum(is.na(b)), 0, .01)
  if(length(c) == 1) c <- rep(c, n_i)
  c[c_free <- is.na(c)] <- rbeta(sum(is.na(c)), 5, 50)

  for(k in 1:iter){
    # max change in parameters
    max_absh <- 0
    
    # t parameters
    if(any(t_free)){
      for(m in 1:nr_iter){
        dv_t <- model_3pl_dv_jmle(model_3pl_dv_Pt, u, t, a, b, c, D)
        dv_t$dv1 <- colSums(dv_t$dv1, na.rm=TRUE)
        dv_t$dv2 <- colSums(dv_t$dv2, na.rm=TRUE)
        if(!is.null(priors$t)){
          dv_t$dv1 <- dv_t$dv1 - (t - priors$t[1]) / priors$t[2]^2
          dv_t$dv2 <- dv_t$dv2 - 1 / priors$t[2]^2
        }
        nr_t <- nr_iteration(t, t_free, dv_t, h_max, decay, bounds_t)
        t <- nr_t$param
        if(max(abs(nr_t$h[t_free])) < conv) break
      }
      max_absh <- max(abs(nr_t$h[t_free]), max_absh)
      # rescale thetas
      if(!is.null(scale))
        t <- (t - mean(t)) / sd(t) * scale[2] + scale[1]
    }
    
    # b parameters
    if(any(b_free)){
      for(m in 1:nr_iter){
        dv_b <- model_3pl_dv_jmle(model_3pl_dv_Pb, u, t, a, b, c, D)
        dv_b$dv1 <- rowSums(dv_b$dv1, na.rm=TRUE)
        dv_b$dv2 <- rowSums(dv_b$dv2, na.rm=TRUE)
        if(!is.null(priors$b)){
          dv_b$dv1 <- dv_b$dv1 - (b - priors$b[1]) / priors$b[2]^2
          dv_b$dv2 <- dv_b$dv2 - 1 / priors$b[2]^2
        }
        nr_b <- nr_iteration(b, b_free, dv_b, h_max, decay, bounds_b)
        b <- nr_b$param
        if(max(abs(nr_b$h[b_free])) < conv) break
      }
      max_absh <- max(abs(nr_b$h[b_free]), max_absh)
    }

    # a parameters
    if(any(a_free)){
      for(m in 1:nr_iter){
        dv_a <- model_3pl_dv_jmle(model_3pl_dv_Pa, u, t, a, b, c, D)
        dv_a$dv1 <- rowSums(dv_a$dv1, na.rm=TRUE)
        dv_a$dv2 <- rowSums(dv_a$dv2, na.rm=TRUE)
        if(!is.null(priors$a)){
          dv_a$dv1 <- dv_a$dv1 - 1/a * (1 + (log(a)-priors$a[1])/priors$a[2]^2)
          dv_a$dv2 <- dv_a$dv2 - 1/a^2 * (1/priors$a[2]^2 - (1 + (log(a)-priors$a[1])/priors$a[2]^2))
        }
        nr_a <- nr_iteration(a, a_free, dv_a, h_max, decay, bounds_a)
        a <- nr_a$param
        if(max(abs(nr_a$h[a_free])) < conv) break
      }
      max_absh <- max(abs(nr_a$h[a_free]), max_absh)
    }
    
    # estimate c parameters
    if(any(c_free)){
      for(m in 1:nr_iter){
        dv_c <- model_3pl_dv_jmle(model_3pl_dv_Pc, u, t, a, b, c, D)
        dv_c$dv1 <- rowSums(dv_c$dv1, na.rm=TRUE)
        dv_c$dv2 <- rowSums(dv_c$dv2, na.rm=TRUE)
        if(!is.null(priors$c)){
          dv_c$dv1 <- dv_c$dv1 - ((priors$c[2]-1)/(1-c) - (priors$c[1]-1)/c)
          dv_c$dv2 <- dv_c$dv2 - ((priors$c[1]-1)/c^2 + (priors$c[2]-1)/(1-c)^2)
        }
        nr_c <- nr_iteration(c, c_free, dv_c, h_max, decay, bounds_c)
        c <- nr_c$param
        if(max(abs(nr_c$h[c_free])) < conv) break
      }
      max_absh <- max(abs(nr_c$h[c_free]), max_absh)
    }
    
    decay <- decay * decay
    
    # model fit
    if(verbose) {
      loglh <- -2 * sum(model_3pl_lh(u, t, a, b, c, D, log=TRUE), na.rm=TRUE)
      cat('iter #', k, ': -2 log-likelihood = ', round(loglh, 2), ', max_change = ', round(max_absh, 4), '\n', sep='')
      tracking$fit[k] <- loglh
      if(any(t_free)) tracking$t[k] <- mean(abs(nr_t$h[t_free]))
      if(any(a_free)) tracking$a[k] <- mean(abs(nr_a$h[a_free]))
      if(any(b_free)) tracking$b[k] <- mean(abs(nr_b$h[b_free]))
      if(any(c_free)) tracking$c[k] <- mean(abs(nr_c$h[c_free]))
    }
    if(max_absh < conv) break
  }
  
  # debugging
  if(verbose)
    estimate_3pl_debug(tracking, k)
  
  # compare with true parameters
  if(!is.null(true_params))
    estimate_3pl_eval(true_params, t, a, b, c, t_free, a_free, b_free, c_free)

  list(t=t, a=a, b=b, c=c)
}


#' @rdname estimate_3pl
#' @param pdv_fn the function to compute derivatives of P w.r.t the estimating parameters 
#' @keywords internal
model_3pl_dv_mmle <- function(pdv_fn, u, quad, a, b, c, D){
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_q <- length(quad$t)
  pdv <- pdv_fn(quad$t, a, b, c, D)
  
  x <- model_3pl_dvC(u, quad, a, b, c, D)
  p0 <- array(unlist(x[[1]]), dim=c(n_p, n_i, n_q))
  p1 <- x[[2]]
  p2 <- x[[3]]
  
  dv1 <- dv2 <- matrix(0, nrow=n_i, ncol=n_p)
  for(k in 1:n_q) {
    dv0 <- t(quad$w[k] * p1[,k] / p2 / p0[,,k] * (-1)^(u+1))
    dv1 <- dv1 + dv0 * pdv$dv1[,k]
    dv2 <- dv2 + dv0 * pdv$dv2[,k]
  }
  dv2 <- rowSums(dv2 - dv1^2, na.rm=TRUE)
  dv1 <- rowSums(dv1, na.rm=TRUE)
  list(dv1=dv1, dv2=dv2)
}

#' @rdname estimate_3pl
#' @description \code{model_3pl_mmle} estimates the item parameters using the 
#' marginal maximum likelihood estimation (MMLE) method
#' @param quad the number of quadrature points
#' @param score_fn the scoring function: 'eap' or 'map'
#' @return \code{model_3pl_mmle} returns estimated t, a, b, c parameters in a list
#' @examples
#' \donttest{
#' # generate data
#' x <- model_3pl_gendata(2000, 40)
#' # free estimation, 40 iterations
#' y <- model_3pl_mmle(x$u, true_params=x, iter=40, verbose=TRUE)
#' }
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_3pl_mmle <- function(u, t=NA, a=NA, b=NA, c=NA, D=1.702, iter=100, conv=1e-3, nr_iter=10, bounds_t=c(-4, 4), bounds_a=c(.01, 2.5), bounds_b=c(-4, 4), bounds_c=c(0, .4), priors=list(t=c(0, 1)), decay=1, quad='11', score_fn=c('eap', 'map'), verbose=FALSE, true_params=NULL){
  # internal config
  h_max <- 1
  if(is.null(priors$t)) priors$t <- c(0, 1)
  quad <- hermite_gauss(quad)
  quad$w <- quad$w * exp(quad$t^2) * dnorm(quad$t, priors$t[1], priors$t[2])
  score_fn <- switch(match.arg(score_fn, score_fn), 'eap'=model_3pl_eap, 'map'=model_3pl_map)
  tracking <- list(fit=rep(NA, iter), t=rep(NA, iter), a=rep(NA, iter), b=rep(NA, iter), c=rep(NA, iter))
  
  # initial values
  n_p <- nrow(u)
  n_i <- ncol(u)
  if(length(t) == 1) t <- rep(t, n_p)
  t[t_free <- is.na(t)] <- rnorm(sum(is.na(t)), 0, .01)
  if(length(a) == 1) a <- rep(a, n_i)
  a[a_free <- is.na(a)] <- rlnorm(sum(is.na(a)), -.1, .01)
  if(length(b) == 1) b <- rep(b, n_i)
  b[b_free <- is.na(b)] <- rnorm(sum(is.na(b)), 0, .01)
  if(length(c) == 1) c <- rep(c, n_i)
  c[c_free <- is.na(c)] <- rbeta(sum(is.na(c)), 5, 50)
  
  for(k in 1:iter){
    # max change in parameters
    max_absh <- 0
    
    # b parameters
    if(any(b_free)){
      for(m in 1:nr_iter){
        dv_b <- model_3pl_dv_mmle(model_3pl_dv_Pb, u, quad, a, b, c, D)
        if(!is.null(priors$b)){
          dv_b$dv1 <- dv_b$dv1 - (b - priors$b[1]) / priors$b[2]^2
          dv_b$dv2 <- dv_b$dv2 - 1 / priors$b[2]^2
        }
        nr_b <- nr_iteration(b, b_free, dv_b, h_max, decay, bounds_b)
        b <- nr_b$param
        if(max(abs(nr_b$h[b_free])) < conv) break
      }
      max_absh <- max(abs(nr_b$h[b_free]), max_absh)
    }
    
    # a parameters
    if(any(a_free)){
      for(m in 1:nr_iter){
        dv_a <- model_3pl_dv_mmle(model_3pl_dv_Pa, u, quad, a, b, c, D)
        if(!is.null(priors$a)){
          dv_a$dv1 <- dv_a$dv1 - 1/a * (1 + (log(a)-priors$a[1])/priors$a[2]^2)
          dv_a$dv2 <- dv_a$dv2 - 1/a^2 * (1/priors$a[2]^2 - (1 + (log(a)-priors$a[1])/priors$a[2]^2))
        }
        nr_a <- nr_iteration(a, a_free, dv_a, h_max, decay * .2, bounds_a)
        a <- nr_a$param
        if(max(abs(nr_a$h[a_free])) < conv) break
      }
      max_absh <- max(abs(nr_a$h[a_free]), max_absh)
    }
    
    # estimate c parameters
    if(any(c_free)){
      for(m in 1:nr_iter){
        dv_c <- model_3pl_dv_mmle(model_3pl_dv_Pc, u, quad, a, b, c, D)
        if(!is.null(priors$c)){
          dv_c$dv1 <- dv_c$dv1 - ((priors$c[2]-1)/(1-c) - (priors$c[1]-1)/c)
          dv_c$dv2 <- dv_c$dv2 - ((priors$c[1]-1)/c^2 + (priors$c[2]-1)/(1-c)^2)
        }
        nr_c <- nr_iteration(c, c_free, dv_c, h_max, decay, bounds_c)
        c <- nr_c$param
        if(max(abs(nr_c$h[c_free])) < conv) break
      }
      max_absh <- max(abs(nr_c$h[c_free]), max_absh)
    }
    
    decay <- decay * decay
    
    # scoring
    if(any(t_free))
      t[t_free] <- score_fn(u, a, b, c, D, priors=priors$t, bounds_t=bounds_t)$t[t_free]

    # model fit
    if(verbose) {
      loglik <- -2 * sum(model_3pl_lh(u, t, a, b, c, D, log=TRUE), na.rm=TRUE)
      cat('iter #', k, ': -2 log-likelihood = ', round(loglik, 2), ', max_change = ', round(max_absh, 4), '\n', sep='')
      tracking$fit[k] <- loglik
      if(any(a_free)) tracking$a[k] <- mean(abs(nr_a$h[a_free]))
      if(any(b_free)) tracking$b[k] <- mean(abs(nr_b$h[b_free]))
      if(any(c_free)) tracking$c[k] <- mean(abs(nr_c$h[c_free]))
    }
    if(max_absh < conv)
      break
  }
  
  # debugging
  if(verbose)
    estimate_3pl_debug(tracking, k)

  # compare with true parameters
  if(!is.null(true_params))
    estimate_3pl_eval(true_params, t, a, b, c, t_free, a_free, b_free, c_free)
  
  list(t=t, a=a, b=b, c=c)
}


#' @rdname estimate_3pl
#' @param index the indices of items being plotted
#' @param intervals intervals on the x-axis
#' @return \code{model_3pl_fitplot} returns a \code{ggplot} object
#' @examples 
#' with(model_3pl_gendata(1000, 20), 
#'      model_3pl_fitplot(u, t, a, b, c, index=c(1, 3, 5)))
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_3pl_fitplot <- function(u, t, a, b, c, D=1.702, index=NULL, intervals=seq(-3, 3, .5)){
  if(is.null(index)) index <- seq(b)
  labels <- (intervals[-length(intervals)] + intervals[-1]) / 2
  groups <- cut(t, intervals, labels=labels)
  
  obs <- aggregate(u, by=list(intervals=groups), mean, na.rm=TRUE)[, c(1, index+1)]
  obs <- melt(obs, id.vars='intervals', variable.name='items')
  obs[, 'type'] <- 'Observed'
  p <- model_3pl_prob(t, a, b, c, D)
  exp <- aggregate(p, by=list(intervals=groups), mean, na.rm=TRUE)[, c(1, index+1)]
  exp <- melt(exp, id.vars='intervals', variable.name='items')
  exp[, 'type'] <- 'Expected'
  data <- rbind(obs, exp)
  data$intervals <- as.numeric(levels(data$intervals)[data$intervals])
  levels(data$items) <- gsub('V', 'Item ', levels(data$items))
  
  ggplot(data, aes_string('intervals', 'value', color='type', group='type')) + 
    geom_line() + facet_wrap(~items) + xlab(expression(theta)) + ylab('Probability') + 
    scale_color_discrete(guide=guide_legend("")) + theme_bw()
}
