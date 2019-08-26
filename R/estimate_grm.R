#' Estimation of the Graded Response Model
#' @description Estimate the GRM using the joint or marginal
#' maximum likelihood estimation 
#' @name estimate_grm
NULL

#' @rdname estimate_grm
#' @description \code{model_grm_eap} scores response vectors using the EAP method
#' @return \code{model_grm_eap} returns theta estimates and standard errors in a list
#' @examples 
#' with(model_grm_gendata(10, 50, 3), 
#'      cbind(true=t, est=model_grm_eap(u, a, b)$t))
#' @importFrom stats dnorm 
#' @export
model_grm_eap <- function(u, a, b, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4)){
  quad <- hermite_gauss('11')
  quad$w <- quad$w * exp(quad$t^2) * dnorm(quad$t, priors[1], priors[2])
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_q <- length(quad$t)
  
  p <- model_grm_prob(quad$t, a, b, D)
  ix <- model_polytomous_3dindex(u)
  lh <- array(NA, c(n_p, n_i, n_q))
  for(q in 1:n_q)
    lh[,,q] <- array(p[q,,][ix[,-1]], c(n_p, n_i))
  lh <- apply(lh, c(1, 3), prod, na.rm=T)
  t <- ((lh / (lh %*% quad$w)[,1]) %*% (quad$w * quad$t))[,1]
  t[t < bounds_t[1]] <- bounds_t[1]
  t[t > bounds_t[2]] <- bounds_t[2]
  t_sd <- ((lh / (lh %*% quad$w)[,1] * outer(t, quad$t, '-')^2) %*% quad$w)[,1]
  t_sd <- sqrt(t_sd)
  list(t=t, sd=t_sd)
}

#' @rdname estimate_grm
#' @description \code{model_grm_map} scores response vectors using the MAP method
#' @return \code{model_grm_map} returns theta estimates in a list
#' @examples 
#' with(model_grm_gendata(10, 50, 3), 
#'      cbind(true=t, est=model_grm_map(u, a, b)$t))
#' @export
model_grm_map <- function(u, a, b, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4), iter=30, conv=1e-3){
  ix <- model_polytomous_3dindex(u)
  t <- rnorm(dim(u)[1], 0, .01)
  t_free <- rep(T, length(t))
  for(m in 1:iter){
    dv <- model_grm_dv_jmle(ix, model_grm_dv_Pt(t, a, b, D))
    dv$dv1 <- rowSums(dv$dv1, na.rm=T)
    dv$dv2 <- rowSums(dv$dv2, na.rm=T)
    if(!is.null(priors)){
      dv$dv1 <- dv$dv1 - (t - priors[1]) / priors[2]^2
      dv$dv2 <- dv$dv2 - 1 / priors[2]^2
    }
    nr <- nr_iteration(t, t_free, dv, 1.0, 1.0, bounds_t)
    t <- nr$param
    if(max(abs(nr$h)) < conv) break
  }
  list(t=t)
}

#' @rdname estimate_grm
#' @keywords internal
model_grm_dv_Pt <- function(t, a, b, D){
  n_c <- ncol(b) + 1
  p <- model_grm_prob(t, a, b, D, raw=T)
  dv1 <- aperm(aperm(p*(1-p), c(2,1,3)) * D * a, c(2,1,3))
  dv1 <- dv1[,,1:n_c] - dv1[,,-1]
  dv2 <- aperm(aperm(p*(1-p)*(1-2*p), c(2,1,3)) * (D * a)^2, c(2,1,3))
  dv2 <- dv2[,,1:n_c] - dv2[,,-1]
  p <- p[,,1:n_c] - p[,,-1]
  list(p=p, dv1=dv1, dv2=dv2)
}

#' @rdname estimate_grm
#' @keywords internal
model_grm_dv_Pa <- function(t, a, b, D){
  n_c <- ncol(b) + 1
  p <- model_grm_prob(t, a, b, D, raw=T)
  term0 <- D * outer(t, cbind(0, b, 0), '-')
  dv1 <- p * (1-p) * term0
  dv1 <- dv1[,,1:n_c] - dv1[,,-1]
  dv2 <- p * (1-p) * (1-2*p) * term0^2
  dv2 <- dv2[,,1:n_c] - dv2[,,-1]
  p <- p[,,1:n_c] - p[,,-1]
  list(p=p, dv1=dv1, dv2=dv2)
}

#' @rdname estimate_grm
#' @keywords internal
model_grm_dv_Pb <- function(t, a, b, D){
  n_p <- length(t)
  n_i <- nrow(b)
  n_c <- ncol(b) + 1
  p <- model_grm_prob(t, a, b, D, raw=T)
  dv1 <- dv2 <- array(0, c(n_p, n_i, n_c, n_c-1))
  for(k in 1:(n_c-1)){
    term0 <- t(t(p[,,k+1]*(1-p[,,k+1])) * (-D * a))
    dv1[,,k,k] <- -1 * term0
    dv1[,,k+1,k] <- term0
    term1<- t(t(p[,,k+1]*(1-p[,,k+1])*(1-2*p[,,k+1])) * (D*a)^2)
    dv2[,,k,k] <- -1 * term1
    dv2[,,k+1,k] <- term1
  }
  p <- p[,,1:n_c] - p[,,-1]
  list(p=p, dv1=dv1, dv2=dv2)
}

#' @rdname estimate_grm
#' @param u_ix the 3d indices
#' @param dvp the derivatives of P
#' @keywords internal
model_grm_dv_jmle <- function(u_ix, dvp){
  n_p <- max(u_ix[,1])
  n_i <- max(u_ix[,2])
  dv1 <- array(with(dvp, dv1[u_ix]/p[u_ix]), c(n_p, n_i))
  dv2 <- array(with(dvp, dv2[u_ix]/p[u_ix]), c(n_p, n_i)) - dv1^2
  list(dv1=dv1, dv2=dv2)  
}

#' @rdname estimate_grm
#' @description \code{model_grm_jmle} estimates the parameters using the 
#' joint maximum likelihood estimation (JMLE) method
#' @param u the observed response matrix, 2d matrix
#' @param t ability parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param a discrimination parameters, 1d vector (fixed value) or NA (freely estimate)
#' @param b difficulty parameters, 2d matrix (fixed value) or NA (freely estimate)
#' @param D the scaling constant, 1.702 by default
#' @param iter the maximum iterations
#' @param conv the convergence criterion for the -2 log-likelihood
#' @param nr_iter the maximum newton-raphson iterations, default=10
#' @param scale the scale of theta parameters
#' @param bounds_t bounds of ability parameters
#' @param bounds_a bounds of discrimination parameters
#' @param bounds_b bounds of location parameters
#' @param priors a list of prior distributions
#' @param decay decay rate
#' @param verbose TRUE to print debuggin information
#' @param true_params a list of true parameters for evaluating the estimation accuracy
#' @return \code{model_grm_jmle} returns estimated t, a, b parameters in a list
#' @examples
#' \donttest{
#' # generate data
#' x <- model_grm_gendata(1000, 40, 3)
#' # free calibration, 40 iterations
#' y <- model_grm_jmle(x$u, true_params=x, iter=40, verbose=TRUE)
#' }
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_grm_jmle <- function(u, t=NA, a=NA, b=NA, D=1.702, iter=100, nr_iter=10, conv=1e-3, scale=c(0, 1), bounds_t=c(-4, 4), bounds_a=c(.01, 2.5), bounds_b=c(-4, 4), priors=list(t=c(0, 1)), decay=1, verbose=FALSE, true_params=NULL){
  # configuration
  h_max <- 1.0
  tracking <- list(fit=rep(NA, iter), t=rep(NA, iter), a=rep(NA, iter), b=rep(NA, iter))
  
  # initial values
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_c <- max(u, na.rm=T) + 1
  u_ix <- model_polytomous_3dindex(u)
  if(length(t) == 1) t <- rep(t, n_p)
  t[t_free <- is.na(t)] <- rnorm(sum(is.na(t)), 0, .01)
  if(length(a) == 1) a <- rep(a, n_i)
  a[a_free <- is.na(a)] <- rlnorm(sum(is.na(a)), -.1, .01)
  if(length(b) == 1) b <- array(b, c(n_i, n_c-1))
  b[b_free <- is.na(b)] <- rnorm(sum(is.na(b)), 0, .01)
  b <- t(apply(b, 1, sort))

  # estimate parameters
  for (k in 1:iter){
    # max change in parameters
     max_absh <- 0
    
    # t parameters
    if(any(t_free)){
      for(j in 1:nr_iter){
        dv_t <- model_grm_dv_jmle(u_ix, model_grm_dv_Pt(t, a, b, D))
        dv_t$dv1 <- rowSums(dv_t$dv1, na.rm=T)
        dv_t$dv2 <- rowSums(dv_t$dv2, na.rm=T)
        if(!is.null(priors$t)){
          dv_t$dv1 <- dv_t$dv1 - (t - priors$t[1]) / priors$t[2]^2
          dv_t$dv2 <- dv_t$dv2 - 1 / priors$t[2]^2
        }
        nr_t <- nr_iteration(t, t_free, dv_t, h_max, decay, bounds_t)
        t <- nr_t$param
        if(max(abs(nr_t$h[t_free])) < conv) break
      }
      max_absh <- max(max_absh, abs(nr_t$h[t_free]))
      # rescale theta
      if(!is.null(scale)) t <- (t - mean(t)) / sd(t) * scale[2] + scale[1]
    }
    
    # b parameters
    if(any(b_free)){
      for(j in 1:nr_iter){
        dv_b <- model_grm_dv_Pb(t, a, b, D)
        dv_bh <- array(0, c(n_i, n_c-1))
        for(m in 1:(n_c-1)){
          dv <- model_grm_dv_jmle(u_ix, with(dv_b, list(p=p, dv1=dv1[,,,m], dv2=dv2[,,,m])))
          dv$dv1 <- colSums(dv$dv1, na.rm=T)
          dv$dv2 <- colSums(dv$dv2, na.rm=T)
          if(!is.null(priors$b)){
            dv$dv1 <- dv$dv1 - (b[,m] - priors$b[1]) / priors$b[2]^2
            dv$dv2 <- dv$dv2 - 1 / priors$b[2]^2
          }
          nr <- nr_iteration(b[,m], b_free[,m], dv, h_max, decay, bounds_b)
          b[,m] <- nr$param
          dv_bh[,m] <- nr$h
        }
        b <- t(apply(b, 1, sort))
        if(max(abs(dv_bh[b_free])) < conv) break
      }
      max_absh <- max(max_absh, abs(dv_bh[b_free]))
    }
    
    # a parameters
    if(any(a_free)){
      for(j in 1:nr_iter){
        dv_a <- model_grm_dv_jmle(u_ix, model_grm_dv_Pa(t, a, b, D))
        dv_a$dv1 <- colSums(dv_a$dv1, na.rm=T)
        dv_a$dv2 <- colSums(dv_a$dv2, na.rm=T)
        if(!is.null(priors$a)){
          dv_a$dv1 <- dv_a$dv1 - 1/a * (1 + (log(a)-priors$a[1])/priors$a[2]^2)
          dv_a$dv2 <- dv_a$dv2 - 1/a^2 * (1/priors$a[2]^2 - (1 + (log(a)-priors$a[1])/priors$a[2]^2))
        }
        nr_a <- nr_iteration(a, a_free, dv_a, h_max, decay, bounds_a)
        a <- nr_a$param
        if(max(abs(nr_a$h[a_free])) < conv) break
      }
      max_absh <- max(max_absh, abs(nr_a$h[a_free]))
    }
    
    decay <- decay * decay
    
    # model fit
    if(verbose) {
      loglh <- -2 * sum(model_grm_lh(u, t, a, b, D, log=T), na.rm=T)
      cat('iter #', k, ': -2 log-likelihood = ', round(loglh, 2), ', max_change = ', round(max_absh, 4), '\n', sep='')
      tracking$fit[k] <- loglh
      if(any(t_free)) tracking$t[k] <- mean(abs(nr_t$h[t_free]))
      if(any(a_free)) tracking$a[k] <- mean(abs(nr_a$h[a_free]))
      if(any(b_free)) tracking$b[k] <- mean(abs(dv_bh[b_free]))      
    }
    if(max_absh < conv) 
      break
  }
  
  # debugging
  if(verbose)
    estimate_grm_debug(tracking, k)
  
  # compare with true parameters
  if(!is.null(true_params))
    estimate_grm_eval(true_params, n_c, t, a, b, t_free, a_free, b_free)
  
  list(t=t, a=a, b=b)
}

#' @rdname estimate_grm
#' @keywords internal
model_grm_dv_mmle <- function(u_ix, quad,  pdv){
  n_p <- max(u_ix[,1])
  n_i <- max(u_ix[,2])
  n_q <- length(quad$t)
  p0 <- array(NA, c(n_p, n_i, n_q))
  for(q in 1:n_q)
    p0[,,q] <- array(pdv$p[q,,][u_ix[,-1]], c(n_p, n_i))
  p1 <- apply(p0, c(1, 3), prod, na.rm=T)
  p2 <- (p1 %*% quad$w)[,1]
  
  dv1 <- dv2 <- array(0, c(n_p, n_i))
  dv_common <- t(t(p1 / p2) * quad$w)
  for(q in 1:n_q) {
    dv1 <- dv1 + dv_common[,q] / p0[,,q] * array(pdv$dv1[q,,][u_ix[,-1]], c(n_p, n_i))
    dv2 <- dv2 + dv_common[,q] / p0[,,q] * array(pdv$dv2[q,,][u_ix[,-1]], c(n_p, n_i))
  }
  dv2 <- dv2 - dv1^2
  list(dv1=dv1, dv2=dv2)          
}


#' @rdname estimate_grm
#' @description \code{model_grm_mmle} estimates the parameters using the 
#' marginal maximum likelihood estimation (MMLE) method
#' @param quad_degree the number of quadrature points
#' @param score_fn the scoring method: 'eap' or 'map'
#' @return \code{model_grm_mmle} returns estimated t, a, b parameters in a list
#' @examples
#' \donttest{
#' # generate data
#' x <- model_grm_gendata(1000, 40, 3)
#' # free estimation, 40 iterations
#' y <- model_grm_mmle(x$u, true_params=x, iter=40, verbose=TRUE)
#' }
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_grm_mmle <- function(u, t=NA, a=NA, b=NA, d=NA, D=1.702, iter=100, nr_iter=10, conv=1e-3, bounds_t=c(-4, 4), bounds_a=c(.01, 2.5), bounds_b=c(-4, 4), priors=list(t=c(0, 1)), decay=1, quad_degree='11', score_fn=c('eap', 'map'), verbose=FALSE, true_params=NULL){
  # configuration
  h_max <- 1.0
  score_fn <- switch(match.arg(score_fn, score_fn), 'eap'=model_grm_eap, 'map'=model_grm_map)
  if(is.null(priors$t)) priors$t <- c(0, 1)
  quad <- hermite_gauss(quad_degree)
  quad$w <- quad$w * exp(quad$t^2) * dnorm(quad$t, priors$t[1], priors$t[2])
  tracking <- list(fit=rep(NA, iter), t=rep(NA, iter), a=rep(NA, iter), b=rep(NA, iter), d=rep(NA, iter))
  
  # initial values
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_c <- max(u, na.rm=T) + 1
  u_ix <- model_polytomous_3dindex(u)
  if(length(t) == 1) t <- rep(t, n_p)
  t[t_free <- is.na(t)] <- rnorm(sum(is.na(t)), 0, .01)
  if(length(a) == 1) a <- rep(a, n_i)
  a[a_free <- is.na(a)] <- rlnorm(sum(is.na(a)), -.1, .01)
  if(length(b) == 1) b <- array(b, c(n_i, n_c-1))
  b[b_free <- is.na(b)] <- rnorm(sum(is.na(b)), 0, .1)
  b <- t(apply(b, 1, sort))
  
  # estimate parameters
  for (k in 1:iter){
    # max change in parameters
    max_absh <- 0
    
    # b parameters
    if(any(b_free)){
      for(j in 1:nr_iter){
        dv_b <- model_grm_dv_Pb(quad$t, a, b, D)
        dv_bh <- array(0, c(n_i, n_c-1))
        for(m in 1:(n_c-1)){
          dv <- model_grm_dv_mmle(u_ix, quad, with(dv_b, list(p=p, dv1=dv1[,,,m], dv2=dv2[,,,m])))
          dv$dv1 <- colSums(dv$dv1, na.rm=T)
          dv$dv2 <- colSums(dv$dv2, na.rm=T)
          if(!is.null(priors$b)){
            dv$dv1 <- dv$dv1 - (b[,m] - priors$b[1]) / priors$b[2]^2
            dv$dv2 <- dv$dv2 - 1 / priors$b[2]^2
          }
          nr <- nr_iteration(b[,m], b_free[,m], dv, h_max, decay, bounds_b)
          b[,m] <- nr$param
          dv_bh[,m] <- nr$h
        }
        b <- t(apply(b, 1, sort))
        if(max(abs(dv_bh)) < conv) break
      }
      max_absh <- max(max_absh, abs(dv_bh[b_free]))
    }
    
    # a parameters
    if(any(a_free)){
      for(j in 1:nr_iter){
        dv_a <- model_grm_dv_mmle(u_ix, quad, model_grm_dv_Pa(quad$t, a, b, D))
        dv_a$dv1 <- colSums(dv_a$dv1, na.rm=T)
        dv_a$dv2 <- colSums(dv_a$dv2, na.rm=T)
        if(!is.null(priors$a)){
          dv_a$dv1 <- dv_a$dv1 - 1/a * (1 + (log(a)-priors$a[1])/priors$a[2]^2)
          dv_a$dv2 <- dv_a$dv2 - 1/a^2 * (1/priors$a[2]^2 - (1 + (log(a)-priors$a[1])/priors$a[2]^2))
        }
        nr_a <- nr_iteration(a, a_free, dv_a, h_max, decay, bounds_a)
        a <- nr_a$param
        if(max(abs(nr_a$h[a_free])) < conv) break
      }
      max_absh <- max(max_absh, abs(nr_a$h[a_free]))
    }
    
    # scoring
    if(any(t_free))
      t[t_free] <- score_fn(u, a, b, D, priors=priors$t, bounds_t=bounds_t)$t[t_free]
    
    decay <- decay * decay
    
    # model fit
    if(verbose) {
      loglh <- -2 * sum(model_grm_lh(u, t, a, b, D, log=T), na.rm=T)
      cat('iter #', k, ': -2 log-likelihood = ', round(loglh, 2), ', max_change = ', round(max_absh, 4), '\n', sep='')
      tracking$fit[k] <- loglh
      if(any(a_free)) tracking$a[k] <- mean(abs(nr_a$h[a_free]))
      if(any(b_free)) tracking$d[k] <- mean(abs(dv_bh[b_free]))
    }
    if(max_absh < conv)
      break
  }
  
  # debugging
  if(verbose)
    estimate_grm_debug(tracking, k)
  
  # compare with true parameters
  if(!is.null(true_params))
    estimate_grm_eval(true_params, n_c, t, a, b, t_free, a_free, b_free)
  
  list(t=t, a=a, b=b)
}

#' @rdname estimate_grm
#' @param index the indices of items being plotted
#' @param intervals intervals on the x-axis
#' @return \code{model_grm_fitplot} returns a \code{ggplot} object
#' @examples 
#' with(model_grm_gendata(1000, 20, 3), 
#'      model_grm_fitplot(u, t, a, b, index=c(1, 3, 5)))
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
model_grm_fitplot <- function(u, t, a, b, D=1.702, index=NULL, intervals=seq(-3, 3, .5)){
  if(is.null(index)) index <- seq(b)
  groups <- cut(t, intervals, labels=(intervals[-length(intervals)] + intervals[-1]) / 2)
  
  obs <- aggregate(u, by=list(intervals=groups), mean, na.rm=TRUE)[, c(1, index+1)]
  obs <- melt(obs, id.vars='intervals', variable.name='items')
  obs[, 'type'] <- 'Observed'
  p <- model_grm_prob(t, a, b, D)
  p <- apply(p, 1:2, function(x) sum(x * (seq(x)-1), na.rm=T))
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
