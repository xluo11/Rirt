#' Estimation of the Mixed Format Model
#' @description Estimate the mixed format model
#' @name estimate_mixed
NULL


#' @rdname estimate_mixed
#' @param u the response data, 2d marix
#' @param items a list of 3pl, gpcm, grm items
#' @param D the scaling constant
#' @param priors the prior distribution
#' @param bounds_t the lower- and upper-bound of the parameter
#' @return \code{model_mixed_eap} returns a list of point estimates and
#' standard error of the ability parameters
#' @examples 
#' x <- model_mixed_gendata(200, 30, 5, 5, 3)
#' y <- model_mixed_eap(x$u, x$items)
#' c('corr'=cor(x$t, y$t), 'rmse'=rmse(x$t, y$t))
#' @export
model_mixed_eap <- function(u, items, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4)) {
  n_p <- dim(u)[1]
  n_i <- sapply(names(items), function(x) ifelse(is.null(items[[x]]), 0, nrow(items[[x]])))
  
  u_ix1 <- cumsum(n_i)
  u_ix0 <- c(0, u_ix1[-length(u_ix1)]) + 1
  u_ix1[n_i == 0] <- 0
  u_ix0[n_i == 0] <- 0
  u <- Map(function(x, y) u[, x:y], u_ix0, u_ix1)
  names(u) <- names(items)
  
  quad <- hermite_gauss('11')
  quad$w <- quad$w * exp(quad$t^2) * dnorm(quad$t, priors[1], priors[2])
  n_q <- length(quad$t)
  
  p <- model_mixed_prob(quad$t, items, D)
  lh_3pl <- lh_gpcm <- lh_grm <- 1
  
  if(!is.null(items$'3pl')) {
    lh_3pl <- ifelse(is.na(u$'3pl'), 0, u$'3pl') %*% t(log(p$'3pl')) + (1 - ifelse(is.na(u$'3pl'), 1, u$'3pl')) %*% t(log(1 - p$'3pl'))
    lh_3pl <- exp(lh_3pl)  
  }
  
  if(!is.null(items$'gpcm')) {
    uix_gpcm <- model_polytomous_3dindex(u$'gpcm')[,-1]
    lh_gpcm <- array(NA, dim=c(n_p, n_i['gpcm'], n_q))
    for(i in 1:n_q)
      lh_gpcm[,,i] <- array(p$'gpcm'[i,,][uix_gpcm], dim=c(n_p, n_i['gpcm']))
    lh_gpcm <- exp(apply(log(lh_gpcm), 3, rowSums, na.rm=TRUE))
  }
  
  if(!is.null(items$'grm')) {
    uix_grm <- model_polytomous_3dindex(u$'grm')[,-1]
    lh_grm <- array(NA, dim=c(n_p, n_i['grm'], n_q))
    for(i in 1:n_q)
      lh_grm[,,i] <- array(p$'grm'[i,,][uix_grm], dim=c(n_p, n_i['grm']))
    lh_grm <- exp(apply(log(lh_grm), 3, rowSums, na.rm=TRUE))
  }
  
  lh <- lh_3pl * lh_gpcm * lh_grm
  t <- ((lh / (lh %*% quad$w)[,1]) %*% (quad$w * quad$t))[,1]
  t[t < bounds_t[1]] <- bounds_t[1]
  t[t > bounds_t[2]] <- bounds_t[2]
  t_sd <- ((lh / (lh %*% quad$w)[,1] * outer(t, quad$t, '-')^2) %*% quad$w)[,1]
  t_sd <- sqrt(t_sd)
  list(t=t, sd=t_sd)
}


#' @rdname estimate_mixed
#' @param iter the maximum number of newton-raphson iterations
#' @param conv the convergence criterion
#' @return \code{model_mixed_map} returns a list of point estimates of the ability parameters
#' @examples 
#' x <- model_mixed_gendata(200, 30, 5, 5, 3)
#' y <- model_mixed_map(x$u, x$items)
#' c('corr'=cor(x$t, y$t), 'rmse'=rmse(x$t, y$t))
#' @export
model_mixed_map <- function(u, items, D=1.702, priors=c(0, 1), bounds_t=c(-4, 4), iter=30, conv=1e-3) {
  n_p <- dim(u)[1]
  n_i <- sapply(names(items), function(x) ifelse(is.null(items[[x]]), 0, nrow(items[[x]])))
  
  u_ix1 <- cumsum(n_i)
  u_ix0 <- c(0, u_ix1[-length(u_ix1)]) + 1
  u_ix1[n_i == 0] <- 0
  u_ix0[n_i == 0] <- 0
  u <- Map(function(x, y) u[, x:y], u_ix0, u_ix1)
  names(u) <- names(items)
  
  if(!is.null(items$'gpcm'))
    uix_gpcm <- model_polytomous_3dindex(u$'gpcm')
  if(!is.null(items$'grm'))
    uix_grm <- model_polytomous_3dindex(u$'grm')
  
  t <- rnorm(n_p, 0, .01)
  t_free <- rep(TRUE, n_p)
  for(i in 1:iter){
    dv_3pl <- dv_gpcm <- dv_grm <- numeric(0)
    if(!is.null(items$'3pl'))
      dv_3pl <- with(items$'3pl', model_3pl_dv_jmle(model_3pl_dv_Pt, u$'3pl', t, a, b, c, D))
    if(!is.null(items$'gpcm')){
      d <- as.matrix(items$'gpcm'[, grep('^d[0-9]+$', colnames(items$'gpcm'), value=TRUE)])
      dv_gpcm <- with(items$'gpcm', model_gpcm_dv_jmle(uix_gpcm, model_gpcm_dv_Pt(t, a, b, d, D)))
    }
    if(!is.null(items$'grm')){
      b <- as.matrix(items$'grm'[, grep('^b[0-9]+$', colnames(items$'grm'), value=TRUE)])
      dv_grm <- with(items$'grm', model_grm_dv_jmle(uix_grm, model_grm_dv_Pt(t, a, b, D)))
    }
    dv1 <- rowSums(cbind(t(dv_3pl$dv1), dv_gpcm$dv1, dv_grm$dv1), na.rm=TRUE)
    dv2 <- rowSums(cbind(t(dv_3pl$dv2), dv_gpcm$dv2, dv_grm$dv2), na.rm=TRUE)
    if(!is.null(priors)){
      dv1 <- dv1 - (t - priors[1]) / priors[2]^2
      dv2 <- dv2 - 1 / priors[2]^2      
    }
    nr <- nr_iteration(t, t_free, list(dv1=dv1, dv2=dv2), 1.0, 1.0, bounds_t)
    t <- nr$param
    if(max(abs(nr$h)) < conv) break
  }
  
  list(t=t)
}
