#' Mixed-format model
#' @description Common computations and operations for the mixed format model
#' @name model_mixed
NULL


#' @rdname model_mixed
#' @param n_p the number of test takers
#' @param n_3pl the number of 3pl items
#' @param n_gpcm the number of gpcm items
#' @param n_grm the number of grm items
#' @param n_c the number of score categories for polytomous items
#' @param ... additional arguments
#' @return \code{model_mixed_gendata} returns a list of generated responses, ability paramters and items
#' @examples 
#' # generate 10 3pl items, 5 gpcm items and 5 grm items
#' model_mixed_gendata(10, n_3pl=10, n_gpcm=5, n_grm=5, n_c=3)
#' # generate 5 gpcm items and 5 grm items, 4 score categories each
#' model_mixed_gendata(10, n_gpcm=5, n_grm=5, n_c=4)
#' @export
model_mixed_gendata <- function(n_p, n_3pl=0, n_gpcm=0, n_grm=0, n_c, ...) {
  u_3pl <- u_gpcm <- u_grm <- NULL
  i_3pl <- i_gpcm <- i_grm <- NULL
  
  opts <- list(...)
  t_dist <- if(is.null(opts$t_dist)) c(0, 1) else opts$t_dist
  t_bounds <- if(is.null(opts$t_bounds)) c(-3, 3) else opts$t_bounds
  t <- rnorm(n_p, mean=t_dist[1], sd=t_dist[2])
  t[t < t_bounds[1]] <- t_bounds[1]
  t[t > t_bounds[2]] <- t_bounds[2]
  
  
  if(!is.na(n_3pl) && n_3pl > 0) {
    x_3pl <- model_3pl_gendata(n_p, n_3pl, t=t, ...)
    u_3pl <- x_3pl$u
    i_3pl <- with(x_3pl, data.frame(a=a, b=b, c=c))
  }
  
  if(!is.na(n_gpcm) && n_gpcm > 0) {
    x_gpcm <- model_gpcm_gendata(n_p, n_gpcm, n_c, t=t, ...)
    u_gpcm <- x_gpcm$u
    i_gpcm <- with(x_gpcm, data.frame(a=a, b=b, d))
    colnames(i_gpcm) <- gsub('^X([0-9]+)$', 'd\\1', colnames(i_gpcm))
  }
  
  if(!is.na(n_grm) && n_grm > 0) {
    x_grm <- model_grm_gendata(n_p, n_grm, n_c, t=t, ...)
    u_grm <- x_grm$u
    i_grm <- with(x_grm, data.frame(a=a, b))
    colnames(i_grm) <- gsub('^X([0-9]+)$', 'b\\1', colnames(i_grm))
  }
  
  list(u=cbind(u_3pl, u_gpcm, u_grm), t=t, items=list('3pl'=i_3pl, 'gpcm'=i_gpcm, 'grm'=i_grm))
}

#' @rdname model_mixed
#' @param t ability parameters, a vector
#' @param items a list of '3pl', 'gpcm', and 'grm' items
#' @param D the scaling constant, default=1.702
#' @return \code{model_mixed_prob} returns a list of probabilities for '3pl', 'gpcm', and 'grm' items
#' @examples 
#' # generate 5 people and 4 items of each type
#' with(model_mixed_gendata(n_p=5, n_3pl=4, n_gpcm=4, n_grm=4, n_c=3), 
#'      model_mixed_prob(t, items))
#'      
#' # generate 10 people and 5 gpcm and 5 grm items
#' with(model_mixed_gendata(n_p=10, n_gpcm=4, n_grm=4, n_c=3), 
#'      model_mixed_prob(t, items))
#' @export
model_mixed_prob <- function(t, items, D=1.702){
  rs_3pl <- rs_gpcm <- rs_grm <- NULL
  
  if(!is.null(items$'3pl'))
    rs_3pl <- with(items$'3pl', model_3pl_prob(t, a, b, c, D))
  
  if(!is.null(items$'gpcm')) {
    d <- as.matrix(with(items, gpcm[, grep('^d[0-9]+$', colnames(gpcm), value=TRUE)]))
    rs_gpcm <- with(items, model_gpcm_prob(t, gpcm[, 'a'], gpcm[, 'b'], d, D))
  }
  
  if(!is.null(items$'grm')) {
    b <- as.matrix(with(items, grm[, grep('^b[0-9]+$', colnames(grm), value=TRUE)]))
    rs_grm <- with(items, model_grm_prob(t, grm[, 'a'], b, D))
  }
  
  rs <- list('3pl'=rs_3pl, 'gpcm'=rs_gpcm, 'grm'=rs_grm)
  rs <- rs[match(names(items), names(rs))]
  rs
}

#' @rdname model_mixed
#' @param combine \code{TRUE} to combine results from list to matrix
#' @return \code{model_mixed_info} returns a list or matrix of information 
#' @examples
#' with(model_mixed_gendata(10, 4, 4, 4, 3), model_mixed_info(t, items)) 
#' with(model_mixed_gendata(10, 0, 4, 4, 3), model_mixed_info(t, items)) 
#' @export
model_mixed_info <- function(t, items, D=1.702, combine=TRUE){
  rs_3pl <- rs_gpcm <- rs_grm <- NULL
  
  if(!is.null(items$'3pl'))
    rs_3pl <- with(items$'3pl', model_3pl_info(t, a, b, c, D))
  
  if(!is.null(items$'gpcm')) {
    d <- as.matrix(with(items, gpcm[, grep('^d[0-9]+$', colnames(gpcm), value=TRUE)]))
    rs_gpcm <- with(items, model_gpcm_info(t, gpcm[, 'a'], gpcm[, 'b'], d, D))
    rs_gpcm <- apply(rs_gpcm, 2, rowSums, na.rm=TRUE)
    if(is.null(dim(rs_gpcm)))
      rs_gpcm <- matrix(rs_gpcm, nrow=1)
  }
  
  if(!is.null(items$'grm')) {
    b <- as.matrix(with(items, grm[, grep('^b[0-9]+$', colnames(grm), value=TRUE)]))
    rs_grm <- with(items, model_grm_info(t, grm[, 'a'], b, D))
    rs_grm <- apply(rs_grm, 2, rowSums, na.rm=TRUE)
    if(is.null(dim(rs_grm)))
      rs_grm <- matrix(rs_grm, nrow=1)
  }
  
  rs <- list('3pl'=rs_3pl, 'gpcm'=rs_gpcm, 'grm'=rs_grm)
  rs <- rs[match(names(items), names(rs))]
  if(combine)
    rs <- Reduce(cbind, rs)
  rs
}

#' @rdname model_mixed
#' @param u the response data, a 2d matrix
#' @param log \code{TRUE} to return log-likelihood
#' @examples 
#' with(model_mixed_gendata(10, 4, 4, 4, 3), model_mixed_lh(u, t, items))
#' @export
model_mixed_lh <- function(u, t, items, D=1.702, log=FALSE, combine=TRUE) {
  rs_3pl <- rs_gpcm <- rs_grm <- NULL
  u_ix <- sapply(names(items), function(x) ifelse(is.null(items[[x]]), 0, nrow(items[[x]])))
  u_ix1 <- cumsum(u_ix)
  u_ix0 <- c(0, u_ix1[-length(u_ix)]) + 1
  u_ix1[u_ix == 0] <- 0
  u_ix0[u_ix == 0] <- 0
  u <- Map(function(x, y) u[, x:y], u_ix0, u_ix1)
  names(u) <- names(items)
  
  if(!is.null(items$'3pl'))
    rs_3pl <- with(items$'3pl', model_3pl_lh(u$'3pl', t, a, b, c, D=D, log=log))
  
  if(!is.null(items$'gpcm')) {
    d <- as.matrix(with(items, gpcm[, grep('^d[0-9]+$', colnames(gpcm), value=TRUE)]))
    rs_gpcm <- with(items, model_gpcm_lh(u$'gpcm', t, gpcm[, 'a'], gpcm[,'b'], d, D=D, log=log))
  }
  
  if(!is.null(items$'grm')) {
    b <- as.matrix(with(items, grm[, grep('^b[0-9]+$', colnames(grm), value=TRUE)]))
    rs_grm <- with(items, model_grm_lh(u$'grm', t, grm[, 'a'], b, D=D, log=log))
  }
  
  rs <- list('3pl'=rs_3pl, 'gpcm'=rs_gpcm, 'grm'=rs_grm)
  rs <- rs[match(names(items), names(rs))]
  if(combine)
    rs <- Reduce(cbind, rs)
  rs
}

