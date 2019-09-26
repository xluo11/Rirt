#' Helper Functions
#' @name helpers
NULL

#' @rdname helpers
#' @description \code{model_polytomous_3dindex} creates indices extracting 3D stats
#' @param u the observed response, 2d matrix, values start from 0
#' @keywords internal
model_polytomous_3dindex <- function(u){
  n_p <- dim(u)[1]
  n_i <- dim(u)[2]
  n_c <- max(u) + 1
  cbind(rep(1:n_p, n_i), rep(1:n_i, each=n_p), as.vector(u+1))
}

#' @rdname helpers
#' @description \code{model_polytomous_3dresponse} converts responses from 2D to 3D
#' @keywords internal
model_polytomous_3dresponse <- function(u){
  n_c <- max(u) + 1
  x <- array(0, dim=c(dim(u), n_c))
  x[model_polytomous_3dindex(u)] <- 1
  x
}

#' @rdname helpers
#' @description \code{hermite_gauss} stores pre-computed hermite gaussian 
#' quadratures points and weights
#' @param degree the degree of hermite-gauss quadrature: '20', '11', '7'
#' @keywords internal
hermite_gauss <- function(degree=c('20', '11', '7')){
  switch(match.arg(degree),
         '20'=list(t=c(-5.38748089001123,-4.60368244955074,-3.94476404011562,-3.34785456738321,-2.78880605842813,-2.25497400208927,-1.73853771211658,-1.23407621539532,-0.737473728545394,-0.245340708300901,0.245340708300901,0.737473728545394,1.23407621539532,1.73853771211658,2.25497400208927,2.78880605842813,3.34785456738321,3.94476404011562,4.60368244955074,5.38748089001123),
                   w=c(2.22939364553415E-13,4.39934099227318E-10,1.08606937076928E-07,7.80255647853206E-06,0.000228338636016353,0.00324377334223786,0.0248105208874636,0.109017206020023,0.286675505362834,0.46224366960061,0.46224366960061,0.286675505362834,0.109017206020023,0.0248105208874636,0.00324377334223786,0.000228338636016353,7.80255647853206E-06,1.08606937076928E-07,4.39934099227318E-10,2.22939364553415E-13)),
         '11'=list(t=c(-3.66847084655958,-2.78329009978165,-2.02594801582575,-1.32655708449493,-0.656809566882099,0,0.656809566882099,1.32655708449493,2.02594801582575,2.78329009978165,3.66847084655958), 
                   w=c(1.43956039371425E-06,0.000346819466323345,0.0119113954449115,0.117227875167708,0.429359752356125,0.654759286914591,0.429359752356125,0.117227875167708,0.0119113954449115,0.000346819466323345,1.43956039371425E-06)),
         '7'=list(t=c(-2.651961356835233492447,-1.673551628767471445032,-0.8162878828589646630387,0,0.8162878828589646630387,1.673551628767471445032,2.651961356835233492447),
                  w=c(9.71781245099519154149E-4,0.05451558281912703059218,0.4256072526101278005203,0.810264617556807326765,0.4256072526101278005203,0.0545155828191270305922,9.71781245099519154149E-4)))
}

#' @rdname helpers
#' @description \code{nr_iteration} updates the parameters using the Newton-Raphson method
#' @param param the parameter being estimated
#' @param free TRUE to free parameters, otherwise fix parameters
#' @param dv the first and second derivatives
#' @param h_max the maximum value of h
#' @param lr the learning rate
#' @param bound the lower and upper bounds of the parameter
#' @keywords internal
nr_iteration <- function(param, free, dv, h_max, lr, bound){
  h <- dv$dv1 / dv$dv2
  h[is.na(h)] <- 0
  h <- ifelse(abs(h) > h_max, sign(h) * h_max, h) * lr
  h[!free] <- 0
  param <- param - h
  param[param < bound[1]] <- bound[1]
  param[param > bound[2]] <- bound[2]
  list(param=param, h=h)
}


#' @rdname helpers
#' @param tracking estimation tracking information
#' @param k the number of iterations in estimation
#' @importFrom reshape2 melt
#' @import ggplot2
#' @keywords internal
estimate_3pl_debug <- function(tracking, k){
  xx <- with(tracking, data.frame(fit=fit, t=t, a=a, b=b, c=c))[1:k, ]
  xx$iteration <- 1:k
  xx <- melt(xx, id.vars='iteration')
  xx <- xx[!is.na(xx$value), ]
  g <- ggplot(xx, aes_string(x="iteration", y="value", color="variable")) + 
    geom_line() + facet_wrap(~variable, scales="free") + guides(color=F) + 
    xlab('Iterations') + ylab('') + theme_bw()
  print(g)
  invisible(g)
}

#' @rdname helpers
#' @param true_params a list of true parameters
#' @param t estimated ability parameters
#' @param a estimated discrimination parameters
#' @param b estimated difficulty parameters
#' @param c estimated guessing parameters
#' @param t_free TRUE to estimate ability parameters, otherwise fix
#' @param a_free TRUE to estimate discrimination parameters, otherwise fix
#' @param b_free TRUE to estimate difficulty parameters, otherwise fix
#' @param c_free TRUE to estimate guessing parameters, otherwise fix
#' @import ggplot2
#' @keywords internal
estimate_3pl_eval <- function(true_params, t, a, b, c, t_free, a_free, b_free, c_free) {
  xx <- rbind(data.frame(true=true_params$t, est=t, params='t'),
              data.frame(true=true_params$a, est=a, params='a'),
              data.frame(true=true_params$b, est=b, params='b'),
              data.frame(true=true_params$c, est=c, params='c'))
  g <- ggplot(xx, aes_string(x="true", y="est", color="params")) + 
    geom_point(alpha=.3) + geom_smooth(method='gam', se=F) + 
    facet_wrap(~params, nrow=1, scales='free') + guides(color=F) +
    xlab('True Parameters') + ylab('Est. Parameters') + theme_bw()
  print(g)
  if(any(t_free)) 
    cat('t: corr = ', round(cor(t, true_params$t), 3), ', rmse = ', round(rmse(t, true_params$t), 3),'\n', sep='')
  if(any(a_free)) 
    cat('a: corr = ', round(cor(a, true_params$a), 3), ', rmse = ', round(rmse(a, true_params$a), 3),'\n', sep='')
  if(any(b_free)) 
    cat('b: corr = ', round(cor(b, true_params$b), 3), ', rmse = ', round(rmse(b, true_params$b), 3),'\n', sep='')
  if(any(c_free)) 
    cat('c: corr = ', round(cor(c, true_params$c), 3), ', rmse = ', round(rmse(c, true_params$c), 3),'\n', sep='')
}


#' @rdname helpers
#' @keywords internal
estimate_gpcm_debug <- function(tracking, k) {
  xx <- with(tracking, data.frame(fit=fit, t=t, a=a, b=b, d=d))[1:k, ]
  xx$iteration <- 1:k
  xx <- melt(xx, id.vars='iteration')
  xx <- xx[!is.na(xx$value),]
  g <- ggplot(xx, aes_string(x="iteration", y="value", color="variable")) + 
    geom_line() + facet_wrap(~variable, scales="free") + guides(color=F) + 
    xlab('Iterations') + ylab('') + theme_bw()
  print(g)
  invisible(g)
}

#' @rdname helpers
#' @keywords internal
estimate_gpcm_eval <- function(true_params, n_c, t, a, b, d, t_free, a_free, b_free, d_free) {
  xx <- rbind(data.frame(true=true_params$t, est=t, params='t'),
              data.frame(true=true_params$a, est=a, params='a'),
              data.frame(true=true_params$b, est=b, params='b'))
  for(i in 2:n_c)
    xx <- rbind(xx, data.frame(true=true_params$d[,i], est=d[,i], params=paste('d',i,sep='')))
  g <- ggplot(xx, aes_string(x="true", y="est", color="params")) + 
    geom_point(alpha=.3) + geom_smooth(method='gam', se=F) + 
    facet_wrap(~params, nrow=1, scales='free') + guides(color=F) +
    xlab('True Parameters') + ylab('Est. Parameters') + theme_bw()
  print(g)
  if(any(t_free)) 
    cat('t: corr = ', round(cor(t, true_params$t), 3), ', rmse = ', round(rmse(t, true_params$t), 3),'\n', sep='')
  if(any(a_free)) 
    cat('a: corr = ', round(cor(a, true_params$a), 3), ', rmse = ', round(rmse(a, true_params$a), 3),'\n', sep='')
  if(any(b_free)) 
    cat('b: corr = ', round(cor(b, true_params$b), 3), ', rmse = ', round(rmse(b, true_params$b), 3),'\n', sep='')
  for(i in 2:n_c) 
    if(any(d_free[,i])) 
      cat('d_', i, ': corr = ', round(cor(d[,i], true_params$d[,i]), 3), ', rmse = ', round(rmse(d[,i], true_params$d[,i]), 3),'\n', sep='')
}


#' @rdname helpers
#' @keywords internal
estimate_grm_debug <- function(tracking, k) {
  xx <- with(tracking, data.frame(fit=fit, t=t, a=a, b=b))[1:k,]
  xx$iteration <- 1:k
  xx <- melt(xx, id.vars='iteration')
  xx <- xx[!is.na(xx$value), ]
  g <- ggplot(xx, aes_string(x="iteration", y="value", color="variable")) + 
    geom_line() + facet_wrap(~variable, scales="free") + guides(color=F) + 
    xlab('Iterations') + ylab('') + theme_bw()
  print(g)
  invisible(g)
}

#' @rdname helpers
#' @keywords internal
estimate_grm_eval <- function(true_params, n_c, t, a, b, t_free, a_free, b_free) {
  xx <- rbind(data.frame(true=true_params$t, est=t, params='t'),
              data.frame(true=true_params$a, est=a, params='a'))
  for(i in 1:(n_c-1))
    xx <- rbind(xx, data.frame(true=true_params$b[,i], est=b[,i], params=paste('b',i,sep='')))
  g <- ggplot(xx, aes_string(x="true", y="est", color="params")) + 
    geom_point(alpha=.3) + geom_smooth(method='gam', se=F) + 
    facet_wrap(~params, nrow=1, scales='free') + guides(color=F) +
    xlab('True Parameters') + ylab('Est. Parameters') + theme_bw()
  print(g)
  if(any(t_free)) 
    cat('t: corr = ', round(cor(t, true_params$t), 3), ', rmse = ', round(rmse(t, true_params$t), 3),'\n', sep='')
  if(any(a_free)) 
    cat('a: corr = ', round(cor(a, true_params$a), 3), ', rmse = ', round(rmse(a, true_params$a), 3),'\n', sep='')
  for(i in 1:(n_c-1)) 
    if(any(b_free[,i])) 
      cat('b_', i, ': corr = ', round(cor(b[,i], true_params$b[,i]), 3), ', rmse = ', round(rmse(b[,i], true_params$b[,i]), 3),'\n', sep='')
}

