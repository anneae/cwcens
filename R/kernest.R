#' Estimate the probability in state and restricted mean time in an illness-death
#' model with component-wise censoring
#'
#'This function implements the
#'
#' @param dat a dataframe with one row per individual with the variables
#' \code{t1-tm} and \code{x1-xm} where \code{m}
#' is the number of visits and \code{ti} and \code{xi} are the time and status at
#' each visit, and the variables \code{dtime} and \code{dstatus} which are the time
#' and event indicator for death, and the variable \code{nvisits} is the number of visits.
#' \code{xi=1} represents alive and event free
#' (state 1 in a multistate illness-death model), and \code{xi=2} represents alive with event.
#' @param bandwidth specifies the bandwidth to be used for the kernel estimator.
#' This can be selected data-adaptively using the \code{dab()} function.
#' @param tau2 the maximum time that individuals were at risk for a visit,
#' based on the study design.
#' @param prob.times a vector of times at which probability in state will be estimated
#' @param mu.times a vector of restriction times at which restricted mean time
#' in state will be estimated
#' @param boundary specifies how kernel estimation is done in the left boundary region
#' from zero to the bandwidth. The default is \code{boundary.kernel}, meaning a boundary
#'  kernel is used in the left boundary region. Set \code{boundary = 'interpolation'}
#'  to use linear interpolation through the points (0,1) and (h, r) where h is the
#'  bandwidth and r is the kernel estimate at time h.
#' @param kfun specifies the kernel function to be used for estimation. The default
#' is \code{epanechnikov}; other possible values are \code{triweight} and \code{uniform}
#' @param std.err If \code{std.err= 'asymptotic'}  or \code{'boot'}, the function
#' calculates the standard error estimates and 95% confidence intervals for each
#' quantity using the asymptotic or bootstrap estimators. (\code{std.err= 'none'} is the default.)
#' @param B the number of bootstrap samples; the default value of 50 for the sake
#' of computation time, but we recommend increasing it
#' @param boot.seed If \code{boot.seed} is specified, \code{set.seed(boot.seed)}
#'  will be run before generating bootstrap samples, so the samples can be reproduced.
#' @param scale a scaling factor for the restricted mean time in state output. For example,
#' if times are in days and you want the output to reflect restricted mean years in state,
#' set \code{scale = 365.25}.
#'
#' @return A list with up to two elements, \code{prob.info}  if \code{prob.times} was non-null,
#' and \code{mu.info} if \code{mu.times} was non-null. \code{prob.info} contains
#'  probability in state estimates and \code{mu.info} contains restricted mean time
#'  in state estimates.
#'  The columns in \code{prob.info} are \code{t, p1, p2, p3} for time and
#'  probability in state 1, 2 and 3, respectively. If \code{std.err ='boot'} or
#'  \code{'asymptotic'}, additional columns are added with standard error estimates
#'  and lower and upper limits of the 95% confidence interval for each estimate.
#'  The columns in \code{mu.info} are analogous.
#' @export
#'
#' @examples mydat <- simdat(50, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
#' vital.lfu=c(30.4*36, 30.4*48),
#' visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48), scatter.sd=10)
#' kernel.est(mydat, bandwidth=12*30.4, tau2=30.4*48, prob.times=30.4*48, mu.times=30.4*48,
#' boundary = 'boundary.kernel', kfun='epanechnikov',
#' std.err='none',  scale=12*30.4)
kernel.est <- function(dat, bandwidth, tau2, prob.times=NULL, mu.times=NULL,
                       boundary = 'boundary.kernel', kfun='epanechnikov',
                       std.err='none', B=50, boot.seed = NA,
                       scale=1){

  utils::globalVariables(c("Kq", "K1"))

  if (!is.null(prob.times)) if (all.equal(sort(prob.times),prob.times)!=T) stop('prob.times must be in ascending order.')
  if (!is.null(mu.times)) if (all.equal(sort(mu.times),mu.times)!=T) stop('mu.times must be in ascending order.')
  if (!(boundary %in% c('boundary.kernel','interpolation'))) stop('The only allowable entries for boundary are "boundary.kernel" and "interpolation".')

  if (kfun == 'epanechnikov') {
    Kq <<- function(x, q) {
      sigk1 <- sqrt(0.2)
      2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
        (1 + ((q - 1) / (q + 1) / sigk1) ^ 2 + 2 / sigk1 ^ 2 * (1 - q) / (1 + q) ^ 2 * x)
      # Could also use  6*(1+x)*(q-x)*((1+q)^(-3))*(1+5*((1-q)/(1+q))^2+10*(1-q)*((1+q)^(-2))*x)
      # Equivalent when x<=q which is the only time we'd use it!
    }
    K1 <<- function(u) {
      0.75 * (1 - u ^ 2) * (abs(u) < 1)
    }
  }
  else if (kfun == 'uniform') {
    Kq <<- function(x, q) {
      (1+q)^(-1)*(1+3*((1-q)/(1+q))^2+6*(1-q)*(1+q)^(-2)*x)
    }
    K1 <<- function(u) 0.5* (abs(u) < 1)
  }
  else if (kfun == 'triweight'){
    Kq <<- function(x, q) {
      30*(1+x)^2*(q-x)^2*(1+q)^(-5)*(1+7*((1-q)/(1+q))^2+14*(1-q)*((1+q)^(-2))*x)
    }
    K1 <<- function(u) (15/16)*(1-u^2)^2 * (abs(u) < 1)
  }

  # Prep data
  nvisits<-dat$nvisits[1]
  ID<- 1:nrow(dat)
  X <- dat$dtime
  E <- dat$dstatus
  V.all <- c(t(as.matrix(dat[,paste('t',1:nvisits, sep='')])))
  Y.all <- as.numeric(c(t(as.matrix(dat[,paste('x',1:nvisits, sep='')])))==1)
  ID.all <- rep(1:nrow(dat), each = nvisits)
  missing <- is.na(V.all)
  V.all <- V.all[!missing]
  Y.all <- Y.all[!missing]
  ID.all <- ID.all[!missing]
  N<-length(X)

  #####################################################
  # generate survival data
  #####################################################
  fit <- survival::survfit(survival::Surv(X, E)~1)
  Left <- c(0,summary(fit)$time)
  Right <- c(summary(fit)$time,Inf)
  Value <- c(1,summary(fit)$surv)

  #####################################################
  # calculate values at the times of interest
  #####################################################
  if (length(prob.times)>0){
    prob.matrix<-matrix(NA, nrow = length(prob.times), ncol = 3)
    prob.matrix[,1]<-f2(prob.times, v.all = V.all, y.all = Y.all, h = bandwidth, N = N, tau2 = tau2, boundary = boundary,
                        Left=Left, Right= Right, Value=Value)
    prob.matrix[,3]<-1-sapply(1:length(prob.times), function (x)
      Sd(prob.times[x], Value, Left, Right))
    prob.matrix[,2]<-1-prob.matrix[,1]-prob.matrix[,3]
  }

  if (length(mu.times)>0){
    rm.matrix<-matrix(NA, nrow = length(mu.times), ncol = 3)
    mu.times2<-c(0, mu.times)
    rm.diff<-sapply(1:length(mu.times), function(x) int_f2(mu.times2[x], mu.times2[x+1],bandwidth, V.all,Y.all,N,tau2, boundary, Left, Right, Value))
    rm.matrix[,1]<-cumsum(rm.diff)/scale
    rm.matrix[,3]<-mu.times/scale-sapply(1:length(mu.times), function(i)
      summary(fit, rmean =mu.times[i], scale =scale, extend = T)$table[5])
    rm.matrix[,2]<-mu.times/scale-rm.matrix[,1]-rm.matrix[,3]
  }

  if (std.err=='none')  {
    if (length(prob.times)>0){
      prob.matrix<-cbind(prob.times/scale, prob.matrix)
      colnames(prob.matrix)<-c('time','p1','p2','p3')
    }
    if (length(mu.times)>0){
      rm.matrix<-cbind(mu.times/scale, rm.matrix)
      colnames(rm.matrix)<-c('time','mu1','mu2','mu3')
    }
  }
  else if (std.err == 'asymptotic'){
    #####################################################
    # generate longitudinal data
    #####################################################
    xi.all <- xi(V.all,bandwidth, V.all, Y.all, N, tau2 = tau2)
    xi.all2 <- xi(V.all,bandwidth, V.all, 1-Y.all, N, tau2 = tau2)
    lambda.all <- lambda(V.all,bandwidth, V.all, N, tau2 = tau2)
    sd.all <- xi.all
    for(i in 1:length(sd.all))
    {
      sd.all[i] <- Sd(V.all[i], Value, Left, Right)
    }
    #####################################################
    # calculate mu, S_X at terminal time
    #####################################################
    sx2.all <- rep(0,N)
    for(i in 1:N)
    {
      sx2.all[i] <- sum(X>=X[i])/N
    }
    mu2.all <- rep(0,N)
    mu2.2.all <- rep(0,N)
    for(i in 1:N)
    {
      tt <- X[i]
      mu2.all[i] <- sum(sd.all[V.all<tt]*Y.all[V.all<tt]/lambda.all[V.all<tt])/N
      mu2.2.all[i] <- sum(sd.all[V.all<tt]*(1-Y.all)[V.all<tt]/lambda.all[V.all<tt])/N
    }
    mu2.3.all<-mu.times-sapply(1:length(X), function(i)
      summary(fit, rmean =X[i], extend = T)$table[5])
    # The mu2s are NOT scaled!

    #####################################################
    # calculate Phi_i
    #####################################################
    Phi1.1<-Phi2.1<-Phi3.1<-matrix(NA, nrow = N, ncol = length(mu.times))
    Phi1.2<-Phi2.2<-Phi3.2<-matrix(NA, nrow = N, ncol = length(mu.times))
    Phi3<-matrix(NA, nrow = N, ncol = length(mu.times))
    phi3<-matrix(NA, nrow = N, ncol = length(prob.times))
    for(i in 1:N)
    {
      V.i <- V.all[ID.all == i]
      Y.i <- Y.all[ID.all == i]
      sd.i <- sd.all[ID.all == i]
      lambda.i <- lambda.all[ID.all == i]
      xi.i <- xi.all[ID.all == i]
      xi.i2<-xi.all2[ID.all == i]

      if (length(mu.times)>0){
        for (j in 1:length(mu.times)){
          # Variance of RMTIS1
          Phi1.1[i,j] <- (sum(1/(sx2.all[X<X[i] & X<mu.times[j] & E == TRUE])^2)/N - (1/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1))*rm.matrix[j,1]*scale +
            (mu2.all[i]/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1) - sum((mu2.all/(sx2.all)^2)[X<X[i] & X<mu.times[j] & E == TRUE])/N
          Phi2.1[i,j] <- sum((sd.i*Y.i/lambda.i)[V.i<mu.times[j]]) #- out
          Phi3.1[i,j] <- sum((sd.i*xi.i/lambda.i^2)[V.i<mu.times[j]]) #- out

          # Variance of RMTIS2
          Phi1.2[i,j] <- (sum(1/(sx2.all[X<X[i] & X<mu.times[j] & E == TRUE])^2)/N - (1/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1))*rm.matrix[j,2]*scale +
            (mu2.2.all[i]/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1) - sum((mu2.2.all/(sx2.all)^2)[X<X[i] & X<mu.times[j] & E == TRUE])/N
          Phi2.2[i,j] <- sum((sd.i*(1-Y.i)/lambda.i)[V.i<mu.times[j]]) #- out
          Phi3.2[i,j] <- sum((sd.i*(xi.i2)/lambda.i^2)[V.i<mu.times[j]]) #- out

          #Variance of RMTIS3
          Phi3[i,j] <- (sum(1/(sx2.all[X<X[i] & X<mu.times[j] & E == TRUE])^2)/N - (1/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1))*rm.matrix[j,3]*scale +
            (mu2.3.all[i]/sx2.all[i])*(X[i]<mu.times[j] & E[i] == 1) - sum((mu2.3.all/(sx2.all)^2)[X<X[i] & X<mu.times[j] & E == TRUE])/N

        }
      }

      if (length(prob.times)>0){
        # variance of sqrt(n)*(probability in state 3)
        for (j in 1:length(prob.times)){
          phi3[i,j]<- -Sd(prob.times[j], Value, Left, Right)*(sum(1/(sx2.all[X<X[i] & X<prob.times[j] & E == TRUE])^2)/N - (1/sx2.all[i])*(X[i]<prob.times[j] & E[i] == 1))
        }
      }
    }

    Phi <-  (Phi1.1 + Phi2.1 - Phi3.1)/scale
    Phi.2 <-  (Phi1.2 + Phi2.2 - Phi3.2)/scale
    Phi3 <- Phi3/scale

    if (length(prob.times)>0){
      # variance of sqrt(n)*(prob in state 1). Same as state 2.
      var1 <- (3/5*((1-prob.matrix[,3])*prob.matrix[,1]-prob.matrix[,1]^2)/lambda(prob.times,bandwidth, V.all, N, tau2 = tau2))/bandwidth

      prob.matrix<-cbind(prob.matrix, sqrt(var1)/sqrt(N), sqrt(var1)/sqrt(N), apply(phi3, 2, sd)/sqrt(N))
      prob.matrix<-cbind(prob.times/scale, prob.matrix, prob.matrix[,1:3]-1.96*prob.matrix[,4:6], prob.matrix[,1:3]+1.96*prob.matrix[,4:6])
      colnames(prob.matrix)<-c('time','p1','p2','p3','p1se', 'p2se','p3se',
                               'p1lower','p2lower','p3lower','p1upper','p2upper','p3upper')
    }
    if (length(mu.times)>0){
      rm.matrix<-cbind(rm.matrix, apply(Phi, 2, sd)/sqrt(N), apply(Phi.2, 2, sd)/sqrt(N), apply(Phi3, 2, sd)/sqrt(N))
      rm.matrix<-cbind(mu.times/scale, rm.matrix, rm.matrix[,1:3]-1.96*rm.matrix[,4:6], rm.matrix[,1:3]+1.96*rm.matrix[,4:6])
      colnames(rm.matrix)<-c('time','mu1','mu2','mu3','mu1se', 'mu2se','mu3se',
                             'mu1lower','mu2lower','mu3lower','mu1upper','mu2upper','mu3upper')
    }
  }

  else if (std.err == 'boot'){
    if (length(prob.times)>0) prob.boot<-matrix(NA, nrow = B, ncol = 3*length(prob.times))
    if (length(mu.times)>0) rm.boot<-matrix(NA, nrow = B, ncol = 3*length(mu.times))
    if(!is.na(boot.seed)) set.seed(boot.seed)
    indices<-matrix(sample(nrow(dat), nrow(dat)*B, replace = T), nrow = B)

    for (i in 1:B) {
      dat.bs <- dat[indices[i,],]
      # Prep data
      nvisits<-dat.bs$nvisits[1]
      ID<- 1:nrow(dat.bs)
      X <- dat.bs$dtime
      E <- dat.bs$dstatus
      V.all <- c(t(as.matrix(dat.bs[,paste('t',1:nvisits, sep='')])))
      Y.all <- as.numeric(c(t(as.matrix(dat.bs[,paste('x',1:nvisits, sep='')])))==1)
      ID.all <- rep(1:nrow(dat.bs), each = nvisits)
      missing <- is.na(V.all)
      V.all <- V.all[!missing]
      Y.all <- Y.all[!missing]
      ID.all <- ID.all[!missing]
      N<-length(X)

      #####################################################
      # generate survival data
      #####################################################
      fit <- survival::survfit(survival::Surv(X, E)~1)
      Left <- c(0,summary(fit)$time)
      Right <- c(summary(fit)$time,Inf)
      Value <- c(1,summary(fit)$surv)

      #####################################################
      # calculate values at the times of interest
      #####################################################
      if (length(prob.times)>0){
        prob.matrix.boot<-matrix(NA, nrow = length(prob.times), ncol = 3)
        prob.matrix.boot[,1]<-f2(prob.times, v.all = V.all, y.all = Y.all, h = bandwidth, N = N, tau2 = tau2, boundary = boundary,
                            Left=Left, Right= Right, Value=Value)
        prob.matrix.boot[,3]<-1-sapply(1:length(prob.times), function (x)
          Sd(prob.times[x], Value, Left, Right))
        prob.matrix.boot[,2]<-1-prob.matrix.boot[,1]-prob.matrix.boot[,3]
        prob.boot[i,]<-t(prob.matrix.boot)
      }
      if (length(mu.times)>0){
        rm.matrix.boot<-matrix(NA, nrow = length(mu.times), ncol = 3)
        mu.times2<-c(0, mu.times)
        rm.diff<-sapply(1:length(mu.times), function(x) int_f2(mu.times2[x], mu.times2[x+1],bandwidth, V.all,Y.all,N,tau2, boundary, Left, Right, Value, warn = F))
        rm.matrix.boot[,1]<-cumsum(rm.diff)/scale
        rm.matrix.boot[,3]<-mu.times/scale-sapply(1:length(mu.times), function(i)
          summary(fit, rmean =mu.times[i], scale =scale, extend = T)$table[5])
        rm.matrix.boot[,2]<-mu.times/scale-rm.matrix.boot[,1]-rm.matrix.boot[,3]
        rm.boot[i,]<-t(rm.matrix.boot)
      }
    }
    if (length(prob.times)>0) {
      prob.se<-cbind(matrix(apply(prob.boot,2,sd), ncol = 3, byrow = T),
                matrix(apply(prob.boot,2,function(x) quantile(x,probs = 0.025,na.rm = T)), ncol=3, byrow = T),
                matrix(apply(prob.boot,2,function(x) quantile(x,probs = 0.975,na.rm = T)), ncol=3, byrow = T))
      prob.matrix<-cbind(prob.times/scale, prob.matrix, prob.se)
      colnames(prob.matrix)<-c('time','p1','p2','p3','p1se', 'p2se','p3se',
                               'p1lower','p2lower','p3lower','p1upper','p2upper','p3upper')
    }
    if (length(mu.times)>0) {
      rm.se<-cbind(matrix(apply(rm.boot,2,sd), ncol = 3, byrow = T),
                matrix(apply(rm.boot,2,function(x) quantile(x,probs = 0.025,na.rm = T)), ncol=3, byrow = T),
                matrix(apply(rm.boot,2,function(x) quantile(x,probs = 0.975,na.rm = T)), ncol=3, byrow = T))
      rm.matrix<-cbind(mu.times/scale, rm.matrix, rm.se)
      colnames(rm.matrix)<-c('time','mu1','mu2','mu3','mu1se', 'mu2se','mu3se',
                           'mu1lower','mu2lower','mu3lower','mu1upper','mu2upper','mu3upper')
    }
  }
  if (!is.null(prob.times)&!is.null(mu.times)) return(list(prob.info = prob.matrix, mu.info =rm.matrix))
  else if(!is.null(prob.times)) return(list(prob.info = prob.matrix))
  else if(!is.null(mu.times)) return(list(mu.info = rm.matrix))
}
