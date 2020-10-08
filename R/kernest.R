# This function calculates the probability in each state at prob.times and the
# restricted mean time in each state at mu.times.

# The dataset specified in dat should have one row per person with the variables
# dtime, dstatus, t1-tm and x1-xm. After the last visit t_ and x_ are NA.
# x=1 means alive and event free, x=2 means alive with event

# bandwidth specifies the bandwidth to be used for the kernel estimator. This can
# be selected data-adaptively using the dab() function.
# tau2

# boundary specifies how kernel estimation is done in the left boundary region.
# The default is 'boundary.kernel', meaning a boundary kernel is used in the
# left boundary region, but if you set boundary = 'interpolation',
# linear interpolation through the points (0,1) and (h, \hat r_h(h))
# is used to estimate r(t) in the left boundary region.

# kfun can take values of 'epanechnikov','triweight' or 'uniform'

# If std.err= 'asymptotic'  or 'boot', it also calculates the standard errors and 95% CIs using
# the asymptotic or bootstrap estimators. (std.err= 'none' is the default.)
# B is the number of bootstrap samples; the default value of B is 50 for the sake
# of computation time, but we recommend increasing it.
# If boot.seed is specified, set.seed(boot.seed) will be run before generating bootstrap
# samples so the samples are reproducible.

# scale is a scaling factor for the restricted mean time in state output. For example.
# if times are in days and you want the output to reflect restricted mean years in state,
# set scale = 365.25.

### The output is a list with up to two elements, prob.info for probability in state estimates
### and mu.info for restricted mean time in state estimates.
### The columns in prob.info will be t, p1, p2, p3 for time and probability in the three states,
### and the columns in mu.info will be analogous.
### If std.err ='boot' or 'asymptotic', additional columns will contain the standard error
### and lower and upper limits of the 95% confidence interval for each estimate.

kernel.est <- function(dat, bandwidth, tau2, prob.times=NULL, mu.times=NULL,
                       boundary = 'boundary.kernel', kfun='epanechnikov',
                       std.err='none', B=50, boot.seed = NA,
                       scale=1){

  if (!is.null(prob.times)) if (all.equal(sort(prob.times),prob.times)!=T) stop('prob.times must be in ascending order.')
  if (!is.null(mu.times)) if (all.equal(sort(mu.times),mu.times)!=T) stop('mu.times must be in ascending order.')
  if (!(boundary %in% c('boundary.kernel','interpolation'))) stop('The only allowable entries for boundary are "boundary.kernel" and "interpolation".')

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
  fit <- survfit(Surv(X, E)~1)
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
      fit <- survfit(Surv(X, E)~1)
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


################ ################ ################
################ EXAMPLES
################ ################ ################
#source('~/Dropbox/Ann Eaton PhD research not shared/RMST/computing code/cwcens/R/internal_fns.R')
#source('~/Dropbox/Ann Eaton PhD research not shared/RMST/computing code/cwcens/R/simdat.R')
#library(msm); library(pracma); library(survival)

# Plot estimates from two visit processes
#dat<-simdat(55, visit.schedule = 30.4*c(12,24,36,48,60), vital.lfu = c(30.4*66, 30.4*78))
#test<-kernel.est(dat, 12*30.4, NULL, NULL, 1:1200, T, std.err = 'formula',
#                 tau2 = 78*30.4, scale = 12*30.4, boundary = 'interpolation')
#plot(1:1200, test[1,1:1200])
#dat<-simdat(55, visit.schedule = NA, visitrate = function(x) sapply(x, function(t) 1/365.25),
#            vital.lfu = c(30.4*66, 30.4*78))
#test<-kernel.est(dat, 12*30.4, NULL, NULL, 1:1200, T, std.err = 'formula',
#                 tau2 = 78*30.4, scale = 12*30.4, boundary = 'interpolation')
#plot(1:1200, test[1,1:1200])

# Estimate RMST with boundary kernel and formula
#test<-kernel.est(dat, 12*30.4, c(1, 12*30.4, 60*30.4), 12*30.4, 1:5, T, std.err = 'formula',
#                 tau2 = 78*30.4, scale = 12*30.4)
# Estimate RMST with interpolation and bootstrap
#test<-kernel.est(dat, 12*30.4, c(1, 12*30.4, 60*30.4), 12*30.4, 1:5, T, boundary = 'interpolation',
#                 std.err = 'boot',
#                 tau2 = 78*30.4, scale = 12*30.4, B=500)

# normal visit results - compare to old results
# to check whether boundary kernel gives reasonable results
# and whether formula based results are reasonable
#start<-Sys.time()
#set.seed(22)
#M<-500
#simres.formula<-matrix(NA, ncol = 6*4, nrow = M)
#for (m in 1:M){
#    dat<-simdat(m, visit.schedule = 30.4*c(12,24,36,48,60), vital.lfu = c(30.4*66, 30.4*78))
#    simres.formula[m,]<-t(kernel.est(dat, 12*30.4, 5*12*30.4, 5*12*30.4, NULL, T, std.err = 'formula',
#                     tau2 = 78*30.4, scale = 12*30.4))
#}
#print(Sys.time()-start)
## Time difference of 24.32166 mins

#start<-Sys.time()
#set.seed(22)
#M<-3
#simres<-matrix(NA, ncol = 6*4, nrow = M)
#for (m in 1:M){
#    dat<-simdat(sample(-2^25:2^25,1), visit.schedule = 30.4*c(12,24,36,48,60), vital.lfu = c(30.4*66, 30.4*78))
#    simres[m,]<-t(kernel.est(dat, 12*30.4, 5*12*30.4, 5*12*30.4, NULL, T, std.err = 'boot',
#                             tau2 = 78*30.4, scale = 12*30.4, B=500))
#}
#print(Sys.time()-start)
