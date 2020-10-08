# This function calculates the probability in each state at prob.times and the
# restricted mean time in each state at mu.times.
# It also calculates either the formula based standard errors and 95% CIs,
# or the bootstrap versions. (std.err= 'boot' or 'formula')

### TO DO! replace SE with
# CI which can take values of 'none', 'boot' or 'asymptotic'

# It uses the boundary kernel on the left, but if you set boundary = 'interpolation'
# it will linearly interpolate from the point (0,1) to (h, \hat r_h(h)).

# The dataset dat should have one row per person with the variables
# dtime, dstatus, t1-tm and x1-xm. After the last visit t_ and x_ are NA.
# x=1 means alive and event free, x=2 means alive with event

### UPDATE: THE OUTPUT IS GOING TO BE A LIST WITH UP TO 2 ELEMENTS
### result$prob.info and result$mu.info
### The columns in prob.info will be t, p1, p2, p3,
### And if CI ='boot' or 'asymptotic':
### p1.se, p2.se, p3.se, p1.lower, p2.lower, p3.lower, p1.upper, p2.upper, p3.upper
### And analogous columns for mu.info
### OLD:
### The output consists of a matrix with 4 rows and
### 3*length(prob.times) + 3*length(mu.times) + 3*length(curve.times) columns.
### The order of the elements in the first row is, for the first
### 3*length(prob.times) + 3*length(mu.times) estimates,
### probability or time in state 1,2,3 at each time, then move to the next time.
### Then for the last 3*length(curve.times), it's probability in state 1 at
### all times, followed by probability in state 2 at all times, followed by
### probability in state 3 at all times.
### The second row contains SE estimate for each element in the first row,
### except SE is not estimated for curve.times, the 3rd and 4th rows contain
### lower and upper limits of the 95% confidence interval.
### If SE = F, the elements in the second, third and fourth row will all be NA.

### If incl.dis = T, the output is a LIST with two elements, res and discon.
### discon is a vector that gives the integral excluded from the integration
### If no interval was excluded, it's c(0,0).

kernel.est <- function(dat, bandwidth,  tau2,
                       prob.times=NULL, mu.times=NULL,
                       boundary = 'boundary.kernel', kfun='epanechnikov',
                       SE=F, std.err='boot', B=500, boot.seed = NA,
                       scale=1, incl.discon=F){
  if (!is.null(prob.times)) if (all.equal(sort(prob.times),prob.times)!=T) stop('prob.times must be in ascending order.')
  if (!is.null(mu.times)) if (all.equal(sort(mu.times),mu.times)!=T) stop('mu.times must be in ascending order.')
  if (!(boundary %in% c('boundary.kernel','interpolation'))) stop('The only allowable entries for boundary are "boundary.kernel" and "interpolation".')

  if(is.null(mu.times) & incl.discon==T) stop('No integration is done; incl.discon must be F.')

  # Prep data for Yifei's method
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

  #plot(curve.times,xi(curve.times, bandwidth,V.all, Y.all, N, tau2))
  #plot(curve.times,lambda(curve.times, bandwidth,V.all, N, tau2))
  #hist(V.all,breaks = 20)
  #####################################################
  # calculate values at the times of interest
  #####################################################
  rm.matrix<-NULL
  if (length(mu.times)>0){
    if (boundary=='boundary.kernel' & lambda(0, bandwidth, V.all, N,tau2)<0) {
      dis<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
                   interval = c(0,bandwidth))$root
    }
    else dis<-0
    discon<-c(0,0)
    rm.matrix<-matrix(NA, nrow = length(mu.times), ncol = 3)
    mu.times2<-c(0, mu.times)
    rm.diff<-rep(NA, length(mu.times))
    for (j in 1:length(mu.times)){
      if (dis<=mu.times2[j] | dis>=mu.times2[j+1] ) rm.diff[j]<-integral(f2, xmin = mu.times2[j], xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                                                         tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                                                         reltol = 1e-5)
      else {
        part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                            tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                            reltol = 1e-5))
        discon[1]<-dis
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-5, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1e-5
        }
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-4, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1e-4
        }
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-3, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1e-3
        }
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-2, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1e-2
        }
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-1, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1e-1
        }
        if (inherits(part1, "try-error")){
          part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[1]<-dis-1
        }

        part2<-try(integral(f2, xmin = dis, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                            tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                            reltol = 1e-5))
        discon[2]<-dis
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1e-5, xmax =mu.times2[j+1] , v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1e-5
        }
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1e-4, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1e-4
        }
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1e-3, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1e-3
        }
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1e-2, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1e-2
        }
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1e-1, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1e-1
        }
        if (inherits(part2, "try-error")){
          part2<-try(integral(f2, xmin = dis+1, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                              tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                              reltol = 1e-5))
          discon[2]<-dis+1
        }
        rm.diff[j]<-part1+part2
      }
    }
    if (discon[2]-discon[1]>mu.times[1]/1000) warning('The interval being excluded from the integration has width>(earliest RM time)/1000')
    discon.to.ret <-  discon

    #   discon<-c(0,0)
    #   if (boundary=='boundary.kernel' & xi(0, bandwidth, V.all, Y.all, N,tau2)<0) {
    #       discon[1]<- uniroot(function (x) xi(x, bandwidth, V.all, Y.all, N,tau2),
    #                        interval = c(0,bandwidth))$root
    #       discon[2]<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
    #                        interval = c(0,bandwidth))$root
    #       discon<-sort(discon)
    #       if (discon[2]-discon[1]>mu.times[1]/200) warning('The interval being excluded from the integration has width>(earliest RM time)/200')
    #       discon.to.ret <-  discon
    #   }
    #   rm.matrix<-matrix(NA, nrow = length(mu.times), ncol = 3)
    #   mu.times2<-c(0, mu.times)
    #   rm.diff<-sapply(1:length(mu.times), function(j)
    #       if (discon[2]<=mu.times2[j] | discon[1]>=mu.times2[j+1] ) integral(f2, xmin = mu.times2[j], xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
    #                tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
    #                reltol = 1e-5)
    #       else integral(f2, xmin = mu.times2[j], xmax = discon[1]-1e-4, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
    #                     tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
    #                     reltol = 1e-5)+
    #            integral(f2, xmin = discon[2]+1e-4, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
    #             tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
    #             reltol = 1e-5))
    rm.matrix[,1]<-cumsum(rm.diff)/scale
    rm.matrix[,3]<-mu.times/scale-sapply(1:length(mu.times), function(i)
      summary(fit, rmean =mu.times[i], scale =scale, extend = T)$table[5])
    rm.matrix[,2]<-mu.times/scale-rm.matrix[,1]-rm.matrix[,3]
  }
  prob.matrix<-NULL
  if (length(prob.times)>0){
    prob.matrix<-matrix(NA, nrow = length(prob.times), ncol = 3)
    prob.matrix[,1]<-f2(prob.times, v.all = V.all, y.all = Y.all, h = bandwidth, N = N, tau2 = tau2, boundary = boundary,
                        Left=Left, Right= Right, Value=Value)
    prob.matrix[,3]<-1-sapply(1:length(prob.times), function (x)
      Sd(prob.times[x], Value, Left, Right))
    prob.matrix[,2]<-1-prob.matrix[,1]-prob.matrix[,3]
  }

  if (SE == F)  {
    prob.matrix<-cbind(prob.times/scale, prob.matrix)
    rm.matrix<-cbind(mu.times/scale, rm.matrix)
    colnames(prob.matrix)<-c('time','p1','p2','p3')
    colnames(rm.matrix)<-c('time','mu1','mu2','mu3')
  }
  else if (SE == T & std.err == 'formula'){
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
    phi1<-phi2<-phi3<-matrix(NA, nrow = N, ncol = length(prob.times))
    var1<-var2<-rep(NA, length(prob.times))
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
        # variance of sqrt(n)*(prob in state 1). Same as state 2.
        var1 <- (3/5*((1-prob.matrix[,3])*prob.matrix[,1]-prob.matrix[,1]^2)/lambda(prob.times,bandwidth, V.all, N, tau2 = tau2))/bandwidth
        var2 <- var1
      }
    }
    Phi <-  (Phi1.1 + Phi2.1 - Phi3.1)/scale
    Phi.2 <-  (Phi1.2 + Phi2.2 - Phi3.2)/scale
    Phi3 <- Phi3/scale

    prob.matrix<-cbind(prob.matrix, sqrt(var1)/sqrt(N), sqrt(var2)/sqrt(N), apply(phi3, 2, sd)/sqrt(N))
    rm.matrix<-cbind(rm.matrix, apply(Phi, 2, sd)/sqrt(N), apply(Phi.2, 2, sd)/sqrt(N), apply(Phi3, 2, sd)/sqrt(N))

    prob.matrix<-cbind(prob.times/scale, prob.matrix, prob.matrix[,1:3]-1.96*prob.matrix[,4:6], prob.matrix[,1:3]+1.96*prob.matrix[,4:6])
    rm.matrix<-cbind(mu.times/scale, rm.matrix, rm.matrix[,1:3]-1.96*rm.matrix[,4:6], rm.matrix[,1:3]+1.96*rm.matrix[,4:6])
    colnames(prob.matrix)<-c('time','p1','p2','p3','p1se', 'p2se','p3se',
                             'p1lower','p2lower','p3lower','p1upper','p2upper','p3upper')
    colnames(rm.matrix)<-c('time','mu1','mu2','mu3','mu1se', 'mu2se','mu3se',
                             'mu1lower','mu2lower','mu3lower','mu1upper','mu2upper','mu3upper')
  }

  else if (SE == T & std.err == 'boot'){
    prob.boot<-matrix(NA, nrow = B, ncol = 3*length(prob.times))
    rm.boot<-matrix(NA, nrow = B, ncol = 3*length(mu.times))
    if(!is.na(boot.seed)) set.seed(boot.seed)
    indices<-matrix(sample(nrow(dat), nrow(dat)*B, replace = T), nrow = B)
    # we only estimate SE for prob.times, mu.times
    for (i in 1:B) {
      dat.bs <- dat[indices[i,],]
      # Prep data for Yifei's method
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
      rm.matrix<-NULL
      if (length(mu.times)>0){

        if (boundary=='boundary.kernel' & lambda(0, bandwidth, V.all, N,tau2)<0) {
          dis<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
                       interval = c(0,bandwidth))$root
        }
        else dis<-0
        discon<-c(0,0)
        rm.matrix<-matrix(NA, nrow = length(mu.times), ncol = 3)
        mu.times2<-c(0, mu.times)
        rm.diff<-rep(NA, length(mu.times))
        for (j in 1:length(mu.times)){
          if (dis<=mu.times2[j] | dis>=mu.times2[j+1] ) rm.diff[j]<-integral(f2, xmin = mu.times2[j], xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                                                             tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                                                             reltol = 1e-5)
          else {
            part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                reltol = 1e-5))
            discon[1]<-dis
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-5, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1e-5
            }
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-4, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1e-4
            }
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-3, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1e-3
            }
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-2, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1e-2
            }
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1e-1, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1e-1
            }
            if (inherits(part1, "try-error")){
              part1<-try(integral(f2, xmin = mu.times2[j], xmax = dis-1, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[1]<-dis-1
            }

            part2<-try(integral(f2, xmin = dis, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                reltol = 1e-5))
            discon[2]<-dis
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1e-5, xmax =mu.times2[j+1] , v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1e-5
            }
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1e-4, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1e-4
            }
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1e-3, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1e-3
            }
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1e-2, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1e-2
            }
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1e-1, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1e-1
            }
            if (inherits(part2, "try-error")){
              part2<-try(integral(f2, xmin = dis+1, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5))
              discon[2]<-dis+1
            }

            rm.diff[j]<-part1+part2
          }
        }

        #discon<-c(0,0)
        #if (boundary=='boundary.kernel' & xi(0, bandwidth, V.all, Y.all, N,tau2)<0) {
        #    discon[1]<- uniroot(function (x) xi(x, bandwidth, V.all, Y.all, N,tau2),
        #                        interval = c(0,bandwidth))$root
        #    discon[2]<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
        #                       interval = c(0,bandwidth))$root
        #    discon<-sort(discon)
        #}
        #rm.matrix<-matrix(NA, nrow = length(mu.times), ncol = 3)
        #mu.times2<-c(0, mu.times)
        #rm.diff<-sapply(1:length(mu.times), function(j)
        #    if (discon[2]<=mu.times2[j] | discon[1]>=mu.times2[j+1] ) integral(f2, xmin = mu.times2[j], xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
        #                tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
        #                reltol = 1e-5)
        #    else integral(f2, xmin = mu.times2[j], xmax = discon[1]-1e-4, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
        #                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
        #                  reltol = 1e-5)+
        #        integral(f2, xmin = discon[2]+1e-4, xmax = mu.times2[j+1], v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
        #                 tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
        #                 reltol = 1e-5))
        rm.matrix[,1]<-cumsum(rm.diff)/scale
        rm.matrix[,3]<-mu.times/scale-sapply(1:length(mu.times), function(i)
          summary(fit, rmean =mu.times[i], scale =scale, extend = T)$table[5])
        rm.matrix[,2]<-mu.times/scale-rm.matrix[,1]-rm.matrix[,3]
#        rm.matrix<-t(rm.matrix)
      }
      prob.matrix<-NULL
      if (length(prob.times)>0){
        prob.matrix<-matrix(NA, nrow = length(prob.times), ncol = 3)
        prob.matrix[,1]<-f2(prob.times, v.all = V.all, y.all = Y.all, h = bandwidth, N = N, tau2 = tau2, boundary = boundary,
                            Left=Left, Right= Right, Value=Value)
        prob.matrix[,3]<-1-sapply(1:length(prob.times), function (x)
          Sd(prob.times[x], Value, Left, Right))
        prob.matrix[,2]<-1-prob.matrix[,1]-prob.matrix[,3]
#        prob.matrix<-t(prob.matrix)
      }
      prob.boot[i,]<-t(prob.matrix)
      rm.boot[i,]<-t(rm.matrix)

#      bootsam[i,1:se_length]<-c(prob.matrix,rm.matrix)
    }
    prob.se<-cbind(matrix(apply(prob.boot,2,sd), ncol = 3, byrow = T),
                matrix(apply(prob.boot,2,function(x) quantile(x,probs = 0.025,na.rm = T)), ncol=3, byrow = T),
                matrix(apply(prob.boot,2,function(x) quantile(x,probs = 0.975,na.rm = T)), ncol=3, byrow = T))
    rm.se<-cbind(matrix(apply(rm.boot,2,sd), ncol = 3, byrow = T),
                matrix(apply(rm.boot,2,function(x) quantile(x,probs = 0.025,na.rm = T)), ncol=3, byrow = T),
                matrix(apply(rm.boot,2,function(x) quantile(x,probs = 0.975,na.rm = T)), ncol=3, byrow = T))

    #    res<-rbind(res,apply(bootsam,2,sd),
#               apply(bootsam,2,function(x) quantile(x,probs = 0.025,na.rm = T)),
#               apply(bootsam,2,function(x) quantile(x,probs = 0.975, na.rm = T)))

    prob.matrix<-cbind(prob.times/scale, prob.matrix, prob.se)
    rm.matrix<-cbind(mu.times/scale, rm.matrix, rm.se)

    colnames(prob.matrix)<-c('time','p1','p2','p3','p1se', 'p2se','p3se',
                             'p1lower','p2lower','p3lower','p1upper','p2upper','p3upper')
    colnames(rm.matrix)<-c('time','mu1','mu2','mu3','mu1se', 'mu2se','mu3se',
                           'mu1lower','mu2lower','mu3lower','mu1upper','mu2upper','mu3upper')
  }

  if (incl.discon==T) return(list(prob.info = prob.matrix,
                                  mu.info =rm.matrix, discon = discon.to.ret))
  else return(list(prob.info = prob.matrix, mu.info =rm.matrix))
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
