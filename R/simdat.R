#' Simulate component-wise censored data from an illness-death model
#'
#' This function simulates component-wise censored data from an irreversible or
#' reversible illness-death model. The data is component-wise censored because
#' illness status is intermittently observed at visits, while vital status (alive/dead)
#' is fully observed until right censoring time. All individuals start alive and illness-free.
#' Time from state entry to transition is simulated from a Weibull distribution with
#' user-specified parameters, and the user specifies the visit process.

#' @param n sample size
#' @param scale12 The model has states 1, 2 and 3 representing illness free,
#' alive with illness, and death. \code{scale12} is the scale parameter for the Weibull
#' distribution used to generate the time from entry in state 1 until transition to state 2.
#' Default is 1/0.0008.
#' @param scale13 the scale parameter for the Weibull
#' distribution used to generate the time from entry in state 1 until transition to state 3.
#' Default is 1/0.0002.
#' @param scale23 the scale parameter for the Weibull
#' distribution used to generate the time from entry in state 2 until transition to state 3.
#' Default is 1/0.0016.
#' @param shape12 the shape parameter for the Weibull
#' distribution used to generate the time from entry in state 1 until transition to state 2.
#' Default is 1.
#' @param shape13 the shape parameter for the Weibull
#' distribution used to generate the time from entry in state 1 until transition to state 3.
#' Default is 1.
#' @param shape23 the shape parameter for the Weibull
#' distribution used to generate the time from entry in state 2 until transition to state 3.
#' Default is 1.
#' @param scale21 the scale parameter for the Weibull
#' distribution used to generate the time from entry in state 2 until transition to state 1.
#' Default is \code{NULL}, which implies an irreversible illness-death model. If \code{scale21} is
#' a positive number, the model is reversible, i.e. transitions from state 2 back
#' to state 1 are allowed. This version of \code{simdat} can only handle reversible models
#' where all scale parameters are equal to 1.
#' @param shape21 the shape parameter for the Weibull
#' distribution used to generate the time from entry in state 2 until transition to state 1.
#' Default is 1.
#' @param vital.lfu defines the right censoring time or distribution. If \code{vital.lfu}
#' is a vector, right censoring will occur according to a uniform distribution on
#' the interval defined by the vector. If \code{vital.lfu} is a scalar, all observations
#' will be right-censored at \code{vital.lfu}.
#' @param visit.schedule exactly one of \code{visit.schedule}, \code{visit.rate} and \code{renew.param}
#' should be non-missing. If \code{visit.schedule} is non-missing, visits will be generated
#' according to a truncated normal distribution around each time in \code{visit.schedule},
#' truncated at 3 standard deviations.
#' @param scatter.sd the standard deviation of the truncated normal distribution used to
#' generate visit times
#' @param missing.rate a vector of the same length as \code{visit.schedule},
#' which specifies the proportion of  individuals with missing data at each time
#' in \code{visit.schedule}. If \code{missing.rate=0} (the default), no visits are missed.
#' @param visit.rate exactly one of \code{visit.schedule}, \code{visit.rate} and \code{renew.param}
#' should be non-missing. If \code{visit.rate} is non-missing,  visits will be generated
#' according to a Poisson process with rate function \code{visit.rate} as a function of time,
#' for example, \code{visit.rate = function(t) .0001*t}. For a constant
#' visit rate, use \code{visit.rate = function(x) sapply(x, function(t) C)} where \code{C}
#' is the constant rate.
#' @param renew.param exactly one of \code{visit.schedule}, \code{visit.rate} and \code{renew.param}
#' should be non-missing. If \code{renew.param} is non-missing,  visits will be generated
#' according to a Weibull renewal visit process with shape and scale parameters equal to
#' the first and second entries of \code{renew.param}.
#' @param visit.postprog if equal to 0, that once illness is observed at a visit,
#' no more visits occur for that individual. The default is 1.
#' @param seed If \code{seed} is specified, \code{set.seed(seed)} will be run before
#' generating the data, so the data can be reproduced.
#'
#' @return A dataframe with one row per individual with the variables \code{dtime},
#' \code{dstatus}, \code{state2obs}, \code{laststate1}, \code{t1}-\code{tm}, \code{x1}-
#' \code{xm}, and \code{nvisits}.
#' \code{nvisits} and \code{m} are the largest number of visits observed in the dataset.
#' \code{dstatus=1} represents death and \code{dstatus=0} represents right censoring.
#' \code{ti} and \code{xi} are the time and observed status at the \code{i}-th visit.
#' \code{xi=1} means alive and event free (state 1 in a multistate illness-death model),
#' and \code{xi=2} means alive with event. \code{state2obs} is the time of the first visit
#' in state 2 (\code{Inf} if state 2 was never observed) and \code{laststate1} is the time
#' of the last observation in state 1 before state 2 was observed (\code{0} if an
#' individual was already in state 2 at the first visit). Visits are numbered by
#'  ascending time and missing visits do not get a placeholder. That is, if
#'  \code{visit.schedule} was 30, 60, 90 and 120 days, but the actual visit times
#'  for an individual were \code{45,40, NA, 125} because of normal scatter and missing visits,
#'  the entries \code{t1-t4} will be \code{40, 45, 125, NA}.
#' @export
#' @import stats
#' @examples simdat(50, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
#' vital.lfu=c(30.4*36, 30.4*48),
#' visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48),
#' scatter.sd=10)
simdat<-function(n, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
                 shape12=1, shape13=1, shape23=1,
                 scale21=NULL, shape21=1, vital.lfu=c(30.4*36, 30.4*48),
                 visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48),
                 scatter.sd=10, missing.rate = 0,
                 visit.rate=NA, renew.param = NA,
                 visit.postprog = 1, seed = NULL){
  if (!is.function(visit.rate)+is.na(visit.schedule[1])+is.na(renew.param[1])!=2) stop('You must specify exactly one of visit.schedule, visit.rate and renew.param.')
  if (max(missing.rate)>0 & is.na(visit.schedule[1])) stop('missing.rate must be zero unless you are using the visit.schedule option.')
  if (!is.null(seed)) set.seed(seed)

  # Generate true transition times for an irreversible model
  # Can have weibull distributions
  if (is.null(scale21)){
    sim1<- data.frame(one_to_two=rweibull(n, shape = shape12, scale = scale12),
                      one_to_three=rweibull(n, shape = shape13, scale = scale13))
    sim1$two_to_three<- sim1$one_to_two+rweibull(n, shape = shape23, scale = scale23)
    sim1$dtime <- ifelse(sim1$one_to_three<=sim1$one_to_two,
                         sim1$one_to_three, sim1$two_to_three)
    truetime <- matrix(ncol = 1, nrow = n)
    truestate <- matrix(ncol = 1, nrow = n)
    truetime[,1]<-sim1$one_to_two
    truestate[,1]<-2
  }
  # Reversible must have constant hazards for 1->3 and 2->3 transitions
  # - otherwise computing was too slow
  else if (!is.null(scale21)){
    if (shape13 !=1 | shape23!=1) stop('shape13 and shape23 must be 1.')
    # Make matrices representing true transition times and true states
    # Every time a person enters state 1 or 2, their time from entering that state
    # until transition to the other state comes from the weibull distribution
    truetime <- matrix(rweibull(n, shape = shape12, scale = scale12), ncol = 1, nrow = n)
    truestate <- matrix(2, ncol = 1, nrow = n)
    ntransitions<-1
    while (min(truetime[,ntransitions])<max(vital.lfu)){
      truetime<-cbind(truetime, truetime[,ntransitions]+rweibull(n, shape = shape21, scale = scale21))
      truestate<-cbind(truestate, rep(1,n))
      ntransitions<-ntransitions+1
      truetime<-cbind(truetime, truetime[,ntransitions]+rweibull(n, shape = shape12, scale = scale12))
      truestate<-cbind(truestate, rep(2,n))
      ntransitions<-ntransitions+1
    }
    # Generate each person's death time from a non stationary poisson process where
    # rates come from which state they were in over time (1 or 2)
    poisson_incr<-rexp(n=n)
    sim1<-data.frame(dtime=rep(NA, n))
    for (i in 1:n){
      times<-truetime[i,]
      solvefort<-function(t) {
        gridpts<-c(0, times[times<t], t)
        sum(diff(gridpts)*stepfun(truetime[i,], rep(c(1/scale13, 1/scale23), (ntransitions+1)/2))(gridpts[1:(length(gridpts)-1)]))- poisson_incr[i]
      }
      sim1$dtime[i]<- uniroot(solvefort, interval = c(0, 3000), extendInt = 'upX')$root
    }
  }

  # Right censoring for death
  if (length(vital.lfu)==1) sim1$C_D <- vital.lfu
  else if (length(vital.lfu==2)) sim1$C_D <- runif(n, min = vital.lfu[1], max = vital.lfu[2])
  sim1$dstatus<- as.numeric(sim1$dtime<sim1$C_D)
  sim1$dtime[sim1$dstatus == 0] <- sim1$C_D[sim1$dstatus == 0]

  # Generate visit times from truncated normal distribution
  if (!is.na(visit.schedule[1])){
    for (j in 1:length(visit.schedule)){
      sim1[,paste('t',j, sep='')]  <- visit.schedule[j] + msm::rtnorm(n, 0, scatter.sd, -3*scatter.sd, 3*scatter.sd)
    }
    nvisits<-length(visit.schedule)
    # Make some visits missing, if needed
    if (length(missing.rate)==1 & missing.rate[1] == 0) missing.rate<-rep(0, length(visit.schedule))
    if (length(missing.rate) != length(visit.schedule)) stop ('missing.rate and visit.schedule must have same length.')
    if (max(missing.rate)>0){
      for (j in 1:nvisits){
        missing<- (rbinom(n = n, size = 1, prob = missing.rate[j])==1)
        sim1[missing, paste('t',j, sep='')] <- NA
      }
    }
    # Make visits before time zero missing
    sim1[,paste('t',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]<0]<-NA
    # Sort visits in ascending order
    sim1[,paste('t',1:nvisits, sep='')]<-
      t(apply(sim1[,paste('t',1:nvisits, sep='')], 1, sort, na.last=T))
  }
  # Or, generates visits from a poisson process
  else if (is.function(visit.rate)){
    next_event_time<-function(last_ev_time, visit.rate, poisson_incr){
      solvefort<-function(t) integrate(visit.rate, lower = last_ev_time, upper = t)$value - poisson_incr
      uniroot(solvefort, interval = c(last_ev_time, last_ev_time+6*30.4), extendInt = 'upX')$root
    }
    # Generate each person's first visit
    poisson_incr<-rexp(n=n)
    times<-matrix(ncol = 1, nrow = n)
    times[,1]<-sapply(1:n, function(i) next_event_time(0, visit.rate, poisson_incr[i]))
    nvisits<-1
    while (min(times[,nvisits])<max(vital.lfu)){
      poisson_incr<-rexp(n=n)
      times<-cbind(times, sapply(1:n, function(i) next_event_time(times[i,nvisits],visit.rate,poisson_incr[i])))
      nvisits<-nvisits+1
    }
    sim1[,paste('t', 1:nvisits, sep = '')]<-times
  }
  # Or, generate visits from a weibull renewal process
  else if (!is.na(renew.param[1])){
    times<-matrix(ncol = 1, nrow = n)
    times[,1]<-rweibull(n, shape = renew.param[1], scale = renew.param[2])
    nvisits<-1
    while (min(times[,nvisits])<max(vital.lfu)){
      times<-cbind(times, times[,nvisits]+rweibull(n, shape = renew.param[1], scale = renew.param[2]))
      nvisits<-nvisits+1
    }
    sim1[,paste('t', 1:nvisits, sep = '')]<-times
  }

  # Make visits after end of followup missing
  sim1[,paste('t',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]>sim1$dtime]<-NA

  # Record current state at each visit
  sim1[,paste('x',1:nvisits, sep='')] <-NA
  for (i in 1:n){
    truei<-stepfun(truetime[i,], c(1, truestate[i,]))
    sim1[i,paste('x',1:nvisits, sep='')] <-truei(sim1[i,paste('t',1:nvisits, sep='')] )
  }
  # Record the first time state 2 was observed (observed progression time) and the last time state 1 observed before the first time state 2 was observed
  sim1$state2obs<-sapply(1:n, function(i) min(c(Inf, sim1[i,paste('t',1:nvisits,sep ='')][sim1[i,paste('x',1:nvisits,sep ='')]==2]), na.rm = T))
  sim1$laststate1<-sapply(1:n, function(i) max(c(0, sim1[i,paste('t',1:nvisits,sep ='')][sim1[i,paste('x',1:nvisits,sep ='')]==1 & sim1[i,paste('t',1:nvisits,sep ='')]<sim1$state2obs[i]]), na.rm = T))

  # If visit.postprog=0, remove all visits after progression first observed
  if (visit.postprog == 0){
    sim1[,paste('x',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]>sim1$state2obs]<-NA
    sim1[,paste('t',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]>sim1$state2obs]<-NA
  }
  # Delete variables that represent visits that never happened
  while(sum(!is.na(sim1[,paste('t',nvisits, sep='')]))==0){
    nvisits<-nvisits-1
    if (nvisits == 0) stop('No visits ocurred; change parameters to get component-wise censored data.')
  }
  sim1<-sim1[,c('dtime','dstatus','state2obs','laststate1',
                paste('t',1:nvisits, sep=''),paste('x',1:nvisits, sep=''))]
  sim1$nvisits<-nvisits

  return(sim1)
}
