# This function simulates data from an illness-death model. All patients start
# alive and illness-free, and time to illness, time to death without illness, and
# time from illness to death follow a weibull distributions.
# Illness status can be oberved only at visits, which occur according to either
# a truncated normal visit process, a (optionally non-stationary) poisson
# process or a weibull renewal process.
# There is an option to have visits stop after illness is observed.
# Vital status (Alive/dead) is observed until a right censoring time, which
# can be fixed or random.
# There is also an option to make the transition health to illness reversible
# by specifying a value for scale21 greater than zero.

# Mssing visits does not work with the aj funciton, it ignores all data after a missing visit
# msm may do this as well?

########################## Function parameters #################################
# seed is the random seed
# Specify one of visit.schedule (for truncated normal), visitrate (for stationary
# or non-stationary poisson visit process) and renew_param (for weibull renewal process).
# visit.schedule gives the times in days when visits are planned; they occur
# at the planned times +/- N(mean = 0, SD=scattersd) days. The normal
# distribution is truncated at -3*scattersd and 3*scattersd.
# If scattersd = 'nogaps', the SD will be half the time to the first visit,
# and the first visit time distribution will be truncated at 0, t1+1SD.
# Other visit time distributions will be truncated halfway between visits.
# Alternatively, you can specify a rate function for the visit process and the
# visits will occur according to that process. If it's constant, visits will occur
# according to a poisson process. If you want a fixed rate, you specify visitrate
# as function(x) sapply(x, function(t) C) where C is the constant rate.
# Othewise, you can do, for example, visit.rate = function(x) .0001*x.
# Finally, you can specify renew_param, which is a vector of the shape, scale
# parameters for a weibull renewal visit process.
# The default shape and scale parameters mean that there are constant
# transition intensities of .0008 (1->2), .0002 (1->3), .0016 (2->3).
# Setting visitpostprog=0 means that once a person has illness observed at a visit,
# no more visits will occur.
# missing_1 specifies what percent of people missed the first visit. The rate
# doubles at each subsequent visit
# vital.lfu allows you to set a common right censoring time for death.
# If vital.lfu is a vector, right censoring will occur according to a
# uniform distribution on the interval defined by the vector.
# Visits cannot occur after right censoring time.

simdat<-function(seed, visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48), n=250,
                 missing_1 = 0, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
                 shape12=1, shape13=1, shape23=1, scattersd=10,
                 shape21=1, scale21=NULL,
                 visitrate=NA, renew_param = NA, visitpostprog = 1, vital.lfu=30.4*60){
  if (!is.function(visitrate)+is.na(visit.schedule[1])+is.na(renew_param[1])!=2) stop('You must specify exactly one of visit.schedule, visitrate and renew_param.')
  set.seed(seed)

  # Irreversible can have weibull distributions
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

  # Generate visit times
  if (!is.na(visit.schedule[1])){
    if (scattersd=='nogaps'){
      # nogaps means the SD will be half the time to the first visits,
      # and the first visit time distribution will be truncated at 0, (midpoint of visit 1 and 2)
      # The upper truncation for the last visit will be visit time + 1SD.
      # Other visit time distributions will be truncated halfway between visits.
      scattersd<-visit.schedule[1]/2
      sim1[,'t1']<-rtnorm(n, visit.schedule[1], scattersd, 0,
                          ifelse(length(visit.schedule)>1, (visit.schedule[1]+visit.schedule[2])/2, visit.schedule[1]+scattersd))
      if (length(visit.schedule)>1){
        for (j in 2:length(visit.schedule)){
          sim1[,paste('t',j, sep='')]  <- rtnorm(n, visit.schedule[j], scattersd, (visit.schedule[j-1]+visit.schedule[j])/2,
                                                 ifelse(length(visit.schedule)>j, (visit.schedule[j]+visit.schedule[j+1])/2, visit.schedule[j]+scattersd))
        }
      }
    }
    else{ # If scattersd is a number, truncate at -3*SD, 3*SD
      for (j in 1:length(visit.schedule)){
        sim1[,paste('t',j, sep='')]  <- visit.schedule[j] + rtnorm(n, 0, scattersd, -3*scattersd, 3*scattersd)
      }
    }
    nvisits<-length(visit.schedule)
    # sort them in ascending order
    sim1[,paste('t',1:nvisits, sep='')]<-
      t(apply(sim1[,paste('t',1:nvisits, sep='')], 1, sort))
  }
  if (is.function(visitrate)){
    # Function that generates from non-stationary poisson process
    next_event_time<-function(last_ev_time, visitrate, poisson_incr){
      solvefort<-function(t) integrate(visitrate, lower = last_ev_time, upper = t)$value - poisson_incr
      uniroot(solvefort, interval = c(last_ev_time, last_ev_time+6*30.4), extendInt = 'upX')$root
    }
    # Generate each person's first visit
    poisson_incr<-rexp(n=n)
    times<-matrix(ncol = 1, nrow = n)
    times[,1]<-sapply(1:n, function(i) next_event_time(0, visitrate, poisson_incr[i]))
    nvisits<-1
    while (min(times[,nvisits])<max(vital.lfu)){
      poisson_incr<-rexp(n=n)
      times<-cbind(times, sapply(1:n, function(i) next_event_time(times[i,nvisits],visitrate,poisson_incr[i])))
      nvisits<-nvisits+1
    }
    sim1[,paste('t', 1:nvisits, sep = '')]<-times
  }
  if (!is.na(renew_param[1])){
    times<-matrix(ncol = 1, nrow = n)
    times[,1]<-rweibull(n, shape = renew_param[1], scale = renew_param[2])
    nvisits<-1
    while (min(times[,nvisits])<max(vital.lfu)){
      times<-cbind(times, times[,nvisits]+rweibull(n, shape = renew_param[1], scale = renew_param[2]))
      nvisits<-nvisits+1
    }
    sim1[,paste('t', 1:nvisits, sep = '')]<-times
  }

  # record current state at each visit
  sim1[,paste('x',1:nvisits, sep='')] <-NA
  for (i in 1:n){
    truei<-stepfun(truetime[i,], c(1, truestate[i,]))
    sim1[i,paste('x',1:nvisits, sep='')] <-truei(sim1[i,paste('t',1:nvisits, sep='')] )
  }
  # Make visits after end of followup, before 0 years missing
  for (j in 1:nvisits){
    sim1[,paste('x',j, sep='')][sim1[,paste('t',j, sep='')]<0] <- NA
    sim1[,paste('t',j, sep='')][sim1[,paste('t',j, sep='')]<0] <- NA
    sim1[,paste('x',j, sep='')][sim1[,paste('t',j, sep='')]>sim1$dtime] <- NA
    sim1[,paste('t',j, sep='')][sim1[,paste('t',j, sep='')]>sim1$dtime] <- NA
  }
  # make some random visits missing, if needed
  if (missing_1!= 0 & is.function(visitrate)) stop('Missing visit rate must be zero for visit.rate specification.')
  if (missing_1>0){
    for (j in 1:nvisits){
      missing<- (rbinom(n = n, size = 1, prob = min(1,(missing_1/100)*(2^(j-1))))==1)
      sim1[missing, paste('x',j, sep='')] <- NA
      sim1[missing, paste('t',j, sep='')] <- NA
    }
  }
  # The first time state 2 was observed (observed progression time), last time state 1 observed
  sim1$state2obs<-sapply(1:n, function(i) min(c(Inf, sim1[i,paste('t',1:nvisits,sep ='')][sim1[i,paste('x',1:nvisits,sep ='')]==2]), na.rm = T))
  sim1$laststate1<-sapply(1:n, function(i) max(c(0, sim1[i,paste('t',1:nvisits,sep ='')][sim1[i,paste('x',1:nvisits,sep ='')]==1 & sim1[i,paste('t',1:nvisits,sep ='')]<sim1$state2obs[i]]), na.rm = T))

  # If visitpostprog=0, remove all visits after progression
  if (visitpostprog == 0){
    sim1[,paste('x',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]>sim1$state2obs]<-NA
    sim1[,paste('t',1:nvisits,sep ='')][sim1[,paste('t',1:nvisits,sep ='')]>sim1$state2obs]<-NA
  }
  # delete variables that represent visits that never happened
  while(sum(!is.na(sim1[,paste('t',nvisits, sep='')]))==0){
    nvisits<-nvisits-1
  }

  sim1<-sim1[,c('dtime','dstatus','state2obs','laststate1',
                paste('t',1:nvisits, sep=''),paste('x',1:nvisits, sep=''))]
  sim1$nvisits<-nvisits
  return(sim1)
}

#dat<-simdat(9538, n=5, scale21=NULL)
#dat2<-simdat(9538, n=25, scale21=1/.0016)


