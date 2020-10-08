##### Implements the kernel method with the boundary kernel of Muller 1991 on the LEFT ONLY.

Kq <- function(x, q) {
  sigk1 <- sqrt(0.2)
  2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
    (1 + ((q - 1) / (q + 1) / sigk1) ^ 2 + 2 / sigk1 ^ 2 * (1 - q) / (1 + q) ^ 2 * x)
}

K1 <- function(u) {
  0.75 * (1 - u ^ 2) * (abs(u) < 1)
}

K.unif<-function(u,x,b){ # implements the left boundary ONLY
  #    if (b<=x & x<=(1-b)) return(K1((x-u)/b))
  if (b<=x) return(K1((x-u)/b))
  if (b>x) return(Kq((x-u)/b, x/b))
  #    if (x>(1-b)) return(Kq(-(x-u)/b, (1-x)/b))
}

lambda <- function(v,h,v.all,N,tau2)
{
  # Change visit time scale to the interval 0,1
  v<-v/tau2; h<-h/tau2; v.all<-v.all/tau2

  lambda.v <- rep(0,length(v))
  for(i in 1:length(v))
  {
    lambda.v[i] <- sum(K.unif(v.all,v[i],h))/(N*h*tau2)
  }
  lambda.v
}

xi <- function(v,h, v.all, y.all, N,tau2)
{
  # Change visit time scale to the interval 0,1
  v<-v/tau2; h<-h/tau2; v.all<-v.all/tau2
  xi.v <- rep(0,length(v))
  for(i in 1:length(v))
  {
    xi.v[i] <- sum(K.unif(v.all,v[i],h)*y.all)/(N*h*tau2)
  }
  xi.v
}

Sd <- function(u, value, left, right)
{
  if(u == 0)
  {1}
  else{
    value[which(u<=right & u>left)]
  }
}

f2 <- function(v,h,v.all,y.all,N,tau2, boundary, Left, Right, Value)
{
  lambda.all <- lambda(v, h, v.all, N,tau2)
  xi.all <- xi(v, h, v.all, y.all, N,tau2)
  sd.all <- v
  for(i in 1:length(v))
  {
    sd.all[i] <- Sd(v[i], Value, Left, Right)
  }
  lambda.all[lambda.all == 0] <- Inf

  rh.all <-xi.all/lambda.all
  if (boundary=='interpolation'){
    lambda.h <- lambda(h, h, v.all, N,tau2)
    lambda.h[lambda.h==0] <- Inf
    xi.h <- xi(h, h, v.all, y.all, N,tau2)
    rh.all[v <= h] <- 1+((xi.h/lambda.h-1)/h)*v[v <= h]
  }
  sd.all*rh.all
}

int_f2<-function(lower, upper, bandwidth, V.all,Y.all,N,tau2, boundary, Left, Right, Value, warn = T){
  if (boundary =='interpolation' | lower>bandwidth|
      lambda(lower, bandwidth, V.all, N,tau2)>=0|lambda(upper, bandwidth, V.all, N,tau2)<0){
    return(integral(f2, xmin = lower, xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
             tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
             reltol = 1e-5))
  }
  dis<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
                 interval = c(lower,upper))$root
  leftpow <- -5
  rightpow <- -5
  while (inherits(integral(f2, xmin =lower, xmax = dis-10^(leftpow), v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
           tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
           reltol = 1e-5), "try-error")){
    leftpow <- leftpow+1
  }
  while (inherits(integral(f2, xmin =dis+10^(rightpow), xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                           tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                           reltol = 1e-5), "try-error")){
    rightpow <- rightpow+1
  }
  if (warn == T) warning('Due to a discontinuity in hat r_h(t), the interval[', dis-10^(leftpow),',',dis+10^(rightpow),'] was excluded from integration.')

  return(integral(f2, xmin = lower, xmax = dis-10^(leftpow), v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                  reltol = 1e-5)+
          integral(f2, xmin = dis+10^(rightpow), xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
           tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
           reltol = 1e-5))
}


######## dab function
### Stands for Data Adapative Bandwidth and it finds the bandwidth using the
### method described in Sun, Huang and Wang 2017
dab<-function(lower, upper, V.all, Y.all, ID.all, tau, tau2, length.out = 5000){
  # Find h hat that minimizes cross validation error. If you use a first-order kernel, the
  # resulting h hat satisfies the regularity condition in Sun, Huang and Wang 2017.
  # The resulting h hat can then be used in the estimator, with the kernel of your choice.
  # (We've used the Epanechnikov kernel and the boundary kernel extension of Epanechnikov).

  hhh <- seq(from = lower, to = upper, length.out = length.out)
  UCV.nw <- hhh
  for(i in 1:length(hhh))
  {
    UCV.nw[i] <- sum((Y.all[V.all<tau]-mLeaveone(V.all[V.all<tau],ID.all[V.all<tau],hhh[i],V.all,Y.all,ID.all, tau2))^2)
  }
  hh1 <- hhh[which(UCV.nw  == min(UCV.nw,na.rm = TRUE))]
  hh1 <- max(hh1)
  # plot(hhh,UCV.nw)
  return(hh1)
}

###### dab calls the mLeaveone function (originally from Yifei's data.R file)
mLeaveone <- function(v,id,h,v.all,y.all,id.all, tau2)
{
  K.unif <- function(v,h)
  {
    (v<h & v > 0)/h
  }
  m.v <- rep(0,length(v))
  for(i in 1:length(v))
  {
    if(v[i]<h){
      m.v[i] <- sum(K.unif(h-v.all[id.all != id[i]],h)*y.all[id.all != id[i]])/sum(K.unif(h-v.all[id.all != id[i]],h))
    }else if(v[i]<tau2-h){
      m.v[i] <- sum(K.unif(v[i]-v.all[id.all != id[i]],h)*y.all[id.all != id[i]])/sum(K.unif(v[i]-v.all[id.all != id[i]],h))
    }else{
      m.v[i] <- sum(K.unif(tau2-h-v.all[id.all != id[i]],h)*y.all[id.all != id[i]])/sum(K.unif(tau2-h-v.all[id.all != id[i]],h))
    }
  }
  m.v[is.na(m.v)] <- 0
  m.v
}
