K.unif<-function(u,x,b){
  # Implements the kernel method with the boundary kernel of Muller 1991 on the LEFT AND RIGHT.
  if (b<=x & x<=(1-b)) return(K1((x-u)/b))
  #if (b<=x) return(K1((x-u)/b))
  if (b>x) return(Kq((x-u)/b, x/b))
  if (x>(1-b)) return(Kq(-(x-u)/b, (1-x)/b))
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
  ##### Don't need any boundary correction on the RHS
  if (upper < tau2 - bandwidth & (boundary =='interpolation' | lower>bandwidth|
                                  lambda(lower, bandwidth, V.all, N,tau2)>=0|lambda(upper, bandwidth, V.all, N,tau2)<0)){
    myint <- try(pracma::integral(f2, xmin =lower, xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5), silent = T)
    if (inherits(myint, "try-error")){
      brkpt <- optimize(function(x) lambda(x, bandwidth, V.all, N, tau2),
                        interval=c(lower, bandwidth), maximum=F)$minimum
      myint <- int_w_discon(lower, brkpt, bandwidth, V.all,Y.all,N,tau2,
                            boundary, Left, Right, Value, warn = warn)+
        int_w_discon(brkpt, upper, bandwidth, V.all,Y.all,N,tau2,
                     boundary, Left, Right, Value, warn = warn)
    }
    return(myint)
  }
  if (upper < tau2 - bandwidth)   return(int_w_discon(lower, upper, bandwidth, V.all,Y.all,N,tau2,
                                                      boundary, Left, Right, Value, warn = warn))
  ##### Don't need boundary correction on the LHS
  if (lower>bandwidth & (lambda(upper, bandwidth, V.all, N,tau2)>=0 |lambda(lower, bandwidth, V.all, N,tau2)<0 )){
    myint <- try(pracma::integral(f2, xmin =lower, xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5), silent = T)
    if (inherits(myint, "try-error")){
      brkpt <- optimize(function(x) lambda(x, bandwidth, V.all, N, tau2),
                        interval=c(tau2-bandwidth, upper), maximum=F)$minimum
      myint <- int_w_discon(lower, brkpt, bandwidth, V.all,Y.all,N,tau2,
                            boundary, Left, Right, Value, warn = warn)+
        int_w_discon(brkpt, upper, bandwidth, V.all,Y.all,N,tau2,
                     boundary, Left, Right, Value, warn = warn)
    }
    return(myint)
  }
  if (lower>bandwidth)   return(int_w_discon(lower, upper, bandwidth, V.all,Y.all,N,tau2,
                                                      boundary, Left, Right, Value, warn = warn))
  ####### If you need boundary correction on both sides, break into two pieces
  # First, the piece from lower to tau2/2 which needs no correction on RHS
  if (boundary =='interpolation' | lambda(lower, bandwidth, V.all, N,tau2)>=0 | lambda(tau2/2, bandwidth, V.all, N,tau2)<0){
    myint <- try(pracma::integral(f2, xmin =lower, xmax = tau2/2, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5), silent = T)
    if (inherits(myint, "try-error")){
      brkpt <- optimize(function(x) lambda(x, bandwidth, V.all, N, tau2),
                        interval=c(lower, tau2/2), maximum=F)$minimum
      myint <- int_w_discon(lower, brkpt, bandwidth, V.all,Y.all,N,tau2,
                            boundary, Left, Right, Value, warn = warn)+
        int_w_discon(brkpt, tau2/2, bandwidth, V.all,Y.all,N,tau2,
                     boundary, Left, Right, Value, warn = warn)
    }
    part1 <- myint
  }
  else part1 <- int_w_discon(lower, tau2/2, bandwidth, V.all,Y.all,N,tau2,
                                                      boundary, Left, Right, Value, warn = warn)
  # Next, the piece from tau2/2 to upper which needs no correction on RHS
  if (lambda(upper, bandwidth, V.all, N,tau2)>=0 | lambda(tau2/2, bandwidth, V.all, N,tau2)<0){
    myint <- try(pracma::integral(f2, xmin =tau2/2, xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5), silent = T)
    if (inherits(myint, "try-error")){
      brkpt <- optimize(function(x) lambda(x, bandwidth, V.all, N, tau2),
                        interval=c(tau2/2, upper), maximum=F)$minimum
      myint <- int_w_discon(tau2/2, brkpt, bandwidth, V.all,Y.all,N,tau2,
                            boundary, Left, Right, Value, warn = warn)+
        int_w_discon(brkpt, upper, bandwidth, V.all,Y.all,N,tau2,
                     boundary, Left, Right, Value, warn = warn)
    }
    part2 <- myint
  }
  else part2 <- int_w_discon(tau2/2, upper, bandwidth, V.all,Y.all,N,tau2,
                             boundary, Left, Right, Value, warn = warn)
  return(part1 + part2)
}


int_w_discon<-function(lower, upper, bandwidth, V.all,Y.all,N,tau2, boundary, Left, Right, Value, warn = T){
  dis<-uniroot(function (x) lambda(x, bandwidth, V.all, N,tau2),
               interval = c(lower,upper))$root
  leftpow <- -5
  rightpow <- -5
  myint1<-try(pracma::integral(f2, xmin =lower, xmax = dis-10^(leftpow), v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                               tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                               reltol = 1e-5), silent = T)
  while (inherits(myint1, "try-error")){
    leftpow <- leftpow+1
    myint1<-try(pracma::integral(f2, xmin =lower, xmax = dis-10^(leftpow), v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                 tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                 reltol = 1e-5), silent = T)
  }
  myint2<-try(pracma::integral(f2, xmin =dis+10^(rightpow), xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                               tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                               reltol = 1e-5), silent = T)
  while (inherits(myint2, "try-error")){
    rightpow <- rightpow+1
    myint2<-try(pracma::integral(f2, xmin =dis+10^(rightpow), xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                 tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                 reltol = 1e-5), silent = T)
  }
  if (warn == T) warning('Due to a discontinuity in hat r_h(t), the interval[', dis-10^(leftpow),',',dis+10^(rightpow),'] was excluded from integration.')
  return(myint1+myint2)
}

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

#### Old functions BEFORE I updated them to incorporate a boundary correction on the RHS
int_f2_old<-function(lower, upper, bandwidth, V.all,Y.all,N,tau2, boundary, Left, Right, Value, warn = T){
  if (boundary =='interpolation' | lower>bandwidth|
      lambda(lower, bandwidth, V.all, N,tau2)>=0|lambda(upper, bandwidth, V.all, N,tau2)<0){
    myint <- try(pracma::integral(f2, xmin =lower, xmax = upper, v.all = V.all, y.all = Y.all, h = bandwidth, N = N,
                                  tau2 = tau2, boundary = boundary, Left=Left, Right= Right, Value=Value,
                                  reltol = 1e-5), silent = T)
    if (inherits(myint, "try-error")){
      brkpt <- optimize(function(x) lambda(x, bandwidth, V.all, N, tau2),
                        interval=c(lower, bandwidth), maximum=F)$minimum
      myint <- int_w_discon(lower, brkpt, bandwidth, V.all,Y.all,N,tau2,
                            boundary, Left, Right, Value, warn = warn)+
        int_w_discon(brkpt, upper, bandwidth, V.all,Y.all,N,tau2,
                     boundary, Left, Right, Value, warn = warn)
    }
    return(myint)
  }
  return(int_w_discon(lower, upper, bandwidth, V.all,Y.all,N,tau2,
                      boundary, Left, Right, Value, warn = warn))
}
