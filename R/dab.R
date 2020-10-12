#' Find bandwidth that minimizes cross validation error
#'
#' This function searches a grid to find the bandwidth that
#' minimizes mean squared leave-one-out cross validation error.
#' The procedure uses a first-order kernel, and thus
#' yields a bandwidth on the order of n^(-1/3). We recommend using a
#' bandwidth on the order of n^(-7/24) in the kernel estimator, which can be
#' obtained by multiplying the bandwidth returned by the \code{dab} function
#' by n^(1/24).
#'
#' @param dat a dataframe with one row per individual with the variables
#' \code{t1-tm} and \code{x1-xm} where \code{m}
#' is the number of visits and \code{ti} and \code{xi} are the time and status at
#' each visit. \code{xi=1} represents alive and event free
#' (state 1 in a multistate illness-death model), and \code{xi=2} represents alive with event.
#' @param lower the lower end of the interval that will be searched
#' @param upper the upper end of the interval that will be searched
#' @param tau the maximum visit time that will be used in the cross-validation
#' @param tau2  the maximum time that individuals were at risk for a visit,
#' based on the study design
#' @param length.out the number of grid points in the interval [lower, upper]
#' that will be considered. Default is 5000.
#'
#' @return the bandwidth minimizes cross-validation error. We recommend multiplying by
#' n^(1/24) to get a bandwidth on the order of n^(-7/24) for the kernel estiamtor.

#' @export
#'
#' @examples mydat <- simdat(50, scale12=1/.0008, scale13=1/.0002, scale23=1/.0016,
#' vital.lfu=c(30.4*36, 30.4*48),
#' visit.schedule = 30.4*c(6, 12, 18, 24, 30, 36, 42, 48), scatter.sd=10)
#' dab(mydat, lower = 30.4, upper = 36*30.4, tau = 30.4*48, tau2=30.4*48,length.out = 100)
dab<-function(dat, lower, upper, tau, tau2, length.out = 5000){

  # Prep data
  nvisits<-dat$nvisits[1]
  ID<- 1:nrow(dat)
  V.all <- c(t(as.matrix(dat[,paste('t',1:nvisits, sep='')])))
  Y.all <- as.numeric(c(t(as.matrix(dat[,paste('x',1:nvisits, sep='')])))==1)
  ID.all <- rep(1:nrow(dat), each = nvisits)
  missing <- is.na(V.all)
  V.all <- V.all[!missing]
  Y.all <- Y.all[!missing]
  ID.all <- ID.all[!missing]

  hhh <- seq(from = lower, to = upper, length.out = length.out)
  UCV.nw <- hhh
  for(i in 1:length(hhh))
  {
    UCV.nw[i] <- sum((Y.all[V.all<tau]-mLeaveone(V.all[V.all<tau],ID.all[V.all<tau],hhh[i],V.all,Y.all,ID.all, tau2))^2)
  }
  hh1 <- hhh[which(UCV.nw  == min(UCV.nw,na.rm = TRUE))]
  hh1 <- max(hh1)
  return(hh1)
}
