
#' @importFrom stats ecdf median
LRfitone=function(id,time,status,wH,wD){
  o<-order(id,time)
  id<-id[o]
  time<-time[o]
  status<-status[o]
  # unique event times
  t<-sort(unique(time[status!=0]))
  t<-t[t>0]
  m<-length(t)

  # compute at-risk prob function pai()
  X=time[status %in% c(0,2)]

  ########Descriptive stats##########
  ## total obs, #events, median FU
  desc<-c(length(unique(id)),sum(status==1),sum(status==2),median(X))
  names(desc)<-c("N","N_Rec","N_D","MedFU")
  ###################################

  # prompt error if length(X)!=n
  pai=ecdf(-X)(-t)
  #dL and S
  uid=unique(id)
  n=length(uid)
  dL=matrix(0,n,m)
  dN=matrix(0,n,m)
  # at-risk process
  Y=matrix(0,n,m)

  # i=1
  # system.time({
  for (i in 1:n){
    Y[i,]=(X[i]>=t)

    ind_i=(id==uid[i])
    time_i=time[ind_i]
    status_i=status[ind_i]

    ind_i_rec=(status_i==1)
    mi=sum(ind_i_rec)
    ind_i_death=(status_i==2)

    if (mi>0){
      dLi=wH(time_i[ind_i_rec],1:mi)
      rank_i_rec=round(m*ecdf(t)(time_i[ind_i_rec]))
      # sum(t<=time_i[ind_i_rec])
      dL[i,rank_i_rec]=dLi
    }
    if (any(ind_i_death)){
      dNi=wD(time_i[ind_i_death])
      rank_i_death=round(m*ecdf(t)(time_i[ind_i_death]))
      dN[i,rank_i_death]=1
      dL[i,rank_i_death]=dL[i,rank_i_death]+dNi

    }
  }
  # })

  # cbind(id,time,status)[id==uid[3],]
  # t[which(dL[3,]==1)]
  # t[which(dN[3,]==1)]

  ### KM estimator for S(.) ####
  # dLambda=dN/(n*matrix(rep(pai,each=n),n,m))
  # dLambda[3,]
  # 1/sum(X>=0.027378)

  dLambda=colMeans(dN)/pai
  St=cumprod(1-dLambda)
  # plot(stepfun(t,c(1,St)),do.points=F,col="red",lwd=2)
  # lines(survfit(Surv(time[status!=1],status[status!=1]>0)~1),conf.int = F)

  ## Conditional cumulative loss function R(.) ####
  dR=colMeans(dL)/pai
  #################
  # Loss rate
  #################
  dt=t-c(0,t[-m])
  # numerator
  muR=cumsum(St*dR)
  muD=cumsum(St*dt)

  #log-loss rate (LR)
  llr=log(muR/muD)
  # llr1=llr
  # plot(t,exp(llr1),type='l')
  # lines(t,exp(llr))

  ############ variance estimation #####################
  paiM<-matrix(rep(pai,each=n),n,m)
  dRM<-matrix(rep(dR,each=n),n,m)
  dLambdaM<-matrix(rep(dLambda,each=n),n,m)
  StM<-matrix(rep(St,each=n),n,m)
  dtM<-matrix(rep(dt,each=n),n,m)

  muRM<-matrix(rep(muR,each=n),n,m)
  muDM<-matrix(rep(muD,each=n),n,m)

  #IF for dR
  dpsiR<-(dL-Y*dRM)/paiM
  #IF for S (without -S(t))
  dpsiD<-(dN-Y*dLambdaM)/paiM
  psi<-t(apply(dpsiD,1,cumsum))
  # colMeans(psi) should be 0

  #IF for log-muR
  IFlmuR<- t(apply(StM*(dpsiR-psi*dRM),1,cumsum))/muRM
  # IF for log-muD
  IFlmuD<- -t(apply(StM*psi*dtM,1,cumsum))/muDM
  #IF for log-LR
  IFllr<- IFlmuR - IFlmuD

  se_llr <- sqrt(colMeans(IFllr^2)/(n-1))
  se_lmuR <- NULL
  se_lmuD <- NULL

    se_lmuR <- sqrt(colMeans(IFlmuR^2)/(n-1))
    se_lmuD <- sqrt(colMeans(IFlmuD^2)/(n-1))


  return(list(llr=llr,t=t,se_llr=se_llr,lmuR=log(muR), se_lmuR=se_lmuR,
              lmuD=log(muD), se_lmuD=se_lmuD,St=St,IFllr=IFllr,
              IFlmuR=IFlmuR,IFlmuD=IFlmuD,desc=desc))
}


#' Estimate the while-alive loss (event) rate
#'
#' Estimate and make inference on the while-alive loss (or event) rate
#' across \eqn{J} groups under a user-specified loss function
#'
#' @param id A vector of id variable.
#' @param time A vector of follow-up times.
#' @param status A vector of event type variable; 1 = recurrent event, 2 = death,
#' and 0 = censoring.
#' @param trt A vector of categorical (binary or multiclass)
#'  variable for treatment group.
#' @param Dweight A non-negative weight for death relative to the
#' recurrent event; Default is 0.
#' @param wH A function of \eqn{t} and \eqn{m} to weight recurrent event;
#' \eqn{t}: time; \eqn{m}: existing number of recurrent event; Default is
#' the constant function of 1.
#' @param wD A function of \eqn{t} and \eqn{m} to weight death;
#' \eqn{t}: time; \eqn{m}: existing number of recurrent event; Default is
#' the constant function of 0.
#' @return An object of class \code{LRfit}. See \code{\link{LRfit.object}}
#' for details.
#' @examples
#'# load the HF-ACTION trial data
#'head(hfaction_cpx12)
#'# fit the data
#'dat<-hfaction_cpx12
#'obj<-LRfit(dat$id,dat$time,dat$status,dat$trt)
#'# print the event numbers by group
#'obj
#'# summarize the inference results for tau=3.5 years
#'# with joint test with RMST
#'summary(obj,tau=3.5,joint.test=TRUE)
#'# plot the estimated survival-completed cumulative loss
#'# by group, with 95% confidence intervals
#'plot(obj,conf=TRUE,xlab="Time (years)",xlim=c(0, 3.5),ylim=c(0,3),
#'     ylab="Survival-completed cumulative frequency")
#' @keywords LRfit
#' @export
#' @aliases LRfit
#' @seealso \code{\link{LRfit.object}},
#' \code{\link{summary.LRfit}}, \code{\link{plot.LRfit}}.
LRfit=function(id,time,status,trt,Dweight=0,wH=NULL,wD=NULL){

  if (is.null(wH)){
    wH=function(t,k){
      return(1)
    }
  }
  if (is.null(wD)){
    wD=function(t){
      return(Dweight)
    }
  }


  group<-sort(unique(trt))

  LRfit_g=function(g){
    return(LRfitone(id[trt==g],time[trt==g],status[trt==g],wH,wD))
  }
  obj<-sapply(group,LRfit_g)

  colnames(obj)<-group

  obj<-list(obj,match.call())
  names(obj)<-c("content","call")

  class(obj)<-"LRfit"
  return(obj)

}


