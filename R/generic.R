

#' Print a short summary of LRfit objects
#'
#' Print the results for the restricted mean times in favor of treatment.
#'
#' @param x An object returned by \code{\link{LRfit}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
print.LRfit=function(x,...){
  cat("Call:\n")
  print(x$call)
  # cat("\n")
  # get descriptive statistics
  x<-x$content
  J<-ncol(x)
  desc<-matrix(unlist(x["desc",]),J,4,byrow=TRUE)
  rownames(desc)<-colnames(x)
  colnames(desc)<-c("N","Rec. event","Death","Med. Follow-up")
  print(desc)
}




#' Summary of the analysis results
#'
#' Summarize the inferential results for group-specific while-alive loss rates
#' and the loss rate ratios at a user-specified length of follow-up.
#'
#'
#' @param object An object returned by \code{\link{LRfit}}.
#' @param tau A positive real number for the follow-up time; Default is the maximum
#' event time in the data.
#' @param ref The label of the group to be set as the reference; Default is the
#' first level by numerical or alphabetical order
#' @param joint.test If TRUE, joint analysis with the restricted mean
#' survival time (RMST) will be performed; Default is FALSE.
#' @param ...  Additional arguments affecting the summary produced.
#'
#' @return An object of class \code{summary.LRfit} with components
#'
#' \item{LRtab}{A \eqn{(J\times 4)}-dimensional matrix containing the inference
#' results for the log-loss rate; Columns include
#' \code{Estimate}, \code{Std.Err}, \code{Z value}, and \code{Pr(>|z|)}.}
#'
#' \item{Rtab}{A \eqn{(J\times 4)}-dimensional matrix containing the inference
#' results for the log-raw cumulative loss if \code{joint.test=TRUE}; Columns include
#' \code{Estimate}, \code{Std.Err}, \code{Z value}, and \code{Pr(>|z|)}.}
#'
#' \item{Dtab}{A \eqn{(J\times 4)}-dimensional matrix containing the inference
#' results for the log-RMST if \code{joint.test=TRUE}; Columns include
#' \code{Estimate}, \code{Std.Err}, \code{Z value}, and \code{Pr(>|z|)}.}
#'
#' \item{LRpval}{\eqn{p}-value for the \eqn{(J-1)}-df chi-square test of group difference in the
#' loss rate.}
#'
#' \item{LRDpval}{\eqn{p}-value for the \eqn{2(J-1)}-df joint chi-square test of group difference in the
#' loss rate and RMST.}
#'
#' \item{...}{}
#' @seealso \code{\link{LRfit}}, \code{\link{plot.LRfit}}.
#' @keywords LRfit
#' @importFrom stats pchisq pnorm
#' @examples
#' #See examples for LRfit().
#' @export
summary.LRfit=function(object,tau=NULL,ref=NULL,joint.test=FALSE,...){

  x<-object$content
  call<-object$call

  # number of groups
  J<-ncol(x)
  # minimum and maximum taus allowed
  t.min<-max(sapply(x["t",],min))
  t.max<-max(sapply(x["t",],max))


  #### check on the input value tau ##########
  # If null, set tau=t.max
  if (is.null(tau)){
    tau=t.max
  }else{
     if (tau<t.min||tau>t.max){
       tau=t.max
    #   stop(
    #     paste0("tau value is outside the range of empirical observations.\n"),
    #     "Specify a tau within [",t.min,", ",t.max,"].\n")
     }
  }
  ###########################################

  #### check on the reference level ##########
  groups<-colnames(x)

  if (is.null(ref)){
    ref<- groups[1]
  }else{
    ref<-as.character(ref)
    # the specified ref must be one of the group levels
    if (!ref %in% groups){
      stop("Reference level specified is not a valid level of the
          treatment variable.\n")
    }else{
      # regroup x with data from the ref level
      # in the first column
      x<-cbind(x[,ref],x[,colnames(x)!=ref])
      colnames(x)<-c(ref,groups[groups!=ref])
      groups<-colnames(x)
    }
  }

  # ref=1
  # colnames(x)-ref
  ####################


  ### Output table for log-LRR, log-raw loss rate, and log-RMST ratio ####
  LRtab=matrix(NA,J,4)
  Rtab=matrix(NA,J,4)
  Dtab=matrix(NA,J,4)

  rownames(LRtab)=c(paste0("Ref (Group ",groups[1],")"),
                    paste0("Group ",groups[2:J]," vs ",ref))
  colnames(LRtab)=c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")


  rownames(Rtab)=c(paste0("Ref (Group ",groups[1],")"),
                   paste0("Group ",groups[2:J]," vs ",ref))
  colnames(Rtab)=c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")
  rownames(Dtab)=c(paste0("Ref (Group ",groups[1],")"),
                   paste0("Group ",groups[2:J]," vs ",ref))
  colnames(Dtab)=c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")

  #####################################
  ## Fill in the output table for log-LRR
  #####################################

  ### (J-1)-degree of freedom test ####
  llr<-rep(NA,J) ## J-vector of log-LR
  se_llr<-rep(NA,J) ## J-vector of se of log-LR

  ### Fill in reference group ###
  i1<-sum(x[,1]$t<=tau) # index for tau
  # output table
  LRtab[1,1]<-x[,1]$llr[i1]
  LRtab[1,2]<-x[,1]$se_llr[i1]
  # vectors of log-LR and se
  llr[1]<-x[,1]$llr[i1]
  se_llr[1]<-x[,1]$se_llr[i1]


  ## 1. Coefficient matrix for joint test ###
  # M: (J-1) x J
  # -1  1
  # -1    1
  # ...     1
  # -1        1
  #
  ## 2. Fill vectors of log-LR and se

  M<-matrix(0,J-1,J)
  M[,1]=rep(-1,J-1)
   for (j in 2:J){
     M[j-1,j]=1
     ij<-sum(x[,j]$t<=tau)
     llr[j]<-x[,j]$llr[ij]
     se_llr[j]<-x[,j]$se_llr[ij]
   }

  LRtab[2:J,1]=M%*%llr
  LRtab[2:J,2]=sqrt(diag(M%*%diag(se_llr^2)%*%t(M)))

  ###### output table z score and p-value ###
  LRtab[,3]<-LRtab[,1]/LRtab[,2]
  LRtab[,4]<-2*(1-pnorm(abs(LRtab[,3])))
  # LRtab[1,4]=NA
  ### (J-1)-d.f. joint test on LRR
  LRchisq<-t(M%*%llr)%*%solve(M%*%diag(se_llr^2)%*%t(M))%*%(M%*%llr)
  LRpval<-1-pchisq(LRchisq,J-1)
  ##########################################

  #########################
  # Joint test with RMST  #
  #########################
  LRDchisq<-NULL
  LRDpval<-NULL

  if (joint.test){
    lmuR<-rep(NA,J)
    lmuD<-rep(NA,J)
    var_lmuR<-rep(NA,J)
    var_lmuD<-rep(NA,J)

    # vector of log-LR and log-RMST
    # llr_lmuD=c(llr[1],lmuD[1],llr[2],lmuD[2],...,llr[J],lmuD[J])
    llr_lmuD<-rep(NA,2*J)

    i1<-sum(x[,1]$t<=tau)
    lmuR[1]<-x[,1]$lmuR[i1]
    lmuD[1]<-x[,1]$lmuD[i1]



    # fill in the reference group
    llr_lmuD[1:2]<-c(llr[1],lmuD[1])

    # Influence function for llr and lmuD (ref group)
    IF1<-cbind(x[,1]$IFllr[,i1], x[,1]$IFlmuR[,i1],x[,1]$IFlmuD[,i1])
    n1<-nrow(IF1)
    # variance of llr, lmuR, and lmuD for ref group
    V1<-t(IF1)%*%IF1/(n1*(n1-1))

    # coefficient matrix for joint tests
    # M2: 2(J-1) x 2J
    # -1    1
    #    -1   1
    # -1        1
    #    -1       1
    # ...
    # -1              1
    #    -1             1
    #############################

    M2<-matrix(0,2*(J-1),2*J)
    M2[seq(1,2*J-3,by=2),1]<- -1
    M2[seq(2,2*J-2,by=2),2]<- -1

    # variance matrix for (llr, lmuD)
    S_LRD<-matrix(0,2*J,2*J)
    S_LRD[1:2,1:2]<-V1[c(1,3),c(1,3)]
    var_lmuR[1]<-V1[2,2]
    var_lmuD[1]<-V1[3,3]

    # output table for raw loss and RMST
    Rtab[1,1]<-lmuR[1]
    Rtab[1,2]<-sqrt(var_lmuR[1])
    Dtab[1,1]<-lmuD[1]
    Dtab[1,2]<-sqrt(var_lmuD[1])


    for (j in 2:J){
      ij<-sum(x[,j]$t<=tau)
      lmuR[j]<-x[,j]$lmuR[ij]
      lmuD[j]<-x[,j]$lmuD[ij]

      # fill in llr and lmuD
      llr_lmuD[(2*j-1):(2*j)]<-c(llr[j],lmuD[j])
      # influence functions for group j
      IFj<-cbind(x[,j]$IFllr[,ij],x[,j]$IFlmuR[,ij],x[,j]$IFlmuD[,ij])
      nj<-nrow(IFj)
      # variance of llr, lmuR, lmuD for sample j
      Vj<-t(IFj)%*%IFj/(nj*(nj-1))

      # coefficient matrix for joint tests
      M2[(2*j-3):(2*j-2),(2*j-1):(2*j)]<-diag(rep(1,2))
      # variance matrix for (llr, lmuD)
      S_LRD[(2*j-1):(2*j),(2*j-1):(2*j)]<-Vj[c(1,3),c(1,3)]
      var_lmuR[j]<-Vj[2,2]
      var_lmuD[j]<-Vj[3,3]

      Rtab[j,1]<-lmuR[j]
      Rtab[j,2]<-sqrt(var_lmuR[j])
      Dtab[j,1]<-lmuD[j]
      Dtab[j,2]<-sqrt(var_lmuD[j])

    }
    LRDchisq <-t(M2%*%llr_lmuD)%*%solve(M2%*%S_LRD%*%t(M2))%*%(M2%*%llr_lmuD)
    LRDpval <- 1-pchisq(LRDchisq,2*(J-1))

    # complete the output tables for raw loss and RMST
    Rtab[2:J,1]=M%*%lmuR
    Rtab[2:J,2]=sqrt(diag(M%*%diag(var_lmuR)%*%t(M)))
    Rtab[,3]<-Rtab[,1]/Rtab[,2]
    Rtab[,4]<-2*(1-pnorm(abs(Rtab[,3])))
    # Rtab[1,4]=NA

    Dtab[2:J,1]=M%*%lmuD
    Dtab[2:J,2]=sqrt(diag(M%*%diag(var_lmuD)%*%t(M)))
    Dtab[,3]<-Dtab[,1]/Dtab[,2]
    Dtab[,4]<-2*(1-pnorm(abs(Dtab[,3])))
    # Dtab[1,4]=NA
    }

  # printCoefmat(Rtab, P.values=TRUE, has.Pvalue=TRUE)
  result<-list(tau=tau,LRtab=LRtab,Rtab=Rtab,Dtab=Dtab,J=J,
               LRchisq=LRchisq,LRpval=LRpval,
               LRDchisq=LRDchisq,LRDpval=LRDpval,call=call,
               joint.test=joint.test)
  class(result)<-"summary.LRfit"
  return(result)



}





#' Print method for summary.LRfit objects
#'
#' Produces a printed summary of the results for the while-alive loss rate
#'
#' @param x An object returned by \code{\link{summary.LRfit}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
#' @importFrom stats printCoefmat
print.summary.LRfit=function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  tau=x$tau
  joint.test=x$joint.test
  J<-x$J
  # p-value for the (J-1)-d.f. test on loss rate
  LRchisq=x$LRchisq
  LRpval=x$LRpval

  # table for loss rate
  LRtab<-x$LRtab
  cat("Analysis of log loss rate (LR) by tau = ",tau, ":\n",sep="")
  printCoefmat(LRtab, P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")

  cat("Test of group difference in while-alive LR\n")
  cat("X-squared = ", LRchisq, ", df = ",J-1,", p = ",LRpval, sep="")

  ## exponentiated table for loss rate ratio
  za<-qnorm(0.975)
  beta<-LRtab[2:J,1]
  se<-LRtab[2:J,2]
  LRR<-cbind(exp(beta),exp(beta-za*se),exp(beta+za*se))
  colnames(LRR)<-c("LR ratio","95% lower CL","95% higher CL")
  rownames(LRR)<-rownames(LRtab)[2:J]
  cat("\n\n")
  cat("Point and interval estimates for the LR ratio:\n")
  print(LRR)
  if (joint.test){
    # table for RMST
    Dtab<-x$Dtab
    LRDchisq<-x$LRDchisq
    LRDpval<-x$LRDpval
    cat("\n\n")
    cat("Analysis of log RMST (restricted mean survival time) by tau = ",tau, ":\n",sep="")
    printCoefmat(Dtab, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n\n")
    cat("Test of group difference in while-alive LR and RMST\n")
    cat("X-squared = ", LRDchisq, ", df = ",2*(J-1),", p = ",LRDpval, sep="")

  }




}




#' Plot the estimated survival-completed cumulative loss curve
#'
#' Plot the estimated survival-completed cumulative loss (while-alive
#' loss rate times the length of follow-up) as a function of
#' the time horizon.
#'
#' @param x An object returned by \code{\link{LRfit}}.
#' @param group Specifies the group to be plotted.
#' @param conf If TRUE, 95\% confidence limits for the target curve are overlaid.
#' @param group.col A vector of colors for the group-specific curves; must be commensurate
#' with the number of groups.
#' @param main A main title for the plot.
#' @param xlim The x limits of the plot.
#' @param ylim The y limits of the plot.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of y.
#' @param conf.lty Line type for the confidence limits if \code{conf=TRUE}.
#' @param lwd Line width.
#' @param legend If TRUE, a crude legend for the group-specific curves will appear
#' on the bottom right corner of the graph.
#' @param ... Other arguments that can be passed to the underlying \code{plot} method.
#' @return No return value, called for side effects.
#' @seealso \code{\link{LRfit}}, \code{\link{summary.LRfit}}.
#' @keywords LRfit
#' @importFrom stats qnorm
#' @importFrom graphics lines
#' @export
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
plot.LRfit=function(x,group=NULL,conf=FALSE,main=NULL,xlim=NULL, ylim=NULL,xlab="Follow-up time",
                     ylab="Survival-completed cumulative loss",group.col=NULL,conf.lty=3,lwd=2,
                    legend=TRUE,...){
  x<-x$content
  if (is.null(group)){
    group<-colnames(x)
  }
  group<-as.character(group)
  # number of groups to be plotted
  J<-length(group)

  if (is.null(group.col)){
    group.col<-c("red","blue","purple","yellow","green","brown")
  }

  if (length(group.col)<J){
    stop(paste0("Please specify a ", J,"-vector of colors in the group.col argument."))
  }


  # get the data for the plot
  results<-list()
  xmax<-0
  ymax<-0
  za<-qnorm(0.975)

  for (j in 1:J){
    xj<-x[,group[j]]
    tj<-xj$t
    Lj<-tj*exp(xj$llr)
    sej<-xj$se_llr
    Lj_lo<-Lj*exp(-za*sej)
    Lj_hi<-Lj*exp(za*sej)

    results[[j]]=rbind(tj,Lj,Lj_lo,Lj_hi)
    xmax<-max(xmax,tj)
    if (conf){
      ymax<-max(ymax,Lj_hi)
    }else{
      ymax<-max(ymax,Lj)
    }
  }

  if(is.null(xlim)){
    xlim<-c(0,xmax)
  }
  if(is.null(ylim)){
    ylim<-c(0,ymax)
  }

  plot(NULL,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,...)

  # plot the curves
  for (j in 1:J){
    rj<-results[[j]]
    tj<-c(0,rj["tj",])
    Lj<-c(0,rj["Lj",])
    Lj_lo<-c(0,rj["Lj_lo",])
    Lj_hi<-c(0,rj["Lj_hi",])

    lines(tj,Lj,lty=1,lwd=lwd,col=group.col[j])
    if (conf){
      lines(tj,Lj_lo,lty=conf.lty,lwd=lwd,col=group.col[j])
      lines(tj,Lj_hi,lty=conf.lty,lwd=lwd,col=group.col[j])
    }
  }
  # if legend is true, print legend on the bottom right corner
if (legend){
  legend("bottomright",lty=1,lwd=lwd,col=group.col[1:J],group)
}
}

