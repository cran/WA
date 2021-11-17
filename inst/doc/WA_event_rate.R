## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=F-------------------------------------------------------------------
#  obj<-LRfit(id,time,status,trt)

## ----eval=F-------------------------------------------------------------------
#  obj<-LRfit(id,time,status,trt,Dweight=0)

## ----eval=F-------------------------------------------------------------------
#  # define w_m(t)
#  wH<-function(t,m) return (1)
#  # define w_m^D(t)
#  wD<-function(t,m) return (0)
#  obj<-LRfit(id,time,status,trt,wH=wH,wD=wD)

## ----eval=F-------------------------------------------------------------------
#  summary(obj,tau)

## ----eval=F-------------------------------------------------------------------
#  plot(obj,conf=TRUE)

## ----setup--------------------------------------------------------------------
library(WA)
head(hfaction_cpx12)

## -----------------------------------------------------------------------------
dat<-hfaction_cpx12
obj<-LRfit(dat$id,dat$time,dat$status,dat$trt)
## print some descriptive information 
obj

## ----fig.align='center',fig.width=7.2,fig.height=4.5--------------------------
plot(obj,conf=T,xlab="Time (years)",xlim=c(0, 3.5),ylim=c(0,3))

## -----------------------------------------------------------------------------
## summarize the inference results at tau=3.5 years
summary(obj,tau=3.5,joint.test=T)

## -----------------------------------------------------------------------------
## fit the data using the weighted composite loss
obj2<-LRfit(dat$id,dat$time,dat$status,dat$trt,Dweight=2)
## summarize the results at tau=3.5 years
summary(obj2,tau=3.5)

