senmCI<-function(y, z, mset, gamma=1, inner=0, trim=3, lambda=1/2,
                  alpha=0.05, twosided=TRUE, upper=TRUE, TonT = FALSE){
  stopifnot(gamma>=1)
  stopifnot((inner>=0)&(inner<=trim))
  stopifnot((lambda>0)&(lambda<1))
  stopifnot(is.vector(y)&is.vector(z)&is.vector(mset))
  stopifnot((length(z)==length(y)))
  stopifnot((length(z)==length(mset)))
  stopifnot(all(!is.na(y)))
  stopifnot(all((z==0)|(z==1))) #z is 1 for treated, 0 for control
  tbcheck<-table(z,mset)
  ck<-all(tbcheck[2,]==1)&all(tbcheck[1,]>=1)
  if (!ck){
    warning("Every matched set must contain one treated subject and at least one control.")
    stopifnot(ck)
  }
  stopifnot((alpha>0)&(alpha<1))
  stopifnot(is.logical(twosided))
  stopifnot(is.logical(upper))

  #Convert y to matrix
  mset<-as.integer(mset)
  o<-order(mset,1-z)
  y<-y[o]
  z<-z[o]
  mset<-mset[o]
  tb<-table(mset) #need to check
  nset<-length(tb)
  setsize<-max(tb)

  makeymat<-function(yj){
    ymat<-matrix(NA,nset,setsize)
    m<-0
    for (i in 1:nset){
      ymat[i,1:tb[i]] <- yj[(m+1):(m+tb[i])]
      m<-m+tb[i]
    }
    ymat
  }

  ymat<-makeymat(y)

#Determine a finite interval, (-mx,mx) that will be searched for CI limits
  mx1<-max(ymat[,1],na.rm=TRUE)-min(ymat[,-1],na.rm=TRUE)
  mx2<-max(ymat[,-1],na.rm=TRUE)-min(ymat[,1],na.rm=TRUE)
  mx<-max(mx1,mx2)

# senmInternal is a version of senm() that does not repeat unneeded steps
senmInternal <- function (tau=0, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
                    alternative="greater", TonT = FALSE){
  if (alternative=="less"){
    ymat<-(-ymat)
    tau<-(-tau)
  }
  yadj<-ymat
  yadj[, 1] <- yadj[, 1] - tau
  ms <- mscorev(yadj, inner = inner, trim = trim, qu = lambda, TonT = TonT)
  separable1v(ms, gamma = gamma)
}

#Internal function returning the tau that corresponds with a specified Normal quantile
  solve4tau<-function(dev,alternative="greater"){
    f<-function(taus){
      ntaus<-length(taus)
      o<-rep(NA,ntaus)
      for (i in 1:ntaus){
        d<-senmInternal(tau=taus[i],gamma=gamma,inner=inner,trim=trim,lambda=lambda,alternative=alternative,TonT=TonT)$deviate
        o[i]<-as.vector(d-dev)
      }
      o
    }
    rt<-stats::uniroot(f,c(-mx,mx))
    unlist(rt$root)
  }

  if (twosided) nq<-(-stats::qnorm(alpha/2))
  else nq<-(-stats::qnorm(alpha))
  estL<-solve4tau(0,alternative="greater")
  if (gamma>1) estH<-solve4tau(0,alternative="less")
  else estH<-estL
  if (twosided|upper) ciL<-solve4tau(nq,alternative="greater")
  else ciL<-(-Inf)
  if (twosided|(!upper)) ciH<-solve4tau(nq,alternative="less")
  else ciH<-Inf

  CI<-"Lower"
  if(twosided) CI<-"Two-sided"
  if((!twosided)&upper) CI<-"Upper"
  desc<-c(1-alpha,gamma,CI)
  names(desc)<-c("Coverage","Gamma","Confidence Interval")

  list(PointEstimates=c(estL,estH),ConfidenceInterval=c(ciL,ciH),description=desc)
}
