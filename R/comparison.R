comparison<-function(y,z,mset,w,gamma=1,inner=0,trim=3,lambda=0.5,
                     TonT=FALSE,apriori=FALSE,Scheffe=FALSE){
  #Check input
  stopifnot((0 <= inner) & (inner <= trim))
  stopifnot((inner == 0)|(trim<Inf))
  stopifnot((lambda > 0) & (lambda < 1))
  stopifnot(gamma >= 1)
  stopifnot((TonT==TRUE)|(TonT==FALSE))
  stopifnot((apriori==TRUE)|(apriori==FALSE))
  stopifnot((Scheffe==TRUE)|(Scheffe==FALSE))
  stopifnot(all(as.vector(!is.na(y))))
  stopifnot(all((z==0)|(z==1))) #z is 1 for treated, 0 for control
  tbcheck<-table(z,mset)
  ck<-all(tbcheck[2,]==1)&all(tbcheck[1,]>=1)
  if (!ck){
    warning("Every matched set must contain one treated subject and at least one control.")
    stopifnot(ck)
  }
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(length(z)==(dim(y)[1]))
  stopifnot(length(mset)==(dim(y)[1]))
  nvars<-dim(y)[2]
  if (nvars<2){
    warning("y must have at least two outcomes in two columns.")
    stopifnot(nvars>=2)
  }
  stopifnot(length(w)==(dim(y)[2]))
  stopifnot(sum(abs(w))>0)
  if (!is.null(colnames(y))) names(w)<-colnames(y)


  mset<-as.integer(mset)
  o<-order(mset,1-z)
  y<-y[o,]
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

  scoremat<-w[1]*mscorev(makeymat(as.vector(y[,1])),inner=inner,trim=trim,qu=lambda,TonT=TonT)
  for (j in 2:nvars){
    scoremat<-scoremat+w[j]*mscorev(makeymat(as.vector(y[,j])),inner=inner,trim=trim,qu=lambda,TonT=TonT)
  }

  deviate<-separable1v(scoremat,gamma=gamma)$deviate
  if (Scheffe){
    ScheffePVal<-1-stats::pchisq(max(0,deviate)^2,nvars)
    list(deviate=deviate,ScheffePVal=ScheffePVal,weights=w)
  }
  else if (apriori){
    aprioriPVal<-1-stats::pnorm(deviate)
    list(deviate=deviate,aprioriPVal=aprioriPVal,weights=w)
  }
  else list(deviate=deviate,weights=w)
}
