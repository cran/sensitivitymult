principal<-function(y,z,mset,w=NULL,gamma=1,inner=0,trim=3,lambda=0.5,
                     TonT=FALSE,apriori=FALSE,Scheffe=FALSE,detail=FALSE,
                     cor=FALSE){
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
  if (trim==Inf){
    warning("Function principal() requires trim<Inf.  See documentation of trim for discussion of trim==Inf.")
    stopifnot(trim<Inf)
  }
  stopifnot(is.matrix(y)|is.data.frame(y))
  stopifnot(length(z)==(dim(y)[1]))
  stopifnot(length(mset)==(dim(y)[1]))
  nvars<-dim(y)[2]
  if (nvars<2){
    warning("y must have at least two outcomes in two columns.")
    stopifnot(nvars>=2)
  }
  if (!is.null(w)) stopifnot(sum(abs(w))>0)

  mset<-as.integer(mset)
  nset<-length(unique(mset))
  ck<-nset==(dim(tbcheck)[2])
  if (!ck){
    warning("mset should either contain integer values or should be a factor.")
    stopifnot(ck)
  }

  o<-order(mset,1-z)
  y<-y[o,]
  z<-z[o]
  mset<-mset[o]
  tb<-table(mset)
  setsize<-max(tb)

  mscrs<-array(NA,c(nset,setsize,nvars))

  makeymat<-function(yj){
    ymat<-matrix(NA,nset,setsize)
    m<-0
    for (i in 1:nset){
      ymat[i,1:tb[i]] <- yj[(m+1):(m+tb[i])]
      m<-m+tb[i]
    }
    ymat
  }

  for (j in 1:nvars) mscrs[,,j]<-mscorev(makeymat(as.vector(y[,j])),inner=inner,trim=trim,qu=lambda,TonT=TonT)

  mscorsmat<-matrix(as.vector(mscrs),nset*setsize,nvars)
  mscorsmat<-mscorsmat[!is.na(mscorsmat[,1]),] #elimiate NA padding
  prn<-stats::princomp(mscorsmat,cor=cor)
  if (detail) dout<-list(sdev=prn$sdev,center=round(prn$center,7),scale=prn$scale)
  ld<-as.matrix(prn$loadings[1:nvars,1:nvars])
  for (j in 1:nvars) if (ld[1,j]!=0) ld[,j]<-ld[,j]*sign(ld[1,j]) #make loading + for first variable
  pscrs<-array(0,c(nset,setsize,nvars))
  #Use loadings to compute prinicpal component scores
  for (j in 1:nvars){
    for (i in 1:nvars){
      pscrs[,,j]<-pscrs[,,j]+mscrs[,,i]*ld[i,j]
    }
  }

  if(is.null(w)) scoremat<-as.matrix(pscrs[,,1])
  else if (length(w)==1) scoremat<-w*as.matrix(pscrs[,,1])
  else{
    scoremat<-w[1]*as.matrix(pscrs[,,1])
    for (j in 2:(length(w))){
      scoremat<-scoremat+w[j]*as.matrix(pscrs[,,j])
    }
  }

  if (!is.null(colnames(y))) rownames(ld)<-colnames(y)
  if (!is.null(w)) names(w)<-colnames(ld)[1:length(w)]
  if (is.null(w)) scheffe.dimension<-1
  else scheffe.dimension<-length(w)

  deviate<-separable1v(scoremat,gamma=gamma)$deviate

  if (!detail){
    if (Scheffe){
      ScheffePVal<-1-stats::pchisq(max(0,deviate)^2,scheffe.dimension)
      list(deviate=deviate,ScheffePVal=ScheffePVal,scheffe.dimension=scheffe.dimension,weights=w,loadings=ld)
    }
    else if (apriori){
      aprioriPVal<-1-stats::pnorm(deviate)
      list(deviate=deviate,aprioriPVal=aprioriPVal,weights=w,loadings=ld)
    }
    else list(deviate=deviate,weights=w,loadings=ld)
  }
  else{
    if (Scheffe){
      ScheffePVal<-1-stats::pchisq(max(0,deviate)^2,scheffe.dimension)
      list(deviate=deviate,ScheffePVal=ScheffePVal,weights=w,scheffe.dimension=scheffe.dimension,loadings=ld,princomp.detail=dout)
    }
    else if (apriori){
      aprioriPVal<-1-stats::pnorm(deviate)
      list(deviate=deviate,aprioriPVal=aprioriPVal,weights=w,loadings=ld,princomp.detail=dout)
    }
    else list(deviate=deviate,weights=w,loadings=ld,princomp.detail=dout)
  }
}
