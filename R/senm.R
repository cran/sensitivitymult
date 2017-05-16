senm <- function (y, z, mset, gamma = 1, inner = 0, trim = 3, lambda = 1/2,
              tau = 0, alternative="greater", TonT = FALSE)
    {
        #Check input
        stopifnot((alternative=="greater")|(alternative=="less"))
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

        if (alternative=="less"){
          ymat<-(-ymat)
          tau<-(-tau)
        }
        if (!(tau == 0)) ymat[, 1] <- ymat[, 1] - tau

        ms <- mscorev(ymat, inner = inner, trim = trim, qu = lambda, TonT = TonT)
        separable1v(ms, gamma = gamma)
    }
