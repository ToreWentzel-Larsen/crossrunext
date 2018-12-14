library(Rmpfr)

# auxiliary functions:
# matrix of cumulative probabilities for L for each value of C:
cumsumm <- function(mtrx) {
  mtrs <- mtrx
  nrw <- nrow(mtrx)
  for (rw in 1:nrw) mtrs[rw,] <- cumsum(mtrx[rw,])
  mtrs
} # end function cumsumm
cumsummcol <- function(mtrx) {
  mtrs <- mtrx
  ncl <- ncol(mtrx)
  for (cl in 1:ncl) mtrs[,cl] <- cumsum(mtrx[,cl])
  mtrs
} # end function cumsummcol
# matrix of "box" probabilities ("times" scale), P(C>=c, L<=l):
boxprobt <- function(mtrx) {
  mtrs <- mtrx
  rw <- nrow(mtrx)
  cl <- ncol(mtrx)
  for (row in 1:rw) mtrs[row,] <- cumsum(mtrx[row,])
  for (col in 1:cl) mtrs[,col] <- rev(cumsum(rev(mtrs[,col])))
  return(mtrs)
} # end function boxprobt

# function for computing the simultaqneous distribution for number of
# crossings and longest run in n independent binomial observations with
# probability 1/2 (the symmetrixc case). Probabilities are multiplied 
# by m^(n-1) where the multiplier m by default is 2. The computations 
# use the package Rmpfr to enhance accuracy.
crossrunsymm <- function(nmax=100, mult=2, check=0, prec=120, printn=FALSE) {
  nill <- mpfr(0,prec)
  one <- mpfr(1,prec)
  two <- mpfr(2,prec)
  multm <- mpfr(mult,prec)
  pmultm <- multm/two
  pt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  qt <- pt
  for (nn in 2:nmax) {
    pt[[nn]] <- mpfr2array(rep(nill,nn*nn), dim=c(nn,nn))
    rownames(pt[[nn]]) <- c(0:(nn-1))
    colnames(pt[[nn]]) <- c(1:nn)
    pt[[nn]][1,nn] <- pmultm^(nn-1)
    for (ff in 2:nn) {
      if (nn-ff+1<=ff-1) {
        f1 <- ff
        pt[[nn]][(1+1):(nn-f1+2),f1-1] <-
          pt[[nn]][2:(nn-f1+2),f1-1] +
          (pmultm^(f1-2)) * qt[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
      } # end if last part shortest
      if (nn-ff+1>ff-1) {
        f2 <- ff
        pt[[nn]][2:(nn-f2+2),f2-1] <-
          pt[[nn]][2:(nn-f2+2),f2-1] +
          (pmultm^(f1-2)) * qt[[nn-f2+1]][1:(nn-f2+1),f2-1]
        pt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] <-
          pt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] +
          (pmultm^(f1-2)) * pt[[nn-f2+1]][1:(nn-f2+1),f2:(nn-f2+1)]
      } # end if last part longest
    } # end for ff
    qt[[nn]] <- cumsumm(pt[[nn]])
    rownames(qt[[nn]]) <- c(0:(nn-1))
    colnames(qt[[nn]]) <- c(1:nn)
    if (printn) {
      print(nn)
      print(Sys.time())
    } # end optional timing information
  } # end for nn
  names(pt) <-paste("pt", 1:nmax, sep="")
  names(qt) <-paste("qt", 1:nmax, sep="")
  return(list(pt=pt,qt=qt))
} # end function crossrunsymm

# function for computing the simultaneous distribution for number of
# crossings and longest run in n independent binomial observations with
# probability p. Probabilities are multiplied by m^(n-1) where the 
# multiplier m by default is 2. The computations use the package Rmpfr 
#to enhance accuracy.
crossrunbin <- function(nmax=100, prob=.5, mult=2, prec=120, printn=FALSE) {
  nill <- mpfr(0,prec)
  one <- mpfr(1,prec)
  multm <- mpfr(mult,prec)
  pm <- mpfr(prob,prec)
  qm <- one - pm
  pmultm <- pm*multm
  qmultm <- qm*multm
  # conditioning of S= first value, pat: above 0, pbt: below 0
  # suffix t: probabilities times multm^(n-1). 
  # n=1:
  pat <- list(pt1=mpfr2array(one, dim=c(1,1)))
  pbt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  pt <- list(pt1=mpfr2array(one, dim=c(1,1))) 
  qat <- list(pt1=mpfr2array(one, dim=c(1,1))) 
  qbt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  qt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  for (nn in 2:nmax) {
    pat[[nn]] <- mpfr2array(rep(nill,nn*nn), dim=c(nn,nn))
    pbt[[nn]] <- mpfr2array(rep(nill,nn*nn), dim=c(nn,nn))
    rownames(pat[[nn]]) <- c(0:(nn-1))
    rownames(pbt[[nn]]) <- c(0:(nn-1))
    colnames(pat[[nn]]) <- c(1:nn)
    colnames(pbt[[nn]]) <- c(1:nn)
    pat[[nn]][1,nn] <- (pmultm^(nn-1)) # from cond on no crossing
    pbt[[nn]][1,nn] <- (qmultm^(nn-1)) # from cond on no crossing
    for (ff in 2:nn) { # from cond on first crossing at ff
      if (nn-ff+1<=ff-1) { # if last part shortest:
        f1 <- ff # unnecessary, but makes code checking easier
        pat[[nn]][2:(nn-f1+2),f1-1] <-
          pat[[nn]][2:(nn-f1+2),f1-1] +
          (pmultm^(f1-2)) * qmultm * qbt[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
        pbt[[nn]][2:(nn-f1+2),f1-1] <-
          pbt[[nn]][2:(nn-f1+2),f1-1] +
          (qmultm^(f1-2)) * pmultm * qat[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
      } # end if last part shortest
      if (nn-ff+1>ff-1) {# if last part longest
        f2 <- ff # unnecessary, but makes code checking easier
        pat[[nn]][2:(nn-f2+2),f2-1] <-
          pat[[nn]][2:(nn-f2+2),f2-1] +
          (pmultm^(f2-2)) * qmultm * qbt[[nn-f2+1]][1:(nn-f2+1),f2-1]
        pat[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] <-
          pat[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] +
          (pmultm^(f2-2)) * qmultm * pbt[[nn-f2+1]][1:(nn-f2+1),f2:(nn-f2+1)]
        pbt[[nn]][2:(nn-f2+2),f2-1] <-
          pbt[[nn]][2:(nn-f2+2),f2-1] +
          (qmultm^(f2-2)) * pmultm * qat[[nn-f2+1]][1:(nn-f2+1),f2-1]
        pbt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] <-
          pbt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] +
          (qmultm^(f2-2)) * pmultm * pat[[nn-f2+1]][1:(nn-f2+1),f2:(nn-f2+1)]
      } # end if last part longest
    } # end for ff
    pt[[nn]] <- pm*pat[[nn]] + qm*pbt[[nn]]
    qat[[nn]] <- cumsumm(pat[[nn]])
    qbt[[nn]] <- cumsumm(pbt[[nn]])
    qt[[nn]] <- pm*qat[[nn]] + qm*qbt[[nn]]
    rownames(pt[[nn]]) <- c(0:(nn-1))
    colnames(pt[[nn]]) <- c(1:nn)
    rownames(qat[[nn]]) <- c(0:(nn-1))
    colnames(qat[[nn]]) <- c(1:nn)
    rownames(qbt[[nn]]) <- c(0:(nn-1))
    rownames(qat[[nn]]) <- c(0:(nn-1))
    colnames(qt[[nn]]) <- c(1:nn)
    colnames(qt[[nn]]) <- c(1:nn)
    if (printn) {
      print(nn)
      print(Sys.time())
    } # end optional timing information
  } # end for nn
  names(pat) <-paste("pat", 1:nmax, sep="")
  names(pbt) <-paste("pbt", 1:nmax, sep="")
  names(pt) <-paste("pt", 1:nmax, sep="")
  names(qat) <-paste("qat", 1:nmax, sep="")
  names(qbt) <-paste("qbt", 1:nmax, sep="")
  names(qt) <-paste("qt", 1:nmax, sep="")
  return(list(pat=pat,pbt=pbt,pt=pt,qat=qat,qbt=qbt,qt=qt))
} # end function crossrunbin

# wrapper for crossrunbin, success probability defined as pnorm(shift)
# for a shift in a standard normal variable:
crossrunshift <- function(nmax=100, shift=0, mult=2, prec=120, printn=FALSE) {
  prob <- pnorm(shift)
  return(crossrunbin(nmax=nmax, prob=prob, mult=mult, prec=prec, printn=printn))
} # end function crossrunshift

# function for computing the simultaneous distribution for number of
# crossings and longest run in n independent binomial observations with
# possibly varying probability p. Probabilities are multiplied by m^(n-1) 
# where the multiplier m by default is 2. The computations use the package 
# Rmpfr to enhance accuracy.
crossrunchange <- function(nmax=100, prob=rep(.5,100), mult=2, prec=120, printn=FALSE) {
  nill <- mpfr(0,prec)
  one <- mpfr(1,prec)
  multm <- mpfr(mult,prec)
  pm <- mpfr(prob,prec)
  qm <- one - pm
  pmultm <- pm*multm
  qmultm <- qm*multm
  # conditioning of S= first value, pat: above 0, pbt: below 0
  # suffix t: probabilities times multm^(n-1). 
  # n=1:
  pat <- list(pt1=mpfr2array(one, dim=c(1,1)))
  pbt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  pt <- list(pt1=mpfr2array(one, dim=c(1,1))) 
  qat <- list(pt1=mpfr2array(one, dim=c(1,1))) 
  qbt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  qt <- list(pt1=mpfr2array(one, dim=c(1,1)))
  for (nn in 2:nmax) {
    pat[[nn]] <- mpfr2array(rep(nill,nn*nn), dim=c(nn,nn))
    pbt[[nn]] <- mpfr2array(rep(nill,nn*nn), dim=c(nn,nn))
    rownames(pat[[nn]]) <- c(0:(nn-1))
    rownames(pbt[[nn]]) <- c(0:(nn-1))
    colnames(pat[[nn]]) <- c(1:nn)
    colnames(pbt[[nn]]) <- c(1:nn)
    pat[[nn]][1,nn] <- prod(pmultm[(nmax+2-nn):nmax]) # from cond on no crossing
    pbt[[nn]][1,nn] <- prod(qmultm[(nmax+2-nn):nmax]) # from cond on no crossing
    for (ff in 2:nn) { # from cond on first crossing at ff
      if (nn-ff+1<=ff-1) { # if last part shortest:
        f1 <- ff # unnecessary, but makes code checking easier
        if (f1==2) {
          prodmulta <- one
          prodmultb <- one
        } else {
          prodmulta <- prod(pmultm[(nmax+2-nn):(nmax-nn+f1-1)])
          prodmultb <- prod(qmultm[(nmax+2-nn):(nmax-nn+f1-1)])
        }
        pat[[nn]][2:(nn-f1+2),f1-1] <-
          pat[[nn]][2:(nn-f1+2),f1-1] +
          prodmulta * qmultm[nmax-nn+f1] * qbt[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
        pbt[[nn]][2:(nn-f1+2),f1-1] <-
          pbt[[nn]][2:(nn-f1+2),f1-1] +
          prodmultb * pmultm[nmax-nn+f1] * qat[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
      } # end if last part shortest
      if (nn-ff+1>ff-1) {# if last part longest
        f2 <- ff # unnecessary, but makes code checking easier
        if (f2==2) {
          prodmulta <- one
          prodmultb <- one
        } else {
          prodmulta <- prod(pmultm[(nmax+2-nn):(nmax-nn+f2-1)])
          prodmultb <- prod(qmultm[(nmax+2-nn):(nmax-nn+f2-1)])
        }
        pat[[nn]][2:(nn-f2+2),f2-1] <-
          pat[[nn]][2:(nn-f2+2),f2-1] +
          prodmulta * qmultm[nmax-nn+f2] * qbt[[nn-f2+1]][1:(nn-f2+1),f2-1]
        pat[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] <-
          pat[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] +
          prodmulta * qmultm[nmax-nn+f2] * pbt[[nn-f2+1]][1:(nn-f2+1),f2:(nn-f2+1)]
        pbt[[nn]][2:(nn-f2+2),f2-1] <-
          pbt[[nn]][2:(nn-f2+2),f2-1] +
          prodmultb * pmultm[nmax-nn+f2] * qat[[nn-f2+1]][1:(nn-f2+1),f2-1]
        pbt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] <-
          pbt[[nn]][2:(nn-f2+2),f2:(nn-f2+1)] +
          prodmultb * pmultm[nmax-nn+f2] * pat[[nn-f2+1]][1:(nn-f2+1),f2:(nn-f2+1)]
      } # end if last part longest
    } # end for ff
    pt[[nn]] <- pm[nmax-nn+1]*pat[[nn]] + qm[nmax-nn+1]*pbt[[nn]]
    qat[[nn]] <- cumsumm(pat[[nn]])
    qbt[[nn]] <- cumsumm(pbt[[nn]])
    qt[[nn]] <- pm[nmax-nn+1]*qat[[nn]] + qm[nmax-nn+1]*qbt[[nn]]
    rownames(pt[[nn]]) <- c(0:(nn-1))
    colnames(pt[[nn]]) <- c(1:nn)
    rownames(qat[[nn]]) <- c(0:(nn-1))
    colnames(qat[[nn]]) <- c(1:nn)
    rownames(qbt[[nn]]) <- c(0:(nn-1))
    rownames(qat[[nn]]) <- c(0:(nn-1))
    colnames(qt[[nn]]) <- c(1:nn)
    colnames(qt[[nn]]) <- c(1:nn)
    if (printn) {
      print(nn)
      print(Sys.time())
    } # end optional timing information
  } # end for nn
  names(pat) <-paste("pat", 1:nmax, sep="")
  names(pbt) <-paste("pbt", 1:nmax, sep="")
  names(pt) <-paste("pt", 1:nmax, sep="")
  names(qat) <-paste("qat", 1:nmax, sep="")
  names(qbt) <-paste("qbt", 1:nmax, sep="")
  names(qt) <-paste("qt", 1:nmax, sep="")
  return(list(pat=pat,pbt=pbt,pt=pt,qat=qat,qbt=qbt,qt=qt))
} # end function crossrunchange

# compute simultaneous distributions for n=1, ..., 100 in the symmetric case:
Sys.time()
cr100 <- crossrunsymm(100)
Sys.time() # about 3 minutes (R3.4.4, Rmpfr 0-7.0)


# checks of crossrunsymm: 
# shows the simultaneous distributions for n=1, ... 7 for 
# comparisons with the result of brute force calculations:
for (nn in 1:7) print(cr100$pt[[nn]]) 
# comparison with result of brute force calculations
# n=1: 1, n=2, reverse diagonal:
cr100$pt[[1]]
cr100$pt[[2]]
matrix(c(0,0,1, 0,2,0, 1,0,0), ncol=3) -
  matrix(as.numeric(cr100$pt[[3]]),ncol=3)
matrix(c(0,0,0,1, 0,1,3,0, 0,2,0,0, 1,0,0,0), ncol=4) -
  matrix(as.numeric(cr100$pt[[4]]),ncol=4)
matrix(c(0,0,0,0,1, 0,0,3,4,0, 0,2,3,0,0, 0,2,0,0,0, 1,0,0,0,0), ncol=5) -
  matrix(as.numeric(cr100$pt[[5]]),ncol=5)
matrix(c(0,0,0,0,0,1, 0,0,1,6,5,0,  0,1,6,4,0,0,
         0,2,3,0,0,0, 0,2,0,0,0,0, 1,0,0,0,0,0), ncol=6) -
  matrix(as.numeric(cr100$pt[[6]]),ncol=6)
matrix(c(0,0,0,0,0,0,1, 0,0,0,4,10,6,0,  0,0,6,12,5,0,0, 0,2,6,4,0,0,0,
         0,2,3,0,0,0,0, 0,2,0,0,0,0,0, 1,0,0,0,0,0,0), ncol=7) -
  matrix(as.numeric(cr100$pt[[7]]),ncol=7)
# ok

# check of marginal distributions of C from the simultaneous
# distribution compared with binomial formulas in the 
# symmetric case (printn=TRUE keeps track of the progress 
# in the time-consuming computation):
for (nn in 2:100) {
  print(nn)
  print(range(as.numeric(cumsumm(cr100$pt[[nn]])[-1,nn] - chooseMpfr.all(nn-1))))
}
# ok

# crossrunbin distributions, check with crossrunsymm for p=0.5:
crb100.0 <- crossrunbin(printn=TRUE)$pt
for (nn in 1:100) {
  print(nn)
  print(sum(abs(crb100.0[[nn]] - cr100$pt[[nn]])))
} # no differences, ok

# binomial probabilities .6, .7, .8, .9:
Sys.time()
crb100..6 <- crossrunbin(prob=.6, printn=TRUE)$pt
crb100..7 <- crossrunbin(prob=.7, printn=TRUE)$pt
crb100..8 <- crossrunbin(prob=.8, printn=TRUE)$pt
crb100..9 <- crossrunbin(prob=.9, printn=TRUE)$pt
Sys.time()


# comparison with explicit computations for low n:
exactbin <- function(n, p=.5, prec=120) {
  nill <- mpfr(0,prec)
  one <- mpfr(1,prec)
  two <- mpfr(2,prec)
  pm <- mpfr(p,prec)
  qm <- one - pm
  res <- mpfr2array(one, dim=c(1,1))
  if (n==2) res <- mpfr2array(c(nill,two*pm*qm,pm^2+qm^2,nill), dim=c(2,2))
  if (n==3) res <- mpfr2array(c(nill,nill,pm*qm,nill,two*pm*qm,
                                nill,pm^3+qm^3,nill,nill), dim=c(3,3))
  if (n==4) res <- mpfr2array(
    c(rep(nill,3),2*pm^2*qm^2,  nill,2*pm^2*qm^2,two*pm*qm*(1-pm*qm),nill,
      nill,two*pm*qm*(pm^2+qm^2),nill,nill,  pm^4+qm^4,rep(nill,3)), dim=c(4,4))
  if (n==5) res <- mpfr2array(
    c(rep(nill,4),pm^2*qm^2,   nill,nill,pm*qm*(1-pm*qm),4*pm^2*qm^2,nill,
      nill,2*pm^2*qm^2,pm*qm*(2*pm^3+pm*qm+2*qm^3),nill,nill,
      nill,2*pm*qm*(pm^3+qm^3),rep(nill,3),    pm^5+qm^5,rep(nill,4)), dim=c(5,5))
  if (n==6) res <- mpfr2array(
    c(rep(nill,5),2*pm^3*qm^3, 
      nill,nill,pm^4*qm^2+pm^2*qm^4,2*pm^4*qm^2+8*pm^3*qm^3+2*pm^2*qm^4,
        3*pm^4*qm^2+4*pm^3*qm^3+3*pm^2*qm^4,nill,
      nill,2*pm^3*qm^3,2*pm^5*qm+2*pm^4*qm^2+4*pm^3*qm^3+2*pm^2*qm^4+2*pm*qm^5,
        4*pm^4*qm^2+4*pm^2*qm^4,nill,nill,
      nill,2*pm^4*qm^2+2*pm^2*qm^4,2*pm^5*qm+pm^4*qm^2+pm^2*qm^4+2*pm*qm^5,
        rep(nill,3),
      nill,2*pm^5*qm+2*pm*qm^5,rep(nill,4),   pm^6+qm^6,rep(nill,5)), dim=c(6,6))
  return(res)
} # end function exactbin
# n=2:
sum(abs(crb100..6[[2]]/sum(crb100..6[[2]]) - exactbin(n=2, p=.6)))
sum(abs(crb100..7[[2]]/sum(crb100..7[[2]]) - exactbin(n=2, p=.7)))
sum(abs(crb100..8[[2]]/sum(crb100..8[[2]]) - exactbin(n=2, p=.8)))
sum(abs(crb100..9[[2]]/sum(crb100..9[[2]]) - exactbin(n=2, p=.9)))
# equal
# n=3:
sum(abs(crb100..6[[3]]/sum(crb100..6[[3]]) - exactbin(n=3, p=.6)))
sum(abs(crb100..7[[3]]/sum(crb100..7[[3]]) - exactbin(n=3, p=.7)))
sum(abs(crb100..8[[3]]/sum(crb100..8[[3]]) - exactbin(n=3, p=.8)))
sum(abs(crb100..9[[3]]/sum(crb100..9[[3]]) - exactbin(n=3, p=.9)))
# equal
# n=4:
sum(abs(crb100..6[[4]]/sum(crb100..6[[4]]) - exactbin(n=4, p=.6)))
sum(abs(crb100..7[[4]]/sum(crb100..7[[4]]) - exactbin(n=4, p=.7)))
sum(abs(crb100..8[[4]]/sum(crb100..8[[4]]) - exactbin(n=4, p=.8)))
sum(abs(crb100..9[[4]]/sum(crb100..9[[4]]) - exactbin(n=4, p=.9)))
# in 37. decimal
# n=5:
sum(abs(crb100..6[[5]]/sum(crb100..6[[5]]) - exactbin(n=5, p=.6)))
sum(abs(crb100..7[[5]]/sum(crb100..7[[5]]) - exactbin(n=5, p=.7)))
sum(abs(crb100..8[[5]]/sum(crb100..8[[5]]) - exactbin(n=5, p=.8)))
sum(abs(crb100..9[[5]]/sum(crb100..9[[5]]) - exactbin(n=5, p=.9)))
# in 16. decimal
# n=6:
sum(abs(crb100..6[[6]]/sum(crb100..6[[6]]) - exactbin(n=6, p=.6)))
sum(abs(crb100..7[[6]]/sum(crb100..7[[6]]) - exactbin(n=6, p=.7)))
sum(abs(crb100..8[[6]]/sum(crb100..8[[6]]) - exactbin(n=6, p=.8)))
sum(abs(crb100..9[[6]]/sum(crb100..9[[6]]) - exactbin(n=6, p=.9)))
# in 16. decimal

# check by simulations for n=100, random binomial draws
# by dichotomizing a standard normal variable, to use the
# same random variable for several binomial probabilities:
# auxiliary function:
clshift <- function(seri,shift=0,type=0) {
  rle.sh <- rle(seri+shift>0)$lengths
  if (type==0) res <- length(rle.sh) - 1
  if (type==1) res <- max(rle.sh)
  return(res)
  } # end function clshift
# main function:
simclbin <- function(nser=100, nsim=100000, probs=c(.5,.6,.7,.8,.9)) {
  nprob <- length(probs)
  shifts <- qnorm(probs)
  series <- data.frame(matrix(rnorm(nser*nsim), nrow=nsim))
  res <- data.frame(matrix(rep(NA,2*nsim*nprob),nrow=nsim))
  names(res) <- paste(rep(c("nc","lr"), nprob), rep(probs,rep(2,nprob)), sep="")
  for (shnr in 1:nprob) {
    res[,1+2*(shnr-1)] <- apply(series,1,clshift,shift=shifts[shnr],type=0)
    res[,2*shnr] <- apply(series,1,clshift,shift=shifts[shnr],type=1)
  }
  return(res)
} # end function simclbin
set.seed(83938487)
Sys.time()
cl100simbin <- simclbin() 
Sys.time()
# about 3 minutes
# expectation of C for n=100:
sum(cumsumm(cr100$pt[[100]])[,100]*(0:99))/(2^99)
mean(cl100simbin$nc0.5) # 49.500 vs 49.502, ok
sum(cumsumm(crb100..6[[100]])[,100]*(0:99))/sum(cumsumm(crb100..6[[100]])[,100])
mean(cl100simbin$nc0.6) # 47.52 vs 47.51, ok
sum(cumsumm(crb100..7[[100]])[,100]*(0:99))/sum(cumsumm(crb100..7[[100]])[,100])
mean(cl100simbin$nc0.7) # 41.580 vs 47.575, ok
sum(cumsumm(crb100..8[[100]])[,100]*(0:99))/sum(cumsumm(crb100..8[[100]])[,100])
mean(cl100simbin$nc0.8) # 31.68 vs 31.69, ok
sum(cumsumm(crb100..9[[100]])[,100]*(0:99))/sum(cumsumm(crb100..9[[100]])[,100])
mean(cl100simbin$nc0.9) # 17.82 vs 17.84, ok
# expectation of L for n=100:
sum(cumsummcol(cr100$pt[[100]])[100,]*(1:100))/(2^99)
mean(cl100simbin$lr0.5) # 6.98 vs 6.97, ok
sum(cumsummcol(crb100..6[[100]])[100,]*(1:100))/sum(cumsummcol(crb100..6[[100]])[100,])
mean(cl100simbin$lr0.6) # 7.989 vs 7.988, ok
sum(cumsummcol(crb100..7[[100]])[100,]*(1:100))/sum(cumsummcol(crb100..7[[100]])[100,])
mean(cl100simbin$lr0.7) # 10.69 vs 10.66, ok
sum(cumsummcol(crb100..8[[100]])[100,]*(1:100))/sum(cumsummcol(crb100..8[[100]])[100,])
mean(cl100simbin$lr0.8) # 15.56 vs 15.52, ok
sum(cumsummcol(crb100..9[[100]])[100,]*(1:100))/sum(cumsummcol(crb100..9[[100]])[100,])
mean(cl100simbin$lr0.9) # 26.93 vs 26.89, ok
# standard deviation of C for n=100:
sqrt(sum(cumsumm(cr100$pt[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(cr100$pt[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simbin$nc0.5) # 4.97 vs 4.96, ok
sqrt(sum(cumsumm(crb100..6[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crb100..6[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simbin$nc0.6) # 5.16 vs 5.15, ok
sqrt(sum(cumsumm(crb100..7[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crb100..7[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simbin$nc0.7) # 5.54 vs 5.53, ok
sqrt(sum(cumsumm(crb100..8[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crb100..8[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simbin$nc0.8) # 5.730 vs 5.726, ok
sqrt(sum(cumsumm(crb100..9[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crb100..9[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simbin$nc0.9) # 5.089 vs 5.085, ok
# standard deviation of L for n=100:
sqrt(sum(cumsummcol(cr100$pt[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(cr100$pt[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simbin$lr0.5) # 1.7926 vs 1.7930, ok
sqrt(sum(cumsummcol(crb100..6[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crb100..6[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simbin$lr0.6) # 2.321 vs 2.322, ok
sqrt(sum(cumsummcol(crb100..7[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crb100..7[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simbin$lr0.7) # 3.34 vs 3.33, ok
sqrt(sum(cumsummcol(crb100..8[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crb100..8[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simbin$lr0.8) # 5.129 vs 5.133, ok
sqrt(sum(cumsummcol(crb100..9[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crb100..9[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simbin$lr0.9) # 9.736 vs 9.727, ok
# expectation of C*L for n=100:
(matrix(0:99,nrow=1) %*% cr100$pt[[100]] %*% matrix(1:100,ncol=1))/(sum(cr100$pt[[100]]))
mean(cl100simbin$nc0.5*cl100simbin$lr0.5) # 341.43 vs 341.38, ok
(matrix(0:99,nrow=1) %*% crb100..6[[100]] %*% matrix(1:100,ncol=1))/(sum(crb100..6[[100]]))
mean(cl100simbin$nc0.6*cl100simbin$lr0.6) # 373.89 vs 373.83, ok
(matrix(0:99,nrow=1) %*% crb100..7[[100]] %*% matrix(1:100,ncol=1))/(sum(crb100..7[[100]]))
mean(cl100simbin$nc0.7*cl100simbin$lr0.7) # 434.4 vs 433.4, ok
(matrix(0:99,nrow=1) %*% crb100..8[[100]] %*% matrix(1:100,ncol=1))/(sum(crb100..8[[100]]))
mean(cl100simbin$nc0.8*cl100simbin$lr0.8) # 475.9 vs 474.8, ok
(matrix(0:99,nrow=1) %*% crb100..9[[100]] %*% matrix(1:100,ncol=1))/(sum(crb100..9[[100]]))
mean(cl100simbin$nc0.9*cl100simbin$lr0.9) # 448.56 vs 448.58, ok
# compare cdf of C with simulations, all plots indistinguishible:
plot(x=as.numeric(names(table(cl100simbin$nc0.5))), y=(cumsum(cumsumm(cr100$pt[[100]])[,100])/(2^99))[
  as.numeric(names(table(cl100simbin$nc0.5)))+1], type="l")
points(x=as.numeric(names(table(cl100simbin$nc0.5))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$nc0.5))/sum(table(cl100simbin$nc0.5)))
plot(x=as.numeric(names(table(cl100simbin$nc0.6))), y=(cumsum(cumsumm(crb100..6[[100]])[,100])/(2^99))[
  as.numeric(names(table(cl100simbin$nc0.6)))+1], type="l")
points(x=as.numeric(names(table(cl100simbin$nc0.6))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$nc0.6))/sum(table(cl100simbin$nc0.6)))
plot(x=as.numeric(names(table(cl100simbin$nc0.7))), y=(cumsum(cumsumm(crb100..7[[100]])[,100])/(2^99))[
  as.numeric(names(table(cl100simbin$nc0.7)))+1], type="l")
points(x=as.numeric(names(table(cl100simbin$nc0.7))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$nc0.7))/sum(table(cl100simbin$nc0.7)))
plot(x=as.numeric(names(table(cl100simbin$nc0.8))), y=(cumsum(cumsumm(crb100..8[[100]])[,100])/(2^99))[
  as.numeric(names(table(cl100simbin$nc0.8)))+1], type="l")
points(x=as.numeric(names(table(cl100simbin$nc0.8))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$nc0.8))/sum(table(cl100simbin$nc0.8)))
plot(x=as.numeric(names(table(cl100simbin$nc0.9))), y=(cumsum(cumsumm(crb100..9[[100]])[,100])/(2^99))[
  as.numeric(names(table(cl100simbin$nc0.9)))+1], type="l")
points(x=as.numeric(names(table(cl100simbin$nc0.9))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$nc0.9))/sum(table(cl100simbin$nc0.9)))
# compare cdf of L with simulations, all plots indistinguishible:
plot(x=as.numeric(names(table(cl100simbin$lr0.5))), type="l",
     y=as.numeric(cumsum(cumsummcol(cr100$pt[[100]])[100,])/
                    sum(cumsummcol(cr100$pt[[100]])[100,]))[as.numeric(names(table(cl100simbin$lr0.5)))])
points(x=as.numeric(names(table(cl100simbin$lr0.5))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$lr0.5))/sum(table(cl100simbin$lr0.5)))
plot(x=as.numeric(names(table(cl100simbin$lr0.6))), type="l",
     y=as.numeric(cumsum(cumsummcol(crb100..6[[100]])[100,])/
                    sum(cumsummcol(crb100..6[[100]])[100,]))[as.numeric(names(table(cl100simbin$lr0.6)))])
points(x=as.numeric(names(table(cl100simbin$lr0.6))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$lr0.6))/sum(table(cl100simbin$lr0.6)))
plot(x=as.numeric(names(table(cl100simbin$lr0.7))), type="l",
     y=as.numeric(cumsum(cumsummcol(crb100..7[[100]])[100,])/
                    sum(cumsummcol(crb100..7[[100]])[100,]))[as.numeric(names(table(cl100simbin$lr0.7)))])
points(x=as.numeric(names(table(cl100simbin$lr0.7))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$lr0.7))/sum(table(cl100simbin$lr0.7)))
plot(x=as.numeric(names(table(cl100simbin$lr0.8))), type="l",
     y=as.numeric(cumsum(cumsummcol(crb100..8[[100]])[100,])/
                    sum(cumsummcol(crb100..8[[100]])[100,]))[as.numeric(names(table(cl100simbin$lr0.8)))])
points(x=as.numeric(names(table(cl100simbin$lr0.8))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$lr0.8))/sum(table(cl100simbin$lr0.8)))
plot(x=as.numeric(names(table(cl100simbin$lr0.9))), type="l",
     y=as.numeric(cumsum(cumsummcol(crb100..9[[100]])[100,])/
                    sum(cumsummcol(crb100..9[[100]])[100,]))[as.numeric(names(table(cl100simbin$lr0.9)))])
points(x=as.numeric(names(table(cl100simbin$lr0.9))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simbin$lr0.9))/sum(table(cl100simbin$lr0.9)))

# check of shift distributions where 0 for shift=0, only for n=100:
max(abs(matrix(as.numeric(crb100..6[[100]]), ncol=100)*
          (matrix(as.numeric(cr100$pt[[100]]), ncol=100)==0))) # 0, ok
max(abs(matrix(as.numeric(crb100..7[[100]]), ncol=100)*
          (matrix(as.numeric(cr100$pt[[100]]), ncol=100)==0))) # 0, ok
max(abs(matrix(as.numeric(crb100..8[[100]]), ncol=100)*
          (matrix(as.numeric(cr100$pt[[100]]), ncol=100)==0))) # 0, ok
max(abs(matrix(as.numeric(crb100..9[[100]]), ncol=100)*
          (matrix(as.numeric(cr100$pt[[100]]), ncol=100)==0))) # 0, ok

# check of crossrunchange for constant probability:
crc100..6 <- crossrunchange(nmax=100,prob=rep(.6,100), printn=TRUE)$pt
for (nn in 1:100) {
  print(nn)
  print(sum(abs(crb100..6[[nn]] - crc100..6[[nn]])))
} # differences larger for higher n, in 7. decimal for n=100
# tried with 250 bits precision, not better (in fact 6. decimal):

# check of crossrunchange for two probabilities:
crc100.2p <- crossrunchange(nmax=100,
                            prob=c(rep(.5,40),rep(.6,60)), printn=TRUE)$pt
for (nn in 1:61) {
  print(nn)
  print(sum(abs(crb100..6[[nn]] - crc100.2p[[nn]])))
} # difference in <= 19, decimal for n<=60, large difference for n=61, ok

# example with changing binomial probabilities:
prob.change <- 1/(1+(1/9)^((0:99)/99)) # logistic change from 0.5 to 0.9, n=100
plot(x=1:100, y=prob.change, type="l")
crc100.probchange <- crossrunchange(nmax=100, prob=prob.change, printn=TRUE)$pt

# compare with simulations:
simclchange <- function(nsim=100000, probs=prob.change) {
  nser <- length(probs)
  shifts <- qnorm(probs)
  series <- data.frame(matrix(rnorm(nser*nsim), nrow=nsim))
  res <- data.frame(matrix(rep(NA,2*nsim),nrow=nsim))
  names(res) <- c("nc","lr")
  res[,1] <- apply(series,1,clshift,shift=shifts,type=0)
  res[,2] <- apply(series,1,clshift,shift=shifts,type=1)
  return(res)
} # end function simclchange
set.seed(83938487)
Sys.time()
cl100simchange <- simclchange() 
Sys.time() # about 20 seconds
# expectation of C for n=100:
sum(cumsumm(crc100.probchange[[100]])[,100]*(0:99))/(2^99)
mean(cl100simchange$nc) # 36.047 vs 36.053, ok
# expectation of L for n=100:
sum(cumsummcol(crc100.probchange[[100]])[100,]*(1:100))/(2^99)
mean(cl100simchange$lr) # 15.61 vs 15.58, ok
# standard deviation of C for n=100:
sqrt(sum(cumsumm(crc100.probchange[[100]])[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crc100.probchange[[100]])[,100]*(0:99))/(2^99))^2)
sd(cl100simchange$nc) # 5.43 vs 5.44, ok
# standard deviation of L for n=100:
sqrt(sum(cumsummcol(crc100.probchange[[100]])[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crc100.probchange[[100]])[100,]*(1:100))/(2^99))^2)
sd(cl100simchange$lr) # 5.68 vs 5.66, ok
# expectation of C*L for n=100:
(matrix(0:99,nrow=1) %*% crc100.probchange[[100]] %*% 
    matrix(1:100,ncol=1))/(sum(crc100.probchange[[100]]))
mean(cl100simchange$nc*cl100simchange$lr) # 546.8 vs 545.8, ok
# compare cdf of C with simulations (indistinguishable):
plot(x=as.numeric(names(table(cl100simchange$nc))), 
     y=(cumsum(cumsumm(crc100.probchange[[100]])[,100])/(2^99))[
       as.numeric(names(table(cl100simchange$nc)))+1], type="l")
points(x=as.numeric(names(table(cl100simchange$nc))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simchange$nc))/sum(table(cl100simchange$nc)))
# compare cdf of L with simulations:
plot(x=as.numeric(names(table(cl100simchange$lr))), type="l",
     y=as.numeric(cumsum(cumsummcol(crc100.probchange[[100]])[100,])/
                    sum(cumsummcol(cr100$pt[[100]])[100,]))[as.numeric(names(table(cl100simchange$lr)))])
points(x=as.numeric(names(table(cl100simchange$lr))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100simchange$lr))/sum(table(cl100simchange$lr)))

# plot to include in article for n=24, f=16 (case 1):
cr24matr <- matrix(as.numeric(cr100$pt[[24]]), ncol=24)
rownames(cr24matr) <- 0:23
colnames(cr24matr) <- 1:24
cr24matr
# runs plot for n=24:
set.seed(83938487)
values <- c(runif(15,.5,1),runif(1,-1,-.5),2*rbinom(8,1,.5)-1+runif(8,-.3,.3))
xright <- 70
ytop <- 7
ybottom <- -1.2
ytimes <- -1.4
cextimes <- .8
cexc <- .5
x24 <- 30
xincr <- 1.5
y24 <- ytop -.5
yincr <- .2
x9 <- 5
y9 <- 6
par(mar=c(bottom=3,left=3,top=0,right=0)+.1)
plot(x=1:24, y=values, axes=FALSE, type="o", pch=19, xlab="", ylab="", 
     xlim=c(1,xright), ylim=c(-2,ytop))
lines(x=c(1,24), y=c(0,0), col="blue")
lines(x=c(1,24), y=c(ybottom,ybottom))
points(x=1, y=values[1], col="red", pch=19)
points(x=15:16, y=values[15:16], col="red", pch=19, type="o")
points(x=17:24, y=values[17:24], col="blue", pch=19, type="o")
text(x=c(1,16,24),y=ytimes, labels=c(1,16,24), cex=cextimes)
text(x=12,y=ytimes-.4, labels="Time", cex=cextimes)
lines(x=c(1,1), y=c(ybottom,max(values)), lty="dotted")
lines(x=c(16,16), y=c(ybottom,max(values)), lty="dotted")
lines(x=c(24,24), y=c(ybottom,values[24]), lty="dotted")
text(x=1,y=max(values), labels="s=1", col="red", pos=3)
text(x=16,y=max(values), labels="f=16", col="red", pos=3)
for (c1 in 0:23) {
  print(lines(x=c(x24,x24+24*xincr), 
              y=c(y24-c1*yincr,y24-c1*yincr), lty="dotted"))
}
for (c1 in c(0,10:23))
  print(text(x=x24+xincr,y=y24-c1*yincr-yincr/2, labels=c1, cex=cexc, pos=2))
for (c1 in 1:9)
  print(text(x=x24+xincr,y=y24-c1*yincr-yincr/2, labels=c1, 
             cex=cexc, pos=2, col="red"))
lines(x=c(x24,x24+24*xincr), y=c(y24-24*yincr,y24-24*yincr), lty="dotted")
text(x=x24-xincr,y=y24-12*yincr, labels="C", pos=2)
for (l1 in 1:24) {
  print(lines(x=c(x24+l1*xincr,x24+l1*xincr), 
              y=c(y24-24*yincr,y24), lty="dotted"))
}
for (l1 in c(1:7,9:12)*2)
  print(text(x=x24+l1*xincr-xincr/2,y=y24-yincr, labels=l1, cex=cexc, pos=3))
text(x=x24+16*xincr-xincr/2,y=y24-yincr, labels=16, 
     cex=cexc, pos=3, col="red")
for (l1 in c(c(1:12)*2-1))
  print(text(x=x24+l1*xincr-xincr/2,
             y=y24-24*yincr+yincr, labels=l1, cex=cexc, pos=1))
text(x=x24+12*xincr,y=y24+yincr, labels="L", pos=3)
text(x=x24+12*xincr,y=y24-24*yincr-yincr/2, labels="L", pos=1)
rect(xleft=x24+15*xincr, xright=x24+16*xincr, border=NA,
     ybottom=y24-10*yincr, ytop=y24-1*yincr, col="red")
for (c1 in 0:9) {
  print(lines(x=c(x9,x9+9*xincr), 
              y=c(y9-c1*yincr,y9-c1*yincr), lty="dotted"))
}
for (c1 in 0:8) {
  print(text(x=x9+xincr,y=y9-c1*yincr-yincr/2, 
             labels=c1, cex=cexc, pos=2))
}
for (l1 in 0:9) {
  print(lines(x=c(x9+l1*xincr,x9+l1*xincr), 
              y=c(y9-9*yincr,y9), lty="dotted"))
}
for (l1 in 1:9) {
  print(text(x=x9+l1*xincr-xincr/2,
             y=y9-yincr, labels=l1, cex=cexc, pos=3))
}
text(x=x9+5*xincr,y=y9+yincr, labels="Last 9 (time 16-24)",  
     pos=3, cex=.7)
arrows(x0=x9, x1=x9+9*xincr, y0=y9-11*yincr,y1=y9-11*yincr,
       col="red")
text(x=x9-xincr,y=y9-14*yincr, labels="Cumutative sums\nin each row",  
     pos=4, cex=.6, col="red")
rect(xleft=x9+8*xincr, xright=x9+9*xincr, border=NA,
     ybottom=y9-9*yincr, ytop=y9, col="red")
arrows(x0=x9+9*xincr, x1=x24+15*xincr, col="red",
       y0=y9-4*yincr,y1=y9-4*yincr)
par(mar=c(bottom=5,left=4,top=4,right=2)+.1)

# plot to include in article for n=24, f=9 (case 2):
cr24matr
# runs plot for n=24, f=9:
set.seed(83938487)
values2 <- c(runif(8,.5,1),runif(1,-1,-.5),2*rbinom(15,1,.5)-1+runif(15,-.3,.3))
xright <- 70
ytop <- 7
ybottom2 <- -1.4
ytimes2 <- -1.7
cextimes2 <- .7
cexc <- .5
x24 <- 35.5
xincr <- 1.5
y24 <- ytop -.5
yincr <- .2
x16 <- 1
y16 <- 6
par(mar=c(bottom=3,left=3,top=0,right=0)+.1)
plot(x=1:24, y=values2, axes=FALSE, type="o", pch=19, xlab="", ylab="", 
     xlim=c(1,xright), ylim=c(-2.2,ytop))
points(x=10:24, y=values2[10:24], type="o", pch=19, col="blue")
lines(x=c(1,24), y=c(0,0), col="blue")
lines(x=c(1,24), y=c(ybottom2,ybottom2))
points(x=1, y=values[1], col="red", pch=19)
points(x=8:9, y=values2[8:9], col="red", pch=19, type="o")
text(x=c(1,9,24),y=ytimes2, labels=c(1,9,24), cex=cextimes2)
text(x=12,y=ytimes2-.5, labels="Time", cex=cextimes)
lines(x=c(1,1), y=c(ybottom2,max(values2)+.15), lty="dotted")
lines(x=c(9,9), y=c(ybottom2,max(values2)+.4), lty="dotted")
lines(x=c(24,24), y=c(ybottom2,values2[24]), lty="dotted")
text(x=3,y=max(values2[1:8]), labels="s=1", 
     col="red", cex=cextimes2 ,pos=3)
text(x=9,y=max(values2), labels="f=9", 
     col="red", cex=cextimes2, pos=3)
rect(xleft=x24+7*xincr, xright=x24+8*xincr, border=NA,
     ybottom=y24-17*yincr, ytop=y24-1*yincr, col="red")
rect(xleft=x24+8*xincr, xright=x24+16*xincr, border=NA,
     ybottom=y24-17*yincr, ytop=y24-1*yincr, col="blue")
for (c1 in 0:23) {
  print(lines(x=c(x24,x24+24*xincr), 
              y=c(y24-c1*yincr,y24-c1*yincr), lty="dotted"))
}
for (c1 in c(0,17:23))
  print(text(x=x24+xincr,y=y24-c1*yincr-yincr/2, labels=c1, cex=cexc, pos=2))
for (c1 in 1:16)
  print(text(x=x24+xincr,y=y24-c1*yincr-yincr/2, labels=c1, 
             cex=cexc, pos=2, col="red"))
lines(x=c(x24,x24+24*xincr), y=c(y24-24*yincr,y24-24*yincr), lty="dotted")
text(x=x24-xincr,y=y24-yincr, cex=.6, labels="C", pos=2)
text(x=x24-xincr,y=y24-22*yincr, cex=.6, labels="C", pos=2)
for (l1 in 1:24) {
  print(lines(x=c(x24+l1*xincr,x24+l1*xincr), 
              y=c(y24-24*yincr,y24), lty="dotted"))
}
for (l1 in c(1:3,9:12)*2)
  print(text(x=x24+l1*xincr-xincr/2,y=y24-yincr, labels=l1, cex=cexc, pos=3))
for (l1 in c(5:8)*2)
  print(text(x=x24+l1*xincr-xincr/2,y=y24-yincr, 
             labels=l1, cex=cexc, pos=3, col="blue"))
text(x=x24+8*xincr-xincr/2,y=y24-yincr, labels=8, 
     cex=cexc, pos=3, col="red")
for (l1 in c(c(1:4,9:12)*2-1))
  print(text(x=x24+l1*xincr-xincr/2,
             y=y24-24*yincr+yincr, labels=l1, cex=cexc, pos=1))
for (l1 in c(c(5:8)*2-1))
  print(text(x=x24+l1*xincr-xincr/2, col="blue",
             y=y24-24*yincr+yincr, labels=l1, cex=cexc, pos=1))
text(x=x24+12*xincr,y=y24+yincr, labels="L", pos=3)
text(x=x24+12*xincr,y=y24-24*yincr-yincr/2, labels="L", pos=1)
rect(xleft=x16+7*xincr, xright=x16+8*xincr, border=NA,
     ybottom=y16, ytop=y16-16*yincr, col="red")
rect(xleft=x16+8*xincr, xright=x16+16*xincr, border=NA,
     ybottom=y16-16*yincr, ytop=y16, col="blue")
for (c1 in 0:15) {
  print(lines(x=c(x16,x16+16*xincr), 
              y=c(y16-c1*yincr,y16-c1*yincr), lty="dotted"))
}
for (c1 in (c(1:8)*2-1)) {
  print(text(x=x16+xincr+xincr/2,y=y16-c1*yincr+yincr/2, 
             labels=c1-1, cex=cexc, pos=2))
}
for (c1 in (c(1:8)*2)) {
  print(text(x=x16+15*xincr-xincr/2,y=y16-c1*yincr+yincr/2, 
             labels=c1-1, cex=cexc, pos=4))
}
for (l1 in 0:15) {
  print(lines(x=c(x16+l1*xincr,x16+l1*xincr), 
              y=c(y16-16*yincr,y16), lty="dotted"))
}
for (l1 in c(1:3)*2) {
  print(text(x=x16+l1*xincr-xincr/2,
             y=y16-yincr, labels=l1, cex=cexc, pos=3))
}
for (l1 in (c(1:4)*2-1)) {
  print(text(x=x16+l1*xincr-xincr/2,
             y=y16-15*yincr+yincr/2, labels=l1, cex=cexc, pos=1))
}
print(text(x=x16+8*xincr-xincr/2, col="red",
           y=y16-yincr, labels=8, cex=cexc, pos=3))
for (l1 in (c(5:8)*2)) {
  print(text(x=x16+l1*xincr-xincr/2, col="blue",
             y=y16-yincr, labels=l1, cex=cexc, pos=3))
}
for (l1 in c(5:8)*2-1) {
  print(text(x=x16+l1*xincr-xincr/2, col="blue",
             y=y16-20*yincr, labels=l1, cex=cexc, pos=3))
}
text(x=x16+5*xincr,y=y9+yincr, labels="Last 16 (time 9-24)",  
     pos=3, cex=.5)
arrows(x0=x16, x1=x16+7*xincr, y0=y16-7*yincr,y1=y16-7*yincr,
       col="red", lwd=2)
par(mar=c(bottom=5,left=4,top=4,right=2)+.1)

# joint probabilities for n=15, symmetric case:
symm15 <- matrix(as.numeric(cr100$pt[[15]]),ncol=15)
rownames(symm15) <- paste0("c=",0:14)
colnames(symm15) <- paste0("l=",1:15)
symm15
sum(symm15)
sum(symm15) - 2^14

bin15.7 <- matrix(as.numeric(crb100..6[[15]]),ncol=15)
rownames(bin15.7) <- paste0("c=",0:14)
colnames(bin15.7) <- paste0("l=",1:15)
round(bin15.7,1)
sum(bin15.7) - 2^14

# commands for saving some important objects:
saveRDS(object=cr100$pt, file="crsymm.Rdata")
saveRDS(object=crb100..6, file="cr.6.Rdata")
saveRDS(object=crb100..7, file="cr.7.Rdata")
saveRDS(object=crb100..8, file="cr.8.Rdata")
saveRDS(object=crb100..9, file="cr.9.Rdata")

# tries getting the data from GitHub:
crsymm <- repmis::source_data(rdata=TRUE,
                              url="https://github.com/ToreWentzel-Larsen/crossrun/blob/master/crsymm.Rdata")
crsymm <- readr::read_rds(path="https://github.com/ToreWentzel-Larsen/crossrun/blob/master/crsymm.Rdata")

# vise matrisene for gitt n:
# det symmetriske tilfellet:
asNumeric(cr100$pt[[20]])
# har det også for p=0,6, 0,7, 0,8, 0,9
round(asNumeric(crb100..6$pt20),1) # p=0,6, vises med en desimal

# objekt for 2 sd:
crb100.2sd <- crossrunbin(prob=pnorm(2), printn=TRUE)$pt
# check by simulations:
set.seed(83938487)
Sys.time() # less than 1 minute
cl100.2sdsimbin <- simclbin(probs=pnorm(2)) 
Sys.time()
# about 3 minutes
# expectation of C for n=100:
sum(cumsumm(crb100.2sd$pt100)[,100]*(0:99))/(2^99)
mean(cl100.2sdsimbin$nc2) # 4.40 vs 4.41, ok
# expectation of L for n=100:
sum(cumsummcol(crb100.2sd$pt100)[100,]*(1:100))/(2^99)
mean(cl100.2sdsimbin$lr2) # 62.10 vs 62.0, ok
# standard deviation of C for n=100:
sqrt(sum(cumsumm(crb100.2sd$pt100)[,100]*((0:99)^2))/(2^99) - 
       (sum(cumsumm(crb100.2sd$pt100)[,100]*(0:99))/(2^99))^2)
sd(cl100.2sdsimbin$nc2) # 2.859 vs 2.865, ok
# standard deviation of L for n=100:
sqrt(sum(cumsummcol(crb100.2sd$pt100)[100,]*((1:100)^2))/(2^99) - 
       (sum(cumsummcol(crb100.2sd$pt100)[100,]*(1:100))/(2^99))^2)
sd(cl100.2sdsimbin$lr2) # 21.17 vs 21.20, ok
# expectation of C*L for n=100:
(matrix(0:99,nrow=1) %*% crb100.2sd$pt100 %*% matrix(1:100,ncol=1))/(sum(crb100.2sd$pt100))
mean(cl100.2sdsimbin$nc2*cl100.2sdsimbin$lr2) # 227.2 vs 227.0, ok
(matrix(0:99,nrow=1) %*% crb100.2sd$pt100 %*% matrix(1:100,ncol=1))/(sum(crb100.2sd$pt100))
# compare cdf of C with simulations, all plots indistinguishible:
plot(x=as.numeric(names(table(cl100.2sdsimbin$nc2))), y=(cumsum(cumsumm(crb100.2sd$pt100)[,100])/(2^99))[
  as.numeric(names(table(cl100.2sdsimbin$nc2)))+1], type="l")
points(x=as.numeric(names(table(cl100.2sdsimbin$nc2))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100.2sdsimbin$nc2))/sum(table(cl100.2sdsimbin$nc2)))
# compare cdf of L with simulations, all plots indistinguishible:
plot(x=as.numeric(names(table(cl100.2sdsimbin$lr2))), type="l",
     y=as.numeric(cumsum(cumsummcol(crb100.2sd$pt100)[100,])/
                    sum(cumsummcol(crb100.2sd$pt100)[100,]))[as.numeric(names(table(cl100.2sdsimbin$lr2)))])
points(x=as.numeric(names(table(cl100.2sdsimbin$lr2))), type="l", col="red", lty="dotted",
       y=cumsum(table(cl100.2sdsimbin$lr2))/sum(table(cl100.2sdsimbin$lr2)))







