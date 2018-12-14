library(Rmpfr)
library(crossrun)

crossrunauto <- function(nmax = 100, prob = 0.5, changeprob = 0.5,
                         mult = 2, prec = 120, printn = FALSE) {
  nill <- Rmpfr::mpfr(0, prec)
  one <- Rmpfr::mpfr(1, prec)
  half <- Rmpfr::mpfr(0.5, prec)
  multm <- Rmpfr::mpfr(mult, prec)
  pm <- Rmpfr::mpfr(prob, prec)
  qm <- one - pm
  changeprobm <- Rmpfr::mpfr(changeprob, prec)
  if (pm>=0.5) {
    upprobm <- changeprobm
    downprobm <- (one-pm)*upprobm/pm
  } else {
    downprobm <- changeprobm
    upprob <- pm*upprobm/(one-pm)
  }
  corrm <- one-upprobm/pm
  pmultm <- pm * multm
  qmultm <- qm * multm
  umultm <- upprobm * multm # multiplied probability for going up
  dmultm <- downprobm * multm  # multiplied probability for going down
  topmultm <- (one-downprobm) * mult  # multiplied probability for staying at top
  botmultm <- (one-upprobm) * mult  # multiplied probability for staying at bottom
  # conditioning of S= first value, pat: above 0, pbt: below 0 suffix
  # t: probabilities times multm^(n-1).  n=1:
  pat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  pbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  pt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  for (nn in 2:nmax) {
    pat[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
    pbt[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
    rownames(pat[[nn]]) <- c(0:(nn - 1))
    rownames(pbt[[nn]]) <- c(0:(nn - 1))
    colnames(pat[[nn]]) <- c(1:nn)
    colnames(pbt[[nn]]) <- c(1:nn)
    pat[[nn]][1, nn] <- (topmultm^(nn - 1))  # from cond on no crossing
    pbt[[nn]][1, nn] <- (botmultm^(nn - 1))  # from cond on no crossing
    for (ff in 2:nn) {
      # from cond on first crossing at ff if last part shortest:
      if (nn - ff + 1 <= ff - 1)
      {
        f1 <- ff  # unnecessary, but makes code checking easier
        pat[[nn]][2:(nn-f1+2), f1-1] <- pat[[nn]][2:(nn-f1+2),f1-1] +
          (topmultm^(f1-2))*dmultm*qbt[[nn-f1+1]][1:(nn-f1+1), nn-f1+1]
        pbt[[nn]][2:(nn-f1+2), f1-1] <- pbt[[nn]][2:(nn-f1+2), f1-1] +
          (botmultm^(f1-2))*umultm*qat[[nn-f1+1]][1:(nn-f1+1), nn-f1+1]
      }  # end if last part shortest
      if (nn - ff + 1 > ff - 1)
      {
        # if last part longest
        f2 <- ff  # unnecessary, but makes code checking easier
        pat[[nn]][2:(nn-f2+2), f2-1] <- pat[[nn]][2:(nn-f2+2), f2-1] +
          (topmultm^(f2-2))*dmultm*qbt[[nn-f2+1]][1:(nn-f2+1), f2-1]
        pat[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] <- pat[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] +
          (topmultm^(f2-2))*dmultm*pbt[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
        pbt[[nn]][2:(nn-f2+2), f2-1] <- pbt[[nn]][2:(nn-f2+2), f2-1] +
          (botmultm^(f2-2))*umultm*qat[[nn-f2+1]][1:(nn-f2+1), f2-1]
        pbt[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] <- pbt[[nn]][2:(nn-f2+2), f2:(nn-f2+1)] +
          (botmultm^(f2-2))*umultm*pat[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
      }  # end if last part longest
    }  # end for ff
    pt[[nn]] <- pm * pat[[nn]] + qm * pbt[[nn]]
    qat[[nn]] <- cumsumm(pat[[nn]])
    qbt[[nn]] <- cumsumm(pbt[[nn]])
    qt[[nn]] <- pm * qat[[nn]] + qm * qbt[[nn]]
    rownames(pt[[nn]]) <- c(0:(nn - 1))
    colnames(pt[[nn]]) <- c(1:nn)
    rownames(qat[[nn]]) <- c(0:(nn - 1))
    colnames(qat[[nn]]) <- c(1:nn)
    rownames(qbt[[nn]]) <- c(0:(nn - 1))
    rownames(qat[[nn]]) <- c(0:(nn - 1))
    colnames(qt[[nn]]) <- c(1:nn)
    colnames(qt[[nn]]) <- c(1:nn)
    if (printn)
    {
      print(nn)
      print(Sys.time())
    }  # end optional timing information
  }  # end for nn
  names(pat) <- paste("pat", 1:nmax, sep = "")
  names(pbt) <- paste("pbt", 1:nmax, sep = "")
  names(pt) <- paste("pt", 1:nmax, sep = "")
  names(qat) <- paste("qat", 1:nmax, sep = "")
  names(qbt) <- paste("qbt", 1:nmax, sep = "")
  names(qt) <- paste("qt", 1:nmax, sep = "")
  return(list(pat = pat, pbt = pbt, pt = pt, qat = qat, qbt = qbt,
              qt = qt))
} # end function crossrunauto

# check, symmetric case, independence:
cr20 <- crossrunbin(nmax=20, prob=0.5, printn=TRUE)
cra20 <- crossrunauto(nmax=20, prob=0.5, changeprob=.5,
                      printn=TRUE)
asNumeric(cr20$pt[[20]])
asNumeric(cra20$pt[[20]])
asNumeric(cr20$pt[[20]]) - asNumeric(cra20$pt[[20]])
# the same
# check, p=0.6, independence:
cr20.6 <- crossrunbin(nmax=20, prob=0.6, printn=TRUE)
cra20.6 <- crossrunauto(nmax=20, prob=0.6, changeprob=.6,
                        printn=TRUE)
asNumeric(cr20.6$pt[[20]])
asNumeric(cr20.6$pt[[20]])
asNumeric(cr20.6$pt[[20]]) - asNumeric(cra20.6$pt[[20]])
# the same
# check, symmetric case, some dependence:
cra20.u.4 <- crossrunauto(nmax=20, prob=0.5, changeprob=.4,
                          printn=TRUE)
asNumeric(cr20$pt[[20]])
round(asNumeric(cra20.u.4$pt[[20]]),1)
# check, p=0.6, some dependence:
cra20.6.u.5 <- crossrunauto(nmax=20, prob=0.6, changeprob=.5,
                          printn=TRUE)
round(asNumeric(cr20.6$pt[[20]]),1)
round(asNumeric(cra20.6.u.5$pt[[20]]),1)
# delete small series:
rm(cr20,cr20.6,cra20,cra20.6,cra20.6.u.5,cra20.u.4)

# p=0.6, some dependence, nmax=100:
cra100.6.u.5 <- crossrunauto(nmax=100, prob=0.6, changeprob=.5,
                            printn=TRUE)$pt
save(cra100.6.u.5, file="./cra100.6.u.5.Rdata")

# simulation with autocorrelation, auxiliary function:
clf <- function(seri, type=0) {
  rleser <- rle(seri)$lengths
  if (type==0) res <- length(rleser) - 1
  if (type==1) res <- max(rleser)
  return(res)
} # end auxiliary function clf
# function for simulation with autocorrelation:
simclauto <- function (nser = 100, nsim = 1e+05, prob=0.5, changeprob=0.4) {
  if (prob>=0.5) {
    upprob <- changeprob
    downprob <- (1-prob)*upprob/prob
  } else {
    downprob <- changeprobm
    upprob <- prob*upprob/(1-prob)
  } # end setting downprob and upprob
  un <- data.frame(matrix(stats::runif(nser*nsim), nrow=nsim))
  series <- data.frame(matrix(rep(0,nser*nsim), nrow=nsim))
  series[,1] <- as.numeric(un[,1]<=prob)
  for (nr in 2:nser) series[,nr] <- series[,nr-1]*(un[,nr]>=downprob) +
    (1-series[,nr-1])*(un[,nr]<=upprob)
  res <- data.frame(matrix(rep(NA, 2*nsim), nrow = nsim))
  names(res) <- c("nc", "lr")
  res[,1] <- apply(series, 1, clf, type=0)
  res[,2] <- apply(series, 1, clf, type=1)
  return(res)
} # end function simclauto

# simulations for  p=0.6, some dependence, n=100:
set.seed(83938487)
sim.6.u.5 <- simclauto(prob=0.6, changeprob=0.5)
# joint distribution for n=100, times representation:
joint100.6.u.5 <- cra100.6.u.5[[100]]

# expectation of C:
sum(cumsumm(joint100.6.u.5)[,100]*(0:99))/sum(joint100.6.u.5)
mean(sim.6.u.5$nc) # 39.60 vs 39.62, ok
# expectation of L:
sum(cumsummcol(joint100.6.u.5)[100,]*(1:100))/sum(joint100.6.u.5)
mean(sim.6.u.5$lr) # 9.551 vs 9.543, ok
# standard deviation of C:
sqrt(sum(cumsumm(joint100.6.u.5)[,100]*(0:99)^2)/sum(joint100.6.u.5) - 
       (sum(cumsumm(joint100.6.u.5)[,100]*(0:99))/sum(joint100.6.u.5))^2)
# standard deviation of L:
sqrt(sum(cumsummcol(joint100.6.u.5)[100,]*(1:100)^2)/sum(joint100.6.u.5) - 
       (sum(cumsummcol(joint100.6.u.5)[100,]*(1:100))/sum(joint100.6.u.5))^2)
sd(sim.6.u.5$lr) # 2.858 vs 2.860, ok
# expectation of C*L:
(matrix(0:99,nrow=1) %*% joint100.6.u.5 %*% matrix(1:100,ncol=1))/sum(joint100.6.u.5)
mean(sim.6.u.5$nc*sim.6.u.5$lr) # 371.12 vs 371.0, ok
# cdf of C compared with simulations:
plot(x=as.numeric(names(table(sim.6.u.5$nc))), 
     y=(cumsum(cumsumm(joint100.6.u.5)[,100])/sum(joint100.6.u.5))[
       as.numeric(names(table(sim.6.u.5$nc)))+1], type="l", 
     xlab="C", ylab="CDF C")
points(x=as.numeric(names(table(sim.6.u.5$nc))), type="l", col="red", lty="dotted",
       y=cumsum(table(sim.6.u.5$nc))/sum(table(sim.6.u.5$nc)))
# cdf of L compared with simulations:
plot(x=as.numeric(names(table(sim.6.u.5$lr))), 
     y=(cumsum(cumsummcol(joint100.6.u.5)[100,])/sum(joint100.6.u.5))[
       as.numeric(names(table(sim.6.u.5$lr)))], type="l", 
     xlab="L", ylab="CDF L")
points(x=as.numeric(names(table(sim.6.u.5$lr))), type="l", col="red", lty="dotted",
       y=cumsum(table(sim.6.u.5$lr))/sum(table(sim.6.u.5$lr)))



