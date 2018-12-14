library(crossrun)
library(Rmpfr)

# variant of crossrunbin to show details for n=nmax:
crossrunbinshow <- function (nmax = 100, prob = 0.5, mult = 2, prec = 120, printn = FALSE) {
  nill <- Rmpfr::mpfr(0, prec)
  one <- Rmpfr::mpfr(1, prec)
  multm <- Rmpfr::mpfr(mult, prec)
  pm <- Rmpfr::mpfr(prob, prec)
  qm <- one - pm
  pmultm <- pm * multm
  qmultm <- qm * multm
  pat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  pbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  pt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  qt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
  for (nn in 2:nmax) {
    pat[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, 
                                                               nn))
    pbt[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, 
                                                               nn))
    rownames(pat[[nn]]) <- c(0:(nn - 1))
    rownames(pbt[[nn]]) <- c(0:(nn - 1))
    colnames(pat[[nn]]) <- c(1:nn)
    colnames(pbt[[nn]]) <- c(1:nn)
    pat[[nn]][1, nn] <- (pmultm^(nn - 1))
    pbt[[nn]][1, nn] <- (qmultm^(nn - 1))
    if (nn==nmax) {
      print("start above, contribution from no crossing:")
      suma <- pat[[nn]]
      print(asNumeric(pat[[nn]]))
      print(c("start below, contribution from no crossing:"))
      sumb <- pbt[[nn]]
      print(asNumeric(pbt[[nn]]))
    } # end show contribution for f=1 (no crossing)
    for (ff in 2:nn) {
      if (nn - ff + 1 <= ff - 1) {
        f1 <- ff
        pat[[nn]][2:(nn - f1 + 2), f1 - 1] <- pat[[nn]][2:(nn - 
                                                             f1 + 2), f1 - 1] + (pmultm^(f1 - 2)) * qmultm * 
          qbt[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 + 
                               1]
        pbt[[nn]][2:(nn - f1 + 2), f1 - 1] <- pbt[[nn]][2:(nn - 
                                                             f1 + 2), f1 - 1] + (qmultm^(f1 - 2)) * pmultm * 
          qat[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 + 
                               1]
        if (nn==nmax) {
          print(paste0("case 1, f=",f1))
          print(paste0("start below, last ",nn-f1+1," part"))
          print(asNumeric(pbt[[nn-f1+1]]))
          print("Corresponding contribution to entire series, start above:")
          contriba <- 0*pat[[nn]]
          contriba[2:(nn-f1+2),f1-1] <- (pmultm^(f1-2))*
            qmultm*qbt[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
          suma <- suma+contriba
          print(asNumeric(contriba))
          print(paste0("start above, last ",nn-f1+1," part"))
          print(asNumeric(pat[[nn-f1+1]]))
          print("Corresponding contribution to entire series, start below:")
          contribb <- 0*pbt[[nn]]
          contribb[2:(nn-f1+2),f1-1] <- (qmultm^(f1-2))*
            pmultm*qat[[nn-f1+1]][1:(nn-f1+1),nn-f1+1]
          sumb <- sumb+contribb
          print(asNumeric(contribb))
        } # end show contribution for each f in case 1
      }
      if (nn - ff + 1 > ff - 1) {
        f2 <- ff
        pat[[nn]][2:(nn - f2 + 2), f2 - 1] <- pat[[nn]][2:(nn - 
                                                             f2 + 2), f2 - 1] + (pmultm^(f2 - 2)) * qmultm * 
          qbt[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
        pat[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pat[[nn]][2:(nn - 
                                                                       f2 + 2), f2:(nn - f2 + 1)] + (pmultm^(f2 - 
                                                                                                               2)) * qmultm * pbt[[nn - f2 + 1]][1:(nn - f2 + 
                                                                                                                                                      1), f2:(nn - f2 + 1)]
        pbt[[nn]][2:(nn - f2 + 2), f2 - 1] <- pbt[[nn]][2:(nn - 
                                                             f2 + 2), f2 - 1] + (qmultm^(f2 - 2)) * pmultm * 
          qat[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
        pbt[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pbt[[nn]][2:(nn - 
                                                                       f2 + 2), f2:(nn - f2 + 1)] + (qmultm^(f2 - 
                                                                                                               2)) * pmultm * pat[[nn - f2 + 1]][1:(nn - f2 + 
                                                                                                                                                      1), f2:(nn - f2 + 1)]
        if (nn==nmax) {
          print(paste0("case 2, f=",f2))
          print(paste0("start below, last ",nn-f2+1," part"))
          print(asNumeric(pbt[[nn-f2+1]]))
          print("Corresponding contribution to entire series, start above:")
          contriba <- 0*pat[[nn]]
          contriba[2:(nn-f2+2),f2-1] <- (pmultm^(f2-2))*
            qmultm*qbt[[nn-f2+1]][1:(nn-f2+1),f2-1]
          contriba[2:(nn-f2+2),f2:(nn-f2+1)] <- (pmultm^(f2-2))*
            qmultm*pbt[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
          suma <- suma+contriba
          print(asNumeric(contriba))
          print(paste0("start above, last ",nn-f2+1," part"))
          print(asNumeric(pat[[nn-f2+1]]))
          print("Corresponding contribution to entire series, start below:")
          contribb <- 0*pbt[[nn]]
          contribb[2:(nn-f2+2),f2-1] <- (qmultm^(f2-2))*
            pmultm*qat[[nn-f2+1]][1:(nn-f2+1),f2-1]
          contribb[2:(nn-f2+2),f2:(nn-f2+1)] <- (qmultm^(f2-2))*
            pmultm*pat[[nn-f2+1]][1:(nn-f2+1), f2:(nn-f2+1)]
          sumb <- sumb+contribb
          print(asNumeric(contribb))
        } # end show contribution for each f in case 1
        
              }
    }
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
    if (printn) {
      print(nn)
      print(Sys.time())
    }
  }
  names(pat) <- paste("pat", 1:nmax, sep = "")
  names(pbt) <- paste("pbt", 1:nmax, sep = "")
  names(pt) <- paste("pt", 1:nmax, sep = "")
  names(qat) <- paste("qat", 1:nmax, sep = "")
  names(qbt) <- paste("qbt", 1:nmax, sep = "")
  names(qt) <- paste("qt", 1:nmax, sep = "")
  print("sum of contributions, start above:")
  print(asNumeric(suma))
  print("sum of contributions, start below:")
  print(asNumeric(sumb))
  return(list(pat = pat, pbt = pbt, pt = pt, qat = qat, qbt = qbt, 
              qt = qt))
} # end function crossrunbinshow

c7 <- crossrunbinshow(nmax=7)
asNumeric(c7$pat[[7]])
asNumeric(c7$pbt[[7]])

c7.6 <- crossrunbinshow(nmax=7, prob=0.6)
asNumeric(c7.6$pat[[7]])
asNumeric(c7.6$pbt[[7]])
