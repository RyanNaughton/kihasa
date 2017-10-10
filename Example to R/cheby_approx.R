cheby_approx <- function(alpha,ncheby,h1,a1,dh,da,h,a) {
 
  Base <- kronecker(chebpoly_base(ncheby+1, dh*(h - h1) - 1),chebpoly_base(ncheby+1, da*(a - a1) - 1))
  EVapp <- sum(alpha*Base) #CHECK
  
}


# alpha
# ncheby <- G$ncheby  
# h1 <- ss$S$extmin1[t+1]
# a1 <- ss$S$extmin2[t+1]
# dh <- ss$S$dh[t+1]
# da <- ss$S$da[t+1]
# h <- hpr
# a <- apr