obj_fun <- function(alpha,params,t,P,G,S,h,a,inc_min,inc_wage,c,minv,tinv) {

### Borrar ###
    #alpha=alp2(t,:)
    #a=S.a_(j,t)
    #h=S.h_(j,t)
    #inc=S.inc(i,t)
    #cc=S.c(j,k)
    #minv=S.minv(j,l)
    #tinv=S.tinv(m)
    #inc_min=S.inc_min(i,t)
    #inc_wage=S.inc_wage(i,t)
    
# ## in R
#   t <- 1
#   j <- 1
#   i <- 1
#   m <- 1
#   l <- 1
#   k <- 1
#   alpha <- alp2[t,]
#   a <- ss$S$a_[j,t]
#   h <- ss$S$h_[j,t]
#   inc <- ss$S$inc[i,t]
#   c <- ss$S$c[j,k]
#   minv <- ss$S$minv[j,l]
#   tinv <- ss$S$tinv[m]
#   inc_min <- ss$S$inc_min[i,t]
#   inc_wage <- ss$S$inc_wage[i,t]
  
      u <- utilfun(ss$S$c,P)
      gamma <- matrix(params[1:6], nrow = 3, ncol = 2) # check if it matches 3x2 matrix OR reshape(params([1:6]),3,2);
      phi <- params[7]
      rho <- params[8]
      delta <- params[9]
      inc <- inc_min + inc_wage # 9x1 matrix
      apr <- (1+P$r)*a + inc*(1-tinv) - P$p*minv - c
      hpr <- delta*CES(h,minv,tinv,gamma[t,],rho,phi,0)
      EV <- cheby_approx(alpha,G$ncheby,ss$S$extmin1[t+1],ss$S$extmin2[t+1],ss$S$dh[t+1],ss$S$da[t+1],hpr,apr)
      
      objf <-  u + P$beta*EV
      
  return <- list('objf'=objf,'apr'=apr,'hpr'=hpr)
}
