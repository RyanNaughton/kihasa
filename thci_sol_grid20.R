thci_sol_grid20 <- function(params,P,G,B,S) {
  
  gamma <- params[1:6]
  phi <- params[7]
  rho <- params[8]
  delta <- params[9]
  
  W <- matrix(0, nrow = G$nss, ncol = G$ntime); # EValue function
  #V=zeros(G.nss,G.Ne,G.nper); # Value function
  ### split into 3 matrices by G.nper and join
  V1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V <- list('V1'=V1, 'V2'=V2, 'V3'=V3)
  
  W[,G$ntime] <- P$eta*ss$S$h_[,G$ntime]^(1-P$sigma)/(1-P$sigma)+P$tau1*(1-exp(-ss$S$a_[,G$ntime]))
  
  alp2 <- matrix(0, nrow = G$ntime-1, ncol = (G$ncheby+1)^2) #empty matrix of 3x29^2 = 3x1841
  
  ### start with empty matrices
  VV1 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  VV2 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  VV3 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  bc1_val <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  bc2_val <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  bc3_val <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  auxV <- matrix(0, nrow = G$Ntinv, ncol = 1)
  Minv1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  Minv2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  Minv3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  Tinv1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  Tinv2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  Tinv3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  C1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  C2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  C3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  A1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  A2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  A3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  A <- list('A1'=A1, 'A2'=A2, 'A3'=A3)
  
  for (t in G$ntime-1:1)
  {
    ## A.Polynomial Approximation of Expected Value Function
    # 1. Input is W(t+1)
    # 2. Output is vector of polynomial coefficients alp2(t)
    
    coef1 <- W[,t+1]%*%kronecker(ss$B$B,ss$B$B)  #matrix of 1x29^2 = 1x1841
    aux <- diag(t(ss$B$B)%*%ss$B$B)              #matrix of 1x29
    aux2 <- matrix(kronecker(aux,aux), nrow = 1) #matrix of 1x29^2 = 1x1841
    alp2[t,(1:(G$ncheby+1)^2)] <- coef1/aux2     #fill out one row of matrix of 3x29^2 = 3x1841

    ## B. Optimization using polynomial tensor product  
    for (i in 1:G$Ne)
    {
      for (j in 1:G$nss)
      {
        
        for (k in 1:G$Nc)
        {
          for (l in 1:G$Nc)
          {
            
            objf1 <- obj_fun(alp2[t,], params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[1])
            objf2 <- obj_fun(alp2[t,], params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[2])
            objf3 <- obj_fun(alp2[t,], params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[3])
            
            VV1[k,l] <- objf1$objf[k,l] 
            VV2[k,l] <- objf2$objf[k,l]
            VV3[k,l] <- objf3$objf[k,l]
            
            bc1_val[k,l] <- objf1$apr>=ss$S$a1[t+1] && objf1$apr<=ss$S$a2[t+1] && objf1$hpr>=ss$S$h1[t+1] && objf1$hpr<=ss$S$h2[t+1]
            bc2_val[k,l] <- objf2$apr>=ss$S$a1[t+1] && objf2$apr<=ss$S$a2[t+1] && objf2$hpr>=ss$S$h1[t+1] && objf2$hpr<=ss$S$h2[t+1]
            bc3_val[k,l] <- objf3$apr>=ss$S$a1[t+1] && objf3$apr<=ss$S$a2[t+1] && objf3$hpr>=ss$S$h1[t+1] && objf3$hpr<=ss$S$h2[t+1]
            
          }
        }
        
        #bc1_val[bc1_val < 1] <- NA
        #bc2_val[bc2_val < 1] <- NA
        #bc3_val[bc3_val < 1] <- NA
        
        Ind_bc1 <- VV1*bc1_val
        Ind_bc2 <- VV2*bc2_val
        Ind_bc3 <- VV3*bc3_val
        
        auxV[1] <- max(Ind_bc1) # V_tinv1
        auxV[2] <- max(Ind_bc2) # V_tinv2
        auxV[3] <- max(Ind_bc3) # V_tinv3
        
        Indaux1 <- ind2sub(dim(Ind_bc1),which.max(Ind_bc1))
        Indaux2 <- ind2sub(dim(Ind_bc2),which.max(Ind_bc2))
        Indaux3 <- ind2sub(dim(Ind_bc3),which.max(Ind_bc3))
        
        I_c1 <- Indaux1[1]
        I_c2 <- Indaux2[1]
        I_c3 <- Indaux3[1]
        
        I_minv1 <- Indaux1[2]
        I_minv2 <- Indaux2[2]
        I_minv3 <- Indaux3[2]
        
        if(t==1) V$V1[j,i] <- max(auxV)
        if(t==2) V$V2[j,i] <- max(auxV)
        if(t==3) V$V3[j,i] <- max(auxV)
        
        
        Minv1[j,i] <- ss$S$minv[j,I_minv1]
        Minv2[j,i] <- ss$S$minv[j,I_minv2]
        Minv3[j,i] <- ss$S$minv[j,I_minv3]
        
        C1[j,i] <- ss$S$c[j,I_c1]
        C1[j,i] <- ss$S$c[j,I_c2]
        C1[j,i] <- ss$S$c[j,I_c3]
        
        Tinv1[j,i] <- ss$S$tinv[which.max(auxV)]
        Tinv2[j,i] <- ss$S$tinv[which.max(auxV)]
        Tinv3[j,i] <- ss$S$tinv[which.max(auxV)]
        
        if(t==1) A$A1[j,i] <- (1+P$r)*ss$S$a_[j,t] + ss$S$inc[i,t]*(1-Tinv1[j,i]) - P$p*Minv1[j,i] - C1[j,i]
        if(t==2) A$A2[j,i] <- (1+P$r)*ss$S$a_[j,t] + ss$S$inc[i,t]*(1-Tinv2[j,i]) - P$p*Minv2[j,i] - C2[j,i]
        if(t==3) A$A3[j,i] <- (1+P$r)*ss$S$a_[j,t] + ss$S$inc[i,t]*(1-Tinv3[j,i]) - P$p*Minv3[j,i] - C3[j,i]

        print(j)
      }
      
      print(i)
    }
    
    W[,t] <- pi^(-1/2)*V$V3*G$wt
    #W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt;
    # check integration of cons and inv
    
    print(t)
  }
  
  return <- list('W'=W,'A'=A,'Minv'=Minv,'Tinv'=Tinv,'C'=C)
}
