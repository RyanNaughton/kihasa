thci_sol_grid20 <- function(params,P,G,B,S) {
  
  gamma <- params[1:6]
  phi <- params[7]
  rho <- params[8]
  delta <- params[9]
  
  #V=zeros(G.nss,G.Ne,G.nper); # Value function
  ### split into 3 matrices by G.nper and join
  V1 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V2 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V3 <- matrix(0, nrow = G$nss, ncol = G$Ne)
  V <- list('V1'=V1, 'V2'=V2, 'V3'=V3) ### check: as list?
  W <- matrix(0, nrow = G$nss, ncol = G$ntime); # EValue function
  
  W[,G$ntime] <- P$eta*ss$S$h_[,G$ntime]^(1-P$sigma)/(1-P$sigma)+P$tau1*(1-exp(-ss$S$a_[,G$ntime]))
  
  
  alp2 <- matrix(0, nrow = G$ntime-1, ncol = (G$ncheby+1)^2) #empty matrix of 3x29^2 = 3x1841
  
  ### start empty matrix
  VV1 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  VV2 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  VV3 <- matrix(0, nrow = G$Nc, ncol = G$Nc)
  
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
    for i=1:1:G.Ne       
    for j=1:1:G.nss  
    
    for k=1:G.Nc
    for l=1:G.Nc
    #for m=1:G.Ntinv
    
    objf1 <- obj_fun(alp2[t,],params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[1])
    objf2 <- obj_fun(alp2[t,],params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[2])
    objf3 <- obj_fun(alp2[t,],params,t,P,G,S, ss$S$h_[j,t], ss$S$a_[j,t], ss$S$inc_min[i,t], ss$S$inc_wage[i,t], ss$S$c[j,k], ss$S$minv[j,l], ss$S$tinv[3])
    
    VV1[k,l] <- objf1$objf[k,l]
    VV2[k,l] <- objf2$objf[k,l]
    VV3[k,l] <- objf3$objf[k,l]
    
    #[VV(k,l,m),apr(k,l,m), hpr(k,l,m)]=obj_fun(alp2(t,:),params,t,P,G,S,S.h_(j,t),S.a_(j,t),S.inc_min(i,t),S.inc_wage(i,t),S.c(j,k),S.minv(j,l),S.tinv(m));
    #bc_val(k,l,m)=(apr(k,l,m)>=S.a1(t+1) && apr(k,l,m)<=S.a2(t+1) && hpr(k,l,m)>=S.h1(t+1) && hpr(k,l,m)<=S.h2(t+1));                                                   
    
    #end
    end 
    end
    
    bc_val = double(bc_val);
    bc_val(bc_val == 0)=NaN;
    auxV=zeros(G.Ntinv,1);
    
    Ind_bc_check=VV.*bc_val;
    for m=1:G.Ntinv
    Ind_bc=VV(:,:,m).*bc_val(:,:,m);
    [V_tinv,Indaux] = max(Ind_bc(:));
    auxV(m)=V_tinv;
    [I_c(m), I_minv(m)] = ind2sub(size(Ind_bc),Indaux);
    end
    [V(j,i,t), Ind]=max(auxV(:));                  
    Minv(j,i,t)=S.minv(j,I_minv(Ind));
    C(j,i,t)=S.c(j,I_c(Ind));
    Tinv(j,i,t)=S.tinv(Ind);
    
    A(j,i,t)=(1+P.r)*S.a_(j,t) + S.inc(i,t)*(1-Tinv(j,i,t)) - P.p*Minv(j,i,t) - C(j,i,t);
    #H(j,i,t)=delta*CES(S.h_(j,t),opty(1),gamma(t),rho,phi);
    j    
    end
    i
    end
    
  }

  
  
W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt;
# check integration of cons and inv


print(t)
  
  return <- list('W'=W,'A'=A,'Minv'=Minv,'Tinv'=Tinv,'C'=C)
}
