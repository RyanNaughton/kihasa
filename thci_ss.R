### Function to Create State Space

thci_ss <- function(params,G,P,epsy_sim) {
  
  w <- P$w
  inc_min  <- matrix(P$w_min, nrow = length(P$w_min), ncol = G$Ne)  # 3x3 matrix (3 periods x 3 shocks)
  inc_wage <- matrix(w, nrow = length(w), ncol = G$Ne)*exp(G$eps_y) # 3x3 matrix (3 periods x 3 shocks)
  inc <- inc_min + inc_wage # 3x3 matrix (3 periods x 3 shocks) = inc_min + inc_wage (inc_min = 0.5)
 
  inc_min_sim  <- matrix(P$w_min, nrow = G$npop, ncol = length(P$w_min))   # 3xNpop matrix (3 periods x Npop people)
  inc_wage_sim <- matrix(w, nrow = G$npop, ncol = length(w))*exp(epsy_sim) # 3xNpop matrix (3 periods x Npop people)
  inc_sim <- inc_min_sim + inc_wage_sim 
  
  INC <- matrix(c(inc_sim[,1], inc_sim[,2]/(1+P$r),inc_sim[,3]/(1+P$r)^2), nrow = nrow(inc_sim), ncol = G$Ne)
  PI  <- apply(INC,1,sum)

  ###### Bounds State Vectors #######
  
  h1 <- c(0,0,0,0)
  h2 <- c(2.5,2.5,2.5,2.5)
  #a1=[0,0,0,0]; # with credit constraints
  a1 <- c(-inc[1,1]/(1+P$r)-inc[1,2]/(1+P$r)^2-inc[1,3]/(1+P$r)^3,-inc[1,2]/(1+P$r)-inc[1,3]/(1+P$r)^2,-inc[1,3]/(1+P$r),0)
  a2 <- c(10,10,10,10) 
  
  ### Chebyshev nodes and re-scaled state space vector
  
  z <- rev(cos(((2*(1:G$M)-1)*pi)/(2*G$M))) # For Objective function 
         
                  ext1 <- -(z[1]+1)*(h2-h1)/2/z[1]   # parameter for expanded Chebyshev polynomial approximation 
                  ext2 <- -(z[1]+1)*(a2-a1)/2/z[1]   # parameter for expanded Chebyshev polynomial approximation 
                  extmin1 <- h1 - ext1               # expanded minimal state value
                  extmin2 <- a1 - ext2               # expanded minimal state value
                  extmax1 <- h2 + ext1               # expanded minimal state value
                  extmax2 <- a2 + ext2               # expanded minimal state value
                  dh <- 2/(extmax1-extmin1)          # dz/dh     
                  da <- 2/(extmax2-extmin2)          # dz/da
                  
                  hext <- matrix(0, nrow = G$M, ncol = G$ntime)
                  aext <- matrix(0, nrow = G$M, ncol = G$ntime)
                  hext[1,] <- h1
                  aext[1,] <- a1
                  hext[G$M,] <- h2
                  aext[G$M,] <- a2
                  for (j in 2:G$M-1)
                  {
                    hext[j,] <- extmin1 + (1+z[j])/dh
                    aext[j,] <- extmin2 + (1+z[j])/da
                  }
                  
                  ## State vector
                  h_ <- kronecker(hext, matrix(1,G$M,1)) # Asset vector
                  a_ <- kronecker(matrix(1,G$M,1), aext) # HC vector
                  
                  ## RHS Budget Constraint
                  lb <- matrix(0, nrow = 2, ncol = 3)
                  b1 <- matrix(inc[,1], nrow = G$nss, ncol = length(inc[,1]), byrow = TRUE)+(1+P$r)*matrix(a_[,1], nrow =length(a_[,1]), ncol = G$Ne)
                  b2 <- matrix(inc[,2], nrow = G$nss, ncol = length(inc[,2]), byrow = TRUE)+(1+P$r)*matrix(a_[,2], nrow =length(a_[,2]), ncol = G$Ne)
                  b3 <- matrix(inc[,3], nrow = G$nss, ncol = length(inc[,3]), byrow = TRUE)+(1+P$r)*matrix(a_[,3], nrow =length(a_[,3]), ncol = G$Ne)
                  b <- cbind(b1, b2, b3)

## Grid for choice variables 

  ubT <- inc[1,3] + (1+P$r)*a_[,3]

  c <- matrix(0, nrow = G$nss, ncol = G$Nc)
  minv <- matrix(0, nrow = G$nss, ncol = G$Nc)
  for (j in 1:G$nss)
  {
    c[j,] <- seq(0,ubT[j],length=G$Nc)
    minv[j,] <- seq(0,ubT[j],length=G$Nc)
  }
  
  tinv <- c(0,0.5,1)
  
  S <- list('inc'=inc,'inc_min'=inc_min,'inc_wage'=inc_wage,'h1'=h1,'a1'=a1,'extmin1'=extmin1,'extmin2'=extmin2,
            'c'=c,'minv'=minv,'tinv'=tinv,'dh'=dh,'da'=da,'hext'=hext,'aext'=aext,'h_'=h_,'a_'=a_,'b'=b,'lb'=lb,
            'ubT'=ubT,'a2'=a2,'h2'=h2)
  
  
  ### Polynomial Bases and Derivatives #### 
  
  B <- chebpoly_base(G$ncheby+1,z) #Polynomial base for Objective function
  T <- kronecker(B,B)                
  
  for t=1:1:G.ntime
  [dBh(:,:,t),dBa(:,:,t),ddBh(:,:,t),ddBa(:,:,t)]=Dchebpoly_deriv(G.ncheby+1,z,B,dh(t),da(1,t));
  end
  
  
  ### Polynomial Bases Policy Function ####
  # Basis for states
  B_sim <- chebpoly_base(G$ncheb_pol,z)
  
  # Basis for Income Shocks
  z_epsw= 2*(G.eps_y(:,1)-G.eps_y(1,1))/(G.eps_y(G.Ne,1)-G.eps_y(1,1))-1; 
  B_ew=chebpoly_base(G.Ne-1,z_epsw);
  
  
  # Multidimensional base
  T_sim=kron(B_sim,kron(B_sim,B_ew));
  
  B  = struct('B',B,'dBh',dBh,'dBa',dBa,'ddBh',ddBh,'ddBa',ddBa,'B_sim',B_sim,'B_ew',B_ew,'T_sim',T_sim);
  
return(S,B,inc_min_sim,inc_wage_sim,inc_sim,INC,PI)
}