############################################
##### Parental Investments in Children #####
############################################

rm(list=ls())
#install.packages('statmod')
require("statmod")
source("ind2sub.R")

#load functions
source("CES.R")
source("utilfun.R")
source("thci_ss.R")
source("obj_fun.R")
source("cheby_approx.R")
source("chebpoly_base.R")
source("Dchebpoly_deriv.R")

##### 1. PARAMETERS #####

### Structural Estimates ###

gamma   <-  c(0.5991, 0.4602, 0.6194, 0.20045, 0.2699, 0.1903)
phi     <- 0.4282 
rho     <- 0.7487 
delta   <- 2.23

params0 <- c(gamma, phi, rho, delta)

gamma_b <- matrix(c(0.1,0.8,0.3,0.9,0.25,0.9,0.0,0.5,0.0,0.5,0.0,0.5), nrow=6, byrow=TRUE)
phi_b   <- c(-1.2,0.85)
rho_b   <- c(0.2,0.85)
delta_b <- c(1,4)
bounds  <- c(gamma_b,phi_b,rho_b,delta_b)

params <- params0 # temporary

### Calibrated Estimates ###
  
sigma     <- 0.5            #Elasticity of Substitution
beta      <- 0.96
sigma_eps <- sqrt(0.05)     #shocks to income
sigma_v   <- sqrt(0)
eta       <- 12             #terminal value parameter (human capital)
tau1      <- 12             #terminal value parameter (assets) - If tau/eta too small V_t<-1 is not concave in some cases
r         <- 1/beta-1       #interest rate
p         <- 1              #relative price of investment and consumption
w_min     <- c(0.5,0.5,0.5) ## CHANGE THIS TO TEST TIMING OF INCOME
mu        <- 0              #mean of shock(s)
cmin      <- w_min
u0        <- 0
u1        <- 0.5
u2        <- -0.5
vcv       <- matrix(c(sigma_eps^2,0,0,sigma_v^2), nrow=2, byrow=TRUE)
detV      <- det(vcv)
w         <- c(1,1,1)        
  
P  <- list('beta'=beta,'sigma'=sigma,'sigma_eps'=sigma_eps,'eta'=eta,'tau1'=tau1,'r'=r,'p'=p,'delta'=delta,
           'w_min'=w_min,'mu'=mu,'cmin'=cmin,'u0'=u0,'u1'=u1,'u2'=u2,'w'=w,'sigma_v'=sigma_v)
  
##### 2. STATE SPACE #####
    
# Total states: A=15; H:15; eps_y=3; eps_h1=3 
  
M  <- 30 #Chevyshev points to evaluate objective function
M2 <- M
Ne <- 3  #Gauss-Hermite Points
nss <- M^2 #state space A and H1
ncheby    <- M-2
npop      <- 40 #40000
ntime     <- 4
nper      <- 3
Nc        <- 20
Ntinv     <- 3
ncheb_pol <- 6
neps <- 2
out <- gauss.quad(Ne) # Delivers integration nodes and weights for epsilon
e <- out$nodes
wt <- out$weights
eps_y <- sqrt(2)*e*sigma_eps # error vector
eps_h <- sqrt(2)*e*sigma_v
  
G  <- list('M'=M,'M2'=M2,'Ne'=Ne,'nss'=nss,'ncheby'=ncheby,'npop'=npop,'ntime'=ntime,'nper'=nper,
           'ncheb_pol'=ncheb_pol,'Nc'=Nc,'Ntinv'=Ntinv,'eps_y'=eps_y,'eps_h'=eps_h,'wt'=wt)
  
##### 3. SHOCKS #####
    
epsy_sim <- kronecker(matrix(rnorm(npop*nper),nrow = npop,ncol = nper),sigma_eps)
epsh_sim <- kronecker(matrix(rnorm(npop*nper),nrow = npop,ncol = nper),sigma_v)

##### 4. TEST FUNCTIONS #####

ss <- thci_ss(params0,G,P,epsy_sim)
# [W,A,Minv,Tinv,C]=thci_sol_grid20(params0,P,G,B,S);  
