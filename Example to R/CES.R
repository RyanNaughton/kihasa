CES <- function(h,minv,tinv,gamma,rho,phi,eps_h) {
 
  CES <- (gamma[1]*h^phi + gamma[2]*minv^phi + (1-gamma[1]-gamma[2])*tinv^phi)^(rho/phi)+eps_h 
  
}
