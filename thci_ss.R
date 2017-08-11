### Function to Create State Space

thci_ss <- function(params,G,P,epsy_sim) {
  
  #return?
  
  w <- P$w
  inc_min  <- matrix(P$w_min, nrow = length(P$w_min), ncol = G$Ne)  # 3x3 matrix (3 periods x 3 shocks)
  inc_wage <- matrix(w, nrow = length(w), ncol = G$Ne)*exp(G$eps_y) # 3x3 matrix (3 periods x 3 shocks)
  inc <- inc_min + inc_wage # 3x3 matrix (3 periods x 3 shocks) = inc_min + inc_wage (inc_min = 0.5)
 
  inc_min_sim  <- matrix(P$w_min, nrow = G$npop, ncol = length(P$w_min))   # 3xNpop matrix (3 periods x Npop people)
  inc_wage_sim <- matrix(w, nrow = G$npop, ncol = length(w))*exp(epsy_sim) # 3xNpop matrix (3 periods x Npop people)
  inc_sim <- inc_min_sim + inc_wage_sim 
  
 INC=[inc_sim(:,1),inc_sim(:,2)/(1+P.r),inc_sim(:,3)/(1+P.r)^2];
 PI=sum(INC,2);