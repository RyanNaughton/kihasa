### Function to Create State Space

thci_ss <- function(params,G,P,epsy_sim) {
  
  #return?
  
  w <- P$w
  inc_min <- kronecker(matrix(1,G$Ne,1),P$w_min) # 9x1 vector of w_min = 0.5
  inc_wage <- w*exp(G$eps_y) # vector of 1x3
  inc <- inc_min + inc_wage # 9x1 vector of inc = inc_min + inc_wage (inc_min = 0.5)