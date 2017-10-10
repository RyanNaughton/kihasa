chebpoly_base <- function(n,x) {
#this function returns the Chebyshev Polynomial Base
  
  # Chebyshev Polynomials are defined by:
  # T_0(x)=1
  # T_1(x)=x
  # T_n+1(x)= 2*x*T_n(x) - T_n-1(x)
    
  ## Computing Matrices
    T <- matrix(0, nrow = length(x), ncol = n) # we can use the same base for both state variables
    
    # First-column: polynomial of degree 0
    T[,1] <- matrix(1, nrow = length(x), ncol = 1)
    #No need to replace first column of dTi and ddTi
    
    # Second-column: polynomial of degree 1
    T[,2] <- x
### tell R to skip if i < 3
    for (i in 3:n)
    {
      T[,i] <- 2*x*T[,i-1]-T[,i-2]
    }
    
 return(T)
}