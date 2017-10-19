function u = utilfun(c, sigma, psi, kappa, theta1, theta2, theta3, m, n, K)

u=(c^(1-sigma))/(1-sigma) + psi + kappa*(1+theta1*m+theta2*n+theta3*K);
