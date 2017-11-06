function [f]=chebcoef_obj(fun,alpha,B)
% Function specifies the OLS minimization between the real Expected Value function
% to be approximated and its approximation using the product between chebyshev coefficients 
% and the tensor product of the polynomial bases.

% Inputs
% fun    = function to be approximated.
% m      = number of nodes
% coef   = objective vector to be computed
% Ti     = polynomial base of state_var i

% Outputs
% f      = objective function



%% Objective Function

for i=1:nss
    aux(i,1)= sum(alpha.*T(i,:));
end
f=sum((fun(:,t)- aux).^2);

