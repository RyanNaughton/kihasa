
function [dT1,dT2,ddT1,ddT2]=Dchebpoly_deriv(n,x,T,dx1,dx2)

%this function returns the matrix of first and second derivatives of the
%polynomial base T.

% Let dx=dx/dstate_var: change in polynomial domain due to a change in
% state variable. In general dx/dstate_var= 2/(state_UB - state_LB). Comes from
% change in variables for rescaling
% state_var  = (x+1)*(state_UB - state_LB)/2 + state_LB;
% (state_var - state_LB)*2/(state_UB - state_LB) - 1 = x

% Chebyshev Polynomials of degree n and its derivatives defined by:
% T_0(x)=1
% T_1(x)=x
% T_n+1(x)= 2*x*T_n(x) - T_n-1(x)

% dT_0(x)=0
% dT_1(x)=dx 
% dT_n+1(x)=2*T_n(x)*dx + 2*x*dT_n(x) - dT_n-1(x);

% ddT_0(x)=0;
% ddT_1(x)=0;
% ddT_n+1(x)=4*dT_n(x)*dx + 2*x*ddT_n(x) - ddT_n-1(x);


%% Computing Matrices

% Set initial space
dT1=zeros(length(x),n);
dT2=zeros(length(x),n);
ddT1=zeros(length(x),n);
ddT2=zeros(length(x),n);

% First-column: polynomial of degree 0: No need to replace first column of dTi and ddTi

% Second-column: polynomial of degree 1

dT1(:,2)=dx1;
dT2(:,2)=dx2;

%No need to replace second column of ddTi

for i=3:n
    dT1(:,i)=2*dx1.*T(:,i-1) + 2*x.*dT1(:,i-1) - dT1(:,i-2);
    dT2(:,i)=2*dx2.*T(:,i-1) + 2*x.*dT2(:,i-1) - dT2(:,i-2);
    ddT1(:,i)=4*dx1.*dT1(:,i-1) + 2*x.*ddT1(:,i-1)- ddT1(:,i-2);
    ddT2(:,i)=4*dx2.*dT2(:,i-1) + 2*x.*ddT2(:,i-1)- ddT2(:,i-2);
end;
end


