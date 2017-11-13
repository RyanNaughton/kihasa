clear all; clc;

% Testing:
a1=0;
a2=2;
b1=1;
b2=10; 
c1=5;
c2=8;

% Chevyshev Approximation
M = 10; % points to evaluate objective function
ncheby = M-2; %degrees

%%% Chebyshev nodes and re-scaled state space vector
z = flipud(cos(((2*(1:M)'-1)*pi)/(2*M))); % For Objective function

% parameter for expanded Chebyshev polynomial approximation 
extA = -(z(1)+1)*(a2-a1)/2/z(1);
extB = -(z(1)+1)*(b2-b1)/2/z(1);
extC = -(z(1)+1)*(c2-c1)/2/z(1);

%bounds
extminA = a1 - extA;
extminB = b1 - extB;
extminC = c1 - extC;
extmaxA = a2 + extA;
extmaxB = b2 + extB;
extmaxC = c2 + extC;

%gradient
dA = 2./(extmaxA-extminA);
dB = 2./(extmaxB-extminB);
dC = 2./(extmaxC-extminC);

%vector
Aext(1)=a1;
Bext(1)=b1;
Cext(1)=c1;
Aext(M)=a2;
Bext(M)=b2;
Cext(M)=c2;
for j=2:1:M-1
   Aext(j)= extminA + (1+z(j))./dA;
   Bext(j)= extminB + (1+z(j))./dB;
   Cext(j)= extminC + (1+z(j))./dC;
end

% Linear Vectors
A_vector = Aext; %1:10;
B_vector = Bext; %1:10;
C_vector = Cext; %1:10;

% Expanded Bounds (for the linear interpolation)
A_vector_wide = linspace(0,11,length(A_vector));
B_vector_wide = linspace(0,11,length(B_vector));
C_vector_wide = linspace(0,11,length(C_vector));

% Expand Vector (to calculate function)
SS_A = repmat(A_vector',[length(B_vector)*length(C_vector) 1]);
SS_B = repmat(kron(B_vector',ones(length(C_vector),1)),[length(A_vector) 1]);
SS_C = kron(C_vector',ones([length(A_vector)*length(B_vector),1]));

% Function
ABC_func = SS_A + SS_B + SS_C;
%reshaped for linear interpolation
rsp_func = reshape(ABC_func,[length(A_vector),length(B_vector),length(C_vector)]);

%%% Polynomial Bases and Derivatives %%%% 
T = chebpoly_base(ncheby+1,z); % base for Objective function
B = kron(T, kron(T,T)); % half of the numerator of the coeff
Num = ABC_func'*B; % numerator (bases*function) 
aux = diag(T'*T)'; % square of T
Den = kron(aux, kron(aux,aux));
alpha(1,(1:(ncheby+1)^3)) = Num./Den; % (9X9x9) = 729 coefficients

%algorithm
alp0=alpha*0;
nss = length(ABC_func);
alpha2=fmincon(@(alpha2) chebcoef_obj(ABC_func,alpha2,nss,B),alp0);

tic
% Linear Interpolation
for i = 1:1:length(ABC_func)
    A_next = SS_A(i) + 0.1;
    B_next = SS_B(i) + 0.5;
    C_next = SS_C(i) + 0.9;
    linear(i) = interpn(A_vector_wide,B_vector_wide,C_vector_wide, rsp_func, A_next,B_next,C_next);
end
toc %Elapsed time is 1 seconds.

tic
% Chebyshev Approximation
for i = 1:1:length(ABC_func)
    A_next = SS_A(i) + 0.1;
    B_next = SS_B(i) + 0.5;
    C_next = SS_C(i) + 0.9;
    EV(i) = cheby_approx(alpha,ncheby,extminA,extminB,extminC,dA,dB,dC,A_next,B_next,C_next); %Elapsed time is 0.032402 seconds.
    EV2(i) = cheby_approx(alpha2,ncheby,extminA,extminB,extminC,dA,dB,dC,A_next,B_next,C_next);
end
toc %Elapsed time is 0.15 seconds.

%% compare
x = 1:M^3;
plot(x,ABC_func)
plot(x,ABC_func,x,linear)
plot(x,ABC_func,x,linear,x,EV)
plot(x,ABC_func,x,linear,x,EV,x,EV2)
