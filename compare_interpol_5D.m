clear all; clc;

% Testing:
a1=0;
a2=2;
b1=1;
b2=10; 
c1=5;
c2=8;
d1=1;
d2=3;
e1=4;
e2=6;

% Chevyshev Approximation
M = 7; % points to evaluate objective function
ncheby = M-2; %degrees

%%% Chebyshev nodes and re-scaled state space vector
z = flipud(cos(((2*(1:M)'-1)*pi)/(2*M))); % For Objective function

% parameter for expanded Chebyshev polynomial approximation 
extA = -(z(1)+1)*(a2-a1)/2/z(1);
extB = -(z(1)+1)*(b2-b1)/2/z(1);
extC = -(z(1)+1)*(c2-c1)/2/z(1);
extD = -(z(1)+1)*(d2-d1)/2/z(1);
extE = -(z(1)+1)*(e2-e1)/2/z(1);

%bounds
extminA = a1 - extA;
extminB = b1 - extB;
extminC = c1 - extC;
extminD = d1 - extD;
extminE = e1 - extE;
extmaxA = a2 + extA;
extmaxB = b2 + extB;
extmaxC = c2 + extC;
extmaxD = d2 + extD; 
extmaxE = e2 + extE;

%gradient
dA = 2./(extmaxA-extminA);
dB = 2./(extmaxB-extminB);
dC = 2./(extmaxC-extminC);
dD = 2./(extmaxD-extminD);
dE = 2./(extmaxE-extminE);

%vector
Aext(1)=a1;
Bext(1)=b1;
Cext(1)=c1;
Dext(1)=d1;
Eext(1)=e1;
Aext(M)=a2;
Bext(M)=b2;
Cext(M)=c2;
Dext(M)=d2;
Eext(M)=e2;
for j=2:1:M-1
   Aext(j)= extminA + (1+z(j))./dA;
   Bext(j)= extminB + (1+z(j))./dB;
   Cext(j)= extminC + (1+z(j))./dC;
   Dext(j)= extminD + (1+z(j))./dD;
   Eext(j)= extminE + (1+z(j))./dE;
end

% Linear Vectors
A_vector = Aext; %1:10;
B_vector = Bext; %1:10;
C_vector = Cext; %1:10;
D_vector = Dext;
E_vector = Eext;

% Bounds
A_vector_wide = linspace(0,11,length(A_vector));
B_vector_wide = linspace(0,11,length(B_vector));
C_vector_wide = linspace(0,11,length(C_vector));
D_vector_wide = linspace(0,11,length(D_vector));
E_vector_wide = linspace(0,11,length(E_vector));

% Expand Vector
SS_A = repmat(A_vector',[length(B_vector)*length(C_vector)*length(D_vector)*length(E_vector) 1]);
SS_B = repmat(kron(B_vector',ones(length(A_vector),1)),[length(C_vector)*length(D_vector)*length(E_vector) 1]);
SS_C = repmat(kron(C_vector',ones(length(B_vector)*length(A_vector),1)),[length(D_vector)*length(E_vector) 1]);
SS_D = repmat(kron(D_vector',ones(length(C_vector)*length(B_vector)*length(A_vector),1)),[length(E_vector) 1]);
SS_E = kron(E_vector',ones([length(A_vector)*length(B_vector)*length(C_vector)*length(D_vector),1]));

% Function
ABC_func = SS_A + SS_B + SS_C + SS_D + SS_E;
% reshaped for linear interpolation
rsp_func = reshape(ABC_func,[length(A_vector),length(B_vector),length(C_vector),length(D_vector),length(E_vector)]);

tic
% Linear Interpolation
for i = 1:1:length(ABC_func)
    A_next = SS_A(i) + 0.1;
    B_next = SS_B(i) + 0.5;
    C_next = SS_C(i) + 0.9;
    D_next = SS_D(i) + 0.3;
    E_next = SS_E(i) + 0.7;
    linear(i) = interpn(A_vector_wide,B_vector_wide,C_vector_wide,D_vector_wide,E_vector_wide, rsp_func, A_next,B_next,C_next,D_next,E_next);
end
toc %Elapsed time is 0.998017 seconds.

%%% Polynomial Bases and Derivatives %%%% 
T = chebpoly_base(ncheby+1,z); % base for Objective function
B = kron(T, kron(T, kron(T, kron(T,T)))); % half of the numerator of the coeff
Num = ABC_func'*B; % numerator (bases*function) 
aux = diag(T'*T)'; % square of T
Den = kron(aux, kron(aux, kron(aux, kron(aux,aux))));
alpha(1,(1:(ncheby+1)^5)) = Num./Den; % (9x9x9x9x9) = 59049 coefficients

tic
% Chebyshev Approximation
for i = 1:1:length(ABC_func)
    A_next = SS_A(i) + 0.1;
    B_next = SS_B(i) + 0.5;
    C_next = SS_C(i) + 0.9;
    D_next = SS_D(i) + 0.3;
    E_next = SS_E(i) + 0.7;
    EV(i) = cheby_approx_5D(alpha,ncheby,extminA,extminB,extminC,extminD,extminE,dA,dB,dC,dD,dE,A_next,B_next,C_next,D_next,E_next);
end
toc %Elapsed time is 0.121253 seconds.

%algorithm
alp0 = alpha*0;
nss = length(ABC_func);
tic
alpha2 = fmincon(@(alpha2) chebcoef_obj(ABC_func,alpha2,nss,B), alp0);
toc

tic
for i = 1:1:length(ABC_func)
    A_next = SS_A(i) + 0.1;
    B_next = SS_B(i) + 0.5;
    C_next = SS_C(i) + 0.9;
    D_next = SS_D(i) + 0.3;
    E_next = SS_E(i) + 0.7;
    EV2(i) =  cheby_approx_5D(alpha2,ncheby,extminA,extminB,extminC,extminD,extminE,dA,dB,dC,dD,dE,A_next,B_next,C_next,D_next,E_next);
end
toc %Elapsed time is 4.155648 seconds.

%% compare this using 5 variables (10)

x = 1:M^5;
plot(x,ABC_func)
plot(x,ABC_func,x,linear)
plot(x,ABC_func,x,linear,x,EV)
plot(x,ABC_func,x,linear,x,EV,x,EV2)
