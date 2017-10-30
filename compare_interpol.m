clear all; clc;

% Testing:

A_vector = 1:5;
B_vector = 1:5;
C_vector = 1:5;

A_vector_wide = linspace(-1,10,length(A_vector));
B_vector_wide = linspace(-1,10,length(B_vector));
C_vector_wide = linspace(-1,10,length(C_vector));

%%%%%%%% or adjust vector using this:
% ubT= inc(1,3) + (1+P.r)*a_(:,3);
% for j=1:G.nss
%   c(j,:)= linspace(0,ubT(j),G.Nc);
%   minv(j,:)= linspace(0,ubT(j),G.Nc);
% end
% tinv= [0,0.5,1];

% Expand Vector
SS_A = repmat(A_vector',[length(B_vector)*length(C_vector) 1]);
SS_B = repmat(kron(B_vector',ones(length(C_vector),1)),[length(A_vector) 1]);
SS_C = kron(C_vector',ones([length(A_vector)*length(B_vector),1]));

% Function
ABC_func = SS_A + SS_B + SS_C;
rsp_func = reshape(ABC_func,[5,5,5]); %reshaped for linear interpolation

%linear interpolation
for i = 1:1:length(ABC_func)
    
    A_next = A_vector(i) + 0.1;
    B_next = B_vector(i) + 0.5;
    C_next = C_vector(i) + 0.9;

    linear(i) = interpn(A_vector_wide,B_vector_wide,C_vector_wide, rsp_func, A_next,B_next,C_next);
end

%Chevyshev
M= 5; % points to evaluate objective function
ncheby = M-2;
%ncheb_pol=6;

%[S,B,inc_min_sim,inc_wage_sim,inc_sim,~,PI]=thci_ss(params0,G,P,epsy_sim);
%%%%%% Bounds State Vectors %%%%%%% use same as above
% h1=[0,0,0,0];
% h2=[2.5,2.5,2.5,2.5];
% %a1=[0,0,0,0]; % with credit constraints
% a1(:,:)=[-inc(1,1)/(1+P.r)-inc(1,2)/(1+P.r)^2-inc(1,3)/(1+P.r)^3,-inc(1,2)/(1+P.r)-inc(1,3)/(1+P.r)^2,-inc(1,3)/(1+P.r),0];
% a2=[10,10,10,10]; 
%%% Chebyshev nodes and re-scaled state space vector
z= flipud(cos(((2*(1:G.M)'-1)*pi)/(2*G.M))); % For Objective function 
ext1 = -(z(1)+1)*(h2-h1)/2/z(1);   % parameter for expanded Chebyshev polynomial approximation 
ext2 = -(z(1)+1)*(a2-a1)/2/z(1);   % parameter for expanded Chebyshev polynomial approximation 
extmin1 = h1 - ext1;               % expanded minimal state value
extmin2 = a1 - ext2;               % expanded minimal state value
extmax1 = h2 + ext1;               % expanded minimal state value
extmax2 = a2 + ext2;               % expanded minimal state value
dh = 2./(extmax1-extmin1);          % dz/dh     
da = 2./(extmax2-extmin2);          % dz/da
hext(1,[1:G.ntime])=h1;
aext(1,[1:G.ntime])=a1;
hext(G.M,[1:G.ntime])=h2;
aext(G.M,[1:G.ntime])=a2;
for j=2:1:G.M-1
   hext(j,[1:G.ntime])= extmin1 + (1+z(j))./dh; %a_next_poly
   aext(j,[1:G.ntime])= extmin2 + (1+z(j))./da;
end
%%%%%%%% state vector       
h_ = kron(hext , ones(G.M,1)); % Asset vector
a_ = kron(ones(G.M,1), aext); % HC vector
S = struct('inc',inc,'inc_min',inc_min,'inc_wage',inc_wage,'h1',h1,'a1',a1,'extmin1',extmin1,'extmin2',extmin2,'c',c,'minv',minv,'tinv',tinv,...
            'dh',dh,'da',da,'hext',hext,'aext',aext,'h_',h_,'a_',a_,'b',b,'lb',lb,'ubT',ubT,'a2',a2,'h2',h2);
%%% Polynomial Bases and Derivatives %%%% 
[B]=chebpoly_base(G.ncheby+1,z); %Polynomial base for Objective function
T=kron(B,kron(B,B));
% for t=1:1:G.ntime
% [dBh(:,:,t),dBa(:,:,t),ddBh(:,:,t),ddBa(:,:,t)]=Dchebpoly_deriv(G.ncheby+1,z,B,dh(t),da(1,t));
% end
%%% Polynomial Bases Policy Function %%%%
% Basis for states
% B_sim=chebpoly_base(G.ncheb_pol,z);
% % Basis for Income Shocks
% z_epsw= 2*(G.eps_y(:,1)-G.eps_y(1,1))/(G.eps_y(G.Ne,1)-G.eps_y(1,1))-1; 
% B_ew=chebpoly_base(G.Ne-1,z_epsw);
% % Multidimensional base
% T_sim=kron(B_sim,kron(B_sim,B_ew));
% B  = struct('B',B,'dBh',dBh,'dBa',dBa,'ddBh',ddBh,'ddBa',ddBa,...
%             'B_sim',B_sim,'B_ew',B_ew,'T_sim',T_sim);
        
%[W,A,Minv,Tinv,C]=thci_sol_grid20(params0,P,G,B,S)        
% W(:,G.ntime)=P.eta*S.h_(:,G.ntime).^(1-P.sigma)/(1-P.sigma)+ P.tau1*(1-exp(-S.a_(:,G.ntime))); 

%for t=G.ntime-1:-1:1
t=1;
       coef1=W(:,t+1)'*T %kron(B.B,B.B);
       aux=diag(B.B'*B.B)'; %T
       aux2=kron(aux,aux);
       alp2(t,(1:(G.ncheby+1)^2))=coef1./aux2;
       
%        for i=1:1:G.Ne       
%         for j=1:1:G.nss  
% 
%             for k=1:G.Nc
%                 for l=1:G.Nc
%                     for m=1:G.Ntinv
    
    %[VV(k,l,m),apr(k,l,m), hpr(k,l,m)]=obj_fun(alp2(t,:),params,t,P,G,S,S.h_(j,t),S.a_(j,t),S.inc_min(i,t),S.inc_wage(i,t),S.c(j,k),S.minv(j,l),S.tinv(m));
    
    %%%% these are the next values:
%     apr=(1+P.r)*a + inc*(1-tinv) - P.p*minv - c;
%     hpr=delta*CES(h,minv,tinv,gamma(t,:),rho,phi,0);  
    EV=cheby_approx(alpha,G.ncheby,S.extmin1(t+1),S.extmin2(1,t+1),S.dh(t+1),S.da(1,t+1),hpr,apr)                       
