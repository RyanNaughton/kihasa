function [S,B,inc_min_sim,inc_wage_sim,inc_sim,INC,PI]=thci_ss(params,G,P,epsy_sim)

w=P.w;

inc_min=repmat(P.w_min,G.Ne,1);
inc_wage=w.*exp(G.eps_y); 
inc=inc_min+inc_wage;

inc_min_sim=repmat(P.w_min,G.npop,1);
inc_wage_sim=repmat(w,G.npop,1).*exp(epsy_sim); 
inc_sim=inc_min_sim + inc_wage_sim; 

INC=[inc_sim(:,1),inc_sim(:,2)/(1+P.r),inc_sim(:,3)/(1+P.r)^2];
PI=sum(INC,2);
%pc_inc=prctile(INC,G.pc);
%pc_pi=prctile(PI,G.pc);
       
%%%%%% Bounds State Vectors %%%%%%%

h1=[0,0,0,0];
h2=[2.5,2.5,2.5,2.5];
%a1=[0,0,0,0]; % with credit constraints
a1(:,:)=[-inc(1,1)/(1+P.r)-inc(1,2)/(1+P.r)^2-inc(1,3)/(1+P.r)^3,-inc(1,2)/(1+P.r)-inc(1,3)/(1+P.r)^2,-inc(1,3)/(1+P.r),0];
a2=[10,10,10,10]; 

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
   hext(j,[1:G.ntime])= extmin1 + (1+z(j))./dh;
   aext(j,[1:G.ntime])= extmin2 + (1+z(j))./da;
end


%% State vector
       
h_ = kron(hext , ones(G.M,1)); % Asset vector
a_ = kron(ones(G.M,1), aext); % HC vector

%% RHS Budget Constraint
lb=[0,0,0;0,0,0];
for t=1:G.nper
b(:,:,t)=repmat(inc(:,t)',G.nss,1) +(1+P.r)*repmat(a_(:,t),1,G.Ne);    
end

%% Grid for choice variables 
% Note 1 (lower bounds): Lower bounds of continous choices (c, minv) are set to zero (in
% practice it could be larger than zero because individuals have a minimum
% income. Need to try another version with C_lb, Minv_lb = w_min

% Note 2 (previous upper bounds): Previously we were setting the upper bound for
% C_t=3 and Minv_t=3 = inc(eps_low,t=3) +  inc(eps_low,t=3)/(1+r).
% Reasoning is households can consume: a) per period income evaluated at
% the lower shock not to allow negative consumption; + b) the max amount of
% debt they allowed to hold in period 3, given by the lower bound of A_t=3.
% But the second term only applies for borrowers, so it's wrong. 

% Note 3 (new upper bounds): Should be literally computed from the budget
% constraint. Assuming the lower bounds are zero (see Note 1), the upper
% bound for consumption in period 3 is: 
%
%        C_3,ub= inc(eps_3,lb) + (1+r)A_3,ub - p*Minv_3,lb - A_4,lb
%
% a) Minv_3,lb=0 for now. If C_3opt + MInv_3opt> resources, this is
%    sorted out with BC restrictions in our algorithm.
% b) A_4,lb=0. Households can choose to leave bequests but if they have
% strong preferences for own consumption/investments in the present they
% can do it as long as the final asset position is non-negaive.
% c) A_3,ub=70, but here I have doubts. We could use 70 to compute a
% general upper bound for consumption. But alternatively, since the policy
% function of consumption is conputed for each value of the grid of A_3,
% maybe the correct thing to do is to set A_3,ub(j)=A_3(j), where
% j=1,...,M is the grid point. So alterntives are:
%
%       C_3,ub   = inc(eps_3,lb) + (1+r)*70
%       C_3,ub(j)= inc(eps_3,lb) + (1+r)*A_3(j)

% Old one: ubT=inc(1,3)+ inc(1,3)/(1+P.r); 

ubT= inc(1,3) + (1+P.r)*a_(:,3);
for j=1:G.nss
  c(j,:)= linspace(0,ubT(j),G.Nc);
  minv(j,:)= linspace(0,ubT(j),G.Nc);
end
tinv= [0,0.5,1];

S = struct('inc',inc,'inc_min',inc_min,'inc_wage',inc_wage,'h1',h1,'a1',a1,'extmin1',extmin1,'extmin2',extmin2,'c',c,'minv',minv,'tinv',tinv,...
            'dh',dh,'da',da,'hext',hext,'aext',aext,'h_',h_,'a_',a_,'b',b,'lb',lb,'ubT',ubT,'a2',a2,'h2',h2);

        
       
       
%%% Polynomial Bases and Derivatives %%%% 

[B]=chebpoly_base(G.ncheby+1,z); %Polynomial base for Objective function
T=kron(B,B);                


for t=1:1:G.ntime
[dBh(:,:,t),dBa(:,:,t),ddBh(:,:,t),ddBa(:,:,t)]=Dchebpoly_deriv(G.ncheby+1,z,B,dh(t),da(1,t));
end




%%% Polynomial Bases Policy Function %%%%
% Basis for states
B_sim=chebpoly_base(G.ncheb_pol,z);

% Basis for Income Shocks
z_epsw= 2*(G.eps_y(:,1)-G.eps_y(1,1))/(G.eps_y(G.Ne,1)-G.eps_y(1,1))-1; 
B_ew=chebpoly_base(G.Ne-1,z_epsw);
                                                                                                                     

% Multidimensional base
T_sim=kron(B_sim,kron(B_sim,B_ew));


B  = struct('B',B,'dBh',dBh,'dBa',dBa,'ddBh',ddBh,'ddBa',ddBa,...
            'B_sim',B_sim,'B_ew',B_ew,'T_sim',T_sim);
        


       

             
        


