function [S] = sspace(params,G)

%% Shocks
[e, wt] = GaussHermite(G.Ne);
eps_r = sqrt(2)*e*params(19); % error vector
eps_n = sqrt(2)*e*params(20); % error vector
eps_i = sqrt(2)*e*params(21); % error vector

% DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% vcv = [sigma_eps^2,0;0,sigma_v^2];
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

shocks_i = kron(eps_i,ones(length(eps_n)*length(eps_r),1));
shocks_r = repmat(kron(eps_r,ones(length(eps_n),1)),[length(eps_i) 1]);
shocks_n = repmat(eps_n,[length(eps_i)*length(eps_r) 1]);

shocks = [shocks_i shocks_r shocks_n];
weight = kron(wt, kron(wt,wt)); % 27x1

%% State Space

%Endogenous

assets_lb = 1;
assets_ub = 5;
n_assets = 20;
assets = linspace(assets_lb,assets_ub,n_assets);

matstat = [1 0];
n_matstat = length(matstat);

workexp = [1:10];
workexp_r = [1:3];
workexp_n = [1:3];
n_wrkexp = length(workexp);

sector = [1:3];

% c_min = 0;
% c_max = 5;
c_n = 20;
% c_vector = linspace(c_min,c_max,c_n);

%Exogenous

children = [1 0];

hwages_lb = 1;
hwages_ub = 5;
n_hwages = 5;
hwages = linspace(hwages_lb,hwages_ub,n_hwages);

childK_lb = 1.1;
childK_ub = 4.5;
n_childK = 5;
childK = linspace(childK_lb,childK_ub,n_childK);

%% SS for linear

% SS_K = repmat(childHC',[length(assets)*length(hearnings)*length(workexp)*length(matstat) 1]);
% SS_A = repmat(kron(assets',ones(length(childHC),1)),[length(hearnings)*length(workexp)*length(matstat) 1]);
% SS_H = repmat(kron(hearnings',ones(length(assets)*length(childHC),1)),[length(workexp)*length(matstat) 1]);
% SS_X = repmat(kron(workexp',ones(length(hearnings)*length(assets)*length(childHC),1)),[length(matstat) 1]);
% SS_N = kron(children',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));
% SS_M = kron(matstat',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));

%SS = [SS_M SS_N SS_X SS_H SS_A SS_K];

%% SS for chevyshev

SS_K = repmat(childK',[length(assets)*length(hwages) 1]);
SS_A = repmat(kron(assets',ones(length(childK),1)),[length(hwages) 1]);
SS_H = repmat(kron(hwages',ones(length(assets)*length(childK),1)), 1 );
SS_rows = [SS_H SS_A SS_K]; % rows

SS_X = repmat(workexp, [1 length(matstat)]);
SS_M = kron(matstat, ones([1, length(workexp)]));
SS_N = kron(children, ones([1, length(workexp)]));
SS_cols = [SS_M' SS_N' SS_X']; % columns

%% Chevyshev Approximation

[z_A,ext_A,extmin_A,extmax_A,d_A,vector_A,T_A,T2_A] = cheby_values(n_assets,assets_ub,assets_lb);
[z_H,ext_H,extmin_H,extmax_H,d_H,vector_H,T_H,T2_H] = cheby_values(n_hwages,hwages_ub,hwages_lb);
[z_K,ext_K,extmin_K,extmax_K,d_K,vector_K,T_K,T2_K] = cheby_values(n_childK,assets_ub,assets_lb);

% Basis for Income Shocks
% zeps_r= 2*(eps_r-eps_r(1))/(eps_r(Ne,1)-eps_r(1))-1; 
% zeps_n= 2*(eps_n-eps_n(1))/(eps_n(Ne,1)-eps_n(1))-1; 
% zeps_i= 2*(eps_i-eps_i(1))/(eps_i(Ne,1)-eps_i(1))-1; 
% Teps_r=chebpoly_base(Ne-1,zeps_r);
% Teps_n=chebpoly_base(Ne-1,zeps_n);
% Teps_i=chebpoly_base(Ne-1,zeps_u);
% T2eps_r = diag(Teps_r'*Teps_r);
% T2eps_n = diag(Teps_n'*Teps_n);
% T2eps_i = diag(Teps_i'*Teps_i);   

S = struct('Ne',Ne,'sigma',sigma,'beta',beta,'r',r,'inv',inv,'n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop);



end
