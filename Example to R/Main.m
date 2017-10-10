%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Parental Investments in Children  n%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model v1: Model with investments of time and money
% - Assumptions: 
%   * h(0)=0; a(0) random uniformly distributed
%   * flat income profile with iid shocks
%   * technology shocks
% - Output: c_sim, t_sim, m_sim, h_sim, a_sim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load mm_est1.txt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 1. PARAMETERS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Structural Estimates %%%

gamma     =  [0.5991;0.4602;0.6194;0.20045;0.2699;0.1903];
phi       = 0.4282; 
rho       = 0.7487; 
delta     = 2.23;
     
params0=[gamma;phi;rho;delta];


gamma_b = [0.1,0.8;0.3,0.9;0.25,0.9;0.0,0.5;0.0,0.5;0.0,0.5];
phi_b   = [-1.2,0.85];
rho_b   = [0.2,0.85];
delta_b     = [1,4];
bounds  = [gamma_b;phi_b;rho_b;delta_b];
   
%%% Calibrated %%%

sigma     = 0.5; %Elasticity of Substitution
beta      = 0.96;
sigma_eps  = sqrt(0.05); %shocks to income
sigma_v   = sqrt(0);
eta       = 12; %terminal value parameter (human capital)
tau1      = 12; % terminal value parameter (assets) - If tau/eta too small V_t=1 is not concave in some cases
r         = 1/beta-1; %interest rate
p         = 1; %relative price of investment and consumption
w_min     = [0.5,0.5,0.5];  %% CHANGE THIS TO TEST TIMING OF INCOME
mu        = 0; %mean of shock(s)
cmin      = w_min;
u0        = 0;
u1        = 0.5;
u2        = -0.5;
vcv       =[sigma_eps^2,0;0,sigma_v^2];
detV      =det(vcv);
w         = [1,1,1];        


P  = struct('beta',beta,'sigma',sigma,'sigma_eps',sigma_eps,'eta',eta,'tau1',tau1,'r',r,'p',p,...
           'delta',delta,'w_min',w_min,'mu',mu,'cmin',cmin,'u0',u0,'u1',u1,'u2',u2,'w',w,'sigma_v',sigma_v);
     


%%% State space %%%

% Total states: A=15; H:15; eps_y=3; eps_h1=3 

M= 30; %Chevyshev points to evaluate objective function
M2=M;
Ne = 3; % Gauss-Hermite Points
nss = M^2; %state space A and H1
ncheby    = M-2;
npop      = 40000;
ntime     = 4;
nper      = 3;
Nc        = 20;
Ntinv     = 3;
ncheb_pol=6;
%pc=[5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95];
neps=2;
[e, wt]= GaussHermite_2(Ne); % Delivers integration nodes and weights for epsilon
eps_y=sqrt(2)*e*sigma_eps; % error vector
eps_h=sqrt(2)*e*sigma_v;


G  = struct('M',M,'M2',M2,'Ne',Ne,'nss',nss,'ncheby',ncheby,'npop',npop,...
           'ntime',ntime,'nper',nper,'ncheb_pol',ncheb_pol,'Nc',Nc,'Ntinv',Ntinv,...
           'eps_y',eps_y,'eps_h',eps_h,'wt',wt);

       
%%% Shocks %%%

epsy_sim=randn(npop,nper).*sigma_eps; 
epsh_sim=randn(npop,nper).*sigma_v;        
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4. Test Functions   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S,B,inc_min_sim,inc_wage_sim,inc_sim,~,PI]=thci_ss(params0,G,P,epsy_sim);

tic;   
[W,A,Minv,Tinv,C]=thci_sol_grid20(params0,P,G,B,S);   
toc;

tic;
[alp_C,alp_Minv,alp_Tinv]=thci_polfunc_new(G,P,B,S,C,Minv,Tinv);
toc;

tic;
[c_s,minv_s,tinv_s,m_s,t_s,a_s,h_s] = thci_sim(params0,alp_C,alp_Minv,alp_Tinv,G,P,S,epsy_sim,inc_min_sim,inc_wage_sim,epsh_sim);
toc;

simula_data = thci_plots(c_s,m_s,t_s,a_s,h_s,inc_sim,P,G);

save solJun17.mat


%%
%% load('optCMT_may15.mat')
%%

copt_t3_e1=reshape(C(:,1,3),30,30)';
copt_t3_e2=reshape(C(:,2,3),30,30)';
copt_t3_e3=reshape(C(:,3,3),30,30)';
copt_t2_e1=reshape(C(:,1,2),30,30)';
copt_t2_e2=reshape(C(:,1,2),30,30)';
copt_t2_e3=reshape(C(:,3,2),30,30)';
copt_t1_e1=reshape(C(:,1,1),30,30)';
copt_t1_e2=reshape(C(:,1,1),30,30)';
copt_t1_e3=reshape(C(:,3,1),30,30)';
x=(1:30);
plot(x,copt_t3_e1(1,:),x,copt_t3_e1(30,:));
plot(x,copt_t3_e2(1,:),x,copt_t3_e2(30,:));
plot(x,copt_t3_e3(1,:),x,copt_t3_e3(30,:));
plot(x,copt_t2_e1(1,:),x,copt_t2_e1(30,:));
plot(x,copt_t2_e2(1,:),x,copt_t2_e2(30,:));
plot(x,copt_t2_e3(1,:),x,copt_t2_e3(30,:));
plot(x,copt_t1_e1(1,:),x,copt_t1_e1(30,:));
plot(x,copt_t1_e2(1,:),x,copt_t1_e2(30,:));
plot(x,copt_t1_e3(1,:),x,copt_t1_e3(30,:));







[alp_A,alp_I]=thci_polfunc(G,P,B,S,H,A,I);
[c_s,i_s,a_s,h_s] = thci_sim(params0,alp_A,alp_I,G,P,S,epsy_sim,inc_sim,epsh_sim);
toc;    




  
matlabpool close



% 
% [theta,fval,exit,output,lambda]=ktrlink(@(theta) thci_mdf(theta,P,G,eps_sim,mm_est1),params0,[],[],[],[],bounds(:,1),bounds(:,2),[],[],'knitro_est.opt');
% theta
% save thci.mat;
% 
% matlabpool close
exit;



