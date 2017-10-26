
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Structural Parameters

%Family Background
famb1=0.3; %delta1
famb2=0.3; %delta2

%Parental Income Shocks by Age 14
pincshocks=0.3; %delta3

%Private Expenditures in Schooling
pexpschool=4; %delta4

%SD of Shocks
sd_shock=sqrt(0.5);

%Utility of Consumption
sigma=0.4;

%Disutility of work by Sector
psi_r=0.4;
psi_n=0.4;
%unemployed?

%MRS of HH production wrt consumption/leisure
kappa=0.5;

%Value of Marriage in HH Production
theta1=0.7;

%Value of Number of Children in HH Production
theta2=0.3;

%Value of Child HC in HH Production
theta3=0.4;

%Child HC Production Function
gamma1=0.3;
gamma2=0.2;
gamma3=0.2;
phi=0.4;
rho=0.3;

%Female Share of Consumption
delta1=0.3;
delta2=0.3;

%Variance of HH Income Shocks
var_hhincshock=0.2;

%Family Background Types
alpha_f1=4;
alpha_f2=4;

%Return to College
alpha1=0.3;

%Return to General Experiencce
alpha2=0.3;

%Return to Recent Sector Experience
alpha3_r=0.4;
alpha3_n=0.4;
%unemployed?

%Interest rate
r=0.07;

%Discount rate
beta=0.05;

%Investment in Children
inv=1;

%Vector of Initial Parameters
Params0 = [famb1,famb2,pincshocks,pexpschool,sd_shock,...
           sigma,psi_r,psi_n,kappa,theta1,theta2,theta3,...
           gamma1,gamma2,gamma3,phi,rho,delta1,delta2,...
           var_hhincshock,alpha_f1,alpha_f2,alpha1,alpha2,...
           alpha3_r,alpha3_n];

%Bounds
alpha_b=[0.1 0.8];
delta_b=[0.1 0.8];
gamma_b=[0.2 0.9];
theta_b=[0.1 0.9];
bounds = [alpha_b,delta_b,gamma_b,theta_b];

P = struct();

%% Initial Conditions

edu_levels = [1:3];
abi_levels = [1 2];
types = [kron(abi_levels',ones(length(edu_levels),1)) repmat(edu_levels',[length(abi_levels) 1])];

%% Shocks

Ne = 3; % Gauss-Hermite Points
[e, wt] = GaussHermite(Ne);
sigma_r = sqrt(0.05); %shocks to regular
sigma_n = sqrt(0.05);
sigma_i = sqrt(0.05);

eps_r = sqrt(2)*e*sigma_r; % error vector
eps_n = sqrt(2)*e*sigma_v;
eps_i =

% DON'T NEED THIS NOW, THEY'RE INDEPENDENT
% vcv = [sigma_eps^2,0;0,sigma_v^2];
% detV = det(vcv);
% detV = det(vcv);
% detR = det(R);

% shock_r = [-1:1]; % sector R wages shock
% shock_n = [-1:1]; % sector N wages shock
% shock_hh= [-1:1]; % household income shock
shocks = [kron(shock_hh',ones(length(shock_n)*length(shock_r),1))...
                repmat(kron(shock_r',ones(length(shock_n),1)),[length(shock_hh) 1])...
                repmat(shock_n',[length(shock_hh)*length(shock_r) 1])]; % CHECK THIS w/Italo's board

%% General Parameters

n_incond = length(types);
n_shocks = length(shocks);
n_period = 25;
n_pop = 1000;

G = struct('n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop);

%% QUESTIONS:
%PARAMETERS BY AGE AND EDUCATION
%PROBABILITY FUNCTIONS

%% State Space

%Endogenous

ass_lb = 1;
ass_up = 5;
n_ass = 3;
assets = linspace(ass_lb,ass_up,n_ass);

matstat = [0 1];
n_matstat = length(matstat);

workexp = [1:4];
workexp_r = [1:3];
workexp_n = [1:3];
n_wrkexp = length(workexp);

sector = [1:3];

c_vector = [1:3];

%Exogenous

children = [0 1];

hearn_lb = 1;
hearn_ub = 5;
n_hearn = 3;
hearnings = linspace(hearn_lb,hearn_ub,n_hearn);

childHC_lb = 1.1;
childHC_ub = 4.5;
n_childHC = 3;
childHC = linspace(childHC_lb,childHC_ub,n_childHC);

SS_K = repmat(childHC',[length(assets)*length(hearnings)*length(workexp)*length(matstat) 1]);
SS_A = repmat(kron(assets',ones(length(childHC),1)),[length(hearnings)*length(workexp)*length(matstat) 1]);
SS_H = repmat(kron(hearnings',ones(length(assets)*length(childHC),1)),[length(workexp)*length(matstat) 1]);
SS_X = repmat(kron(workexp',ones(length(hearnings)*length(assets)*length(childHC),1)),[length(matstat) 1]);
SS_N = kron(children',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));
SS_M = kron(matstat',ones([length(childHC)*length(assets)*length(hearnings)*length(workexp),1]));

SS = [SS_M SS_N SS_X SS_H SS_A SS_K];
