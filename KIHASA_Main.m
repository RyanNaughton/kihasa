
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
r=0.7;

%Ivestment in Children
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

%PARAMETERS BY AGE AND EDUCATION
%PROBABILITY FUNCTIONS

P = struct();

%% Initial Conditions

edu = [1:3];
abi = [1 2];
types = [kron(abi',ones(length(edu),1)) repmat(edu',[length(abi) 1])];

%% Shocks
shock_r = [1:3];
shock_n = [1:3];
shocks = [kron(shock_r',ones(length(shock_n),1)) repmat(shock_n',[length(shock_r) 1])];

shock_vector = [-1:1]; % household income shock

%% General Parameters

n_incond = length(types);
n_shocks = length(shocks); %should be 27, missing 1 (unemployment) shock?
n_period = 25;
n_pop = 1000;

G = struct('n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop);

%% State Space

%Endogenous

ass_lb = 1;
ass_up = 5;
n_ass = 10;
assets = linspace(ass_lb,ass_up,n_ass);

matstat = [0 1];

workexp = [1:15];
workexp_r = [1:3];
workexp_n = [1:3];

sector = [1:3];

c_vector = [1:3];

%Exogenous

children = [0 1];

hearn_lb = 1;
hearn_ub = 5;
n_hearn = 5;
hearnings = linspace(hearn_lb,hearn_ub,n_hearn);

childHC_lb = 1.1;
childHC_ub = 4.5;
n_childHC = 5;
childHC = linspace(childHC_lb,childHC_ub,n_childHC);

SS = [kron(matstat',ones(length(sector),1)) repmat(sector',[length(matstat) 1])];

%% Temporarily, import SS from excel

[~, ~, raw] = xlsread('C:\Users\jchenper\Documents\GitHub\kihasa\ini_ss.xlsx','Sheet1','A2:F217');
iniss = reshape([raw{:}],size(raw));
clearvars raw;
SS = iniss;
