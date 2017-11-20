%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB CODE MARRIAGE AND FERTILITY MODEL %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Structural Parameters

% EDUCATION STAGE
% %Family Background
% famb1=0.3; %delta1
% famb2=0.3; %delta2
% 
% %Parental Income Shocks by Age 14
% pincshocks=0.3; %delta3
% 
% %Private Expenditures in Schooling
% pexpschool=4; %delta4
% 
% %SD of Shocks
% sd_shock=sqrt(0.5);

% ADULT STAGE

%Disutility of work by Sector (1-2)
psi_r=-0.6;
psi_n=-0.8;

%Value of Marriage in HH Production (3-5)
theta1_r=0.2;
theta1_n=0.25;
theta1_u=0.4;

% %Value of Child HC in HH Production (6-8)
theta3_r=0.3;
theta3_n=0.3;
theta3_u=0.3;
theta=[theta1_r;theta1_n;theta1_u;theta3_r;theta3_n;theta3_u];

%Child HC Production Function (table 3, CGST_Oct82015)  (9-10)
gamma1=0.5595;
phi=0.4282;

% %Female Share of Consumption
% delta1=0.3;
% delta2=0.3;

%Family Background Types (11-12)
alpha0_r=0.0995;
alpha0_n=0.0115;

%Return to 2yr College (13-14)
alpha11_n=0.116;
alpha11_r=0.016;

%Return to General Experiencce (15-16)
alpha12_n=0.474;
alpha12_r=0.174;

%Return to Recent Sector Experience (17-18)
alpha2_r=0.437;
alpha2_n=0.302;
%unemployed?
alpha=[alpha0_r;alpha0_n;alpha11_r;alpha11_n;alpha12_r;alpha12_n;alpha2_r;alpha2_n];

% Shocks (19-21)
sigma_r = 0.43; %shock to regular
sigma_n = 0.73; %shock to non-regular
sigma_i = 245; %shock to hh income 245

% Probability of marriage (22-26)

omega0_w  =  0.3349; 
omega0_u  =  0.401; 
omega11 = - 0.0581;
omega12 = - 0.0904;
omega2 = 0.0152;
omega=[omega0_w;omega0_u;omega11;omega12;omega2];

% Terminal Value function (27-29)
lambda1=1;
lambda2=1;
lambda3=1;
lambda=[lambda1;lambda2;lambda3];
%Vector of Initial Parameters
params0 = [psi_r;psi_n;gamma1;phi;theta;alpha;sigma_r;sigma_n;sigma_i;omega;lambda];

% %Bounds
% alpha_b=[0.1 0.8];
% delta_b=[0.1 0.8];
% gamma_b=[0.2 0.9];
% theta_b=[0.1 0.9];
% bounds = [alpha_b,delta_b,gamma_b,theta_b];

%% Initial Conditions

edu_levels = [1:3];
abi_levels = [1 2];
types = [kron(abi_levels',ones(length(edu_levels),1)) repmat(edu_levels',[length(abi_levels) 1])];

%% General Parameters
Ne = 3; % Gauss-Hermite Points
%Utility of Consumption
sigma=0.4;
%Discount rate
beta=0.95;
%Interest rate
r=0.07;
%Investment in Children
inv=1;
% state parameters
n_incond = length(types);
n_shocks = 27;
n_period = 20;
n_pop = 1000;

G = struct('Ne',Ne,'sigma',sigma,'beta',beta,'r',r,'inv',inv,'n_incond',n_incond,'n_period',n_period,'n_shocks',n_shocks,'n_pop',n_pop);


%% Test Functions
tic;
S = sspace(params,G);
for z=1:n_incond
    [C(:,:,:,z),R(:,:,:,z),N(:,:,:,z),U(:,:,:,z),M(:,:,:,z)] = solution(G,types(z,1),types(z,2),S); 
end
toc;
[alpC(:,:,:,z),alpR(:,:,:,z),alpT(:,:,:,z),alpU(:,:,:,z),alpM(:,:,:,z)]=thci_polfunc(C(:,:,:,z),R(:,:,:,z),N(:,:,:,z),U(:,:,:,z),M(:,:,:,z),S,G);



 
%% save output
%save solutiontest.mat;
