% loop for initial conditions:
% for n = 1:1:G.n_incond
n=1;  
     abi = types(n,1);
     edu = types(n,2);

% Terminal Value Function:
% TVF = A_T + W_T + Q_T = assets + wages + HH_production

assets = SS_A; % just the vector of all possible assets, not a function
wages = exp(alpha1*edu + alpha2*SS_X + abi);
hhprod = (SS_M+1).^theta1 .* (SS_N+1).^theta2 .* SS_K;

TVF = assets + wages + hhprod;

% loop for time (25):
for t = G.n_period-1:-1:1
    t
    
    % loop for shocks (27):
    for i = 1:1:G.n_shocks
        i
        
        shock_hh= shock_vector(1,i);
        shock_r = shock_vector(2,i);
        shock_n = shock_vector(3,i);
        
        % loop over states (216):
        for j = 1:1:length(SS)
            j
            
            m_j = SS_M(j);  % marital status
            n_j = SS_N(j);  % children
            X_j = SS_X(j);  % experience
            wh_j = SS_H(j); % husband's wage
            A_j = SS_A(j);  % HH's assets
            K_j = SS(j,6);  % HC of child
                
            % transitions for exogenous variables
            wh_0 = wh_j;
            K_0 = K_j + inv;
            m_0 = m_j; % prob of marriage
            n_0 = n_j; % prob of marriage
              
            % probabilities:      
            lambda_r = 0.7; %needs to be function of education & ability
            pi_r = 0.2;
            
            % loop over consumption (5):
            for k = 1:1:length(c_vector) % this is HH consumption
                k
                
                chh = c_vector(k); % HH consumption
                cw = (1+delta1*m_j+delta2*n_j)*chh; % woman's consumption
                
                % transitions:
                    
                    % regular job:
                    X_0 = X_j + 1;
                    w_j = exp(alpha1*edu + alpha2*X_j + abi + shock_r);
                    A_0 = A_j/(1+r) - (w_j + wh_j*m_j + shock_hh) - chh - n_j*inv; % eq. 8 CHECK for NEXT Period
                     
                    V_r_next = interpn(SS_M,SS_N,SS_X,SS_H,SS_A,SS_K,TVF,m_0,n_0,X_0,wh_0,A_0,K_0);
                   
                    % non-regular job:
                    X_0 = X_j + 1;
                    w_j = exp(alpha1*edu + alpha2*SS(j,3) + abi + shock_n);
                    A_0 = A_j/(1+r) - (w_j + wh_j*m_j + shock_hh) - chh - n_j*inv; % CHECK for NEXT Period
                    
                    V_n_next = interpn(SS_M,SS_N,SS_X,SS_H,SS_A,SS_K,TVF,m_0,n_0,X_0,wh_0,A_0,K_0);
                    
                    % unemployment:
                    X_0 = X_j;
                    w_j = w_min;
                    A_0 = A_j/(1+r) - (w_j + wh_j*m_j + shock_hh) - chh - n_j*inv;
                    
                    V_u_next = interpn(SS_M,SS_N,SS_X,SS_H,SS_A,SS_K,TVF,m_0,n_0,X_0,wh_0,A_0,K_0);
                     
                % calulcate the wage in regular and non-regular (using the
                % experience) and the budget constraint (using the assets
                % transtion)
                % start with sector specific value functions
                % value functions:
                    [I_r, Emax_r] = max(V_r_next,V_n_next);
                    [I_n, Emax_n] = max(V_r_next,V_n_next,V_u_next);
                    
                    % calculate utility
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + psi_r + kappa*(1+theta1*SS(j,1)+theta2*SS(j,2)+theta3*SS(j,6));
                    V_r(k) = u_r(k) + beta*lambda_r*Emax_r +;
                    % calculate utility
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + phi_n + kappa*(1+theta1*SS(j,1)+theta2*SS(j,2)+theta3*SS(j,6));
                    V_r(k) = u_r(k) + 

            end
                % save optimal U* & c*
                [U_star,c_star] = max(V);
        end   
    end
end