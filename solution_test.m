% loop for initial conditions:
% for n = 1:1:G.n_incond
n=1;  
     abi = types(n,1);
     edu = types(n,2);

% Terminal Value Function:
% TVF = A_T + W_T + Q_T = assets + wages + HH_production

assets = SS(:,5); % just the vector of all possible assets, not a function
wages = exp(alpha1*edu + alpha2*SS(:,3) + abi);
hhprod = (SS(:,1)+1).^theta1 .* (SS(:,2)+1).^theta2 .* SS(:,6);

TVF = assets + wages + hhprod;

% loop for time (25):
for t = G.n_period-1:-1:1
    t
    % compute coefficients:
    
    % loop for shocks (9, should be 27):
    for i = 1:1:G.n_shocks
        i
        
        shock = shock_vector(i);
        
        % loop over states (216):
        for j = 1:1:length(SS)
            j
            
            % marital status:
            m_j = SS(j,1);
            
            % children:
            n_j = SS(j,2);
            
            % experience:
            X_j = SS(j,3);
            
            % husband wage:
            wh_j = SS(j,4);
            
            % assets:
            A_j = SS(j,5);
            
            % HC of child:
            K_j = SS(j,6);
                  
            % transitions for exogenous variables
            wh_0 = wh_j;
            K_0 = K_j - inv;
              
            % probabilities:      
            lambda_r = 0.7;
            pi_r = 0.2;
            
            % loop over consumption (5):
            for k = 1:1:length(c_vector) % this is HH consumption
                k
                
                chh = c_vector(k);
                cw = (1+delta1*m_j+delta2*n_j)*chh;
                
                % value functions:
                
                    % transitions for endogenous variables (assets and work
                    % experience)
                    
                    % regular job:
                     X_0 = X_j - 1;
                     w_j = exp(alpha1*edu + alpha2*SS(j,3) + abi);
                     A_0 = A_j/(1+r) - (w_j + wh_j*m_j + shock) - chh - n_j*inv; % eq. 8
                     
                     V_r_next = ;
                     V_n_next = ;
                     V_u_next = ;
                     [I_r, Emax_r] = max(V_r_next,V_n_next);
                     [I_n, Emax_n] = max(V_r_next,V_n_next,V_u_next);
                    
                    % calculate utility
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + phi + kappa*(1+theta1*SS(j,1)+theta2*SS(j,2)+theta3*SS(j,6));
                    V_r(k) = u_r(k) + beta*lambda_r*Emax_r +;
                    
                    % non-regular job:
                    X_0 = X_j - 1;
                    w_j = exp(alpha1*edu + alpha2*SS(j,3) + abi); % how is this different from regular?
                    A_0 = A_j/(1+r) - (w_j + wh_j*m_j + shock) - chh - n_j*inv;
                    
                    % calculate utility
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + phi + kappa*(1+theta1*SS(j,1)+theta2*SS(j,2)+theta3*SS(j,6));
                    V_r(k) = u_r(k) + 
                    
                    
                    % unemployment:
                    X_0 = X_j;
                    w_j = w_min;
                     
                % calulcate the wage in regular and non-regular (using the
                % experience) and the budget constraint (using the assets
                % transtion)
                
                % start with sector specific value functions
                

            end
                % save optimal U* & c*
                [U_star,c_star] = max(V);
        end   
    end
end
