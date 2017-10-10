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
            
            % household income:
            w_j = exp(alpha1*edu + alpha2*SS(j,3) + abi);
            hhinc = w_j + wh_j*m_j + shock;
            
            % transitions:
            X_0 = X_j - 1;
            wh_0 = wh_j;
            K_0 = K_j - inv;
            A_0 = A_j/(1+r) - hhinc - c_hh - n_j*inv; % eq. 8
            
            % probabilities:      
            
            % loop over consumption (5):
            for k = 1:1:length(c_vector)
                k
                c = c_vector(k);
                
                % value functions:
                
                % calculate utility
                u(k) = (c^(1-sigma))/(1-sigma) + phi + kappa*(1+theta1*SS(j,1)+theta2*SS(j,2)+theta3*SS(j,6));
                
                V(k) = u(k) +
            end
                % save optimal U* & c*
                [U_star,c_star] = max(V);
        end   
    end
end
