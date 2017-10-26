% loop for initial conditions:
% for n = 1:1:G.n_incond
n=1;
     abi = types(n,1);
     edu = types(n,2);

% Terminal Value Function:
% TVF = A_T + W_T + Q_T = assets + wages + HH_production

assets = SS_A; % just the vector of all possible assets, not a function
wages = exp(alpha1*SS_X + alpha2*edu + abi);
hhprod = (SS_M+1).^theta1 .* (SS_N+1).^theta2 .* SS_K;

TVF = assets + wages + hhprod;
%W(:,G.ntime)=P.eta*S.h_(:,G.ntime).^(1-P.sigma)/(1-P.sigma)+
%P.tau1*(1-exp(-S.a_(:,G.ntime)));

% loop for time (25):
%for t = G.n_period-1:-1:1
    t=1;
    
    if t==G.n_period-1
        Emax = TVF;
    end
    
    assets_wide = linspace(-5,10,n_ass);
    childHC_wide = linspace(-5,10,n_childHC);
    
    %%%%%%%%% this is how you adjust the linepsace
%     ubT= inc(1,3) + (1+P.r)*a_(:,3);
%     for j=1:G.nss
%       c(j,:)= linspace(0,ubT(j),G.Nc);
%       minv(j,:)= linspace(0,ubT(j),G.Nc);
%     end
%     tinv= [0,0.5,1];
    
    % loop for shocks (27):
    %for i = 1:1:G.n_shocks % 27 x 3
        i=1;
        
        shock_hh= shocks(1,i);
        shock_r = shocks(2,i);
        shock_n = shocks(3,i);
    tic    
        % loop over states (216):
        for j = 1:1:length(SS)
            j
            
            % current state variables:
            m_j = SS_M(j);  % marital status
            n_j = SS_N(j);  % children
            X_j = SS_X(j);  % experience
            wh_j = SS_H(j); % husband's wage
            A_j = SS_A(j);  % HH's assets
            K_j = SS(j,6);  % HC of child
            
            % sector-specific state variables:
            w_j_r = exp(alpha1*X_j + alpha2*edu + abi + shock_r); % same
            w_j_n = exp(alpha1*X_j + alpha2*edu + abi + shock_n); % same
            w_j_u = 0; % unemployed women don't have earnings?
              
            % probabilities:      
            lambda_r = (edu + X_j + abi)/10; % currently has abi, shouldn't
            pi_r = 0.2;
            
            % sector-specific probabilities:
            prob_marr_r = (edu + abi + 3 + t)/10; % only change with time & type
            prob_marr_n = (edu + abi + 2 + t)/10; % only change with time & type
            prob_marr_u = (edu + abi + 1 + t)/10; % only change with time & type
            
            % transitions for exogenous variables:
            wh_next = wh_j;
            K_next = K_j + inv;
            
            % loop over consumption (5):
            for k = 1:1:length(c_vector)
                k
                
                chh = c_vector(k); % HH consumption
                cw = (1+delta1*m_j+delta2*n_j)*chh; % woman's consumption - can this be larger than HH consumption?
                
                % regular job:
                
                    % transitions:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_hh) - chh - n_j*inv); % eq. 8
                    if prob_marr_r > 0.5
                        m_next = m_j + 1;
                        n_next = n_j + 1;
                    else
                        m_next = m_j;
                        n_next = n_j;
                    end
                    if m_next > 1 % limit to less than 1
                        m_next = 1;
                        n_next = 1;
                    end
                    
                    % interpolated t+1:
                    tmp = reshape(TVF,[n_childHC,n_ass,n_hearn,n_wrkexp,n_matstat]);
                    V_r_next = interpn(childHC_wide,assets_wide,hearnings,workexp,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);
                    V_r_next
                % non-regular job:
                
                    % transitions:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_hh) - chh - n_j*inv);
                    if prob_marr_n > 0.5
                        m_next = m_j + 1;
                        n_next = n_j + 1;
                    else
                        m_next = m_j;
                        n_next = n_j;
                    end
                    if m_next > 1 % limit to less than 1
                        m_next = 1;
                        n_next = 1;
                    end
                    
                    % interpolated t+1:
                    tmp = reshape(TVF,[n_childHC,n_ass,n_hearn,n_wrkexp,n_matstat]);
                    V_n_next = interpn(childHC_wide,assets_wide,hearnings,workexp,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);
                    V_n_next
                % unemployed:
                    
                    % transitions:
                    X_next = X_j + 0;
                    A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_hh) - chh - n_j*inv);
                    if prob_marr_u > 0.5
                        m_next = m_j + 1;
                        n_next = n_j + 1;
                    else
                        m_next = m_j;
                        n_next = n_j;
                    end
                    if m_next > 1 % limit to less than 1
                        m_next = 1;
                        n_next = 1;
                    end
                    
                    % interpolated t+1:
                    tmp = reshape(TVF,[n_childHC,n_ass,n_hearn,n_wrkexp,n_matstat]);
                    V_u_next = interpn(childHC_wide,assets_wide,hearnings,workexp,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);
                    V_u_next
                % Value Functions:
                    
                    % Sector-Specific Utility:
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + psi_r + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_n(k) = (cw^(1-sigma))/(1-sigma) + psi_n + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_u(k) = (cw^(1-sigma))/(1-sigma) + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    
                    % Expected Utility:
                    Emax = max([V_r_next,V_n_next,V_u_next]);
                    
                    % Sector-Specific Value Functions:
                    V_r(k) = u_r(k) + beta * Emax;
                    V_n(k) = u_n(k) + beta * Emax;
                    V_u(k) = u_u(k) + beta * Emax;
                    
            end
                % save optimal V* & c*
                [V_r_star, Index_r_k] = max(V_r);
                [V_n_star, Index_n_k] = max(V_n);
                [V_u_star, Index_u_k] = max(V_u);
                c_r_star = c_vector(Index_r_k);
                c_n_star = c_vector(Index_n_k);
                c_u_star = c_vector(Index_u_k);
                c_star_vector = [c_r_star, c_n_star, c_u_star];
                
                % save labor choice:
                [V_star, Index_l] = max([V_r_star,V_n_star,V_u_star]);
                c_star(j) = c_star_vector(Index_l); % add i to have 216 x 3?
                l_star(j) = Index_l;
                V_star_aux(j) = V_star;
        end
     toc   
        % Integrate over shocks
        W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt;
        Em(:,t,k) = detR*detV^(-1/2)*pi^(-(size(R,1))/2)*G.w(i)*G.w(j)*V + Em(:,t,k);
    %end
%end