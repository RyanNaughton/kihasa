% loop for initial conditions:
% for n = 1:1:G.n_incond
n=1;
     abi = types(n,1);
     edu = types(n,2);

% Terminal Value Function:
% TVF = A_T + W_T + Q_T = assets + wages + HH_production

assets = SS_A; 
wages = exp(alpha1*SS_X + alpha2*edu + abi);
hhprod = (SS_M+1).^theta1 .* (SS_N+1).^theta2 .* SS_K;
% matrix of J=10x5x5=500 rows and 10x2=20 cols
TVF = assets + wages + hhprod;

tic
% loop for time (25):
for t = G.n_period-1:-1:1
    t
    toc
    if t==G.n_period-1
        Emax = TVF;
    else
        Emax = W;
    end
    
    % Chevyshev Approximation
    B = kron(T_A, kron(T_H,T_K)); % half of the numerator of the coeff
    Num = TVF'*B; % TVF'*kron(T_A, kron(T_H,T_K)); % numerator (bases*function) 
    Den = kron(T2_A, kron(T2_H,T2_K)); % square of T multiplied
    for x = 1:1:(n_matstat*n_wrkexp)
        alpha(x,:) = Num(x,:)./Den';
    end
    % alpha contains 20 rows of 19x4x4 = 304 coefficients
    
    % maybe won't need this
    assets_wide = linspace(-5,100,n_assets); % probably going to change by period
    childHC_wide = linspace(-5,10,n_childK);
    workexp_wide = linspace(-1,10,n_wrkexp);    

    % loop for shocks (27):
    for i = 1:1:G.n_shocks % 27 x 3
        i
        
        shock_i = shocks_i(i);
        shock_r = shocks_r(i);
        shock_n = shocks_n(i);
        
        % loop for work experience and marital status (20):
        for x = 1:1:(n_matstat*n_wrkexp)
            x;
            coeffs = alpha(1,:);
            % current state variables:
            m_j = SS_M(x);  % marital status
            n_j = SS_N(x);  % children
            X_j = SS_X(x);  % experience
         
        % loop over states (20 assets x 5 child HC x 5 husband wages = 500):
        for j = 1:1:length(SS_rows)
            j;
            
            % current state variables:
            wh_j = SS_H(j); % husband's wage
            A_j = SS_A(j);  % HH's assets
            K_j = SS_K(j);  % HC of child
            
            % sector-specific state variables:
            w_j_r = exp(alpha1*X_j + alpha2*edu + abi + shock_r); % same
            w_j_n = exp(alpha1*X_j + alpha2*edu + abi + shock_n); % same
            w_j_u = 0; % unemployed women don't have earnings?
              
            % probabilities:      
            lambda_r = (edu + X_j + abi)/10; % currently has abi, shouldn't
            pi_r = 0.2; % one Emax for temporary/unemployed and one Emax for all
            
            % sector-specific probabilities: % CHECK THIS WHEN AGE VARIES
            prob_marr_r = (edu + abi + 35 + t)/100; % only change with time & type
            prob_marr_n = (edu + abi + 30 + t)/100; % only change with time & type
            prob_marr_u = (edu + abi + 25 + t)/100; % only change with time & type
            probably = a0 + a1*col4 + a2*col2 + a3*work + a4*t;
            
            % transitions for exogenous variables:
            wh_next = wh_j;
            K_next = K_j + inv; % make a function (CES)
            
            % save wages in R and N for all states and all shocks
            % wages_r_allstates = 216 x 23
            
            % make consumption vector HERE (for each time-shock-state)
            % from 0 to max potential assets
            %c_vector = ;
    %%%%%%%%% this is how you adjust the linepsace
%     ubT= inc(1,3) + (1+P.r)*a_(:,3);
%     for j=1:G.nss
%       c(j,:)= linspace(0,ubT(j),G.Nc);
%       minv(j,:)= linspace(0,ubT(j),G.Nc);
%     end
%     tinv= [0,0.5,1];

            % loop over consumption (5):
            for k = 1:1:length(c_vector)
                k;
                
                chh = c_vector(k); % HH consumption
                cw = (1+delta1*m_j+delta2*n_j)*chh; % woman's consumption - can this be larger than HH consumption?
                
                % single:
                if x <= 10
                    
                    % regular job:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh - n_j*inv); % eq. 8
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
                    
                    % value function:
                    V_r_next = interpn(childHC_wide,assets_wide,hearnings,workexp_wide,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);

                    % non-regular job:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh - n_j*inv);
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
                    
                    % value function:
                    V_n_next = interpn(childHC_wide,assets_wide,hearnings,workexp_wide,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);
                    
                    % unemployed:
                    X_next = X_j + 0;
                    A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh - n_j*inv);
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
                    
                    % value function:
                    V_u_next = interpn(childHC_wide,assets_wide,hearnings,workexp_wide,matstat, tmp, K_next,A_next,wh_next,X_next,m_next);
                    
                else   
                % Married:
                
                    % Regular:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh - n_j*inv);
                    m_next = m_j;
                    n_next = n_j;
                    Vm_r_next = ;
                    
                    % Non-regular:
                    X_next = X_j + 1;
                    A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh - n_j*inv);
                    m_next = m_j;
                    n_next = n_j;
                    Vm_n_next = ;
                    
                    % Unemployed:
                    X_next = X_j + 0;
                    A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh - n_j*inv);
                    m_next = m_j;
                    n_next = n_j;
                    Vm_u_next = ;
                    
                end
                
                % Together:
                    
                    % Sector-Specific Utility:
                    u_r(k) = (cw^(1-sigma))/(1-sigma) + psi_r + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_n(k) = (cw^(1-sigma))/(1-sigma) + psi_n + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_u(k) = (cw^(1-sigma))/(1-sigma) + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    
                    % Expected Utility:
                    EU = max([V_r_next,V_n_next,V_u_next,Vm_r_next,Vm_n_next, Vm_u_next]);
                    
                    % Sector-Specific Value Functions:
                    V_r(k) = u_r(k) + beta * EU;
                    V_n(k) = u_n(k) + beta * EU;
                    V_u(k) = u_u(k) + beta * EU;         
            end
            
                % save optimal V* & c*
                %%%validate that consumption is feasible
                [V_r_star, Index_r_k] = max(V_r);
                [V_n_star, Index_n_k] = max(V_n);
                [V_u_star, Index_u_k] = max(V_u);
                c_r_star = c_vector(Index_r_k);
                c_n_star = c_vector(Index_n_k);
                c_u_star = c_vector(Index_u_k);
                c_star_vector = [c_r_star, c_n_star, c_u_star];
                
                % save labor choice:
                [V_star, Index_l] = max([V_r_star,V_n_star,V_u_star]);
                c_star(j, i, t) = c_star_vector(Index_l);
                l_star(j, i, t) = Index_l;
                V_star_aux(j, i, t) = V_star; 
                
        end
        
    end
    
    % Integrate over shocks
%     W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt; 
%     Em(:,t,k) = detR*detV^(-1/2)*pi^(-(size(R,1))/2)*G.w(i)*G.w(j)*V + Em(:,t,k);
    W = pi^(-1/2)*V_star_aux(:,:,t)*weight;   
end
