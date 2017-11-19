function [c_func, lr_func, ln_func, lu_func, m_func] = solution(G,abi,edu,S,shocks)

% loop for initial conditions:
% for n = 1:1:G.n_incond
% n=1;
%      abi = types(n,1);
%      edu = types(n,2);

% Terminal Value Function: TVF = A_T + W_T + Q_T = assets + wages + HH_prod

    assets = SS_A; 
    wages = exp(alpha1*SS_X + alpha2*edu + abi);
    hhprod = (SS_M+1).^theta1 .* (SS_N+1).^theta2 .* SS_K;
    % matrix of J=10x5x5=500 rows and 10x2=20 cols
    TVF = assets + wages + hhprod;
    
tic
% loop for time (20):
for t = G.n_period-1:-10:1
    t
    toc
    if t==G.n_period-1
        Emax = TVF;
        % Chevyshev Approximation - alpha contains 20 rows of 19x4x4 = 304 coefficients  
        Num = Emax'*kron(T_A, kron(T_H,T_K)); % numerator (bases*function) 
        Den = kron(T2_A, kron(T2_H,T2_K)); % square of T multiplied
        for x = 1:1:(n_matstat*n_wrkexp)
            alpha(x,:) = Num(x,:)./Den';
        end
    else
        Emax = W(t+1);
        % use 20 new VF (W) to get 20 new coefficients
        Num = Emax'*kron(T_A, kron(T_H,T_K));
        Den = kron(T2_A, kron(T2_H,T2_K));
        for x = 1:1:(n_matstat*n_wrkexp)
            alpha(x,:) = Num(x,:)./Den';
        end
    end
    
    % loop for work experience and marital status (20):
    for x = 1:10:length(SS_cols)
        x;
        
        % current state variables:
        m_j = SS_M(x);  % marital status
        n_j = SS_N(x);  % children
        X_j = SS_X(x);  % experience
    
        % loop for shocks (27):
        for i = 1:1:G.n_shocks % 27 x 3
            i;
        
            shock_i = shocks_i(i);
            shock_r = shocks_r(i);
            shock_n = shocks_n(i);
            
            % loop over continuous states (20 assets x 5 child HC x 5 hwages = 500):
            for j = 1:1:5 %00:length(SS_rows)
                j;
            
                % current state variables:
                wh_j = SS_H(j); % husband's wage
                A_j = SS_A(j);  % HH's assets
                K_j = SS_K(j);  % HC of child
            
                % sector-specific state variables:
                w_j_r = exp(alpha1*X_j + alpha2*edu + abi + shock_r); % same
                w_j_n = exp(alpha1*X_j + alpha2*edu + abi + shock_n); % same
                w_j_u = 0; % unemployed women don't have earnings?
              
                % sector-specific probabilities: % CHECK THIS WHEN AGE VARIES
                prob_marr_r = (edu + abi + 35 + t)/100; % only change with time & type
                prob_marr_n = (edu + abi + 30 + t)/100; % only change with time & type
                prob_marr_u = (edu + abi + 25 + t)/100; % only change with time & type
                %probability = a0*abi + a1*col4 + a2*col2 + a3*work + a4*t;
                %will be different
            
                % transitions for exogenous variables:
                K_next = (gamma1*K_j^phi + (1-gamma1)*inv^phi)^(1/phi); % make a function (CES)
                wh_next = wh_j; % no transition
            
                % consumption vector
                chh_r_min = A_j + (w_j_r + wh_j*m_j + shock_i) - n_j*inv - extmax_A/(1+r);
                chh_r_max = A_j + (w_j_r + wh_j*m_j + shock_i) - n_j*inv - extmin_A/(1+r);
                cr_vector = linspace(chh_r_min,chh_r_max,c_n);
                cr_vector(cr_vector < 0) = 0;
                chh_n_min = A_j + (w_j_n + wh_j*m_j + shock_i) - n_j*inv - extmax_A/(1+r);
                chh_n_max = A_j + (w_j_n + wh_j*m_j + shock_i) - n_j*inv - extmin_A/(1+r);
                cn_vector = linspace(chh_n_min,chh_n_max,c_n);
                cn_vector(cn_vector < 0) = 0;
                chh_u_min = A_j + (w_j_u + wh_j*m_j + shock_i) - n_j*inv - extmax_A/(1+r);
                chh_u_max = A_j + (w_j_u + wh_j*m_j + shock_i) - n_j*inv - extmin_A/(1+r);
                cu_vector = linspace(chh_u_min,chh_u_max,c_n);
                cu_vector(cu_vector < 0) = 0;
                
                % loop over consumption
                for k = 1:1:c_n
                    k;
                    
                    % HH consumption
                    chh_r = cr_vector(k);
                    chh_n = cn_vector(k);
                    chh_u = cu_vector(k);
                    % woman's consumption
                    cw_r = 0.5*chh_r;
                    cw_n = 0.5*chh_n;
                    cw_u = 0.5*chh_u;
                    
                    % Sector-Specific Utility:
                    u_r(k) = (cw_r^(1-sigma))/(1-sigma) + psi_r + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_n(k) = (cw_n^(1-sigma))/(1-sigma) + psi_n + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    u_u(k) = (cw_u^(1-sigma))/(1-sigma) + kappa*(1+theta1*m_j+theta2*n_j+theta3*K_j);
                    
                    % married:
                    if x <= 10
                    
                        % regular job:
                        A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh_r - n_j*inv); % eq. 8
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vm_r_next = sum(alpha(x_next,:).*Base,2); %cheby_approx
                        Amr_next(k)=A_next;
                        Vmr_next(k,x)=Vm_r_next;
                        % non-regular job:
                        A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh_n - n_j*inv);
                        x_next = x + 1;
                        if x_next == 11
                            x_next = 10;
                        end
                        % value function:
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vm_n_next = sum(alpha(x_next,:).*Base,2);
                        Amn_next(k)=A_next;
                        Vmn_next(k,x)=Vm_n_next;
                        % unemployed:
                        A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh_u - n_j*inv);
                        x_next = x;
                        % value function:
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vm_u_next = sum(alpha(x_next,:).*Base,2);
                        Amu_next(k)=A_next;
                        Vmu_next(k,x)=Vm_u_next;
                        % Sector-Specific Value Functions
                        Vm_r(k) = u_r(k) + beta * Vm_r_next;
                        Vm_n(k) = u_n(k) + beta * Vm_n_next;
                        Vm_u(k) = u_u(k) + beta * Vm_u_next;
                        % save marriage values (for marriage decision)
                        Vm_r_aux(k,x) = Vm_r(k);
                        Vm_n_aux(k,x) = Vm_n(k);
                        Vm_u_aux(k,x) = Vm_u(k);
                        
                    % Single:
                    else
                        
                        % Regular:
                        A_next = (1+r) * (A_j + (w_j_r + wh_j*m_j + shock_i) - chh_r - n_j*inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vs_r_next = sum(alpha(x_next,:).*Base,2);
                        Asr_next(k)=A_next;
                        Vsr_next(k,x)=Vs_r_next;
                        % Non-regular:
                        A_next = (1+r) * (A_j + (w_j_n + wh_j*m_j + shock_i) - chh_n - n_j*inv);
                        x_next = x + 1;
                        if x_next == 21
                            x_next = 20;
                        end
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vs_n_next = sum(alpha(x_next,:).*Base,2);
                        Asn_next(k)=A_next;
                        Vsn_next(k,x)=Vs_n_next;
                        % Unemployed:
                        A_next = (1+r) * (A_j + (w_j_u + wh_j*m_j + shock_i) - chh_u - n_j*inv);
                        x_next = x;
                        Base=kron(chebpoly_base(18+1, d_A*(A_next - extmin_A) - 1),kron(chebpoly_base(3+1, d_H*(wh_next - extmin_H) - 1),chebpoly_base(3+1, d_K*(K_next - extmin_K) - 1)));
                        Vs_u_next = sum(alpha(x_next,:).*Base,2);
                        Asu_next(k)=A_next;
                        Vsu_next(k,x)=Vs_u_next;
                        % Sector-Specific Value Functions
                        Vs_r(k) = u_r(k) + beta * Vs_r_next;
                        Vs_n(k) = u_n(k) + beta * Vs_n_next;
                        Vs_u(k) = u_u(k) + beta * Vs_u_next;
                        Vsm_r(k) = prob_marr_r*Vm_r_aux(k,x-10) + (1-prob_marr_r)*Vs_r(k);
                        Vsm_n(k) = prob_marr_n*Vm_n_aux(k,x-10) + (1-prob_marr_n)*Vs_n(k);
                        Vsm_u(k) = prob_marr_u*Vm_u_aux(k,x-10) + (1-prob_marr_u)*Vs_u(k);
                    end
                    end

                % optimal consumption and max Emax
                if x <= 10
                    % check
                    Vm_r(Amr_next < extmin_A) = NaN;
                    Vm_r(Amr_next > extmax_A) = NaN;
                    Vm_n(Amn_next < extmin_A) = NaN;
                    Vm_n(Amn_next > extmax_A) = NaN;
                    Vm_u(Amu_next < extmin_A) = NaN;
                    Vm_u(Amu_next > extmax_A) = NaN;
                    % save optimal
                    [Vm_r_star, Index_mr_k] = max(Vm_r);
                    [Vm_n_star, Index_mn_k] = max(Vm_n);
                    [Vm_u_star, Index_mu_k] = max(Vm_u);
                    cm_r_star = cr_vector(Index_mr_k);
                    cm_n_star = cn_vector(Index_mn_k);
                    cm_u_star = cu_vector(Index_mu_k);
                    cm_star_aux = [cm_r_star, cm_n_star, cm_u_star];
                    [Vm_star, Index_lm] = max([Vm_r_star,Vm_n_star,Vm_u_star]);
                else
                    % check
                    Vs_r(Asr_next < extmin_A) = NaN;
                    Vs_r(Asr_next > extmax_A) = NaN;
                    Vs_n(Asn_next < extmin_A) = NaN;
                    Vs_n(Asn_next > extmax_A) = NaN;
                    Vs_u(Asu_next < extmin_A) = NaN;
                    Vs_u(Asu_next > extmax_A) = NaN;
                    % save optimal
                    [Vs_r_star, Index_sr_k] = max(Vs_r);
                    [Vs_n_star, Index_sn_k] = max(Vs_n);
                    [Vs_u_star, Index_su_k] = max(Vs_u);
                    [Vsm_r_star, Index_smr_k] = max(Vsm_r);
                    [Vsm_n_star, Index_smn_k] = max(Vsm_n);
                    [Vsm_u_star, Index_smu_k] = max(Vsm_u);
                    cs_r_star = cr_vector(Index_sr_k);
                    cs_n_star = cn_vector(Index_sn_k);
                    cs_u_star = cu_vector(Index_su_k);
                    csm_r_star = cr_vector(Index_smr_k);
                    csm_n_star = cn_vector(Index_smn_k);
                    csm_u_star = cu_vector(Index_smu_k);
                    cs_star_aux = [csm_r_star, csm_n_star, csm_u_star, cs_r_star, cs_n_star, cs_u_star];
                    [Vs_star, Index_ls] = max([Vsm_r_star, Vsm_n_star, Vsm_u_star, Vs_r_star, Vs_n_star, Vs_u_star]);
                end
            % save choice:
            if x <= 10
                c_star(j, i, x, t) = cm_star_aux(Index_lm);
                l_star(j, i, x, t) = Index_lm;
                V_star(j, i, x, t) = Vm_star;
            else
                c_star(j, i, x, t) = cs_star_aux(Index_ls);
                l_star(j, i, x, t) = Index_ls;
                V_star(j, i, x, t) = Vs_star;
            end
            
            % save the number assets outside grid
            if x <= 10
                Ar_out(j,i,x,t) = sum(Amr_next < extmin_A) + sum(Amr_next > extmax_A);
                An_out(j,i,x,t) = sum(Amn_next < extmin_A) + sum(Amn_next > extmax_A);
                Au_out(j,i,x,t) = sum(Amu_next < extmin_A) + sum(Amu_next > extmax_A);
            else
                Ar_out(j,i,x,t) = sum(Asr_next < extmin_A) + sum(Asr_next > extmax_A);
                An_out(j,i,x,t) = sum(Asn_next < extmin_A) + sum(Asn_next > extmax_A);
                Au_out(j,i,x,t) = sum(Asu_next < extmin_A) + sum(Asu_next > extmax_A);
            end
            end
        end
        
        % Integrate over shocks
        W(:,x,t) = pi^(-1/2)*V_star(:,:,x,t)*weight;
         
        % reshape policy func
        c_func(:,x,t) = reshape(c_star(:,:,x,t),[],1);
        l_func(:,x,t) = reshape(l_star(:,:,x,t),[],1);
    end %this becomes end of loop x
end

% three labor functions (as 0 or 1)
lr_func = l_func == 1 | l_func == 4;
ln_func = l_func == 2 | l_func == 5;
lu_func = l_func == 3 | l_func == 6;

% marriage function for single women: equals 1 if labor choice is 1, 2, or 3
m_func = l_func == 1 | l_func == 2 | l_func == 3;

end
%% 5 funcs: consumption, regular, non-regular, unemployed, and marriage
