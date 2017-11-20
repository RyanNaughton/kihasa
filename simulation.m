function [c_s,i_s,a_s,h_s] = thci_sim(params,alp_A,alp_I,G,P,S,epsy_sim,inc_sim,k,epsh_sim)

% Initial Conditions

% States 
a_s = zeros(G.npop,G.nper); % save assets according to law of motion
h_s= ones(G.npop,G.nper); % save HC according to law of motion

%a_s(:,1)=a1(1)+(a2(1)-a1(1))*rand(npop,1);



for t=1:1:G.nper
    for n=1:1:G.npop
        
    T_ew(t,:)=chebpoly_base(G.Ne-1, 2*(epsy_sim(n,t) - G.eps_y(1,1))/(G.eps_y(G.Ne,1)-G.eps_y(1,1)) - 1);
    T_h(t,:)=chebpoly_base(G.ncheb_pol, S.dh(t)*(h_s(n,t) - S.h1(t)) - 1);
    T_a(t,:)=chebpoly_base(G.ncheb_pol, S.da(t)*(a_s(n,t) - S.a1(1,t,k)) - 1);
    T_s=kron(T_h(t,:),kron(T_a(t,:),T_ew(t,:)));
    a_s(n,t+1)=sum(alp_A(t,:).*T_s,2);
    i_s(n,t)=sum(alp_I(t,:).*T_s,2);
    c_s(n,t)= inc_sim(n,t) + (1+P.r)*a_s(n,t) - P.p*i_s(n,t) - a_s(n,t+1);
    if i_s(n,t)<0
    h_s(n,t+1)=delta*CES(h_s(n,t),0,gamma(t),rho,phi,epsh_sim(n,t));
    else
    h_s(n,t+1)=delta*CES(h_s(n,t),i_s(n,t),gamma(t),rho,phi,epsh_sim(n,t));
    end
    n
    
    end
    t
end



end

