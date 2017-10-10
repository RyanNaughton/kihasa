function [W,A,Minv,Tinv,C]=thci_sol_grid20(params,P,G,B,S)

gamma=params([1:6]);
phi=params(7);
rho=params(8);
delta=params(9);

V=zeros(G.nss,G.Ne,G.nper); % Value function 
W=zeros(G.nss,G.ntime); % EValue function

W(:,G.ntime)=P.eta*S.h_(:,G.ntime).^(1-P.sigma)/(1-P.sigma)+ P.tau1*(1-exp(-S.a_(:,G.ntime)));   


for t=G.ntime-1:-1:1
    
    %% A.Polynomial Approximation of Expected Value Function
         % 1. Input is W(t+1)
         % 2. Output is vector of polynomial coefficients alp2(t)
     
       coef1=W(:,t+1)'*kron(B.B,B.B);
       aux=diag(B.B'*B.B)';
       aux2=kron(aux,aux);
       alp2(t,(1:(G.ncheby+1)^2))=coef1./aux2;
         
    %% B. Optimization using polynomial tensor product  
    for i=1:1:G.Ne       
        for j=1:1:G.nss  

            for k=1:G.Nc
                for l=1:G.Nc
                    for m=1:G.Ntinv
                            [VV(k,l,m),apr(k,l,m), hpr(k,l,m)]=obj_fun(alp2(t,:),params,t,P,G,S,S.h_(j,t),S.a_(j,t),S.inc_min(i,t),S.inc_wage(i,t),S.c(j,k),S.minv(j,l),S.tinv(m));
                            bc_val(k,l,m)=(apr(k,l,m)>=S.a1(t+1) && apr(k,l,m)<=S.a2(t+1) && hpr(k,l,m)>=S.h1(t+1) && hpr(k,l,m)<=S.h2(t+1));                                                   
                    end
                end 
             end
           
             bc_val = double(bc_val);
             bc_val(bc_val == 0)=NaN;
             auxV=zeros(G.Ntinv,1);

             Ind_bc_check=VV.*bc_val;
             for m=1:1%G.Ntinv
                Ind_bc=VV(:,:,m).*bc_val(:,:,m);
                [V_tinv,Indaux] = max(Ind_bc(:));
                auxV(m)=V_tinv;
                [I_c(m), I_minv(m)] = ind2sub(size(Ind_bc),Indaux);
             end
             [V(j,i,t), Ind]=max(auxV(:));                  
             Minv(j,i,t)=S.minv(j,I_minv(Ind));
             C(j,i,t)=S.c(j,I_c(Ind));
             Tinv(j,i,t)=S.tinv(Ind);
                    
            A(j,i,t)=(1+P.r)*S.a_(j,t) + S.inc(i,t)*(1-Tinv(j,i,t)) - P.p*Minv(j,i,t) - C(j,i,t);
            %H(j,i,t)=delta*CES(S.h_(j,t),opty(1),gamma(t),rho,phi);
           j    
        end
       i
    end
    
    W(:,t)=pi^(-1/2)*V(:,:,t)*G.wt;
    % check integration of cons and inv

 
t   

end


