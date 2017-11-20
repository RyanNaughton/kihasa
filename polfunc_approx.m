
function [alpC,alpR,alpT,alpU,alpM]=thci_polfunc(C,R,T,U,M,S,G)


T_sim=kron(S.T_A, kron(S.T_H,kron(S.T_K,kron(S.Teps_r,kron(S.Teps_n,S.Teps_u)))));
Den= kron(S.T2_A, kron(S.T2_H,kron(S.T2_K,kron(S.T2eps_r,kron(S.T2eps_n,S.T2eps_u)))));

%%% Polynomial Bases and Derivatives %%%% 

 % square of T multiplied
for t=1:1:G.nper
	for x = 1:1:(n_matstat*n_wrkexp)
	NumC(:,x,t) = C(:,x,t)'*T_sim;  
	alpC(:,x,t) = NumC(:,x,t)./Den';
	NumR(:,x,t) = R(:,x,t)'*T_sim; 
	alpR(:,x,t) = NumR(:,x,t)./Den';
	NumT(:,x,t) = T(:,x,t)'*T_sim; 
	alpT(:,x,t) = NumT(:,x,t)./Den';
	NumU(:,x,t) = U(:,x,t)'*T_sim; 
	alpU(:,x,t) = NumU(:,x,t)./Den';
	NumM(:,x,t) = M(:,x,t)'*T_sim; 
	alpM(:,x,t) = NumM(:,x,t)./Den';
	end
end


