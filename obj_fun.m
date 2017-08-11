function [objf,apr,hpr] = obj_fun(alpha,params,t,P,G,S,h,a,inc_min,inc_wage,c,minv,tinv)

%% Borrar %%
%alpha=alp2(t,:)
%a=S.a_(j,t)
%h=S.h_(j,t)
%inc=S.inc(i,t)
%cc=S.c(j,k)
%minv=S.minv(j,l)
%tinv=S.tinv(m)
%inc_min=S.inc_min(i,t)
%inc_wage=S.inc_wage(i,t)

%% 
u=utilfun(c,P);
gamma=reshape(params([1:6]),3,2);
phi=params(7);
rho=params(8);
delta=params(9);
inc=inc_min+inc_wage;
apr=(1+P.r)*a + inc*(1-tinv) - P.p*minv - c;

        hpr=delta*CES(h,minv,tinv,gamma(t,:),rho,phi,0);
        EV=cheby_approx(alpha,G.ncheby,S.extmin1(t+1),S.extmin2(1,t+1),S.dh(t+1),S.da(1,t+1),hpr,apr);
   
    
objf= u + P.beta*EV;

