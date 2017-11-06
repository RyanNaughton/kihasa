function EVapp=cheby_approx(alpha,ncheby,a1,b1,c1,da,db,dc,a,b,c)

Base=kron(chebpoly_base(ncheby+1, da*(a - a1) - 1),kron(chebpoly_base(ncheby+1, db*(b - b1) - 1),chebpoly_base(ncheby+1, dc*(c - c1) - 1)));
EVapp=sum(alpha.*Base,2); 


 