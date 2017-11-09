function EVapp=cheby_approx_5D(alpha,ncheby,a1,b1,c1,d1,e1,da,db,dc,dd,de,a,b,c,d,e)

Base=kron(chebpoly_base(ncheby+1, da*(a - a1) - 1),kron(chebpoly_base(ncheby+1, db*(b - b1) - 1),kron(chebpoly_base(ncheby+1, dc*(c - c1) - 1),kron(chebpoly_base(ncheby+1, dd*(d - d1) - 1),chebpoly_base(ncheby+1, de*(e - e1) - 1)))));
EVapp=sum(alpha.*Base,2); 
