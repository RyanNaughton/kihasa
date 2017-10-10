function EVapp=cheby_approx(alpha,ncheby,h1,a1,dh,da,h,a)

Base=kron(chebpoly_base(ncheby+1, dh*(h - h1) - 1),chebpoly_base(ncheby+1, da*(a - a1) - 1));
EVapp=sum(alpha.*Base,2); 


 