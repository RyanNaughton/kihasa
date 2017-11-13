function [z,ext,extmin,extmax,d,vector,T,T2] = cheby_values(n,ub,lb)
    
    % points
    M = n;
    
    % degrees
    ncheby = M-2;
    
    % Chebyshev
    z = flipud(cos(((2*(1:M)'-1)*pi)/(2*M)));

    % parameter for expanded Chebyshev polynomial approximation 
    ext = -(z(1)+1)*(ub-lb)/2/z(1);
    
    % bounds
    extmin = lb - ext;
    extmax = ub + ext;

    % gradient
    d = 2./(extmax-extmin);

    % vector
    vector(1)=lb;
    vector(M)=ub;
    for j=2:1:M-1
       vector(j)= extmin + (1+z(j))./d;
    end
    
    % Base 
    T = chebpoly_base(ncheby+1,z);
    % Base Squared
    T2 = diag(T'*T);
    
end
