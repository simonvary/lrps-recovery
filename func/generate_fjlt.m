function [A_fjlt, aA_fjlt] = generate_fjlt(n,p)
%GENERATE_FJLT Fast Johnson-Lindenstrauss transform 
    S_fjtl = speye(n,n);
    S_fjtl = S_fjtl(randperm(n, p),:);
    D_fjlt = 2*((rand(n,1) > 0.5) - 0.5);
    D_fjlt = spdiags(D_fjlt(:),0,n,n);

    A_fjlt = @(x) S_fjtl*(dct(D_fjlt*x));
    aA_fjlt = @(y) D_fjlt'*(idct(S_fjtl'*y));
    
end

