function [A_fjlt, aA_fjlt] = generate_fjlt_mat(n,p)
%GENERATE_FJLT Fast Johnson-Lindenstrauss transform 
    S_fjtl = speye(n,n);
    S_fjtl = S_fjtl(randperm(n, p),:);
    D_fjlt = 2*((rand(n,1) > 0.5) - 0.5);
    D_fjlt = spdiags(D_fjlt(:),0,n,n);
    
    fjlt_mat = S_fjtl*(dct(D_fjlt*eye(n,n)));
    A_fjlt = @(x) fjlt_mat*x;
    aA_fjlt = @(y)fjlt_mat'*y;
    
end

