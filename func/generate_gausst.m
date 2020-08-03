function [A, aA] = generate_gausst(n,p)
%GENERATE_GAUSST Gaussian transform  
    Phi = randn(p, n) / sqrt(n);
    A = @(x) Phi*x;
    aA = @(y) Phi'*y;
end