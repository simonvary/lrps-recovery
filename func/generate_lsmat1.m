function [L, S, X] = generate_lsmat1(m, n, r, s, c)
%GENERATE_PROBLEM Generate low-rank plus sparse problem in line with 
%               Netrapalli et al 14, Yi et al 14, etc
%           The only exception is for r = 0, we take mean(abs(L(:))) to be 0.2
%           this modification allows for solution of CS only problems

    P = randn(m,r);
    Q = randn(n,r);
    L = P*Q';
    
    S = zeros(m, n);
    S_ind = randperm(m*n,s);
    S(S_ind) = (2*rand(s,1)-1)*c* max(mean(abs(L(:))), 0.2);
    X = L + S;
end

