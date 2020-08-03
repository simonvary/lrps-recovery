function [success_rate, errs_X, errs_L, errs_S, errs_res] = phase_single_general(mat_size, delta, rho_r, rho_s, num_trials, generate_matrix, generate_sensing, lsrec, tol)
%PHASE_SINGLE_GENERAL General function to perform a single experiment with the setup 
%   [m n] = mat_size;
%   delta = p/(m*n)
%   rho_r = r(m+n-r)/p
%   rho_s = s/p
%   num_trials - number of trials
%   generate_matrix - function handle to LS(r,s) generator
%   generate_sensing - function handle to A, aA generator
%   lsrec - recovery algorithm, takes (b, A, aA, r, s),
%           all additional info should be setup before

    if length(mat_size) == 1
        m = mat_size;
        n = mat_size;
    else
        m = mat_size(1);
        n = mat_size(2);
    end
    
    p = ceil(delta*m*n);
    r = ceil(((m+n)-sqrt((m+n)^2 - 4*rho_r*p))/2);
    s = ceil(rho_s*p);
    
    errs_X = zeros(num_trials, 1);
    errs_L = zeros(num_trials, 1);
    errs_S = zeros(num_trials, 1);
    errs_res = zeros(num_trials, 1);
    
    for i = 1:num_trials      
        % Setup the problem
        [L_true, S_true, X_true] = generate_matrix(r,s);
        [A, aA] = generate_sensing(m*n, p);
        b = A(X_true(:));
        
        % Run recovery
        [L, S, X, ~] = lsrec(b, A, aA, r, s);
        
        errs_X(i) = relerr(X, X_true);
        errs_L(i) = relerr(L, L_true);
        errs_S(i) = relerr(S, S_true);
        errs_res(i) = relerr(A(X(:)), b);
    end
    success_rate = sum(errs_X <= tol)/num_trials;
end

