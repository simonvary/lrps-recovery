function [L, S, X, out] = lsrec_naht(b, A, aA, r, s, opts)
%LSREC_NAHT Normalized Alternating Hard Thresholding for low-rank plus 
%           sparse recovery

%   Input parameters:
%   b               - vector of observed entries
%   opts.mask       - indices of observed entries
%   A               - forward operator - R^mn to R^p, linear, tight frame
%   aA              - adjoint operator - R^p to R^mn
%   r, s            - rank and sparsity
%   opts.m, opts.n  - matrix size
%   opts.alpha_L    - stepsize in direction of L (set = [] for normalized)
%   opts.alpha_S    - stepsize in direction of S (set = [] for normalized)
%   opts.verb       - verbosity
%   opts.time_iter  - timings of individual iterations
%                     (accurate when each iteration takes longer than ~ 1/10 sec. )
 

    err_res = zeros(opts.MAX_ITER, 1);
    err_rel = zeros(opts.MAX_ITER, 1);
    
    is_true_given = ~isempty(opts.L_true) && ~isempty(opts.S_true);
    if is_true_given     % true error tracking
        opts.X_true = opts.L_true + opts.S_true;
        err_L = zeros(opts.MAX_ITER, 1);
        err_S = zeros(opts.MAX_ITER, 1);
        err_X = zeros(opts.MAX_ITER, 1);
        ind_true = find(opts.S_true ~= 0);          % for support tracking
    end

    if opts.time_iter
        iter_timings = zeros(opts.MAX_ITER, 1);
    end
    
    tic;
    if isempty(opts.init) % make an initialization
        % L, S, A(l), A(s) are zero
        grad_L = reshape(aA(-b), [opts.m, opts.n]);
        L = -grad_L;
        [L, U, ~, ~] = matproj_rank(L, r);        
        Al = A(L(:));
        
        grad_S = reshape(aA(Al - b), [opts.m, opts.n]);
        S = -grad_S;
        [S, T_ind] = matproj_sparse(S, s, []);
        As = A(S(:));
        
        X = L + S;
    else
        [L, S, X, U, ~, T_ind] = matproj_ls_accaltproj(opts.init.X, r, s, [], [], []);
        Al = A(L(:));
        As = A(S(:));
    end
    
    % error tracking for initial guess
    err_res(1) = relerr(A(X(:)), b);
    if is_true_given 
        err_L(1) = relerr(L, opts.L_true);
        err_S(1) = relerr(S, opts.S_true);
        err_X(1) = relerr(X, opts.X_true);
    end
    
    if opts.time_iter
        iter_timings(1) = toc;
    end
    
    done = 0;
    l = 2;
    while (~done)
        % Gradient step for L
        grad_L = reshape(aA(Al + As - b), [opts.m, opts.n]);
        if isempty(opts.alpha_L)
            %proj_grad_L = U*(U'*L);
            proj_grad_L = matproj_lssub(grad_L, U, T_ind, 4);
            c1 = norm(proj_grad_L,'fro')^2;
            c2 = norm(A(proj_grad_L(:)))^2;
            L = L - (c1/c2) * grad_L;
        else
            L = L - opts.alpha_L * grad_L;
        end       
        % Thresholding for L
        [L, U, ~, ~] = matproj_rank(L, r);
        Al = A(L(:));

        % Gradient step for S
        grad_S = reshape(aA(Al + As - b), [opts.m, opts.n]);
        if isempty(opts.alpha_S)
            %proj_grad_S = proj_grad_S; %zeros(opts.m, opts.n);
            proj_grad_S = matproj_lssub(grad_S, U, T_ind, 4);
            proj_grad_S(T_ind) = grad_S(T_ind);
            c1 = norm(proj_grad_S,'fro')^2;
            c2 = norm(A(proj_grad_S(:)))^2;
            S = S - (c1/c2) * grad_S;
        else
            S = S - opts.alpha_S * grad_S;
        end
        % Thresholding for S
        [S, T_ind] = matproj_sparse(S, s, []);      
        As = A(S(:));

        X = L + S;
        if opts.verb
            fprintf('Iter%d: %1.5f\n', l, relerr(X, opts.X_true))
        end
    
        err_res(l) = relerr(Al + As, b);
        
        if l > (opts.shift_rel+1)   % err_res(1) == 0 is the initial guess !
            err_rel(l) = (err_res(l)/err_res(l-opts.shift_rel))^(1/opts.shift_rel);
        else
            err_rel(l) = 0;
        end

        if is_true_given 
            err_L(l) = relerr(L, opts.L_true);
            err_S(l) = relerr(S, opts.S_true);
            err_X(l) = relerr(X, opts.X_true);
        end
        
        if opts.time_iter
            iter_timings(l) = toc;
        end
        
        l = l + 1;
        done = (l > opts.MAX_ITER || err_res(l-1) <= opts.tol_res || err_rel(l-1) >= opts.tol_rel);
    end
    
    out.iter = l-1;
    if opts.time_iter
        out.iter_timings = iter_timings(1:out.iter);
        out.time = out.iter_timings(end);
    else
        out.time = toc;
    end
    out.err_res = err_res(1:out.iter);
    out.err_rel = err_rel(1:out.iter);
    if is_true_given 
        out.err_L   = err_L(1:out.iter);
        out.err_S   = err_S(1:out.iter);    
        out.err_X   = err_X(1:out.iter);
        ind = find(S ~= 0);       % support checking
        out.err_S_ind = length(intersect(ind_true,ind));
       %fprintf('Finished (r=%d, s=%d) in %d iterations, err_res=%1.2e, err_X=%1.2e\n', r, s, l, out.err_res(end), out.err_X(end))
    else       
       % fprintf('Finished (r=%d, s=%d) in %d iterations, err_res=%1.2e\n', r, s, l, out.err_res(end))
    end
end


function [X] = matproj_lssub(M, U, T, k)
    % Project matrix M on P_lr, P_sp
    for i = 1:k
        X = U*(U'*M);
        X(T) = M(T);
    end
end
