function [L, S, X, out] = lsrec_niht(b, A, aA, r, s, matproj_ls, opts)
%LSREC_NIHT (Normalized) Iterative Hard Thresholding for low-rank plus sparse recovery

%   Input parameters:
%   b               - vector of observed entries
%   opts.mask       - indices of observed entries
%   A               - forward operator - R^mn to R^p, linear, tight frame
%   aA              - adjoint operator - R^p to R^mn
%   r, s            - rank and sparsity
%   matproj_ls      - thresholding operation
%               function handle
%                   [L, S, X] = matproj_LS(X, r, s)
%   opts.m, opts.n  - matrix size
%   opts.alpha      - stepsize


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
        X = reshape(aA(b), [opts.m, opts.n]);
       % fprintf('Iter%d (P1): %1.5f\n', 1, relerr(X, opts.X_true))
        [L, S, X, U, ~, T] = matproj_ls_accaltproj(X, r, s, [], [], []);
        %fprintf('Iter%d (P2): %1.5f\n', 1, relerr(X, opts.X_true))
    else
        X = opts.init.X;
        [L, S, X, U, ~, T] = matproj_ls(X, r, s, [], [], []);
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
        % Gradient step
        Ax = A(X(:));
        grad = reshape(aA(Ax - b), [opts.m, opts.n]);

        if isempty(opts.alpha)
            proj_grad_ls = matproj_lssub(grad, U, T, 4);
            c1 = norm(proj_grad_ls,'fro')^2;
            c2 = norm(A(proj_grad_ls(:)))^2;
            X = X - (c1/c2) * grad;
        else
            X = X - opts.alpha * grad;
        end
       % fprintf('Iter%d (P1): %1.5f\n', l, relerr(X, opts.X_true))

        % Thresholding step
        [L, S, X, U, ~, T] = matproj_ls(X, r, s, min(0.1*err_res(l-1), 1e-5), L, S);

        err_res(l) = relerr(Ax, b);

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
        
        if opts.verb
            fprintf('Iter%d: %1.5f\n', l, relerr(Ax, b))
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
        %fprintf('Finished (r=%d, s=%d, p=%d) in %d iterations, err_res=%1.2e, err_X=%1.2e\n', r, s, p, l, out.err_res(end), out.err_X(end))
    else
        %fprintf('Finished (r=%d, s=%d, p=%d) in %d iterations, err_res=%1.2e\n', r, s, p, l, out.err_res(end))
    end
end


function [X] = matproj_lssub(M, U, T, k)
    % Project matrix M on P_lr, P_sp
    for i = 1:k
        X = U*(U'*M);
        X(T) = M(T);
    end
end
