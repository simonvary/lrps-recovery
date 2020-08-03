function [L, S, X, out] = lsrec_cvx(b, A, r, s, opts)
%LSREC_CVX Summary of this function goes here
%   Detailed explanation goes here

    % low rank approximation using nuclear norm
    lambda = sqrt(r/s);
    tic;
    warning('off')
    cvx_begin quiet
        cvx_precision default
        variable L(opts.m,opts.n)
        variable S(opts.m,opts.n)
        minimize norm_nuc(L) + lambda*norm(S(:),1)
        subject to
        A(L(:)+S(:))==b
    cvx_end
    warning('on')
    [~ , ind] = matproj_sparse(S, s, []);
  %  [L , ~, ~, ~] = matproj_rank(L, r);
    X = L + S;

    out.time = toc;
    out.lambda = lambda;
    out.err_res = relerr(A(X(:)), b);
    is_true_given = ~isempty(opts.L_true) && ~isempty(opts.S_true);
    if is_true_given     % true error tracking
        opts.X_true = opts.L_true + opts.S_true;
        out.err_L = relerr(L, opts.L_true);
        out.err_S = relerr(S, opts.S_true);
        out.err_X = relerr(X, opts.X_true);
        ind_true = find(opts.S_true ~= 0);          % for support tracking
        out.err_S_ind = length(intersect(ind_true,ind));
        if opts.verb == 1
            fprintf('Finished (r=%d, s=%d), err_res=%1.2e, err_X=%1.2e\n', r, s, out.err_res, out.err_X)
        end
    else
        if opts.verb == 1
            fprintf('Finished (r=%d, s=%d), err_res=%1.2e\n', r, s, out.err_res)
        end
    end

end
