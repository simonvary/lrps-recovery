function [L, S, X, out] = lsrec_sparcs(b, A, aA, r, s, opts)
%LSREC_SPARCS SPArse and low Rank decomposition via Compressive Sensing 
%           Implementation of SpaRCS by Waters, Sankaranarayanan, Baraniuk
%           adapted from code: https://github.com/image-science-lab/SpaRCS
%           in the paper: SpaRCS: Recovering Low-Rank and Sparse Matrices 
%                         from Compressive Measurements, NIPS 2011
%
%           Should be rewritten in a more concise form, messy code, 
%           To be finished. Use the original sparcs() function

%   Input parameters:
%   b               - vector of observed entries
%   opts.mask       - indices of observed entries
%   A               - forward operator - R^mn to R^p, linear, tight frame
%   aA              - adjoint operator - R^p to R^mn
%   r, s            - rank and sparsity
%   opts.m, opts.n  - matrix size
%   opts.alpha_L    - stepsize in direction of L (set = [] for normalized)
%   opts.alpha_S    - stepsize in direction of S (set = [] for normalized)
%
%

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
    
    done = 0;
    l = 2;
    while (~done) ((err>tol)&&(iterCnt<=maxIter))
        iterCnt = iterCnt + 1;
    
        %%Do greedy matrix steps
        %Compute new basis for residual, store in Psiprime
        rt = reshape(aA(b-A(Xhat(:))-A(s_cosamp(:))), [opts.m, opts.n]);
        % switch for different SVDs
        switch lower(svdMethod)
            case 'svdlibc'
                [U diagS V] = svdlibc(rt, selectAtom*r);
            case 'propack'
                [U,S,V] = lansvd(rt,selectAtom*r,'L');
                diagS = diag(S);
            case 'svds'
                [U, S, V] = svds(rt, selectAtom*r);
                diagS = diag(S);
            case 'spsvd'
                [U, S, V] = spsvd(rt, selectAtom*r, 1e-2);
                diagS = diag(S);
            otherwise
                [U,S,V] = svd(rt,0);
                diagS = diag(S);
        end
    
        Uprime = U(:, 1:min(selectAtom*r,size(U,2)));
        Vprime = V(:, 1:min(selectAtom*r,size(V,2)));
    
    
        %Merge supports
        PsitildeU = [PsihatU Uprime];
        PsitildeV = [PsihatV Vprime];
    
        %Do least squares estimation to get Xhat
        AP = @(z) APsiUV(z,A,PsitildeU, PsitildeV);
        APt = @(z) APsitUV(z,At,PsitildeU, PsitildeV);
        ALS = @(z) APt(AP(z));
        
        [alpha, res, iter] = cgsolve(ALS, APt(b-A(s_cosamp(:))), 1e-6, 100, 0);
        Xtilde = PsitildeU*diag(alpha)*PsitildeV'; 
    
        %Update Psihat
        switch lower(svdMethod)
            case 'svdlibc'
                [U diagS V] = svdlibc(Xtilde, r);
            case 'propack'
                [U,S,V] = lansvd(Xtilde,r,'L');
                diagS = diag(S);
            case 'svds'
                [U, S, V] = svds(Xtilde, r, 'L');
                diagS = diag(S);
            otherwise
                [U,S,V] = svd(Xtilde,0);
                diagS = diag(S);
        end
        PsihatU = U; PsihatV = V;
        Xhat = U*S*V';
    
        %%Do sparse matrix steps
        rr = b - A(s_cosamp(:)) - A(Xhat(:));
        proxy = aA(rr);% does not need to reshape proxy = proxy(:);
        
        %---Estimate support
        [tmp,ww]= sort(abs(proxy),'descend');
        tt= union(find(ne(s_cosamp,0)),ww(1:(selectAtom*K)));
    
        % Preparation for cg_solve
        PP_tt = @(z) A_I(A,z,tt,N);
        PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
        PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
        
        qq = PP_transpose_tt(b-A(Xhat(:)));   
        
        %Pseudo-inverse
        [w, res, iter2] = cgsolve(PPtranspose_PP_tt, qq,1e-6, 100, 0);
        bb= 0*s_cosamp; bb(tt)= w;
        
        %---Prune
        kk = iterCnt;
    	[tmp,ww2]= sort(abs(bb),'descend');
        s_cosamp=0*bb;
        s_cosamp(ww2(1:K)) = bb(ww2(1:K));
        aa = s_cosamp;
    
        errTemp = norm(A(Xhat)+A(s_cosamp)-b)/norm(b);
        if (VERBOSE)
            fprintf('Iter: %d Err: %f\n',iterCnt,errTemp);
        end
    
        if ((err-errTemp)/err <= tol_early) %Early termination condition
            err = errTemp;
            if (VERBOSE)
                disp(sprintf('Terminating.'));
            end
            break 
        end


        err = errTemp;
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
        
        l = l + 1;
        done = (l > opts.MAX_ITER || err_res(l-1) <= opts.tol_res || err_rel(l-1) >= opts.tol_rel);
    end
    out.time = toc;
    out.iter = l-1;
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
