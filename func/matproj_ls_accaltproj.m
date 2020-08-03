function [L, S, X, U, V, T] = matproj_ls_accaltproj(M, r, s, tol, L0, S0)
    % Robust PCA
    clear para
    para.verb = 0;
    if ~isempty(tol)
        para.tol = tol;
    end
    
    if (r == 0) && (s == 0)
        U = zeros(size(M,1),1);
        V = zeros(size(M,2),1);
        L = zeros(size(M));
        S = zeros(size(M));
        T = [];
    elseif (r == 0) && (s > 0)
        % Low-rank is zero
        U = zeros(size(M,1),1);
        V = zeros(size(M,2),1);
        L = zeros(size(M));
        % Sparse is hard thresholded
        [ S, ~ ] = matproj_sparse( M, s, [] );
        T = (S~=0);
    elseif (r > 0) && (s == 0)
        S = zeros(size(M));
        T = [];
        [L, U, ~, V] = matproj_rank( M, r );
    else
        [L, S, U, V, ~] = AccAltProj( M, r, para, L0, S0);
        % hard threshold the sparse component
        [ S, ~ ] = matproj_sparse( S, s, [] );
        T = (S~=0);
    end
    X = L + S;
end
