function [ M, U, Sin, V] = matproj_rank( M, r )
%MATPROJ_RANK Euclidean projection of space of rank r matrices

    try
        [U, Sin, V] = svd(M,'econ');
        U = U(:,1:r);
        Sin = Sin(1:r,1:r);
        V = V(:,1:r);
        %[U, Sin, V] = svds(M,r);
        M = U*Sin*V';
    catch ME
        fprintf('Failed to converge svd. Projection by pivoted QR.\n')
        [U,R1,P1] = qr(M);
        [V,R2,P2] = qr(P1*R1');
        Sin = P2*R2';
        
        U = U(:,1:r);
        V = V(:,1:r);
        Sin = Sin(1:r,1:r);
        try
            [U, Sin, V] = svds(U*Sin*V', r);
        catch ME2
            fprintf('Failed to make svd on pivoted QR, keeping only QR.\n')
        end
       % U = U(:,1:r);
        %Sin = Sin(1:r,1:r);        %V = V(:,1:r);
        M = U*Sin*V';
    end
end