function [ M, ind ] = matproj_sparse( M, s, xi )
%MATPROJ_SPARSE Euclidean projection of space of sparsity k matrices
%   Detailed explanation goes here

    if ~isempty(s)
        [~, ind2] = sort(abs(M(:)));
        ind = ind2(numel(M)-s+1:end);
        M(ind2(1:(numel(M)-s))) = 0;
    elseif ~isempty(xi)
        ind = (M(:) >= xi);
        M(~ind) = 0;
    end
end

