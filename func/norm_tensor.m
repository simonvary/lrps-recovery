function [ out ] = norm_tensor(Y)
%NORM_TENSOR Frobenius norm of a tenson
    out = sqrt(sum(sum(sum(Y.^2))));
end

