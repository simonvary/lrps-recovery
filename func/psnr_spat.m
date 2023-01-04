function [out] = psnr_spat(Y_approx, Y_true)
%PSNR_SPAT Compute spatial PSNR of a tensor approximation
% Note: Very slow, implement in faster way
    out = zeros(size(Y_true, 1), size(Y_true, 2));
    for i = 1:size(Y_true, 1)
        for j = 1:size(Y_true, 2)
            out(i, j) = 10*log10( max(max( Y_approx(i,j,:) ))^2/ immse(Y_approx(i,j,:),Y_true(i,j,:)));
        end
    end
end