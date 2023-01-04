function [out] = psnr_spec(Y_approx, Y_true)
%PSNR_spec Compute spectral PSNR of a tensor approximation
    out = zeros(1, size(Y_true, 3));
    max_val = max(max(max(Y_true)));
    for b = 1:size(Y_true, 3)
      %  out(b) = 10*log10( max_val^2/ immse(Y_approx(:,:,b),Y_true(:,:,b))); the same thing
        out(b) = psnr(Y_approx(:,:,b), Y_true(:,:,b), max_val);
    end
end