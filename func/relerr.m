function [ out ] = relerr(Y1, Y2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    out = norm_tensor(Y1-Y2)/norm_tensor(Y2);

end

