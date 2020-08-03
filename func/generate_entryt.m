function [A_entryt, aA_entryt] = generate_entryt(n,P_omega)
%GENERATE_ENTRYT Entry-wise subsample transform
    if max(P_omega) > n
        fprintf('Error! Dimension n not large enough!\n')
    else
        S = speye(n,n);
        S = S(P_omega,:);

        A_entryt = @(x) S*x;
        aA_entryt = @(y) S'*y;
    end
end

