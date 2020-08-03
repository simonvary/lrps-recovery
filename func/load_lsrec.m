function [opts] = load_lsrec(lsrec_name, matsize)
%LOAD_LSREC Loads default param settings for ls recovery algorithms

    opts.m = matsize(1);
    opts.n = matsize(2);
    
    opts.MAX_ITER = 500;
        
    opts.tol_res = 1e-6;
    opts.tol_chg = 1e-5;
    opts.tol_rel = 1-(1e-3);        % 10 = turned off
    opts.shift_rel = 15;
    
    opts.init = [];                 % initial guess
    opts.verb = 0;
    opts.time_iter = 0;
    
    opts.L_true = [];
    opts.S_true = [];

    if strcmp(lsrec_name, 'lsrec_iht') 
        opts.alpha = 1;             % step-size
    elseif strcmp(lsrec_name, 'lsrec_niht')
        opts.alpha = [];
    elseif strcmp(lsrec_name, 'lsrec_aht')
        opts.alpha_L = 1;
        opts.alpha_S = 1;
    elseif strcmp(lsrec_name, 'lsrec_naht')
        opts.alpha_L = [];
        opts.alpha_S = [];
    elseif  strcmp(lsrec_name, 'lsrec_cvx')
        
    end
end

