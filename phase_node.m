function [finished] = phase_node(node_id, num_cores, num_trials, max_success, tol_success, deltas, rho_rs, rho_ss, htc_single)
%PHASE_NODE Run expriment for a range of deltas and pairs of (rho_rs, rho_ss)
%
%   INPUT:
%       node_id (int)                   -- ID for the computing node
%       num_cores (int)                 -- Number of cores to utilize
%       num_trials (int)                -- Number of independent trials
%       max_success (int)               -- Successful recoveries until stopping
%       tol_success (float)             -- Probability tolerance for successful recovery
%       deltas, rho_rs, rho_ss (array)  -- Lists of parameters
%       htc_single (fun_handle)         -- Function of (delta, rho_r, rho_s)
%
%   For every rho_rs(i), rho_ss(i), compute through a list of deltas(j) 
%   in an increasing order until recovery occurrs at least max_success times. 
%   Jobs are spread over num_cores via parallel loop on a (num_node)^th 
%   computing node.
%
%   Each experiment consists of running htc_single(delta, rho_r, rho_s)
%   function that is a specific setup of the htc_single_general function.
    
    errs_X = zeros(length(rho_rs),length(deltas), num_trials);
    errs_L = zeros(length(rho_rs),length(deltas), num_trials);
    errs_S = zeros(length(rho_rs),length(deltas), num_trials);
    errs_res = zeros(length(rho_rs),length(deltas), num_trials);
    success_rates = zeros(length(rho_rs),length(deltas));
    
    delete(gcp('nocreate'))
    poolobj = parpool('local', num_cores)
    parfor (i = 1:numel(rho_rs), num_cores)
        rng(node_id*num_cores + i)
        
        rho_r = rho_rs(i);
        rho_s = rho_ss(i);
        
        tmp_errs_X = NaN(length(deltas), num_trials);
        tmp_errs_L = NaN(length(deltas), num_trials);
        tmp_errs_S = NaN(length(deltas), num_trials);
        tmp_errs_res = NaN(length(deltas), num_trials);
        success_rate_tmp = NaN(length(deltas),1);
        nsuccess = 0;
        j = 1;
        while (nsuccess < max_success) && (j <= length(deltas))
            delta = deltas(j);
            
            %fprintf('rho_r = %1.2e, rho_s = %1.2e', rho_r, rho_s);
            [success_rate_tmp_scalar, tmp_errs_X_vec, tmp_errs_L_vec, ...
                tmp_errs_S_vec, tmp_errs_res_vec] = htc_single(delta, rho_r, rho_s);
            
            success_rate_tmp(j) = success_rate_tmp_scalar;
            tmp_errs_X(j, :) = tmp_errs_X_vec;
            tmp_errs_L(j, :) = tmp_errs_L_vec;
            tmp_errs_S(j, :) = tmp_errs_S_vec;
            tmp_errs_res(j,:) = tmp_errs_res_vec;
        
            fprintf('Finished (delta = %1.2f, rho_r=%1.2f, rho_s=%1.2f) with %1.2f ratio.\n', ...
                            delta, rho_r, rho_s, success_rate_tmp_scalar)
            nsuccess = nsuccess + floor(success_rate_tmp_scalar + tol_success)
            j = j + 1;
        end
        success_rates(i,:) = success_rate_tmp;
        errs_X(i,:,:) = tmp_errs_X;
        errs_L(i,:,:) = tmp_errs_L;
        errs_S(i,:,:) = tmp_errs_S;
        errs_res(i,:,:) = tmp_errs_res;
    end
    
    save(sprintf('experiment.mat', node_id), ...
          'deltas', ...
          'rho_rs', ...
          'rho_ss', ...
          'success_rates', ...
          'errs_X', ...
          'errs_L', ...
          'errs_S', ...
          'errs_res', '-v7.3')
      
    delete(poolobj)
    finished = 1;
end

