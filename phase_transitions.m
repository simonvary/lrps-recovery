%% Code for Fig. 1 and Fig. A.6

addpath('func_rpca/AccAltProj_for_RPCA/')
addpath('func_rpca/')
addpath('func/')
addpath(genpath('func_rpca/SpaRCS/'))

rng(1)

%%
    
% EXPERIMENT SETUP
num_trials = 20;
mat_size = [100 100]; 
max_success = 5;
tol_success = 1e-2;
    
% Function to generate subsampling operators 
generate_sensing = @(n, p) generate_fjlt(n,p);
    
% Function to generate low-rank plus sparse matrices
c = 1;
generate_matrix = @(r,s) generate_lsmat1(mat_size(1), mat_size(2), r, s, c);

% Projection operator
matproj_ls = @(M, r, s, tol, L0, S0) matproj_ls_accaltproj(M, r, s, [], [], []);

% Settings for the recovery algorithm
opts = load_lsrec( 'lsrec_niht', [mat_size(1) mat_size(2)]);
opts.MAX_ITER = 300;
opts.tol_res = 1e-6;
opts.alpha = [];

% Function handle of the recovery algorithm with the settings above
lsrec = @ (b, A, aA, r, s) lsrec_niht(b, A, aA, r, s, matproj_ls, opts);

% Function handle for performing experiment for a specific (delta, rho_r, rho_s)
% using generate_matrix, generate sensing, lsrec and num_trials of times.
phase_single = @(delta, rho_r, rho_s) phase_single_general( ...
                      mat_size, delta, rho_r, rho_s, ...
                      num_trials, generate_matrix, generate_sensing, ...
                      lsrec, tol_success);
% PARAMETER RANGE
delta_from = 0.02;
delta_to = 1;
delta_by = 0.05;

rho_r_from = 0;
rho_r_to = 1;
rho_r_by = 0.1;

rho_s_from = 0;
rho_s_to = 1;
rho_s_by = 0.1;

deltas = delta_from:delta_by:delta_to;
[rho_ss, rho_rs] = meshgrid( rho_s_from:rho_s_by:rho_s_to, ...
                             rho_r_from:rho_r_by:rho_r_to);

ind_viable = (rho_rs + rho_ss <= 1) & (rho_rs > 0 | rho_ss > 0);
rho_rs = rho_rs(ind_viable);
rho_ss = rho_ss(ind_viable);

% Run the experiment for pairs in (rho_rs, rho_ss) with increasing delta
% until 5 success are reached.
finished = phase_node(num_cores, num_trials, max_success, tol_success, ...
    deltas, rho_rs, rho_ss, htc_single);

