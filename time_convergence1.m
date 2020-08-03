%% Code for Fig. 2

addpath('func_rpca/AccAltProj_for_RPCA/')
addpath('func_rpca/')
addpath('func/')
addpath(genpath('func_rpca/SpaRCS/'))

rng(1)

%%
m = 100;
n = 100;

rng(0)

delta = 0.5;
rho_r = 0.1; % rho_r = 0.05, rho_r = 0.2
rho_s = rho_r;
rho = rho_r + rho_s;

p = round(delta*m*n);
s = round(rho_s*p);
r = round(((m+n)-sqrt( (m+n)^2 -4*rho_r*p))/2);

c = 1;
[L_true, S_true, X_true] = generate_lsmat1(m, n, r, s, c);

opts = load_lsrec( 'lsrec_naht', [m n]);
opts.L_true = L_true;
opts.S_true = S_true;
opts.MAX_ITER = 300;
opts.tol_res = 1e-4;
opts.tol_rel = 1-(1e-3);
opts.alpha_L = [];
opts.alpha_S = [];
opts.alpha = [];
opts.verb = 1;
opts.time_iter = 1;

[A, aA] = generate_gausst(m*n,p);

A_sparcs = @(Z) A(Z(:));
aA_sparcs = @(z) reshape(aA(z), m, n);

b = A(X_true(:));

matproj_ls1 = @(M, r, s, tol, L0, S0) matproj_ls_accaltproj(M, r, s, [], [], []);
matproj_ls2 = @(M, r, s, tol, L0, S0) matproj_ls_godec(M, r, s, [], [], []);

%%
[L1, S1, err_res, sparcs_err_X, sparcs_iter_timings] = sparcs(b, r, s, A_sparcs, aA_sparcs, 'svds', opts.tol_res , 30, 1,-Inf, X_true); 

[L, S, X, out_naht] = lsrec_naht(b, A, aA, r, s, opts);
[L, S, X, out_niht1] = lsrec_niht(b, A, aA, r, s, matproj_ls1, opts);
[L, S, X, out_niht2] = lsrec_niht(b, A, aA, r, s, matproj_ls2, opts);

%%
labfsz = 25;
legfsz = 20;
lwd = 1.5;

fig = figure;

semilogy(out_niht1.iter_timings, out_niht1.err_X, 'o-', 'LineWidth', lwd); hold on;
semilogy(out_niht2.iter_timings, out_niht2.err_X, 'o-', 'LineWidth', lwd); 
semilogy(out_naht.iter_timings, out_naht.err_X, 'o-', 'LineWidth', lwd);
semilogy([0; sparcs_iter_timings], [0.7; sparcs_err_X], 'o-', 'LineWidth', lwd); 
hold off;

ax = gca;
ax.FontSize = 18; 

xlim([0, 30])
ylim([0.5*0.0001, 1])
xlabel('Time (sec.)', 'FontSize', labfsz, 'Interpreter', 'LaTex');
ylabel('$$\| X^\ell - X_0\|_F\big/\|X_0\|_F$$', 'FontSize', labfsz, 'Interpreter', 'LaTex');

leg = legend('NIHT (AccAltProj)', 'NIHT (GoDec)', 'NAHT', 'SpaRCS');
set(leg, 'FontSize', legfsz, 'Interpreter', 'LaTex')

pbaspect([1.4 1 1])
