%% Code for Fig. 3

addpath('func_rpca/AccAltProj_for_RPCA/')
addpath('func_rpca/')
addpath('func/')
addpath(genpath('func_rpca/SpaRCS/'))

rng(1)

%% Setup and recovery

m = 100;
n = 100;

delta = 0.5;
rho_r = 0.05;
rho_s = 0.25; % rho_s = 0.15;
rho = rho_r + rho_s;

p = round(delta*m*n);
s = round(rho_s*p);
r = round(((m+n)-sqrt( (m+n)^2 -4*rho_r*p))/2);

c = 1;
[L_true, S_true, X_true] = generate_lsmat1(m, n, r, s, c);

opts = load_lsrec( 'lsrec_naht', [m n]);
opts.L_true = L_true;
opts.S_true = S_true;
opts.MAX_ITER = 500;
opts.tol_res = 1e-8;
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

[~, ~, ~, out_naht] = lsrec_naht(b, A, aA, r, s, opts);
[~, ~, ~, out_niht] = lsrec_niht(b, A, aA, r, s, matproj_ls1, opts);

%% Plot


ratio = 0.8;
labfsz = 25*ratio;
legfsz = 20*ratio;
ticfsz = 18*ratio;
lwd = 1.8;

fig = figure;

semilogy(out_niht.iter_timings, out_niht.err_L, '-o', 'LineWidth', lwd, 'color', [0, 0.4470, 0.7410]); hold on;
semilogy(out_niht.iter_timings, out_niht.err_S, '-o', 'LineWidth', lwd, 'color',[0.4940, 0.1840, 0.5560]);
semilogy(out_naht.iter_timings, out_naht.err_L, '-o', 'LineWidth', lwd, 'color', [0.8500, 0.3250, 0.0980]);
semilogy(out_naht.iter_timings, out_naht.err_S, '-o', 'LineWidth', lwd, 'color', [0.9290, 0.6940, 0.1250]);
hold off;

ax = gca;
ax.FontSize = ticfsz; 

%xlim([0, 30])
ylim([0.5*0.0001, 5])
%ylim([0.5*opts.tol_res, 1])
xlabel('Time (sec.)', 'FontSize', labfsz, 'Interpreter', 'LaTex');
ylabel('Rel. Frob. error', 'FontSize', labfsz, 'Interpreter', 'LaTex');


leg = legend('$$L$$, NIHT (AccAltProj)', '$$S$$, NIHT (AccAltProj)', '$$L$$, NAHT', '$$S$$, NAHT');
set(leg, 'FontSize', legfsz, 'Interpreter', 'LaTex')

pbaspect([1.6 1 1])

