%% Code for Fig. 5

addpath('func_rpca/AccAltProj_for_RPCA/')
addpath('func_rpca/')
addpath('func/')

rng(1)

%% Load data

load('data/grss18-building512.mat')

Y_true = image;

image_size = size(Y_true);

Y_true = tens2mat(Y_true, [], 3);
Y_true = Y_true ./ vecnorm(Y_true);

alpha_s = 0.013;
s = round(alpha_s*numel(Y_true));
r = 3;

m = size(Y_true,1);
n = size(Y_true,2);


[L, S, Y_ls, U, V, T] = matproj_ls_accaltproj(Y_true, r, s, 1e-5, [], []);

%% NIHT

delta = 0.33;
p = round(prod(image_size)*delta);

rho_r = r*(m+n - r)/p
rho_s = s/p


[A, aA] = generate_fjlt(prod(image_size), round(prod(image_size)*delta));
b = A(Y_true(:));

opts = load_lsrec( 'lsrec_niht', size(Y_true));
opts.L_true = [];
opts.S_true = [];
opts.MAX_ITER = 300;
opts.tol_res = 1e-3;
opts.verb = 1;

matproj_ls = @(M, r, s, tol, L0, S0) matproj_ls_accaltproj(M, r, s, tol, [], []);

% NIHT using low-rank plus sparse
[L_niht, S_niht, Y_niht, out_niht] = lsrec_niht(b, A, aA, r, s, matproj_ls, opts);

% NIHT using only rank
[L_mc, S_mc, Y_mc, out_mc] = lsrec_niht(b, A, aA, r, 0, matproj_ls, opts);



%%
Y_niht_tens = mat2tens(Y_niht, image_size, [], 3);
Y_mc_tens = mat2tens(Y_mc, image_size, [], 3);
Y_true_tens = mat2tens(Y_true, image_size, [], 3);

% Mean PSNR
mean(psnr_spec(Y_mc_tens, Y_true_tens))
mean(psnr_spec(Y_niht_tens, Y_true_tens))

% Spatial PSNR
psnr_niht = psnr_spat(Y_niht_tens, Y_true_tens);
psnr_mc = psnr_spat(Y_mc_tens, Y_true_tens);

%% Draw reconstruction

r_band = 15;
g_band = 12;
b_band = 5 ;
fprintf('RGB rendering using wavelengths: %3.1f, %3.1f, %3.1f nm.\n', wavelengths([r_band, g_band, b_band]))


clim = [0, 1];
img_true = Y_true_tens(:,:,[r_band g_band b_band]);
img_niht = Y_niht_tens(:,:,[r_band g_band b_band]);
img_mc = Y_mc_tens(:,:,[r_band g_band b_band]);

% Normalize
max_true = 0.9*reshape(max(max(img_true(:,:,:))),[],1);
for i=1:3
    img_true(:,:,i) = img_true(:,:,i) ./ max_true(i);
    img_niht(:,:,i) = img_niht(:,:,i) ./ max_true(i);
    img_mc(:,:,i)   = img_mc(:,:,i)   ./ max_true(i);
end

margins = [0.01, 0.01];
title_fsz = 14;
fig = figure;

subplot_tight(1, 3, 1, margins);
imagesc(img_true(:,:,:),clim)
title('Groundtruth',  'FontSize', title_fsz, 'Interpreter', 'LaTex')
axis equal
axis off

subplot_tight(1, 3, 2, margins);
imagesc(img_niht(:,:,:),clim)
title('Low-rank plus sparse',  'FontSize', title_fsz, 'Interpreter', 'LaTex')
axis equal
axis off

subplot_tight(1, 3, 3, margins);
imagesc(img_mc(:,:,:),clim)
title('Low-rank',  'FontSize', title_fsz, 'Interpreter', 'LaTex')
axis equal
axis off

%% Draw spatial PSNR
sel1_x = (129:257)+128;
sel1_y = 385:512;
sel2_x = (129:256)-32;
sel2_y = 1:128;
lwd = 3;
clim = [13, 38];

fig = figure;

subplot_tight(2, 1, 1, margins)
imagesc(psnr_mc, clim)
hold on;
rectangle('Position',[min(sel1_y) min(sel1_x) 127 127], 'LineWidth', lwd);
rectangle('Position',[min(sel2_y) min(sel2_x) 127 127], 'LineWidth', lwd);
hold off;
colormap(flipud(jet))
c = colorbar;
c.TickLabelInterpreter = 'LaTex';
c.Ticks = [15, 20, 25, 30, 35];
c.FontSize = 18;
axis off
axis equal
box on


subplot_tight(2, 1, 2, margins)
imagesc(psnr_niht, clim)
hold on;
rectangle('Position',[min(sel1_y) min(sel1_x) 127 127],  'LineWidth', lwd);
rectangle('Position',[min(sel2_y) min(sel2_x) 127 127], 'LineWidth', lwd);
hold off;
colormap(flipud(jet))
c = colorbar;
c.TickLabelInterpreter = 'LaTex';
c.Ticks = [15, 20, 25, 30, 35];
c.FontSize = 18;
axis off
axis equal
box on


%% Draw Detail 1 & 2

sel1_x = (129:256)-32;
sel1_y = 1:128;

sel2_x = (129:257)+128;
sel2_y = 385:512;

color_pal = flipud(brewermap(30,'Spectral'));

sel_band = 23;
wavelengths(sel_band)

max_true = 0.8*max(max(Y_true_tens(sel1_x,sel1_y,sel_band)));
min_true = 1.1*min(min(Y_true_tens(sel1_x,sel1_y,sel_band)));
clim = [min_true, max_true];

fig = figure;

subplot_tight(3, 2, 1, margins)
imagesc(Y_true_tens(sel1_x,sel1_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on

subplot_tight(3, 2, 2, margins)
imagesc(Y_true_tens(sel2_x,sel2_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on

subplot_tight(3, 2, 3, margins)
imagesc(Y_niht_tens(sel1_x,sel1_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on

subplot_tight(3, 2, 4, margins)
imagesc(Y_niht_tens(sel2_x,sel2_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on

subplot_tight(3, 2, 5, margins)
imagesc(Y_mc_tens(sel1_x,sel1_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on

subplot_tight(3, 2, 6, margins)
imagesc(Y_mc_tens(sel2_x,sel2_y, sel_band),clim)
axis off;axis equal;
colormap(color_pal)
box on
