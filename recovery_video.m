%% Code for Fig. 4

addpath('func_rpca/AccAltProj_for_RPCA/')
addpath('func_rpca/')
addpath('func/')

rng(1)

%% Load data


video_name = 'data/shop.m4v';

video = VideoReader(video_name);
video_tens = video.read();


selx = 1:190; % x-indices
sely = 7:146; % y-indices
selt = 1:150; % time-indices
video_tens = video_tens(sely,selx,1,selt);

video_size = size(video_tens);
video_size = video_size([1,2,4]);
clear video

video_mat = tens2mat(video_tens,4)';
video_mat = im2double(video_mat);

% implay(video_tens) % Check to play the video

[m n] = size(video_mat);


%%

% Set up for the NIHT solver
opts = load_lsrec( 'lsrec_niht', [m n]);
opts.X_true = video_mat;
opts.MAX_ITER = 500;
opts.tol_res = 1e-6;
opts.verb = 1;
opts.L_true = [];
opts.S_true = [];
opts.verb = 0;

% Take 0.33*m*n measurements
delta = 0.33;
p = round(delta*m*n);
alpha_s = 0.0495;
s = round(alpha_s*m*n);
r = 1;
rho_s = s/p;
rho_r = r*(m+n-r)/p;

[A, aA] = generate_fjlt(m*n,p);
b = A(video_mat(:));

% Function wrapper for Robust PCA via AccAltProj algorithm
matproj_ls = @(M, r, s, tol, L0, S0) matproj_ls_accaltproj(M, r, s, tol, [], []);
% Another option for Robust PCA is GoDec decomposition
% matproj_ls = @(M, r, s, tol, L0, S0) matproj_ls_godec(M, r, s, tol, [], []);

% Perform Robust PCA
[L_rpca, S_rpca, X_rpca, ~, ~, ~] = matproj_ls(video_mat, r, s, 1e-9, [], []);

% Perform NIHT
[L_niht, S_niht, X_niht, out_niht] = lsrec_niht(b, A, aA, r, s, matproj_ls, opts);

%% subsampled information

for i = 1:size(X_niht,2)
    psnrs_niht(i) = psnr(X_niht(:,i), video_mat(:,i));
end
mean(psnrs_niht)


for i = 1:size(X_rpca,2)
    psnrs_rpca(i) = psnr(X_rpca(:,i), video_mat(:,i));
end
mean(psnrs_rpca)

%% show results
L_tens = mat2tens(L_niht,size(video_tens),[1,2,3],4);
S_tens = mat2tens(S_niht,size(video_tens),[1,2,3],4);
X_niht_tens = mat2tens(X_niht,size(video_tens),[1,2,3],4);

L_rpca_tens = mat2tens(L_rpca,size(video_tens),[1,2,3],4);
S_rpca_tens = mat2tens(S_rpca,size(video_tens),[1,2,3],4);
X_rpca_tens = mat2tens(X_rpca,size(video_tens),[1,2,3],4);
dynamic_rpca_tens = X_rpca_tens;
dynamic_rpca_tens(S_rpca_tens==0) = 1;

dynamic_niht_tens = X_niht_tens;
dynamic_niht_tens(S_tens==0) = 1;

%% Plot

ind_frames = [15 60 80 130 150];
n_frames = length(ind_frames);

margins = [0.01, 0.01];
clim = [0, 1];
title_fsz = 14;

fig = figure();
colormap gray
for i = 1:n_frames
    subplot_tight(4, n_frames, i, margins);
    imagesc(X_rpca_tens(:,:,ind_frames(i)), clim)
    title(sprintf('Frame %d', ind_frames(i)), 'FontSize', title_fsz, 'Interpreter', 'LaTex')
    axis square
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([length(selx)/length(sely) 1 1])
end

for i = 1:n_frames
    subplot_tight(4, n_frames, i+n_frames, margins);
    imagesc(X_niht_tens(:,:,ind_frames(i)), clim)
    axis square
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([length(selx)/length(sely) 1 1])
end

for i = 1:n_frames
    subplot_tight(4, n_frames, i+2*n_frames, margins);
    imagesc(dynamic_rpca_tens(:,:,ind_frames(i)), clim)
    axis square
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([length(selx)/length(sely) 1 1])
end

for i = 1:n_frames
    subplot_tight(4, n_frames, i+3*n_frames, margins);
    imagesc(dynamic_niht_tens(:,:,ind_frames(i)), clim)
    axis square
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([length(selx)/length(sely) 1 1])
end


