clc, close all, clear all, addpath helper_functions;
%% load data 
[time, x] = load_sunspot_numbers();
%% train/test split
T = length(x);
x_test = x(time>=2020); test_time = time(time>=2020);
test_length = length(x_test);
train_length = T-test_length;
x_train = x(1:train_length); train_time = time(1:train_length);
%% figure 1
figure(1)
subplot(1,2,2)
plot(train_time, x_train, 'Displayname', 'train', 'color', 'k'), hold on;
plot(test_time, x_test, 'Displayname', 'train', 'linestyle', ':', 'color', 'k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = 12*11;    % hypothesized AR order 
L = 12*11;    % hypothesized MA order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sPCA
acf = estimate_acf(x_train, L);
covmat_spca = toeplitz(acf(1:L)/acf(1));
[U, val] = eig(covmat_spca);
singval = sqrt(diag(val));
threshold = optimal_SVHT_coef(L/(2*L-1), 0)*median(singval);
Ur = U(:, singval>threshold);
U_spca_weighted = U*sqrt(val);
%% estimate LRR
lrr_spca = estimate_lrr(Ur);
%% prediction
[pred_spca, ~] = predict_ts(x_train, lrr_spca, test_length);
%% figure 1
figure(1)
subplot(1,2,1)
plot(lrr_spca, 'Displayname', 'sPCA', 'color', 'r'), hold on;
subplot(1,2,2)
plot(test_time, pred_spca, 'Displayname', 'sPCA', 'color', 'r'), hold on;
%% figure 2
figure(2)
subplot(1,4,1)
imagesc(covmat_spca), axis equal, axis off;
subplot(1,4,2)
plot(U_spca_weighted(:, end-7:end))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SSA
X = hankel(x_train(1:L), x_train(L:end));
[U, S, ~] = svd(hankel(x_train(1:L), x_train(L:end)));
threshold = optimal_SVHT_coef(L/(2*L-1), 0)*median(diag(S));
Ur_ssa = U(:, diag(S)>threshold);
%% frame projection
ssa_proj = ssa_projection(x_train, Ur_ssa);                                % SSA, following [6]
% basic_ssa_rec = basic_ssa(x_train, L, 1:sum(diag(S)>threshold));         % vanilla SSA 
%% estimate LRR
lrr_ssa = estimate_lrr(Ur_ssa);
%% prediction
[~, ts_complete_ssa] = predict_ts(ssa_proj, lrr_ssa, test_length);
%% figure 1
figure(1)
subplot(1,2,1)
plot(lrr_ssa, 'Displayname', 'SSA', 'color', 'g', 'linestyle', '--');
subplot(1,2,2)
plot(time, ts_complete_ssa, 'Displayname', 'SSA', 'color', 'g', 'linestyle', '--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wiener
lrr_wiener = estimate_wiener(x_train(1:end-1), x_train(P+1:end), P);
[pred_wiener, ~] = predict_ts(x_train, lrr_wiener, test_length);
%% figure 1
figure(1)
subplot(1,2,1)
plot(lrr_wiener, 'Displayname', 'Wiener', 'color', 'b'), hold on;
subplot(1,2,2)
plot(test_time, pred_wiener, 'Displayname', 'Wiener', 'color', 'b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AR
acf = estimate_acf(x_train, P+1);
lrr_ar = estimate_ar(acf, P+1);
[pred_ar, ~] = predict_ts(x_train, lrr_ar, test_length);
%% figure 1
figure(1)
subplot(1,2,1)
plot(lrr_ar, 'Displayname', 'AR', 'color', 	'#EDB120', 'linestyle', ':', 'linewidth', 2);
subplot(1,2,2)
plot(test_time, pred_ar, 'Displayname', 'AR', 'color', '#EDB120', 'linestyle', ':', 'linewidth', 2);
xlim([2015, 2025]);
legend()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classical PCA
pca_data = reshape(x_train(1:L*floor(length(x_train)/L)), [L, floor(length(x_train)/L)]);
covmat_pca = pca_data*pca_data';
[U, val] = eig(covmat_pca);
U_pca_weighted = U*sqrt(val);
%% figure 2
figure(2)
subplot(1,4,4)
imagesc(covmat_pca), axis equal, axis off;
subplot(1,4,3)
plot(U_pca_weighted(:, end-7:end))
