clear all; warning off; clc; 
addpath(genpath(pwd));

%% Settings
dataset = 'Yale64';     % dataset
metric = 'Gaussian';    % kernel type (Gaussian or Exponential)
r = 0.8;                % missing ratio (missing completely at random)
niter = 5;              % number of repeated experiments
topk = 10;              % top-k items for recall calculation
seed = 2023;            % random seed
% Baselines of data imputation
baseline_impute = {'zero','mean','knn','em','svt','gr','kfmc'};
% Baselines of distance calibration
baseline_calibrate = {'dc','trf','ee','kc'};

%% Load a dataset
load([dataset, '_', num2str(r), 'miss.mat']);

%% Kernel and Distance Estimation
fprintf(['\n[Spectral Clustering]_', dataset, '_', metric, ', r = %1.1f, niter = %1.0f:'], r, niter);
fprintf('\nStage-I. Kernel and Distance Estimation');
for iter = 1 : niter
    fprintf('\n  iter = %1.0f: ', iter);
    X0 = X_list{iter};
    Y0 = Y_list{iter};
    idx = Idx_list{iter};

    %% Data Imputation 
    X_list{iter} = generate_x(X0, idx, baseline_impute);
    
    %% Distance Calibration & Kernel Correction
    [K_list{iter}, D_list{iter}] = generate_k(X_list{iter}, metric, baseline_calibrate);
end

%% Spectral Clustering
fprintf('\nStage-II. Spectral Clustering');
clear Score
for iter = 1 : niter
    fprintf('\n  iter = %1.0f: ', iter);
    num_cluster = length(unique(Y_list{iter}));
    
    %% SC: Standard Spectral Clustering
    rng(seed + iter);
    method = 'SC'; parameter.knn = 10;
    [S_list.SC{iter}, Label.SC{iter}] = spectral_cluster(K_list{iter}, num_cluster, method, parameter, iter);
    Score.SC{iter} = eval_cluster(Label.SC{iter}, Y_list{iter});
    
    %% KSSC: Kernel Sparse Subspace Clustering
    rng(seed + iter);
    method = 'KSSC'; parameter.lambda = 0.08; parameter.knn = 10;
    [S_list.KSSC{1,iter}, Label.KSSC{1,iter}] = spectral_cluster(K_list{1,iter}, num_cluster, method, parameter);
    Score.KSSC{1,iter} = eval_cluster(Label.KSSC{1,iter}, Y_list{1,iter});

    %% KLSR: Kernel Least-Squares Representation
    rng(seed + iter);
    method = 'KLSR'; parameter.lambda = 8; parameter.knn = 10;
    [S_list.KLSR{1,iter}, Label.KLSR{1,iter}] = spectral_cluster(K_list{1,iter}, num_cluster, method, parameter);
    Score.KLSR{1,iter} = eval_cluster(Label.KLSR{1,iter}, Y_list{1,iter});

   %% KSL-Sp: Kernel Self-expressive Learning with Schatten p-norm
    rng(seed + iter);
    method = 'KSL-Sp'; parameter.lambda = 0.6; parameter.knn = 10; 
    [S_list.KSLsp{1,iter}, Label.KSLsp{1,iter}] = spectral_cluster(K_list{1,iter}, num_cluster, method, parameter, iter);
    Score.KSLsp{1,iter} = eval_cluster(Label.KSLsp{1,iter}, Y_list{1,iter});

   %% AKLSR: Adaptive Kernel Least-Squares Representation
    rng(seed + iter);
    method = 'AKLSR'; parameter.lambda = 0.5; parameter.knn = 10; 
    [S_list.AKLSR{1,iter}, Label.AKLSR{1,iter}] = spectral_cluster(K_list{1,iter}, num_cluster, method, parameter, iter);
    Score.AKLSR{1,iter} = eval_cluster(Label.AKLSR{1,iter}, Y_list{1,iter});
end

%% Estimation Accuracy & Clustering Performance
[recall, recall_std] = eval_recall(D_list, 'distance', topk);
[error, error_std] = eval_error(D_list, K_list);
[Stat, Stat_std] = statistic_cluster(Score, error, error_std, recall, recall_std);
fprintf(['\n', '[Spectral Clustering]_', dataset, '_', metric, ', r = %1.1f, niter = %1.0f:\n'], r, niter);
disp(Stat([1:4,6:4:18,21:end], :));
    