function [A, label] = cluster_ksl_sp(K, num_cluster, lambda, knn, eval_tol, alpha_init)
% function [A, label] = cluster_ksl_sp(K, num_cluster, lambda, knn, eval_tol, alpha_init)
%
% Perform Kernel Self-expressive Learning with Schatten p-norm (KSL-Sp)
% on a given kernel matrix. (see Section 4.1)
%
% @param K              Kernel matrix
% @param num_cluster    Number of clusters
% @param lambda         Penalty coefficient (default 0.01)
% @param knn            K-nearest neighbors of kNN graph (default 10)
% @param eval_tol       Stop criterion for getting Connectivity matrix (default 1e-4)
% @param alpha_init     Step size
%  
% @return A             Learned affinity matrix
% @return label         Cluster labels for all samples
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

% Default values for the input arguments
if nargin < 2 || isempty(num_cluster), num_cluster = 10; end
if nargin < 3 || isempty(lambda), lambda = 0.01; end
if nargin < 4 || isempty(knn), knn = 10; end
if nargin < 5 || isempty(eval_tol), eval_tol = 1e-4; end
if nargin < 6 || isempty(alpha_init), alpha_init = 0.01; end

%% Step 1. Learn the affinity matrix A
A = get_affinity_matrix(K, eval_tol, alpha_init, lambda);

%% Step 2. Construct a kNN graph on A
n = size(K, 1);
[~, idx] = sort(A, 2, 'descend');
pos = zeros(n, n);
for i = 1:n
    pos(i, idx(i,1:knn)) = 1;
end
pos = ((pos + pos') > 0);
A(~pos) = 0;
A(1:n+1:n^2) = 0;

%% Step 3. Perform normalized cut
label = spectralcluster(A, num_cluster, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric', 'ClusterMethod', 'kmeans');

end


%%
function A = get_affinity_matrix(K, eval_tol, alpha_init, lambda)
    n = size(K, 1);
    C_old = zeros(n, n);
    [A, s, B] = svd(C_old);
    U = A * s;
    V = B;
    U_old = U;
    iter = 0;

    gradient_U = (2 * K * (U_old * V' - eye(n)) * V + lambda * U);
    flag1 = norm(gradient_U, 'fro');
    alpha = alpha_init;

    while flag1 > eval_tol && iter < 500
        iter = iter + 1;
        U = U_old - alpha * gradient_U;
        V = 2 * K * U / (2 * U' * K * U + lambda * eye(n));
        U_old = U;
        gradient_U = (2 * K * (U * V' - eye(n)) * V + lambda * U);
        flag1 = norm(gradient_U, 'fro');
    end

    A = U * V';
    A = real(A);
    A = (A + A') / 2;
end