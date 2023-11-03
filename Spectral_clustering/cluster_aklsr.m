function [A, label] = cluster_aklsr(K, num_cluster, lambda, knn, rho, maxiter, tol)
% function [A, label] = cluster_aklsr(K, num_cluster, lambda, knn, rho, maxiter, tol)
%
% Perform Adaptive Kernel Least-Squares Representation (AKLSR)
% on a given kernel matrix. (see Section 4.2)
%
% @param K             Initial kernel matrix
% @param num_cluster   Number of clusters
% @param lambda        Regularization parameter
% @param knn           K-nearest neighbors of kNN graph (default 10)
% @param rho           Penalty parameter for ADMM
% @param maxiter       Maximum number of iterations
% @param tol           Tolerance for convergence
% 
% @return C            Learned affinity matrix
% @return label        Cluster labels for all samples
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

% Default values for the input arguments
if nargin < 2 || isempty(num_cluster), num_cluster = 10; end
if nargin < 3 || isempty(lambda), lambda = 0.1; end
if nargin < 4 || isempty(knn), knn = 10; end
if nargin < 5 || isempty(rho), rho = 0.01; end
if nargin < 6 || isempty(maxiter), maxiter = 200; end
if nargin < 7 || isempty(tol), tol = 1e-3; end

%% Step 1.Learn the affinity matrix C
% Get the size of K0
n = size(K, 1);

% Initialize variables
K0 = K;
A = zeros(n, n);
C = zeros(n, n);
U = zeros(n, n);
I = eye(n, n);

% Main loop for ADMM iterations
for iter = 1:maxiter
    % Update K
    K = (2 * K0 - U + rho * A) / (rho + 2);
    
    % Update K to PSD
    [Q, D] = eig(K);
    d = diag(D);
    p = d > 1e-5*d(n);
    K = Q(:,p) * D(p,p) * Q(:,p)';

    % Update A
    A = (rho * K + U + 2 * C' - C * C' - I) / rho;

    % Update C
    C = real(inv(2 * lambda * I + A + A')) * A'*2;

    % Update U
    U = U + rho * (K - A);

    % Calculate primal and dual residuals to check convergence
    primal_residual = norm(K - A, 'fro')^2;
    dual_residual = norm(rho * (K - A), 'fro')^2;

    % Check if the residuals are lower than the tolerance
    if primal_residual < tol && dual_residual < tol
        break;
    end
end

C = real(C);
C = (abs(C) + abs(C)')/2;

%% Step 2. Construct a kNN graph on C
[~, idx] = sort(C, 2, 'descend');
pos = zeros(n, n);
for i = 1:n
    pos(i, idx(i,1:knn)) = 1;
end
pos = ((pos + pos') > 0);
C(~pos) = 0;
C(1:n+1:n^2) = 0;

%% Step 3. Perform normalized cut
label = spectralcluster(C, num_cluster, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric', 'ClusterMethod', 'kmeans');

end