function [A, label] = cluster_ksl_pp(K, num_cluster, lambda, knn, p, rho, gamma, alpha, maxiter, tol)
% function [A, label] = cluster_ksl_pp(K, num_cluster, lambda, knn, p, rho, gamma, alpha, maxiter, tol)
%
% Perform Kernel Self-expressive Learning with Proximal p-norm (KSL-Pp)
% on a given kernel matrix. (see Section 4.1)
%
% @param K                              Kernel matrix
% @param num_cluster                    Number of clusters
% @param lambda                         Penalty coefficient (default 0.01)
% @param knn                            K-nearest neighbors of kNN graph (default 10)
% @param p, rho, gamma, alpha           Hyperparameters
% @param maxiter                        Maximum number of iterations
% @param tol                            Tolerance threshold
% 
% @return A                             Learned affinity matrix
% @return label                         Cluster labels for all samples
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

% Default values for the input arguments
if nargin < 2 || isempty(num_cluster), num_cluster = 10; end
if nargin < 3 || isempty(lambda), lambda = 1; end
if nargin < 4 || isempty(knn), knn = 10; end
if nargin < 5 || isempty(p), p = 0.5; end
if nargin < 6 || isempty(rho), rho = 1; end
if nargin < 7 || isempty(gamma), gamma = 0.01; end
if nargin < 8 || isempty(alpha), alpha = 0.01; end
if nargin < 9 || isempty(maxiter), maxiter = 200; end
if nargin < 10 || isempty(tol), tol = 1e-3; end

%% Step 1. Learn the affinity matrix A
% Initialize variables
n = size(K, 1);
C = zeros(n);
Z = zeros(n);
U = zeros(n);

% Algorithm iterations
for t = 1:maxiter
    % C-update
    C = inv(2*K + (1/rho)*eye(n)) * (2*K + U + (1/rho)*Z);

    % Z-update
    for i = 1:n
        for j = 1:n
            c_ij = C(i, j);
            u_ij = U(i, j);

            % Solve for z_ij using solve_z function
            z_ij = solve_z_golden_section(lambda, p, rho, gamma, c_ij, u_ij, alpha, 100, 1e-6, Z(i, j));
            % z_ij = solve_z_bisection(lambda, p, rho, gamma, c_ij, u_ij, alpha, 100, 1e-6, Z(i, j));

            % Update Z
            Z(i, j) = z_ij;
        end
    end

    % U-update
    U = U + rho * (Z - C);

    % Check convergence
    if norm(C - Z, 'fro') <= tol && norm(Z - Z_prev, 'fro') <= tol
        break;
    end

    % Store previous Z for convergence check
    Z_prev = Z;
end

% Compute A
C = real(C);
A = (abs(C) + abs(C')) / 2;

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