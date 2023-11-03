function [A, label] = cluster_kssc(K, num_cluster, lambda, knn, maxiter)
% function [A, label] = cluster_kssc(K, num_cluster, lambda, knn, maxiter)
%
% Perform Kernel Sparse Subspace Clustering (KSSC) on a given kernel matrix.
%
% @param K              Kernel matrix
% @param num_cluster    Number of clusters
% @param lambda         Hyperparameter (default 0.1)
% @param knn            K-nearest neighbors of kNN graph (default 10)
% @param maxiter        Maximum number of iterations (default 100)
% 
% @return A             Learned affinity matrix
% @return label         Cluster labels for all samples 
%
% <Reference>
% Vishal M Patel and René Vidal. "Kernel sparse subspace clustering." 
% IEEE International Conference on Image Processing, pages 2849–2853. IEEE, 2014.

if (nargin < 5)
    maxiter = 100;
end
if (nargin < 4)
    knn = 10;
end
if (nargin < 3)
    lambda = 0.1;
end

%% Step 1. Learn the affinity matrix C and A
n = size(K, 1);
I = eye(n, n);
C = real(inv(K + I) * K);  % C: affinity matrix
A = C;
Y = 0;
u = 0.1; % 0.01 for Gaussian kernel; 0.1 for other kernels
e = 1e-4;
iter = 0;
while iter < maxiter
    iter = iter + 1;
    % A
    A_new = inv(K+u*eye(n))*(K+u*C-Y);
    % C
    temp = A_new + Y/u;
    C_new = max(0,temp-lambda/u) + min(0,temp+lambda/u);
    C_new = C_new - diag(diag(C_new));
    Y = Y + u*(A_new-C_new);
    % Stopping criteria
    stopC = max(norm(C_new-C,'fro')/norm(C,'fro'), norm(A_new-A,'fro')/norm(A,'fro'));
    isstopC = stopC < e;
    if isstopC
        break;
    end
    % Update
    C = C_new;    
    A = A_new;
end
C(isnan(C)) = 0;
C(1:n+1:n^2) = 0;
C = abs(real(C));
A = C + C';  % symmetric

%% Step 2. Construct a kNN graph on A
[~, idx] = sort(A, 2, 'descend');
pos = zeros(n, n);
for i = 1:n
    pos(i, idx(i,1:knn)) = 1;
end
pos = ((pos + pos') > 0);
A(~pos) = 0;
A(1:n+1:n^2) = 0;

%% Step 3. Perform normalized cut
try
    label = spectralcluster(real(A), num_cluster, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric', 'ClusterMethod', 'kmeans');
catch
    label = spectralcluster(real(A), num_cluster);
end

end