function [A, label] = cluster_klsr(K, num_cluster, lambda, knn)
% function [A, label] = cluster_klsr(K, num_cluster, lambda, knn)
%
% Perform Kernel Least Squares Representation (KLSR) Clustering 
% on a given kernel matrix. (see reference)
%
% @param K              Kernel matrix
% @param num_cluster    Number of clusters
% @param lambda         Hyperparameter (default 8)
% @param knn            K-nearest neighbors of kNN graph (default 10)
% 
% @return A             Learned affinity matrix
% @return label         Cluster labels for all samples 
%
% <Reference>
% Fan Jicong, et al. "A simple approach to automated spectral clustering." 
% Advances in Neural Information Processing Systems, 2022, 35: 9907-9921.

if (nargin < 4)
    knn = 10;
end
if (nargin < 3)
    lambda = 8;
end

%% Step 1. Learn the affinity matrix C and A
n = size(K, 1);
I = eye(n, n);
C = real(inv(K + lambda*I) * K);  % C: affinity matrix
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
label = spectralcluster(A, num_cluster, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric', 'ClusterMethod', 'kmeans');

end