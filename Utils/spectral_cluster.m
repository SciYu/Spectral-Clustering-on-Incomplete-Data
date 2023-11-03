function [S_list, Label] = spectral_cluster(K_list, num_cluster, method, parameter, iter)
% function [S_list, Label] = spectral_cluster(K_list, num_cluster, method, parameter, iter)
%
% @param K_list          Approximate kernel matrices using baselines
% @param num_cluster     Number of clusters
% @param method          Spectral clustering method
% @param parameter       Parameters for clustering methods
% 
% @return label          Predicted labels using clustering method

if (nargin < 4)
    parameter.knn = 10;
end
if (nargin < 3)
    method = 'SC';
end

fprintf(['\n  ',method,': ']);
baseline = fieldnames(K_list);
for i = 1 : length(baseline)
    name = baseline{i};
    fprintf([name, ', ']);
    for j = 1 : 2
        eval(['[S_list.',name,', Label.',name, '(:,j)] = cluster_kernel(K_list.', name,', num_cluster, method, parameter);']); 
    end
end

end


%%
function [S, label] = cluster_kernel(K, num_cluster, method, parameter)

if (nargin < 4)
    parameter.knn = 10;
end
if (nargin < 3)
    method = 'SC';
end

K = real(K);
if strcmp(method, 'SC')
    % Standard Spectral Clustering
    S = make_kNN_graph(K, parameter.knn);
    try
        label = spectralcluster(S, num_cluster, 'Distance', 'precomputed', 'LaplacianNormalization', 'symmetric', 'ClusterMethod', 'kmeans');
    catch
        label = spectralcluster(S, num_cluster);
    end
elseif strcmp(method, 'KSSC')
    % Kernel Sparse Subspace Clustering
    lambda = parameter.lambda; % 0.1
    knn = parameter.knn; % 10
    [S, label] = cluster_kssc(K, num_cluster, lambda, knn);
elseif strcmp(method, 'KLSR')
    % Kernel Least Squares Representation
    lambda = parameter.lambda; % 8 or 25
    knn = parameter.knn; % 10
    [S, label] = cluster_klsr(K, num_cluster, lambda, knn);
elseif strcmp(method, 'KSL-Sp')
    % Kernel Self-expressive Learning with Schatten p-norm
    lambda = parameter.lambda; % 1
    knn = parameter.knn; % 10
    [S, label] = cluster_ksl_sp(K, num_cluster, lambda, knn);
elseif strcmp(method, 'AKLSR')
    % Adaptive Kernel Least-Squares Representation
    lambda = parameter.lambda; % 0.5
    knn = parameter.knn; % 10
    [S, label] = cluster_aklsr(K, num_cluster, lambda, knn);
end

end