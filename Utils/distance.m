function [D] = distance(X, metric, type)
% function [D] = distance(X, metric, type)
%
% Approximate a distance matrix for samples with NaN values. (see Section 2.2)
%
% @param X        d*n, each column is a sample
% @param metric   Default "euclidean"
% @param type     Default "true" (calculate the true distance matrix)
% @return D       n*n distance matrix
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

if (nargin < 3)
    type = 'true';
end
if (nargin < 2)
    metric = 'euclidean';
end

if strcmp(metric, 'euclidean')
    switch type
        case 'true'
            D = pdist2(X', X', 'euclidean');
        case 'miss'
            [d, n] = size(X);
            D = zeros(n);
            for i = 1 : n
                x = X(:,i);
                xx = X - x;
                O = ~isnan(xx);
                D(i,:) = sqrt(nansum(xx.^2)) .* sqrt(d./(sum(O)+eps));
            end
            D(1:n+1:n^2) = 0;
        case 'basic'
            [d, n] = size(X);
            O = ~isnan(X);
            D = zeros(n);
            for i = 1 : n-1
                for j = i+1 : n
                    ii = O(:,i) & O(:,j);
                    D(i,j) = norm(X(ii,i)-X(ii,j)) * sqrt(d/(nnz(ii)+eps));
                end
            end
            D(1:n+1:n^2) = 0;
            D = real(D);
            D = D + D';
    end
end
D = max(D, 0);
end