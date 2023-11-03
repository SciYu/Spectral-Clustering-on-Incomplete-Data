function [X_norm] = normalize_x(X, preproess)
% function [X_norm] = normalize_x(X, low, high)
%
% To normalize data in X
%
% @param X         d * n (each column is a sample)
% @param low       Default "0"
% @param high      Default "1"
% 
% @return X_norm   Normalized data matrix

if (nargin < 2)
    preproess = 'nonormalize';
end

X = double(X);
n = size(X,2);

switch preproess
    case 'nonormalize'
        X_norm = X;
    case 'normalize1' 
        % normalize to [0, 1]
        X = X - repmat(min(X,[],2), 1, n);
        X_norm = X ./ (repmat(max(X,[],2), 1, n) + eps);
    case 'normalize2'
        % normalize to [-1, 1]
        low = -1; high = 1;
        X = X - repmat(min(X,[],2), 1, n);
        X = X ./ (repmat(max(X,[],2), 1, n) + eps);
        diff = 0.5 - 0.5*(low+high);
        X_norm = (X-diff) * (high-low);
end
end