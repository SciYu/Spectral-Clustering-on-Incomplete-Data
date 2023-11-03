function S_knn = make_kNN_graph(S, knn)
% function S_knn = make_kNN_graph(S, knn)
%
% Construct a kNN graph of a given affinity matrix.
%
% @param S       Affinity matrix
% @param knn     K-nearest neighbors (default 10)
%
% @return S_knn  kNN affinity matrix

if nargin < 2
    knn = 10;
end

n = size(S, 1);
[~, idx] = sort(S, 2, 'descend');
pos = zeros(n, n);
for i = 1:n
    pos(i, idx(i,1:knn+1)) = 1;
end
pos = ((pos + pos') > 0);
S(~pos) = 0;
S(1:n+1:n^2) = 0;
S_knn = S;

end