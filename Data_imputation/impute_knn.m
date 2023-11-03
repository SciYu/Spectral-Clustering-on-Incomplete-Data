function [imputedX] = impute_knn(Xmiss, k)
% function [imputedX] = impute_knn(Xmiss, k)
%
% Impute a data matrix. Each column is a sample. Each NaN value is replaced
% by the mean of the sample's k-nearest neighbors with known values. If all 
% k-nearest samples' corresponding features are NaN, then replaced by zero.
%
% @param Xmiss      Incomplete data matrix, each column is a sample
% @param k          Default 10 (k-nearest neighbors)
% 
% @return imputedX  Imputed matrix with all data samples
%
% <Reference>
% Lorenzo Beretta and Alessandro Santaniello. "Nearest neighbor imputation algorithms: a critical
% evaluation." BMC Medical Informatics and Decision Making, 16(3):197–208, 2016.

if (nargin < 2)
    k = 10;
end
low = min(Xmiss(:)); high = max(Xmiss(:));

n = size(Xmiss, 2);
idx = isnan(Xmiss);
Xmiss(idx) = 0;

Dmiss = distance(Xmiss, 'euclidean');
Dmiss(1:n+1:n^2) = inf;
[~, Imiss] = sort(Dmiss, 2, 'ascend');
Ximp = zeros(size(Xmiss));
for i = 1 : n
    Ximp(:,i) = nanmean(Xmiss(:, Imiss(i,1:k)), 2); % row mean
end
Ximp(isnan(Ximp)) = 0;
Ximp(~idx) = 0;
imputedX = Xmiss + Ximp;

imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;

end