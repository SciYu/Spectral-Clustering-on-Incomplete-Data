function [imputedX] = impute_kfmc(Xmiss, method)
% function [imputedX] = impute_kfmc(Xmiss, method)
%
% Impute a data matrix. Each column is a sample. The complete matrix is obtained
% by the Kernelized Factorization Matrix Completion (KFMC) method (see reference).
%
% @param Xmiss     Incomplete data matrix, each column is a sample
% @param method    Default 'poly' (Polynomial kernel function)
% 
% @return imputedX Imputed matrix with all data samples
%
% <Reference>
% Fan, Jicong, and Madeleine Udell. "Online high rank matrix completion." 
% Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition. 2019.

if (nargin < 2)
    method = 'poly';
end
low = min(Xmiss(:)); high = max(Xmiss(:));

Xzero = Xmiss;
Xzero(isnan(Xzero)) = 0;
M = double(~isnan(Xmiss));

if strcmp(method, 'poly')
    ker.type = 'poly'; ker.par = [1 2];
    alpha = 1; beta = 1; d = 100;
    [imputedX, ~] = KFMC(Xzero, M, d, alpha, beta, ker);
elseif strcmp(method, 'rbf')
    ker.type = 'rbf'; ker.par = []; ker.par_c = 1;
    alpha = 0; beta = 0.001; d = 100;
    [imputedX, ~] = KFMC(Xzero, M, d, alpha, beta, ker);
end

imputedX(imputedX<low) = low;
imputedX(imputedX>high) = high;

end