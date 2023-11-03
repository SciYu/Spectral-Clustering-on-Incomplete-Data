function [D_new] = correct_distance(D0, maxiter, type, k)
% function [D_new] = correct_distance(D0, maxiter, type, k)
%
% Correct a distance matrix by a kernel correction method (see Sections 3,2 and 3.4)
%
% @param  D0        Initial distance matrix (estimated from incomplete data)
% @param  maxiter   Maximum iterations (default 100)
% @param  type      Default 'eig' (eigen-decomposition)
% @param  k         Hyperparameter in rSVD (randomized SVD)
%
% @return D_new     Corrected distance matrix
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

if (nargin < 4)
    k = 10;
end
if (nargin < 3)
    type = 'eig'; % eigen-decomposition
end
if (nargin < 2)
    maxiter = 100;
end

n = size(D0,1);
low = 0; high = 1;

% Calculate a Gaussian kernel K0 from D0
sigma = median(D0(:));
K0 = exp(-D0.^2/(sigma^2));

% Correct K0 to K_new
K_new = correct_kernel(K0, maxiter, low, high, type, k);

% Convert K_new back to D_new
D_new = sqrt(log(K_new) .* (-sigma^2));

D_new = (D_new + D_new') / 2;
D_new(1:n+1:n^2) = 0;
D_new = max(real(D_new), 0);

end