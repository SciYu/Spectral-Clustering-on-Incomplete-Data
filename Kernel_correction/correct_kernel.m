function [K_new] = correct_kernel(K0, maxiter, low, high, type, k)
% function [K_new] = correct_kernel(K0, type, k, maxiter)
%
% Correct a kernel matrix by a kernel correction (KC) method (see Section 3.2)
%
% @param  K0        Initial kernel matrix (estimated from incomplete data)
% @param  maxiter   Maximum iterations (default 100)
% @param  low       Default 0 (lower bound of kernel values) 
% @param  high      Default 1 (upper bound of kernel values)
% @param  type      Default 'eig' (eigen-decomposition)
% @param  k         Hyperparameter in rSVD (randomized SVD)
%
% @return K_new     Corrected kernel matrix
%
% <Reference>
% Fangchen Yu, Runze Zhao, et al. "Boosting Spectral Clustering on Incomplete Data 
% via Kernel Correction and Affinity Learning", NeurIPS, 2023.

if (nargin < 2)
    maxiter = 100;
end
if (nargin < 3)
    low = 0;
end
if (nargin < 4)
    high = 1;
end
if (nargin < 5)
    type = 'eig'; % eigen-decomposition
end
if (nargin < 6)
    k = 10;
end

n = size(K0, 1);
K_new = nearpsd(K0, maxiter, low, high, type, k);

K_new = (K_new + K_new') / 2;
K_new(1:n+1:n^2) = 1;
K_new = max(real(K_new), 0);

end


%%
function [X, iter] = nearpsd(A, maxits, low, high, type, k)
% function [X, iter] = nearpsd(A, maxits, low, high, type, k)
%
% Computes the nearest positive semi-definite matrix 
% for a given square matrix.
%
% @param A        A square matrix to be calibrated
% @param maxits   Max num of iters allowed, default 100
% @param low      Default 0 
% @param high     Default 1
% @param type     Default 'eig' (eigen-decomposition)
% @param k        Hyperparameter in rSVD (randomized SVD)
%
% @return X       Nearest psd matrix to A
% @return iter    Number of iterations taken

if  ~isequal(A,A')
    A = (A + A') / 2;
end
if nargin < 6
    k = 10;
end
if nargin < 5
    type = 'eig';
end
if nargin < 4
    high = 1;
end
if nargin < 3
    low = 0; 
end
if nargin < 2
    maxits = 100;
end

% threshold for convergence & eigs
tolconv = 1.0e-5;
toleigs = 1.0e-5;

n = size(A,1);

U = zeros(n);
Y = A;

if strcmp(type, 'eig')
    [~, D, ~] = eig(Y);
elseif strcmp(type, 'svd')
    [~, D, ~] = svd(Y);
elseif strcmp(type, 'rsvd')
    [~, D, ~] = rsvd(Y, k); 
end
d = diag(D);

iter = 0;
while 1
    T = Y - U;

    % project onto psd matrices
    if strcmp(type, 'eig')
        [Q, D, ~] = eig(T);
        d = diag(D);
        p = d > toleigs*d(n);
        X = Q(:,p) * D(p,p) * Q(:,p)';
    elseif strcmp(type, 'svd')
        [Q, D, V] = svd(T);
        X = Q * max(D,0) * V';
    elseif strcmp(type, 'rsvd')
        [Q, D, V] = rsvd(T, k);
        X = Q * max(D,0) * V';
    end

    % update correction
    U = X - T;

    % maximum iteration & convergence test
    iter = iter + 1;
    if iter == maxits
        %fprintf('Max iterations reached. ');
        break; 
    end
    if norm(Y-X,'inf')/norm(Y,'inf') <= tolconv 
    	break;
    end
    
    % problem-dependent here
    Y = X;
    Y(1:n+1:n*n) = 1;
    Y(Y<low) = low;
    Y(Y>high) = high;
end

Y(1:n+1:n*n) = 1;
Y(Y<low) = low;
Y(Y>high) = high;

end
